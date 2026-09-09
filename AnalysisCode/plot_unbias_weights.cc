// plot_unbias_weights.cc
// Quick ROOT macro to verify unbiasing by comparing unweighted vs weighted pT.
//
// Builds: c++ -O2 plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights
//
// Example:
//   ./plot_unbias_weights \
//     --input unbias_weights.root \
//     --out unbias_weights_plots.root \
//     --pt-min 50 --pt-max 200 --nbins 60

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TMath.h>
#include <TLine.h>

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

static void die(const std::string &msg) {
  std::cerr << "error: " << msg << std::endl;
  std::exit(1);
}

static std::string get_arg(int argc, char* argv[], const std::string &flag, const std::string &def = "") {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == flag && i + 1 < argc) return argv[i + 1];
  }
  return def;
}

// =====================================================================
// BASIS VECTORS
// =====================================================================
// This macro no longer re-derives the basis.  unbias_weights.cc writes the
// per-jet basis values g_k directly into its output — source jets in
// tweights.g_basis, target jets in the tbasis_target tree, with per-function
// titles in meta.basis_labels — so the basis-vector plots below always match
// exactly the basis that was fit (EEC terms, several subjet radii, custom
// powers, ...), with no reclustering and no dependence on the old fixed basis.


static void fill_eec(TH1D *h, const std::vector<double> &cpt,
                     const std::vector<double> &ceta,
                     const std::vector<double> &cphi,
                     double wjet) {
  const size_t nconst = cpt.size();
  if (nconst < 2) return;
  const double norm = 1.0 / (120.0 * 120.0);   // z = pt_i pt_k / 120^2
  for (size_t i = 0; i < nconst; ++i) {
    if (cpt[i] <= 0.0) continue;
    for (size_t k = i + 1; k < nconst; ++k) {
      if (cpt[k] <= 0.0) continue;
      double dphi = TVector2::Phi_mpi_pi(cphi[i] - cphi[k]);
      double deta = ceta[i] - ceta[k];
      double theta = std::sqrt(deta * deta + dphi * dphi);
      double w = (cpt[i] * cpt[k]) * norm;
      h->Fill(theta, wjet * w);
    }
  }
}

int main(int argc, char* argv[]) {
  const std::string input = get_arg(argc, argv, "--input");
  const std::string out_name = get_arg(argc, argv, "--out", "unbias_weights_plots.root");
  const std::string target_input = get_arg(argc, argv, "--target-input");
  const std::string target_tree = get_arg(argc, argv, "--target-tree", "tX");
  const double pt_min = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
  const double pt_max = std::stod(get_arg(argc, argv, "--pt-max", "200.0"));
  const int nbins = std::stoi(get_arg(argc, argv, "--nbins", "80"));
  const int eec_bins = std::stoi(get_arg(argc, argv, "--eec-bins", "60"));
  const double eec_max = std::stod(get_arg(argc, argv, "--eec-max", "0.4"));

  if (input.empty()) die("--input is required");
  if (target_input.empty()) die("--target-input is required to compare against unbiased target");


  TFile in(input.c_str(), "READ");
  if (in.IsZombie()) die("failed to open input file: " + input);
  TTree *t = dynamic_cast<TTree*>(in.Get("tweights"));
  if (!t) die("missing tree 'tweights' in input file");

  
  
  float pt = 0.0f;
  double w_unbias = 1.0;
  double w_base = 1.0;
  double w_total = 1.0;
  std::vector<double> *const_pt = nullptr;
  std::vector<double> *const_eta = nullptr;
  std::vector<double> *const_phi = nullptr;

  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("w_unbias", &w_unbias);
  t->SetBranchAddress("w_base", &w_base);
  t->SetBranchAddress("w_total", &w_total);
  const bool has_const = (t->GetBranch("const_pt") != nullptr &&
                          t->GetBranch("const_eta") != nullptr &&
                          t->GetBranch("const_phi") != nullptr);
  if (has_const) {
    t->SetBranchAddress("const_pt", &const_pt);
    t->SetBranchAddress("const_eta", &const_eta);
    t->SetBranchAddress("const_phi", &const_phi);
  }

  TH1D *h_base = new TH1D("h_base", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_total = new TH1D("h_total", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_target = new TH1D("h_target", "pT; p_{T} [GeV]; entries", nbins, pt_min, pt_max);
  TH1D *h_w_unbias = new TH1D("h_w_unbias", "unbias weight; w; entries", 80, 0.0, 2.0);
  // log-spaced theta binning for the EEC histograms
  const double eec_logmin = TMath::Log10(0.004);
  const double eec_logmax = TMath::Log10(eec_max);
  const double eec_binwidth = (eec_logmax - eec_logmin) / eec_bins;
  std::vector<double> eec_edges(eec_bins + 1);
  for (int i = 0; i <= eec_bins; ++i)
    eec_edges[i] = std::pow(10.0, eec_logmin + i * eec_binwidth);

  TH1D *h_eec_base   = new TH1D("h_eec_base",   "EEC; #theta; EEC", eec_bins, eec_edges.data());
  TH1D *h_eec_total  = new TH1D("h_eec_total",  "EEC; #theta; EEC", eec_bins, eec_edges.data());
  TH1D *h_eec_target = new TH1D("h_eec_target", "EEC; #theta; EEC", eec_bins, eec_edges.data());

  const Long64_t nentries = t->GetEntries();
  double sumw_base = 0.0;
  double sumw_total = 0.0;
  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    h_base->Fill(pt, w_base);
    h_total->Fill(pt, w_total);
    h_w_unbias->Fill(w_unbias);
    if (pt >= pt_min && pt <= pt_max && has_const && const_pt && const_eta && const_phi) {
      if (const_pt->size() == const_eta->size() && const_pt->size() == const_phi->size()) {
        fill_eec(h_eec_base, *const_pt, *const_eta, *const_phi, w_base);
        fill_eec(h_eec_total, *const_pt, *const_eta, *const_phi, w_total);
        sumw_base += w_base;
        sumw_total += w_total;
      }
    }
  }

  // Load unbiased target tree (tX) for comparison
  TFile tfin(target_input.c_str(), "READ");
  if (tfin.IsZombie()) die("failed to open target input file: " + target_input);
  TTree *tt = dynamic_cast<TTree*>(tfin.Get(target_tree.c_str()));
  if (!tt) die("missing target tree: " + target_tree);

  float tpt = 0.0f;
  float tw = 1.0f;
  std::vector<double> *tconst_pt = nullptr;
  std::vector<double> *tconst_eta = nullptr;
  std::vector<double> *tconst_phi = nullptr;
  tt->SetBranchAddress("pt", &tpt);
  tt->SetBranchAddress("weight", &tw);
  const bool target_has_const = (tt->GetBranch("const_pt") != nullptr &&
                                 tt->GetBranch("const_eta") != nullptr &&
                                 tt->GetBranch("const_phi") != nullptr);
  if (target_has_const) {
    tt->SetBranchAddress("const_pt", &tconst_pt);
    tt->SetBranchAddress("const_eta", &tconst_eta);
    tt->SetBranchAddress("const_phi", &tconst_phi);
  }
  const Long64_t tentries = tt->GetEntries();
  double sumw_target = 0.0;
  for (Long64_t i = 0; i < tentries; ++i) {
    tt->GetEntry(i);
    h_target->Fill(tpt, tw);
    if (tpt >= pt_min && tpt <= pt_max && target_has_const && tconst_pt && tconst_eta && tconst_phi) {
      if (tconst_pt->size() == tconst_eta->size() && tconst_pt->size() == tconst_phi->size()) {
        fill_eec(h_eec_target, *tconst_pt, *tconst_eta, *tconst_phi, tw);
        sumw_target += tw;
      }
    }
  }

  h_base->SetLineColor(kBlue + 1);
  h_total->SetLineColor(kRed + 1);
  h_target->SetLineColor(kGreen + 2);


  h_base->SetLineWidth(2);
  h_total->SetLineWidth(3);
  h_target->SetLineWidth(3);
  h_eec_base->SetLineColor(kBlue + 1);
  h_eec_total->SetLineColor(kRed + 1);
  h_eec_target->SetLineColor(kGreen + 2);
  h_eec_base->SetLineWidth(2);
  h_eec_total->SetLineWidth(3);
  h_eec_target->SetLineWidth(3);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");

  TCanvas *c = new TCanvas("c_unbias", "unbias check", 900, 900);
  TPad *p_top = new TPad("p_pt_top", "", 0.0, 0.30, 1.0, 1.0);
  TPad *p_bot = new TPad("p_pt_bot", "", 0.0, 0.0, 1.0, 0.30);
  p_top->SetBottomMargin(0.0);
  p_top->SetLeftMargin(0.12);
  p_top->SetRightMargin(0.04);
  p_bot->SetTopMargin(0.0);
  p_bot->SetBottomMargin(0.30);
  p_bot->SetLeftMargin(0.12);
  p_bot->SetRightMargin(0.04);
  p_top->Draw();
  p_bot->Draw();

  p_top->cd();
  gPad->SetLogy();
  gPad->SetTicks(1,1); 
  h_base->Draw("hist");
  h_total->Draw("hist same");
  h_target->Draw("hist same");

  TLegend *leg = new TLegend(0.52, 0.7, 0.88, 0.9);
  leg->AddEntry(h_base, "biased sample (base weight)", "l");
  leg->AddEntry(h_total, "weighted sample (base x unbias)", "l");
  leg->AddEntry(h_target, Form("unbiased target %s", target_tree.c_str()), "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  p_bot->cd();
  gPad->SetTicks(1,1);
  TH1D *h_ratio_base = (TH1D*)h_base->Clone("h_ratio_base_pt");
  TH1D *h_ratio_total = (TH1D*)h_total->Clone("h_ratio_total_pt");
  h_ratio_base->Divide(h_target);
  h_ratio_total->Divide(h_target);
  h_ratio_total->SetTitle("; p_{T} [GeV]; ratio to target");
  h_ratio_total->SetLineColor(kRed + 1);
  h_ratio_base->SetLineColor(kBlue + 1);
  h_ratio_total->SetLineWidth(2);
  h_ratio_base->SetLineWidth(2);
  h_ratio_total->SetMinimum(0.5);
  h_ratio_total->SetMaximum(1.5);
  h_ratio_total->GetYaxis()->SetNdivisions(505);
  h_ratio_total->GetYaxis()->SetTitleSize(0.10);
  h_ratio_total->GetYaxis()->SetTitleOffset(0.5);
  h_ratio_total->GetYaxis()->SetLabelSize(0.08);
  h_ratio_total->GetXaxis()->SetTitleSize(0.10);
  h_ratio_total->GetXaxis()->SetLabelSize(0.08);
  h_ratio_total->Draw("hist");
  h_ratio_base->Draw("hist same");
  TLine *lr_pt = new TLine(h_ratio_total->GetXaxis()->GetXmin(), 1.0,
                           h_ratio_total->GetXaxis()->GetXmax(), 1.0);
  lr_pt->SetLineStyle(2);
  lr_pt->SetLineColor(kGray + 2);
  lr_pt->Draw("same");

  // Normalize EEC by total jet weight (per-jet average)
  if (sumw_base > 0.0) h_eec_base->Scale(1.0 / sumw_base);
  if (sumw_total > 0.0) h_eec_total->Scale(1.0 / sumw_total);
  if (sumw_target > 0.0) h_eec_target->Scale(1.0 / sumw_target);

  TFile out(out_name.c_str(), "RECREATE");
  h_base->Write();
  h_total->Write();
  h_target->Write();
  h_w_unbias->Write();
  h_eec_base->Write();
  h_eec_total->Write();
  h_eec_target->Write();
  c->Write();

  c->SaveAs("plot_unbias_pT_check.pdf");


  TCanvas *cw = new TCanvas("c_weights", "weight check", 700, 600);
  cw->SetTicks(1,1);
  cw->SetLogy(); 
  h_w_unbias->SetLineColor(kRed + 1);
  h_w_unbias->SetLineWidth(2);
  h_w_unbias->Draw("hist");
  TLegend *lw = new TLegend(0.55, 0.75, 0.88, 0.88);
  lw->AddEntry(h_w_unbias, "unbias weights", "l");
  lw->SetBorderSize(0);
  lw->SetFillStyle(0);
  lw->Draw();
  cw->SaveAs("plot_unbias_weights_weights.pdf");
  delete cw;

  // EEC comparison plot
  TCanvas *ce = new TCanvas("c_eec", "EEC comparison", 800, 800);
  TPad *p_eec_top = new TPad("p_eec_top", "", 0.0, 0.30, 1.0, 1.0);
  TPad *p_eec_bot = new TPad("p_eec_bot", "", 0.0, 0.0, 1.0, 0.30);
  p_eec_top->SetBottomMargin(0.0);
  p_eec_top->SetLeftMargin(0.12);
  p_eec_top->SetRightMargin(0.04);
  p_eec_bot->SetTopMargin(0.0);
  p_eec_bot->SetBottomMargin(0.30);
  p_eec_bot->SetLeftMargin(0.12);
  p_eec_bot->SetRightMargin(0.04);
  p_eec_top->Draw();
  p_eec_bot->Draw();

  p_eec_top->cd();
  gPad->SetLogy(); 
  gPad->SetLogx();
  gPad->SetTicks(1,1); 
  h_eec_target->Draw("hist");
  h_eec_base->Draw("hist same");
  h_eec_total->Draw("hist same");
  TLegend *lege = new TLegend(0.52, 0.3, 0.88, 0.5);
  lege->AddEntry(h_eec_base, "biased sample (base weight)", "l");
  lege->AddEntry(h_eec_total, "weighted sample (base x unbias)", "l");
  lege->AddEntry(h_eec_target, Form("unbiased target %s", target_tree.c_str()), "l");
  lege->SetBorderSize(0);
  lege->SetFillStyle(0);
  lege->Draw();

  p_eec_bot->cd();
  gPad->SetLogx();
  gPad->SetTicks(1,1);
  TH1D *h_ratio_base_eec = (TH1D*)h_eec_base->Clone("h_ratio_base_eec");
  TH1D *h_ratio_total_eec = (TH1D*)h_eec_total->Clone("h_ratio_total_eec");
  h_ratio_base_eec->Divide(h_eec_target);
  h_ratio_total_eec->Divide(h_eec_target);
  h_ratio_total_eec->SetTitle("; #theta; ratio to target");
  h_ratio_total_eec->SetLineColor(kRed + 1);
  h_ratio_base_eec->SetLineColor(kBlue + 1);
  h_ratio_total_eec->SetLineWidth(2);
  h_ratio_base_eec->SetLineWidth(2);
  h_ratio_total_eec->SetMinimum(0.5);
  h_ratio_total_eec->SetMaximum(1.5);
  h_ratio_total_eec->GetYaxis()->SetNdivisions(505);
  h_ratio_total_eec->GetYaxis()->SetTitleSize(0.10);
  h_ratio_total_eec->GetYaxis()->SetTitleOffset(0.5);
  h_ratio_total_eec->GetYaxis()->SetLabelSize(0.08);
  h_ratio_total_eec->GetXaxis()->SetTitleSize(0.10);
  h_ratio_total_eec->GetXaxis()->SetLabelSize(0.08);
  h_ratio_total_eec->Draw("hist");
  h_ratio_base_eec->Draw("hist same");
  TLine *lr_eec = new TLine(h_ratio_total_eec->GetXaxis()->GetXmin(), 1.0,
                            h_ratio_total_eec->GetXaxis()->GetXmax(), 1.0);
  lr_eec->SetLineStyle(2);
  lr_eec->SetLineColor(kGray + 2);
  lr_eec->Draw("same");

  ce->SaveAs("plot_eec_compare.pdf");
  ce->Write();

  out.Close();

  in.Close();
  tfin.Close();
  {

    // ---------------------------------------------------------------
    // Basis-vector closure plots, read STRAIGHT from the fit output:
    //   source jets : tweights.g_basis        (weighted by w_base, w_total)
    //   target jets : tbasis_target.g_basis    (weighted by weight)
    //   titles      : meta.basis_labels
    // No basis is re-derived here, so these plots always match the fit
    // exactly (EEC terms, multiple subjet radii, custom powers, ...).
    // ---------------------------------------------------------------
    TH1::AddDirectory(kFALSE);   // histos below are transient, not file-owned
    TFile inb(input.c_str(), "READ");
    TTree *tsrc = inb.IsZombie() ? nullptr : dynamic_cast<TTree*>(inb.Get("tweights"));
    TTree *ttar = inb.IsZombie() ? nullptr : dynamic_cast<TTree*>(inb.Get("tbasis_target"));
    TTree *tmet = inb.IsZombie() ? nullptr : dynamic_cast<TTree*>(inb.Get("meta"));

    const bool have_src = tsrc && tsrc->GetBranch("g_basis");
    const bool have_tar = ttar && ttar->GetBranch("g_basis");
    if (!have_src || !have_tar) {
      std::cout << "warn: basis-vector plot skipped -- fit output is missing "
                << (!have_src ? "tweights.g_basis" : "tbasis_target.g_basis")
                << ".  Re-run unbias_weights.cc (it stores per-jet basis values)."
                << std::endl;
    } else {
      // per-function titles (optional)
      std::vector<std::string> *labels = nullptr;
      if (tmet && tmet->GetBranch("basis_labels")) {
        tmet->SetBranchAddress("basis_labels", &labels);
        tmet->GetEntry(0);
      }

      // source (tweights) branches
      std::vector<double> *sg = nullptr;
      float  spt = 0.0f; double sw_base = 1.0, sw_total = 1.0;
      tsrc->SetBranchAddress("g_basis", &sg);
      tsrc->SetBranchAddress("pt", &spt);
      tsrc->SetBranchAddress("w_base", &sw_base);
      tsrc->SetBranchAddress("w_total", &sw_total);

      // target (tbasis_target) branches
      std::vector<double> *tg = nullptr;
      double tgw = 1.0;
      ttar->SetBranchAddress("g_basis", &tg);
      ttar->SetBranchAddress("weight", &tgw);

      // number of basis functions from the first source entry
      int nphys = 0;
      if (tsrc->GetEntries() > 0) { tsrc->GetEntry(0); nphys = sg ? (int)sg->size() : 0; }

      if (nphys == 0) {
        std::cout << "warn: basis-vector plot skipped -- g_basis is empty." << std::endl;
      } else {
        auto label_of = [&](int k) -> std::string {
          if (labels && k < (int)labels->size() && !(*labels)[k].empty()) return (*labels)[k];
          return std::string(Form("g_%d", k));
        };

        const Long64_t ns = tsrc->GetEntries();
        const Long64_t nt = ttar->GetEntries();

        // pass 1: per-function positive min/max (source total + target) for log ranges
        std::vector<double> gmin(nphys, std::numeric_limits<double>::infinity());
        std::vector<double> gmax(nphys, 0.0);
        auto scan = [&](std::vector<double> *g) {
          if (!g) return;
          for (int k = 0; k < nphys && k < (int)g->size(); ++k) {
            const double v = std::abs((*g)[k]);
            if (v > 0.0) { gmin[k] = std::min(gmin[k], v); gmax[k] = std::max(gmax[k], v); }
          }
        };
        for (Long64_t i = 0; i < ns; ++i) { tsrc->GetEntry(i); if (spt < pt_min || spt > pt_max) continue; scan(sg); }
        for (Long64_t i = 0; i < nt; ++i) { ttar->GetEntry(i); scan(tg); }

        // build histograms with data-driven log binning (handles tiny EEC and huge pT^n alike)
        const int nbins_basis = 20;
        std::vector<TH1D*> hbase(nphys), htot(nphys), htar(nphys);
        for (int k = 0; k < nphys; ++k) {
          double lo = std::isfinite(gmin[k]) ? gmin[k] : 1e-3;
          double hi = (gmax[k] > lo) ? gmax[k] : lo * 10.0;
          double logmin = std::log10(lo) - 0.10;
          double logmax = std::log10(hi) + 0.10;
          if (!(logmax > logmin)) { logmin = -3.0; logmax = 3.0; }
          std::vector<double> edges(nbins_basis + 1);
          for (int i = 0; i <= nbins_basis; ++i)
            edges[i] = std::pow(10.0, logmin + (logmax - logmin) * i / nbins_basis);
          const std::string ttl = label_of(k) + "; " + label_of(k) + "; entries";
          hbase[k] = new TH1D(Form("hg_base_%d",   k), ttl.c_str(), nbins_basis, edges.data());
          htot[k]  = new TH1D(Form("hg_total_%d",  k), ttl.c_str(), nbins_basis, edges.data());
          htar[k]  = new TH1D(Form("hg_target_%d", k), ttl.c_str(), nbins_basis, edges.data());
          for (TH1D *h : { hbase[k], htot[k], htar[k] }) { h->Sumw2(); h->SetLineWidth(1); h->SetMarkerSize(0.5); }
          hbase[k]->SetLineColor(kBlue + 1);  hbase[k]->SetMarkerColor(kBlue + 1);  hbase[k]->SetMarkerStyle(20);
          htot[k] ->SetLineColor(kRed + 1);   htot[k] ->SetMarkerColor(kRed + 1);   htot[k] ->SetMarkerStyle(21);
          htar[k] ->SetLineColor(kGreen + 2); htar[k] ->SetMarkerColor(kGreen + 2); htar[k] ->SetMarkerStyle(22);
        }

        // pass 2: fill
        for (Long64_t i = 0; i < ns; ++i) {
          tsrc->GetEntry(i);
          if (spt < pt_min || spt > pt_max) continue;
          if (!sg) continue;
          for (int k = 0; k < nphys && k < (int)sg->size(); ++k) {
            const double v = std::abs((*sg)[k]);
            hbase[k]->Fill(v, sw_base);
            htot[k] ->Fill(v, sw_total);
          }
        }
        for (Long64_t i = 0; i < nt; ++i) {
          ttar->GetEntry(i);
          if (!tg) continue;
          for (int k = 0; k < nphys && k < (int)tg->size(); ++k)
            htar[k]->Fill(std::abs((*tg)[k]), tgw);
        }

        // canvas grid adapts to the number of basis functions
        const int ncol = std::max(1, (int)std::ceil(std::sqrt((double)nphys)));
        const int nrow = (nphys + ncol - 1) / ncol;
        TCanvas *ctw = new TCanvas("c_basis_weighted", "basis vectors (weighted)", 480 * ncol, 430 * nrow);
        ctw->Divide(ncol, nrow);
        TLegend *legw = new TLegend(0.12, 0.72, 0.88, 0.88);
        legw->AddEntry(htar[0],  Form("unbiased target %s", target_tree.c_str()), "lp");
        legw->AddEntry(hbase[0], "biased sample (base weight)", "lp");
        legw->AddEntry(htot[0],  "weighted sample (base x unbias)", "lp");
        legw->SetBorderSize(0); legw->SetFillStyle(0);

        for (int k = 0; k < nphys; ++k) {
          ctw->cd(k + 1);
          TPad *p_top = new TPad(Form("p_bt_%d", k), "", 0.0, 0.30, 1.0, 1.0);
          TPad *p_bot = new TPad(Form("p_bb_%d", k), "", 0.0, 0.0, 1.0, 0.30);
          p_top->SetBottomMargin(0.0); p_top->SetLeftMargin(0.14); p_top->SetRightMargin(0.04);
          p_bot->SetTopMargin(0.0); p_bot->SetBottomMargin(0.30); p_bot->SetLeftMargin(0.14); p_bot->SetRightMargin(0.04);
          p_top->Draw(); p_bot->Draw();

          p_top->cd(); gPad->SetLogy(); gPad->SetLogx(); gPad->SetTicks(1, 1);
          htar[k]->Draw("E1"); hbase[k]->Draw("E1 same"); htot[k]->Draw("E1 same");
          if (k == 0) legw->Draw();

          p_bot->cd(); gPad->SetLogx(); gPad->SetTicks(1, 1);
          TH1D *rb = (TH1D*)hbase[k]->Clone(Form("hg_ratio_base_%d",  k));
          TH1D *rt = (TH1D*)htot[k] ->Clone(Form("hg_ratio_total_%d", k));
          rb->Divide(htar[k]); rt->Divide(htar[k]);
          rt->SetTitle((std::string("; ") + label_of(k) + "; ratio to target").c_str());
          rt->SetLineColor(kRed + 1); rb->SetLineColor(kBlue + 1);
          rt->SetMinimum(0.5); rt->SetMaximum(1.5);
          rt->GetYaxis()->SetNdivisions(505);
          rt->GetYaxis()->SetTitleSize(0.10); rt->GetYaxis()->SetTitleOffset(0.5); rt->GetYaxis()->SetLabelSize(0.08);
          rt->GetXaxis()->SetTitleSize(0.10); rt->GetXaxis()->SetLabelSize(0.08);
          rt->Draw("E1"); rb->Draw("E1 same");
          TLine *lr = new TLine(rt->GetXaxis()->GetXmin(), 1.0, rt->GetXaxis()->GetXmax(), 1.0);
          lr->SetLineStyle(2); lr->SetLineColor(kGray + 2); lr->Draw("same");
        }
        ctw->SaveAs("plot_theta_basis_weighted_compare.pdf");
      }
    }
    inb.Close();
  }

  std::cout << "wrote " << out_name << ", plot_unbias_pT_check.pdf, "
            << "plot_unbias_weights_weights.pdf, plot_eec_compare.pdf and "
            << "plot_theta_basis_weighted_compare.pdf" << std::endl;
  return 0;
}
