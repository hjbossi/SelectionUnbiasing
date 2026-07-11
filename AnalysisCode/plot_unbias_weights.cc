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

#include "subjet_basis.h"

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
// BASIS FUNCTION DEFINITIONS + EVALUATION
// =====================================================================
// The basis (Cambridge/Aachen R=0.1 subjet pT power sums, g_n = Σ pT^n
// for n = 3..10) and its self-contained reclustering live in
// subjet_basis.h, shared with unbias_weights.cc so the fit and these
// validation plots use an identical definition.

static void set_logx_range_from_content(TH1D *h1, TH1D *h2, TH1D *h3) {
  double xmin = std::numeric_limits<double>::infinity();
  double xmax = 0.0;
  auto scan = [&](TH1D *h) {
    const int n = h->GetNbinsX();
    for (int i = 1; i <= n; ++i) {
      const double c = h->GetBinContent(i);
      if (c <= 0.0) continue;
      const double lo = h->GetBinLowEdge(i);
      const double hi = lo + h->GetBinWidth(i);
      if (lo > 0.0) xmin = std::min(xmin, lo);
      if (hi > 0.0) xmax = std::max(xmax, hi);
    }
  };
  scan(h1);
  scan(h2);
  scan(h3);
  if (!std::isfinite(xmin) || xmax <= xmin) return;
  const double logmin = std::log10(xmin);
  const double logmax = std::log10(xmax);
  const double span = std::max(1e-6, logmax - logmin);
  const double newmin = std::pow(10.0, logmin - 0.05 * span);
  const double newmax = std::pow(10.0, logmax + 0.05 * span);
  h1->GetXaxis()->SetRangeUser(newmin, newmax);
  h2->GetXaxis()->SetRangeUser(newmin, newmax);
  h3->GetXaxis()->SetRangeUser(newmin, newmax);
}

static void fill_eec(TH1D *h, const std::vector<double> &cpt,
                     const std::vector<double> &ceta,
                     const std::vector<double> &cphi,
                     double wjet) {
  const size_t nconst = cpt.size();
  if (nconst < 2) return;
  // double sumpt = 0.0;
  // for (size_t i = 0; i < nconst; ++i) {
  //   if (cpt[i] > 0.0) sumpt += cpt[i];
  // }
  // if (sumpt <= 0.0) return;
  // change the norm to be 120
  const double norm = 1.0 / (120.0 * 120.0);
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
  double w_analytic = 1.0;
  int source = -1;
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
  TH1D *h_w_analytic = new TH1D("h_w_analytic", "analytic weight; w; entries", 80, 0.0, 2.0);
  double eec_logmin = TMath::Log10(0.004);
  double eec_logmax = TMath::Log10(eec_max);
  double eec_binwidth = (eec_logmax - eec_logmin) /eec_bins;
  double *eec_edges = new double[eec_bins + 1];

  for (int i = 0; i <= eec_bins; i++) {
      eec_edges[i] = pow(10, eec_logmin + i * eec_binwidth);
  }
       
  TH1D *h_eec_base = new TH1D("h_eec_base", "EEC; #theta; EEC", eec_bins,eec_edges );
  TH1D *h_eec_total = new TH1D("h_eec_total", "EEC; #theta; EEC",  eec_bins,eec_edges);
  TH1D *h_eec_target = new TH1D("h_eec_target", "EEC; #theta; EEC",  eec_bins,eec_edges);

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
    // Weighted theta-basis distributions from tweights vs target tree (if constituents available)
    {
      
      std::vector<BasisFuncDef> basis = get_default_basis();
      const int nphys = (int)basis.size();

      int nbins = 20; // Desired number of log bins

       
      TH1D *hgw_base[nphys];
      TH1D *hgw_total[nphys];
      TH1D *hgw_target[nphys];
      std::vector<double> gmin_base(nphys, std::numeric_limits<double>::infinity());
      std::vector<double> gmax_base(nphys, 0.0);
      std::vector<double> gmin_total(nphys, std::numeric_limits<double>::infinity());
      std::vector<double> gmax_total(nphys, 0.0);
      std::vector<double> gmin_target(nphys, std::numeric_limits<double>::infinity());
      std::vector<double> gmax_target(nphys, 0.0);
      for (int k = 0; k < nphys; ++k) {
        // g_n = Σ pT^n grows steeply with n; size the log range so the
        // edges always cover the data (display auto-ranges below).
        const int nexp = basis[k].n;
        double logmin = -2.0;
        double logmax = nexp * TMath::Log10(300.0) + 2.0;
        double binwidth = (logmax - logmin) / nbins;
        double *edges = new double[nbins + 1];

        for (int i = 0; i <= nbins; i++) {
            edges[i] = pow(10, logmin + i * binwidth);
        }

        hgw_base[k] = new TH1D(Form("hThetaW_base_%d", k), Form("basis g_%d; g_%d; entries", nexp, nexp),  nbins, edges);
        hgw_total[k] = new TH1D(Form("hThetaW_total_%d", k), Form("basis g_%d; g_%d; entries", nexp, nexp),  nbins, edges);
        hgw_target[k] = new TH1D(Form("hThetaW_target_%d", k), Form("basis g_%d; g_%d; entries", nexp, nexp),  nbins, edges);
        hgw_base[k]->Sumw2();
        hgw_total[k]->Sumw2();
        hgw_target[k]->Sumw2();
        hgw_base[k]->SetLineColor(kBlue + 1);
        hgw_total[k]->SetLineColor(kRed + 1);
        hgw_target[k]->SetLineColor(kGreen + 2);
        hgw_base[k]->SetLineWidth(1);
        hgw_total[k]->SetLineWidth(1);
        hgw_target[k]->SetLineWidth(1);
        hgw_base[k]->SetMarkerStyle(20);
        hgw_total[k]->SetMarkerStyle(21);
        hgw_target[k]->SetMarkerStyle(22);
        hgw_base[k]->SetMarkerSize(0.5);
        hgw_total[k]->SetMarkerSize(0.5);
        hgw_target[k]->SetMarkerSize(0.5);
        hgw_base[k]->SetMarkerColor(kBlue + 1);
        hgw_total[k]->SetMarkerColor(kRed + 1);
        hgw_target[k]->SetMarkerColor(kGreen + 2);
      }

      // Fill weighted distributions from tweights (use a fresh file/tree to avoid branch conflicts)
      TFile in2(input.c_str(), "READ");
      TTree *tw2 = in2.IsZombie() ? nullptr : dynamic_cast<TTree*>(in2.Get("tweights"));
      if (!tw2 || !tw2->GetBranch("const_pt") || !tw2->GetBranch("const_eta") || !tw2->GetBranch("const_phi")) {
        std::cout << "warn: missing const_* branches in tweights for weighted theta compare" << std::endl;
      } else {
        float wpt = 0.0f;
        double ww_base = 1.0;
        double ww_total = 1.0;
        std::vector<double> *wconst_pt = nullptr;
        std::vector<double> *wconst_eta = nullptr;
        std::vector<double> *wconst_phi = nullptr;
        tw2->SetBranchAddress("pt", &wpt);
        tw2->SetBranchAddress("w_base", &ww_base);
        tw2->SetBranchAddress("w_total", &ww_total);
        tw2->SetBranchAddress("const_pt", &wconst_pt);
        tw2->SetBranchAddress("const_eta", &wconst_eta);
        tw2->SetBranchAddress("const_phi", &wconst_phi);
        const Long64_t wn = tw2->GetEntries();
        for (Long64_t i = 0; i < wn; ++i) {
          tw2->GetEntry(i);
          if (wpt < pt_min || wpt > pt_max) continue;
          if (!wconst_pt || !wconst_eta || !wconst_phi) continue;
          if (wconst_pt->size() != wconst_eta->size() || wconst_pt->size() != wconst_phi->size()) continue;
          // (basis, x, tcp->at(j), tce->at(j), tcf->at(j), dR_min);
          auto g = evaluate_basis(basis, wpt, *wconst_pt, *wconst_eta, *wconst_phi, kSubjetR);
          for (int k = 0; k < nphys; ++k) {
            const double gv = std::abs(g[k]);
            hgw_base[k]->Fill(gv, ww_base);
            hgw_total[k]->Fill(gv, ww_total);
            if (gv > 0.0) {
              gmin_base[k] = std::min(gmin_base[k], gv);
              gmax_base[k] = std::max(gmax_base[k], gv);
              gmin_total[k] = std::min(gmin_total[k], gv);
              gmax_total[k] = std::max(gmax_total[k], gv);
            }
          }
        }
      }
      in2.Close();

      // Fill target distributions (use a fresh file/tree to avoid branch conflicts)
      TFile tfin2(target_input.c_str(), "READ");
      TTree *tt2 = tfin2.IsZombie() ? nullptr : dynamic_cast<TTree*>(tfin2.Get(target_tree.c_str()));
      if (!tt2 || !tt2->GetBranch("const_pt") || !tt2->GetBranch("const_eta") || !tt2->GetBranch("const_phi")) {
        std::cout << "warn: missing const_* branches in target for weighted theta compare" << std::endl;
      } else {
        float tpt2 = 0.0f;
        float tw2 = 1.0f;
        std::vector<double> *tconst_pt2 = nullptr;
        std::vector<double> *tconst_eta2 = nullptr;
        std::vector<double> *tconst_phi2 = nullptr;
        tt2->SetBranchAddress("pt", &tpt2);
        tt2->SetBranchAddress("weight", &tw2);
        tt2->SetBranchAddress("const_pt", &tconst_pt2);
        tt2->SetBranchAddress("const_eta", &tconst_eta2);
        tt2->SetBranchAddress("const_phi", &tconst_phi2);
        const Long64_t tn2 = tt2->GetEntries();
        for (Long64_t i = 0; i < tn2; ++i) {
          tt2->GetEntry(i);
          if (tpt2 < pt_min || tpt2 > pt_max) continue;
          if (!tconst_pt2 || !tconst_eta2 || !tconst_phi2) continue;
          if (tconst_pt2->size() != tconst_eta2->size() || tconst_pt2->size() != tconst_phi2->size()) continue;
          //const std::vector<double> g = compute_theta_basis(tpt2, *tconst_pt2, *tconst_eta2, *tconst_phi2);
          auto g = evaluate_basis(basis, tpt2, *tconst_pt2, *tconst_eta2, *tconst_phi2, kSubjetR);

          for (int k = 0; k < nphys; ++k) {
            if(i == 0 ){
              std::cout << "g[" << k << "]: " << g[k] << " weight: " << tw2 <<  std::endl;
            }
            const double gv = std::abs(g[k]);
            hgw_target[k]->Fill(gv, tw2);
            if (gv > 0.0) {
              gmin_target[k] = std::min(gmin_target[k], gv);
              gmax_target[k] = std::max(gmax_target[k], gv);
            }
          }
        }
      }
      tfin2.Close();
      
      // for(int j = 0; j < 5 ; j++){
        
      // int size = hgw_base[j]->GetNbinsX(); 
      //   // 4. Fill the new histogram with old data
      //   for (int i = 1; i <=  size; i++) {
      //       hgw_base_rebin[j]->Fill(hgw_base[j]->GetBinCenter(i), hgw_base[j]->GetBinContent(i));
      //       hgw_total_rebin[j]->Fill(hgw_total[j]->GetBinCenter(i), hgw_total[j]->GetBinContent(i));
      //       hgw_target_rebin[j]->Fill(hgw_target[j]->GetBinCenter(i), hgw_target[j]->GetBinContent(i));
      //   }
    
      // }
      // Overlay plots with ratio panel
      TCanvas *ctw = new TCanvas("c_theta_basis_weighted", "theta basis weighted", 6500, 5000);
      ctw->Divide(4, 2);
      TLegend *legw = new TLegend(0.12, 0.72, 0.88, 0.88);
      legw->AddEntry(hgw_target[0], Form("unbiased target %s", target_tree.c_str()), "l");
      legw->AddEntry(hgw_base[0], "biased sample (base weight)", "l");
      legw->AddEntry(hgw_total[0], "weighted sample (base x unbias)", "l");
      legw->SetBorderSize(0);
      legw->SetFillStyle(0);
      for (int k = 0; k < nphys; ++k) {
        ctw->cd(k + 1);
        TPad *p_top = new TPad(Form("p_theta_top_%d", k), "", 0.0, 0.30, 1.0, 1.0);
        TPad *p_bot = new TPad(Form("p_theta_bot_%d", k), "", 0.0, 0.0, 1.0, 0.30);
        p_top->SetBottomMargin(0.0);
        p_top->SetLeftMargin(0.14);
        p_top->SetRightMargin(0.04);
        p_bot->SetTopMargin(0.0);
        p_bot->SetBottomMargin(0.30);
        p_bot->SetLeftMargin(0.14);
        p_bot->SetRightMargin(0.04);
        p_top->Draw();
        p_bot->Draw();

        p_top->cd();
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetTicks(1,1);
        std::cout << "On canvas k=" << k << " then drawing hist with integral " << hgw_target[k]->Integral() << std::endl;
        // Use per-basis min/max from the actual g values to set display range.
        double xmin = std::min({gmin_base[k], gmin_total[k], gmin_target[k]});
        double xmax = std::max({gmax_base[k], gmax_total[k], gmax_target[k]});
        if (std::isfinite(xmin) && xmax > xmin) {
          const double logmin = std::log10(xmin);
          const double logmax = std::log10(xmax);
          const double span = std::max(1e-6, logmax - logmin);
          const double newmin = std::pow(10.0, logmin - 0.05 * span);
          const double newmax = std::pow(10.0, logmax + 0.05 * span);
          hgw_target[k]->GetXaxis()->SetRangeUser(newmin, newmax);
          hgw_base[k]->GetXaxis()->SetRangeUser(newmin, newmax);
          hgw_total[k]->GetXaxis()->SetRangeUser(newmin, newmax);
        } else {
          set_logx_range_from_content(hgw_target[k], hgw_base[k], hgw_total[k]);
        }
        hgw_target[k]->Draw("E1");
        hgw_base[k]->Draw("E1 same");
        hgw_total[k]->Draw("E1 same");
        if (k == 0) legw->Draw();

        // Ratio to target (preserve binning)
        p_bot->cd();
        gPad->SetLogx();
        gPad->SetTicks(1,1); 
        TH1D *h_ratio_base = (TH1D*)hgw_base[k]->Clone(Form("hThetaW_ratio_base_%d", k));
        TH1D *h_ratio_total = (TH1D*)hgw_total[k]->Clone(Form("hThetaW_ratio_total_%d", k));
        h_ratio_base->Divide(hgw_target[k]);
        h_ratio_total->Divide(hgw_target[k]);
        h_ratio_total->SetTitle("; g_{k}; ratio to target");
        h_ratio_total->SetLineColor(kRed + 1);
        h_ratio_base->SetLineColor(kBlue + 1);
        h_ratio_total->SetLineWidth(1);
        h_ratio_base->SetLineWidth(1);
        h_ratio_total->SetMarkerStyle(21);
        h_ratio_base->SetMarkerStyle(20);
        h_ratio_total->SetMarkerSize(0.5);
        h_ratio_base->SetMarkerSize(0.5);
        h_ratio_total->SetMinimum(0.5);
        h_ratio_total->SetMaximum(1.5);
        h_ratio_total->GetYaxis()->SetNdivisions(505);
        h_ratio_total->GetYaxis()->SetTitleSize(0.10);
        h_ratio_total->GetYaxis()->SetTitleOffset(0.5);
        h_ratio_total->GetYaxis()->SetLabelSize(0.08);
        h_ratio_total->GetXaxis()->SetTitleSize(0.10);
        h_ratio_total->GetXaxis()->SetLabelSize(0.08);
        h_ratio_total->Draw("E1");
        h_ratio_base->Draw("E1 same");
        TLine *lr = new TLine(h_ratio_total->GetXaxis()->GetXmin(), 1.0,
                              h_ratio_total->GetXaxis()->GetXmax(), 1.0);
        lr->SetLineStyle(2);
        lr->SetLineColor(kGray + 2);
        lr->Draw("same");
      }
      ctw->SaveAs("plot_theta_basis_weighted_compare.pdf");
    }
  }

  std::cout << "wrote " << out_name << " and plot_unbias_weights_check.pdf" << std::endl;
  return 0;
}