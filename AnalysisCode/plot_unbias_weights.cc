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
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "subjet_basis.h"     // get_default_basis(), kSubjetR, recluster_ca_subjet_pts()
#include "basis_functions.h"  // generic configurable Basis framework (JetData, evaluate_basis, ...)

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
// The basis-vector plots below use exactly the basis that was fit.  This macro
// reconstructs the runtime basis from the fit metadata (meta.basis_labels +
// meta.basis) via rebuild_basis_from_meta(), then re-evaluates it from each
// jet's constituents (reclustering the stored subjet radii on the fly), and
// titles each panel with the corresponding meta.basis_labels entry.  The plots
// therefore follow the fitted basis (EEC terms, several subjet radii, custom
// powers, ...) rather than any fixed default basis.

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

// =====================================================================
// Basis reconstruction from fit metadata
// =====================================================================
// Defensive placeholder used when a basis_labels entry cannot be parsed, so
// the reconstructed basis stays index-aligned with lambda / basis_labels.
class ZeroBasisFunction : public BasisFunction {
public:
  explicit ZeroBasisFunction(std::string label) : label_(std::move(label)) {}
  double evaluate(const JetData &) const override { return 0.0; }
  std::string label() const override { return label_; }
private:
  std::string label_;
};

// Reconstruct the runtime basis that was actually fit, straight from the ROOT
// metadata written by unbias_weights.cc:
//   * meta.basis_labels (vector<string>) -- authoritative order + panel titles
//   * meta.basis        (string)         -- machine-readable signature (norm=..)
// Each label is parsed back into the matching BasisFunction, in order.  On an
// older file without a meta / basis_labels branch we warn and fall back to the
// fixed default basis (get_default_basis()).  labels_out receives the labels
// used (either the metadata labels or the default-basis labels).
static Basis rebuild_basis_from_meta(const std::string &input,
                                     std::vector<std::string> &labels_out) {
  Basis basis;
  labels_out.clear();

  auto fallback = [&]() {
    std::cout << "warn: meta / basis_labels not found in " << input
              << " -- falling back to get_default_basis()." << std::endl;
    basis.clear();
    labels_out.clear();
    for (const BasisFuncDef &d : get_default_basis()) {
      labels_out.push_back(d.label);
      basis.push_back(std::make_unique<SubjetMomentBasisFunction>(
          kSubjetR, static_cast<double>(d.n), d.label));
    }
  };

  TFile f(input.c_str(), "READ");
  TTree *meta = f.IsZombie() ? nullptr : dynamic_cast<TTree *>(f.Get("meta"));
  if (!meta || !meta->GetBranch("basis_labels")) {
    fallback();
    return basis;
  }

  std::vector<std::string> *labels = nullptr;
  std::string *sig = nullptr;
  meta->SetBranchAddress("basis_labels", &labels);
  if (meta->GetBranch("basis")) meta->SetBranchAddress("basis", &sig);
  meta->GetEntry(0);

  if (!labels) {  // branch present but could not be read back
    fallback();
    return basis;
  }

  // EEC normalisation from the machine-readable signature (norm=%lf); the
  // historical default is 120 when the signature carries no norm.
  double eec_norm = 120.0;
  if (sig) {
    const auto pos = sig->find("norm=");
    if (pos != std::string::npos) {
      double v = 0.0;
      if (std::sscanf(sig->c_str() + pos, "norm=%lf", &v) == 1) eec_norm = v;
    }
  }

  for (const std::string &lab : *labels) {
    labels_out.push_back(lab);
    double R = 0.0, n = 0.0, A = 0.0, B = 0.0, E = 0.0;
    int    m = 0;
    if (std::sscanf(lab.c_str(),
                    "sjmom: sum pT^{%lf} (stored subjets R=%lf)", &n, &R) == 2) {
      basis.push_back(std::make_unique<SubjetMomentBasisFunction>(R, n, lab));
    } else if (std::sscanf(lab.c_str(),
                           "eec: dR^{-%lf} ln^{%lf} z^{%d} [dR<%lf]",
                           &A, &B, &m, &E) == 4) {
      basis.push_back(std::make_unique<EECBasisFunction>(A, B, m, E, eec_norm, lab));
    } else {
      std::cout << "warn: unrecognised basis label \"" << lab
                << "\" -- inserting a zero-valued placeholder to keep index"
                   " alignment with lambda / basis_labels." << std::endl;
      basis.push_back(std::make_unique<ZeroBasisFunction>(lab));
    }
  }
  return basis;
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
    // ---------------------------------------------------------------
    // Weighted theta-basis distributions.
    //
    // The basis is reconstructed from the fit metadata (meta.basis /
    // meta.basis_labels) and RE-EVALUATED from each jet's constituents, so
    // the panels follow exactly the basis that was fit and are titled by the
    // matching meta.basis_labels entry.  Source jets come from tweights
    // (weighted by w_base, w_total); target jets from the unbiased target
    // tree (weighted by weight).  We do not depend on any stored g_basis /
    // subjet branches.
    // ---------------------------------------------------------------
    TH1::AddDirectory(kFALSE);   // histos below are transient, not file-owned

    // Reconstruct the fitted basis + its per-function labels from metadata.
    std::vector<std::string> basis_labels;
    Basis basis = rebuild_basis_from_meta(input, basis_labels);
    const int nphys = (int)basis.size();

    std::cout << "theta-basis: reconstructed " << nphys
              << " basis function(s) from metadata:" << std::endl;
    for (int k = 0; k < nphys; ++k)
      std::cout << "  [" << k << "] " << basis_labels[k] << std::endl;

    if (nphys == 0) {
      std::cout << "warn: basis-vector plot skipped -- reconstructed basis is empty." << std::endl;
    } else {
      const BasisInputs req = basis_inputs(basis);

      // Evaluate the generic basis for one jet from its constituents: recluster
      // the stored subjet radii on the fly, build the shared EEC pair table if
      // any EEC term is present, then evaluate every basis function.
      auto eval_jet = [&](double jpt,
                          const std::vector<double> &cpt,
                          const std::vector<double> &ceta,
                          const std::vector<double> &cphi) -> std::vector<double> {
        JetData jd;
        jd.pt = jpt;
        jd.const_pt = &cpt; jd.const_eta = &ceta; jd.const_phi = &cphi;
        std::map<int, std::vector<double>> subj;   // stable storage for the pointers below
        for (int rtag : req.subjet_rtags) {
          subj[rtag] = recluster_ca_subjet_pts(cpt, ceta, cphi, radius_from_tag(rtag));
          jd.subjet_pt[rtag] = &subj[rtag];
        }
        if (req.constituents) build_constituent_pairs(jd, 0.001);
        return evaluate_basis(basis, jd);
      };

      // pass 1 (single evaluation): store per-jet g vectors + weights and
      // accumulate per-function positive min/max of |g| for the log binning.
      std::vector<double> gmin(nphys, std::numeric_limits<double>::infinity());
      std::vector<double> gmax(nphys, 0.0);
      auto accumulate = [&](const std::vector<double> &g) {
        for (int k = 0; k < nphys && k < (int)g.size(); ++k) {
          const double v = std::abs(g[k]);
          if (v > 0.0) { gmin[k] = std::min(gmin[k], v); gmax[k] = std::max(gmax[k], v); }
        }
      };

      std::vector<std::vector<double>> src_g;      // per source jet g vector
      std::vector<double> src_wbase, src_wtotal;
      std::vector<std::vector<double>> tar_g;      // per target jet g vector
      std::vector<double> tar_w;

      // source jets: tweights (weighted by w_base, w_total)
      TFile insrc(input.c_str(), "READ");
      TTree *tsrc = insrc.IsZombie() ? nullptr : dynamic_cast<TTree*>(insrc.Get("tweights"));
      if (!tsrc || !tsrc->GetBranch("const_pt") || !tsrc->GetBranch("const_eta") ||
          !tsrc->GetBranch("const_phi")) {
        std::cout << "warn: missing const_* branches in tweights for weighted theta compare" << std::endl;
      } else {
        float  spt = 0.0f; double sw_base = 1.0, sw_total = 1.0;
        std::vector<double> *scp = nullptr, *sce = nullptr, *scf = nullptr;
        tsrc->SetBranchAddress("pt", &spt);
        tsrc->SetBranchAddress("w_base", &sw_base);
        tsrc->SetBranchAddress("w_total", &sw_total);
        tsrc->SetBranchAddress("const_pt", &scp);
        tsrc->SetBranchAddress("const_eta", &sce);
        tsrc->SetBranchAddress("const_phi", &scf);
        const Long64_t ns = tsrc->GetEntries();
        for (Long64_t i = 0; i < ns; ++i) {
          tsrc->GetEntry(i);
          if (spt < pt_min || spt > pt_max) continue;
          if (!scp || !sce || !scf) continue;
          if (scp->size() != sce->size() || scp->size() != scf->size()) continue;
          std::vector<double> g = eval_jet(spt, *scp, *sce, *scf);
          accumulate(g);
          src_g.push_back(std::move(g));
          src_wbase.push_back(sw_base);
          src_wtotal.push_back(sw_total);
        }
      }
      insrc.Close();

      // target jets: target_tree in target_input (weighted by weight)
      TFile intar(target_input.c_str(), "READ");
      TTree *ttar = intar.IsZombie() ? nullptr : dynamic_cast<TTree*>(intar.Get(target_tree.c_str()));
      if (!ttar || !ttar->GetBranch("const_pt") || !ttar->GetBranch("const_eta") ||
          !ttar->GetBranch("const_phi")) {
        std::cout << "warn: missing const_* branches in target for weighted theta compare" << std::endl;
      } else {
        float tpt = 0.0f, tw = 1.0f;
        std::vector<double> *tcp = nullptr, *tce = nullptr, *tcf = nullptr;
        ttar->SetBranchAddress("pt", &tpt);
        ttar->SetBranchAddress("weight", &tw);
        ttar->SetBranchAddress("const_pt", &tcp);
        ttar->SetBranchAddress("const_eta", &tce);
        ttar->SetBranchAddress("const_phi", &tcf);
        const Long64_t nt = ttar->GetEntries();
        for (Long64_t i = 0; i < nt; ++i) {
          ttar->GetEntry(i);
          if (tpt < pt_min || tpt > pt_max) continue;
          if (!tcp || !tce || !tcf) continue;
          if (tcp->size() != tce->size() || tcp->size() != tcf->size()) continue;
          std::vector<double> g = eval_jet(tpt, *tcp, *tce, *tcf);
          accumulate(g);
          tar_g.push_back(std::move(g));
          tar_w.push_back(tw);
        }
      }
      intar.Close();

      // build histograms with data-driven log binning, using each function's
      // own observed |g| range (handles tiny EEC and huge pT^n alike)
      const int nbins = 20;
      std::vector<TH1D*> hbase(nphys), htot(nphys), htar(nphys);
      for (int k = 0; k < nphys; ++k) {
        double lo = std::isfinite(gmin[k]) ? gmin[k] : 1e-3;
        double hi = (gmax[k] > lo) ? gmax[k] : lo * 10.0;
        double logmin = std::log10(lo) - 0.10;
        double logmax = std::log10(hi) + 0.10;
        if (!(logmax > logmin)) { logmin = -3.0; logmax = 3.0; }
        std::vector<double> edges(nbins + 1);
        for (int i = 0; i <= nbins; ++i)
          edges[i] = std::pow(10.0, logmin + (logmax - logmin) * i / nbins);
        const std::string lab = (k < (int)basis_labels.size())
                                    ? basis_labels[k] : std::string(Form("g_%d", k));
        const std::string ttl = lab + "; " + lab + "; entries";   // x-axis title = basis label
        hbase[k] = new TH1D(Form("hg_base_%d",   k), ttl.c_str(), nbins, edges.data());
        htot[k]  = new TH1D(Form("hg_total_%d",  k), ttl.c_str(), nbins, edges.data());
        htar[k]  = new TH1D(Form("hg_target_%d", k), ttl.c_str(), nbins, edges.data());
        for (TH1D *h : { hbase[k], htot[k], htar[k] }) { h->Sumw2(); h->SetLineWidth(1); h->SetMarkerSize(0.5); }
        hbase[k]->SetLineColor(kBlue + 1);  hbase[k]->SetMarkerColor(kBlue + 1);  hbase[k]->SetMarkerStyle(20);
        htot[k] ->SetLineColor(kRed + 1);   htot[k] ->SetMarkerColor(kRed + 1);   htot[k] ->SetMarkerStyle(21);
        htar[k] ->SetLineColor(kGreen + 2); htar[k] ->SetMarkerColor(kGreen + 2); htar[k] ->SetMarkerStyle(22);
      }

      // pass 2: fill from the stored g vectors (basis is NOT evaluated twice)
      for (size_t i = 0; i < src_g.size(); ++i)
        for (int k = 0; k < nphys && k < (int)src_g[i].size(); ++k) {
          const double v = std::abs(src_g[i][k]);
          hbase[k]->Fill(v, src_wbase[i]);
          htot[k] ->Fill(v, src_wtotal[i]);
        }
      for (size_t i = 0; i < tar_g.size(); ++i)
        for (int k = 0; k < nphys && k < (int)tar_g[i].size(); ++k)
          htar[k]->Fill(std::abs(tar_g[i][k]), tar_w[i]);

      // canvas grid adapts to the number of basis functions
      const int cols = std::max(1, (int)std::ceil(std::sqrt((double)nphys)));
      const int rows = (int)std::ceil((double)nphys / cols);
      TCanvas *ctw = new TCanvas("c_basis_weighted", "basis vectors (weighted)", 480 * cols, 430 * rows);
      ctw->Divide(cols, rows);
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
        const std::string rlab = (k < (int)basis_labels.size())
                                     ? basis_labels[k] : std::string(Form("g_%d", k));
        rt->SetTitle((std::string("; ") + rlab + "; ratio to target").c_str());
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

  std::cout << "wrote " << out_name << " and plot_unbias_weights_check.pdf" << std::endl;
  return 0;
}