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
#include <iostream>
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
// BASIS FUNCTION DEFINITIONS
// =====================================================================

struct BasisFuncDef {
    double A;               // ΔR^{−A}
    double B;               // [ln(ΔR)]^B
    int    m;               // (p_T^i p_T^k / p_T²)^m
    double E;               // angular upper cut: Θ(ΔR < E)
    std::string label;
};

static std::vector<BasisFuncDef> get_default_basis() {
    return {
        {1.0, 4.0, 1, 0.4, "g1: dR^{-3/2} ln^4 z^1  [dR<0.1]"},   // eq 4.15
        {1.0, 3.0, 1, 0.4, "g2: dR^{-1} z^2  [dR<0.2]"},           // eq 4.16
        {1.0, 2.0, 1, 0.4, "g3: dR^{-3/2} z^2  [dR<0.2]"},         // eq 4.17
        {1.0, 1.0, 1, 0.4, "g4: dR^{-1} ln^4 z^2  [dR<0.2]"},      // eq 4.18
        {1.0, 0.0, 1, 0.4, "g5: dR^{-3/2} ln^4 z^2  [dR<0.2]"},    // eq 4.19
        
        {0.0, 4.0, 1, 0.4, "g6: "},   // eq 4.15
        {0.0, 3.0, 1, 0.4, "g7: "},           // eq 4.16
        {0.0, 2.0, 1, 0.4, "g8:"},         // eq 4.17
        {0.0, 1.0, 1, 0.4, "g9: ]"},      // eq 4.18
        {0.0, 0.0, 1, 0.4, "g10"},    // eq 4.19
        
        {-1.0, 4.0, 1, 0.4, "g11: "},   // eq 4.15
        {-1.0, 3.0, 1, 0.4, "g12: "},           // eq 4.16
        {-1.0, 2.0, 1, 0.4, "g13:"},         // eq 4.17
        {-1.0, 1.0, 1, 0.4, "g14: ]"},      // eq 4.18
        {-1.0, 0.0, 1, 0.4, "g15"},    // eq 4.19
    };
}

// =====================================================================
// BASIS FUNCTION EVALUATION
// =====================================================================

/// Evaluate all basis functions for one jet.
/// dR_min screens the small-angle divergence in ΔR^{−A}.
static std::vector<double> evaluate_basis(
        const std::vector<BasisFuncDef> &basis,
        double ptjet,
        const std::vector<double> &cpt,
        const std::vector<double> &ceta,
        const std::vector<double> &cphi,
        double dR_min)
{
    const int nb = (int)basis.size();
    std::vector<double> g(nb, 0.0);
    if (ptjet <= 0.0) return g;

    const double pt2 = ptjet * ptjet;
    const size_t nc  = cpt.size();

    for (size_t i = 0; i < nc; ++i) {
        if (cpt[i] <= 0.0) continue;
        for (size_t k = i + 1; k < nc; ++k) {
            if (cpt[k] <= 0.0) continue;

            const double dphi = TVector2::Phi_mpi_pi(cphi[i] - cphi[k]);
            const double deta = ceta[i] - ceta[k];
            const double dR   = std::sqrt(deta * deta + dphi * dphi);
            if (dR <= dR_min) continue;   // screen small-angle divergence

            const double lndR = std::log(dR);
            const double z    = (cpt[i] * cpt[k]) / pt2;

            for (int j = 0; j < nb; ++j) {
                const auto &bf = basis[j];
                if (dR >= bf.E) continue;
                double val = std::pow(dR, -bf.A);
                if (bf.B > 0.0) val *= std::pow(lndR, bf.B);
                val *= std::pow(z, bf.m);
                g[j] += val;
            }
        }
    }
    return g;
}

// static std::vector<double> compute_theta_basis(double ptjet,
//                                                const std::vector<double> &cpt,
//                                                const std::vector<double> &ceta,
//                                                const std::vector<double> &cphi) {
//   std::vector<double> g(5, 0.0);
//   if (ptjet <= 0.0) return g;
//   const size_t nconst = cpt.size();
//   for (size_t i = 0; i < nconst; ++i) {
//     if (cpt[i] <= 0.0) continue;
//     for (size_t k = i + 1; k < nconst; ++k) {
//       if (cpt[k] <= 0.0) continue;
//       double dphi = TVector2::Phi_mpi_pi(cphi[i] - cphi[k]);
//       double deta = ceta[i] - ceta[k];
//       double theta = std::sqrt(deta * deta + dphi * dphi);
//       if (theta <= 0.0) continue;
//       double lnth = std::log(theta);
//       double ln4 = std::pow(lnth, 4);
//       double ratio = (cpt[i] * cpt[k]) / (ptjet * ptjet);
//       double ratio2 = ratio * ratio;
//       double t_m1 = std::pow(theta, -1.0);
//       double t_m32 = std::pow(theta, -1.5);
//       if (theta < 0.1) {
//         g[0] += t_m32 * ln4 * ratio;
//       }
//       g[1] += t_m1 * ratio2;
//       if (theta < 0.2) {
//         g[2] += t_m32 * ratio2;
//         g[3] += t_m1 * ln4 * ratio2;
//         g[4] += t_m32 * ln4 * ratio2;
//       }
//     }
//   }
//   return g;
// }

static void fill_eec(TH1D *h, const std::vector<double> &cpt,
                     const std::vector<double> &ceta,
                     const std::vector<double> &cphi,
                     double wjet) {
  const size_t nconst = cpt.size();
  if (nconst < 2) return;
  double sumpt = 0.0;
  for (size_t i = 0; i < nconst; ++i) {
    if (cpt[i] > 0.0) sumpt += cpt[i];
  }
  if (sumpt <= 0.0) return;
  const double norm = 1.0 / (sumpt * sumpt);
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
  const bool has_w_analytic = (t->GetBranch("w_analytic") != nullptr);
  if (has_w_analytic) {
    t->SetBranchAddress("w_analytic", &w_analytic);
  }
  if (t->GetBranch("source")) {
    t->SetBranchAddress("source", &source);
  }
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
    if (has_w_analytic) h_w_analytic->Fill(w_analytic);
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
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");

  TCanvas *c = new TCanvas("c_unbias", "unbias check", 900, 700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.04);
  h_base->Draw("hist");
  h_total->Draw("hist same");
  h_target->Draw("hist same");

  TLegend *leg = new TLegend(0.52, 0.70, 0.88, 0.88);
  leg->AddEntry(h_base, "biased sample (base weight)", "l");
  leg->AddEntry(h_total, "weighted sample (base x unbias)", "l");
  leg->AddEntry(h_target, Form("unbiased target %s", target_tree.c_str()), "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  // Normalize EEC by total jet weight (per-jet average)
  if (sumw_base > 0.0) h_eec_base->Scale(1.0 / sumw_base);
  if (sumw_total > 0.0) h_eec_total->Scale(1.0 / sumw_total);
  if (sumw_target > 0.0) h_eec_target->Scale(1.0 / sumw_target);

  TFile out(out_name.c_str(), "RECREATE");
  h_base->Write();
  h_total->Write();
  h_target->Write();
  h_w_unbias->Write();
  if (has_w_analytic) h_w_analytic->Write();
  h_eec_base->Write();
  h_eec_total->Write();
  h_eec_target->Write();
  c->Write();

  c->SaveAs("plot_unbias_weights_check.pdf");

  if (has_w_analytic) {
    TCanvas *cw = new TCanvas("c_weights", "weight check", 700, 600);
    h_w_unbias->SetLineColor(kRed + 1);
    h_w_unbias->SetLineWidth(2);
    h_w_analytic->SetLineColor(kBlack);
    h_w_analytic->SetLineWidth(2);
    h_w_analytic->SetLineStyle(2);
    h_w_analytic->Draw("hist");
    h_w_unbias->Draw("hist same");
    TLegend *lw = new TLegend(0.55, 0.75, 0.88, 0.88);
    lw->AddEntry(h_w_analytic, "analytic weights", "l");
    lw->AddEntry(h_w_unbias, "unbias weights", "l");
    lw->SetBorderSize(0);
    lw->SetFillStyle(0);
    lw->Draw();
    cw->SaveAs("plot_unbias_weights_weights.pdf");
    delete cw;
  }

  // EEC comparison plot
  TCanvas *ce = new TCanvas("c_eec", "EEC comparison", 800, 600);
  ce->SetLogy(); 
  ce->SetLogx(); 
  h_eec_target->Draw("hist");
  h_eec_base->Draw("hist same");
  h_eec_total->Draw("hist same");
  TLegend *lege = new TLegend(0.52, 0.70, 0.88, 0.88);
  lege->AddEntry(h_eec_base, "biased sample (base weight)", "l");
  lege->AddEntry(h_eec_total, "weighted sample (base x unbias)", "l");
  lege->AddEntry(h_eec_target, Form("unbiased target %s", target_tree.c_str()), "l");
  lege->SetBorderSize(0);
  lege->SetFillStyle(0);
  lege->Draw();
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

      std::vector<double> min = {6.396393e04*0.001, 5.151869e-01*0.001,5.468804*0.001,2.700238e02*0.001, 4.075668e03*0.001 };
      std::vector<double> max = {6.396393e04*1000, 5.151869e-01*1000,5.468804*1000,2.700238e02*1000, 4.075668e03*1000 };
      int nbins = 100; // Desired number of log bins

       
      TH1D *hgw_base[nphys];
      TH1D *hgw_total[nphys];
      TH1D *hgw_target[nphys];
      for (int k = 0; k < nphys; ++k) { 
        double logmin = TMath::Log10(10e-5);
        double logmax = TMath::Log10(10e5);
        double binwidth = (logmax - logmin) / nbins;
        double *edges = new double[nbins + 1];

        for (int i = 0; i <= nbins; i++) {
            edges[i] = pow(10, logmin + i * binwidth);
        }
       
        hgw_base[k] = new TH1D(Form("hThetaW_base_%d", k), Form("theta basis %d; g_%d; entries", k+1, k+1),  nbins, edges);
        hgw_total[k] = new TH1D(Form("hThetaW_total_%d", k), Form("theta basis %d; g_%d; entries", k+1, k+1),  nbins, edges);
        hgw_target[k] = new TH1D(Form("hThetaW_target_%d", k), Form("theta basis %d; g_%d; entries", k+1, k+1),  nbins, edges);             
        hgw_base[k]->SetLineColor(kBlue + 1);
        hgw_total[k]->SetLineColor(kRed + 1);
        hgw_target[k]->SetLineColor(kGreen + 2);        hgw_base[k]->SetLineWidth(2);
        hgw_total[k]->SetLineWidth(3);
        hgw_target[k]->SetLineWidth(3);      
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
          auto g = evaluate_basis(basis, wpt, *wconst_pt, *wconst_eta, *wconst_phi, 0.001);
          for (int k = 0; k < nphys; ++k) {
            hgw_base[k]->Fill(abs(g[k]), ww_base);
            hgw_total[k]->Fill(abs(g[k]), ww_total);
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
          auto g = evaluate_basis(basis, tpt2, *tconst_pt2, *tconst_eta2, *tconst_phi2, 0.001);

          for (int k = 0; k < nphys; ++k) {
            if(i == 0 ){
              std::cout << "g[" << k << "]: " << g[k] << " weight: " << tw2 <<  std::endl;
            }
            hgw_target[k]->Fill(abs(g[k]), tw2);
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
      TCanvas *ctw = new TCanvas("c_theta_basis_weighted", "theta basis weighted", 1000, 800);
      ctw->Divide(5, 3);
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
        hgw_target[k]->Draw("hist");
        hgw_base[k]->Draw("hist same");
        hgw_total[k]->Draw("hist same");
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
