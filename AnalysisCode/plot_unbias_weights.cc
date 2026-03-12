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

#include <cstdlib>
#include <iostream>
#include <string>

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

int main(int argc, char* argv[]) {
  const std::string input = get_arg(argc, argv, "--input");
  const std::string out_name = get_arg(argc, argv, "--out", "unbias_weights_plots.root");
  const std::string target_input = get_arg(argc, argv, "--target-input");
  const std::string target_tree = get_arg(argc, argv, "--target-tree", "tX");
  const double pt_min = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
  const double pt_max = std::stod(get_arg(argc, argv, "--pt-max", "200.0"));
  const int nbins = std::stoi(get_arg(argc, argv, "--nbins", "80"));

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

  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("w_unbias", &w_unbias);
  t->SetBranchAddress("w_base", &w_base);
  t->SetBranchAddress("w_total", &w_total);

  TH1D *h_raw = new TH1D("h_raw", "pT; p_{T} [GeV]; entries", nbins, pt_min, pt_max);
  TH1D *h_base = new TH1D("h_base", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_unbias = new TH1D("h_unbias", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_total = new TH1D("h_total", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_target = new TH1D("h_target", "pT; p_{T} [GeV]; entries", nbins, pt_min, pt_max);

  const Long64_t nentries = t->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    h_raw->Fill(pt, 1.0);
    h_base->Fill(pt, w_base);
    h_unbias->Fill(pt, w_unbias);
    h_total->Fill(pt, w_total);
  }

  // Load unbiased target tree (tX) for comparison
  TFile tfin(target_input.c_str(), "READ");
  if (tfin.IsZombie()) die("failed to open target input file: " + target_input);
  TTree *tt = dynamic_cast<TTree*>(tfin.Get(target_tree.c_str()));
  if (!tt) die("missing target tree: " + target_tree);

  float tpt = 0.0f;
  float tw = 1.0f;
  tt->SetBranchAddress("pt", &tpt);
  tt->SetBranchAddress("weight", &tw);
  const Long64_t tentries = tt->GetEntries();
  for (Long64_t i = 0; i < tentries; ++i) {
    tt->GetEntry(i);
    h_target->Fill(tpt, tw);
  }

  h_raw->SetLineColor(kBlack);
  h_base->SetLineColor(kBlue + 1);
  h_unbias->SetLineColor(kOrange + 1);
  h_total->SetLineColor(kRed + 1);
  h_target->SetLineColor(kGreen + 2);

  h_raw->SetLineWidth(2);
  h_base->SetLineWidth(2);
  h_unbias->SetLineWidth(2);
  h_total->SetLineWidth(2);
  h_target->SetLineWidth(2);

  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c_unbias", "unbias check", 900, 800);
  c->Divide(1, 2);
  c->cd(1);
  gPad->SetPad(0.0, 0.30, 1.0, 1.0);
  gPad->SetLogy();
  gPad->SetBottomMargin(0.02);
  h_raw->Draw("hist");
  h_base->Draw("hist same");
  h_unbias->Draw("hist same");
  h_total->Draw("hist same");
  h_target->Draw("hist same");

  TLegend *leg = new TLegend(0.52, 0.68, 0.88, 0.88);
  leg->AddEntry(h_raw, "raw (w=1)", "l");
  leg->AddEntry(h_base, "base weight", "l");
  leg->AddEntry(h_unbias, "unbias weight", "l");
  leg->AddEntry(h_total, "total weight", "l");
  leg->AddEntry(h_target, Form("target %s", target_tree.c_str()), "l");
  leg->Draw();

  // Ratio panel: total / target
  c->cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.30);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.30);
  TH1D *h_ratio = (TH1D*)h_total->Clone("h_ratio");
  h_ratio->SetTitle("; p_{T} [GeV]; total / target");
  h_ratio->Divide(h_target);
  h_ratio->SetLineColor(kRed + 1);
  h_ratio->SetLineWidth(2);
  h_ratio->SetMinimum(0.5);
  h_ratio->SetMaximum(1.5);
  h_ratio->Draw("hist");

  TFile out(out_name.c_str(), "RECREATE");
  h_raw->Write();
  h_base->Write();
  h_unbias->Write();
  h_total->Write();
  h_target->Write();
  h_ratio->Write();
  c->Write();
  out.Close();

  c->SaveAs("unbias_weights_check.pdf");

  in.Close();
  tfin.Close();
  std::cout << "wrote " << out_name << " and unbias_weights_check.pdf" << std::endl;
  return 0;
}
