// startBasis.C: file to plot the initial basis and see the effect of a pT shift on these distirbutions
// Hannah Bossi, <hannah.bossi@cern.ch>
// December 19th, 2025

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <memory>

using namespace std;

// ----------------------------------------------------
// get the files from a specified input directory
void GetFiles(char const *input, vector<string> &files) {
  TSystemDirectory dir(input, input);
  TList *list = dir.GetListOfFiles();

  if (list) {
    TSystemFile *file;
    string fname;
    TIter next(list);
    while ((file = (TSystemFile *)next())) {
      fname = file->GetName();

      if (file->IsDirectory() && (fname.find(".") == string::npos)) {
        string newDir = string(input) + fname + "/";
        GetFiles(newDir.c_str(), files);
      } else if ((fname.find(".root") != string::npos)) {
        files.push_back(string(input) + fname);
        cout << files.back() << endl;
      }
    }
  }
}

// ----------------------------------------------------
// Fill a TChain from the specified files
void FillChain(TChain &chain, vector<string> &files) {
  for (auto &file : files) {
    chain.Add(file.c_str());
  }
}

// ----------------------------------------------------

// main function
void startBasis(const char* inputDir = "/home/hbossi/SelectionUnbiasing/MCOutput/v3/",
                const char* outFile  = "startBasis_output.root",
                double pTLow = 100,
                double pTHigh = 140,
                double p = 8.0,
                double d = 0.1,
                unsigned rngSeed = 12345) {

  // config
  const char* treeName = "tgenBefore";

  std::vector<string> files;
  GetFiles(inputDir, files);

  TChain chain(treeName);
  FillChain(chain, files);

  cout << "Total entries in chain: "
       << chain.GetEntries() << endl;

  // -----------------------------
  // TTreeReader setup
  // -----------------------------
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t>    nJets(reader, "nJets");
  TTreeReaderArray<Float_t> pt(reader, "pt");
  TTreeReaderArray<Float_t> eta(reader, "eta");
  TTreeReaderArray<Float_t> phi(reader, "phi");
  TTreeReaderArray<Float_t> mass(reader, "mass");
  TTreeReaderValue<Float_t> weight(reader, "weight");
  TTreeReaderValue<std::vector<std::vector<double>>> const_pt(reader, "const_pt");
  TTreeReaderValue<std::vector<std::vector<double>>> const_eta(reader, "const_eta");
  TTreeReaderValue<std::vector<std::vector<double>>> const_phi(reader, "const_phi");

  // -----------------------------
  // Subjet radius scan: forward the precomputed C/A subjets from ppjets_root.cc
  // -----------------------------
  // ppjets_root.cc stores subjets over R = 0.01..0.20 in branches
  // "subjet_pt_R0pXX" (vector<vector<double>>, indexed [jet][subjet]).  For each
  // such branch present in the input we set up a reader and, below, a matching
  // per-jet "subjet_pt_R0pXX" output branch (vector<double>) on the flat trees
  // — mirroring how const_pt/eta/phi are propagated.  Radii absent from the
  // input (e.g. older ntuples) are skipped, so base behaviour is unchanged when
  // the branches don't exist.
  const int nRadiiScan = 20;
  chain.LoadTree(0);  // make the branch list available for the presence check
  std::vector<std::string> subjet_names;
  std::vector<std::unique_ptr<TTreeReaderValue<std::vector<std::vector<double>>>>> subjet_in;
  for (int iR = 0; iR < nRadiiScan; ++iR) {
    std::string bname = Form("subjet_pt_R0p%02d", iR + 1);
    if (!chain.GetBranch(bname.c_str())) continue;   // radius not in this input
    subjet_names.push_back(bname);
    subjet_in.emplace_back(
      std::make_unique<TTreeReaderValue<std::vector<std::vector<double>>>>(
        reader, subjet_names.back().c_str()));
  }
  const size_t nSubR = subjet_names.size();
  std::cout << "Forwarding " << nSubR << " subjet_pt_R0pXX branch(es)" << std::endl;


  // Histograms
  const int nBins = 100;
  const double ptMin = 90.0;
  const double ptMax = 150.0;
  TH1D* hPtX  = new TH1D("hPtX", "; #it{p}_{T} (GeV); d#it{N}/d#it{p}_{T}", nBins, ptMin, ptMax);
  TH1D* hPtY = new TH1D("hPtY","; #it{p}_{T} (GeV); d#it{N}/d#it{p}_{T}", nBins, ptMin, ptMax);

  hPtX->SetLineColor(kBlack);
  hPtY->SetLineColor(kRed);
  hPtX->SetLineWidth(2);
  hPtY->SetLineWidth(2);

  // -----------------------------
  // Output trees (per-jet)
  // -----------------------------
  TFile fout(outFile, "RECREATE");
  TTree *tX = new TTree("tX", "biased X jets (raw)");
  TTree *tY = new TTree("tY", "biased Y jets (raw, windowed)");
  TTree *tYp = new TTree("tYprime", "Y' sample (oversampled)");
  TTree *tRef = new TTree("tRef", "reference sample (X + Y)");
  TTree *tBiased = new TTree("tBiased", "biased sample (X + Y')");
  TTree *tPP = new TTree("tPP", "unbiased sample (X + Y)"); // backward compat

  float out_pt = 0.0f;
  float out_pt_raw = 0.0f;
  float out_pt_shifted = 0.0f;
  float out_eta = 0.0f;
  float out_phi = 0.0f;
  float out_mass = 0.0f;
  float out_weight = 1.0f;
  int out_source = 0; // 0 = X, 1 = Y, 2 = Y'
  std::vector<double> out_const_pt;
  std::vector<double> out_const_eta;
  std::vector<double> out_const_phi;
  // one per-jet subjet-pT vector per forwarded radius (indices align with
  // subjet_names).  Fixed size => &out_subjet_pt[r] stays valid for the whole
  // run, so the branch addresses below never dangle.
  std::vector<std::vector<double>> out_subjet_pt(nSubR);

  auto setup_tree = [&](TTree *t) {
    t->Branch("pt", &out_pt, "pt/F");
    t->Branch("pt_raw", &out_pt_raw, "pt_raw/F");
    t->Branch("pt_shifted", &out_pt_shifted, "pt_shifted/F");
    t->Branch("eta", &out_eta, "eta/F");
    t->Branch("phi", &out_phi, "phi/F");
    t->Branch("mass", &out_mass, "mass/F");
    t->Branch("weight", &out_weight, "weight/F");
    t->Branch("source", &out_source, "source/I");
    t->Branch("const_pt", &out_const_pt);
    t->Branch("const_eta", &out_const_eta);
    t->Branch("const_phi", &out_const_phi);
    for (size_t r = 0; r < nSubR; ++r)
      t->Branch(subjet_names[r].c_str(), &out_subjet_pt[r]);
  };
  setup_tree(tX);
  setup_tree(tY);
  setup_tree(tYp);
  setup_tree(tRef);
  setup_tree(tBiased);
  setup_tree(tPP);

  struct JetRec {
    float pt;
    float pt_raw;
    float pt_shifted;
    float eta;
    float phi;
    float mass;
    float weight;
    std::vector<double> const_pt;
    std::vector<double> const_eta;
    std::vector<double> const_phi;
    std::vector<std::vector<double>> subjet_pt;  // per forwarded radius
  };
  std::vector<JetRec> y_cache;
  std::mt19937 rng(rngSeed);
  std::uniform_real_distribution<double> uni01(0.0, 1.0);

  // -----------------------------
  // Event loop
  // -----------------------------
  Long64_t iev = 0;
  while (reader.Next()) {

    if (iev % 100000 == 0)
      cout << "On event " << iev << endl;
    
    ++iev;

    for (int j = 0; j < *nJets; ++j) {

      double ptRaw = pt[j];
      if (ptRaw <= 0) continue;
      // realistic eta cut
      if(eta[j] < -2.4 || eta[j] > 2.4) continue;

      // gather this jet's forwarded subjets (all present radii) once per jet.
      // Element-wise assignment keeps the outer vector (and thus the branch
      // addresses set up above) from reallocating.
      for (size_t r = 0; r < nSubR; ++r) {
        const std::vector<std::vector<double>>* sj = subjet_in[r]->Get();
        out_subjet_pt[r] = (sj && j < (int)sj->size())
                           ? (*sj)[j] : std::vector<double>();
      }

      // Every jet in the pT window enters both as an X jet (source 0) and as
      // a Y jet (source 1): X is downsampled into the biased sample by the
      // bias weight, while Y is cached and later oversampled into Y'.
      if ((pTLow < ptRaw) && (ptRaw < pTHigh)) {
        out_pt = ptRaw;
        out_pt_raw = ptRaw;
        out_pt_shifted = ptRaw;
        out_eta = eta[j];
        out_phi = phi[j];
        out_mass = mass[j];
        out_weight = *weight;
        out_const_pt = (*const_pt)[j];
        out_const_eta = (*const_eta)[j];
        out_const_phi = (*const_phi)[j];

        // ---- X jet ----
        hPtX->Fill(ptRaw, *weight);
        out_source = 0;
        tX->Fill();
        tRef->Fill();
        tPP->Fill();

        // Downsample X with probability = biasWeight
        const double ptL_p   = std::pow(pTLow, p);
        const double ptH_p   = std::pow(pTHigh, p);
        const double ptRaw_p = std::pow(ptRaw, p);
        const double biasWeight =
            (ptL_p - ptH_p * d + (-1.0 + d) * ptRaw_p) / (ptL_p - ptH_p);
        if (uni01(rng) < biasWeight) tBiased->Fill();

        // ---- Y jet ----
        hPtY->Fill(ptRaw, *weight);
        out_source = 1;
        tY->Fill();
        tRef->Fill();
        tPP->Fill();

        JetRec rec;
        rec.pt = out_pt;
        rec.pt_raw = out_pt_raw;
        rec.pt_shifted = out_pt_shifted;
        rec.eta = out_eta;
        rec.phi = out_phi;
        rec.mass = out_mass;
        rec.weight = out_weight;
        rec.const_pt = out_const_pt;
        rec.const_eta = out_const_eta;
        rec.const_phi = out_const_phi;
        rec.subjet_pt = out_subjet_pt;   // copy all forwarded radii
        y_cache.push_back(std::move(rec));
      }
    }
  }

  // Build Y' by oversampling Y (bias) while preserving total biased sample size.
  if (!y_cache.empty()) {
    const size_t x_count = static_cast<size_t>(tX->GetEntries());
    const size_t y_count = y_cache.size();
    const size_t x_kept = static_cast<size_t>(tBiased->GetEntries()); // after downsampling X
    const size_t desired_total = x_count + y_count;
    size_t yprime_count = 0;
    if (desired_total > x_kept) {
      yprime_count = desired_total - x_kept;
    }
    std::uniform_int_distribution<size_t> uni(0, y_count - 1);
    for (size_t i = 0; i < yprime_count; ++i) {
      const JetRec &rec = y_cache[uni(rng)];
      out_pt = rec.pt;
      out_pt_raw = rec.pt_raw;
      out_pt_shifted = rec.pt_shifted;
      out_eta = rec.eta;
      out_phi = rec.phi;
      out_mass = rec.mass;
      out_weight = rec.weight;
      out_source = 2;
      out_const_pt = rec.const_pt;
      out_const_eta = rec.const_eta;
      out_const_phi = rec.const_phi;
      // element-wise restore keeps out_subjet_pt (and the branch addresses) put
      for (size_t r = 0; r < nSubR && r < rec.subjet_pt.size(); ++r)
        out_subjet_pt[r] = rec.subjet_pt[r];
      tYp->Fill();
      tBiased->Fill();
    }
  }
  
  // first clone pT X 
  TH1D* hpTpp = (TH1D*)hPtX->Clone("hpTpp");
  // then add pT Y
  hpTpp->Add(hPtY);
  hpTpp->SetLineColor(kGreen+3);

  // build histograms from biased and reference trees to visualize the bias
  TH1D* hRef = (TH1D*)hPtX->Clone("hRef");
  hRef->Reset("ICES");
  hRef->SetTitle("; #it{p}_{T} (GeV); d#it{N}/d#it{p}_{T}");
  TH1D* hBiased = (TH1D*)hRef->Clone("hBiased");

  float tmp_pt = 0.0f;
  float tmp_w = 1.0f;

  tRef->SetBranchAddress("pt_raw", &tmp_pt);
  tRef->SetBranchAddress("weight", &tmp_w);
  for (Long64_t i = 0; i < tRef->GetEntries(); ++i) {
    tRef->GetEntry(i);
    hRef->Fill(tmp_pt, tmp_w);
  }

  tBiased->SetBranchAddress("pt_raw", &tmp_pt);
  tBiased->SetBranchAddress("weight", &tmp_w);
  for (Long64_t i = 0; i < tBiased->GetEntries(); ++i) {
    tBiased->GetEntry(i);
    hBiased->Fill(tmp_pt, tmp_w);
  }
  hRef->SetLineColor(kGreen+2);
  hBiased->SetLineColor(kMagenta+2);
  hRef->SetLineWidth(3);
  hBiased->SetLineWidth(3);

  // -----------------------------
  // Plot
  // -----------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas* c = new TCanvas("c", "Jet pT shift", 600, 600);
  c->SetLogy();
  c->SetTopMargin(0.05);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.05);
  c->SetTickx(1);
  c->SetTicky(1);

  hRef->Draw("hist");
  hBiased->Draw("hist same");
  hpTpp->Draw("hist same");

  TLegend* leg = new TLegend(0.50, 0.78, 0.88, 0.88);
  leg->AddEntry(hRef, "Unbiased (X #cup Y)", "l");
  leg->AddEntry(hBiased, "Biased (X down + Y')", "l");
  leg->SetBorderSize(0);
  leg->Draw();

  c->SaveAs(Form("jetPtShift_TChainReader_%s.pdf", "072226"));

  // Save histograms and trees for downstream unbiasing.
  hPtX->Write();
  hPtY->Write();
  hpTpp->Write();
  c->Write();
  tX->Write();
  tY->Write();
  tYp->Write();
  tRef->Write();
  tBiased->Write();
  tPP->Write();

  // Analytic weights for toy test
  if (!y_cache.empty()) {
    const double y = static_cast<double>(y_cache.size());
    const double yprime = static_cast<double>(tYp->GetEntries());
    TTree tW("tWeights", "analytic class weights");
    double wX_out  = yprime / (yprime + y);
    double wYp_out = y / (yprime + y);
    tW.Branch("wX", &wX_out, "wX/D");
    tW.Branch("wYprime", &wYp_out, "wYprime/D");
    tW.Fill();
    tW.Write();
  }
  // Sanity check: biased vs reference counts
  std::cout << "tRef entries: " << tRef->GetEntries()
            << "  tBiased entries: " << tBiased->GetEntries()
            << "  diff=" << (tBiased->GetEntries() - tRef->GetEntries())
            << std::endl;
  fout.Close();
}
