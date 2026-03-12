// startBasis.C: file to plot the initial basis and see the effect of a pT shift on these distirbutions
// Hannah Bossi, <hannah.bossi@cern.ch>
// December 19th, 2025

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>



#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>

#include <vector>
#include <string>
#include <iostream>

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
void startBasis() {

  // config
  const char* inputDir =  "/home/hbossi/SelectionUnbiasing/MCOutput/";  // <-- directory containing ROOT files
  const char* treeName = "tgenBefore";
  const char* outFile  = "startBasis_output.root";

  // target (unbiased) window
  const double pTWindowXLow  = 80;
  const double pTWindowXHigh = 140;

  // biased window (shifted up)
  const double ptShift = 10.0;     
  
  
  // then the window Y becomes
  const double pTWindowYLow  = pTWindowXLow  + ptShift; 
  const double pTWindowYHigh = pTWindowXHigh + ptShift; 

  
  // get the files and 
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
  TTree *tY = new TTree("tY", "biased Y jets (pT-shifted)");
  TTree *tPP = new TTree("tPP", "unbiased sample (X + Y)");

  float out_pt = 0.0f;
  float out_pt_raw = 0.0f;
  float out_pt_shifted = 0.0f;
  float out_eta = 0.0f;
  float out_phi = 0.0f;
  float out_mass = 0.0f;
  float out_weight = 1.0f;
  int out_source = 0; // 0 = X, 1 = Y
  std::vector<double> out_const_pt;
  std::vector<double> out_const_eta;
  std::vector<double> out_const_phi;

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
  };
  setup_tree(tX);
  setup_tree(tY);
  setup_tree(tPP);

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

      // X histos
      // only fill these histograms if it is within the X window
      if((pTWindowXLow < ptRaw) && (ptRaw < pTWindowXHigh)){
            hPtX->Fill(ptRaw, *weight);
            out_pt = ptRaw;
            out_pt_raw = ptRaw;
            out_pt_shifted = ptRaw + ptShift;
            out_eta = eta[j];
            out_phi = phi[j];
            out_mass = mass[j];
            out_weight = *weight;
            out_source = 0;
            out_const_pt = (*const_pt)[j];
            out_const_eta = (*const_eta)[j];
            out_const_phi = (*const_phi)[j];
            tX->Fill();
            tPP->Fill();
      }
  
      
      // // for each jet, calculate the basis vector
      // int nConst = const_pt->size();
      // for (size_t j = 0; j < nConst; ++j) {

      // }
      

      // Y shift (absolute)
      double ptShifted = ptRaw + ptShift;

      // alternative fractional shift
      // double ptShifted = ptRaw * ptFrac;

      if((pTWindowYLow < ptShifted) && (ptShifted < pTWindowYHigh)){
          hPtY->Fill(ptShifted, *weight);
          out_pt = ptShifted;
          out_pt_raw = ptRaw;
          out_pt_shifted = ptShifted;
          out_eta = eta[j];
          out_phi = phi[j];
          out_mass = mass[j];
          out_weight = *weight;
          out_source = 1;
          out_const_pt = (*const_pt)[j];
          out_const_eta = (*const_eta)[j];
          out_const_phi = (*const_phi)[j];
          tY->Fill();
          tPP->Fill();
      }
    }
    
    // now create the pp and the AA samples by combining X and Y
    // for now, let's just create the pp sample (unbiased)

  
  }
  
  // first clone pT X 
  TH1D* hpTpp = (TH1D*)hPtX->Clone("hpTpp"); 
  // then add pT Y
  hpTpp->Add(hPtY); 
  hpTpp->SetLineColor(kGreen+3);




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

  hPtX->Draw("hist");
  hPtY->Draw("hist same");
  hpTpp->Draw("hist same"); 

  TLegend* leg = new TLegend(0.50, 0.75, 0.88, 0.88);
  leg->AddEntry(hPtX, "Original jets", "l");
  leg->AddEntry(hPtY, Form("Shifted jets (#Delta #it{p}_{T} = %.1f GeV)", ptShift),"l"  );
  leg->AddEntry(hpTpp, "Unbiased sample","l"  );
  leg->SetBorderSize(0);
  leg->Draw();

  c->SaveAs("jetPtShift_TChainReader.pdf");

  // Save histograms and trees for downstream unbiasing.
  hPtX->Write();
  hPtY->Write();
  hpTpp->Write();
  c->Write();
  tX->Write();
  tY->Write();
  tPP->Write();
  fout.Close();
}
