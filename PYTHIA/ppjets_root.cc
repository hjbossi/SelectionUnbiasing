// ppjets_root.cc: code to generate pp jets are the reference com energy for unbiasing studies
// Hannah Bossi, <hannah.bossi@cern.ch>
// December 18th, 2025


// --------------------- default preamble below -----------------------------
// based off of main143.cc, which is one of the standard pythia example codes
// main143.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Rene Brun, Axel Naumann and Bernhard Meirose

// Keywords: analysis; root

// This is a simple test program, based on main101.cc,
// but modified to use ROOT for histogramming.
// It studies the charged multiplicity distribution at the LHC.

// WARNING: for currently unknown reasons it may hang
// with an empty canvas on a Mac.

// ---------------------------------------------------------------------------

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, TTree for writing output
#include "TTree.h"
#include "TMath.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// jet stuff
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

// ROOT, for saving file.
#include "TFile.h"
#define MAXJETS 100
#define MAXCONST 100


using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  int nEvent = 10;//1e5;
  
  // pp beams
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  // set to reference center of mass energy for PbPb LHC run 2
  pythia.readString("Beams:eCM = 5360.");
  pythia.readString("HardQCD:all = on");//jets
  pythia.readString("PhaseSpace:pTHatMin = 50."); // choose a lower pthatmin to avoid edge effects
  pythia.readString("PhaseSpace:pTHatMax = 200.");
  pythia.init();
  
  
  // keeping some options here commented out just in case. 
  
  // vincia = 2 // dire = 3
  //pythia.readString("PartonShowers:Model = 3");


  // SWITCH TO TURN ON/OFF Hadronization, default is on
  // pythia.readString("HadronLevel:Hadronize=  off");
  
  // set alpha_strong value at scale M_Z^2.
  //pythia.readString("SigmaProcess:alphaSvalue= 0.1136"); // value from https://arxiv.org/abs/2412.15164
  //pythia.readString("SigmaProcess:alphaSvalue= 0.1180"); // value from world average https://pdg.lbl.gov/2024/reviews/rpp2024-rev-qcd.pdf

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("UnbiasingTest_PYTHIApp_pthatmin50_071125.root", "RECREATE");
  TTree*    tree = new TTree("tgenBefore","Pythia8 event tree");
  
  // ------------------ Define the tree ----------------
  // jet variables 
  Int_t nJets; // to track the number of this part of the tree
  float pt[MAXJETS];
  float y[MAXJETS];
  float eta[MAXJETS];
  float phi[MAXJETS];
  float mass[MAXJETS];
  Int_t nSubJets[MAXJETS];
  
  // consituent variables - have this be a vector of vectors the size of the number of jets 
  std::vector<std::vector<double>> const_pt; 
  std::vector<std::vector<double>> const_eta; 
  std::vector<std::vector<double>> const_phi; 
  
  // subjet variables - have this be a vector of vectors the size of the number of subjets
  std::vector<std::vector<double>> subjet_pt; 
  std::vector<std::vector<double>> subjet_eta; 
  std::vector<std::vector<double>> subjet_phi; 
  
  // event variables 
  float weight;
  
  // --------------------------------------------------------
  
  
  // ------------------ create the branches ------------------

  // branch for jet variables 
  tree->Branch("nJets",&nJets,"nJets/I");
  tree->Branch("pt",pt,"pt[nJets]/F");
  tree->Branch("eta",eta,"eta[nJets]/F");
  tree->Branch("phi",phi,"phi[nJets]/F");
  tree->Branch("mass",mass, "mass[nJets]/F"); 
  tree->Branch("nSubJets",nSubJets, "nSubJets[nJets]/I"); 

  // branch for constituent variables
  tree->Branch("const_pt", &const_pt); 
  tree->Branch("const_eta", &const_eta); 
  tree->Branch("const_phi", &const_phi); 
  
  // branch for subjet variables
  tree->Branch("subjet_pt", &subjet_pt); 
  tree->Branch("subjet_eta", &subjet_eta); 
  tree->Branch("subjet_phi", &subjet_phi); 

  // branch for event 
  tree->Branch("weight",&weight,"weight/F");
  // -----------------------------------------------------------

  
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
    if (!pythia.next()) continue;
    
    if (iEvent % 100000 == 0) std::cout << "On event " << iEvent << std::endl;

    // Find number of all final charged particles.
    weight = pythia.info.weight(); 
    
    const_pt.clear();
    const_eta.clear();
    const_phi.clear();
    
    subjet_pt.clear();
    subjet_eta.clear();
    subjet_phi.clear();
    
    // use a classic jet definition - antiKT jets - 0.4 
    fastjet::JetDefinition antiKT = fastjet::JetDefinition( fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
    fastjet::JetDefinition CA3    = fastjet::JetDefinition( fastjet::cambridge_algorithm, 0.3, fastjet::E_scheme, fastjet::Best);

    
        // Collect final-state particles for jet stuff
        std::vector<fastjet::PseudoJet> particles;
        for (int i = 0; i < pythia.event.size(); ++i) {
            // make some basic cuts on the particles, may need to modify this
            if (!pythia.event[i].isFinal()) continue;
            if (!pythia.event[i].isVisible()) continue;
            
            // fill the particles pseudojet vector for jet finding in the next step
            particles.emplace_back(
                pythia.event[i].px(),
                pythia.event[i].py(),
                pythia.event[i].pz(),
                pythia.event[i].e()
            );
        }

        // Cluster jets
        fastjet::ClusterSequence cs(particles, antiKT);
        std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(20.0));
        
        nJets = jets.size(); 
        int index = 0; 
        

        
        for (size_t i = 0; i < jets.size(); ++i) {
          const fastjet::PseudoJet &jet = jets[i];
          pt[index]   = jet.pt(); 
          eta[index]  = jet.eta(); 
          phi[index]  = jet.phi(); 
          mass[index] = jet.m(); 
          std::cout << "pt jet: " << jet.pt() << std::endl;
          
          // for each jet loop through the constituents
          std::vector<fastjet::PseudoJet> constituents = jet.constituents();
          
          std::vector<double> jetConst_pt; 
          std::vector<double> jetConst_eta;
          std::vector<double> jetConst_phi;
      
          for (const auto &c : constituents) {
            jetConst_pt.push_back(c.pt()); 
            jetConst_eta.push_back(c.eta()); 
            jetConst_phi.push_back(c.phi()); 
          }
          const_pt.emplace_back(jetConst_pt); 
          const_eta.emplace_back(jetConst_eta); 
          const_phi.emplace_back(jetConst_phi);  
          
          // recluster the constituents into C/A jets
          fastjet::ClusterSequence cs2(constituents, CA);
          std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(cs2.inclusive_jets(0.0));
          nSubJets[index] = subjets.size();
          std::vector<double> subjetConst_pt; 
          std::vector<double> subjetConst_eta;
          std::vector<double> subjetConst_phi;
      
          for (const auto &s : subjets) {
            subjetConst_pt.push_back(s.pt()); 
            subjetConst_eta.push_back(s.eta()); 
            subjetConst_phi.push_back(s.phi()); 
            std::cout << "pt subjet: " << s.pt() << std::endl;
          }
          
          subjet_pt.emplace_back(subjetConst_pt); 
          subjet_eta.emplace_back(subjetConst_eta); 
          subjet_phi.emplace_back(subjetConst_phi); 
                    
          index++;
       } // end loop over the jets
       tree->Fill(); 
        
    } // end the loop over the number of events
  

  // Statistics on event generation.
  //pythia.stat();

  tree->Write(); 
  delete outFile;

  // Done.
  return 0;
}
