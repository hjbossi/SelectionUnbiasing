#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/Pythia8Yoda.h"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;
using namespace fastjet;
using namespace YODA;

int main() {

    // Set up PYTHIA to produce events at pp reference energy for the LHC
    Pythia pythia;
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 5360.");
    pythia.readString("HardQCD:all = on");//jets
    pythia.readString("PhaseSpace:pTHatMin = 50."); // choose a lower pthatmin to avoid edge effects
    pythia.readString("PhaseSpace:pTHatMax = 200.");
    pythia.init();

    // now set uyp the yoda part following pythia example 114
    Pythia8Yoda p8y("PPJETS_YODA", "ppjets_yoda.yoda");


    // Set up jet clustering
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);

    // Booking of normal types, as defined in Pythia8Yoda.
    auto h_jet_pt         = p8y.bookHisto1D(100, 0.0, 200.0, "jet-pt");
    auto h_jet_eta        = p8y.bookHisto1D(20, -5.0, 5.0, "jet-eta");
    auto h_jet_phi        = p8y.bookHisto1D(100, 0.0,6.5, "jet-phi");
    auto h_constituent_pt = p8y.bookHisto1D(100, 0.0, 50.0, "jet-constituent-pt"); 
    auto h_constituent_eta= p8y.bookHisto1D(20, -5.0, 5.0, "jet-constituent-eta");
    auto h_constituent_phi= p8y.bookHisto1D(100, 0.0,6.5, "jet-constituent-phi");
    
    // Event loop
    int nEvents = 1000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        // Collect final-state particles
        std::vector<PseudoJet> particles;
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (!pythia.event[i].isFinal()) continue;
            if (!pythia.event[i].isVisible()) continue;
            particles.emplace_back(
                pythia.event[i].px(),
                pythia.event[i].py(),
                pythia.event[i].pz(),
                pythia.event[i].e()
            );
        }

        // Cluster jets
        ClusterSequence cs(particles, jet_def);
        std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(20.0));

        // Fill histograms
        for (const auto& jet : jets) {
            h_jet_pt->fill(jet.pt());
	    h_jet_eta->fill(jet.eta());
	    h_jet_phi->fill(jet.phi()); 
            for (const auto& consti : jet.constituents()) {
                h_constituent_pt->fill(consti.pt());
		h_constituent_eta->fill(consti.eta());
		h_constituent_phi->fill(consti.phi()); 
            }
        }
    }

    // write to a file
    p8y.write();
    
    // Done
    pythia.stat();
    return 0;
}
