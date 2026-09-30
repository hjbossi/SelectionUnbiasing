// subjet_basis.h
// =====================================================================
// Subjet reclustering: two interchangeable engines.
//
// The basis-function framework itself lives in basis_functions.h, and the
// standard workflow consumes the subjets that ppjets_root.cc already
// clustered at ntuple-production time (branches "subjet_pt_R0pXX").  This
// header supplies the RECLUSTERING path, used by unbias_weights.cc for
// older constituent-only ntuples, or whenever asked for explicitly via
// --subjet-source recluster / --subjet-source fastjet ("auto" uses
// reclustering only when the stored branch is absent, and currently
// falls back to the homemade engine in that case).
//
// Two engines are provided, selected independently of each other so they
// can be run side by side as a crosscheck of implementation:
//
//   * recluster_ca_subjet_pts(...)
//       The original dependency-free Cambridge/Aachen reclustering,
//       implemented here without any FastJet dependency.  Constituents
//       are treated as massless (only pt/eta/phi are stored) and merged
//       with the E-scheme (4-momentum addition).  Kept exactly as before
//       and retained on its own as a crosscheck.
//
//   * recluster_fastjet_subjet_pts(...)
//       The same Cambridge/Aachen algorithm and radius, but clustered by
//       the standalone FastJet library itself (fastjet::ClusterSequence
//       with fastjet::cambridge_algorithm), reading inclusive jets.  Same
//       massless (pt,eta,phi) constituents and E-scheme recombination, so
//       the two engines are a genuine crosscheck of implementation, not
//       of algorithm choice.
//
// Because this file now offers the FastJet engine, any translation unit
// that includes subjet_basis.h requires linking against FastJet, e.g.
//   c++ -std=c++17 -O2 unbias_weights.cc \
//       $(root-config --cflags --libs) \
//       $(/data/ALEPH/MC/mcgen/fastjet-install/bin/fastjet-config --cxxflags --libs) \
//       -o unbias_weights
// =====================================================================
#pragma once

#include <TVector2.h>

#include <cmath>
#include <limits>
#include <vector>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>

// A pseudojet carrying a 4-momentum for E-scheme recombination.
struct SubjetPseudo {
    double px, py, pz, E;
};

static inline SubjetPseudo make_pseudo(double pt, double eta, double phi) {
    SubjetPseudo p;
    p.px = pt * std::cos(phi);
    p.py = pt * std::sin(phi);
    p.pz = pt * std::sinh(eta);
    p.E  = pt * std::cosh(eta);   // massless constituent
    return p;
}

static inline double pseudo_pt(const SubjetPseudo &p) {
    return std::sqrt(p.px * p.px + p.py * p.py);
}

static inline double pseudo_eta(const SubjetPseudo &p) {
    const double pt = pseudo_pt(p);
    if (pt <= 0.0) return (p.pz >= 0.0) ? 1e6 : -1e6;
    return std::asinh(p.pz / pt);   // pseudorapidity
}

static inline double pseudo_phi(const SubjetPseudo &p) {
    return std::atan2(p.py, p.px);
}

/// Recluster jet constituents into inclusive C/A subjets of radius R and
/// return their pT values.  Homemade, dependency-free engine.
///
/// C/A distances: d_ij = ΔR_ij² / R², d_iB = 1.  Since d_iB is fixed at
/// 1, the closest pair merges iff its ΔR < R; once no pair is within R,
/// every remaining pseudojet is a final subjet.
static std::vector<double> recluster_ca_subjet_pts(
        const std::vector<double> &cpt,
        const std::vector<double> &ceta,
        const std::vector<double> &cphi,
        double R)
{
    std::vector<SubjetPseudo> jets;
    jets.reserve(cpt.size());
    const size_t nc = cpt.size();
    for (size_t i = 0; i < nc; ++i) {
        if (cpt[i] <= 0.0) continue;
        if (i >= ceta.size() || i >= cphi.size()) break;
        jets.push_back(make_pseudo(cpt[i], ceta[i], cphi[i]));
    }

    const double R2 = R * R;
    while (jets.size() > 1) {
        // Find the geometrically closest pair (smallest ΔR²).
        double best_d2 = std::numeric_limits<double>::infinity();
        size_t bi = 0, bj = 0;
        for (size_t i = 0; i < jets.size(); ++i) {
            const double eta_i = pseudo_eta(jets[i]);
            const double phi_i = pseudo_phi(jets[i]);
            for (size_t j = i + 1; j < jets.size(); ++j) {
                const double deta = eta_i - pseudo_eta(jets[j]);
                const double dphi = TVector2::Phi_mpi_pi(phi_i - pseudo_phi(jets[j]));
                const double d2   = deta * deta + dphi * dphi;
                if (d2 < best_d2) { best_d2 = d2; bi = i; bj = j; }
            }
        }
        if (best_d2 >= R2) break;   // no pair within R -> all remaining are subjets

        // Merge bi and bj (E-scheme), removing the higher index first.
        SubjetPseudo merged;
        merged.px = jets[bi].px + jets[bj].px;
        merged.py = jets[bi].py + jets[bj].py;
        merged.pz = jets[bi].pz + jets[bj].pz;
        merged.E  = jets[bi].E  + jets[bj].E;
        jets[bi] = merged;
        jets.erase(jets.begin() + bj);
    }

    std::vector<double> pts;
    pts.reserve(jets.size());
    for (const auto &j : jets) pts.push_back(pseudo_pt(j));
    return pts;
}

/// Recluster jet constituents into inclusive C/A subjets of radius R and
/// return their pT values.  Standalone-FastJet engine: same algorithm
/// (Cambridge/Aachen), same radius, same massless (pt,eta,phi) inputs and
/// E-scheme recombination as recluster_ca_subjet_pts(...) above, but the
/// clustering itself is performed by fastjet::ClusterSequence rather than
/// the homemade loop -- a crosscheck of implementation, not of algorithm.
static std::vector<double> recluster_fastjet_subjet_pts(
        const std::vector<double> &cpt,
        const std::vector<double> &ceta,
        const std::vector<double> &cphi,
        double R)
{
    std::vector<fastjet::PseudoJet> constituents;
    constituents.reserve(cpt.size());
    const size_t nc = cpt.size();
    for (size_t i = 0; i < nc; ++i) {
        if (cpt[i] <= 0.0) continue;
        if (i >= ceta.size() || i >= cphi.size()) break;
        const SubjetPseudo p = make_pseudo(cpt[i], ceta[i], cphi[i]);
        constituents.emplace_back(p.px, p.py, p.pz, p.E);   // massless, E-scheme
    }

    const fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, R);
    fastjet::ClusterSequence cs(constituents, jet_def);
    const std::vector<fastjet::PseudoJet> subjets = cs.inclusive_jets();

    std::vector<double> pts;
    pts.reserve(subjets.size());
    for (const auto &j : subjets) pts.push_back(j.pt());
    return pts;
}