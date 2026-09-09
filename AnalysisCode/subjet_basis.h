// subjet_basis.h
// =====================================================================
// Self-contained Cambridge/Aachen subjet reclustering.
//
// The basis-function framework itself now lives in basis_functions.h, and
// the standard workflow consumes the subjets that ppjets_root.cc already
// clustered at ntuple-production time (branches "subjet_pt_R0pXX").  This
// header is only the reclustering FALLBACK, used by unbias_weights.cc for
// older constituent-only ntuples (--subjet-source recluster, or "auto"
// when the stored branch is absent).
//
// The C/A reclustering is implemented here without a FastJet dependency,
// so the tools build with root-config alone.  Constituents are treated as
// massless (only pt/eta/phi are stored) and merged with the E-scheme
// (4-momentum addition).
// =====================================================================
#pragma once

#include <TVector2.h>

#include <cmath>
#include <limits>
#include <vector>

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
/// return their pT values.
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
