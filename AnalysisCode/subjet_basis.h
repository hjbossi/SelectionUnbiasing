// subjet_basis.h
// =====================================================================
// Shared basis definition for the selection-unbiasing chain.
//
// New basis (replaces the old EEC-style ΔR^{-A} ln^B z^m functions):
//
//   For each R=0.4 jet, recluster its constituents into R=0.1 subjets
//   with the Cambridge/Aachen (C/A) algorithm, then form
//
//       g_n = Σ_{subjets} (pT_subjet)^n ,   n = 3, 4, ..., 10
//
//   i.e. one basis function per exponent (8 functions total).
//
// The C/A reclustering is implemented self-contained here (no FastJet
// dependency) so both unbias_weights.cc and plot_unbias_weights.cc can
// be built with root-config alone and stay bit-for-bit identical.
//
// Constituents are treated as massless (only pt/eta/phi are stored),
// and merged with the E-scheme (4-momentum addition).
// =====================================================================
#pragma once

#include <TVector2.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

// ---- Subjet reclustering configuration ------------------------------
static constexpr double kSubjetR = 0.1;   // C/A subjet radius
static constexpr int    kBasisNMin = 3;   // lowest pT exponent
static constexpr int    kBasisNMax = 10;  // highest pT exponent

// =====================================================================
// Basis function bookkeeping
// =====================================================================
// Each basis function is g_n = Σ_subjets pT^n.  We keep the same struct
// name / accessor names the fit and plot code already used so that the
// surrounding machinery (nphys = basis.size(), basis[j].label, ...) is
// unchanged.
struct BasisFuncDef {
    int         n;      // pT exponent
    std::string label;
};

static std::vector<BasisFuncDef> get_default_basis() {
    std::vector<BasisFuncDef> basis;
    for (int n = kBasisNMin; n <= kBasisNMax; ++n) {
        basis.push_back({n, "g" + std::to_string(n) +
                             ": sum pT_subjet^" + std::to_string(n) +
                             " (C/A R=0.1)"});
    }
    return basis;
}

// =====================================================================
// Self-contained Cambridge/Aachen reclustering
// =====================================================================
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

// =====================================================================
// BASIS FUNCTION EVALUATION
// =====================================================================
/// Evaluate all basis functions g_n = Σ_subjets pT^n for one jet.
/// R_subjet is the C/A reclustering radius (default kSubjetR = 0.1).
static std::vector<double> evaluate_basis(
        const std::vector<BasisFuncDef> &basis,
        double ptjet,
        const std::vector<double> &cpt,
        const std::vector<double> &ceta,
        const std::vector<double> &cphi,
        double R_subjet = kSubjetR)
{
    const int nb = (int)basis.size();
    std::vector<double> g(nb, 0.0);
    if (ptjet <= 0.0) return g;

    const std::vector<double> sjpt =
        recluster_ca_subjet_pts(cpt, ceta, cphi, R_subjet);

    for (const double pt : sjpt) {
        if (pt <= 0.0) continue;
        for (int j = 0; j < nb; ++j)
            g[j] += std::pow(pt, (double)basis[j].n);
    }
    return g;
}
