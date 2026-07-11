// basis_functions.h
// =====================================================================
// Configurable basis-function framework for the selection-unbiasing fit.
//
// This header replaces the fixed, single-family basis that used to live
// in `subjet_basis.h` (C/A R=0.1 subjet pT power sums, reconstructed by
// reclustering the jet constituents on the fly).  The analysis chain has
// since been updated so that the subjets are *already clustered* during
// ROOT-file production (`PYTHIA/ppjets_root.cc`) and stored per radius in
// the branches
//
//     subjet_pt_R0pXX , subjet_eta_R0pXX , subjet_phi_R0pXX
//
// with 0pXX == 0.XX (R = 0.01 ... 0.20 in steps of 0.01).  The fit should
// therefore *consume the stored subjets directly* rather than recluster.
//
// The framework supports an arbitrary basis built from any combination of
//
//   (1) Energy-correlator (EEC) basis functions
//         g = Sum_{i<k, dR<E}  dR^{-A} [ln dR]^B  z^m ,   z = pt_i pt_k / norm^2
//       (the historical basis from the older unbias_weights.cc), and
//
//   (2) Subjet-pT moment basis functions
//         g = Sum_{subjets at radius R}  pT_subjet^n
//       i.e. Mellin moments of the subjet pT spectrum at a chosen subjet
//       radius R.  The user picks one or more radii and one or more powers.
//
// DESIGN
// ------
//   * `BasisFunction`  -- abstract base: one scalar moment per jet, plus a
//                         declaration of which inputs it needs (constituents
//                         and/or which subjet radii) so the reader only loads
//                         the branches that are actually used.
//   * `EECBasisFunction`, `SubjetMomentBasisFunction` -- the two families.
//   * `JetData`        -- a per-jet input bundle (constituents + a lazily
//                         built pair table + the stored subjet lists keyed by
//                         radius).  Basis functions read from it; they never
//                         touch ROOT, so this header is pure C++ and unit
//                         testable (define BASIS_NO_ROOT to drop the ROOT dep).
//   * `BasisConfig` / `build_basis()` -- one clearly-marked place to declare
//                         the desired basis.  Adding a term is a one-line edit;
//                         adding a whole new *family* is a new BasisFunction
//                         subclass and one line in `build_basis`.
//
// The pairwise EEC geometry (dR, ln dR, pt_i*pt_k) is computed once per jet in
// `build_constituent_pairs()` and shared by every EEC term, so no expensive
// transcendental is recomputed and there is no duplicated pair-loop logic.
// =====================================================================
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------
// phi wrapping: use ROOT's TVector2 by default, but allow a ROOT-free
// build (for unit tests) via -DBASIS_NO_ROOT.
// ---------------------------------------------------------------------
#ifdef BASIS_NO_ROOT
namespace basis_detail {
inline double phi_mpi_pi(double dphi) {
    while (dphi >   M_PI) dphi -= 2.0 * M_PI;
    while (dphi <= -M_PI) dphi += 2.0 * M_PI;
    return dphi;
}
}  // namespace basis_detail
#else
#include <TVector2.h>
namespace basis_detail {
inline double phi_mpi_pi(double dphi) { return TVector2::Phi_mpi_pi(dphi); }
}  // namespace basis_detail
#endif

// =====================================================================
// Per-jet input bundle
// =====================================================================
// A precomputed constituent pair, shared by all EEC-type basis functions.
struct ConstituentPair {
    double dR;       // angular distance
    double lndR;     // ln(dR)   (precomputed once)
    double pt_prod;  // pt_i * pt_k   (z = pt_prod / norm^2, norm set per term)
};

// Everything a basis function might read for a single jet.  The reader fills
// only the fields the configured basis actually needs (see BasisInputs).
struct JetData {
    double pt = 0.0;  // jet pT

    // Raw constituent kinematics (optional; owned by the caller).
    const std::vector<double>* const_pt  = nullptr;
    const std::vector<double>* const_eta = nullptr;
    const std::vector<double>* const_phi = nullptr;

    // Pairwise geometry, built once per jet if any EEC term is present.
    std::vector<ConstituentPair> pairs;

    // Stored (pre-clustered) subjet pT lists, keyed by radius tag
    // rtag = round(R * 100)  ->  e.g. R=0.10 -> 10 -> branch "..._R0p10".
    std::map<int, const std::vector<double>*> subjet_pt;

    // ---- helpers -----------------------------------------------------
    static int radius_tag(double R) { return static_cast<int>(std::lround(R * 100.0)); }

    const std::vector<double>* subjets_at(double R) const {
        auto it = subjet_pt.find(radius_tag(R));
        return (it == subjet_pt.end()) ? nullptr : it->second;
    }

    void reset() {
        pt = 0.0;
        const_pt = const_eta = const_phi = nullptr;
        pairs.clear();
        subjet_pt.clear();
    }
};

// Branch-name helper matching the encoding used in ppjets_root.cc:
//   subjet_<field>_R0pXX   with XX = round(R*100), zero padded to 2 digits.
inline std::string subjet_branch_name(const std::string& field, double R) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "subjet_%s_R0p%02d", field.c_str(),
                  JetData::radius_tag(R));
    return std::string(buf);
}

// Convert a radius tag (round(R*100)) back to a radius, e.g. 10 -> 0.10.
inline double radius_from_tag(int rtag) { return rtag / 100.0; }

// Build the shared pair table for a jet (dR-screened at dR_min).  Single
// source of truth for constituent-pair geometry.
inline void build_constituent_pairs(JetData& jet, double dR_min) {
    jet.pairs.clear();
    if (!jet.const_pt || !jet.const_eta || !jet.const_phi) return;
    const std::vector<double>& pt  = *jet.const_pt;
    const std::vector<double>& eta = *jet.const_eta;
    const std::vector<double>& phi = *jet.const_phi;
    const size_t nc = pt.size();
    if (eta.size() < nc || phi.size() < nc) return;
    jet.pairs.reserve(nc * (nc - 1) / 2);
    for (size_t i = 0; i < nc; ++i) {
        if (pt[i] <= 0.0) continue;
        for (size_t k = i + 1; k < nc; ++k) {
            if (pt[k] <= 0.0) continue;
            const double dphi = basis_detail::phi_mpi_pi(phi[i] - phi[k]);
            const double deta = eta[i] - eta[k];
            const double dR   = std::sqrt(deta * deta + dphi * dphi);
            if (dR <= dR_min) continue;  // screen small-angle divergence
            jet.pairs.push_back({dR, std::log(dR), pt[i] * pt[k]});
        }
    }
}

// =====================================================================
// Abstract basis function
// =====================================================================
// One scalar moment g(jet) per basis function.  Subclasses also declare
// which inputs they consume so the reader can wire up only those branches.
class BasisFunction {
public:
    virtual ~BasisFunction() = default;

    // The moment value for one jet.
    virtual double evaluate(const JetData& jet) const = 0;

    // Human-readable description (printed in the fit log / stored as metadata).
    virtual std::string label() const = 0;

    // ---- input requirements -----------------------------------------
    // Does this function read the constituent pair table?
    virtual bool needs_constituents() const { return false; }
    // Which stored subjet radii (as rtags) does this function read?
    virtual void collect_radii(std::set<int>& rtags) const { (void)rtags; }
};

using Basis = std::vector<std::unique_ptr<BasisFunction>>;

// =====================================================================
// (1) Energy-correlator basis function
//     g = Sum_{pairs, dR<E}  dR^{-A} [ln dR]^B  (pt_i pt_k / norm^2)^m
// =====================================================================
class EECBasisFunction : public BasisFunction {
public:
    EECBasisFunction(double A, double B, int m, double E,
                     double norm, std::string label)
        : A_(A), B_(B), m_(m), E_(E), norm2_(norm * norm),
          label_(std::move(label)) {}

    bool needs_constituents() const override { return true; }

    double evaluate(const JetData& jet) const override {
        double g = 0.0;
        for (const ConstituentPair& p : jet.pairs) {
            if (p.dR >= E_) continue;              // angular upper cut
            double val = std::pow(p.dR, -A_);
            if (B_ > 0.0) val *= std::pow(p.lndR, B_);
            const double z = p.pt_prod / norm2_;
            val *= std::pow(z, static_cast<double>(m_));
            g += val;
        }
        return g;
    }

    std::string label() const override { return label_; }

private:
    double      A_;      // dR^{-A}
    double      B_;      // [ln dR]^B
    int         m_;      // z^m
    double      E_;      // angular upper cut: keep dR < E
    double      norm2_;  // z = pt_i pt_k / norm2_
    std::string label_;
};

// =====================================================================
// (2) Subjet-pT moment basis function
//     g = Sum_{subjets at radius R}  pT_subjet^n
//     (Mellin moment of the stored subjet-pT spectrum at radius R.)
// =====================================================================
class SubjetMomentBasisFunction : public BasisFunction {
public:
    SubjetMomentBasisFunction(double R, double n, std::string label)
        : R_(R), n_(n), rtag_(JetData::radius_tag(R)), label_(std::move(label)) {}

    void collect_radii(std::set<int>& rtags) const override { rtags.insert(rtag_); }

    double evaluate(const JetData& jet) const override {
        auto it = jet.subjet_pt.find(rtag_);
        if (it == jet.subjet_pt.end() || it->second == nullptr) return 0.0;
        double g = 0.0;
        for (const double pt : *it->second)
            if (pt > 0.0) g += std::pow(pt, n_);
        return g;
    }

    std::string label() const override { return label_; }

    double  radius() const { return R_; }
    int     radius_tag() const { return rtag_; }
    double  power() const { return n_; }

private:
    double      R_;      // subjet radius (selects the stored branch)
    double      n_;      // pT power (double -> allows fractional Mellin moments)
    int         rtag_;   // round(R*100)
    std::string label_;
};

// =====================================================================
// Basis configuration + builder
// =====================================================================
// One EEC term.  (A, B, m, E) match the historical BasisFuncDef fields.
struct EECTermSpec {
    double      A;
    double      B;
    int         m;
    double      E;
    std::string label;   // optional; auto-generated if empty
};

// A group of subjet-pT moments to include at a single radius.
struct SubjetMomentSpec {
    double              R;        // subjet radius (must be a stored radius)
    std::vector<double> powers;   // e.g. {3,4,5,6,7,8,9,10}
};

// -------- label helpers ----------------------------------------------
inline std::string make_eec_label(const EECTermSpec& t) {
    char buf[128];
    std::snprintf(buf, sizeof(buf),
                  "eec: dR^{-%.3g} ln^{%.3g} z^{%d} [dR<%.3g]", t.A, t.B, t.m, t.E);
    return std::string(buf);
}

inline std::string make_subjet_moment_label(double R, double n) {
    char buf[128];
    std::snprintf(buf, sizeof(buf),
                  "sjmom: sum pT^{%.3g} (stored subjets R=%.2f)", n, R);
    return std::string(buf);
}

// -------- the historical EEC basis (ported verbatim from the old
//          unbias_weights.cc get_default_basis) --------------------
inline std::vector<EECTermSpec> default_eec_terms() {
    return {
        {  1.0, 4.0, 1, 0.4, "" }, {  1.0, 3.0, 1, 0.4, "" },
        {  1.0, 2.0, 1, 0.4, "" }, {  1.0, 1.0, 1, 0.4, "" },
        {  1.0, 0.0, 1, 0.4, "" },
        {  0.0, 4.0, 1, 0.4, "" }, {  0.0, 3.0, 1, 0.4, "" },
        {  0.0, 2.0, 1, 0.4, "" }, {  0.0, 1.0, 1, 0.4, "" },
        {  0.0, 0.0, 1, 0.4, "" },
        { -1.0, 4.0, 1, 0.4, "" }, { -1.0, 3.0, 1, 0.4, "" },
        { -1.0, 2.0, 1, 0.4, "" }, { -1.0, 1.0, 1, 0.4, "" },
        { -1.0, 0.0, 1, 0.4, "" },
        {  1.0, 4.0, 2, 0.4, "" }, {  1.0, 3.0, 2, 0.4, "" },
        {  1.0, 2.0, 2, 0.4, "" }, {  1.0, 1.0, 2, 0.4, "" },
        {  1.0, 0.0, 2, 0.4, "" },
        {  0.0, 4.0, 2, 0.4, "" }, {  0.0, 3.0, 2, 0.4, "" },
        {  0.0, 2.0, 2, 0.4, "" }, {  0.0, 1.0, 2, 0.4, "" },
        {  0.0, 0.0, 2, 0.4, "" },
        { -1.0, 4.0, 2, 0.4, "" }, { -1.0, 3.0, 2, 0.4, "" },
        { -1.0, 2.0, 2, 0.4, "" }, { -1.0, 1.0, 2, 0.4, "" },
        { -1.0, 0.0, 2, 0.4, "" },
    };
}

// =====================================================================
//  >>> EDIT HERE to change the basis <<<
// ---------------------------------------------------------------------
// The default reproduces the current analysis basis (C/A R=0.1 subjet pT
// power sums, n = 3..10) but now sourced from the STORED subjets rather
// than on-the-fly reclustering.  Turn on `use_eec` to add the historical
// energy-correlator functions, add more {R, powers} groups for extra
// subjet radii, etc.  unbias_weights.cc can also override these fields
// from the command line.
// =====================================================================
struct BasisConfig {
    // ---- Energy-correlator family ----
    bool                     use_eec    = false;         // off -> subjet-only default
    double                   eec_norm   = 120.0;         // z = pt_i pt_k / eec_norm^2
    double                   eec_dR_min = 0.001;         // small-angle screen (shared)
    std::vector<EECTermSpec> eec_terms  = default_eec_terms();

    // ---- Subjet-pT moment family ----
    // Any number of radii, each with its own list of powers.
    std::vector<SubjetMomentSpec> subjet_moments = {
        { 0.10, { 3, 4, 5, 6, 7, 8, 9, 10 } },
    };
};

// Assemble the basis from a configuration.  Order: EEC terms first (if
// enabled), then subjet moments grouped by radius.
inline Basis build_basis(const BasisConfig& cfg) {
    Basis basis;

    if (cfg.use_eec) {
        for (const EECTermSpec& t : cfg.eec_terms) {
            basis.push_back(std::make_unique<EECBasisFunction>(
                t.A, t.B, t.m, t.E, cfg.eec_norm,
                t.label.empty() ? make_eec_label(t) : t.label));
        }
    }

    for (const SubjetMomentSpec& sm : cfg.subjet_moments) {
        for (const double n : sm.powers) {
            basis.push_back(std::make_unique<SubjetMomentBasisFunction>(
                sm.R, n, make_subjet_moment_label(sm.R, n)));
        }
    }

    return basis;
}

// =====================================================================
// Input-requirement query + evaluation
// =====================================================================
// Which inputs must the reader load for a given basis?
struct BasisInputs {
    bool          constituents = false;  // any EEC term present?
    std::set<int> subjet_rtags;          // stored subjet radii needed (rtags)
};

inline BasisInputs basis_inputs(const Basis& basis) {
    BasisInputs req;
    for (const std::unique_ptr<BasisFunction>& b : basis) {
        if (b->needs_constituents()) req.constituents = true;
        b->collect_radii(req.subjet_rtags);
    }
    return req;
}

// Evaluate every basis function for one jet -> length-nphys moment vector.
inline std::vector<double> evaluate_basis(const Basis& basis, const JetData& jet) {
    std::vector<double> g(basis.size(), 0.0);
    for (size_t j = 0; j < basis.size(); ++j) g[j] = basis[j]->evaluate(jet);
    return g;
}

// A compact machine-readable description of the basis (stored in metadata).
inline std::string basis_signature(const BasisConfig& cfg) {
    std::string s;
    if (cfg.use_eec) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "eec[%zu,norm=%.3g]+",
                      cfg.eec_terms.size(), cfg.eec_norm);
        s += buf;
    }
    for (const SubjetMomentSpec& sm : cfg.subjet_moments) {
        char buf[96];
        std::string pw;
        for (size_t i = 0; i < sm.powers.size(); ++i) {
            char pb[24];
            std::snprintf(pb, sizeof(pb), "%s%.3g",
                          i ? "," : "", sm.powers[i]);
            pw += pb;
        }
        std::snprintf(buf, sizeof(buf), "sjmom_R%.2f[%s]+", sm.R, pw.c_str());
        s += buf;
    }
    if (!s.empty() && s.back() == '+') s.pop_back();
    return s;
}
