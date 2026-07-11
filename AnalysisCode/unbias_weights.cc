// unbias_weights.cc
// =====================================================================
// Event-by-event reweighting to remove selection bias from heavy-ion
// jet samples, following the information-theoretic framework of
// Andres, Bossi, Holguin (arXiv:2501.17219 and "Unbiasing_by_reweighting").
//
// Builds:
//   c++ -std=c++17 -O2 unbias_weights.cc $(root-config --cflags --libs) -o unbias_weights
//
// =====================================================================
// VERSION 4 — numerical fixes (unchanged, see below)
// =====================================================================
//   FIX 1: dR_min screen on the EEC basis (ΔR^{−A} diverges as ΔR→0).
//   FIX 2: Auto learning rate = max_step / max(|g_scaled|).
//   FIX 3: Exponent clamping |λ·g| ≤ max_exp (default 10).
//   FIX 4: Enhanced debug diagnostics.
//
// =====================================================================
// VERSION 6 — fully configurable, multi-family basis
// =====================================================================
//
// The basis vectors used to construct the moments are now COMPLETELY
// CONFIGURABLE and may mix two families in arbitrary combination:
//
//   (1) Energy-correlator (EEC) basis functions
//         g = Σ_{i<k, ΔR<E} ΔR^{-A} [ln ΔR]^B z^m ,  z = pt_i pt_k / norm²
//       (the historical basis from the older unbias_weights.cc), and
//
//   (2) Subjet-pT moment basis functions
//         g = Σ_{subjets at radius R} pT_subjet^n
//       i.e. Mellin moments of the subjet-pT spectrum, for one or more
//       subjet radii R and one or more powers n.
//
// The definitions, the two families, and the builder live in the new
// self-contained header  basis_functions.h.  Adding a term is a one-line
// edit (or a command-line flag); adding a whole new *family* is a new
// BasisFunction subclass — the reader and optimiser below never change.
//
// Subjet input: ppjets_root.cc clusters C/A subjets at ntuple-production
// time over a radius scan R = 0.01..0.20 and stores them in branches
// "subjet_pt_R0pXX".  The subjet-moment functions CONSUME THOSE STORED
// SUBJETS DIRECTLY (the desired workflow).  A reclustering fallback (via
// subjet_basis.h) is retained only for older constituent-only ntuples,
// selectable per radius through --subjet-source.
//
// Key flags (see full list in main):
//   --basis-eec   {off,on}            enable the historical EEC family
//   --eec-norm    <val>               z normalisation (default 120)
//   --subjet-radii  "0.1,0.05"        subjet radii for the moment family
//   --subjet-powers "3,4,...,10"      powers n for the moment family
//   --subjet-R    <R>                 single-radius shortcut (default 0.1)
//   --subjet-source {auto,precomputed,recluster}   how to obtain subjets
//   --dR-min      <val>               EEC small-angle screen (default 1e-3)
//
// The default basis (no basis flags) reproduces the previous analysis:
// C/A R=0.1 subjet pT power sums, n = 3..10 — now sourced from the stored
// subjets rather than on-the-fly reclustering.
// =====================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TVector2.h>

#include "basis_functions.h"   // configurable multi-family basis framework
#include "subjet_basis.h"      // recluster_ca_subjet_pts (reclustering fallback)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#define MAXJ 100  // matches the fixed jet-array size written by ppjets_root.cc

static void die(const std::string &msg) {
    std::cerr << "error: " << msg << std::endl;
    std::exit(1);
}

static std::string get_arg(int argc, char *argv[],
                           const std::string &flag,
                           const std::string &def = "") {
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == flag && i + 1 < argc)
            return argv[i + 1];
    return def;
}

// Parse a comma-separated list of numbers, e.g. "3,4,5" or "0.1, 0.05".
static std::vector<double> parse_doubles(const std::string &s) {
    std::vector<double> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const size_t a = tok.find_first_not_of(" \t");
        if (a == std::string::npos) continue;
        const size_t b = tok.find_last_not_of(" \t");
        out.push_back(std::stod(tok.substr(a, b - a + 1)));
    }
    return out;
}

// =====================================================================
// Per-tree input plan
// =====================================================================
// For each stored subjet radius the basis needs, decide (per tree) whether
// to read the precomputed branch or to recluster from constituents, honouring
// --subjet-source.  Also records whether constituents must be loaded at all
// (either because an EEC term needs them, or because a radius is reclustered).
struct RadiusPlan {
    std::map<int, bool> precomputed;      // rtag -> true=precomputed / false=recluster
    bool                need_constituents = false;
};

static RadiusPlan plan_tree_inputs(TTree *t,
                                   const BasisInputs &req,
                                   const std::string &source_mode,
                                   const std::string &label) {
    RadiusPlan plan;
    plan.need_constituents = req.constituents;

    for (const int rtag : req.subjet_rtags) {
        const double R  = radius_from_tag(rtag);
        const std::string bn = subjet_branch_name("pt", R);
        const bool have = (t->GetBranch(bn.c_str()) != nullptr);

        bool precomp;
        if (source_mode == "precomputed") {
            if (!have)
                die(label + ": --subjet-source precomputed but branch '" + bn +
                    "' was not found");
            precomp = true;
        } else if (source_mode == "recluster") {
            precomp = false;
        } else {  // auto
            precomp = have;
        }
        plan.precomputed[rtag] = precomp;
        if (!precomp) plan.need_constituents = true;

        std::cout << "  [" << label << "] subjet R=" << std::setprecision(2)
                  << std::fixed << R << std::setprecision(6)
                  << " -> " << (precomp ? ("precomputed '" + bn + "'")
                                        : std::string("recluster from const_pt/eta/phi"))
                  << std::endl;
    }
    if (req.constituents)
        std::cout << "  [" << label << "] EEC terms read const_pt/eta/phi" << std::endl;

    if (plan.need_constituents &&
        !(t->GetBranch("const_pt") && t->GetBranch("const_eta") &&
          t->GetBranch("const_phi")))
        die(label + ": constituents required (EEC and/or reclustering) but "
            "const_pt/const_eta/const_phi branches are missing");

    return plan;
}

// =====================================================================
// Generic per-jet reader (shared by the target and source passes)
// =====================================================================
// Iterates every selected jet of `t`, assembles a JetData (constituents +
// stored/reclustered subjets for the required radii), and invokes
//   cb(const JetData& jd, double pt, double weight, int source)
// exactly once per selected jet.  Handles both "array" trees (nJets + arrays,
// e.g. tgenBefore from ppjets_root.cc) and "flat" one-jet-per-entry trees
// (e.g. the toy tX/tY/tRef/tBiased trees).  This single helper replaces the
// four near-identical read blocks the previous version carried.
template <class JetFn>
static void read_tree_jets(TTree *t,
                           const std::string &n_branch,
                           const std::string &pt_branch,
                           const std::string &wt_branch,
                           double pt_min, double pt_max,
                           const BasisInputs &req,
                           const RadiusPlan &plan,
                           double dR_min,
                           JetFn &&cb) {
    const bool array_style = (t->GetBranch(n_branch.c_str()) != nullptr);
    const bool has_source  = (t->GetBranch("source") != nullptr);
    const bool need_const  = plan.need_constituents;

    // Scratch storage for reclustered subjet lists (kept alive across the cb
    // call; JetData holds pointers into it).
    std::map<int, std::vector<double>> recl;

    // ---- assemble one jet's JetData and dispatch --------------------
    // `get_precomp(rtag)` returns the precomputed subjet-pT list for this jet
    // (or nullptr if unavailable); constituent lists are passed explicitly.
    auto emit = [&](double jetpt, double weight, int source,
                    const std::vector<double> *cpt,
                    const std::vector<double> *ceta,
                    const std::vector<double> *cphi,
                    const std::function<const std::vector<double> *(int)> &get_precomp) {
        JetData jd;
        jd.pt = jetpt;
        if (need_const) { jd.const_pt = cpt; jd.const_eta = ceta; jd.const_phi = cphi; }

        bool ok = true;
        for (const auto &kv : plan.precomputed) {
            const int  rtag    = kv.first;
            const bool precomp = kv.second;
            if (precomp) {
                const std::vector<double> *sp = get_precomp(rtag);
                if (!sp) { ok = false; break; }          // aligned entry missing
                jd.subjet_pt[rtag] = sp;
            } else {
                if (!cpt || !ceta || !cphi) { ok = false; break; }
                recl[rtag] = recluster_ca_subjet_pts(*cpt, *ceta, *cphi,
                                                     radius_from_tag(rtag));
                jd.subjet_pt[rtag] = &recl[rtag];
            }
        }
        if (!ok) return;

        if (req.constituents) build_constituent_pairs(jd, dR_min);
        cb(static_cast<const JetData &>(jd), jetpt, weight, source);
    };

    const Long64_t ne = t->GetEntries();

    if (array_style) {
        int   nJ = 0;
        float pt[MAXJ] = {0};
        float wt = 1.0f;
        int   src = -1;
        std::vector<std::vector<double>> *cp = nullptr, *ce = nullptr, *cf = nullptr;
        std::map<int, std::vector<std::vector<double>> *> sjb;  // precomputed subjet branches

        t->SetBranchAddress(n_branch.c_str(), &nJ);
        t->SetBranchAddress(pt_branch.c_str(), pt);
        t->SetBranchAddress(wt_branch.c_str(), &wt);
        if (has_source) t->SetBranchAddress("source", &src);
        if (need_const) {
            t->SetBranchAddress("const_pt",  &cp);
            t->SetBranchAddress("const_eta", &ce);
            t->SetBranchAddress("const_phi", &cf);
        }
        for (const auto &kv : plan.precomputed)
            if (kv.second) {
                sjb[kv.first] = nullptr;
                t->SetBranchAddress(subjet_branch_name("pt", radius_from_tag(kv.first)).c_str(),
                                    &sjb[kv.first]);
            }

        for (Long64_t ev = 0; ev < ne; ++ev) {
            t->GetEntry(ev);
            const int nJc = std::min(nJ, MAXJ);
            for (int j = 0; j < nJc; ++j) {
                const double x = pt[j];
                if (x < pt_min || x > pt_max) continue;

                const std::vector<double> *cpt = nullptr, *ceta = nullptr, *cphi = nullptr;
                if (need_const) {
                    if (!cp || !ce || !cf) die("missing constituent branches");
                    if (j >= (int)cp->size() || j >= (int)ce->size() || j >= (int)cf->size())
                        continue;
                    cpt = &cp->at(j); ceta = &ce->at(j); cphi = &cf->at(j);
                }
                emit(x, (double)wt, src, cpt, ceta, cphi,
                     [&](int rtag) -> const std::vector<double> * {
                         auto it = sjb.find(rtag);
                         if (it == sjb.end() || !it->second ||
                             j >= (int)it->second->size())
                             return nullptr;
                         return &it->second->at(j);
                     });
            }
        }
    } else {
        float pv = 0.0f;
        float wt = 1.0f;
        int   src = -1;
        std::vector<double> *cp = nullptr, *ce = nullptr, *cf = nullptr;
        std::map<int, std::vector<double> *> sjb;

        t->SetBranchAddress(pt_branch.c_str(), &pv);
        t->SetBranchAddress(wt_branch.c_str(), &wt);
        if (has_source) t->SetBranchAddress("source", &src);
        if (need_const) {
            t->SetBranchAddress("const_pt",  &cp);
            t->SetBranchAddress("const_eta", &ce);
            t->SetBranchAddress("const_phi", &cf);
        }
        for (const auto &kv : plan.precomputed)
            if (kv.second) {
                sjb[kv.first] = nullptr;
                t->SetBranchAddress(subjet_branch_name("pt", radius_from_tag(kv.first)).c_str(),
                                    &sjb[kv.first]);
            }

        for (Long64_t ev = 0; ev < ne; ++ev) {
            t->GetEntry(ev);
            const double x = pv;
            if (x < pt_min || x > pt_max) continue;
            emit(x, (double)wt, src, cp, ce, cf,
                 [&](int rtag) -> const std::vector<double> * {
                     auto it = sjb.find(rtag);
                     return (it == sjb.end()) ? nullptr : it->second;
                 });
        }
    }
}

// =====================================================================
// Print dot-product distribution diagnostics
// =====================================================================
static void print_dot_diagnostics(
        const std::vector<double> &lambda,
        const std::vector<std::vector<double>> &gvals,
        int nphys, size_t N)
{
    std::vector<double> dots(N);
    for (size_t i = 0; i < N; ++i) {
        double dot = 0.0;
        for (int j = 0; j < nphys; ++j)
            dot += lambda[j] * gvals[i][j];
        dots[i] = dot;
    }
    std::sort(dots.begin(), dots.end());
    auto pct = [&](double p) -> double {
        size_t idx = std::min((size_t)(p * N), N - 1);
        return dots[idx];
    };
    std::cout << "  dot lam*g distribution:"
              << "  min=" << dots.front()
              << "  p1=" << pct(0.01)
              << "  p10=" << pct(0.10)
              << "  p50=" << pct(0.50)
              << "  p90=" << pct(0.90)
              << "  p99=" << pct(0.99)
              << "  max=" << dots.back()
              << std::endl;
}

// =====================================================================
// Print per-component gradient diagnostics
// =====================================================================
static void print_gradient_diagnostics(
        const std::vector<double> &grad,
        const std::vector<double> &sj,
        const std::vector<double> &d,
        const std::vector<double> &c,
        const std::vector<double> &basis_scale,
        int nphys)
{
    // Physical gradient, and decompose what drives each component.
    std::cout << "  gradient (physical dL/dlam_real):" << std::endl;
    double gnorm2 = 0.0;
    for (int k = 0; k < nphys; ++k) {
        const double gp = grad[k] * basis_scale[k];
        gnorm2 += gp * gp;
        std::cout << "    grad[" << k << "] = " << std::scientific << gp;
        // Show the residual driver s_j for diagonal term.
        std::cout << "   (s[" << k << "]=" << sj[k] << ")";
        std::cout << std::fixed << std::endl;
    }
    std::cout << "  |grad| = " << std::scientific << std::sqrt(gnorm2)
              << std::fixed << std::endl;
}

// =====================================================================
// Print the covariance matrix conditioning
// =====================================================================
static void print_covariance_diagnostics(
        const std::vector<std::vector<double>> &GG,
        const std::vector<double> &Sj_unc,
        const std::vector<double> &d,
        double S, double S_unc, int nphys)
{
    std::cout << "  Covariance diag (Cov[j,j]) and off-diag correlations:"
              << std::endl;
    // Compute full covariance matrix (using unclamped stats).
    std::vector<std::vector<double>> Cov(nphys, std::vector<double>(nphys, 0.0));
    std::vector<double> d_unc(nphys, 0.0);
    if (S_unc > 0) {
        for (int j = 0; j < nphys; ++j)
            d_unc[j] = Sj_unc[j] / S_unc;
    }
    for (int j = 0; j < nphys; ++j)
        for (int k = j; k < nphys; ++k) {
            const double gg = GG[j][k];
            // Cov using the full d (weighted avg over all jets) and
            // GG from unclamped jets.
            Cov[j][k] = gg / S - d[j] * d[k];
            Cov[k][j] = Cov[j][k];
        }

    // Print diagonal (variance).
    for (int j = 0; j < nphys; ++j)
        std::cout << "    Cov[" << j << "," << j << "] = "
                  << std::scientific << Cov[j][j] << std::fixed << std::endl;

    // Print correlation matrix.
    std::cout << "  Correlation matrix:" << std::endl;
    for (int j = 0; j < nphys; ++j) {
        std::cout << "    [" << j << "]:";
        for (int k = 0; k < nphys; ++k) {
            double denom = std::sqrt(std::fabs(Cov[j][j] * Cov[k][k]));
            double corr = (denom > 1e-30) ? Cov[j][k] / denom : 0.0;
            std::cout << " " << std::setw(7) << std::setprecision(3)
                      << std::fixed << corr;
        }
        std::cout << std::setprecision(6) << std::endl;
    }
}

// =====================================================================
// MAIN
// =====================================================================

int main(int argc, char *argv[]) {

    // ---- Command-line arguments ----
    const std::string input        = get_arg(argc, argv, "--input");
    const std::string tree_name    = get_arg(argc, argv, "--tree",          "tgenBefore");
    const std::string tgt_tree_n   = get_arg(argc, argv, "--target-tree",   tree_name);
    const std::string n_branch     = get_arg(argc, argv, "--n-branch",      "nJets");
    const std::string pt_branch    = get_arg(argc, argv, "--pt-branch",     "pt");
    const std::string wt_branch    = get_arg(argc, argv, "--weight-branch", "weight");
    const std::string target_input = get_arg(argc, argv, "--target-input");
    const std::string out_name     = get_arg(argc, argv, "--out",           "unbias_weights.root");
    const std::string run_mode     = get_arg(argc, argv, "--mode",          "run");
    const std::string scale_mode   = get_arg(argc, argv, "--scale-basis",   "target");
    const std::string lr_spec      = get_arg(argc, argv, "--adam-lr",       "auto");
    const std::string subjet_source = get_arg(argc, argv, "--subjet-source", "auto");

    // ---- Basis configuration flags (VERSION 6) ----
    const std::string basis_eec_s  = get_arg(argc, argv, "--basis-eec",     "off");
    const std::string subjet_radii_s  = get_arg(argc, argv, "--subjet-radii",  "");
    const std::string subjet_powers_s = get_arg(argc, argv, "--subjet-powers", "");

    const double pt_min  = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
    const double pt_max  = std::stod(get_arg(argc, argv, "--pt-max", "1e9"));
    const double dR_min  = std::stod(get_arg(argc, argv, "--dR-min", "0.001"));
    const double subjet_R = std::stod(get_arg(argc, argv, "--subjet-R", "0.1"));
    const double eec_norm = std::stod(get_arg(argc, argv, "--eec-norm", "120.0"));

    const int    max_iter        = std::stoi(get_arg(argc, argv, "--max-iter",       "200000"));
    const double loss_tol        = std::stod(get_arg(argc, argv, "--loss-tol",       "1e-10"));
    const double tol             = std::stod(get_arg(argc, argv, "--tol",            "1e-6"));
    const double adam_beta1      = std::stod(get_arg(argc, argv, "--adam-beta1",     "0.9"));
    const double adam_beta2      = std::stod(get_arg(argc, argv, "--adam-beta2",     "0.999"));
    const double adam_eps        = std::stod(get_arg(argc, argv, "--adam-eps",       "1e-8"));
    const int    print_interval  = std::stoi(get_arg(argc, argv, "--print-interval", "1000"));
    const int    debug_steps     = std::stoi(get_arg(argc, argv, "--debug-steps",    "50"));
    const double max_exp         = std::stod(get_arg(argc, argv, "--max-exp",        "10.0"));
    const double max_step_target = std::stod(get_arg(argc, argv, "--max-step-target","0.01"));
    const int    patience        = std::stoi(get_arg(argc, argv, "--patience",       "5000"));
    const double lr_decay        = std::stod(get_arg(argc, argv, "--lr-decay",       "0.5"));
    const double lr_min_val      = std::stod(get_arg(argc, argv, "--lr-min",         "1e-12"));

    if (input.empty())        die("--input is required");
    if (target_input.empty()) die("--target-input is required");
    if (run_mode != "run" && run_mode != "debug") die("--mode must be 'run' or 'debug'");
    if (subjet_source != "auto" && subjet_source != "precomputed" && subjet_source != "recluster")
        die("--subjet-source must be 'auto', 'precomputed', or 'recluster'");

    // =================================================================
    // Build the (configurable) basis.
    // -----------------------------------------------------------------
    // Start from the header defaults (BasisConfig / build_basis in
    // basis_functions.h), then apply the command-line overrides below.
    // To hard-code a bespoke basis (e.g. per-radius power lists, custom
    // EEC terms), edit BasisConfig in basis_functions.h and comment out
    // the CLI-override block marked >>> below.
    // =================================================================
    BasisConfig cfg;
    cfg.use_eec    = (basis_eec_s == "on" || basis_eec_s == "1" || basis_eec_s == "true");
    cfg.eec_norm   = eec_norm;
    cfg.eec_dR_min = dR_min;

    // >>> CLI override of the subjet-moment family (radii × powers) <<<
    {
        std::vector<double> radii = subjet_radii_s.empty()
            ? std::vector<double>{ subjet_R }
            : parse_doubles(subjet_radii_s);
        std::vector<double> powers = subjet_powers_s.empty()
            ? std::vector<double>{ 3, 4, 5, 6, 7, 8, 9, 10 }
            : parse_doubles(subjet_powers_s);
        cfg.subjet_moments.clear();
        for (const double R : radii) {
            if (R <= 0.0) die("--subjet-radii entries must be positive");
            cfg.subjet_moments.push_back({ R, powers });
        }
    }
    // >>> end CLI override <<<

    Basis basis = build_basis(cfg);
    const int nphys = (int)basis.size();
    if (nphys == 0) die("empty basis: enable EEC and/or subjet moments");

    const BasisInputs req = basis_inputs(basis);

    // Primary subjet radius (smallest requested) — used for the metadata
    // and for the optional per-jet subjet_pt output branch.
    const bool   have_subjets = !req.subjet_rtags.empty();
    const double primary_R    = have_subjets ? radius_from_tag(*req.subjet_rtags.begin()) : 0.0;

    std::cout << "=== Basis functions (" << nphys << ") ===" << std::endl;
    for (int j = 0; j < nphys; ++j)
        std::cout << "  [" << j << "] " << basis[j]->label() << std::endl;
    std::cout << "  [N] normalisation (analytic)" << std::endl;
    std::cout << "  signature = " << basis_signature(cfg) << std::endl;
    std::cout << "  EEC family: " << (cfg.use_eec ? "on" : "off")
              << " (norm=" << cfg.eec_norm << ", dR_min=" << cfg.eec_dR_min << ")" << std::endl;
    std::cout << "  --subjet-source = " << subjet_source << std::endl << std::endl;

    // ---- Read TARGET sample ----
    std::vector<double> c(nphys, 0.0);
    double tgt_sumW = 0.0;
    size_t tgt_njets = 0;
    {
        TFile tf(target_input.c_str(), "READ");
        if (tf.IsZombie()) die("failed to open target file: " + target_input);
        TTree *tt = dynamic_cast<TTree*>(tf.Get(tgt_tree_n.c_str()));
        if (!tt) die("failed to find target tree: " + tgt_tree_n);

        const RadiusPlan tplan = plan_tree_inputs(tt, req, subjet_source, "target");
        read_tree_jets(tt, n_branch, pt_branch, wt_branch, pt_min, pt_max,
                       req, tplan, dR_min,
                       [&](const JetData &jd, double /*x*/, double w, int /*src*/) {
                           const std::vector<double> gv = evaluate_basis(basis, jd);
                           tgt_sumW += w; ++tgt_njets;
                           for (int b = 0; b < nphys; ++b) c[b] += w * gv[b];
                       });
        tf.Close();
    }
    if (tgt_sumW <= 0.0) die("target has non-positive weight sum");
    for (int j = 0; j < nphys; ++j) c[j] /= tgt_sumW;

    std::cout << "Target: " << tgt_njets << " jets, sum_w = " << tgt_sumW << std::endl;
    std::cout << "Target moments c[j] (weighted avg, before scaling):" << std::endl;
    for (int j = 0; j < nphys; ++j)
        std::cout << "  c[" << j << "] = " << c[j] << std::endl;
    std::cout << std::endl;

    // ---- Read SOURCE sample ----
    TFile in_file(input.c_str(), "READ");
    if (in_file.IsZombie()) die("failed to open: " + input);
    TTree *tree = dynamic_cast<TTree*>(in_file.Get(tree_name.c_str()));
    if (!tree) die("failed to find tree: " + tree_name);

    double wX_an = 0, wYp_an = 0; bool has_an = false;
    if (TTree *tW = dynamic_cast<TTree*>(in_file.Get("tWeights"))) {
        tW->SetBranchAddress("wX", &wX_an);
        tW->SetBranchAddress("wYprime", &wYp_an);
        if (tW->GetEntries()>0) { tW->GetEntry(0); has_an=true; }
    }

    std::vector<double> pts, base_w;
    std::vector<std::vector<double>> gvals;
    std::vector<int> sources;
    std::vector<std::vector<double>> ocp, oce, ocf;  // constituents (for output, if read)
    std::vector<std::vector<double>> osj;            // primary-radius subjet pT's (for output)

    pts.reserve((size_t)tree->GetEntries()*2); base_w.reserve((size_t)tree->GetEntries()*2);

    const RadiusPlan splan = plan_tree_inputs(tree, req, subjet_source, "source");
    read_tree_jets(tree, n_branch, pt_branch, wt_branch, pt_min, pt_max,
                   req, splan, dR_min,
                   [&](const JetData &jd, double x, double w, int src) {
                       pts.push_back(x); base_w.push_back(w); sources.push_back(src);
                       gvals.push_back(evaluate_basis(basis, jd));
                       if (jd.const_pt && jd.const_eta && jd.const_phi) {
                           ocp.push_back(*jd.const_pt);
                           oce.push_back(*jd.const_eta);
                           ocf.push_back(*jd.const_phi);
                       }
                       if (have_subjets) {
                           const std::vector<double> *sp = jd.subjets_at(primary_R);
                           if (sp) osj.push_back(*sp);
                       }
                   });

    const size_t N = pts.size();
    if (N == 0) die("no jets passed selection");
    const double Nd = (double)N;
    std::cout << "Source: " << N << " jets loaded." << std::endl;

    // ---- Basis scaling ----
    std::vector<double> bscale(nphys, 1.0);
    if (scale_mode == "target") {
        for (int j=0; j<nphys; ++j) {
            double s = std::fabs(c[j]);
            bscale[j] = (s>1e-30)?s:1.0;
            c[j] /= bscale[j];
        }
        for (size_t i=0; i<N; ++i)
            for (int j=0; j<nphys; ++j)
                gvals[i][j] /= bscale[j];
        std::cout << "Basis scaling:" << std::endl;
        for (int j=0; j<nphys; ++j)
            std::cout << "  s[" << j << "]=" << std::setw(12) << bscale[j]
                      << "  c_sc[" << j << "]=" << c[j] << std::endl;
        std::cout << std::endl;
    }

    // ---- Basis statistics & auto lr ----
    double max_abs_g = 0.0;
    {
        std::cout << "Basis stats (scaled):" << std::endl;
        for (int j=0; j<nphys; ++j) {
            double mn=1e30, mx=-1e30, su=0, su2=0;
            int nzero = 0;
            for (size_t i=0; i<N; ++i) {
                double v = gvals[i][j];
                mn=std::min(mn,v); mx=std::max(mx,v);
                su+=v; su2+=v*v;
                if (v == 0.0) ++nzero;
            }
            double mean = su/Nd;
            double var = su2/Nd - mean*mean;
            double amx = std::max(std::fabs(mn), std::fabs(mx));
            max_abs_g = std::max(max_abs_g, amx);
            std::cout << "  g[" << j << "]:"
                      << "  min=" << std::setw(10) << mn
                      << "  max=" << std::setw(10) << mx
                      << "  mean=" << std::setw(10) << mean
                      << "  std=" << std::setw(10) << std::sqrt(std::max(0.,var))
                      << "  nzero=" << nzero
                      << std::endl;
        }
        std::cout << "  max|g| = " << max_abs_g << std::endl << std::endl;
    }

    double adam_lr;
    if (lr_spec == "auto") {
        adam_lr = max_step_target / std::max(1.0, max_abs_g);
        std::cout << "Auto lr = " << max_step_target << " / " << max_abs_g
                  << " = " << std::scientific << adam_lr << std::fixed << std::endl;
    } else {
        adam_lr = std::stod(lr_spec);
        std::cout << "Manual lr = " << std::scientific << adam_lr << std::fixed << std::endl;
        // Safety warning.
        double worst = adam_lr * max_abs_g;
        if (worst > 1.0)
            std::cout << "  WARNING: lr * max|g| = " << worst
                      << " >> 1.  First step will change max exponent by "
                      << worst << ".  Consider --adam-lr auto." << std::endl;
    }
    std::cout << "max_exp=" << max_exp
              << "  patience=" << patience
              << "  lr_decay=" << lr_decay << std::endl << std::endl;

    // ================================================================
    // ADAM OPTIMISATION
    // ================================================================
    const bool dbg = (run_mode == "debug");
    const int  MI  = dbg ? debug_steps : max_iter;
    const int  PI  = dbg ? 1 : print_interval;

    std::vector<double> lam(nphys, 0.0);
    std::vector<double> m_a(nphys, 0.0), v_a(nphys, 0.0);
    double lam_norm = 0.0;

    double best_loss = 1e30;
    std::vector<double> best_lam(nphys, 0.0);
    double best_lnorm = 0.0;
    int best_it = 0;
    double cur_lr = adam_lr;

    std::cout << "Starting Adam (" << nphys << " params + analytic norm)" << std::endl;
    std::cout << "  lr=" << std::scientific << cur_lr << std::fixed
              << "  max_iter=" << MI << std::endl << std::endl;

    for (int it = 1; it <= MI; ++it) {

        // ---- Reweighted statistics ----
        double S = 0, sumW2 = 0, S_unc = 0;
        int ncl = 0;
        std::vector<double> Sj(nphys, 0.0), Sj_unc(nphys, 0.0);
        std::vector<std::vector<double>> GG(nphys, std::vector<double>(nphys, 0.0));

        for (size_t i = 0; i < N; ++i) {
            double dot = 0;
            for (int j=0; j<nphys; ++j) dot += lam[j] * gvals[i][j];
            bool cl = (dot > max_exp || dot < -max_exp);
            double dc = std::max(-max_exp, std::min(max_exp, dot));
            if (cl) ++ncl;
            double Wi = base_w[i] * std::exp(-dc);
            S += Wi; sumW2 += Wi*Wi;
            for (int j=0; j<nphys; ++j) Sj[j] += Wi * gvals[i][j];
            if (!cl) {
                S_unc += Wi;
                for (int j=0; j<nphys; ++j) {
                    double Wg = Wi * gvals[i][j];
                    Sj_unc[j] += Wg;
                    for (int k=j; k<nphys; ++k) GG[j][k] += Wg * gvals[i][k];
                }
            }
        }
        if (S <= 0) die("S=0");

        // ---- Moments, normalisation ----
        std::vector<double> d(nphys);
        for (int j=0; j<nphys; ++j) d[j] = Sj[j] / S;
        lam_norm = std::log(S / Nd);
        double Neff = (sumW2>0) ? S*S/sumW2 : 0;

        // ---- Loss ----
        std::vector<double> r(nphys);
        double loss = 0, mx_rel = 0;
        for (int j=0; j<nphys; ++j) {
            double dn = c[j]+d[j];
            r[j] = (std::fabs(dn)>1e-30) ? (c[j]-d[j])/dn : 0;
            loss += r[j]*r[j];
            double rl = std::fabs(c[j]-d[j])/(std::fabs(c[j])+1e-30);
            if (rl > mx_rel) mx_rel = rl;
        }

        if (loss < best_loss) {
            best_loss = loss; best_lam = lam; best_lnorm = lam_norm; best_it = it;
        }

        // ---- Plateau lr decay ----
        if (!dbg && (it - best_it) >= patience && cur_lr > lr_min_val) {
            cur_lr = std::max(lr_min_val, cur_lr * lr_decay);
            best_it = it;
            std::cout << "  [lr decay] iter " << it << " -> lr="
                      << std::scientific << cur_lr << std::fixed << std::endl;
        }

        // ---- Print ----
        bool last = (it == MI);
        bool prn  = (PI>0 && it%PI==0) || it==1 || last;
        if (prn) {
            std::cout << "iter " << std::setw(6) << it
                      << "  loss=" << std::scientific << std::setprecision(6) << loss
                      << "  mx_rel=" << mx_rel
                      << "  Neff=" << std::fixed << std::setprecision(0) << Neff
                      << "/" << N
                      << "  lnorm=" << std::setprecision(6) << lam_norm
                      << "  cl=" << ncl
                      << "  lr=" << std::scientific << cur_lr << std::fixed
                      << "  best=" << std::scientific << best_loss << std::fixed
                      << std::endl;
            for (int j=0; j<nphys; ++j) {
                std::cout << "    [" << j << "]"
                          << "  d=" << std::scientific << std::setw(13) << d[j]*bscale[j]
                          << "  c=" << std::setw(13) << c[j]*bscale[j]
                          << "  r=" << std::setw(11) << r[j]
                          << "  lam=" << std::setw(11) << lam[j]/bscale[j]
                          << std::fixed << std::endl;
            }
        }

        // ---- Debug extras ----
        if (dbg && prn) {
            print_dot_diagnostics(lam, gvals, nphys, N);

            // Gradient scalar weights s_j (computed below for gradient anyway).
            std::vector<double> sj(nphys, 0.0);
            for (int j=0; j<nphys; ++j) {
                double dn = c[j]+d[j];
                if (std::fabs(dn)>1e-30) sj[j] = c[j]*r[j]/(dn*dn);
            }

            // Print covariance diagnostics every 10 iters in debug.
            if (it == 1 || it % 10 == 0 || last)
                print_covariance_diagnostics(GG, Sj_unc, d, S, S_unc, nphys);
        }

        // ---- Convergence ----
        if (!dbg && (loss < loss_tol || mx_rel < tol)) {
            std::cout << "\n*** Converged at iter " << it
                      << "  loss=" << std::scientific << loss
                      << "  mx_rel=" << mx_rel << std::fixed << std::endl;
            break;
        }

        // ---- Gradient ----
        std::vector<double> sj(nphys, 0.0);
        for (int j=0; j<nphys; ++j) {
            double dn = c[j]+d[j];
            if (std::fabs(dn)>1e-30) sj[j] = c[j]*r[j]/(dn*dn);
        }
        std::vector<double> grad(nphys, 0.0);
        for (int k=0; k<nphys; ++k) {
            double acc = 0;
            for (int j=0; j<nphys; ++j) {
                double gg = (j<=k) ? GG[j][k] : GG[k][j];
                acc += sj[j] * (gg - d[j] * Sj_unc[k]);
            }
            grad[k] = (4.0 / S) * acc;
        }

        if (dbg && prn)
            print_gradient_diagnostics(grad, sj, d, c, bscale, nphys);

        // ---- Adam update ----
        double bc1 = 1.0 - std::pow(adam_beta1, (double)it);
        double bc2 = 1.0 - std::pow(adam_beta2, (double)it);
        for (int k=0; k<nphys; ++k) {
            m_a[k] = adam_beta1*m_a[k] + (1-adam_beta1)*grad[k];
            v_a[k] = adam_beta2*v_a[k] + (1-adam_beta2)*grad[k]*grad[k];
            double mh = m_a[k]/bc1, vh = v_a[k]/bc2;
            lam[k] -= cur_lr * mh / (std::sqrt(vh) + adam_eps);
        }
    }

    lam = best_lam; lam_norm = best_lnorm;
    std::cout << "\nBest at iter " << best_it
              << " (loss=" << std::scientific << best_loss << std::fixed << ")" << std::endl;

    // ---- Un-scale ----
    if (scale_mode == "target") {
        std::cout << "\nUn-scaling:" << std::endl;
        for (int j=0; j<nphys; ++j) {
            lam[j] /= bscale[j];
            std::cout << "  lam[" << j << "]=" << std::setw(14) << lam[j]
                      << "  (s=" << bscale[j] << ")" << std::endl;
        }
        for (size_t i=0; i<N; ++i)
            for (int j=0; j<nphys; ++j) gvals[i][j] *= bscale[j];
    }
    std::cout << "  lam_norm=" << lam_norm << std::endl;

    // ---- Write output ----
    TFile of(out_name.c_str(), "RECREATE");
    TTree tw("tweights", "per-jet unbiasing weights");
    float  op=0; double wu=1,wb=1,wt2=1,wa=1; int os=-1;
    std::vector<double> xcp,xce,xcf,xsj;
    bool sco = (ocp.size()==N);
    bool ssj = (osj.size()==N);
    tw.Branch("pt",&op,"pt/F");
    tw.Branch("w_unbias",&wu,"w_unbias/D");
    tw.Branch("w_base",&wb,"w_base/D");
    tw.Branch("w_total",&wt2,"w_total/D");
    tw.Branch("source",&os,"source/I");
    if (has_an) tw.Branch("w_analytic",&wa,"w_analytic/D");
    if (sco) { tw.Branch("const_pt",&xcp); tw.Branch("const_eta",&xce); tw.Branch("const_phi",&xcf); }
    if (ssj) tw.Branch("subjet_pt",&xsj);

    for (size_t i=0; i<N; ++i) {
        double dot=0;
        for (int j=0; j<nphys; ++j) dot += lam[j]*gvals[i][j];
        double dc = std::max(-max_exp, std::min(max_exp, dot));
        wu = std::exp(-dc - lam_norm);
        wb = base_w[i]; wt2 = wb*wu; op = (float)pts[i];
        os = (i<sources.size()) ? sources[i] : -1;
        if (has_an) wa = (os==0) ? wX_an : wYp_an;
        if (sco) { xcp=ocp[i]; xce=oce[i]; xcf=ocf[i]; }
        if (ssj) xsj=osj[i];
        tw.Fill();
    }

    TTree tm("meta","fit metadata");
    auto la = lam; la.push_back(lam_norm);
    std::string bl = basis_signature(cfg) + "+norm_analytic";
    std::string subjet_src_label = subjet_source;
    tm.Branch("lambda",&la); tm.Branch("lambda_norm",&lam_norm,"lambda_norm/D");
    tm.Branch("basis",&bl);
    tm.Branch("subjet_source",&subjet_src_label);
    tm.Branch("pt_min",const_cast<double*>(&pt_min),"pt_min/D");
    tm.Branch("pt_max",const_cast<double*>(&pt_max),"pt_max/D");
    tm.Branch("best_loss",const_cast<double*>(&best_loss),"best_loss/D");
    double sr_out = primary_R;
    tm.Branch("subjet_R",&sr_out,"subjet_R/D");
    tm.Fill();
    tw.Write(); tm.Write(); of.Close(); in_file.Close();

    // ---- Verification ----
    std::cout << "\nFinal lambdas:" << std::endl;
    for (int j=0; j<nphys; ++j)
        std::cout << "  [" << j << "] " << std::setw(14) << lam[j]
                  << "  (" << basis[j]->label() << ")" << std::endl;
    std::cout << "  lam_norm = " << lam_norm << std::endl;

    double cS=0,cS2=0; std::vector<double> cd(nphys,0);
    for (size_t i=0; i<N; ++i) {
        double dot=0;
        for (int j=0; j<nphys; ++j) dot += lam[j]*gvals[i][j];
        double dc = std::max(-max_exp, std::min(max_exp, dot));
        double w = base_w[i]*std::exp(-dc-lam_norm);
        cS+=w; cS2+=w*w;
        for (int j=0; j<nphys; ++j) cd[j]+=w*gvals[i][j];
    }
    for (int j=0; j<nphys; ++j) cd[j]/=cS;

    std::cout << "\nVerification:" << std::endl;
    std::cout << "  avg(w)=" << cS/Nd << " (target 1.0)" << std::endl;
    std::cout << "  Neff=" << std::setprecision(0) << cS*cS/cS2
              << "/" << N << std::setprecision(6) << std::endl;
    std::cout << "  Moments:" << std::endl;
    for (int j=0; j<nphys; ++j) {
        double cj = c[j]*bscale[j];
        std::cout << "    [" << j << "] d=" << std::setw(12) << cd[j]
                  << " c=" << std::setw(12) << cj
                  << " ratio=" << cd[j]/cj << std::endl;
    }
    std::cout << "\nWrote " << out_name << std::endl;
    return 0;
}
