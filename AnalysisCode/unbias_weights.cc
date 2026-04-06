// unbias_weights.cc
// =====================================================================
// Event-by-event reweighting to remove selection bias from heavy-ion
// jet samples, following the information-theoretic framework of
// Andres, Bossi, Holguin (arXiv:2501.17219 and "Unbiasing_by_reweighting").
//
// Builds:
//   c++ -std=c++17 -O2 unbias_weights.cc \
//       $(root-config --cflags --libs) -o unbias_weights
//
// =====================================================================
// VERSION 4 — fixes and improvements
// =====================================================================
//
// FIX 1: dR_min = 0.005 cut on all basis functions.
//   The basis functions involve ΔR^{−A} which diverges as ΔR → 0.
//   A few jets with nearly-collinear constituent pairs produce extreme
//   values (g_scaled up to 100,000× the mean), which makes the
//   exponential weight form exp(−λ·g) catastrophically unstable.
//   Requiring dR > 0.005 screens these divergences.
//
// FIX 2: Auto learning rate (default) = max_step / max(|g_scaled|).
//   Ensures each Adam step changes the worst-case exponent by at most
//   max_step (default 0.01).  This is the ONLY safe way to set lr
//   for this problem because Adam normalizes steps to ±lr regardless
//   of gradient magnitude.
//
// FIX 3: Exponent clamping |λ·g| ≤ max_exp (default 10).
//   Safety net — no jet weight exceeds exp(10) ≈ 22000× base.
//
// FIX 4: Enhanced debug diagnostics (see debug mode output).
//
// =====================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

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

// =====================================================================
// BASIS FUNCTION DEFINITIONS
// =====================================================================

struct BasisFuncDef {
    double A;               // ΔR^{−A}
    double B;               // [ln(ΔR)]^B
    int    m;               // (p_T^i p_T^k / p_T²)^m
    double E;               // angular upper cut: Θ(ΔR < E)
    std::string label;
};

static std::vector<BasisFuncDef> get_default_basis() {
    return {
        {1.0, 4.0, 1, 0.4, "g1: dR^{-3/2} ln^4 z^1  [dR<0.1]"},   // eq 4.15
        {1.0, 3.0, 1, 0.4, "g2: dR^{-1} z^2  [dR<0.2]"},           // eq 4.16
        {1.0, 2.0, 1, 0.4, "g3: dR^{-3/2} z^2  [dR<0.2]"},         // eq 4.17
        {1.0, 1.0, 1, 0.4, "g4: dR^{-1} ln^4 z^2  [dR<0.2]"},      // eq 4.18
        {1.0, 0.0, 1, 0.4, "g5: dR^{-3/2} ln^4 z^2  [dR<0.2]"},    // eq 4.19
        
        {0.0, 4.0, 1, 0.4, "g6: "},   // eq 4.15
        {0.0, 3.0, 1, 0.4, "g7: "},           // eq 4.16
        {0.0, 2.0, 1, 0.4, "g8:"},         // eq 4.17
        {0.0, 1.0, 1, 0.4, "g9: ]"},      // eq 4.18
        {0.0, 0.0, 1, 0.4, "g10"},    // eq 4.19
        
        {-1.0, 4.0, 1, 0.4, "g11: "},   // eq 4.15
        {-1.0, 3.0, 1, 0.4, "g12: "},           // eq 4.16
        {-1.0, 2.0, 1, 0.4, "g13:"},         // eq 4.17
        {-1.0, 1.0, 1, 0.4, "g14: ]"},      // eq 4.18
        {-1.0, 0.0, 1, 0.4, "g15"},    // eq 4.19
    };
}

// =====================================================================
// BASIS FUNCTION EVALUATION
// =====================================================================

/// Evaluate all basis functions for one jet.
/// dR_min screens the small-angle divergence in ΔR^{−A}.
static std::vector<double> evaluate_basis(
        const std::vector<BasisFuncDef> &basis,
        double ptjet,
        const std::vector<double> &cpt,
        const std::vector<double> &ceta,
        const std::vector<double> &cphi,
        double dR_min)
{
    const int nb = (int)basis.size();
    std::vector<double> g(nb, 0.0);
    if (ptjet <= 0.0) return g;

    const double pt2 = ptjet * ptjet;
    const size_t nc  = cpt.size();

    for (size_t i = 0; i < nc; ++i) {
        if (cpt[i] <= 0.0) continue;
        for (size_t k = i + 1; k < nc; ++k) {
            if (cpt[k] <= 0.0) continue;

            const double dphi = TVector2::Phi_mpi_pi(cphi[i] - cphi[k]);
            const double deta = ceta[i] - ceta[k];
            const double dR   = std::sqrt(deta * deta + dphi * dphi);
            if (dR <= dR_min) continue;   // screen small-angle divergence

            const double lndR = std::log(dR);
            const double z    = (cpt[i] * cpt[k]) / pt2;

            for (int j = 0; j < nb; ++j) {
                const auto &bf = basis[j];
                if (dR >= bf.E) continue;
                double val = std::pow(dR, -bf.A);
                if (bf.B > 0.0) val *= std::pow(lndR, bf.B);
                val *= std::pow(z, bf.m);
                g[j] += val;
            }
        }
    }
    return g;
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

    const double pt_min  = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
    const double pt_max  = std::stod(get_arg(argc, argv, "--pt-max", "1e9"));
    const double dR_min  = std::stod(get_arg(argc, argv, "--dR-min", "0.001"));

    const int    max_iter        = std::stoi(get_arg(argc, argv, "--max-iter",       "100000"));
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

    // ---- Basis functions ----
    std::vector<BasisFuncDef> basis = get_default_basis();
    const int nphys = (int)basis.size();

    std::cout << "=== Physics basis functions (" << nphys << ") ===" << std::endl;
    for (int j = 0; j < nphys; ++j)
        std::cout << "  [" << j << "] " << basis[j].label << std::endl;
    std::cout << "  [N] normalisation (analytic)" << std::endl;
    std::cout << "  dR_min = " << dR_min << std::endl << std::endl;

    // ---- Read TARGET sample ----
    std::vector<double> c(nphys, 0.0);
    double tgt_sumW = 0.0;
    size_t tgt_njets = 0;
    {
        TFile tf(target_input.c_str(), "READ");
        if (tf.IsZombie()) die("failed to open target file: " + target_input);
        TTree *tt = dynamic_cast<TTree*>(tf.Get(tgt_tree_n.c_str()));
        if (!tt) die("failed to find target tree: " + tgt_tree_n);
        const Long64_t te = tt->GetEntries();
        const bool tev = (tt->GetBranch(n_branch.c_str()) != nullptr);

        if (tev) {
            int tnJ = 0; float tpt[100] = {0}; float tw = 1.0f;
            std::vector<std::vector<double>> *tcp=0,*tce=0,*tcf=0;
            tt->SetBranchAddress(n_branch.c_str(), &tnJ);
            tt->SetBranchAddress(pt_branch.c_str(), tpt);
            tt->SetBranchAddress(wt_branch.c_str(), &tw);
            tt->SetBranchAddress("const_pt", &tcp);
            tt->SetBranchAddress("const_eta", &tce);
            tt->SetBranchAddress("const_phi", &tcf);
            for (Long64_t ev = 0; ev < te; ++ev) {
                tt->GetEntry(ev);
                for (int j = 0; j < tnJ; ++j) {
                    double x = tpt[j];
                    if (x < pt_min || x > pt_max) continue;
                    if (!tcp||!tce||!tcf) die("missing target branches");
                    if (j>=(int)tcp->size()||j>=(int)tce->size()||j>=(int)tcf->size()) continue;
                    double w = (double)tw; tgt_sumW += w; ++tgt_njets;
                    auto gv = evaluate_basis(basis, x, tcp->at(j), tce->at(j), tcf->at(j), dR_min);
                    for (int b = 0; b < nphys; ++b) c[b] += w * gv[b];
                }
            }
        } else {
            float tp = 0.0f; float tw = 1.0f;
            std::vector<double> *tcp=0,*tce=0,*tcf=0;
            tt->SetBranchAddress(pt_branch.c_str(), &tp);
            tt->SetBranchAddress(wt_branch.c_str(), &tw);
            tt->SetBranchAddress("const_pt", &tcp);
            tt->SetBranchAddress("const_eta", &tce);
            tt->SetBranchAddress("const_phi", &tcf);
            for (Long64_t ev = 0; ev < te; ++ev) {
                tt->GetEntry(ev);
                double x = tp;
                if (x < pt_min || x > pt_max) continue;
                if (!tcp||!tce||!tcf) die("missing target branches");
                double w = (double)tw; tgt_sumW += w; ++tgt_njets;
                auto gv = evaluate_basis(basis, x, *tcp, *tce, *tcf, dR_min);
                for (int b = 0; b < nphys; ++b) c[b] += w * gv[b];
            }
        }
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
    const Long64_t ne = tree->GetEntries();
    const bool iev = (tree->GetBranch(n_branch.c_str()) != nullptr);
    const bool hsb = (tree->GetBranch("source") != nullptr);

    double wX_an = 0, wYp_an = 0; bool has_an = false;
    if (TTree *tW = dynamic_cast<TTree*>(in_file.Get("tWeights"))) {
        tW->SetBranchAddress("wX", &wX_an);
        tW->SetBranchAddress("wYprime", &wYp_an);
        if (tW->GetEntries()>0) { tW->GetEntry(0); has_an=true; }
    }

    std::vector<double> pts, base_w;
    std::vector<std::vector<double>> gvals;
    std::vector<int> sources;
    std::vector<std::vector<double>> ocp, oce, ocf;  // constituent arrays for output

    pts.reserve(ne*2); base_w.reserve(ne*2);

    auto read_jets = [&](auto get_jet_data) {
        // Generic lambda to avoid duplicating event vs flat logic.
        // get_jet_data fills pts/base_w/gvals/sources/ocp/oce/ocf.
    };

    if (iev) {
        int nJ=0; float pt[100]={0}; float wt=1.0f; int src=-1;
        std::vector<std::vector<double>> *cp=0,*ce=0,*cf=0;
        bool sc = tree->GetBranch("const_pt") && tree->GetBranch("const_eta")
                  && tree->GetBranch("const_phi");
        tree->SetBranchAddress(n_branch.c_str(), &nJ);
        tree->SetBranchAddress(pt_branch.c_str(), pt);
        tree->SetBranchAddress(wt_branch.c_str(), &wt);
        if (sc) { tree->SetBranchAddress("const_pt",&cp);
                  tree->SetBranchAddress("const_eta",&ce);
                  tree->SetBranchAddress("const_phi",&cf); }
        if (hsb) tree->SetBranchAddress("source", &src);
        for (Long64_t ev=0; ev<ne; ++ev) {
            tree->GetEntry(ev);
            for (int j=0; j<nJ; ++j) {
                double x=pt[j]; if (x<pt_min||x>pt_max) continue;
                if (!sc||!cp||!ce||!cf) die("missing constituent branches");
                if (j>=(int)cp->size()||j>=(int)ce->size()||j>=(int)cf->size()) continue;
                pts.push_back(x); base_w.push_back((double)wt); sources.push_back(src);
                gvals.push_back(evaluate_basis(basis,x,cp->at(j),ce->at(j),cf->at(j),dR_min));
                ocp.push_back(cp->at(j)); oce.push_back(ce->at(j)); ocf.push_back(cf->at(j));
            }
        }
    } else {
        float pv=0; float wt=1.0f; int src=-1;
        std::vector<double> *cp=0,*ce=0,*cf=0;
        bool sc = tree->GetBranch("const_pt") && tree->GetBranch("const_eta")
                  && tree->GetBranch("const_phi");
        tree->SetBranchAddress(pt_branch.c_str(), &pv);
        tree->SetBranchAddress(wt_branch.c_str(), &wt);
        if (sc) { tree->SetBranchAddress("const_pt",&cp);
                  tree->SetBranchAddress("const_eta",&ce);
                  tree->SetBranchAddress("const_phi",&cf); }
        if (hsb) tree->SetBranchAddress("source", &src);
        for (Long64_t ev=0; ev<ne; ++ev) {
            tree->GetEntry(ev);
            double x=pv; if (x<pt_min||x>pt_max) continue;
            if (!sc||!cp||!ce||!cf) die("missing constituent branches");
            pts.push_back(x); base_w.push_back((double)wt); sources.push_back(src);
            gvals.push_back(evaluate_basis(basis,x,*cp,*ce,*cf,dR_min));
            ocp.push_back(*cp); oce.push_back(*ce); ocf.push_back(*cf);
        }
    }

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
        std::cout << "Basis stats (scaled, dR_min=" << dR_min << "):" << std::endl;
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
    std::vector<double> xcp,xce,xcf;
    bool sco = (ocp.size()==N);
    tw.Branch("pt",&op,"pt/F");
    tw.Branch("w_unbias",&wu,"w_unbias/D");
    tw.Branch("w_base",&wb,"w_base/D");
    tw.Branch("w_total",&wt2,"w_total/D");
    tw.Branch("source",&os,"source/I");
    if (has_an) tw.Branch("w_analytic",&wa,"w_analytic/D");
    if (sco) { tw.Branch("const_pt",&xcp); tw.Branch("const_eta",&xce); tw.Branch("const_phi",&xcf); }

    for (size_t i=0; i<N; ++i) {
        double dot=0;
        for (int j=0; j<nphys; ++j) dot += lam[j]*gvals[i][j];
        double dc = std::max(-max_exp, std::min(max_exp, dot));
        wu = std::exp(-dc - lam_norm);
        wb = base_w[i]; wt2 = wb*wu; op = (float)pts[i];
        os = (i<sources.size()) ? sources[i] : -1;
        if (has_an) wa = (os==0) ? wX_an : wYp_an;
        if (sco) { xcp=ocp[i]; xce=oce[i]; xcf=ocf[i]; }
        tw.Fill();
    }

    TTree tm("meta","fit metadata");
    auto la = lam; la.push_back(lam_norm);
    std::string bl = "eec_5+norm_analytic";
    tm.Branch("lambda",&la); tm.Branch("lambda_norm",&lam_norm,"lambda_norm/D");
    tm.Branch("basis",&bl);
    tm.Branch("pt_min",const_cast<double*>(&pt_min),"pt_min/D");
    tm.Branch("pt_max",const_cast<double*>(&pt_max),"pt_max/D");
    tm.Branch("best_loss",const_cast<double*>(&best_loss),"best_loss/D");
    double dr_out = dR_min;
    tm.Branch("dR_min",&dr_out,"dR_min/D");
    tm.Fill();
    tw.Write(); tm.Write(); of.Close(); in_file.Close();

    // ---- Verification ----
    std::cout << "\nFinal lambdas:" << std::endl;
    for (int j=0; j<nphys; ++j)
        std::cout << "  [" << j << "] " << std::setw(14) << lam[j]
                  << "  (" << basis[j].label << ")" << std::endl;
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
