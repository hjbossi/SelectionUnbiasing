// plot_unbias_weights.cc
// ROOT macro to verify unbiasing by comparing unweighted vs weighted pT,
// EEC, and arbitrary basis-vector distributions, at paper-figure quality.
//
// Builds: c++ -std=c++17 -O2 plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights
//
// Example (defaults):
//   ./plot_unbias_weights \
//     --input unbias_weights.root \
//     --target-input startBasis_output.root --target-tree tRef \
//     --pt-min 100 --pt-max 140
//
// Example (paper figure: bigger text, PNG+PDF, a curated subset of panels):
//   ./plot_unbias_weights \
//     --input unbias_weights.root --target-input startBasis_output.root \
//     --label-size 0.050 --title-size 0.055 --legend-size 0.042 \
//     --line-width 3 --marker-size 1.2 \
//     --formats pdf,png \
//     --basis-filter sjmom_R0.10 --basis-ncol 3 --basis-nbins 25
//
// =====================================================================
// STYLE / CONFIGURABILITY
// =====================================================================
// Every plot below is built through a small set of shared helpers
// (PlotConfig, apply_paper_style(), style_hist(), make_ratio_pads(),
// save_canvas()) so that:
//   * axis-label size, axis-title size, legend text size, line width,
//     marker size, and pad margins are single numbers set once via CLI
//     flags and applied everywhere (paper-quality defaults are already
//     larger than ROOT's defaults; tune with --label-size etc.);
//   * every histogram that is filled with a physical event weight calls
//     Sumw2() *before* Fill(), so the statistical error bars shown (and
//     propagated into every ratio panel via TH1::Divide) reflect the true
//     sum-of-weights-squared uncertainty rather than a naive sqrt(N) count;
//   * every comparison histogram is drawn as points with error bars
//     ("E1"/"PE"), never as a bare "hist" line that hides the uncertainty;
//   * canvases can be saved in multiple formats at once (--formats
//     pdf,png,...).
//
// BASIS-VECTOR PANEL CONFIGURABILITY
// -----------------------------------
// unbias_weights.cc can fit a basis with dozens of terms (EEC grid x
// several subjet radii x powers). Previously this macro always drew every
// single one in an auto-sized grid with a fixed 20-bin log axis. Now:
//   --basis-select "0,3,7"   only plot these basis-function indices
//                            (indices follow meta.basis_labels order)
//   --basis-filter <substr>  only plot basis functions whose label
//                            contains this substring (e.g. "eec" or
//                            "sjmom_R0.10"); combine with --basis-select
//                            to intersect both conditions
//   --basis-max <N>          hard cap on the number of panels drawn (after
//                            select/filter), so a curated paper figure
//                            never balloons back out to 30+ panels
//   --basis-nbins <N>        bins per panel (default 20)
//   --basis-ncol <N>         grid columns (0 = auto sqrt-based layout)
//   --basis-pad <val>        log-decade padding around each panel's
//                            data-driven axis range (default 0.10)
//   --basis-cellw/--basis-cellh  pixel size of each grid cell
//
// EEC normalization ("bin center"): the EEC diagnostic histograms below
// normalize pt_i*pt_k by a momentum scale that must match whatever
// unbias_weights.cc used when it fit the weights. That number is written
// once by startBasis.C into a shared TEnv config file (default
// "unbiasing_config.env", key "Unbiasing.BinCenter") and read back here via
// basis_functions.h's read_bin_center() -- pass --config to point at a
// different file, or --eec-norm to override the value outright.

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TMath.h>
#include <TLine.h>
#include <TAxis.h>

#include "basis_functions.h"   // read_bin_center() (shared TEnv config helper)

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

static void die(const std::string &msg) {
  std::cerr << "error: " << msg << std::endl;
  std::exit(1);
}

static std::string get_arg(int argc, char* argv[], const std::string &flag, const std::string &def = "") {
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == flag && i + 1 < argc) return argv[i + 1];
  }
  return def;
}

static bool has_arg(int argc, char* argv[], const std::string &flag) {
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == flag) return true;
  return false;
}

// Parse a comma-separated list of non-negative integers, e.g. "0,3,7".
static std::set<int> parse_int_set(const std::string &s) {
  std::set<int> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    const size_t a = tok.find_first_not_of(" \t");
    if (a == std::string::npos) continue;
    const size_t b = tok.find_last_not_of(" \t");
    out.insert(std::stoi(tok.substr(a, b - a + 1)));
  }
  return out;
}

// Parse a comma-separated list of tokens, e.g. "pdf,png".
static std::vector<std::string> parse_string_list(const std::string &s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    const size_t a = tok.find_first_not_of(" \t");
    if (a == std::string::npos) continue;
    const size_t b = tok.find_last_not_of(" \t");
    out.push_back(tok.substr(a, b - a + 1));
  }
  return out;
}

// =====================================================================
// Shared plot configuration ("paper style")
// =====================================================================
struct PlotConfig {
  // Text sizes are in NDC pad-fraction units (ROOT convention), set for a
  // single, full-height pad; make_ratio_pads() below rescales them for the
  // shrunken ratio sub-pad so the printed text is visually the same size.
  double label_size  = 0.045;   // axis tick-label size
  double title_size  = 0.052;   // axis title size
  double legend_size = 0.040;   // legend entry text size
  double line_width  = 2.5;
  double marker_size = 1.1;
  double margin_left   = 0.14;
  double margin_right  = 0.05;
  double margin_top    = 0.06;
  double margin_bottom = 0.14;
  int    font          = 42;    // ROOT font code (42 = Helvetica, precision 2)
  std::vector<std::string> formats = { "pdf" };
};

static PlotConfig parse_plot_config(int argc, char* argv[]) {
  PlotConfig cfg;
  cfg.label_size  = std::stod(get_arg(argc, argv, "--label-size",  "0.045"));
  cfg.title_size  = std::stod(get_arg(argc, argv, "--title-size",  "0.052"));
  cfg.legend_size = std::stod(get_arg(argc, argv, "--legend-size", "0.040"));
  cfg.line_width  = std::stod(get_arg(argc, argv, "--line-width",  "2.5"));
  cfg.marker_size = std::stod(get_arg(argc, argv, "--marker-size", "1.1"));
  cfg.margin_left   = std::stod(get_arg(argc, argv, "--margin-left",   "0.14"));
  cfg.margin_right  = std::stod(get_arg(argc, argv, "--margin-right",  "0.05"));
  cfg.margin_top    = std::stod(get_arg(argc, argv, "--margin-top",    "0.06"));
  cfg.margin_bottom = std::stod(get_arg(argc, argv, "--margin-bottom", "0.14"));
  cfg.font        = std::stoi(get_arg(argc, argv, "--font", "42"));
  cfg.formats     = parse_string_list(get_arg(argc, argv, "--formats", "pdf"));
  if (cfg.formats.empty()) cfg.formats = { "pdf" };
  return cfg;
}

// Global gStyle setup applied once, up front, so every canvas inherits it.
static void apply_paper_style(const PlotConfig &cfg) {
  gStyle->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(cfg.font, "XYZ");
  gStyle->SetLabelFont(cfg.font, "XYZ");
  gStyle->SetTextFont(cfg.font);
  gStyle->SetTitleSize(cfg.title_size, "XYZ");
  gStyle->SetLabelSize(cfg.label_size, "XYZ");
  gStyle->SetTitleOffset(1.25, "X");
  gStyle->SetTitleOffset(1.35, "Y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetEndErrorSize(4);     // visible caps on statistical error bars
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(cfg.font);
}

// Apply consistent line/marker style + paper-size fonts to one histogram.
// Sumw2() is intentionally NOT called here -- it must be called before the
// first Fill() (see call sites below), not after.
static void style_hist(TH1D *h, int color, int marker, const PlotConfig &cfg) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetLineWidth(cfg.line_width);
  h->SetMarkerSize(cfg.marker_size);
  h->GetXaxis()->SetTitleSize(cfg.title_size);
  h->GetXaxis()->SetLabelSize(cfg.label_size);
  h->GetYaxis()->SetTitleSize(cfg.title_size);
  h->GetYaxis()->SetLabelSize(cfg.label_size);
  h->GetXaxis()->SetTitleFont(cfg.font);
  h->GetYaxis()->SetTitleFont(cfg.font);
  h->GetXaxis()->SetLabelFont(cfg.font);
  h->GetYaxis()->SetLabelFont(cfg.font);
}

static TLegend *make_legend(double x1, double y1, double x2, double y2, const PlotConfig &cfg) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(cfg.legend_size);
  leg->SetTextFont(cfg.font);
  return leg;
}

// Build the standard top(70%)/bottom(30%) main+ratio pad pair. Font sizes on
// the bottom pad are rescaled by the pad-height ratio (~0.70/0.30) so text
// drawn there reads at the same physical size as the top pad's text.
static void make_ratio_pads(const char *base_name, TPad *&p_top, TPad *&p_bot,
                            const PlotConfig &cfg) {
  p_top = new TPad(Form("%s_top", base_name), "", 0.0, 0.32, 1.0, 1.0);
  p_bot = new TPad(Form("%s_bot", base_name), "", 0.0, 0.0, 1.0, 0.32);
  p_top->SetBottomMargin(0.02);
  p_top->SetTopMargin(cfg.margin_top / 0.68);
  p_top->SetLeftMargin(cfg.margin_left);
  p_top->SetRightMargin(cfg.margin_right);
  p_bot->SetTopMargin(0.03);
  p_bot->SetBottomMargin(cfg.margin_bottom / 0.32);
  p_bot->SetLeftMargin(cfg.margin_left);
  p_bot->SetRightMargin(cfg.margin_right);
  p_top->Draw();
  p_bot->Draw();
}

// Style the bottom ratio-panel histogram: rescaled fonts (see above), a
// fixed y-range, and a dashed line at ratio=1.
static void style_ratio_hist(TH1D *h, const std::string &xtitle, const PlotConfig &cfg,
                             double ymin = 0.5, double ymax = 1.5) {
  const double scale = 0.68 / 0.32;   // top-pad-height / bottom-pad-height
  h->SetTitle((std::string("; ") + xtitle + "; ratio to target").c_str());
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->GetYaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleSize(cfg.title_size * scale);
  h->GetYaxis()->SetTitleOffset(1.35 / scale);
  h->GetYaxis()->SetLabelSize(cfg.label_size * scale);
  h->GetXaxis()->SetTitleSize(cfg.title_size * scale);
  h->GetXaxis()->SetLabelSize(cfg.label_size * scale);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetXaxis()->SetTitleFont(cfg.font); h->GetYaxis()->SetTitleFont(cfg.font);
  h->GetXaxis()->SetLabelFont(cfg.font); h->GetYaxis()->SetLabelFont(cfg.font);
  h->SetLineWidth(cfg.line_width);
}

static void draw_unity_line(TH1D *ref_axis_hist) {
  TLine *l = new TLine(ref_axis_hist->GetXaxis()->GetXmin(), 1.0,
                       ref_axis_hist->GetXaxis()->GetXmax(), 1.0);
  l->SetLineStyle(2);
  l->SetLineColor(kGray + 2);
  l->Draw("same");
}

// Save a canvas in every configured format, e.g. {"pdf","png"} ->
// base_name.pdf, base_name.png.
static void save_canvas(TCanvas *c, const std::string &base_name, const PlotConfig &cfg) {
  for (const std::string &fmt : cfg.formats)
    c->SaveAs((base_name + "." + fmt).c_str());
}

// =====================================================================
// BASIS VECTORS
// =====================================================================
// This macro does not re-derive the basis.  unbias_weights.cc writes the
// per-jet basis values g_k directly into its output — source jets in
// tweights.g_basis, target jets in the tbasis_target tree, with per-function
// titles in meta.basis_labels — so the basis-vector plots below always match
// exactly the basis that was fit (EEC terms, several subjet radii, custom
// powers, ...), with no reclustering and no dependence on any fixed basis.
// Which of those (possibly dozens of) functions get their own panel is
// controlled by --basis-select / --basis-filter / --basis-max (see top of
// file).

static void fill_eec(TH1D *h, const std::vector<double> &cpt,
                     const std::vector<double> &ceta,
                     const std::vector<double> &cphi,
                     double wjet, double eec_norm) {
  const size_t nconst = cpt.size();
  if (nconst < 2) return;
  const double norm = 1.0 / (eec_norm * eec_norm);   // z = pt_i pt_k / eec_norm^2
  for (size_t i = 0; i < nconst; ++i) {
    if (cpt[i] <= 0.0) continue;
    for (size_t k = i + 1; k < nconst; ++k) {
      if (cpt[k] <= 0.0) continue;
      double dphi = TVector2::Phi_mpi_pi(cphi[i] - cphi[k]);
      double deta = ceta[i] - ceta[k];
      double theta = std::sqrt(deta * deta + dphi * dphi);
      double w = (cpt[i] * cpt[k]) * norm;
      h->Fill(theta, wjet * w);
    }
  }
}

int main(int argc, char* argv[]) {
  const std::string input = get_arg(argc, argv, "--input");
  const std::string out_name = get_arg(argc, argv, "--out", "unbias_weights_plots.root");
  const std::string target_input = get_arg(argc, argv, "--target-input");
  const std::string target_tree = get_arg(argc, argv, "--target-tree", "tX");
  const double pt_min = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
  const double pt_max = std::stod(get_arg(argc, argv, "--pt-max", "200.0"));
  const int nbins = std::stoi(get_arg(argc, argv, "--nbins", "80"));
  const int eec_bins = std::stoi(get_arg(argc, argv, "--eec-bins", "60"));
  const double eec_max = std::stod(get_arg(argc, argv, "--eec-max", "0.4"));
  const double eec_logmin_val = std::stod(get_arg(argc, argv, "--eec-logmin", "0.004"));

  // ---- Basis-panel configurability ----
  const bool has_basis_select = has_arg(argc, argv, "--basis-select");
  const std::set<int> basis_select = has_basis_select
      ? parse_int_set(get_arg(argc, argv, "--basis-select", ""))
      : std::set<int>();
  const std::string basis_filter = get_arg(argc, argv, "--basis-filter", "");
  const int basis_max     = std::stoi(get_arg(argc, argv, "--basis-max", "0"));  // 0 = no cap
  const int basis_nbins   = std::stoi(get_arg(argc, argv, "--basis-nbins", "20"));
  const int basis_ncol_in = std::stoi(get_arg(argc, argv, "--basis-ncol", "0")); // 0 = auto
  const double basis_pad  = std::stod(get_arg(argc, argv, "--basis-pad", "0.10"));
  const int basis_cellw   = std::stoi(get_arg(argc, argv, "--basis-cellw", "520"));
  const int basis_cellh   = std::stoi(get_arg(argc, argv, "--basis-cellh", "480"));

  const PlotConfig cfg = parse_plot_config(argc, argv);
  apply_paper_style(cfg);

  // ---- Shared bin-center config (written by startBasis.C) ----
  // Same TEnv file unbias_weights.cc reads, so the EEC z-normalization used
  // for these diagnostic plots always matches what the fit itself used,
  // unless explicitly overridden with --eec-norm.
  const std::string config_file = get_arg(argc, argv, "--config", "unbiasing_config.env");
  const double bin_center = read_bin_center(config_file, 120.0);
  const double eec_norm = has_arg(argc, argv, "--eec-norm")
                               ? std::stod(get_arg(argc, argv, "--eec-norm", "120.0"))
                               : bin_center;
  std::cout << "Bin-center config: '" << config_file << "'  BinCenter=" << bin_center
            << "  (EEC norm used here = " << eec_norm << ")" << std::endl;

  if (input.empty()) die("--input is required");
  if (target_input.empty()) die("--target-input is required to compare against unbiased target");


  TFile in(input.c_str(), "READ");
  if (in.IsZombie()) die("failed to open input file: " + input);
  TTree *t = dynamic_cast<TTree*>(in.Get("tweights"));
  if (!t) die("missing tree 'tweights' in input file");

  float pt = 0.0f;
  double w_unbias = 1.0;
  double w_base = 1.0;
  double w_total = 1.0;
  std::vector<double> *const_pt = nullptr;
  std::vector<double> *const_eta = nullptr;
  std::vector<double> *const_phi = nullptr;

  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("w_unbias", &w_unbias);
  t->SetBranchAddress("w_base", &w_base);
  t->SetBranchAddress("w_total", &w_total);
  const bool has_const = (t->GetBranch("const_pt") != nullptr &&
                          t->GetBranch("const_eta") != nullptr &&
                          t->GetBranch("const_phi") != nullptr);
  if (has_const) {
    t->SetBranchAddress("const_pt", &const_pt);
    t->SetBranchAddress("const_eta", &const_eta);
    t->SetBranchAddress("const_phi", &const_phi);
  }

  // ---- Histograms ----
  // Sumw2() is called immediately after construction, *before* any Fill(),
  // on every histogram filled with a physical (non-unity) weight. Without
  // this, ROOT reports Poisson sqrt(N)-style errors for weighted content,
  // which understates/misstates the true statistical uncertainty and is
  // silently wrong for anything downstream that reads these bin errors
  // (including every ratio panel below, via TH1::Divide).
  TH1D *h_base = new TH1D("h_base", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_total = new TH1D("h_total", "pT; p_{T} [GeV]; weighted entries", nbins, pt_min, pt_max);
  TH1D *h_target = new TH1D("h_target", "pT; p_{T} [GeV]; entries", nbins, pt_min, pt_max);
  TH1D *h_w_unbias = new TH1D("h_w_unbias", "unbias weight; w; entries", 80, 0.0, 2.0);
  h_base->Sumw2(); h_total->Sumw2(); h_target->Sumw2(); h_w_unbias->Sumw2();

  // log-spaced theta binning for the EEC histograms
  const double eec_logmin = TMath::Log10(eec_logmin_val);
  const double eec_logmax = TMath::Log10(eec_max);
  const double eec_binwidth = (eec_logmax - eec_logmin) / eec_bins;
  std::vector<double> eec_edges(eec_bins + 1);
  for (int i = 0; i <= eec_bins; ++i)
    eec_edges[i] = std::pow(10.0, eec_logmin + i * eec_binwidth);

  TH1D *h_eec_base   = new TH1D("h_eec_base",   "EEC; #theta; EEC", eec_bins, eec_edges.data());
  TH1D *h_eec_total  = new TH1D("h_eec_total",  "EEC; #theta; EEC", eec_bins, eec_edges.data());
  TH1D *h_eec_target = new TH1D("h_eec_target", "EEC; #theta; EEC", eec_bins, eec_edges.data());
  h_eec_base->Sumw2(); h_eec_total->Sumw2(); h_eec_target->Sumw2();

  const Long64_t nentries = t->GetEntries();
  double sumw_base = 0.0;
  double sumw_total = 0.0;
  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    h_base->Fill(pt, w_base);
    h_total->Fill(pt, w_total);
    h_w_unbias->Fill(w_unbias);
    if (pt >= pt_min && pt <= pt_max && has_const && const_pt && const_eta && const_phi) {
      if (const_pt->size() == const_eta->size() && const_pt->size() == const_phi->size()) {
        fill_eec(h_eec_base, *const_pt, *const_eta, *const_phi, w_base, eec_norm);
        fill_eec(h_eec_total, *const_pt, *const_eta, *const_phi, w_total, eec_norm);
        sumw_base += w_base;
        sumw_total += w_total;
      }
    }
  }

  // Load unbiased target tree (tX) for comparison
  TFile tfin(target_input.c_str(), "READ");
  if (tfin.IsZombie()) die("failed to open target input file: " + target_input);
  TTree *tt = dynamic_cast<TTree*>(tfin.Get(target_tree.c_str()));
  if (!tt) die("missing target tree: " + target_tree);

  float tpt = 0.0f;
  float tw = 1.0f;
  std::vector<double> *tconst_pt = nullptr;
  std::vector<double> *tconst_eta = nullptr;
  std::vector<double> *tconst_phi = nullptr;
  tt->SetBranchAddress("pt", &tpt);
  tt->SetBranchAddress("weight", &tw);
  const bool target_has_const = (tt->GetBranch("const_pt") != nullptr &&
                                 tt->GetBranch("const_eta") != nullptr &&
                                 tt->GetBranch("const_phi") != nullptr);
  if (target_has_const) {
    tt->SetBranchAddress("const_pt", &tconst_pt);
    tt->SetBranchAddress("const_eta", &tconst_eta);
    tt->SetBranchAddress("const_phi", &tconst_phi);
  }
  const Long64_t tentries = tt->GetEntries();
  double sumw_target = 0.0;
  for (Long64_t i = 0; i < tentries; ++i) {
    tt->GetEntry(i);
    h_target->Fill(tpt, tw);
    if (tpt >= pt_min && tpt <= pt_max && target_has_const && tconst_pt && tconst_eta && tconst_phi) {
      if (tconst_pt->size() == tconst_eta->size() && tconst_pt->size() == tconst_phi->size()) {
        fill_eec(h_eec_target, *tconst_pt, *tconst_eta, *tconst_phi, tw, eec_norm);
        sumw_target += tw;
      }
    }
  }

  style_hist(h_base,   kBlue + 1,  20, cfg);
  style_hist(h_total,  kRed + 1,   21, cfg);
  style_hist(h_target, kGreen + 2, 22, cfg);
  style_hist(h_eec_base,   kBlue + 1,  20, cfg);
  style_hist(h_eec_total,  kRed + 1,   21, cfg);
  style_hist(h_eec_target, kGreen + 2, 22, cfg);

  // =====================================================================
  // pT check plot (main + ratio-to-target panels)
  // =====================================================================
  TCanvas *c = new TCanvas("c_unbias", "unbias check", 900, 900);
  TPad *p_top = nullptr, *p_bot = nullptr;
  make_ratio_pads("p_pt", p_top, p_bot, cfg);

  p_top->cd();
  gPad->SetLogy();
  h_target->Draw("E1");   // draw target first/underneath so it sets the axis range
  h_base->Draw("E1 same");
  h_total->Draw("E1 same");

  TLegend *leg = make_legend(0.50, 0.68, 0.90, 0.90, cfg);
  leg->AddEntry(h_base, "biased sample (base weight)", "lp");
  leg->AddEntry(h_total, "weighted sample (base x unbias)", "lp");
  leg->AddEntry(h_target, Form("unbiased target %s", target_tree.c_str()), "lp");
  leg->Draw();

  p_bot->cd();
  // Errors on the ratio come from standard uncorrelated TH1::Divide error
  // propagation; strictly, the biased/weighted samples and the target share
  // some parent jets (they are resampled from a common pool), so this is a
  // display-quality approximation, not a fully correlated treatment.
  TH1D *h_ratio_base = (TH1D*)h_base->Clone("h_ratio_base_pt");
  TH1D *h_ratio_total = (TH1D*)h_total->Clone("h_ratio_total_pt");
  h_ratio_base->Divide(h_target);
  h_ratio_total->Divide(h_target);
  style_ratio_hist(h_ratio_total, "p_{T} [GeV]", cfg);
  h_ratio_total->SetLineColor(kRed + 1);   h_ratio_total->SetMarkerColor(kRed + 1);
  h_ratio_base->SetLineColor(kBlue + 1);   h_ratio_base->SetMarkerColor(kBlue + 1);
  h_ratio_total->Draw("E1");
  h_ratio_base->Draw("E1 same");
  draw_unity_line(h_ratio_total);

  // Normalize EEC by total jet weight (per-jet average)
  if (sumw_base > 0.0) h_eec_base->Scale(1.0 / sumw_base);
  if (sumw_total > 0.0) h_eec_total->Scale(1.0 / sumw_total);
  if (sumw_target > 0.0) h_eec_target->Scale(1.0 / sumw_target);

  TFile out(out_name.c_str(), "RECREATE");
  h_base->Write();
  h_total->Write();
  h_target->Write();
  h_w_unbias->Write();
  h_eec_base->Write();
  h_eec_total->Write();
  h_eec_target->Write();
  c->Write();

  save_canvas(c, "plot_unbias_pT_check", cfg);

  // =====================================================================
  // Weight distribution
  // =====================================================================
  TCanvas *cw = new TCanvas("c_weights", "weight check", 750, 650);
  cw->SetLogy();
  cw->SetLeftMargin(cfg.margin_left); cw->SetRightMargin(cfg.margin_right);
  cw->SetTopMargin(cfg.margin_top);   cw->SetBottomMargin(cfg.margin_bottom);
  h_w_unbias->Draw("E1");
  TLegend *lw = make_legend(0.55, 0.78, 0.90, 0.90, cfg);
  lw->AddEntry(h_w_unbias, "unbias weights", "lp");
  lw->Draw();
  save_canvas(cw, "plot_unbias_weights_weights", cfg);
  delete cw;

  // =====================================================================
  // EEC comparison plot (main + ratio-to-target panels)
  // =====================================================================
  TCanvas *ce = new TCanvas("c_eec", "EEC comparison", 850, 850);
  TPad *p_eec_top = nullptr, *p_eec_bot = nullptr;
  make_ratio_pads("p_eec", p_eec_top, p_eec_bot, cfg);

  p_eec_top->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  h_eec_target->Draw("E1");
  h_eec_base->Draw("E1 same");
  h_eec_total->Draw("E1 same");
  TLegend *lege = make_legend(0.50, 0.20, 0.90, 0.42, cfg);
  lege->AddEntry(h_eec_base, "biased sample (base weight)", "lp");
  lege->AddEntry(h_eec_total, "weighted sample (base x unbias)", "lp");
  lege->AddEntry(h_eec_target, Form("unbiased target %s", target_tree.c_str()), "lp");
  lege->Draw();

  p_eec_bot->cd();
  gPad->SetLogx();
  TH1D *h_ratio_base_eec = (TH1D*)h_eec_base->Clone("h_ratio_base_eec");
  TH1D *h_ratio_total_eec = (TH1D*)h_eec_total->Clone("h_ratio_total_eec");
  h_ratio_base_eec->Divide(h_eec_target);
  h_ratio_total_eec->Divide(h_eec_target);
  style_ratio_hist(h_ratio_total_eec, "#theta", cfg);
  h_ratio_total_eec->SetLineColor(kRed + 1);  h_ratio_total_eec->SetMarkerColor(kRed + 1);
  h_ratio_base_eec->SetLineColor(kBlue + 1);  h_ratio_base_eec->SetMarkerColor(kBlue + 1);
  h_ratio_total_eec->Draw("E1");
  h_ratio_base_eec->Draw("E1 same");
  draw_unity_line(h_ratio_total_eec);

  save_canvas(ce, "plot_eec_compare", cfg);
  ce->Write();

  out.Close();

  in.Close();
  tfin.Close();
  {

    // ---------------------------------------------------------------
    // Basis-vector closure plots, read STRAIGHT from the fit output:
    //   source jets : tweights.g_basis        (weighted by w_base, w_total)
    //   target jets : tbasis_target.g_basis    (weighted by weight)
    //   titles      : meta.basis_labels
    // No basis is re-derived here, so these plots always match the fit
    // exactly (EEC terms, multiple subjet radii, custom powers, ...).
    // Which functions get a panel is controlled by --basis-select /
    // --basis-filter / --basis-max (see top-of-file docs).
    // ---------------------------------------------------------------
    TH1::AddDirectory(kFALSE);   // histos below are transient, not file-owned
    TFile inb(input.c_str(), "READ");
    TTree *tsrc = inb.IsZombie() ? nullptr : dynamic_cast<TTree*>(inb.Get("tweights"));
    TTree *ttar = inb.IsZombie() ? nullptr : dynamic_cast<TTree*>(inb.Get("tbasis_target"));
    TTree *tmet = inb.IsZombie() ? nullptr : dynamic_cast<TTree*>(inb.Get("meta"));

    const bool have_src = tsrc && tsrc->GetBranch("g_basis");
    const bool have_tar = ttar && ttar->GetBranch("g_basis");
    if (!have_src || !have_tar) {
      std::cout << "warn: basis-vector plot skipped -- fit output is missing "
                << (!have_src ? "tweights.g_basis" : "tbasis_target.g_basis")
                << ".  Re-run unbias_weights.cc (it stores per-jet basis values)."
                << std::endl;
    } else {
      // per-function titles (optional)
      std::vector<std::string> *labels = nullptr;
      if (tmet && tmet->GetBranch("basis_labels")) {
        tmet->SetBranchAddress("basis_labels", &labels);
        tmet->GetEntry(0);
      }

      // source (tweights) branches
      std::vector<double> *sg = nullptr;
      float  spt = 0.0f; double sw_base = 1.0, sw_total = 1.0;
      tsrc->SetBranchAddress("g_basis", &sg);
      tsrc->SetBranchAddress("pt", &spt);
      tsrc->SetBranchAddress("w_base", &sw_base);
      tsrc->SetBranchAddress("w_total", &sw_total);

      // target (tbasis_target) branches
      std::vector<double> *tg = nullptr;
      double tgw = 1.0;
      ttar->SetBranchAddress("g_basis", &tg);
      ttar->SetBranchAddress("weight", &tgw);

      // number of basis functions from the first source entry
      int nphys_total = 0;
      if (tsrc->GetEntries() > 0) { tsrc->GetEntry(0); nphys_total = sg ? (int)sg->size() : 0; }

      if (nphys_total == 0) {
        std::cout << "warn: basis-vector plot skipped -- g_basis is empty." << std::endl;
      } else {
        auto label_of = [&](int k) -> std::string {
          if (labels && k < (int)labels->size() && !(*labels)[k].empty()) return (*labels)[k];
          return std::string(Form("g_%d", k));
        };

        // ---- Select which basis-function indices get a panel ----
        std::vector<int> panel_idx;
        for (int k = 0; k < nphys_total; ++k) {
          if (has_basis_select && basis_select.find(k) == basis_select.end()) continue;
          if (!basis_filter.empty() && label_of(k).find(basis_filter) == std::string::npos) continue;
          panel_idx.push_back(k);
        }
        if (panel_idx.empty()) {
          std::cout << "warn: --basis-select/--basis-filter matched zero basis functions "
                       "(of " << nphys_total << " available) -- skipping basis-vector plot."
                    << std::endl;
        } else {
          if (basis_max > 0 && (int)panel_idx.size() > basis_max) {
            std::cout << "note: capping basis-vector panels to --basis-max=" << basis_max
                      << " (of " << panel_idx.size() << " selected)" << std::endl;
            panel_idx.resize(basis_max);
          }
          const int nphys = (int)panel_idx.size();
          std::cout << "Basis-vector panels (" << nphys << " of " << nphys_total << "):" << std::endl;
          for (int k : panel_idx) std::cout << "  [" << k << "] " << label_of(k) << std::endl;

          const Long64_t ns = tsrc->GetEntries();
          const Long64_t nt = ttar->GetEntries();

          // pass 1: per-selected-function positive min/max for log ranges
          std::vector<double> gmin(nphys, std::numeric_limits<double>::infinity());
          std::vector<double> gmax(nphys, 0.0);
          auto scan = [&](std::vector<double> *g) {
            if (!g) return;
            for (int p = 0; p < nphys; ++p) {
              const int k = panel_idx[p];
              if (k >= (int)g->size()) continue;
              const double v = std::abs((*g)[k]);
              if (v > 0.0) { gmin[p] = std::min(gmin[p], v); gmax[p] = std::max(gmax[p], v); }
            }
          };
          for (Long64_t i = 0; i < ns; ++i) { tsrc->GetEntry(i); if (spt < pt_min || spt > pt_max) continue; scan(sg); }
          for (Long64_t i = 0; i < nt; ++i) { ttar->GetEntry(i); scan(tg); }

          // build histograms with data-driven log binning (handles tiny EEC and huge pT^n alike)
          std::vector<TH1D*> hbase(nphys), htot(nphys), htar(nphys);
          for (int p = 0; p < nphys; ++p) {
            double lo = std::isfinite(gmin[p]) ? gmin[p] : 1e-3;
            double hi = (gmax[p] > lo) ? gmax[p] : lo * 10.0;
            double logmin = std::log10(lo) - basis_pad;
            double logmax = std::log10(hi) + basis_pad;
            if (!(logmax > logmin)) { logmin = -3.0; logmax = 3.0; }
            std::vector<double> edges(basis_nbins + 1);
            for (int i = 0; i <= basis_nbins; ++i)
              edges[i] = std::pow(10.0, logmin + (logmax - logmin) * i / basis_nbins);
            const std::string lbl = label_of(panel_idx[p]);
            const std::string ttl = lbl + "; " + lbl + "; entries";
            hbase[p] = new TH1D(Form("hg_base_%d",   panel_idx[p]), ttl.c_str(), basis_nbins, edges.data());
            htot[p]  = new TH1D(Form("hg_total_%d",  panel_idx[p]), ttl.c_str(), basis_nbins, edges.data());
            htar[p]  = new TH1D(Form("hg_target_%d", panel_idx[p]), ttl.c_str(), basis_nbins, edges.data());
            for (TH1D *h : { hbase[p], htot[p], htar[p] }) h->Sumw2();   // before Fill()
            style_hist(hbase[p], kBlue + 1,  20, cfg);
            style_hist(htot[p],  kRed + 1,   21, cfg);
            style_hist(htar[p],  kGreen + 2, 22, cfg);
          }

          // pass 2: fill
          for (Long64_t i = 0; i < ns; ++i) {
            tsrc->GetEntry(i);
            if (spt < pt_min || spt > pt_max) continue;
            if (!sg) continue;
            for (int p = 0; p < nphys; ++p) {
              const int k = panel_idx[p];
              if (k >= (int)sg->size()) continue;
              const double v = std::abs((*sg)[k]);
              hbase[p]->Fill(v, sw_base);
              htot[p] ->Fill(v, sw_total);
            }
          }
          for (Long64_t i = 0; i < nt; ++i) {
            ttar->GetEntry(i);
            if (!tg) continue;
            for (int p = 0; p < nphys; ++p) {
              const int k = panel_idx[p];
              if (k >= (int)tg->size()) continue;
              htar[p]->Fill(std::abs((*tg)[k]), tgw);
            }
          }

          // canvas grid: --basis-ncol overrides the auto sqrt-based layout
          const int ncol = (basis_ncol_in > 0) ? basis_ncol_in
                                               : std::max(1, (int)std::ceil(std::sqrt((double)nphys)));
          const int nrow = (nphys + ncol - 1) / ncol;
          TCanvas *ctw = new TCanvas("c_basis_weighted", "basis vectors (weighted)",
                                     basis_cellw * ncol, basis_cellh * nrow);
          ctw->Divide(ncol, nrow);
          TLegend *legw = make_legend(0.10, 0.70, 0.92, 0.90, cfg);
          legw->AddEntry(htar[0],  Form("unbiased target %s", target_tree.c_str()), "lp");
          legw->AddEntry(hbase[0], "biased sample (base weight)", "lp");
          legw->AddEntry(htot[0],  "weighted sample (base x unbias)", "lp");

          for (int p = 0; p < nphys; ++p) {
            ctw->cd(p + 1);
            TPad *pp_top = nullptr, *pp_bot = nullptr;
            make_ratio_pads(Form("p_bv_%d", p), pp_top, pp_bot, cfg);

            pp_top->cd(); gPad->SetLogy(); gPad->SetLogx();
            htar[p]->Draw("E1"); hbase[p]->Draw("E1 same"); htot[p]->Draw("E1 same");
            if (p == 0) legw->Draw();

            pp_bot->cd(); gPad->SetLogx();
            TH1D *rb = (TH1D*)hbase[p]->Clone(Form("hg_ratio_base_%d",  panel_idx[p]));
            TH1D *rt = (TH1D*)htot[p] ->Clone(Form("hg_ratio_total_%d", panel_idx[p]));
            rb->Divide(htar[p]); rt->Divide(htar[p]);
            style_ratio_hist(rt, label_of(panel_idx[p]), cfg);
            rt->SetLineColor(kRed + 1);  rt->SetMarkerColor(kRed + 1);
            rb->SetLineColor(kBlue + 1); rb->SetMarkerColor(kBlue + 1);
            rt->Draw("E1"); rb->Draw("E1 same");
            draw_unity_line(rt);
          }
          save_canvas(ctw, "plot_theta_basis_weighted_compare", cfg);
        }
      }
    }
    inb.Close();
  }

  std::cout << "wrote " << out_name << " and plot_unbias_pT_check / "
            << "plot_unbias_weights_weights / plot_eec_compare / "
            << "plot_theta_basis_weighted_compare in: ";
  for (const std::string &fmt : cfg.formats) std::cout << "." << fmt << " ";
  std::cout << std::endl;
  return 0;
}