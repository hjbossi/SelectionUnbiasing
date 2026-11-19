// plot_unbias_weights.cc
// ROOT macro to verify unbiasing by comparing unweighted vs weighted pT,
// EEC, and arbitrary basis-vector distributions, at paper-figure quality.
// Also validates closure on several jet-substructure observables that are
// NOT part of the fit basis (jet mass, Soft-Drop groomed radius R_g and
// groomed-splitting k_T, N-subjettiness tau_1/2/3 (+ tau21/tau32), and the
// broadening/width form of the generalised angularities) -- see
// jet_observables.h and the "JET-SUBSTRUCTURE CROSSCHECKS" section below.
//
// Builds (NOTE: now requires linking FastJet, because jet_observables.h's
// Soft Drop groomer uses fastjet::ClusterSequence -- unlike before, ROOT
// flags alone are no longer enough):
//   c++ -std=c++17 -O2 plot_unbias_weights.cc \
//       $(root-config --cflags --libs) \
//       $($FJI/bin/fastjet-config --cxxflags --libs) \
//       -o plot_unbias_weights
//
// To also enable the N-subjettiness crosscheck (tau_1/2/3, tau21, tau32),
// which needs fastjet-contrib's separate Nsubjettiness package (NOT part of
// base FastJet), add -DJET_OBS_HAVE_NSUBJETTINESS and link it, e.g.:
//   c++ -std=c++17 -O2 plot_unbias_weights.cc \
//       $(root-config --cflags --libs) \
//       $($FJI/bin/fastjet-config --cxxflags --libs) \
//       -DJET_OBS_HAVE_NSUBJETTINESS \
//       -I/path/to/fjcontrib-install/include -L/path/to/fjcontrib-install/lib -lNsubjettiness \
//       -o plot_unbias_weights
// Without that define, the tau1/tau2/tau3/tau21/tau32 histograms and plots
// are simply skipped (one note printed to stdout), everything else is
// unaffected. See makefile for a build target with this wired up via
// an (empty-by-default) NSUB_FLAGS make variable.
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
//     --formats pdf,png --plot-dir ../Output/plots \
//     --basis-filter sjmom_R0.10 --basis-ncol 3 --basis-nbins 25
//
// =====================================================================
// STYLE / CONFIGURABILITY
// =====================================================================
// Every plot below is built through a small set of shared helpers
// (PlotConfig, apply_paper_style(), style_hist(), make_ratio_pads(),
// save_canvas(), observable_check_plot()) so that:
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
// PT-BIN LABELING (plot text + file names)
// -----------------------------------------
// Every plot produced below is labeled with its pT bin -- both as text
// drawn directly on the canvas (e.g. "100 < p_{T} < 140 GeV", top-left) and
// as a "_<ptmin>_<ptmax>" suffix on every saved file name (e.g.
// plot_unbias_pT_check_100_140.pdf) and on the default --out ROOT file name.
// Both are derived automatically from --pt-min/--pt-max -- there is no
// separate "bin label" flag to keep in sync by hand, so the label can never
// drift from the window that was actually plotted. Integer-valued bin edges
// are formatted without a decimal point (e.g. "100", not "100.0").
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
// unbias_weights.cc used when it fit the weights. There is no hard-coded
// numeric value for it anywhere in this codebase: it is written once by
// startBasis.C into a shared TEnv config file (default
// "unbiasing_config.env", key "Unbiasing.BinCenter") and read back here via
// basis_functions.h's read_bin_center() -- pass --config to point at a
// different file, or --eec-norm to override the value outright. If the
// config file/key is missing, read_bin_center() aborts with an error
// instead of guessing a number.
//
// JET-SUBSTRUCTURE CROSSCHECKS (jet mass, R_g, k_T, tau_N, broadening)
// ---------------------------------------------------------------------
// None of these observables are part of the fit basis (basis_functions.h)
// -- they are closure crosschecks only, exactly as requested: "does
// reweighting derived from EEC/subjet-moment constraints also produce
// closure on other substructure observables?" (see todo.md). They
// are computed directly from the same per-jet constituents already used
// for the EEC plot above (const_pt/eta/phi, massless convention), plus:
//   * jet mass: read directly from a stored "mass" branch (never
//     recomputed from constituents) -- source jets from --input's
//     "tweights" tree (only present if unbias_weights.cc's input tree
//     carried "eta"/"phi"/"mass"; see that file's "Jet-level eta/phi/mass
//     passthrough" note), target jets directly from --target-input's own
//     tree (already carries "mass" when produced by startBasis.C).
//   * jet radius R (--jet-R, default 0.4, matching the README's stated
//     PYTHIA jet radius -- NOT hard-coded, pass --jet-R to override):
//     used both for the broadening normalization (Delta R_i,jet / R) and
//     as the Soft Drop angular-exponent scale for R_g/k_T.
//   * Soft Drop z_cut/beta (--sd-zcut, default 0.1; --sd-beta, default
//     1.0): grooming parameters for R_g and k_T (kT taken at the same
//     groomed splitting that sets R_g). See jet_observables.h for the
//     exact recursive-declustering implementation.
//   * N-subjettiness beta (--nsub-beta, default 1.0): tau_1/2/3 with
//     one-pass winner-take-all kT axes, NORMALIZED measure (R0 = --jet-R).
//     Only computed/plotted if built with -DJET_OBS_HAVE_NSUBJETTINESS
//     (see build instructions above) -- otherwise skipped.
//   * broadening: the generalised-angularity width (beta=1, kappa=1),
//     using the STORED jet axis (same "eta"/"phi" as above), never an
//     axis re-derived from constituents.
// A jet that never satisfies the Soft Drop condition grooms away to a
// single particle; such jets are excluded from the R_g/k_T histograms
// (counted and reported separately, not silently merged into a bin) --
// see jet_observables.h's SoftDropResult::groomed_away.

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TMath.h>
#include <TLine.h>
#include <TAxis.h>
#include <TSystem.h>

#include "basis_functions.h"   // read_bin_center() (shared TEnv config helper)
#include "jet_observables.h"   // soft_drop() / nsubjettiness() / broadening() -- substructure crosschecks

#include <array>
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

// ---- pT-bin formatting (shared by file names and the on-canvas label) ----
// Integer-valued edges print without a decimal point (e.g. "100"); anything
// else prints with one decimal place (e.g. "137.5").
static std::string format_pt_for_name(double v) {
  char buf[32];
  if (std::abs(v - std::round(v)) < 1e-6)
    std::snprintf(buf, sizeof(buf), "%d", (int)std::llround(v));
  else
    std::snprintf(buf, sizeof(buf), "%.1f", v);
  return std::string(buf);
}

// Filename-safe suffix, e.g. "100_140" -- appended to every saved plot name
// and to the default --out ROOT file name.
static std::string pt_bin_suffix(double pt_min, double pt_max) {
  return format_pt_for_name(pt_min) + "_" + format_pt_for_name(pt_max);
}

// Human-readable on-canvas label, e.g. "100 < p_{T} < 140 GeV".
static std::string pt_bin_label(double pt_min, double pt_max) {
  return format_pt_for_name(pt_min) + " < p_{T} < " + format_pt_for_name(pt_max) + " GeV";
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
  std::string plot_dir = ".";   // directory the plot files are saved in (--plot-dir; created if missing)
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
  cfg.plot_dir    = get_arg(argc, argv, "--plot-dir", ".");
  if (cfg.plot_dir.empty()) cfg.plot_dir = ".";
  gSystem->mkdir(cfg.plot_dir.c_str(), kTRUE);
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

// Draw the pT-bin label (e.g. "100 < p_{T} < 140 GeV") in NDC coordinates of
// whichever pad is currently active. Called once per canvas/top-pad below so
// every plot is unambiguously tagged with the bin it was made from. Default
// (x,y) sits just inside the frame's top-left corner (not above/outside the
// axis box, which is where it used to land relative to the old, larger top
// margin) and clear of the legend, which is always placed on the right half
// of the pad.
static void draw_bin_label(const std::string &label, const PlotConfig &cfg,
                           double x = 0.18, double y = 0.86) {
  TLatex *lt = new TLatex();
  lt->SetNDC();
  lt->SetTextFont(cfg.font);
  lt->SetTextSize(cfg.legend_size);
  lt->DrawLatex(x, y, label.c_str());
}

// Build the standard top(70%)/bottom(30%) main+ratio pad pair, flush against
// each other (no gap): the top pad's bottom margin and the bottom pad's top
// margin are both zero, so the two frames share a single boundary line at
// y=0.32. Because of that, the top pad's histogram must never draw its own
// x-axis labels/title (see hide_shared_xaxis() below) -- the bottom pad's
// x-axis is the only one shown, right under the shared boundary.
static void make_ratio_pads(const char *base_name, TPad *&p_top, TPad *&p_bot,
                            const PlotConfig &cfg) {
  p_top = new TPad(Form("%s_top", base_name), "", 0.0, 0.32, 1.0, 1.0);
  p_bot = new TPad(Form("%s_bot", base_name), "", 0.0, 0.0, 1.0, 0.32);
  p_top->SetBottomMargin(0.0);
  p_top->SetTopMargin(cfg.margin_top / 0.68);
  p_top->SetLeftMargin(cfg.margin_left);
  p_top->SetRightMargin(cfg.margin_right);
  p_bot->SetTopMargin(0.0);
  p_bot->SetBottomMargin(cfg.margin_bottom / 0.32);
  p_bot->SetLeftMargin(cfg.margin_left);
  p_bot->SetRightMargin(cfg.margin_right);
  p_top->Draw();
  p_bot->Draw();
}

// Hide the x-axis labels/title on the histogram that sets the TOP pad's
// frame (always the one drawn first, e.g. "h_target" below) -- with the two
// pads now flush (make_ratio_pads above), the bottom pad's x-axis sits
// immediately under the shared boundary, so the top pad must not draw its
// own copy of it.
static void hide_shared_xaxis(TH1D *h_top_axis) {
  h_top_axis->GetXaxis()->SetLabelSize(0);
  h_top_axis->GetXaxis()->SetTitleSize(0);
}

// Give the top panel enough headroom above the tallest of the three
// histograms that the legend (drawn in the upper area of the pad) never
// overlaps the data. `h_axis` must be the histogram drawn first (it is the
// one whose Maximum actually sets the pad's y-range); `logy` picks a
// multiplicative (linear) vs. additive-in-log-decades (log) padding.
static void set_top_headroom(TH1D *h_axis, double max_val, bool logy,
                             double log_decades = 1.0, double linear_factor = 1.8) {
  if (!(max_val > 0.0)) return;
  if (logy) h_axis->SetMaximum(max_val * std::pow(10.0, log_decades));
  else      h_axis->SetMaximum(max_val * linear_factor);
}

// Style the bottom ratio-panel histogram: rescaled fonts (see above), a
// fixed y-range, and a dashed line at ratio=1.
static void style_ratio_hist(TH1D *h, const std::string &xtitle, const PlotConfig &cfg,
                             double ymin = 0.5, double ymax = 1.5) {
  const double scale = 0.68 / 0.32;   // top-pad-height / bottom-pad-height
  // The y-axis TITLE ("ratio to target") is long enough that scaling it up
  // by the full top/bottom pad-height ratio (to visually match the main
  // panel's title size) makes the rotated string taller than the ratio
  // pad itself, so ROOT clips it against the pad boundary (e.g. "ratio to
  // targe[t]"). Axis LABELS (just numbers, e.g. "0.8", "1.0") stay at the
  // full match-scale since they're short and never overflow; the title
  // uses a capped, smaller scale so the whole string always fits within
  // the pad's height.
  const double title_scale = std::min(scale, 1.45);
  h->SetTitle((std::string("; ") + xtitle + "; ratio to target").c_str());
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->GetYaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleSize(cfg.title_size * title_scale);
  h->GetYaxis()->SetTitleOffset(1.35 / title_scale);
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
// <plot_dir>/base_name.pdf, <plot_dir>/base_name.png.
// --formats none saves nothing (the histograms still go to --out), e.g. when
// the plots are drawn with matplotlib instead (plot_mpl.py, UNBIAS_PLOTTER).
static void save_canvas(TCanvas *c, const std::string &base_name, const PlotConfig &cfg) {
  for (const std::string &fmt : cfg.formats)
    if (fmt != "none") c->SaveAs((cfg.plot_dir + "/" + base_name + "." + fmt).c_str());
}

// =====================================================================
// Generic "base/total/target main+ratio" comparison canvas, used by every
// jet-substructure crosscheck below (jet mass, R_g, k_T, tau_N, tau_N
// ratios, broadening) so none of them repeats the ~30-line drawing block
// that the pT-check and EEC-comparison plots already use above. h_base/
// h_total/h_target must already be filled, Sumw2()'d, and styled (their
// axis titles come from however they were constructed, e.g.
// "; M_{jet} [GeV]; weighted entries").
// =====================================================================
static void observable_check_plot(const std::string &base_name,
                                  const std::string &xtitle,
                                  TH1D *h_base, TH1D *h_total, TH1D *h_target,
                                  bool logy, bool logx,
                                  const std::string &target_tree_name,
                                  const std::string &pt_label,
                                  const std::string &pt_suffix,
                                  const PlotConfig &cfg) {
  TCanvas *c = new TCanvas(("c_" + base_name).c_str(), base_name.c_str(), 850, 850);
  TPad *p_top = nullptr, *p_bot = nullptr;
  make_ratio_pads(("p_" + base_name).c_str(), p_top, p_bot, cfg);

  p_top->cd();
  if (logy) gPad->SetLogy();
  if (logx) gPad->SetLogx();
  set_top_headroom(h_target, std::max({h_base->GetMaximum(), h_total->GetMaximum(),
                                       h_target->GetMaximum()}), logy);
  hide_shared_xaxis(h_target);
  h_target->Draw("E1");
  h_base->Draw("E1 same");
  h_total->Draw("E1 same");
  TLegend *leg = make_legend(0.42, 0.68, 0.82, 0.90, cfg);
  leg->AddEntry(h_base, "biased sample (base weight)", "lp");
  leg->AddEntry(h_total, "weighted sample (base x unbias)", "lp");
  leg->AddEntry(h_target, Form("unbiased target %s", target_tree_name.c_str()), "lp");
  leg->Draw();
  draw_bin_label(pt_label, cfg);

  p_bot->cd();
  if (logx) gPad->SetLogx();
  TH1D *rb = (TH1D*)h_base->Clone((base_name + "_ratio_base").c_str());
  TH1D *rt = (TH1D*)h_total->Clone((base_name + "_ratio_total").c_str());
  rb->Divide(h_target);
  rt->Divide(h_target);
  style_ratio_hist(rt, xtitle, cfg);
  rt->SetLineColor(kRed + 1);  rt->SetMarkerColor(kRed + 1);
  rb->SetLineColor(kBlue + 1); rb->SetMarkerColor(kBlue + 1);
  rt->Draw("E1");
  rb->Draw("E1 same");
  draw_unity_line(rt);

  save_canvas(c, "plot_" + base_name + "_" + pt_suffix, cfg);
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
  const std::string target_input = get_arg(argc, argv, "--target-input");
  const std::string target_tree = get_arg(argc, argv, "--target-tree", "tX");
  const double pt_min = std::stod(get_arg(argc, argv, "--pt-min", "0.0"));
  const double pt_max = std::stod(get_arg(argc, argv, "--pt-max", "200.0"));

  // pT-bin suffix/label, derived from --pt-min/--pt-max (see top-of-file
  // docs): used below both for the default --out ROOT file name and for
  // every saved plot's file name and on-canvas label.
  // Optional --tag <name> (e.g. "lognormal") names the bias type; it is
  // appended to the suffix, so every saved plot and the default --out file
  // carry it (e.g. plot_eec_compare_100_140_lognormal.pdf).
  const std::string bias_tag  = get_arg(argc, argv, "--tag", "");
  const std::string pt_suffix = pt_bin_suffix(pt_min, pt_max) + (bias_tag.empty() ? "" : "_" + bias_tag);
  const std::string pt_label  = pt_bin_label(pt_min, pt_max);

  const std::string out_name = get_arg(argc, argv, "--out",
                                       "unbias_weights_plots_" + pt_suffix + ".root");
  const int nbins = std::stoi(get_arg(argc, argv, "--nbins", "80"));
  const int eec_bins = std::stoi(get_arg(argc, argv, "--eec-bins", "60"));
  const double eec_max = std::stod(get_arg(argc, argv, "--eec-max", "0.4"));
  const double eec_logmin_val = std::stod(get_arg(argc, argv, "--eec-logmin", "0.004"));

  // ---- Jet-substructure crosscheck configuration ----
  // --jet-R: the original jet clustering radius. Not hard-coded anywhere
  // else in this codebase -- defaults to 0.4 (the README's stated PYTHIA
  // jet radius) only because you asked for that default; pass --jet-R to
  // override. Used both for the broadening normalization and as the Soft
  // Drop angular-exponent scale (R0) for R_g/k_T.
  const double jet_R    = std::stod(get_arg(argc, argv, "--jet-R",    "0.4"));
  const double sd_zcut  = std::stod(get_arg(argc, argv, "--sd-zcut",  "0.1"));
  const double sd_beta  = std::stod(get_arg(argc, argv, "--sd-beta",  "1.0"));
  const double nsub_beta = std::stod(get_arg(argc, argv, "--nsub-beta", "1.0"));
  // Binning for the new crosscheck histograms -- purely cosmetic axis
  // ranges/bin counts (like --nbins/--eec-bins above), all overridable.
  const int    mass_bins = std::stoi(get_arg(argc, argv, "--mass-bins", "60"));
  const double mass_max  = std::stod(get_arg(argc, argv, "--mass-max",  "60.0"));
  const int    rg_bins   = std::stoi(get_arg(argc, argv, "--rg-bins",   "25"));
  const int    kt_bins   = std::stoi(get_arg(argc, argv, "--kt-bins",   "40"));
  const double kt_max    = std::stod(get_arg(argc, argv, "--kt-max",    "40.0"));
  const int    tau_bins  = std::stoi(get_arg(argc, argv, "--tau-bins",  "50"));
  const double tau_max   = std::stod(get_arg(argc, argv, "--tau-max",   "1.0"));
  const int    broad_bins = std::stoi(get_arg(argc, argv, "--broadening-bins", "50"));
  const double broad_max  = std::stod(get_arg(argc, argv, "--broadening-max",  "0.6"));

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
  // unless explicitly overridden with --eec-norm. There is no hard-coded
  // fallback: if the config file/key is missing, read_bin_center() aborts
  // with an error instead of guessing a number.
  const std::string config_file = get_arg(argc, argv, "--config", "unbiasing_config.env");
  const double bin_center = read_bin_center(config_file);
  // get_arg's own default value is never actually used here -- has_arg()
  // already guarantees --eec-norm is present before get_arg() is called
  // with it -- so there is no literal number hiding in this default either.
  const double eec_norm = has_arg(argc, argv, "--eec-norm")
                               ? std::stod(get_arg(argc, argv, "--eec-norm"))
                               : bin_center;
  std::cout << "Bin-center config: '" << config_file << "'  BinCenter=" << bin_center
            << "  (EEC norm used here = " << eec_norm << ")" << std::endl;
  std::cout << "pT bin: " << pt_label << "  (file suffix '" << pt_suffix << "')" << std::endl;
  std::cout << "Jet-substructure crosscheck config: jet-R=" << jet_R
            << "  SoftDrop(zcut=" << sd_zcut << ", beta=" << sd_beta << ")"
            << "  nsub-beta=" << nsub_beta
#ifdef JET_OBS_HAVE_NSUBJETTINESS
            << "  (N-subjettiness: ENABLED)"
#else
            << "  (N-subjettiness: DISABLED -- build with -DJET_OBS_HAVE_NSUBJETTINESS "
               "and link fastjet-contrib's Nsubjettiness to enable)"
#endif
            << std::endl;

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

  // ---- Jet-level eta/phi/mass for the source sample, if unbias_weights.cc
  // forwarded them (see that file's "Jet-level eta/phi/mass passthrough"
  // note). All-or-nothing on its side, but read independently here so a
  // partially-populated file (e.g. hand-edited) still gets whichever of the
  // mass-only / broadening crosschecks it can support.
  float src_mass = 0.0f, src_eta = 0.0f, src_phi = 0.0f;
  const bool has_src_mass = (t->GetBranch("mass") != nullptr);
  const bool has_src_eta  = (t->GetBranch("eta")  != nullptr);
  const bool has_src_phi  = (t->GetBranch("phi")  != nullptr);
  if (has_src_mass) t->SetBranchAddress("mass", &src_mass);
  if (has_src_eta)  t->SetBranchAddress("eta",  &src_eta);
  if (has_src_phi)  t->SetBranchAddress("phi",  &src_phi);
  const bool has_src_axis = has_src_eta && has_src_phi;
  if (!has_src_mass)
    std::cout << "note: source tree has no 'mass' branch -- jet-mass crosscheck skipped for the source/weighted samples." << std::endl;
  if (!has_src_axis)
    std::cout << "note: source tree has no 'eta'/'phi' branch(es) -- broadening crosscheck skipped for the source/weighted samples." << std::endl;

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

  // ---- Jet-substructure crosscheck histograms (base/total/target) ----
  TH1D *h_mass_base   = new TH1D("h_mass_base",   "jet mass; M_{jet} [GeV]; weighted entries", mass_bins, 0.0, mass_max);
  TH1D *h_mass_total  = new TH1D("h_mass_total",  "jet mass; M_{jet} [GeV]; weighted entries", mass_bins, 0.0, mass_max);
  TH1D *h_mass_target = new TH1D("h_mass_target", "jet mass; M_{jet} [GeV]; weighted entries", mass_bins, 0.0, mass_max);

  TH1D *h_rg_base   = new TH1D("h_rg_base",   "groomed radius; R_{g}; weighted entries", rg_bins, 0.0, jet_R);
  TH1D *h_rg_total  = new TH1D("h_rg_total",  "groomed radius; R_{g}; weighted entries", rg_bins, 0.0, jet_R);
  TH1D *h_rg_target = new TH1D("h_rg_target", "groomed radius; R_{g}; weighted entries", rg_bins, 0.0, jet_R);

  TH1D *h_kt_base   = new TH1D("h_kt_base",   "groomed k_{T}; k_{T} [GeV]; weighted entries", kt_bins, 0.0, kt_max);
  TH1D *h_kt_total  = new TH1D("h_kt_total",  "groomed k_{T}; k_{T} [GeV]; weighted entries", kt_bins, 0.0, kt_max);
  TH1D *h_kt_target = new TH1D("h_kt_target", "groomed k_{T}; k_{T} [GeV]; weighted entries", kt_bins, 0.0, kt_max);

  TH1D *h_broad_base   = new TH1D("h_broad_base",   "broadening; #lambda_{1}^{(#kappa=1)}; weighted entries", broad_bins, 0.0, broad_max);
  TH1D *h_broad_total  = new TH1D("h_broad_total",  "broadening; #lambda_{1}^{(#kappa=1)}; weighted entries", broad_bins, 0.0, broad_max);
  TH1D *h_broad_target = new TH1D("h_broad_target", "broadening; #lambda_{1}^{(#kappa=1)}; weighted entries", broad_bins, 0.0, broad_max);

  for (TH1D *h : { h_mass_base, h_mass_total, h_mass_target,
                   h_rg_base, h_rg_total, h_rg_target,
                   h_kt_base, h_kt_total, h_kt_target,
                   h_broad_base, h_broad_total, h_broad_target })
    h->Sumw2();

#ifdef JET_OBS_HAVE_NSUBJETTINESS
  TH1D *h_tau1_base   = new TH1D("h_tau1_base",   "#tau_{1}; #tau_{1}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau1_total  = new TH1D("h_tau1_total",  "#tau_{1}; #tau_{1}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau1_target = new TH1D("h_tau1_target", "#tau_{1}; #tau_{1}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau2_base   = new TH1D("h_tau2_base",   "#tau_{2}; #tau_{2}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau2_total  = new TH1D("h_tau2_total",  "#tau_{2}; #tau_{2}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau2_target = new TH1D("h_tau2_target", "#tau_{2}; #tau_{2}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau3_base   = new TH1D("h_tau3_base",   "#tau_{3}; #tau_{3}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau3_total  = new TH1D("h_tau3_total",  "#tau_{3}; #tau_{3}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau3_target = new TH1D("h_tau3_target", "#tau_{3}; #tau_{3}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau21_base   = new TH1D("h_tau21_base",   "#tau_{21}; #tau_{21}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau21_total  = new TH1D("h_tau21_total",  "#tau_{21}; #tau_{21}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau21_target = new TH1D("h_tau21_target", "#tau_{21}; #tau_{21}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau32_base   = new TH1D("h_tau32_base",   "#tau_{32}; #tau_{32}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau32_total  = new TH1D("h_tau32_total",  "#tau_{32}; #tau_{32}; weighted entries", tau_bins, 0.0, tau_max);
  TH1D *h_tau32_target = new TH1D("h_tau32_target", "#tau_{32}; #tau_{32}; weighted entries", tau_bins, 0.0, tau_max);
  for (TH1D *h : { h_tau1_base, h_tau1_total, h_tau1_target,
                   h_tau2_base, h_tau2_total, h_tau2_target,
                   h_tau3_base, h_tau3_total, h_tau3_target,
                   h_tau21_base, h_tau21_total, h_tau21_target,
                   h_tau32_base, h_tau32_total, h_tau32_target })
    h->Sumw2();
#endif

  const Long64_t nentries = t->GetEntries();
  double sumw_base = 0.0;
  double sumw_total = 0.0;
  long n_sd_src_total = 0, n_sd_src_groomed_away = 0;
  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    h_base->Fill(pt, w_base);
    h_total->Fill(pt, w_total);
    h_w_unbias->Fill(w_unbias);
    const bool in_window = (pt >= pt_min && pt <= pt_max);

    if (in_window && has_src_mass) {
      h_mass_base->Fill(src_mass, w_base);
      h_mass_total->Fill(src_mass, w_total);
    }

    if (in_window && has_const && const_pt && const_eta && const_phi &&
        const_pt->size() == const_eta->size() && const_pt->size() == const_phi->size()) {
      fill_eec(h_eec_base, *const_pt, *const_eta, *const_phi, w_base, eec_norm);
      fill_eec(h_eec_total, *const_pt, *const_eta, *const_phi, w_total, eec_norm);
      sumw_base += w_base;
      sumw_total += w_total;

      const jetobs::SoftDropResult sd = jetobs::soft_drop(*const_pt, *const_eta, *const_phi,
                                                          jet_R, sd_zcut, sd_beta);
      ++n_sd_src_total;
      if (sd.groomed_away) {
        ++n_sd_src_groomed_away;
      } else {
        h_rg_base->Fill(sd.Rg, w_base);   h_rg_total->Fill(sd.Rg, w_total);
        h_kt_base->Fill(sd.kT, w_base);   h_kt_total->Fill(sd.kT, w_total);
      }

      if (has_src_axis) {
        const double lam = jetobs::broadening(*const_pt, *const_eta, *const_phi,
                                              src_eta, src_phi, jet_R);
        h_broad_base->Fill(lam, w_base);
        h_broad_total->Fill(lam, w_total);
      }

#ifdef JET_OBS_HAVE_NSUBJETTINESS
      const jetobs::NSubResult ns = jetobs::nsubjettiness(*const_pt, *const_eta, *const_phi,
                                                          jet_R, nsub_beta);
      h_tau1_base->Fill(ns.tau1, w_base); h_tau1_total->Fill(ns.tau1, w_total);
      h_tau2_base->Fill(ns.tau2, w_base); h_tau2_total->Fill(ns.tau2, w_total);
      h_tau3_base->Fill(ns.tau3, w_base); h_tau3_total->Fill(ns.tau3, w_total);
      if (ns.tau1 > 0.0) {
        h_tau21_base->Fill(ns.tau2 / ns.tau1, w_base);
        h_tau21_total->Fill(ns.tau2 / ns.tau1, w_total);
      }
      if (ns.tau2 > 0.0) {
        h_tau32_base->Fill(ns.tau3 / ns.tau2, w_base);
        h_tau32_total->Fill(ns.tau3 / ns.tau2, w_total);
      }
#endif
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

  // Jet-level mass/eta/phi read directly from the SAME target tree, exactly
  // like const_pt/eta/phi above -- startBasis.C's output trees already carry
  // scalar "mass"/"eta"/"phi" branches, so no change to unbias_weights.cc's
  // output is needed on the target side.
  float t_mass = 0.0f, t_eta = 0.0f, t_phi = 0.0f;
  const bool target_has_mass = (tt->GetBranch("mass") != nullptr);
  const bool target_has_eta  = (tt->GetBranch("eta")  != nullptr);
  const bool target_has_phi  = (tt->GetBranch("phi")  != nullptr);
  if (target_has_mass) tt->SetBranchAddress("mass", &t_mass);
  if (target_has_eta)  tt->SetBranchAddress("eta",  &t_eta);
  if (target_has_phi)  tt->SetBranchAddress("phi",  &t_phi);
  const bool target_has_axis = target_has_eta && target_has_phi;
  if (!target_has_mass)
    std::cout << "note: target tree has no 'mass' branch -- jet-mass crosscheck skipped for the target sample." << std::endl;
  if (!target_has_axis)
    std::cout << "note: target tree has no 'eta'/'phi' branch(es) -- broadening crosscheck skipped for the target sample." << std::endl;

  const Long64_t tentries = tt->GetEntries();
  double sumw_target = 0.0;
  long n_sd_tgt_total = 0, n_sd_tgt_groomed_away = 0;
  for (Long64_t i = 0; i < tentries; ++i) {
    tt->GetEntry(i);
    h_target->Fill(tpt, tw);
    const bool t_in_window = (tpt >= pt_min && tpt <= pt_max);

    if (t_in_window && target_has_mass) h_mass_target->Fill(t_mass, tw);

    if (t_in_window && target_has_const && tconst_pt && tconst_eta && tconst_phi &&
        tconst_pt->size() == tconst_eta->size() && tconst_pt->size() == tconst_phi->size()) {
      fill_eec(h_eec_target, *tconst_pt, *tconst_eta, *tconst_phi, tw, eec_norm);
      sumw_target += tw;

      const jetobs::SoftDropResult sd = jetobs::soft_drop(*tconst_pt, *tconst_eta, *tconst_phi,
                                                          jet_R, sd_zcut, sd_beta);
      ++n_sd_tgt_total;
      if (sd.groomed_away) {
        ++n_sd_tgt_groomed_away;
      } else {
        h_rg_target->Fill(sd.Rg, tw);
        h_kt_target->Fill(sd.kT, tw);
      }

      if (target_has_axis) {
        const double lam = jetobs::broadening(*tconst_pt, *tconst_eta, *tconst_phi,
                                              t_eta, t_phi, jet_R);
        h_broad_target->Fill(lam, tw);
      }

#ifdef JET_OBS_HAVE_NSUBJETTINESS
      const jetobs::NSubResult ns = jetobs::nsubjettiness(*tconst_pt, *tconst_eta, *tconst_phi,
                                                          jet_R, nsub_beta);
      h_tau1_target->Fill(ns.tau1, tw);
      h_tau2_target->Fill(ns.tau2, tw);
      h_tau3_target->Fill(ns.tau3, tw);
      if (ns.tau1 > 0.0) h_tau21_target->Fill(ns.tau2 / ns.tau1, tw);
      if (ns.tau2 > 0.0) h_tau32_target->Fill(ns.tau3 / ns.tau2, tw);
#endif
    }
  }

  std::cout << "Soft Drop groomed-away fraction: source "
            << n_sd_src_groomed_away << "/" << n_sd_src_total
            << "  target " << n_sd_tgt_groomed_away << "/" << n_sd_tgt_total
            << "  (these jets are excluded from the R_g/k_T histograms, not binned at 0)"
            << std::endl;

  style_hist(h_base,   kBlue + 1,  20, cfg);
  style_hist(h_total,  kRed + 1,   21, cfg);
  style_hist(h_target, kGreen + 2, 22, cfg);
  style_hist(h_eec_base,   kBlue + 1,  20, cfg);
  style_hist(h_eec_total,  kRed + 1,   21, cfg);
  style_hist(h_eec_target, kGreen + 2, 22, cfg);

  for (auto trip : { std::array<TH1D*,3>{h_mass_base, h_mass_total, h_mass_target},
                     std::array<TH1D*,3>{h_rg_base, h_rg_total, h_rg_target},
                     std::array<TH1D*,3>{h_kt_base, h_kt_total, h_kt_target},
                     std::array<TH1D*,3>{h_broad_base, h_broad_total, h_broad_target} }) {
    style_hist(trip[0], kBlue + 1,  20, cfg);
    style_hist(trip[1], kRed + 1,   21, cfg);
    style_hist(trip[2], kGreen + 2, 22, cfg);
  }
#ifdef JET_OBS_HAVE_NSUBJETTINESS
  for (auto trip : { std::array<TH1D*,3>{h_tau1_base, h_tau1_total, h_tau1_target},
                     std::array<TH1D*,3>{h_tau2_base, h_tau2_total, h_tau2_target},
                     std::array<TH1D*,3>{h_tau3_base, h_tau3_total, h_tau3_target},
                     std::array<TH1D*,3>{h_tau21_base, h_tau21_total, h_tau21_target},
                     std::array<TH1D*,3>{h_tau32_base, h_tau32_total, h_tau32_target} }) {
    style_hist(trip[0], kBlue + 1,  20, cfg);
    style_hist(trip[1], kRed + 1,   21, cfg);
    style_hist(trip[2], kGreen + 2, 22, cfg);
  }
#endif

  // =====================================================================
  // pT check plot (main + ratio-to-target panels)
  // =====================================================================
  TCanvas *c = new TCanvas("c_unbias", "unbias check", 900, 900);
  TPad *p_top = nullptr, *p_bot = nullptr;
  make_ratio_pads("p_pt", p_top, p_bot, cfg);

  p_top->cd();
  gPad->SetLogy();
  set_top_headroom(h_target, std::max({h_base->GetMaximum(), h_total->GetMaximum(),
                                       h_target->GetMaximum()}), /*logy=*/true);
  hide_shared_xaxis(h_target);
  h_target->Draw("E1");   // draw target first/underneath so it sets the axis range
  h_base->Draw("E1 same");
  h_total->Draw("E1 same");

  TLegend *leg = make_legend(0.42, 0.68, 0.82, 0.90, cfg);
  leg->AddEntry(h_base, "biased sample (base weight)", "lp");
  leg->AddEntry(h_total, "weighted sample (base x unbias)", "lp");
  leg->AddEntry(h_target, Form("unbiased target %s", target_tree.c_str()), "lp");
  leg->Draw();
  draw_bin_label(pt_label, cfg);

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
  h_mass_base->Write(); h_mass_total->Write(); h_mass_target->Write();
  h_rg_base->Write(); h_rg_total->Write(); h_rg_target->Write();
  h_kt_base->Write(); h_kt_total->Write(); h_kt_target->Write();
  h_broad_base->Write(); h_broad_total->Write(); h_broad_target->Write();
#ifdef JET_OBS_HAVE_NSUBJETTINESS
  h_tau1_base->Write(); h_tau1_total->Write(); h_tau1_target->Write();
  h_tau2_base->Write(); h_tau2_total->Write(); h_tau2_target->Write();
  h_tau3_base->Write(); h_tau3_total->Write(); h_tau3_target->Write();
  h_tau21_base->Write(); h_tau21_total->Write(); h_tau21_target->Write();
  h_tau32_base->Write(); h_tau32_total->Write(); h_tau32_target->Write();
#endif
  c->Write();

  save_canvas(c, "plot_unbias_pT_check_" + pt_suffix, cfg);

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
  draw_bin_label(pt_label, cfg);
  save_canvas(cw, "plot_unbias_weights_weights_" + pt_suffix, cfg);
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
  set_top_headroom(h_eec_target, std::max({h_eec_base->GetMaximum(), h_eec_total->GetMaximum(),
                                           h_eec_target->GetMaximum()}), /*logy=*/true);
  hide_shared_xaxis(h_eec_target);
  h_eec_target->Draw("E1");
  h_eec_base->Draw("E1 same");
  h_eec_total->Draw("E1 same");
  TLegend *lege = make_legend(0.42, 0.20, 0.82, 0.42, cfg);
  lege->AddEntry(h_eec_base, "biased sample (base weight)", "lp");
  lege->AddEntry(h_eec_total, "weighted sample (base x unbias)", "lp");
  lege->AddEntry(h_eec_target, Form("unbiased target %s", target_tree.c_str()), "lp");
  lege->Draw();
  draw_bin_label(pt_label, cfg);

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

  save_canvas(ce, "plot_eec_compare_" + pt_suffix, cfg);
  ce->Write();

  // =====================================================================
  // Jet-substructure crosscheck plots (main + ratio-to-target panels),
  // each via the shared observable_check_plot() helper.
  // =====================================================================
  if (has_src_mass && target_has_mass) {
    observable_check_plot("jet_mass", "M_{jet} [GeV]", h_mass_base, h_mass_total, h_mass_target,
                          /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  } else {
    std::cout << "note: jet-mass crosscheck plot skipped (missing 'mass' branch on source and/or target)." << std::endl;
  }

  observable_check_plot("rg", "R_{g}", h_rg_base, h_rg_total, h_rg_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  observable_check_plot("kt", "k_{T} [GeV]", h_kt_base, h_kt_total, h_kt_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);

  if (has_src_axis && target_has_axis) {
    observable_check_plot("broadening", "#lambda_{1}^{(#kappa=1)}", h_broad_base, h_broad_total, h_broad_target,
                          /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  } else {
    std::cout << "note: broadening crosscheck plot skipped (missing 'eta'/'phi' branch(es) on source and/or target)." << std::endl;
  }

#ifdef JET_OBS_HAVE_NSUBJETTINESS
  observable_check_plot("tau1", "#tau_{1}", h_tau1_base, h_tau1_total, h_tau1_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  observable_check_plot("tau2", "#tau_{2}", h_tau2_base, h_tau2_total, h_tau2_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  observable_check_plot("tau3", "#tau_{3}", h_tau3_base, h_tau3_total, h_tau3_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  observable_check_plot("tau21", "#tau_{21}", h_tau21_base, h_tau21_total, h_tau21_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
  observable_check_plot("tau32", "#tau_{32}", h_tau32_base, h_tau32_total, h_tau32_target,
                        /*logy=*/true, /*logx=*/false, target_tree, pt_label, pt_suffix, cfg);
#else
  std::cout << "note: N-subjettiness (tau1/tau2/tau3/tau21/tau32) plots skipped -- built without "
               "-DJET_OBS_HAVE_NSUBJETTINESS." << std::endl;
#endif

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
            set_top_headroom(htar[p], std::max({hbase[p]->GetMaximum(), htot[p]->GetMaximum(),
                                                htar[p]->GetMaximum()}), /*logy=*/true);
            hide_shared_xaxis(htar[p]);
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
          // One overall pT-bin label for the whole grid (not per-panel, to
          // avoid cluttering every individual basis-function panel).
          ctw->cd();
          draw_bin_label(pt_label, cfg, 0.02, 0.985);
          save_canvas(ctw, "plot_theta_basis_weighted_compare_" + pt_suffix, cfg);
        }
      }
    }
    inb.Close();
  }

  std::cout << "wrote " << out_name << " and (in " << cfg.plot_dir << "/) plot_unbias_pT_check_" << pt_suffix << " / "
            << "plot_unbias_weights_weights_" << pt_suffix << " / plot_eec_compare_" << pt_suffix << " / "
            << "plot_jet_mass_" << pt_suffix << " / plot_rg_" << pt_suffix << " / plot_kt_" << pt_suffix
            << " / plot_broadening_" << pt_suffix
#ifdef JET_OBS_HAVE_NSUBJETTINESS
            << " / plot_tau{1,2,3,21,32}_" << pt_suffix
#endif
            << " / plot_theta_basis_weighted_compare_" << pt_suffix << " in: ";
  for (const std::string &fmt : cfg.formats) std::cout << "." << fmt << " ";
  std::cout << std::endl;
  return 0;
}