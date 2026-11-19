#!/usr/bin/env python3
"""Matplotlib versions of the ROOT plots in Output/plots/.

Same content, axis labels, legends, log scales and ranges as the ROOT plots, drawn
as step histograms with LaTeX text and a colour-blind-safe palette (Okabe-Ito).

Output layout: <plot-dir>/<bias>/<sample>/pt<lo>-<hi>/<plot>.pdf
  bias   = powerlaw or lognormal
  sample = the sample tag the inputs were made with (UNBIAS_TAG, e.g. 10M); untagged files are skipped
e.g. lognormal/10M/pt100-140/closure_pt.pdf. Made from the ROOT files the chain writes:

  plot                       ROOT counterpart                     input
  pt_spectrum_bias           jetPtShift_TChainReader              startBasis_*.root
  closure_pt                 plot_unbias_pT_check                 unbias_weights_plots_*.root
  unbias_weights             plot_unbias_weights_weights            (histograms saved by
  closure_eec                plot_eec_compare                        plot_unbias_weights)
  closure_jet_mass, _rg, _kt, _broadening, _tau*
                             plot_jet_mass, plot_rg, ...
  basis_functions/g<k>_<description>.pdf                          unbias_weights_*.root
                             plot_theta_basis_weighted_compare
                             (one closure plot per basis function)

Subcommands:
  list        which bias / sample combinations exist in --files-dir, and what they contain
  all         plot them; choose with --bias and --sample (default: everything)
  startbasis  pt_spectrum_bias from one startBasis_*.root file
  unbias      the other plots from one unbias_weights_plots_*.root (+ unbias_weights_*.root)

Usage (from AnalysisCode/):
  python plot_mpl.py list
  python plot_mpl.py all --bias lognormal --sample 10M
  python plot_mpl.py all                                   # every bias and sample
  python plot_mpl.py startbasis ../Output/files/startBasis_100_140_10M.root --plot-dir ../Output/plots/mpl
  python plot_mpl.py unbias --hists ../Output/files/unbias_weights_plots_100_140_lognormal.root \
      --weights ../Output/files/unbias_weights_100_140_lognormal.root \
      --pt-min 100 --pt-max 140 --tag lognormal --plot-dir ../Output/plots/mpl

Text is typeset with LaTeX (needs a TeX installation: latex, dvipng); without one, or
with --no-tex, matplotlib's built-in Computer Modern mathtext is used instead.
The chain calls this automatically when UNBIAS_PLOTTER=matplotlib or both (see SETUP.md).
"""
import argparse
import math
import os
import re
import shutil
import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import uproot

# Okabe-Ito colour-blind-safe palette
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#000000"]
C_BASE, C_TOTAL, C_TARGET = "#0072B2", "#D55E00", "#009E73"   # blue, vermillion, green
C_REF, C_BIASED, C_PP = "#009E73", "#CC79A7", "#999999"       # green, reddish purple, grey

LEG_BASE = "Biased"
LEG_TOTAL = "Reweighted"
LEG_TARGET = "Unbiased"


def setup_style(usetex):
    plt.rcParams.update({
        "axes.prop_cycle": matplotlib.cycler(color=OKABE_ITO),
        "font.family": "serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "legend.fontsize": 9.5,
        "legend.frameon": False,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "pdf.fonttype": 42,
    })
    if usetex:
        plt.rcParams.update({
            "text.usetex": True,
            "font.serif": ["Computer Modern Roman"],
            "text.latex.preamble": r"\usepackage{amsmath}",
        })
    else:
        plt.rcParams.update({"mathtext.fontset": "cm", "font.serif": ["cmr10"],
                             "axes.formatter.use_mathtext": True})


def tex_escape(s):
    """Plain text -> LaTeX-safe text (only needed for names coming from files)."""
    return re.sub(r"([_%&#$])", r"\\\1", s)


def fmt_pt(x):
    """Same as format_pt_for_name() in plot_unbias_weights.cc: 100, 140, 97.5, ..."""
    return f"{x:g}"


def bin_label(pt_min, pt_max):
    return rf"${fmt_pt(pt_min)} < p_{{\mathrm{{T}}}} < {fmt_pt(pt_max)}$ GeV"


def pt_subdir(pt_min, pt_max):
    """Per-bin plot folder, e.g. pt100-140."""
    return f"pt{fmt_pt(pt_min)}-{fmt_pt(pt_max)}"


def parse_tag(tag):
    """What follows the pT bin in the ROOT file names ("", "10M", "lognormal", "lognormal_10M")
    -> (bias, sample). No "lognormal" means the default power-law bias; sample is "" if untagged."""
    parts = [p for p in tag.split("_") if p]
    bias = "lognormal" if "lognormal" in parts else "powerlaw"
    sample = "_".join(p for p in parts if p != "lognormal")
    return bias, sample


def out_dir(plot_dir, tag, pt_min, pt_max):
    """<plot_dir>/<bias>/<sample>/pt<lo>-<hi>"""
    bias, sample = parse_tag(tag)
    return os.path.join(plot_dir, bias, sample, pt_subdir(pt_min, pt_max))


def capitalize(text):
    return text[:1].upper() + text[1:]


def save(fig, plot_dir, name, quiet=False):
    os.makedirs(plot_dir, exist_ok=True)
    path = os.path.join(plot_dir, name + ".pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    if not quiet:
        print(f"plot_mpl: wrote {path}")


# ---------------------------------------------------------------------------
# Histogram helpers
# ---------------------------------------------------------------------------
class Hist:
    """Minimal 1D histogram: edges, values, errors (no under/overflow)."""

    def __init__(self, edges, values, errors):
        self.edges = np.asarray(edges, float)
        self.values = np.asarray(values, float)
        self.errors = np.asarray(errors, float)

    @classmethod
    def from_root(cls, h):
        return cls(h.axis().edges(), h.values(), h.errors())

    @classmethod
    def fill(cls, data, weights, edges):
        v, _ = np.histogram(data, bins=edges, weights=weights)
        w2, _ = np.histogram(data, bins=edges, weights=np.asarray(weights, float) ** 2)
        return cls(edges, v, np.sqrt(w2))

    def centers(self, logx=False):
        if logx:
            return np.sqrt(self.edges[1:] * self.edges[:-1])
        return 0.5 * (self.edges[1:] + self.edges[:-1])

    def divide(self, other):
        """Ratio with uncorrelated error propagation (as TH1::Divide); NaN where undefined."""
        a, ea, b, eb = self.values, self.errors, other.values, other.errors
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(b != 0, a / b, np.nan)
            e = np.where(b != 0, np.sqrt((ea * b) ** 2 + (eb * a) ** 2) / b ** 2, np.nan)
        return Hist(self.edges, r, e)


def data_xlim(hists):
    """x range spanned by the non-empty bins of any of the histograms (full range if all empty)."""
    edges = hists[0].edges
    filled = np.zeros(len(edges) - 1, bool)
    for h in hists:
        filled |= np.isfinite(h.values) & (h.values > 0)
    idx = np.flatnonzero(filled)
    if not idx.size:
        return edges[0], edges[-1]
    return edges[idx[0]], edges[idx[-1] + 1]


def draw_band(ax, h, color, alpha=0.2):
    """Shaded +-1 sigma band around a histogram (bins with undefined values are skipped)."""
    lo, hi = h.values - h.errors, h.values + h.errors
    ok = np.isfinite(lo) & np.isfinite(hi)
    ax.stairs(np.where(ok, hi, np.nan), h.edges, baseline=np.where(ok, lo, np.nan),
              fill=True, color=color, alpha=alpha, lw=0)


def draw_hist(ax, h, color, label=None, logx=False, errors=True, band=False):
    """Step histogram with vertical statistical error bars (and optionally a shaded band)."""
    if band:
        draw_band(ax, h, color)
    ax.stairs(h.values, h.edges, color=color, lw=1.4, label=label)
    if errors:
        sel = np.isfinite(h.values) & (h.values > 0)
        ax.errorbar(h.centers(logx)[sel], h.values[sel], yerr=h.errors[sel],
                    fmt="none", ecolor=color, elinewidth=1.0, capsize=0)


def draw_ratio(ax, r, color, errorbars=False, logx=False):
    """Ratio as a step line with a shaded +-1 sigma band (and optionally error bars)."""
    ax.stairs(r.values, r.edges, color=color, lw=1.4)
    draw_band(ax, r, color, alpha=0.25)
    if errorbars:
        ok = np.isfinite(r.values)
        ax.errorbar(r.centers(logx)[ok], r.values[ok], yerr=r.errors[ok],
                    fmt="none", ecolor=color, elinewidth=1.0, capsize=0)


def legend_handles(target_name):
    from matplotlib.lines import Line2D
    return ([Line2D([], [], color=c, lw=1.4) for c in (C_BASE, C_TOTAL, C_TARGET)],
            [LEG_BASE, LEG_TOTAL, LEG_TARGET])


# ---------------------------------------------------------------------------
# Main panel + ratio-to-target panel
# ---------------------------------------------------------------------------
def draw_ratio_plot(ax_top, ax_bot, h_base, h_total, h_target, xtitle, ytitle,
                    logx, logy, target_name, label=None, legend_loc="upper right", bands_only=False):
    """Main panel + ratio panel. Default: error bars in the main panel, shaded bands in the
    ratio panel. bands_only=True: shaded +-1 sigma bands in both panels, no error bars."""
    for h, col in ((h_target, C_TARGET), (h_base, C_BASE), (h_total, C_TOTAL)):
        draw_hist(ax_top, h, col, logx=logx, errors=not bands_only, band=bands_only)
    if logx:
        ax_top.set_xscale("log")
    if logy:
        ax_top.set_yscale("log")
        # y range from the bin contents (not the error bars, which can reach ~0 on a
        # log axis), with one decade of headroom for the legend, as in the ROOT plots
        vals = np.concatenate([h.values for h in (h_base, h_total, h_target)])
        pos = vals[np.isfinite(vals) & (vals > 0)]
        if pos.size:
            ax_top.set_ylim(0.5 * pos.min(), 10 * pos.max())
    ax_top.set_ylabel(ytitle)
    if legend_loc:
        ax_top.legend(*legend_handles(target_name), loc=legend_loc)
    if label:
        ax_top.text(0.04, 0.94, label, transform=ax_top.transAxes, va="top")

    draw_ratio(ax_bot, h_total.divide(h_target), C_TOTAL)
    draw_ratio(ax_bot, h_base.divide(h_target), C_BASE)
    ax_bot.axhline(1.0, color="black", ls="--", lw=0.8)
    ax_bot.set_ylim(0.5, 1.5)
    ax_bot.set_xlim(*data_xlim([h_base, h_total, h_target]))
    if logx:   # no labels on minor log ticks (they overlap on ranges under a decade)
        ax_bot.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax_bot.set_xlabel(xtitle)
    ax_bot.set_ylabel("Ratio")


def ratio_figure():
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, height_ratios=[3, 1],
                                         figsize=(5, 5.2))
    fig.subplots_adjust(hspace=0.05)
    return fig, ax_top, ax_bot


# ---------------------------------------------------------------------------
# startBasis.C: jetPtShift_TChainReader plot
# ---------------------------------------------------------------------------
def plot_startbasis(path, plot_dir):
    m = re.fullmatch(r"startBasis_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)(.*)\.root", os.path.basename(path))
    if m:
        plot_dir = out_dir(plot_dir, m.group(3), float(m.group(1)), float(m.group(2)))
    # binning as in startBasis.C: pT window +- 10 GeV, 0.6 GeV bins
    if m:
        lo_w, hi_w = float(m.group(1)) - 10.0, float(m.group(2)) + 10.0
        edges = np.linspace(lo_w, hi_w, int(round((hi_w - lo_w) / 0.6)) + 1)
    with uproot.open(path) as f:
        if not m:
            edges = f["hPtX"].axis().edges()
        tr = {n: f[n].arrays(["pt_raw", "weight"], library="np") for n in ("tX", "tY", "tRef", "tBiased")}
    h_ref = Hist.fill(tr["tRef"]["pt_raw"], tr["tRef"]["weight"], edges)
    h_bia = Hist.fill(tr["tBiased"]["pt_raw"], tr["tBiased"]["weight"], edges)
    h_pp = Hist.fill(np.concatenate([tr["tX"]["pt_raw"], tr["tY"]["pt_raw"]]),            # X + Y (hpTpp)
                     np.concatenate([tr["tX"]["weight"], tr["tY"]["weight"]]), edges)

    fig, ax = plt.subplots(figsize=(5, 4.6))
    ax.stairs(h_pp.values, edges, color=C_PP, lw=0.8, zorder=1)   # X + Y (hpTpp), unlabelled as in ROOT
    draw_hist(ax, h_ref, C_REF, label="Unbiased")
    draw_hist(ax, h_bia, C_BIASED, label="Biased")
    ax.set_yscale("log")
    pos = h_ref.values[h_ref.values > 0]
    if pos.size:
        ax.set_ylim(0.5 * pos.min(), 2 * pos.max())
    ax.set_xlim(*data_xlim([h_ref, h_bia]))
    ax.set_xlabel(r"$p_{\mathrm{T}}$ (GeV)")
    ax.set_ylabel(r"$\mathrm{d}N/\mathrm{d}p_{\mathrm{T}}$")
    ax.legend()
    save(fig, plot_dir, "pt_spectrum_bias")


# ---------------------------------------------------------------------------
# plot_unbias_weights.cc plots
# ---------------------------------------------------------------------------
# (name in ROOT file, plot name, x label)
OBSERVABLES = [
    ("mass", "jet_mass", r"$M_{\mathrm{jet}}$ [GeV]"),
    ("rg", "rg", r"$R_{\mathrm{g}}$"),
    ("kt", "kt", r"$k_{\mathrm{T}}$ [GeV]"),
    ("broad", "broadening", r"$\lambda_{1}$"),
    ("tau1", "tau1", r"$\tau_{1}$"),
    ("tau2", "tau2", r"$\tau_{2}$"),
    ("tau3", "tau3", r"$\tau_{3}$"),
    ("tau21", "tau21", r"$\tau_{21}$"),
    ("tau32", "tau32", r"$\tau_{32}$"),
]


def y_title(h):
    try:
        return capitalize(tex_escape(h.member("fYaxis").member("fTitle").strip())) or "Entries"
    except Exception:
        return "Entries"


def plot_unbias(hists_path, weights_path, pt_min, pt_max, tag, target_tree, plot_dir, basis_opts):
    plot_dir = out_dir(plot_dir, tag, pt_min, pt_max)
    label = bin_label(pt_min, pt_max)

    with uproot.open(hists_path) as f:
        keys = {k.split(";")[0] for k in f.keys()}

        def trip(stem):
            return [Hist.from_root(f[f"h{stem}_{s}"]) for s in ("base", "total", "target")]

        # pT check
        if {"h_base", "h_total", "h_target"} <= keys:
            fig, a, b = ratio_figure()
            draw_ratio_plot(a, b, *trip(""), r"$p_{\mathrm{T}}$ [GeV]", y_title(f["h_target"]),
                            logx=False, logy=True, target_name=target_tree, label=label)
            save(fig, plot_dir, "closure_pt")

        # weight distribution
        if "h_w_unbias" in keys:
            hw = Hist.from_root(f["h_w_unbias"])
            fig, ax = plt.subplots(figsize=(5, 4.6))
            draw_hist(ax, hw, "black")
            ax.set_yscale("log")
            pos = hw.values[hw.values > 0]
            if pos.size:
                ax.set_ylim(0.5 * pos.min(), 2 * pos.max())
            ax.set_xlim(*data_xlim([hw]))
            ax.set_xlabel(r"Unbiasing weight $w$")
            ax.set_ylabel("Entries")
            ax.text(0.96, 0.94, label, transform=ax.transAxes, va="top", ha="right")   # data peaks on the left
            save(fig, plot_dir, "unbias_weights")

        # EEC (log-log)
        if {"h_eec_base", "h_eec_total", "h_eec_target"} <= keys:
            fig, a, b = ratio_figure()
            draw_ratio_plot(a, b, *trip("_eec"), r"$\theta$", "EEC",
                            logx=True, logy=True, target_name=target_tree, label=label,
                            legend_loc="lower right")
            save(fig, plot_dir, "closure_eec")

        # jet-substructure crosschecks
        for stem, name, xtitle in OBSERVABLES:
            if not {f"h_{stem}_base", f"h_{stem}_total", f"h_{stem}_target"} <= keys:
                continue
            fig, a, b = ratio_figure()
            draw_ratio_plot(a, b, *trip(f"_{stem}"), xtitle, y_title(f[f"h_{stem}_target"]),
                            logx=False, logy=True, target_name=target_tree, label=label)
            save(fig, plot_dir, f"closure_{name}")

    if weights_path:
        plot_basis(weights_path, pt_min, pt_max, label, target_tree, plot_dir, **basis_opts)


# ---------------------------------------------------------------------------
# Basis functions: labels -> LaTeX formulas, grouped into families
# ---------------------------------------------------------------------------
RE_EEC = re.compile(r"eec: dR\^\{-(\S+?)\} ln\^\{(\S+?)\} z\^\{(\S+?)\} \[dR<(\S+?)\]")
RE_SJ = re.compile(r"sjmom: sum \(pT/(\S+?)\)\^\{(\S+?)\} \(stored subjets R=(\S+?)\)")


def _num(x):
    return f"{x:g}"


def basis_tex(label, norm):
    """basis label (basis_functions.h) -> (formula, description, file-name stem)."""
    m = RE_EEC.fullmatch(label)
    if m:
        a = -float(m.group(1))            # label stores dR^{-A}; exponent is -A
        b, mm, cut = float(m.group(2)), int(float(m.group(3))), float(m.group(4))
        dr = "" if a == 0 else (r"\Delta R_{ik}" if a == 1 else rf"\Delta R_{{ik}}^{{{_num(a)}}}")
        ln = "" if b == 0 else (r"\ln\Delta R_{ik}" if b == 1 else rf"\ln^{{{_num(b)}}}\!\Delta R_{{ik}}")
        z = "" if mm == 0 else (r"z_{ik}" if mm == 1 else rf"z_{{ik}}^{{{mm}}}")
        body = r"\,".join(t for t in (dr, ln, z) if t) or "1"
        formula = rf"$\sum_{{i<k,\ \Delta R_{{ik}}<{_num(cut)}}} {body}$"
        zdef = (rf"$z_{{ik}} = p_{{\mathrm{{T}},i}}\,p_{{\mathrm{{T}},k}}/({_num(norm)}\ \mathrm{{GeV}})^2$" if norm
                else r"$z_{ik} = p_{\mathrm{T},i}\,p_{\mathrm{T},k}/p_{\mathrm{T,0}}^2$")
        desc = "EEC moment, " + zdef
        stem = f"eec_dR{_num(a)}_ln{_num(b)}_z{mm}_cut{_num(cut)}"
        return formula, desc, stem
    m = RE_SJ.fullmatch(label)
    if m:
        nrm, n, R = float(m.group(1)), float(m.group(2)), float(m.group(3))
        formula = rf"$\sum_{{j}} \left(p_{{\mathrm{{T}},j}}/{_num(nrm)}\ \mathrm{{GeV}}\right)^{{{_num(n)}}}$"
        desc = rf"Subjet $p_{{\mathrm{{T}}}}$ moment, C/A subjets $R = {R:.2f}$"
        stem = f"subjet_R{R:.2f}_n{_num(n)}"
        return formula, desc, stem
    return tex_escape(label), "", re.sub(r"[^A-Za-z0-9.]+", "_", label).strip("_")


def plot_basis(weights_path, pt_min, pt_max, label, target_tree, plot_dir,
               nbins=20, pad=0.10, filt="", select=None, max_panels=0):
    """Basis-function closure plots, from tweights/tbasis_target/meta (as in the C++).
    One file per basis function, in <plot_dir>/basis_functions/g<k>_<description>.pdf."""
    import awkward as ak

    with uproot.open(weights_path) as f:
        keys = {k.split(";")[0] for k in f.keys()}
        if "tweights" not in keys or "tbasis_target" not in keys:
            print("plot_mpl: basis plots skipped (no tweights/tbasis_target)")
            return
        src = f["tweights"].arrays(["g_basis", "pt", "w_base", "w_total"], library="ak")
        tar = f["tbasis_target"].arrays(["g_basis", "weight"], library="ak")
        labels, norm = [], None
        if "meta" in keys:
            mk = f["meta"].keys()
            if "basis_labels" in mk:
                labels = list(f["meta"]["basis_labels"].array(library="np")[0])
            if "bin_center" in mk:
                norm = float(f["meta"]["bin_center"].array(library="np")[0])

    sel = (src["pt"] >= pt_min) & (src["pt"] <= pt_max)
    gs = np.abs(ak.to_numpy(src["g_basis"][sel]))
    ws_base = ak.to_numpy(src["w_base"][sel])
    ws_tot = ak.to_numpy(src["w_total"][sel])
    gt = np.abs(ak.to_numpy(tar["g_basis"]))
    wt = ak.to_numpy(tar["weight"])
    nfun = gs.shape[1] if gs.ndim == 2 else 0
    if nfun == 0:
        print("plot_mpl: basis plots skipped (g_basis empty)")
        return

    def lab(k):
        return labels[k] if k < len(labels) and labels[k] else f"g_{k}"

    idx = [k for k in range(nfun)
           if (select is None or k in select) and (not filt or filt in lab(k))]
    if max_panels > 0:
        idx = idx[:max_panels]
    if not idx:
        print("plot_mpl: basis plots skipped (selection matched nothing)")
        return

    bdir = os.path.join(plot_dir, "basis_functions")
    width = max(2, len(str(nfun - 1)))
    for k in idx:
        formula, desc, stem = basis_tex(lab(k), norm)
        vals = np.concatenate([gs[:, k], gt[:, k]])
        pos = vals[vals > 0]
        lo = pos.min() if pos.size else 1e-3
        hi = pos.max() if pos.size and pos.max() > lo else lo * 10
        lmin, lmax = math.log10(lo) - pad, math.log10(hi) + pad
        if not lmax > lmin:
            lmin, lmax = -3.0, 3.0
        edges = np.logspace(lmin, lmax, nbins + 1)
        hb = Hist.fill(gs[:, k], ws_base, edges)
        ht = Hist.fill(gs[:, k], ws_tot, edges)
        hg = Hist.fill(gt[:, k], wt, edges)
        fig, a, b = ratio_figure()
        draw_ratio_plot(a, b, hb, ht, hg, formula, "Entries", logx=True, logy=True,
                        target_name=target_tree, label=label)
        a.set_title(rf"$g_{{{k}}}$: {desc}" if desc else rf"$g_{{{k}}}$", loc="left", fontsize=10)
        save(fig, bdir, f"g{k:0{width}d}_{stem}", quiet=True)
    print(f"plot_mpl: wrote {len(idx)} plots in {bdir}")


# ---------------------------------------------------------------------------
# all: scan a files directory
# ---------------------------------------------------------------------------
RE_START = re.compile(r"startBasis_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)(.*)\.root")
RE_HISTS = re.compile(r"unbias_weights_plots_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)(.*)\.root")


def scan(files_dir):
    """{(bias, sample): {"startbasis": [(lo, hi, file)], "unbias": [(lo, hi, tag, hists, weights)]}}"""
    found = {}
    for n in sorted(os.listdir(files_dir)):
        m = RE_START.fullmatch(n) or RE_HISTS.fullmatch(n)
        if m and not parse_tag(m.group(3))[1]:
            print(f"plot_mpl: skipping {n} (no sample tag; rerun with UNBIAS_TAG set)")
            continue
        m = RE_START.fullmatch(n)
        if m:
            entry = found.setdefault(parse_tag(m.group(3)), {"startbasis": [], "unbias": []})
            entry["startbasis"].append((float(m.group(1)), float(m.group(2)), os.path.join(files_dir, n)))
        m = RE_HISTS.fullmatch(n)
        if m:
            w = os.path.join(files_dir, f"unbias_weights_{m.group(1)}_{m.group(2)}{m.group(3)}.root")
            entry = found.setdefault(parse_tag(m.group(3)), {"startbasis": [], "unbias": []})
            entry["unbias"].append((float(m.group(1)), float(m.group(2)), m.group(3).lstrip("_"),
                                    os.path.join(files_dir, n), w if os.path.exists(w) else None))
    return found


def list_inputs(files_dir):
    found = scan(files_dir)
    if not found:
        print(f"no startBasis_* / unbias_weights_plots_* files in {files_dir}")
        return
    print(f"{'bias':<10} {'sample':<10} {'pT bins':<18} available plots")
    for (bias, sample), e in sorted(found.items()):
        bins = sorted({(lo, hi) for lo, hi, *_ in e["startbasis"] + e["unbias"]})
        what = (["pt_spectrum_bias"] if e["startbasis"] else []) + \
               (["closure + weights + basis"] if e["unbias"] else ["(no unbiasing weights yet)"])
        print(f"{bias:<10} {sample:<10} {', '.join(pt_subdir(lo, hi) for lo, hi in bins):<18} {', '.join(what)}")


def plot_all(files_dir, plot_dir, biases, samples, target_tree, basis_opts):
    found = scan(files_dir)
    chosen = {k: v for k, v in found.items()
              if (not biases or k[0] in biases) and (not samples or k[1] in samples)}
    if not chosen:
        print("plot_mpl: nothing matches; available (bias, sample):",
              ", ".join(f"({b}, {s})" for b, s in sorted(found)) or "none")
        return
    for (bias, sample), e in sorted(chosen.items()):
        print(f"plot_mpl: === {bias} / {sample} ===")
        for lo, hi, path in e["startbasis"]:
            plot_startbasis(path, plot_dir)
        for lo, hi, tag, hists, weights in e["unbias"]:
            plot_unbias(hists, weights, lo, hi, tag, target_tree, plot_dir, basis_opts)


def use_tex_if_available(no_tex=False):
    """Apply the plot style; LaTeX if a TeX installation is found. Returns whether LaTeX is used."""
    usetex = not no_tex and shutil.which("latex") is not None and shutil.which("dvipng") is not None
    if not no_tex and not usetex:
        print("plot_mpl: no LaTeX installation found, using matplotlib mathtext")
    setup_style(usetex)
    return usetex


def main():
    plt.switch_backend("Agg")   # files only (when imported, e.g. from a notebook, the backend is left alone)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-tex", action="store_true", help="use mathtext instead of a LaTeX installation")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def basis_args(p):
        p.add_argument("--target-tree", default="tRef", help="unused, kept for compatibility with the run scripts")
        p.add_argument("--basis-nbins", type=int, default=20)
        p.add_argument("--basis-pad", type=float, default=0.10)
        p.add_argument("--basis-filter", default="")
        p.add_argument("--basis-select", default=None, help="comma list of basis indices")
        p.add_argument("--basis-max", type=int, default=0)

    p = sub.add_parser("startbasis", help="jetPtShift plot from a startBasis_*.root file")
    p.add_argument("input")
    p.add_argument("--plot-dir", default="../Output/plots")

    p = sub.add_parser("unbias", help="plot_unbias_weights plots")
    p.add_argument("--hists", required=True, help="unbias_weights_plots_*.root")
    p.add_argument("--weights", default=None, help="unbias_weights_*.root (for the basis panels)")
    p.add_argument("--pt-min", type=float, required=True)
    p.add_argument("--pt-max", type=float, required=True)
    p.add_argument("--tag", default="")
    p.add_argument("--plot-dir", default="../Output/plots")
    basis_args(p)

    p = sub.add_parser("list", help="bias / sample combinations available in --files-dir")
    p.add_argument("--files-dir", default="../Output/files")

    p = sub.add_parser("all", help="all plots for the chosen bias(es) and sample(s)")
    p.add_argument("--files-dir", default="../Output/files")
    p.add_argument("--plot-dir", default="../Output/plots/mpl")
    p.add_argument("--bias", default="", help="powerlaw, lognormal or a comma list (default: all)")
    p.add_argument("--sample", default="",
                   help="sample tag(s), comma list, e.g. 10M,50K (default: all); see 'list'")
    basis_args(p)

    a = ap.parse_args()
    use_tex_if_available(a.no_tex)

    if a.cmd == "list":
        list_inputs(a.files_dir)
        return
    if a.cmd == "startbasis":
        plot_startbasis(a.input, a.plot_dir)
        return
    basis_opts = dict(nbins=a.basis_nbins, pad=a.basis_pad, filt=a.basis_filter,
                      select=None if a.basis_select is None else {int(x) for x in a.basis_select.split(",") if x},
                      max_panels=a.basis_max)
    if a.cmd == "unbias":
        plot_unbias(a.hists, a.weights, a.pt_min, a.pt_max, a.tag, a.target_tree, a.plot_dir, basis_opts)
    else:
        split = lambda v: {x.strip() for x in v.split(",") if x.strip()}  # noqa: E731
        plot_all(a.files_dir, a.plot_dir, split(a.bias), split(a.sample), a.target_tree, basis_opts)


if __name__ == "__main__":
    sys.exit(main())
