# SelectionUnbiasing

Event-by-event reweighting to remove selection bias from heavy-ion jet samples, following the information-theoretic framework (arXiv:2501.17219).

## Overview
This repo provides a small pipeline to:
1. Build a biased vs unbiased toy dataset from MC ROOT files, per pT bin.
2. Learn per-jet unbiasing weights that match target moments of basis functions, with the basis (subjet radii, EEC term grid) configurable independently for each pT bin.
3. Validate the result with pT and EEC comparisons, labeled by pT bin.

PYTHIA only needs to generate jets at a single radius (default R=0.4);
subjet reclustering now happens entirely inside the unbiasing code
(`AnalysisCode/subjet_basis.h`), not at ntuple-production time.

## Basis functions
The basis is configurable and may mix two families in any combination
(`AnalysisCode/basis_functions.h`):

1. **Subjet-pT moments** (the default) — for one or more subjet radii R,

       g_n = Σ_subjets (pT_subjet)^n ,   n = 3, 4, ..., 10   (8 functions per radius)

   These consume the C/A subjets that `PYTHIA/ppjets_root.cc` already clustered
   at ntuple-production time and stored in the `subjet_pt_R0pXX` branches
   (R = 0.01 … 0.20), where available. Otherwise `AnalysisCode/subjet_basis.h`
   reclusters them on the fly from constituents using one of two
   interchangeable engines (see "Subjet reclustering" below).

2. **Energy correlators** (on by default; `--basis-eec off` to disable) —

       g = Σ_{i<k, ΔR<E} ΔR^{-A} [ln ΔR]^B z^m ,   z = pt_i pt_k / norm²

   `A`, `B`, `m`, and the angular cut `E` can all vary per term, and the
   full term grid can vary **per pT bin** — see "Per-bin configuration"
   below.

Change the compiled default basis by editing `BasisConfig` in
`basis_functions.h`, or override per-run via `--subjet-radii` /
`--subjet-R` / `--subjet-powers` / `--basis-eec` flags. **In normal use,
neither the compiled defaults nor these flags are what actually sets the
subjet radii or EEC term grid** — those come from each bin's config file
and CSV (below); the flags are overrides on top of that.

`unbias_weights.cc` writes the per-jet basis values it fit into its output
(`tweights.g_basis`, `tbasis_target.g_basis`, labels in `meta.basis_labels`), so
the validation plots always match the fitted basis without re-deriving it.

## Per-bin configuration (subjet radii, EEC terms, bin center)

Every pT bin has its own small TEnv config file (default name
`unbiasing_config.env`), written once by `startBasis.C` at the point the
bin's pT window is defined, and read back by every downstream tool
(`basis_functions.h`'s `read_bin_center()` / `read_subjet_radii()` /
`read_eec_terms_from_config()`). All three are **strict**: if the config
file or a required key is missing, the read aborts with a clear error
rather than silently falling back to some hard-coded number. There is no
hard-coded bin center, subjet radius, or EEC term grid anywhere in this
codebase outside these files.

The config file carries:
- `Unbiasing.BinCenter` — momentum-scale normalization (defaults to the
  midpoint of `[pTLow,pTHigh]`, or set explicitly via `startBasis.C`'s
  `binCenter` argument).
- `Unbiasing.PtLow` / `Unbiasing.PtHigh` — the bin's pT window.
- `Unbiasing.SubjetRadii` — a comma list of subjet radii for this bin
  (e.g. `"0.10,0.20"`).
- `Unbiasing.EECTermsFile` — path to this bin's EEC term grid, a small CSV
  file with header `A,B,m,E` and one term per row. If `startBasis.C` isn't
  given a path, it auto-generates one next to the config file from a
  compiled default grid, which can then be hand-edited for that bin.

`--subjet-radii` / `--subjet-R` / `--subjet-powers` can still override the
config-file radii/powers on the command line; there is no CLI override for
the EEC term grid — edit that bin's CSV file instead.

## Subjet reclustering

`AnalysisCode/subjet_basis.h` provides two interchangeable Cambridge/Aachen
reclustering engines, both consuming the same massless (pt,eta,phi)
constituents with E-scheme recombination, selected via `--subjet-source`:

- `recluster` — the original dependency-free homemade C/A loop.
- `fastjet` — the same algorithm and radius, clustered by the standalone
  FastJet library (`fastjet::ClusterSequence` + `fastjet::cambridge_algorithm`).

The two are kept side by side as a crosscheck of implementation. Because
`subjet_basis.h` now offers the FastJet engine, **any translation unit that
includes it requires linking against FastJet** — this is a hard dependency
going forward (see "Requirements" and the build commands below).

`--subjet-source` also accepts `precomputed` (read the stored
`subjet_pt_R0pXX` branches) and `auto` (use the stored branch when present,
otherwise recluster with the homemade engine).

## Requirements
- ROOT with C++17 support
- A standalone FastJet install (e.g. `fastjet-3.5.0`, built with
  `--enable-cgal=no` or similar at a non-system prefix), with its
  `fastjet-config` tool available
- A dataset of ROOT files with a tree like `tgenBefore` containing:
  - `nJets` (optional, for array-style trees)
  - `pt`, `eta`, `phi`, `mass`, `weight`
  - `const_pt`, `const_eta`, `const_phi` (constituent arrays)
  - `subjet_pt_R0pXX` (precomputed subjets; optional, reclustered if absent)

## Quick Start

### 1) Build the toy biased/unbiased samples, per pT bin
`AnalysisCode/claude/driveBins.C` drives `startBasis.C` once per pT bin,
each writing its own output ROOT file and its own config file (and, unless
pointed at an existing CSV, its own default EEC terms file). The checked-in
example sets up two bins (100–140 GeV and 140–180 GeV) sharing the same
subjet radii and the same example EEC term grid
(`AnalysisCode/claude/eec_terms_example.csv`, which reproduces the original
compiled default grid exactly):

```bash
root -l -b -q AnalysisCode/claude/driveBins.C
```

Edit `driveBins.C` to change the pT windows, input directory, bias
parameters, or to point a bin at its own subjet radii / EEC terms file once
you want the bins to diverge.

To build a single bin by hand instead:

```bash
root -l -q 'AnalysisCode/startBasis.C("/path/to/MCOutput/", "startBasis_output.root", 100, 140, 8.0, 0.1, 12345, -1.0, "unbiasing_config.env", "0.10,0.20", "")'
```

Outputs, per bin:
- `startBasis_<bin>.root` with trees:
  - `tX`: unbiased X jets
  - `tY`: Y jets
  - `tYprime`: oversampled Y
  - `tRef`: reference (X + Y)
  - `tBiased`: biased sample (downsampled X + Y')
  - `tPP`: same content as `tRef`, kept for backward compatibility
  - `tWeights`: analytic class weights (`wX`, `wYprime`) for the toy test
- `unbiasing_config_<bin>.env`: the shared TEnv config for this bin
  (`BinCenter`/`PtLow`/`PtHigh`/`SubjetRadii`/`EECTermsFile`)
- `jetPtShift_TChainReader_<date>.pdf`

### 2) Compute unbiasing weights, per pT bin
Build (now requires linking against FastJet):

```bash
c++ -std=c++17 -O2 AnalysisCode/unbias_weights.cc \
    $(root-config --cflags --libs) \
    $(/path/to/fastjet-install/bin/fastjet-config --cxxflags --libs) \
    -o unbias_weights
```

or use `AnalysisCode/Makefile` / `AnalysisCode/runUnbiasing.sh`, which also
runs the binary once per pT bin. If FastJet was installed to a non-system
prefix, point `LD_LIBRARY_PATH` at its `lib/` directory before running:

```bash
export LD_LIBRARY_PATH=/path/to/fastjet-install/lib:$LD_LIBRARY_PATH

./unbias_weights \
  --input startBasis_100_140.root \
  --tree tBiased \
  --target-input startBasis_100_140.root \
  --target-tree tRef \
  --config unbiasing_config_100_140.env \
  --pt-min 100 --pt-max 140 \
  --subjet-source fastjet \
  --out unbias_weights_100_140.root \
  --mode run
```

`--target-input` defaults to `--input`, so one file suffices. `--config`
must point at that bin's config file — it is what supplies the bin center,
subjet radii, and EEC term grid (see "Per-bin configuration" above); there
is no built-in fallback if it's missing.

Outputs:
- `unbias_weights_<bin>.root` with tree `tweights` containing `w_unbias`, `w_base`,
  `w_total`, and `g_basis`; plus `meta` and `tbasis_target`.

### 3) Validate with plots, per pT bin
Build and run (or use `AnalysisCode/runPlotUnbiasedWeights.sh`):

```bash
c++ -std=c++17 -O2 AnalysisCode/plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights

./plot_unbias_weights \
  --input unbias_weights_100_140.root \
  --target-input startBasis_100_140.root \
  --target-tree tRef \
  --config unbiasing_config_100_140.env \
  --pt-min 100 --pt-max 140 \
  --out unbias_weights_plots_100_140.root
```

(`plot_unbias_weights.cc` doesn't recluster — it only reads already-computed
output — so it needs no FastJet link.)

Outputs, every plot labeled on-canvas and in the file name with its pT bin
(e.g. "100 < p_T < 140 GeV", `_100_140` suffix), derived from `--pt-min`/
`--pt-max` so the label can never drift from the window actually plotted:
- `plot_unbias_pT_check_<bin>.pdf`
- `plot_unbias_weights_weights_<bin>.pdf`
- `plot_eec_compare_<bin>.pdf`
- `plot_theta_basis_weighted_compare_<bin>.pdf`

## Files
- `PYTHIA/ppjets_root.cc`: generate pp jets at a single jet radius (default R=0.4), store constituents + (optionally) precomputed C/A subjets over R = 0.01…0.20
- `PYTHIA/ppjets_yoda.cc`: same generation with YODA histogram output
- `AnalysisCode/startBasis.C`: build one pT bin's biased/unbiased toy samples; writes that bin's TEnv config (bin center, pT window, subjet radii, EEC terms file path)
- `AnalysisCode/claude/driveBins.C`: driver macro running `startBasis.C` once per pT bin (currently: 100–140 GeV, 140–180 GeV)
- `AnalysisCode/claude/eec_terms_example.csv`: example per-bin EEC term grid (reproduces the original compiled default grid)
- `AnalysisCode/basis_functions.h`: configurable basis framework (subjet moments + EEC), per-bin config-file readers (`read_bin_center`/`read_subjet_radii`/`read_eec_terms_from_config`), EEC terms CSV reader/writer
- `AnalysisCode/subjet_basis.h`: subjet reclustering — two interchangeable C/A engines (homemade, standalone-FastJet); requires linking against FastJet
- `AnalysisCode/unbias_weights.cc`: learn unbiasing weights (Adam optimizer)
- `AnalysisCode/plot_unbias_weights.cc`: validate unbiasing with pT/EEC/basis plots, labeled per pT bin
- `AnalysisCode/plot_loss_from_nohup.py`: plot loss vs iteration from logs
- `AnalysisCode/runUnbiasing.sh`: example unbiasing command, run once per pT bin
- `AnalysisCode/runPlotUnbiasedWeights.sh`: example plotting command, run once per pT bin
- `AnalysisCode/claude/Makefile`: build target for `unbias_weights` (ROOT + FastJet)
- `AnalysisCode/claude/todo.md`: running to-do list for this analysis chain

## Notes
- `--adam-lr auto` (default) computes a safe learning rate based on max |g|.
- `--subjet-R` (default 0.1) selects a single subjet radius; `--subjet-radii`
  takes a comma-separated list. If neither is passed, the subjet radii come
  from `--config`'s `Unbiasing.SubjetRadii` instead.
- `--subjet-source` (`auto` | `precomputed` | `recluster` | `fastjet`)
  chooses whether to read the stored subjet branches, recluster with the
  homemade C/A engine, or recluster with the standalone-FastJet engine.
- `--config <file>` (default `unbiasing_config.env`) selects the per-bin TEnv
  config file; it must carry `Unbiasing.BinCenter`, `Unbiasing.SubjetRadii`,
  and `Unbiasing.EECTermsFile`, or the run aborts.
- `--mode debug` prints per-iteration gradient and covariance diagnostics.
- For array-style trees, `nJets` + arrays are expected; for flat trees, `pt` is a scalar.

## Citation