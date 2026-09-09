# SelectionUnbiasing

Event-by-event reweighting to remove selection bias from heavy-ion jet samples, following the information-theoretic framework (arXiv:2501.17219).

## Overview
This repo provides a small pipeline to:
1. Build a biased vs unbiased toy dataset from MC ROOT files.
2. Learn per-jet unbiasing weights that match target moments of basis functions.
3. Validate the result with pT and EEC comparisons.

## Basis functions
The basis is configurable and may mix two families in any combination
(`AnalysisCode/basis_functions.h`):

1. **Subjet-pT moments** (the default) — for one or more subjet radii R,

       g_n = Σ_subjets (pT_subjet)^n ,   n = 3, 4, ..., 10   (8 functions at R=0.1)

   These consume the C/A subjets that `PYTHIA/ppjets_root.cc` already clustered
   at ntuple-production time and stored in the `subjet_pt_R0pXX` branches
   (R = 0.01 … 0.20). `AnalysisCode/subjet_basis.h` provides an on-the-fly
   reclustering fallback for older constituent-only ntuples.

2. **Energy correlators** (off by default, `--basis-eec on`) —

       g = Σ_{i<k, ΔR<E} ΔR^{-A} [ln ΔR]^B z^m ,   z = pt_i pt_k / norm²

Change the basis either by editing `BasisConfig` in `basis_functions.h` or with
the `--subjet-radii` / `--subjet-powers` / `--basis-eec` flags. The
implementation needs no FastJet, so the tools build with `root-config` alone.

`unbias_weights.cc` writes the per-jet basis values it fit into its output
(`tweights.g_basis`, `tbasis_target.g_basis`, labels in `meta.basis_labels`), so
the validation plots always match the fitted basis without re-deriving it.

## Requirements
- ROOT with C++17 support
- A dataset of ROOT files with a tree like `tgenBefore` containing:
  - `nJets` (optional, for array-style trees)
  - `pt`, `eta`, `phi`, `mass`, `weight`
  - `const_pt`, `const_eta`, `const_phi` (constituent arrays)
  - `subjet_pt_R0pXX` (precomputed subjets; optional, reclustered if absent)

## Quick Start
### 1) Build the toy biased/unbiased samples
Edit `AnalysisCode/startBasis.C` if needed (input directory, pT window, bias parameters), then run:

```bash
root -l -q 'AnalysisCode/startBasis.C("/path/to/MCOutput/", "startBasis_output.root", 100, 140, 8.0, 0.1, 12345)'
```

Outputs:
- `startBasis_output.root` with trees:
  - `tX`: unbiased X jets
  - `tY`: Y jets
  - `tYprime`: oversampled Y
  - `tRef`: reference (X + Y)
  - `tBiased`: biased sample (downsampled X + Y')
  - `tPP`: same content as `tRef`, kept for backward compatibility
  - `tWeights`: analytic class weights (`wX`, `wYprime`) for the toy test
- `jetPtShift_TChainReader.pdf`

### 2) Compute unbiasing weights
Build and run (or use `AnalysisCode/runUnbiasing.sh`):

```bash
c++ -std=c++17 -O2 AnalysisCode/unbias_weights.cc $(root-config --cflags --libs) -o unbias_weights

./unbias_weights \
  --input startBasis_output.root \
  --tree tBiased \
  --target-input startBasis_output.root \
  --target-tree tRef \
  --pt-min 100 --pt-max 140 \
  --out unbias_weights.root \
  --mode run
```

`--target-input` defaults to `--input`, so one file suffices.

Outputs:
- `unbias_weights.root` with tree `tweights` containing `w_unbias`, `w_base`,
  `w_total`, and `g_basis`; plus `meta` and `tbasis_target`.

### 3) Validate with plots
Build and run (or use `AnalysisCode/runPlotUnbiasedWeights.sh`):

```bash
c++ -std=c++17 -O2 AnalysisCode/plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights

./plot_unbias_weights \
  --input unbias_weights.root \
  --target-input startBasis_output.root \
  --target-tree tRef \
  --pt-min 100 --pt-max 140
```

Outputs:
- `plot_unbias_pT_check.pdf`
- `plot_unbias_weights_weights.pdf`
- `plot_eec_compare.pdf`
- `plot_theta_basis_weighted_compare.pdf`

## Files
- `PYTHIA/ppjets_root.cc`: generate pp jets, store constituents + C/A subjets over R = 0.01…0.20
- `PYTHIA/ppjets_yoda.cc`: same generation with YODA histogram output
- `AnalysisCode/startBasis.C`: build biased/unbiased toy samples
- `AnalysisCode/basis_functions.h`: configurable basis framework (subjet moments + EEC)
- `AnalysisCode/subjet_basis.h`: self-contained C/A reclustering (fallback for older ntuples)
- `AnalysisCode/unbias_weights.cc`: learn unbiasing weights (Adam optimizer)
- `AnalysisCode/plot_unbias_weights.cc`: validate unbiasing with pT/EEC/basis plots
- `AnalysisCode/plot_loss_from_nohup.py`: plot loss vs iteration from logs
- `AnalysisCode/runUnbiasing.sh`: example unbiasing command
- `AnalysisCode/runPlotUnbiasedWeights.sh`: example plotting command

## Notes
- `--adam-lr auto` (default) computes a safe learning rate based on max |g|.
- `--subjet-R` (default 0.1) selects a single subjet radius; `--subjet-radii`
  takes a comma-separated list.
- `--subjet-source` (`auto` | `precomputed` | `recluster`) chooses whether to
  read the stored subjet branches or recluster from constituents.
- `--mode debug` prints per-iteration gradient and covariance diagnostics.
- For array-style trees, `nJets` + arrays are expected; for flat trees, `pt` is a scalar.

## Citation
