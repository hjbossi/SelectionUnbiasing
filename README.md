# SelectionUnbiasing

Event-by-event reweighting to remove selection bias from heavy-ion jet samples, following the information-theoretic framework of Andres, Bossi, Holguin (arXiv:2501.17219).

## Overview
This repo provides a small pipeline to:
1. Build a biased vs unbiased toy dataset from MC ROOT files.
2. Learn per-jet unbiasing weights that match target moments of basis functions.
3. Validate the result with pT and EEC comparisons.

## Requirements
- ROOT with C++17 support
- A dataset of ROOT files with a tree like `tgenBefore` containing:
  - `nJets` (optional, for array-style trees)
  - `pt`, `eta`, `phi`, `mass`, `weight`
  - `const_pt`, `const_eta`, `const_phi` (constituent arrays)

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
- `jetPtShift_TChainReader.pdf`

### 2) Compute unbiasing weights
Build and run:

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

Outputs:
- `unbias_weights.root` with tree `tweights` containing `w_unbias`, `w_base`, `w_total`.

### 3) Validate with plots
Build and run:

```bash
c++ -O2 AnalysisCode/plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights

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
- `AnalysisCode/startBasis.C`: build biased/unbiased toy samples
- `AnalysisCode/unbias_weights.cc`: learn unbiasing weights (Adam optimizer)
- `AnalysisCode/plot_unbias_weights.cc`: validate unbiasing with pT/EEC/basis plots
- `AnalysisCode/plot_loss_from_nohup.py`: plot loss vs iteration from logs
- `AnalysisCode/runUnbiasing.sh`: example unbiasing command
- `AnalysisCode/runPlotUnbiasedWeights.sh`: example plotting command

## Notes
- `--adam-lr auto` (default) computes a safe learning rate based on max |g|.
- `--dR-min` screens small-angle divergences in basis functions.
- For array-style trees, `nJets` + arrays are expected; for flat trees, `pt` is a scalar.

## Citation
If you use this method, please cite:
- Andres, Bossi, Holguin, arXiv:2501.17219
