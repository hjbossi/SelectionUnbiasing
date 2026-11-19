# Setup on a new machine

Everything below runs from `AnalysisCode/` unless stated otherwise.

## 1. Dependencies

ROOT 6, FastJet 3 (with `fastjet-config`), a C++17 compiler, PYTHIA 8 (only to generate test input),
and for the Python notebooks `numpy`, `matplotlib`, `uproot` and `ipykernel`.

Example with conda (macOS or Linux):

```bash
conda create -n unbias -c conda-forge --strict-channel-priority \
    root fastjet compilers python=3.12 matplotlib numpy uproot ipykernel pythia8
```

To add the notebook packages to an existing env: `conda install -n unbias -c conda-forge uproot ipykernel`.

The conda FastJet package may not ship `fastjet-config`. If so, build FastJet from source:

```bash
curl -LO https://fastjet.fr/repo/fastjet-3.5.0.tar.gz
tar xzf fastjet-3.5.0.tar.gz && cd fastjet-3.5.0
./configure --prefix=$HOME/fastjet-install --enable-cgal=no --disable-auto-ptr
make -j8 && make install
```

## 2. Environment (once per shell)

```bash
conda activate unbias                     # ROOT, compilers, PYTHIA (sets $CONDA_PREFIX, PYTHIA8DATA)
export FJI=$HOME/fastjet-install          # FastJet install prefix (contains bin/fastjet-config)
export UNBIAS_INPUT_DIR=/path/to/MCOutput/ # directory with the input .root files
cd AnalysisCode
```

If the env was created with `-p /some/path` (or moved), conda does not know it by name: activate it
by path instead, `conda activate /some/path/unbias`.

In a script or other non-interactive shell, `conda activate` may not work. Then put the env's `bin/`
first on `PATH` instead; that is enough for building and running the analysis (not for generating
PYTHIA test input, which needs the activated env):

```bash
export PATH=/path/to/conda-env/bin:$PATH
```

Do not set `LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` for the build: `make` records FastJet's library
directory in the binaries (rpath), so a second FastJet (e.g. inside the conda env) is not picked up
by mistake. That mix-up causes a segfault on macOS.

## 3. Build

```bash
make        # builds unbias_weights and plot_unbias_weights
```

`make` finds `fastjet-config` via `$FJI`, then `PATH`; override with `make FASTJET_CONFIG=/path/to/fastjet-config`.

## 4. Run

```bash
root -l -b -q driveBins.C      # builds the toy biased/reference samples (uses $UNBIAS_INPUT_DIR)
./runUnbiasing.sh              # computes the weights
./runPlotUnbiasedWeights.sh    # validation plots
```

`driveBins.C` aborts with a message if the input directory is missing or contains no entries.
Every `.root` file in that directory is read, so keep unrelated ROOT files out of it.

All outputs go to `Output/` at the top of the repository (created automatically, git-ignored):
ROOT, config and CSV files in `Output/files/`, plots in `Output/plots/root/` (ROOT) and
`Output/plots/mpl/` (matplotlib, see below).

Every run needs a sample tag naming the input sample, e.g. `export UNBIAS_TAG=10M` (the chain
stops if it is not set). Every output name ends in it (`startBasis_100_140_10M.root`,
`startBasis_100_140_lognormal_10M.root`, ...), so outputs from different samples never overwrite each
other. Use the same tag for all steps (`driveBins*.C`, the `run*.sh` scripts) and set `TAG` in the
notebook to match.

### Fit length and convergence

`unbias_weights` runs a fixed 20 iterations by default, which is usually not converged. To run each
fit until it converges instead:

```bash
export UNBIAS_MAX_ITER=auto            # default 20 (fixed); or any number of iterations
export UNBIAS_CONVERGE_ACCURACY=5e-3   # optional: stop once every moment matches its target to 0.5%
```

before `runUnbiasing.sh` / `runLognormal.sh`. With `auto`, a fit stops when the largest relative
mismatch between reweighted and target moments (`mx_rel` in the log) drops below the accuracy, or
when the loss stops improving, or after 5000 iterations. The log ends each fit with
`Best at iter N ... [converged]` or `[NOT converged: stopped at max-iter]`, and the `meta` tree of
`unbias_weights_*.root` records `n_iter` and `converged`. A smaller accuracy (e.g. `2e-3`) is
tighter but takes longer; there is no point going far below the statistical precision of the
moments (roughly 1/sqrt(number of jets in the bin)).

### ROOT or matplotlib plots

`UNBIAS_PLOTTER` chooses who draws the plots in `Output/plots/`:

| value | plots |
|---|---|
| `root` (default) | ROOT only, in `Output/plots/root/` |
| `matplotlib` | matplotlib (`plot_mpl.py`) only, in `Output/plots/mpl/<bias>/<sample>/pt<lo>-<hi>/` |
| `both` | both of the above |

`export UNBIAS_PLOTTER=both` before `driveBins*.C` and the `run*.sh` scripts. The matplotlib plots are
drawn from the ROOT files the chain writes (`plot_unbias_weights` still runs, to fill the histograms),
with the env's `python3` (override with `UNBIAS_PYTHON`). Rebuild once with `make` after updating
(`plot_unbias_weights` gained `--formats none`).

The matplotlib plots are step histograms with LaTeX text (needs a TeX installation with `latex`
and `dvipng` on `PATH`, e.g. MacTeX; otherwise `plot_mpl.py` falls back to matplotlib's mathtext,
or force that with `--no-tex`) and the colour-blind-safe Okabe-Ito palette. Layout:

```
<bias>/<sample>/pt<lo>-<hi>/<plot>.pdf      e.g. lognormal/10M/pt100-140/closure_pt.pdf
```

`<bias>` is `powerlaw` or `lognormal`; `<sample>` is the `UNBIAS_TAG` the inputs were made with.

| matplotlib plot | ROOT plot |
|---|---|
| `pt_spectrum_bias` | `jetPtShift_TChainReader` |
| `closure_pt` | `plot_unbias_pT_check` |
| `unbias_weights` | `plot_unbias_weights_weights` |
| `closure_eec` | `plot_eec_compare` |
| `closure_jet_mass`, `closure_rg`, `closure_kt`, `closure_broadening`, `closure_tau*` | `plot_jet_mass`, `plot_rg`, ... |
| `basis_functions/g<k>_<description>.pdf`, one per basis function, e.g. `g20_subjet_R0.10_n1.pdf` | `plot_theta_basis_weighted_compare` |
| `pt_bias_vs_expected` (made by `plot_pt_bias.ipynb`, pT 100–140 only) | none |

To (re)draw matplotlib plots from existing outputs without rerunning anything (from `AnalysisCode/`;
inputs from `../Output/files`, plots to `../Output/plots/mpl`):

```bash
python plot_mpl.py list                                   # which bias / sample sets exist, and what they contain
python plot_mpl.py all --bias lognormal --sample 10M      # one set
python plot_mpl.py all                                    # everything
```

`--bias` and `--sample` take comma lists. The closure, weight and basis plots need the unbiasing
weights (`runUnbiasing.sh` / `runLognormal.sh`) for that set; until then only `pt_spectrum_bias` exists.

## Test input (no production data)

`PYTHIA/ppjets_root.cc` generates a PYTHIA sample (jets + constituents only, ~1.1 kB/event).
`PYTHIA/generate.sh` compiles it and runs it as parallel jobs with different random seeds.
Run from the repository root, with the conda env activated (for `$CONDA_PREFIX`) and `FJI` set:

```bash
PYTHIA/generate.sh <nEvents> <nJobs> <outDir> [firstSeed]
# e.g. 10M events on 16 cores (~15 min, ~11 GB):
PYTHIA/generate.sh 10000000 16 /path/to/MCOutput_10M
export UNBIAS_INPUT_DIR=/path/to/MCOutput_10M/
```

It writes one `pythia_<i>.root` per job; `startBasis.C` reads all of them, so no merging is needed.
Speed is ~1.5 ms/event per core. Per 1M events there are ~58k jets in 100–140 GeV and ~11k in
140–180 GeV. 50 000 events is enough for a smoke test; use several million for usable statistics.
To add statistics to an existing sample, generate into a new directory with a different `firstSeed`
(e.g. 1000), otherwise the events repeat.

The binary can also be run directly: `ppjets [nEvent] [outFile] [seed]` (seed 0 = PYTHIA default).

## Notebooks (Jupyter / VS Code)

The notebooks (e.g. `AnalysisCode/plot_pt_bias.ipynb`) run in the conda env. Register it once as a
Jupyter kernel:

```bash
conda activate unbias      # or: conda activate /path/to/unbias
python -m ipykernel install --user --name unbias --display-name "Python (unbias)"
```

In VS Code (with the Python and Jupyter extensions): open the notebook, click **Select Kernel**
(top right) → **Jupyter Kernel…** → **Python (unbias)**. In Jupyter Lab/Notebook, pick
**Python (unbias)** from the kernel menu. The notebooks read from `../Output/files/`, so run the
corresponding chain first (`driveBins.C` / `driveBinsLognormal.C`). `plot_pt_bias.ipynb` uses the
style of `plot_mpl.py` and saves into `Output/plots/mpl/<bias>/<TAG>/pt<lo>-<hi>/`; set `TAG` in it.
