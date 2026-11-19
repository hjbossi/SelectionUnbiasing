#!/bin/bash
# Full chain with the log-normal X bias (see BIAS.md), for both pT bins:
#   1. driveBinsLognormal.C  -> biased/reference trees + per-bin configs
#   2. unbias_weights        -> unbiasing weights
#   3. plot_unbias_weights   -> validation plots
# Every output carries "_lognormal" in its name, so it never overwrites the
# power-law outputs of driveBins.C / runUnbiasing.sh / runPlotUnbiasedWeights.sh.
#
# Needs the same environment as the power-law chain (see SETUP.md):
#   ROOT on PATH, FJI (or fastjet-config on PATH), and UNBIAS_INPUT_DIR
#   (or pass the input directory as the first argument).
set -e

# Run from this script's directory (AnalysisCode/), whatever the caller's cwd.
cd "$(dirname "$0")"
# Generated files live in ../Output/ (git-ignored): ROOT files in files/, plots in plots/.
FILES=../Output/files
PLOTS=../Output/plots
mkdir -p "$FILES" "$PLOTS"
# Sample tag (required, e.g. export UNBIAS_TAG=10M): every name ends in _lognormal_<tag>.
if [ -z "$UNBIAS_TAG" ]; then
  echo "$(basename "$0"): UNBIAS_TAG is not set. Name the input sample, e.g. export UNBIAS_TAG=10M" >&2
  exit 1
fi
T=lognormal_$UNBIAS_TAG
# Fit length (see unbias_weights.cc): UNBIAS_MAX_ITER = number of iterations (default 20)
# or "auto": run until every moment matches its target to UNBIAS_CONVERGE_ACCURACY
# (relative, default 5e-3), or the loss stops improving, or 5000 iterations.
FIT_ARGS="--max-iter ${UNBIAS_MAX_ITER:-20} --converge-accuracy ${UNBIAS_CONVERGE_ACCURACY:-5e-3}"
# Plotting backend (UNBIAS_PLOTTER, see SETUP.md): root (default), matplotlib or both.
# ROOT plots go to $PLOTS/root, matplotlib plots (plot_mpl.py) to $PLOTS/mpl.
# plot_unbias_weights always runs: it fills the histograms that plot_mpl.py reads.
PLOTTER=${UNBIAS_PLOTTER:-root}
PY=${UNBIAS_PYTHON:-python3}
ROOT_DIR=$PLOTS/root
case "$PLOTTER" in
  root)       ROOT_FORMATS=pdf;  MPL_DIR= ;;
  matplotlib) ROOT_FORMATS=none; MPL_DIR=$PLOTS/mpl ;;
  both)       ROOT_FORMATS=pdf;  MPL_DIR=$PLOTS/mpl ;;
  *) echo "UNBIAS_PLOTTER=$PLOTTER not understood (root, matplotlib or both)" >&2; exit 1 ;;
esac
mkdir -p "$ROOT_DIR"

INPUT_DIR=${1:-$UNBIAS_INPUT_DIR}
if [ -z "$INPUT_DIR" ]; then
  echo "runLognormal.sh: no input directory. Set UNBIAS_INPUT_DIR or pass it as the first argument." >&2
  exit 1
fi

# 0. build (no-op if up to date; binaries carry FastJet's lib dir as an rpath)
make

# 1. biased/reference trees + configs
root -l -b -q "driveBinsLognormal.C(\"$INPUT_DIR\")"

# 2. unbiasing weights
for B in "100 140" "140 180"; do
  set -- $B; LO=$1; HI=$2
  ./unbias_weights \
    --input $FILES/startBasis_${LO}_${HI}_${T}.root --tree tBiased \
    --target-input $FILES/startBasis_${LO}_${HI}_${T}.root --target-tree tRef \
    --config $FILES/unbiasing_config_${LO}_${HI}_${T}.env \
    --pt-min $LO --pt-max $HI --subjet-source fastjet \
    --out $FILES/unbias_weights_${LO}_${HI}_${T}.root --mode run ${FIT_ARGS}
done

# 3. validation plots (file names get _<LO>_<HI>_lognormal)
SUBSTRUCTURE_BINNING="--mass-bins 30 --rg-bins 13 --kt-bins 20 --broadening-bins 25 --tau-bins 25"
for B in "100 140" "140 180"; do
  set -- $B; LO=$1; HI=$2
  ./plot_unbias_weights \
    --input $FILES/unbias_weights_${LO}_${HI}_${T}.root \
    --target-input $FILES/startBasis_${LO}_${HI}_${T}.root --target-tree tRef \
    --config $FILES/unbiasing_config_${LO}_${HI}_${T}.env \
    --pt-min $LO --pt-max $HI --tag $T \
    --jet-R 0.4 --sd-zcut 0.1 --sd-beta 1.0 --nsub-beta 1.0 \
    ${SUBSTRUCTURE_BINNING} \
    --plot-dir $ROOT_DIR --formats $ROOT_FORMATS \
    --out $FILES/unbias_weights_plots_${LO}_${HI}_${T}.root
  if [ -n "$MPL_DIR" ]; then
    $PY plot_mpl.py unbias \
        --hists $FILES/unbias_weights_plots_${LO}_${HI}_${T}.root \
        --weights $FILES/unbias_weights_${LO}_${HI}_${T}.root \
        --pt-min $LO --pt-max $HI --tag $T --target-tree tRef --plot-dir $MPL_DIR
  fi
done
