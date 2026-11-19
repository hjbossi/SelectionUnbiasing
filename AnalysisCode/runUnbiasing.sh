#!/bin/bash
# Compute unbiasing weights from the startBasis.C output, for both pT bins.
#
# Build:
#   c++ -std=c++17 -O2 unbias_weights.cc \
#       $(root-config --cflags --libs) \
#       $($FJI/bin/fastjet-config --cxxflags --libs) \
#       -o unbias_weights
# or just `make` (see makefile; needs $FJI set or fastjet-config on PATH).
#
# NOTE: run driveBins.C first (root -l -b -q 'driveBins.C("/path/to/MCOutput/")'). It
# calls startBasis.C once per bin, writing each bin's own
# startBasis_<bin>.root, unbiasing_config_<bin>.env (bin center + pT window +
# subjet radii + EEC terms file path), and
# unbiasing_config_<bin>_eec_terms.csv (the EEC term grid).
#
# NOTE: each unbiasing_config_<bin>.env must carry "Unbiasing.SubjetRadii"
# and "Unbiasing.EECTermsFile" (both written by the current startBasis.C /
# driveBins.C) -- unbias_weights.cc aborts if either is missing. If your
# config files predate this, rerun driveBins.C to regenerate them.
#
# NOTE: subjet_basis.h now links against a standalone FastJet install
# (in addition to the original dependency-free homemade C/A reclustering,
# still available via --subjet-source recluster as a crosscheck). FastJet
# was installed to a non-system prefix, so its shared library is not on the
# default linker search path -- point LD_LIBRARY_PATH at it before running
# the built binary.
# (Binaries built with `make` already carry the FastJet lib dir as an rpath,
# so this is only needed for hand-built binaries. Set FJI to the FastJet
# install prefix; on macOS the variable is DYLD_LIBRARY_PATH.)
if [ -n "$FJI" ]; then
  export LD_LIBRARY_PATH=$FJI/lib:$LD_LIBRARY_PATH
  export DYLD_LIBRARY_PATH=$FJI/lib:$DYLD_LIBRARY_PATH
fi

# Run from this script's directory (AnalysisCode/), whatever the caller's cwd.
cd "$(dirname "$0")" || exit 1
# Generated files live in ../Output/files/ (git-ignored), written there by driveBins.C.
FILES=../Output/files
mkdir -p "$FILES"
# Sample tag (required, e.g. export UNBIAS_TAG=10M): must match the one used for
# driveBins.C; every file name below ends in _<tag>.
if [ -z "$UNBIAS_TAG" ]; then
  echo "$(basename "$0"): UNBIAS_TAG is not set. Name the input sample, e.g. export UNBIAS_TAG=10M" >&2
  exit 1
fi
S=_$UNBIAS_TAG
# Fit length (see unbias_weights.cc): UNBIAS_MAX_ITER = number of iterations (default 20)
# or "auto": run until every moment matches its target to UNBIAS_CONVERGE_ACCURACY
# (relative, default 5e-3), or the loss stops improving, or 5000 iterations.
FIT_ARGS="--max-iter ${UNBIAS_MAX_ITER:-20} --converge-accuracy ${UNBIAS_CONVERGE_ACCURACY:-5e-3}"

# --subjet-source fastjet: recluster subjets on the fly from constituents
# using the standalone-FastJet Cambridge/Aachen engine in subjet_basis.h
# (recluster_fastjet_subjet_pts), rather than reading precomputed
# subjet_pt_R0pXX branches. Use --subjet-source recluster instead to run
# the homemade dependency-free C/A engine as a crosscheck.

# ---- Bin 1: 100-140 GeV ----
./unbias_weights \
    --input $FILES/startBasis_100_140${S}.root \
    --tree tBiased \
    --target-input $FILES/startBasis_100_140${S}.root \
    --target-tree tRef \
    --config $FILES/unbiasing_config_100_140${S}.env \
    --pt-min 100 --pt-max 140 \
    --subjet-source fastjet \
    --out $FILES/unbias_weights_100_140${S}.root \
    --mode run ${FIT_ARGS}

# ---- Bin 2: 140-180 GeV ----
./unbias_weights \
    --input $FILES/startBasis_140_180${S}.root \
    --tree tBiased \
    --target-input $FILES/startBasis_140_180${S}.root \
    --target-tree tRef \
    --config $FILES/unbiasing_config_140_180${S}.env \
    --pt-min 140 --pt-max 180 \
    --subjet-source fastjet \
    --out $FILES/unbias_weights_140_180${S}.root \
    --mode run ${FIT_ARGS}