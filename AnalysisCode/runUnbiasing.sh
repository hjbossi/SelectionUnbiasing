#!/bin/bash
# Compute unbiasing weights from the startBasis.C output, for both pT bins.
#
# Build:
#   c++ -std=c++17 -O2 unbias_weights.cc \
#       $(root-config --cflags --libs) \
#       $(/data/ALEPH/MC/mcgen/fastjet-install/bin/fastjet-config --cxxflags --libs) \
#       -o unbias_weights
#
# NOTE: run claude/driveBins.C first (root -l -b -q claude/driveBins.C). It
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
export LD_LIBRARY_PATH=/data/ALEPH/MC/mcgen/fastjet-install/lib:$LD_LIBRARY_PATH

# --subjet-source fastjet: recluster subjets on the fly from constituents
# using the standalone-FastJet Cambridge/Aachen engine in subjet_basis.h
# (recluster_fastjet_subjet_pts), rather than reading precomputed
# subjet_pt_R0pXX branches. Use --subjet-source recluster instead to run
# the homemade dependency-free C/A engine as a crosscheck.

# ---- Bin 1: 100-140 GeV ----
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

# ---- Bin 2: 140-180 GeV ----
./unbias_weights \
    --input startBasis_140_180.root \
    --tree tBiased \
    --target-input startBasis_140_180.root \
    --target-tree tRef \
    --config unbiasing_config_140_180.env \
    --pt-min 140 --pt-max 180 \
    --subjet-source fastjet \
    --out unbias_weights_140_180.root \
    --mode run