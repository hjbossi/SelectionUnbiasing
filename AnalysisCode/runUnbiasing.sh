#!/bin/bash
# Compute unbiasing weights from the startBasis.C output.
#
# Build:
#   c++ -std=c++17 -O2 unbias_weights.cc \
#       $(root-config --cflags --libs) \
#       $(/data/ALEPH/MC/mcgen/fastjet-install/bin/fastjet-config --cxxflags --libs) \
#       -o unbias_weights
#
# NOTE: run startBasis.C first. It writes the shared "unbiasing_config.env"
# TEnv file (bin center + pT window) that unbias_weights.cc reads below via
# --config -- this is the single place the bin-center normalization (used to
# come from a hard-coded 120.0) is set.
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
./unbias_weights \
    --input startBasis_output.root \
    --tree tBiased \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --config unbiasing_config.env \
    --pt-min 100 --pt-max 140 \
    --subjet-source fastjet \
    --out unbias_weights.root \
    --mode run