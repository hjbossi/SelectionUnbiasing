#!/bin/bash
# Compute unbiasing weights from the startBasis.C output.
#
# Build:
#   c++ -std=c++17 -O2 unbias_weights.cc $(root-config --cflags --libs) -o unbias_weights
#
# NOTE: run startBasis.C first. It writes the shared "unbiasing_config.env"
# TEnv file (bin center + pT window) that unbias_weights.cc reads below via
# --config -- this is the single place the bin-center normalization (used to
# come from a hard-coded 120.0) is set.

./unbias_weights \
    --input startBasis_output.root \
    --tree tBiased \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --config unbiasing_config.env \
    --pt-min 100 --pt-max 140 \
    --out unbias_weights.root \
    --mode run