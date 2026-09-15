#!/bin/bash
# Validate the unbiasing with pT / EEC / basis-vector comparison plots.
#
# Build:
#   c++ -std=c++17 -O2 plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights
#
# NOTE: uses the same shared "unbiasing_config.env" TEnv file (written by
# startBasis.C, consumed by unbias_weights.cc) so the EEC normalization in
# these diagnostic plots always matches what the fit itself used.

./plot_unbias_weights \
    --input unbias_weights.root \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --config unbiasing_config.env \
    --pt-min 100 --pt-max 140