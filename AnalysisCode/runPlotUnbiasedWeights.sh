#!/bin/bash
# Validate the unbiasing with pT / EEC / basis-vector comparison plots, at
# paper-figure quality (larger fonts, proper statistical error bars,
# configurable basis-panel selection).
#
# Build:
#   c++ -std=c++17 -O2 plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights
#
# NOTE: uses the same shared "unbiasing_config.env" TEnv file (written by
# startBasis.C, consumed by unbias_weights.cc) so the EEC normalization in
# these diagnostic plots always matches what the fit itself used.
#
# Default run (all basis panels, default paper-style sizing):
./plot_unbias_weights \
    --input unbias_weights.root \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --config unbiasing_config.env \
    --pt-min 100 --pt-max 140

# Example paper-figure run instead (bigger text, PDF+PNG, a curated subset
# of basis panels) -- uncomment / adapt as needed:
#
# ./plot_unbias_weights \
#     --input unbias_weights.root \
#     --target-input startBasis_output.root \
#     --target-tree tRef \
#     --config unbiasing_config.env \
#     --pt-min 100 --pt-max 140 \
#     --label-size 0.050 --title-size 0.058 --legend-size 0.044 \
#     --line-width 3 --marker-size 1.3 \
#     --formats pdf,png \
#     --basis-filter "sjmom_R0.10" --basis-ncol 3 --basis-nbins 25