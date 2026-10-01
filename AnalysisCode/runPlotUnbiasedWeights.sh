#!/bin/bash
# Validate the unbiasing with pT / EEC / basis-vector comparison plots, at
# paper-figure quality (larger fonts, proper statistical error bars,
# configurable basis-panel selection), for both pT bins.
#
# Build:
#   c++ -std=c++17 -O2 plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights
#
# NOTE: uses each bin's own shared "unbiasing_config_<bin>.env" TEnv file
# (written by startBasis.C / claude/driveBins.C, consumed by
# unbias_weights.cc) so the EEC normalization in these diagnostic plots
# always matches what that bin's fit used.
#
# NOTE: every plot below is labeled with its pT bin (e.g. "100 < p_{T} < 140
# GeV", drawn directly on the canvas) and every saved file name carries the
# bin as a "_<ptmin>_<ptmax>" suffix (e.g. plot_unbias_pT_check_100_140.pdf),
# derived automatically from --pt-min/--pt-max -- no extra flag needed, so
# the label can never drift from the window actually plotted.

# ---- Bin 1: 100-140 GeV ----
./plot_unbias_weights \
    --input unbias_weights_100_140.root \
    --target-input startBasis_100_140.root \
    --target-tree tRef \
    --config unbiasing_config_100_140.env \
    --pt-min 100 --pt-max 140 \
    --out unbias_weights_plots_100_140.root

# ---- Bin 2: 140-180 GeV ----
./plot_unbias_weights \
    --input unbias_weights_140_180.root \
    --target-input startBasis_140_180.root \
    --target-tree tRef \
    --config unbiasing_config_140_180.env \
    --pt-min 140 --pt-max 180 \
    --out unbias_weights_plots_140_180.root

# Example paper-figure run instead (bigger text, PDF+PNG, a curated subset
# of basis panels) -- uncomment / adapt as needed, per bin:
#
# ./plot_unbias_weights \
#     --input unbias_weights_100_140.root \
#     --target-input startBasis_100_140.root \
#     --target-tree tRef \
#     --config unbiasing_config_100_140.env \
#     --pt-min 100 --pt-max 140 \
#     --label-size 0.050 --title-size 0.058 --legend-size 0.044 \
#     --line-width 3 --marker-size 1.3 \
#     --formats pdf,png \
#     --basis-filter "sjmom_R0.10" --basis-ncol 3 --basis-nbins 25