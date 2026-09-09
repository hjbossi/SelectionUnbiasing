#!/bin/bash
# Validate the unbiasing with pT / EEC / basis-vector comparison plots.
#
# Build:
#   c++ -std=c++17 -O2 plot_unbias_weights.cc $(root-config --cflags --libs) -o plot_unbias_weights

./plot_unbias_weights \
    --input unbias_weights.root \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --pt-min 100 --pt-max 140
