#!/bin/bash
# Compute unbiasing weights from the startBasis.C output.
#
# Build:
#   c++ -std=c++17 -O2 unbias_weights.cc $(root-config --cflags --libs) -o unbias_weights

./unbias_weights \
    --input startBasis_output.root \
    --tree tBiased \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --pt-min 100 --pt-max 140 \
    --out unbias_weights.root \
    --mode run
