
#   ./unbias_weights \
#     --input startBasis_output.root \
#     --tree tY \
#     --basis moments \
#     --pt-branch pt_raw \
#     --target-input startBasis_output.root \
#     --target-tree tX \
#     --moments 1,2,3,4 --pt0 1.0 \
#     --pt-min 90 --pt-max 140 \
#     --out unbias_weights.root


# --adam-lr 0.001 
    
  # ./unbias_weights \
  #   --input startBasis_output.root \
  #   --tree tBiased \
  #   --basis theta \
  #   --target-input startBasis_output.root \
  #   --target-tree tRef \
  #   --pt-min 100 --pt-max 140 \
  #   --optimizer adam \
  #   --scale-basis target \
  #   --out unbias_weights.root
    
./unbias_weights \
    --input startBasis_output.root \
    --tree tBiased \
    --target-input startBasis_output.root \
    --target-tree tRef \
    --pt-min 100 --pt-max 140 \
    --out unbias_weights.root \
    --mode run