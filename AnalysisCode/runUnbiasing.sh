
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
    
    
  ./unbias_weights \
    --input startBasis_output.root \
    --tree tY \
    --basis combined \
    --pt-branch pt_raw \
    --target-input startBasis_output.root \
    --target-tree tX \
    --moments 1,2,3,4 --pt0 1.0 --pt-min 90 --pt-max 140 \
    --out unbias_weights.root