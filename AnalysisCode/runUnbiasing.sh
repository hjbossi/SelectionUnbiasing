./unbias_weights \
     --input UnbiasingTest_PYTHIApp_pthatmin50_121825.root \
     --tree tgenBefore --n-branch nJets --pt-branch pt --weight-branch weight \
     --target-hist target.root:hpT \
     --moments 1,2,3,4 --pt0 1.0 --pt-min 50 --pt-max 200 \
     --out unbias_weights.root