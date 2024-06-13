#!bin/bash

run_prefix=$1
run1="NC2"
run2="NC3"

folder1="${run_prefix}_${run1}_png"
folder2="${run_prefix}_${run2}_png"

for obs in "123"*"_observation.png"; do
    gname="${obs%_observation.png}"
    model1="$folder1/${gname}_model.png"
    out1="$folder1/${gname}_out.png"
    residual1="$folder1/${gname}_residual.png"
    
    #model2="$folder2/${gname}_model.png"
    out2="$folder2/${gname}_out.png"
    residual2="$folder2/${gname}_residual.png"
    
    montage "$obs" "$model1" "$out1" "$residual1" "$out2" "$residual2" -tile "6x1" -geometry "175x175+2+2" "${gname}_components_combined.png"
done