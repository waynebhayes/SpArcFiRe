#!/bin/bash

in_dir=$1
out_dir=$2

default_in="$(pwd)"/"sparcfire-in"
default_tmp="$(pwd)"/"sparcfire-tmp"
default_out="$(pwd)"/"sparcfire-out"

if [[ ! $in_dir || ! $out_dir ]]; then
    echo "No input or output directories supplied as cmd line arguments."
    echo "Using sparcfire-in and sparcfire-out..."
    in_dir=$default_in
    out_dir=$default_out

    if [[ ! -d $in_dir || ! -d $out_dir ]]; then
        echo "Cannot find either $in_dir and/or $out_dir. Quitting."
        exit
    fi
fi

pre_galfit_in="pre_galfit-in"
pre_galfit_out="pre_galfit-out"

post_galfit_in="post_galfit-in"
post_galfit_out="post_galfit-out"

if [[ -d $post_galfit_in ]]; then
    echo "Deleting already existing $post_galfit_in in 5 seconds"
    sleep 5
    rm -rf $post_galfit_in
fi

if [[ -d $post_galfit_out ]]; then
    echo "Deleting already existing $post_galfit_out in 5 seconds"
    sleep 5
    rm -rf $post_galfit_out
fi

# Be careful about restarting
mv $in_dir $pre_galfit_in
mv $out_dir $pre_galfit_out

cp $pre_galfit_out/"galaxy.csv" $pre_galfit_out/"pre_galfit_galaxy.csv"
cp $pre_galfit_out/"galaxy_arcs.csv" $pre_galfit_out/"pre_galfit_galaxy_arcs.csv"

mkdir -p $default_in $default_tmp $default_out

# Populate input folder with models
python3 "${SPARCFIRE_HOME}/GalfitModule/Utilities/grab_model_from_output.py" $pre_galfit_in $pre_galfit_out $default_in

ext="*.fits"
conv_fits="-convert-FITS "

if [ -x "`/bin/which fitspng 2>/dev/null`" ]; then
    fitspng=$(/bin/which fitspng)
    working_dir=$(pwd)
    
    cd $default_in
    $fitspng *".fits"
    cd $working_dir
    
    ext="*.png"
    conv_fits=""
fi

# Prep for parallel
#input_arr=($(ls "$default_in/"*".fits"))
input_arr=($(find "$default_in" -name "$ext"))
input_count="${#input_arr[@]}"
cpu_count=$(nproc --all)
cpu_count=$(( cpu_count < input_count ? cpu_count : input_count ))

parallel_file="send_to_parallel_sparcfire"

# Ceiling division
per_cpu=$(( (input_count/cpu_count)+(input_count%cpu_count>0) ))
for (( cpu_num=0; cpu_num<$cpu_count; ++cpu_num )); do
    new_dir="${default_in}_${cpu_num}"
    mkdir -p $new_dir

    # Run sparcfire with defaults on assuming it has already been setup
    # Also no need for star masking
    if [[ $arr_start -lt $input_count ]]; then
        echo "${SPARCFIRE_HOME}/scripts/SpArcFiRe ${conv_fits}-compute-starmask false -ignore-starmask $new_dir $default_tmp $default_out -generateFitQuality 0 -writeBulgeMask 1 -allowArcBeyond2pi 0 -unsharpMaskAmt 8 -useDeProjectStretch 0 -fixToCenter 1 -errRatioThres 3"
    fi
    
    arr_start=$(( $cpu_num*$per_cpu  ))

    for (( idx=$arr_start; idx<$(( $arr_start+$per_cpu )); ++idx )); do

        if [[ $idx -ge $input_count ]]; then
            break 2
        fi

        gfits="${input_arr[$idx]}"
        #gfits="${gfits##*/}"
        cp $gfits $new_dir
    done
done > $parallel_file

#parallel_script="${SPARCFIRE_HOME}/GalfitModule/ParallelDrivers/distrib_slurm"
parallel_script="${SPARCFIRE_HOME}/GalfitModule/ParallelDrivers/parallel"

# RUN SPARCFIRE
echo "Running SpArcFiRe (again) with $cpu_count nodes"
cat $parallel_file | $parallel_script $cpu_count

mv $in_dir $post_galfit_in
mv $out_dir $post_galfit_out

cp $post_galfit_out/"galaxy.csv" $post_galfit_out/"post_galfit_galaxy.csv"
cp $post_galfit_out/"galaxy_arcs.csv" $post_galfit_out/"post_galfit_galaxy_arcs.csv"

mv $pre_galfit_in $in_dir
mv $pre_galfit_out $out_dir

cp $post_galfit_out/"post_sparcfire_galaxy.csv" $out_dir/"post_galfit_galaxy.csv"
cp $post_galfit_out/"post_sparcfire_galaxy_arcs.csv" $out_dir/"post_galfit_galaxy_arcs.csv"

# Cleanup
rm -rf $parallel_file "sparcfire-in_"*
# Sparcfire junk
rm -rf "a" "not" "tty" "SpArcFiRe-stdin.stdin.txt" *"_settings.txt"
