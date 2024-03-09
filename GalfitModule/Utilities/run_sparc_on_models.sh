#!/bin/bash

basename=$1
if [[ ! $basename ]]; then
    echo "Please provide a basename which matches the runname used."
    echo "The default is GALFIT_output. *NOT proceeding on that assumption."
    exit
fi

in_dir=$2
tmp_dir=$3
out_dir=$4

default_in="$(pwd)"/"sparcfire-in"
default_tmp="$(pwd)"/"sparcfire-tmp"
default_out="$(pwd)"/"sparcfire-out"

if [[ ! $in_dir || ! $tmp_dir || ! $out_dir ]]; then
    echo "No input, temp, or output directories supplied as cmd line arguments."
    echo "Using sparcfire-in, sparcfire-tmp, and sparcfire-out..."
    in_dir=$default_in
    tmp_dir=$default_tmp
    out_dir=$default_out

    if [[ ! -d $in_dir || ! -d $tmp_dir ||  ! -d $out_dir ]]; then
        echo "Cannot find any of: $in_dir, $tmp_dir, $out_dir. Quitting."
        exit
    fi
fi

pre_galfit_in=$in_dir
pre_galfit_out=$out_dir

basename_dir="$pre_galfit_out"/"$basename"

if [[ ! -d $basename_dir ]]; then 
    echo "Cannot find $basename_dir. Did you specify the right name? Quitting."
    exit
fi

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

# Assume galaxy.csv corresponds to the results that were just used to create models
cp $pre_galfit_out/"galaxy.csv" $pre_galfit_out/"${basename}"/"${basename}_pre_galfit_galaxy.csv"
cp $pre_galfit_out/"galaxy_arcs.csv" $pre_galfit_out/"${basename}"/"${basename}_pre_galfit_galaxy_arcs.csv"

mkdir -p $post_galfit_in $post_galfit_out

# Populate input folder with models
# True indicates to flip u/d the models to keep consistent with SpArcFiRe's processing
#python3 "${SPARCFIRE_HOME}/GalfitModule/Utilities/extract_model_from_galfit_output.py" $pre_galfit_in $pre_galfit_out $default_in "true"
echo "Populating new input folder with arm-only models... (this may take awhile)"
for gfits in "$tmp_dir"/"galfits"/*"_for_sparcfire.fits"; do   
    gfits_base="${gfits##*/}"
    cp $gfits "$post_galfit_in"/"${gfits_base/_for_sparcfire/}"
done

# Set up additional command line inputs to SpArcFiRe
ext="*.fits"
conv_fits="-convert-FITS "
elps_dir="$(pwd)"/"elps-dir"
elps=""
if [[ -d $elps_dir ]]; then
    elps="-elps_dir $elps_dir "
    
    # Since the GALFIT files are already cropped, we don't want 
    # SpArcFiRe to re-crop them.
    # Don't need this after all
    # for efile in "$elps_dir"/*"_elps.txt"; do
    #     sed -i s/"cropRad=[0-9]*"/"cropRad=0"/g "$efile"
    # done
fi
#img_standardize="1"

# if [ -x "`/bin/which fitspng 2>/dev/null`" ]; then
#     fitspng=$(/bin/which fitspng)
#     working_dir=$(pwd)
    
#     cd $default_in
#     #$fitspng *".fits"

#     # Process all image files
#     for f in *".fits"; do

#         png_galaxy="${f%.fits}.png"

#         #convert "$f" "$png_galaxy"
    
#         # Ignore case, and suppress errors if no files 
#         shopt -s nullglob
#         shopt -s nocaseglob

#         # Get image's width and height, in one go
#         # Do I need this???
#         read w h < <(identify -format "%w %h" "$f")
        
#         if [ $((w % 2)) -eq 1 ]
#         then
#             w=$((++w))
#         fi

#         if [ $((w/2 % 2)) -eq 1 ]
#         then
#             w=$((w+2))
#         fi
        
#         convert "$f" "-resize" "${w}x${w}" "$png_galaxy"

#     done

#     cd $working_dir
    
#     ext="*.png"
#     conv_fits=""
#     #img_std="0"
# fi

# Prep for parallel
# Determine how to split up input for processing
echo "Preparing to run SpArcFiRe with distributed computing"
#input_arr=($(ls "$default_in/"*".fits"))
input_arr=($(find "$post_galfit_in" -name "$ext"))
input_count="${#input_arr[@]}"
cpu_count=$(nproc --all)
cpu_count=$((cpu_count / 4)) # Since MATLAB uses multithreading (3)
cpu_count=$(( cpu_count < input_count ? cpu_count : input_count ))

parallel_file="send_to_parallel_sparcfire"

# Make directories for parallelizing and write calls to SpArcFiRe to file
# Ceiling division
per_cpu=$(( (input_count/cpu_count)+(input_count%cpu_count>0) ))
for (( cpu_num=0; cpu_num<$cpu_count; ++cpu_num )); do
    new_dir="${default_in}_${cpu_num}"
    mkdir -p $new_dir

    # Run sparcfire with defaults on assuming it has already been setup
    # Also no need for star masking
    # Pad images to even so that we can turn off image standardization
    # elps and conv_fits have spaces in them already
    if [[ $arr_start -lt $input_count ]]; then
        echo "${SPARCFIRE_HOME}/scripts/SpArcFiRe ${conv_fits}-compute-starmask false -ignore-starmask ${elps}$new_dir $default_tmp $post_galfit_out -generateFitQuality 0 -allowArcBeyond2pi 0 -unsharpMaskAmt 0 -useDeProjectStretch 1 -fixToCenter 0 -medFiltRad 0 -useImageStandardization 1" # -errRatioThres 2.8"
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
#cat "$parallel_file" | "nice" "-19" "$parallel_script" "SPARCFIRE_ON_GALFIT" "-M" "all"
cat "$parallel_file" | "nice" "-19" "$parallel_script" "$cpu_count"

# Rename results to include basename
cp $post_galfit_out/"galaxy.csv" $post_galfit_out/"${basename}_post_galfit_galaxy.csv"
cp $post_galfit_out/"galaxy_arcs.csv" $post_galfit_out/"${basename}_post_galfit_galaxy_arcs.csv"

# Copy results to basename output directory
cp $post_galfit_out/"${basename}_post_galfit_galaxy.csv" "$basename_dir"/"${basename}_post_galfit_galaxy.csv"
cp $post_galfit_out/"${basename}_post_galfit_galaxy_arcs.csv" "$basename_dir"/"${basename}_post_galfit_galaxy_arcs.csv"

# Cleanup
rm -rf "$parallel_file" "${default_in}_"*
# Sparcfire junk
rm -rf "a" "not" "tty" "SpArcFiRe-stdin.stdin.txt" *"_settings.txt"
