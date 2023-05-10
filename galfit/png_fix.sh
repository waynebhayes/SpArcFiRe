#!/bin/bash

# author: Matthew Portman
# For fixing PNGS when there's an error with fitspng options

# Silencing bash output where necessary
# Was initially only catching error output... now whatever
silent() {
	"$@" >/dev/null
}

# ***************************************************************

# Use for when script runs ON TOP of sparcfire - if it gets integrated
# into Sparcfire, then this will have to change... well a lot will.

# Hardcoding in, tmp, out directory names for now
in_dir=$PWD/sparcfire-in
tmp_dir=$PWD/sparcfire-tmp
out_dir=$PWD/sparcfire-out

spout="sparcfire-out" # for convenient use in loop

# Also writing in default variable value for fitspng
# fitspng_param="1,150"
fitspng_param="0.25,1" # Default values

# 1/19 - added spaces after the '(' make sure that works in test, bash is weird
files=( $in_dir/*.fits ) #($( ls $in_dir/*.fits ))

for item in ${files[@]}
do
	# Removing '.fits'
	gal_path_in="${item%.fits}"
	
	# Grabbing galaxy name 
	gal_name="${gal_path_in##*/}"

	#echo $gal_name

	#Replacing directory name
	gal_path_out="${gal_path_in/-in/-out}"
	#echo $gal_path

	# TODO: GENERATE PSF HERE -- TO DO, text necessary in feedme, need to write separate script to download all psfield files
	# Rather check for PSF here and generate it if its not here... implement as a function

	# Changing to tmp to run fitspng
	gal_path_fits="${tmp_dir}/galfits/${gal_name}_out.fits"

	#echo $gal_path_fitspng

	# Converting fits to png for easy viewing and their residuals
	silent fitspng -fr "$fitspng_param" -o "${gal_name}.png" "${gal_path_fits}[1]"
	silent fitspng -fr "$fitspng_param" -o "${gal_name}_out.png" "${gal_path_fits}[2]"
	silent fitspng -fr "$fitspng_param" -o "${gal_name}_residual.png" "${gal_path_fits}[3]"
	
	# Combining the three with Sparcfire's images using ImageMagick
	#silent montage "${gal_name}.png" "${gal_name}_out.png" "${gal_name}_residual.png" "${spout}/${gal_name}/${gal_name}-A_input.png" "${spout}/${gal_name}/${gal_name}-C_preproc.png" "${spout}/${gal_name}/${gal_name}-J_logSpiralArcs-merged.png" -geometry 150x125+2+4 "${gal_name}_combined.png"
    silent montage "${gal_name}.png" "${gal_name}_out.png" "${gal_name}_residual.png" -tile 3x1 -geometry "150x150+2+0<" "${gal_name}_combined.png"

done

# Cleaning up 
silent rm fit.log #rm galfit.*
mv *_combined.png $out_dir/all_galfit_out/galfit_png/
rm *.png
