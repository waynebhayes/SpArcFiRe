#!/bin/bash

# Author: Matthew Portman
# Date (see GITHUB, I may forget to change this): 1/26/21

# 1/26/21
# New usage ---
# Will default to running on top of sparc directory 
# unless otherwise specified
# 
# Command line args
#
# -r, --rerun-galfit
#		Re-run GALFIT using the galfit.01 file
#
# -p, --path
#		Specify path to Sparcfire's in/tmp/out directories in that order;
#		e.g. '-p /home/sparcfire_in /home/sparcfire_tmp /home/sparcfire_out'
#

# Change log - 1/19/21 
# Updated silent to catch everything
# Updated call to feedme gen to include paths
# Changed 'which' to 'type -P' which seems better
# See note for files= and masks=

# For controlling galfitting via sparcfire

# ***************************************************************
USAGE="USAGE:

	bash ./control_script.sh [-r|--rerun-galfit] [-p|--path] IN-DIRECTORY TMP-DIRECTORY OUT-DIRECTORY
	
	This script is the wrapping script for running GALFIT using SpArcFiRe to inform 
	the input. By default, it runs from the directory it is called and checks for the
	existence (verbatim) of the 'sparcfire-in' 'sparcfire-tmp' and 'sparcfire-out' 
	directories in the current directory. If this is not the case, it will fail and exit. 
	
	You may specify	a different in, tmp, and out directory using the '-p' option. Please
	do not specify a symlink, it discomforts the programmer.
"

die() { (echo "$USAGE"; echo "FATAL ERROR: $@") >&2; return 1
}

# Silencing bash output where necessary
# Was initially only catching error output... now whatever
silent() {
	"$@" >/dev/null
}

# Insert code to check for directory structure and exit
# ***************************************************************
# INSERTING CODE FOR GENERALIZING DIRECTORY STRUCTURE - 1/6/21

# Use for when script runs ON TOP of sparcfire - if it gets integrated
# into Sparcfire, then this will have to change... well a lot will.

# Hardcoding in, tmp, out directory names for now
in_dir=$PWD/sparcfire-in
tmp_dir=$PWD/sparcfire-tmp
out_dir=$PWD/sparcfire-out

spout="sparcfire-out" # for convenient use in loop

# ***************************************************************

# Taking command line arguments
# Code nabbed from stackoverflow
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r|--rerun-galfit)
	rerun=1
    shift # past argument
    ;;
    -p|--path)
	# Checking if empty... just in case
	[ -z "$2" ] && echo "No in_dir given, using default..." || in_dir="$2"
	[ -z "$3" ] && echo "No in_dir given, using default..." || tmp_dir="$3"
	[ -z "$4" ] && echo "No in_dir given, using default..." || out_dir="$4"
    shift # past argument
    shift # past in
    shift # past tmp
    shift # past out
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

set -- "${POSITIONAL[@]}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi

# ***************************************************************

# Checking for things necessary to run, GALFIT, directories specified, etc.

if [ ! -d "$in_dir" ] || [ ! -d "$tmp_dir" ] || [ ! -d "$out_dir" ]; then 
    die "Couldn't find SpArcFiRe's in, tmp, out directories. Exitting."
fi

if silent type -P galfit; then
	echo -e "Found GALFIT. Proceeding.\n"
	run_galfit="$( type -P galfit )" #$( which galfit )
else
    die "Can't call GALFIT by 'type -P galfit'. Exitting."
fi

# Updated to use type - 1/19
if silent type -P python3; then
	python=python3
elif silent type -P python2.7; then
	python=python2.7
    echo -e "Python version is 2.7. Feedme gen will NOT be able to run (currently needs >3.6). Proceeding anyway.\n"
elif silent type -P python; then 
	python=python
else
	die "Python couldn't run OR no usable Python (default, 2.7, 3). Exitting."
fi

echo -e "Using $python\n"

if silent type -P fitspng; then
	echo -e "Found fitspng. Proceeding.\n"
else
	echo -e "No fitspng found, pngs will not be automatically generated.\n"
	echo -e "You can find information on fitspng at http://integral.physics.muni.cz/fitspng/\n"
fi

# Also writing in default variable value for fitspng
# fitspng_param="1,150"
fitspng_param="0.25,1" # Default values

# ***************************************************************

# Making some directories and cleaning them

silent mkdir $tmp_dir/galfit_masks \
$tmp_dir/galfits \
$tmp_dir/psf_files \
$out_dir/all_galfit_out \
$out_dir/all_galfit_out/galfit_png 

silent rm -f \ #$tmp_dir/galfit_masks/* \
$tmp_dir/galfits/*.fits \
#$out_dir/all_galfit_out/*.fits \
#$out_dir/all_galfit_out/galfit_png/*.png \

# Defining the max number of iterations
max_it=150

# Running galfit feedme generator
feedme_gen="sparc_to_galfit_feedme_gen.py"
feed="Running feedme generator...\n"
echo -e ${feed}
$python $feedme_gen $in_dir $tmp_dir $out_dir # Now uses the directory as input on command line

# 1/19 - added spaces after the '(' make sure that works in test, bash is weird
files=( $in_dir/*.fits ) #($( ls $in_dir/*.fits ))
masks=( $tmp_dir/galfit_masks/*star-rm.fits ) #($( ls $tmp_dir/galfit_masks/ ))

#echo $files
# Generate Star masks here if not already done so
#echo "masks ${#masks[@]}"
#echo "files ${#files[@]}"

if [ ${#masks[@]} -ne ${#files[@]} ]; then
	
	echo "Star Masking"
	cd star_removal # Unfortunately sextractor does not seem to play nice with a single directory so unfortunately we must cd
	echo "running sextractor"
    # IF THIS FAILS IT'S LIKELY BECAUSE THE DIRECTORIES NEED TO BE ABSOLUTE PATHS
	# TODO: CONFIRM THIS
	$python remove_stars_with_sextractor.py $in_dir/ $tmp_dir/galfit_masks/
	cd ..
else
	echo -e "Star Masks already ran. Proceeding.\n"
fi

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

	#Modifying path to feedme to run GALFIT
	feedme_path="${gal_path_out}/autogen_feedme_galfit.in"

	# TODO: GENERATE PSF HERE -- TO DO, text necessary in feedme, need to write separate script to download all psfield files
	# Rather check for PSF here and generate it if its not here... implement as a function

	echo "Galfitting $feedme_path with" "$run_galfit"
	$run_galfit -imax "$max_it" "$feedme_path"
	
	# Adding re-run if chosen 
	
	if [ ! -z "$rerun" ] && [ -f "./galfit.01" ]; then
	    echo "Re-running GALFIT..."
	    $run_galfit "${gal_path_out}/galfit.01" #$run_galfit -imax "$max_it" "${gal_path_out}/galfit.01"
	fi

	# Changing to tmp to run fitspng
	gal_path_fits="${tmp_dir}/galfits/${gal_name}_out.fits"

	#echo $gal_path_fitspng

	# Converting fits to png for easy viewing and their residuals
    # TODO: Drop these in tmp dir
	silent fitspng -fr "$fitspng_param" -o "${gal_name}.png" "${gal_path_fits}[1]"
	silent fitspng -fr "$fitspng_param" -o "${gal_name}_out.png" "${gal_path_fits}[2]"
	silent fitspng -fr "$fitspng_param" -o "${gal_name}_residual.png" "${gal_path_fits}[3]"
	
	# Adding these lines for easy comment/uncomment if fitspng won't cooperate
	# using the previous settings
	#silent fitspng -o "${gal_name}.png" "${gal_path_fits}[1]"
    #silent fitspng -o "${gal_name}_out.png" "${gal_path_fits}[2]"
    #silent fitspng -o "${gal_name}_residual.png" "${gal_path_fits}[3]"
	
	# Combining the three with Sparcfire's images using ImageMagick
	silent montage "${gal_name}.png" "${gal_name}_out.png" "${gal_name}_residual.png" "${spout}/${gal_name}/${gal_name}-A_input.png" "${spout}/${gal_name}/${gal_name}-C_preproc.png" "${spout}/${gal_name}/${gal_name}-J_logSpiralArcs-merged.png" -geometry 150x125+2+4 "${gal_name}_combined.png"

	# Moving successful fit into output folder
	cp $gal_path_fits "${gal_path_out}/"

	# Moving galfit output into output folder
	mv galfit.* "${gal_path_out}/"

done

# For future use in analysis
cp $tmp_dir/galfits/*.fits $out_dir/all_galfit_out/

# Cleaning up 
silent rm galfit.* fit.log
mv *_combined.png $out_dir/all_galfit_out/galfit_png/
rm *.png

# Running a script to compare input to galfit and output
# See comparison_params.csv for *just* the differences
# Each galaxy folder contains the input, output, and difference in a text file galfit_io_compare
echo "Running in_out_comparison.py"
# TODO: This needs to be updated
$python in_out_comparison.py $in_dir $tmp_dir $out_dir

# For running fitspng on all galfit output assuming all in one folder
# Keep the below line in case of desire to parallelize... which will be strong
# ls $out_dir/all_fits/*.fits | xargs printf -- '%s[2] ' | xargs -n 1 fitspng
# can use xargs -P number where the number is the number of processes
# In fact, I could parallelize the whole thing if I run galfit as a separate process via xargs
# on such small images, it shouldn't be so bad. 

# For resizing images...
# for i in *; do convert $i -resize 200% "./larger_images/${i%.png}_large.png"; done
