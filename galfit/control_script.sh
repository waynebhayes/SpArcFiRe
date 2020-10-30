#!/bin/bash

# Author: Matthew Portman
# Date (see GITHUB, I may forget to change this): 10/22/20

# For controlling galfitting via sparcfire

# Silencing bash output where necessary
silent() {
	"$@" 2>/dev/null
}

# Insert code to check for directory structure and exit

# Making some directories and cleaning them

silent mkdir ./sparcfire-tmp/galfit_masks \
./sparcfire-tmp/galfits \
./sparcfire-temp/PSF_files \
./sparcfire-out/all_galfit_out \
./sparcfire-out/all_galfit_out/galfit_png 

silent rm -f \ #./sparcfire-tmp/galfit_masks/* \
./sparcfire-tmp/galfits/*.fits \
./sparcfire-out/all_galfit_out/*.fits \
./sparcfire-out/all_galfit_out/galfit_png/*.png \


if which galfit; then
	echo -e "Found GALFIT. Proceeding.\n"
	run_galfit=$(which galfit)
else
	echo -e "Can't call GALFIT. Try 'which galfit' " && exit
fi

if which fitspng; then
	echo -e "Found fitspng. Proceeding.\n"
else
	echo -e "No fitspng found, pngs will not be automatically generated.\n"
	echo -e "You can find information on fitspng at http://integral.physics.muni.cz/fitspng/\n"
fi

# Running galfit feedme generator featuring fun error handling
feedme_gen=sparc_to_galfit_feedme_gen.py
feed="Running feedme generator...\n"
if silent python $feedme_gen; then
	echo -e ${feed}
	python=python
elif silent python3 $feedme_gen; then
        echo -e ${feed}
	python=python3
elif silent python2.7 $feedme_gen; then 
	echo -e ${feed}
	python=python2
else
	echo -e "Python couldn't run OR no usable Python (default, 2.7, 3). Proceeding with previous output." # && exit 
fi

# Relative file paths for now

files=(./sparcfire-in/*.fits) #($( ls ./sparcfire-in/*.fits ))
masks=(./sparcfire-tmp/galfit_masks/*star-rm.fits) #($( ls ./sparcfire-tmp/galfit_masks/ ))

# Generate Star masks here if not already done so
if [ ${#masks[@]} -ne ${#files[@]} ]; then
	
	echo "Star Masking"
	cd star_removal # Unfortunately sextractor does not seem to play nice with a single directory so unfortunately we must cd
	echo "running sextractor"
	$python remove_stars_with_sextractor.py ../sparcfire-in/ ../sparcfire-tmp/galfit_masks/
	cd ..
else
	echo -e "Star Masks already ran. Proceeding.\n"
fi

spout="sparcfire-out"

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

	# GENERATE PSF HERE -- TO DO, text necessary in feedme, need to write separate script to download all psfield files
	# Rather check for PSF here and generate it if its not here... implement as a function

	echo "Galfitting" $feedme_path
	$run_galfit $feedme_path

	# Changing to tmp to run fitspng
	gal_path_fits="./sparcfire-tmp/galfits/${gal_name}_out.fits"

	#echo $gal_path_fitspng

	# Converting fits to png for easy viewing and their residuals
	silent fitspng -fr "1,150" -o "${gal_name}.png" "${gal_path_fits}[1]"
	silent fitspng -fr "1,150" -o "${gal_name}_out.png" "${gal_path_fits}[2]"
	silent fitspng -fr "1,150" -o "${gal_name}_residual.png" "${gal_path_fits}[3]"
	
	# Combining the three with Sparcfire's images using ImageMagick
	silent montage "${gal_name}.png" "${gal_name}_out.png" "${gal_name}_residual.png" "${spout}/${gal_name}/${gal_name}-A_input.png" "${spout}/${gal_name}/${gal_name}-C_preproc.png" "${spout}/${gal_name}/${gal_name}-J_logSpiralArcs-merged.png" -geometry 150x125+2+4 "${gal_name}_combined.png"

	# Moving successful fit into output folder
	cp "$gal_path_fits" "${gal_path_out}/"

	# Moving galfit output into output folder
	mv galfit.01 "${gal_path_out}/"

done

# For future use in analysis
cp ./sparcfire-tmp/galfits/*.fits ./sparcfire-out/all_galfit_out/

# Cleaning up 
silent rm galfit.* fit.log
mv *_combined.png ./sparcfire-out/all_galfit_out/galfit_png/
rm *.png

# Running a script to compare input to galfit and output
# See comparison_params.csv for *just* the differences
# Each galaxy folder contains the input, output, and difference in a text file galfit_io_compare
$python in_out_comparison.py

# For running fitspng on all galfit output assuming all in one folder
# Keep the below line in case of desire to parallelize... which will be strong
# ls sparcfire-out/all_fits/*.fits | xargs printf -- '%s[2] ' | xargs -n 1 fitspng
# can use xargs -P number where the number is the number of processes
# In fact, I could parallelize the whole thing if I run galfit as a separate process via xargs
# on such small images, it shouldn't be so bad. 

# For resizing images...
# for i in *; do convert $i -resize 200% "./larger_images/${i%.png}_large.png"; done
