#!/bin/bash

#This is meant to clean up galfit files created when using the generateFitQuality option.
#Galfit will output a series of files, called galfit.01, galfit.02, ... etc
#This script simply moves the file to the corresponding galaxy directory

#TODO move fit.log
#TODO check if the generateFitQuality variable has been set so that this is safer

BASEOUTDIR=$1

echo "Moving Galfit Files"

for file in $(find "$SPARCFIRE_HOME/scripts" -iregex '.*/galfit\.[0-9]+')
do
    secondLine=$(sed '2q;d' $file)
    out=$(echo "$secondLine" | cut -c 1-20 --complement )
    outdir=$(echo "$out" | rev | cut -d '/' -f2- | rev)
    if [[ -d "$outdir" ]]
    then
        outpath=$(echo "$out" | rev | cut -c 1-7 --complement | rev)
        outpath="$outpath-galfit_control_params"
        mv $file $outpath
    fi
done

mv fit.log "$BASEOUTDIR/galfit_fit.log"
