#Description:
# SpArcFiRe to Galfit
#Usage:
# - Run on ics.uci.edu servers to access sparcfire
# - Place png/fits files in sparcfire_input
# - Requires fv (fits viewer) to view output automatically

#Run SpArcFiRe on images within in/ directory
echo "Running SpArcFiRe"
#rm -rf tmp/* out/* 
#/home/sparcfire/bin/SpArcFiRe -convert-FITS in tmp out -clusSizeCutoff 25 -unsharpMaskAmt 8 -unsharpMaskSigma 20 -stopThres 0.0500 -mergeChkMinClusSz 350 -errRatioThres 7.5000 -lookForBulge 1
if [ "$?" != "0" ]; then
	echo "SpArcFiRe error, aborting" 1>&2
	exit 1
fi


#Run python conversion script on each sparfire output
#echo "Running python sparcfire to galfit conversion script"
#rm -rf galfit_in/* galfit_out/*
#python2.7 s2g.py

#activate ? (Will added this on 5/7)

python2.7 compare.py out
#if [ "$?" != "0" ]; then
#	echo "Python error, aborting" 1>&2
#	exit 1
#fi

#Run galfit on each feedme file
#echo "Running GALFIT"
#for f in galfit_in/*
#do
#	galfit "$f"
#	if [ "$?" != "0" ]; then
#		echo "GALFIT error, aborting" 1>&2
#		exit 1
#	fi
#done

#Show original image, sparcfire output with arcs, galfit light curve, and subtraction
#echo "Showing images"
#fv galfit_out/*
