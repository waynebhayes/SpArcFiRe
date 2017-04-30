#Description:
# SpArcFiRe to Galfit
#Usage:
# - Run on ics.uci.edu servers to access sparcfire
# - Place png/fits files in sparcfire_input
# - Requires fv (fits viewer) to view output automatically

#Run SpArcFiRe on images within in/ directory
echo "Running SpArcFiRe"
#rm -rf tmp/* out/* 
#/home/sparcfire/bin/SpArcFiRe -convert-FITS in tmp out
if [ "$?" != "0" ]; then
	echo "SpArcFiRe error, aborting" 1>&2
	exit 1
fi


#Run python conversion script on each sparfire output
#echo "Running python sparcfire to galfit conversion script"
rm -rf galfit_in/* galfit_out/*
python2.7 s2g.py

#Run galfit on each feedme file
echo "Running GALFIT"
for f in galfit_in/*
do
	galfit "$f"
done

#un python comparison and confidence estimation script
python2.7 compare_galfit.py

#Show original image, sparcfire output with arcs, galfit light curve, and subtraction
#echo "Showing images"
#fv galfit_out/*
