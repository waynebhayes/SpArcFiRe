#Description:
# SpArcFiRe to Galfit
#Usage:
# - Run on ics.uci.edu servers to access sparcfire
# - Place png/fits files in 
echo "Running SpArcFiRe"
/home/sparcfire/bin/SpArcFiRe sparcfire_input/ sparcfire_temp/ sparcfire_output/
if [ "$?" != "0" ]; then
	echo "SpArcFiRe error, aborting" 1>&2
	exit 1
fi

#Run python conversion script on each sparfire output
echo "Running python sparcfire to galfit conversion script"
python main.py sparcfire_output/galaxy.csv sparcfire_output galaxy_arcs.csv
if [ "$?" != "0" ]; then
	echo "Python error, aborting" 1>&2
	exit 1
fi

#Run galfit on each feedme file
echo "Running GALFIT"
#galfit galfit_input/*.feedme
if [ "$?" != "0" ]; then
	echo "GALFIT error, aborting" 1>&2
	exit 1
fi

#Show original image, sparcfire output with arcs, galfit light curve, and subtraction
echo "Showing images"

