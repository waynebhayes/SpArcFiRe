#!/bin/bash
die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; return 1;}

#Run test$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCVIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCVIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCVIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1
source $SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1;

#Cleanup log files
rm $cwd/*settings.txt;
rm $cwd/*stdin.txt;

if cmp -s "test_data/test_out/galaxy.csv" "test_data/G.out/galaxy.csv"; then
	printf 'Pass'
else
	printf 'Fail'
fi;

#Cleanup data created by test
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/*;
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp/*;
