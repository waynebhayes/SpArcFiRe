#!/bin/bash
die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}

#Run test$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCVIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCVIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCVIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1
source $SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1;
TEST_RESULT=0

if cmp -s "test_data/test_out/galaxy.csv" "test_data/G.out/galaxy.csv"; then
	TEST_RESULT=1
else
	TEST_RESULT=0
fi;

#Cleanup data created by test
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/*;
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp/*;

exit $TEST_RESULT
