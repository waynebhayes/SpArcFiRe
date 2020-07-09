#!/bin/bash
set -e
die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}
$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1 > fit_quality_test.txt 2> fit_quality_err.txt
TEST_RESULT=0

if cmp -s "$SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/test_out/galaxy.csv" "$SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/galaxy.csv"; then
	TEST_RESULT=1
else
	TEST_RESULT=0
fi;


#Cleanup data created by test
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/*;
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp/*;

diff $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/test_out/galaxy.csv $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/galaxy.csv > fit_quality_out.txt

exit $TEST_RESULT
