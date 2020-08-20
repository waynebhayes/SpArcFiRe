#!/bin/bash

die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}


mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp
mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out
$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1 > fit_quality_test.txt 2> fit_quality_err.txt
TEST_RESULT=1

diff <(cut -f1-39,42-150 -d$','  $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/galaxy.csv) <(cut -f1-39,42-150 -d$','  $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/test.out/galaxy.csv) > comp.txt

if ! [ -s "comp.txt" ]; then
	TEST_RESULT=0
fi;

#Cleanup data created by test
rm comp.txt
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/*;
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp/*;

echo $TEST_RESULT
exit $TEST_RESULT
