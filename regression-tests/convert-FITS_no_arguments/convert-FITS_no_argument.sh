#!/bin/bash

die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}
TEST_RESULT=1

#Run test
mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.tmp
mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.out
$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.out > convertFits_test.txt 2> convertFITS_err.txt

diff <(cut -f1-39,42-150 -d$','  $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.out/galaxy.csv) <(cut -f1-39,42-150 -d$','  $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/test.out/galaxy.csv) > comp.txt

if ! [ -s "comp.txt" ]; then
	TEST_RESULT=0
fi;

#Cleanup data created by test
rm comp.txt
rm -rf  $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.out/*;
rm -rf  $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.tmp/*;

exit $TEST_RESULT
