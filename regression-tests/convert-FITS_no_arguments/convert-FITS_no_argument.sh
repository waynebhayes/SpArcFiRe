#!/bin/bash
die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; return 1;}

#Run test
source $SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.out;

if cmp -s "test_data/test_out/galaxy.csv" "test_data/G.out/galaxy.csv"; then
	printf 'Pass'
else
	printf 'Fail'
fi;

#Cleanup data created by test
rm -rf  $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.out/*;
rm -rf  $SPARCFIRE_HOME/regression-tests/convert-FITS_no_arguments/test_data/G.tmp/*;
