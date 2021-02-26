#!/bin/bash

die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}

TEST_DIR=$SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data
mkdir -p $TEST_DIR/G.tmp
mkdir -p $TEST_DIR/G.out

# Cleanup data created by test
rm -rf $TEST_DIR/G.out/*;
rm -rf $TEST_DIR/G.tmp/*;

$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $TEST_DIR/G.in $TEST_DIR/G.tmp $TEST_DIR/G.out -generateFitQuality 1 > fit_quality_test.txt 2> fit_quality_err.txt
TEST_RESULT=1

diff <(cut -f1-39,42-150 -d$','  $TEST_DIR/G.out/galaxy.csv) <(cut -f1-39,42-150 -d$','  $TEST_DIR/Gcorrect.out/galaxy.csv) > comp.txt

if ! [ -s "comp.txt" ]; then
	TEST_RESULT=0
fi;

echo $TEST_RESULT
exit $TEST_RESULT
