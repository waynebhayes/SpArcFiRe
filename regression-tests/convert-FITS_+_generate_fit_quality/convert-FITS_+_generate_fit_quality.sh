#!/bin/bash

die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}

TEST_DIR=$SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data
mkdir -p $TEST_DIR/G.tmp
mkdir -p $TEST_DIR/G.out

# Cleanup data created by test
rm -rf $TEST_DIR/G.out/*;
rm -rf $TEST_DIR/G.tmp/*;

$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $TEST_DIR/G.in $TEST_DIR/G.tmp $TEST_DIR/G.out -generateFitQuality 1 > fit_quality_test.txt 2> fit_quality_err.txt

csv2tsv $TEST_DIR/G.out/galaxy.csv $TEST_DIR/Gcorrect.out/galaxy.csv
if regression-diff.sh $TEST_DIR/G.out/galaxy.tsv $TEST_DIR/Gcorrect.out/galaxy.tsv; then
    echo SUCCESS
else
    echo FAIL; exit 1
fi
