#!/bin/bash

die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}

abs(){
	return 'scale=8;sqrt($1 - $2) ^ 2)' | bc
}

mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp
mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out
$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1 > fit_quality_test.txt 2> fit_quality_err.txt
TEST_RESULT=1

cut -f11-22,38-40,44-50,51-57,62-83,86-162 -d$','  $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/galaxy.csv | sed '1d' > /tmp/test_numbers1
cut -f11-22,38-40,44-50,51-57,62-83,86-162 -d$','  $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/test.out/galaxy.csv | sed '1d' > /tmp/test_numbers2

sed 's/[^0-9.EefF-]/ /g' /tmp/test_numbers1 > /tmp/sPaRCfIrE_TC1
sed 's/[^0-9.EefF-]/ /g' /tmp/test_numbers2 > /tmp/sPaRCfIrE_TC2

paste /tmp/sPaRCfIrE_TC1 /tmp/sPaRCfIrE_TC2 > /tmp/comp1

awk -F'[-,]' 'function abs(x) {return x < 0 ? -v : v} {for(i=1;i<=NF;i++) if(abs($1-$2) > 1.0e-6) print "abs($1-$2"}' /tmp/comp1 > /tmp/sPaRCfIrE_TCP1


if ! [ -s "/tmp/sPaRCfIrE_TCP1" ]; then
	TEST_RESULT=0
fi;

#Cleanup data created by test
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/*;
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp/*;

echo $TEST_RESULT
exit $TEST_RESULT
