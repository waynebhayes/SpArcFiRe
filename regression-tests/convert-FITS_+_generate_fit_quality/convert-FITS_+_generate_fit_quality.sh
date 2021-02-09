#!/bin/bash

die() { echo "FATAL ERROR IN CONVERT FITS: $@" >&2; exit 1;}

abs(){
	return 'scale=8;sqrt($1 - $2) ^ 2)' | bc
}

mkdir -p $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out
$SPARCFIRE_HOME/scripts/SpArcFiRe -convert-FITS $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out -generateFitQuality 1 > fit_quality_test.txt 2> fit_quality_err.txt
TEST_RESULT=1

sed 's/[^0-9.EefF-]/ /g' $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/galaxy.csv > /tmp/sPaRCfIrE_TC1
sed 's/[^0-9.EefF-]/ /g' $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/test.out/galaxy.csv > /tmp/sPaRCfIrE_TC2

paste /tmp/sPaRCfIrE_TC1 /tmp/sPaRCfIrE_TC2 > /tmp/sPaRCfIrE_TCP1

awk '{for(i=1;i<=NF;i++) if(abs($1 $2) > 1.0e-6) print "different"}' > /tmp/sPaRCfIrE_TCP1 /tmp/comp.txt

if ! [ -s "/tmp/comp.txt" ]; then
	TEST_RESULT=0
fi;

#Cleanup data created by test
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.out/*;
rm -rf $SPARCFIRE_HOME/regression-tests/convert-FITS_+_generate_fit_quality/test_data/G.tmp/*;

echo $TEST_RESULT
exit $TEST_RESULT
