#!/bin/bash
die() { echo "FATAL ERROR IN IMAGE GUIDING: $@" >&2; exit 1;}

abs(){
	return 'scale=8;sqrt($1 - $2) ^ 2)' | bc
}


TEST_RESULT=1
for i in $(seq 0 .25 1)
#for i in $(seq 0 .25 0)
    do
    mkdir -p $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.tmp/
    mkdir -p $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.out/$i
    $SPARCFIRE_HOME/scripts/SpArcFiRe -guide_dir $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.guide $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.tmp/ $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.out/$i -imageGuidingThreshold $i > tempImageGuiding.txt;

    cut -f11-22,38-40,44-50,51-57,62-83,86-162 -d$','  $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.out/"$i"/galaxy.csv | sed '1d' > /tmp/test_numbers1
    cut -f11-22,38-40,44-50,51-57,62-83,86-162 -d$','  $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/test.out/"$i"/galaxy.csv | sed '1d' > /tmp/test_numbers2

    sed 's/[^0-9.EefF-]/ /g' /tmp/test_numbers1 > /tmp/sPaRCfIrE_TC1
    sed 's/[^0-9.EefF-]/ /g' /tmp/test_numbers2 > /tmp/sPaRCfIrE_TC2

    paste /tmp/sPaRCfIrE_TC1 /tmp/sPaRCfIrE_TC2 > /tmp/comp1

    awk -F'[-,]' 'function abs(x) {return x < 0 ? -v : v} {for(i=1;i<=NF;i++) if(abs($1-$2) > 1.0e-6) print "abs($1-$2"}' /tmp/comp1 > /tmp/sPaRCfIrE_TCP1

    if ! [ -s "comp.txt" ]; then
	    
    #if cmp -s "test_data/test.out/$i/galaxy.csv" "test_data/G.out/$i/galaxy.csv"; then
        TEST_RESULT=0
    else
        TEST_RESULT=1
    fi;

    #Cleanup data created by test
    rm -rf $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.out/$i/*;
    rm -rf $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.tmp/*;

done
echo $TEST_RESULT
rm -f *_settings.txt
exit $TEST_RESULT

