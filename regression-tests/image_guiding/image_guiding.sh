#!/bin/bash
set -e
die() { echo "FATAL ERROR IN IMAGE GUIDING: $@" >&2; exit 1;}

TEST_RESULT=0
for i in $(seq 0 .25 1)
#for i in $(seq 0 .25 0)
    do
    $SPARCFIRE_HOME/scripts/SpArcFiRe -guide_dir $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.guide $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.in/ $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.tmp/ $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.out/$i -imageGuidingThreshold $i > tempImageGuiding.txt;


    if cmp -s "test_data/test.out/$i/galaxy.csv" "test_data/G.out/$i/galaxy.csv"; then
        TEST_RESULT=1
    #else
        #TEST_RESULT=0
    fi;

    #Cleanup data created by test
    rm -rf $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.out/$i/*;
    rm -rf $SPARCFIRE_HOME/regression-tests/image_guiding/test_data/G.tmp/*;

    diff "tempImageGuiding.txt" "$SPARCFIRE_HOME/regression-tests/image_guiding/image_guiding_$i.txt"
done

rm -f *_settings.txt
exit $TEST_RESULT

