#!/bin/bash
die() { echo "FATAL ERROR IN IMAGE GUIDING: $@" >&2; exit 1;}

TEST_RESULT=1
TEST_DIR=$SPARCFIRE_HOME/regression-tests/image_guiding/test_data
for i in $(seq 0 .25 1)
#for i in $(seq 0 .25 0)
    do
    mkdir -p $TEST_DIR/G.tmp/
    mkdir -p $TEST_DIR/G.out/$i

    # Cleanup before test runs
    rm -rf $TEST_DIR/G.out/$i/*;
    rm -rf $TEST_DIR/G.tmp/*;

    $SPARCFIRE_HOME/scripts/SpArcFiRe -guide_dir $TEST_DIR/G.guide $TEST_DIR/$TEST_DIR/G.tmp $TEST_DIR/G.out/$i -imageGuidingThreshold $i > tempImageGuiding.txt;

    diff <(cut -f1-39,42-150 -d$','  $TEST_DIR/G.out/"$i"/galaxy.csv) <(cut -f1-39,42-150 -d$','  $TEST_DIR/Gcorrect.out/$i/galaxy.csv") > comp.txt
    if ! [ -s "comp.txt" ]; then
	    
    #if cmp -s "test_data/Gcorrect.out/$i/galaxy.csv" "test_data/G.out/$i/galaxy.csv"; then
        TEST_RESULT=0
    else
        TEST_RESULT=1
    fi;

done
echo $TEST_RESULT
rm -f *_settings.txt
exit $TEST_RESULT

