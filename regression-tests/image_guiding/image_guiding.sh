#!/bin/bash
die() { echo "FATAL ERROR IN IMAGE GUIDING: $@" >&2; exit 1;}

FAIL_NUM=0
TEST_DIR=$SPARCFIRE_HOME/regression-tests/image_guiding/test_data
for i in $(seq 0 .25 1)
    do
    mkdir -p $TEST_DIR/G.tmp/
    mkdir -p $TEST_DIR/G.out/$i

    # Cleanup before test runs
    rm -rf $TEST_DIR/G.out/$i/*;
    rm -rf $TEST_DIR/G.tmp/*;

    $SPARCFIRE_HOME/scripts/SpArcFiRe -guide_dir $TEST_DIR/G.guide $TEST_DIR/$TEST_DIR/G.tmp $TEST_DIR/G.out/$i -imageGuidingThreshold $i > tempImageGuiding.txt;

	       csv2tsv $TEST_DIR/G.out/"$i"/galaxy.csv $TEST_DIR/Gcorrect.out/$i/galaxy.csv
    regression-diff.sh $TEST_DIR/G.out/"$i"/galaxy.tsv $TEST_DIR/Gcorrect.out/$i/galaxy.tsv || ((++FAIL_NUM))
done
rm -f *_settings.txt
exit $FAIL_NUM

