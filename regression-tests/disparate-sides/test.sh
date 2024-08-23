#!/bin/bash
TMP=/tmp/disparate-sides.$$
#trap "/bin/rm -f $TMP" 0 1 2 3 15
cd regression-tests/disparate-sides
../../scripts/most-disparate-side.sh -ebm `ls -Sr *.tsv` 2>&1 | tee $TMP
fgrep -v 'TMPDIR is ' $TMP > $TMP.2 && mv $TMP.2 $TMP

echo "Checking numerical differences..."
../../scripts/ndiff <(fgrep -v p-value correct.out) <(fgrep -v p-value $TMP)
