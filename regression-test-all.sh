#!/bin/bash
case "$1" in
-use-git-at)
    if [ -f git-at ] && [ `wc -l < git-at` -eq 2 -a `git log -1 --format=%at` -eq `tail -1 git-at` ]; then
	echo -n "Repo unchanged; returning same status code as "
	tail -1 git-at | xargs -I{} date -d @{} +%Y-%m-%d-%H:%M:%S
	exit `head -1 git-at`
    fi
    shift
    ;;
esac

USAGE="USAGE: $0 [ list of tests to run, defaults to regression-tests/*/*.sh ]"
source ~/bin/misc.sh
PATH=`pwd`:`pwd`/scripts:$PATH
export PATH

#SpArcFiRe ONLY!!
source setup.sh

NUM_FAILS=0
STDBUF=''
if which stdbuf >/dev/null; then
    STDBUF='stdbuf -oL -eL'
fi
if [ $# -eq 0 ]; then
    set regression-tests/*/*.sh
fi
for r
do
    REG_DIR=`dirname "$r"`
    NEW_FAILS=0
    export REG_DIR
    echo --- running test $r ---
    if eval $STDBUF "$r"; then # force output and error to be line buffered
	:
    else
	NEW_FAILS=$?
	(( NUM_FAILS+=$NEW_FAILS ))
    fi
    echo --- test $r incurred $NEW_FAILS failures, cumulative failures is $NUM_FAILS ---
done
echo Total number of failures: $NUM_FAILS
(echo $NUM_FAILS; git log -1 --format=%at) > git-at
exit $NUM_FAILS
