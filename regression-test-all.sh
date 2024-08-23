#!/bin/bash
USAGE="USAGE: $0 [-use-git-at]  [ list of tests to run, defaults to regression-tests/*/*.sh ]"
NL='
'
die() { echo "$USAGE${NL}FATAL ERROR: $@" >&2; exit 1;}
TestGawk(){ echo hello | gawk '{++array[$1][$1]}' || die "Upgrade to a [g]awk that supports multi-dimensional arrays"; }
TestGawk

while [ "X$1" != X ]; do
    case "$1" in
    -use-git-at)
	if [ -f git-at ] && [ `wc -l < git-at` -eq 2 -a `git log -1 --format=%at` -eq `tail -1 git-at` ]; then
	    echo -n "Repo unchanged; returning same status code as "
	    tail -1 git-at | xargs -I{} date -d @{} +%Y-%m-%d-%H:%M:%S
	    exit `head -1 git-at`
	fi
	shift
	;;
    -make) echo "nothing to make; ignoring this option" >&2; shift;;
    *) [ -f "$1" ] || die "unknown option '$1'"
	break ;;
    esac
done

PATH=`pwd`:`pwd`/scripts:$PATH
export PATH

if [ ! -x scripts/delete-commas-inside-quotes ]; then
    compile="gcc -o scripts/delete-commas-inside-quotes scripts/delete-commas-inside-quotes.c" || die "please compile scripts/delete-commas-inside-quotes.c"
    eval $compile || die "please do whatever's needed to compile: '$compile'"
fi

#SpArcFiRe ONLY!!
source setup.bash || die "setup failed"

# This now happens in the SpArcFiRe script
# make all

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
    if eval time $STDBUF "$r"; then # force output and error to be line buffered
	:
    else
	NEW_FAILS=$?
	(( NUM_FAILS+=$NEW_FAILS ))
    fi
    echo --- test $r incurred $NEW_FAILS failures, cumulative failures is $NUM_FAILS ---
done

echo --- running GALFIT unit and reg tests ---
echo
if [ -n "${PYTHON3}" ]; then
    GALFIT_REG_OUTPUT="$($PYTHON3 ${SPARCFIRE_HOME}/GalfitModule/RegTest/RegTest.py)"
    echo "$GALFIT_REG_OUTPUT"
    len=$((${#GALFIT_REG_OUTPUT}))
    GALFIT_FAILS="${GALFIT_REG_OUTPUT:(-2):$len}"
else
    echo "Python3 installation either not found or too old... tallying a failure."
    GALFIT_FAILS=1
fi


(( NUM_FAILS+=GALFIT_FAILS ))
echo
echo --- GALFIT unit and reg tests incurred $GALFIT_FAILS failures, cumulative failures is $NUM_FAILS ---

echo Total number of failures: $NUM_FAILS
(echo $NUM_FAILS; git log -1 --format=%at) > git-at
exit $NUM_FAILS
