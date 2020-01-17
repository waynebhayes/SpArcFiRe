#!/bin/sh
USAGE="$0 [SpArcFiRe Repo Directory]
This script does not change any files, but is run each time SpArcFiRe wishes to determine
the base of the SpArcFiRe repo.
This script sets the SpArcFiRe repo directory---which is usually just the same directory this
file resides. If you give it an argument, it'll use that directory instead. Otherwise it'll
assume that it's own directory is the SpArcFiRe directory, and all it'll do is echo that
directory."

MYDIR=`dirname $0`
MYDIR=`cd $MYDIR; /bin/pwd`
case $# in
0) echo $MYDIR ;;
*) echo "$1";;
esac
