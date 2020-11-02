#!/bin/sh
USAGE="USAGE:

    source ./setup.sh [SpArcFiRe Repo Directory]

This script does not change any files, but sets the environment variable SPARCFIRE_HOME to
the base of the SpArcFiRe repo---which is usually just the same directory this
file resides. If you give it an argument, it'll use that directory instead. Otherwise it'll
assume that it's own directory is the SpArcFiRe directory, and all it'll do is echo that
directory.
Note: the 'source' command assumes you're using bash as your shell; if not, you'll need
to figure out how to set the environment variable SPARCFIRE_HOME yourself.
"
die() { (echo "$USAGE"; echo "FATAL ERROR: $@") >&2; exit 1
}

if [ `basename $0` == setup.sh ]; then
    die "You've run this script; source it instead by typing:
    source $0"
fi

NL='
' # newline as a variable, for ease later

# Source of following line: https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
case $# in
1) export SPARCFIRE_HOME="$1";;
*) export SPARCFIRE_HOME="$MYDIR" ;;
esac
if [ "$SPARCFIRE_HOME" = . ]; then SPARCFIRE_HOME=`/bin/pwd`; fi

echo -n "Checking you have Python 2.7 installed"
if python --version | fgrep 2.7 2>/dev/null; then
    PYTHON=python
elif python2.7 --version 2>/dev/null; then
    PYTHON=python2.7
else
    die "You need to install Python 2.7, and have the executable called python2.7"
fi
echo "Your python2.7 is called $PYTHON:"
$PYTHON --version

PIP_NEED='numpy|Pillow|scipy|astropy'
PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`
if [ `echo "$PIP_HAVE" | wc -l` -eq 4 ]; then
    echo "SUCCESS! SPARCFIRE_HOME set to $SPARCFIRE_HOME. Now adding $SPARCFIRE_HOME/scripts to PATH."
    export PATH="$SPARCFIRE_HOME/scripts:$PATH"
else
    echo "$USAGE"
    echo "We need all of the following Python packages: `echo "$PIP_NEED" | sed 's/|/ /g'`"
    echo But you only have the following:
    echo "$PIP_HAVE"
    echo "To get the missing packages, please execute following commands, and note you may need to specify '--user':"
    echo ""
    echo "$PIP_NEED" | tr '|' "$NL" | fgrep -v -f <(echo "$PIP_HAVE" | tr ' ' "$NL") |
	awk '{printf "\t'$PYTHON' -m pip install [--user] %s\n",$0}'
    echo ""
    echo SpArcFiRe repo is in "$SPARCFIRE_HOME"
    echo "If you ran this script without the word 'source' before it, you messed up. Try again."
    die "missing pip packages"
fi
