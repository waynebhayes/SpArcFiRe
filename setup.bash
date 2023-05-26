#!/bin/bash
SETUP_USAGE="USAGE:

    source setup.bash [SpArcFiRe Repo Directory]

This script does not change any files, but sets the environment variable SPARCFIRE_HOME to
the base of the SpArcFiRe repo---which is usually just the same directory this
file resides. If you give it an argument, it'll use that directory instead. Otherwise it'll
assume that it's own directory is the SpArcFiRe directory, and all it'll do is echo that
directory.
Note: the 'source' command assumes you're using bash as your shell; if not, you'll need
to figure out how to set the environment variable SPARCFIRE_HOME yourself.
"
fail() { echo "$@" >&2; return 1
}

[ "$SPARCFIRE_HOME" != "" ] && return

if [ "`echo $0 | sed 's,.*/,,'`" = setup.sh ]; then
    fail "$SETUP_USAGE
    SETUP ERROR: You've run this script; source it instead by typing:
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

PIP_NEED='numpy|Pillow|scipy|astropy|pandas'

echo -n "Checking you have Python 2.7 or 3 installed"
if python --version | fgrep 2.7 2>/dev/null; then
    PYTHON=python
    PYTHON_SUFFIX=".py2"
    PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

elif python --version | fgrep 3 2>/dev/null; then
    PYTHON=python
    PYTHON_SUFFIX=".py3"
    PIP_HAVE=`(pip3 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

elif python2.7 --version 2>/dev/null; then
    PYTHON=python2.7
    PYTHON_SUFFIX=".py2"
    PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

elif python3 --version 2>/dev/null; then
    PYTHON=python3
    PYTHON_SUFFIX=".py3"
    PIP_HAVE=`(pip3 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

else
    fail "$SETUP_USAGE${NL} SETUP ERROR: You need to install Python 2.7 or 3, and have the executable called python2.7 or python3"
fi
echo "Your python is called $PYTHON:"
$PYTHON --version
export SPARCFIRE_PYTHON="$PYTHON"
export PYTHON_SUFFIX="$PYTHON_SUFFIX"

#PIP_NEED='numpy|Pillow|scipy|astropy|pandas'
#PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`
if [ `echo "$PIP_HAVE" | wc -l` -eq 4 ]; then
    echo "SUCCESS! SPARCFIRE_HOME set to $SPARCFIRE_HOME. Now adding $SPARCFIRE_HOME/scripts to PATH."
    export PATH="$SPARCFIRE_HOME/scripts:$PATH"
else
    (echo "We need all of the following Python packages: `echo "$PIP_NEED" | sed 's/|/ /g'`"
    echo But you only have the following:
    if echo $PIP_HAVE | grep . >/dev/null; then echo "$PIP_HAVE"; else echo "    (none)"; fi
    echo "To get the missing packages, please execute following commands, and note you may need to specify '--user':"
    echo ""
    echo "$PIP_NEED" | tr '|' "$NL" | fgrep -v -f <(echo "$PIP_HAVE" | tr ' ' "$NL") |
	awk '{printf "\t'$PYTHON' -m pip install [--user] %s\n",$0}'
    echo ""
    echo SpArcFiRe repo is in "$SPARCFIRE_HOME"
    echo "If you ran this script without the word 'source' before it, you messed up. Try again.") >&2
    fail "SETUP ERROR: missing pip packages"
fi
