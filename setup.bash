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

PIP_NEED='numpy|Pillow|scipy|astropy'

echo "Checking you have Python 2.7 or 3 installed."
if python --version 2>&1 | grep -q '2.7.5'; then
    PYTHON=python
    PYTHON_SUFFIX=".py2"
    PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

elif python --version 2>&1 | fgrep -q "3"; then
    PYTHON=python
    PYTHON_SUFFIX=".py3"
    PIP_HAVE=`(pip3 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

elif python2.7 --version 2>&1 | grep -q '2.7.5'; then
    PYTHON=python2.7
    PYTHON_SUFFIX=".py2"
    PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

elif python3.7 --version 2>&1; then
    PYTHON=python3.7
    PYTHON_SUFFIX=".py3"
    PIP_HAVE=`(pip3 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`

else
    fail "$SETUP_USAGE${NL} SETUP ERROR: You need to install Python 2.7 or 3, and have the executable called python2.7 or python3"
fi
echo "Your python is called by $PYTHON and is version:"
$PYTHON --version
export SPARCFIRE_PYTHON="${PYTHON}"
export PYTHON_SUFFIX="${PYTHON_SUFFIX}"

# Regex to check if python >=3.7 (includes 3.10+)
# NOTE: Would succeed for 3.1
#if python3 --version 2>&1 | grep -q '3\.[1,7-9]\?[0-9]\.[0-9]\?[0-9]'; then
echo ""
echo "Now checking for python3(.7 or greater) for GalfitModule."
PIP_NEED3='numpy|Pillow|scipy|astropy|pandas|IPython'
if python3 -c 'import sys; assert sys.version_info >= (3,7), "Python3.7 or newer needed."'; then
    export PYTHON3=python3
    PIP_HAVE3=`(pip3 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED3"`
    echo "Your python is called by $PYTHON3 and is version:"
    $PYTHON3 --version
else
    echo "python3 currently calls $(python3 --version)"
    echo "Python3.7 or newer must be called by 'python3' to use GalfitModule."
fi
echo ""

#PIP_NEED='numpy|Pillow|scipy|astropy|pandas'
#PIP_HAVE=`(pip2 list; $PYTHON -m pip list) 2>/dev/null | awk '{print $1}' | sort -u | egrep "$PIP_NEED"`
if [ `echo "$PIP_HAVE" | wc -l` -eq 4 ]; then
    echo "SUCCESS! SPARCFIRE_HOME set to $SPARCFIRE_HOME. Now adding $SPARCFIRE_HOME/scripts to PATH."
    echo ""
    export PATH="$SPARCFIRE_HOME/scripts:$PATH"
else
    (echo "We need all of the following Python packages for $PYTHON: `echo "$PIP_NEED" | sed 's/|/ /g'`"
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


if [ `echo "$PIP_HAVE3" | wc -l` -eq 5 ]; then
    true
elif [ -z "${PYTHON3}" ]; then
    true
else
    (echo "We need all of the following Python packages for $PYTHON3: `echo "$PIP_NEED3" | sed 's/|/ /g'`"
    echo But you only have the following:
    if echo $PIP_HAVE3 | grep . >/dev/null; then echo "$PIP_HAVE3"; else echo "    (none)"; fi
    echo "To get the missing packages, please execute following commands, and note you may need to specify '--user':"
    echo ""
    echo "$PIP_NEED3" | tr '|' "$NL" | fgrep -v -f <(echo "$PIP_HAVE3" | tr ' ' "$NL") |
	awk '{printf "\t'$PYTHON' -m pip install [--user] %s\n",$0}'
    ) >&2
fi
