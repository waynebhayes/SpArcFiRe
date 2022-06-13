#!/bin/bash

#1)Create symbolic link to sparcfire:
##1.A) Get SpArcFiRe dir path:
CURRENTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PARENTDIR="$(dirname "$CURRENTDIR")"
SPARCFIREDIR="$(dirname "$PARENTDIR")"

##1.B) Create symbolic link:
ln -s "$SPARCFIREDIR"/scripts/ "$HOME"/bin
