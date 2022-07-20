#!/bin/bash

#1) Compile delete-commas-inside-quotes script for csvtotsv
##1.A) Get SpArcFiRe dir path:
CURRENTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PARENTDIR="$(dirname "$CURRENTDIR")"
SPARCFIREDIR="$(dirname "$PARENTDIR")"

##1.B) Compile delete-commas-inside-quotes script for csvtotsv
gcc -o "$SPARCFIREDIR"/scripts/delete-commas-inside-quotes "$SPARCFIREDIR"/scripts/delete-commas-inside-quotes.c
