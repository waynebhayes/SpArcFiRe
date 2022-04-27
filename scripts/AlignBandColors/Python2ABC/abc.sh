#!/bin/sh

rm -rf "tmp"
mkdir "tmp"
python2 -W "ignore" shift_gal.py $*
