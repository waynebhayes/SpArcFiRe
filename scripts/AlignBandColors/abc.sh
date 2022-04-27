#!/bin/sh

rm -rf "tmp"
mkdir "tmp"
python3 -W "ignore" shift_gal.py $*
