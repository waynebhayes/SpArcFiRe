#!/bin/bash

if [[ $# -ne 3 ]]; then echo "USAGE: downloadFields ra dec outdir"; exit; fi

# download and unzip the field images to outdir
wget -O - -o /dev/null "https://dr12.sdss.org/fields/raDec?ra=$1&dec=$2" | awk '/bz2/ {print}' | grep -o '".*"' | sed 's;^";https://dr12.sdss.org;g' | sed 's;";;g' | xargs -n1 wget -P $3

echo "Unzipping field images..."
bzip2 -d $3/*
