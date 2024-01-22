#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#

USAGE=

die() { echo "$@" >&2; exit 1
}
# Variables you may need to change
MCRROOT=$HOME/bin/ArcFinder/MCRlib/v717
exe_dir=$HOME/bin/ArcPNGcreator

LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR;
#echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
# batchConvertFromFits FITS_IN_DIR '\d+' r.fits .png UPPER_QUANTILE NUM_ASINHs LOWER_QUANTILE MASK_IN_DIR r_starmask.png 255 OUT_DIR
USAGE="$0 [-p blackQlvl qLvl nReps ] inFitsFile {starMaskFile|NONE} outPNGfile"
PARAMS="0.75 2 0.25" # on the command line this would be "0.25 0.75 2"
case "$1" in
-p) PARAMS="$3 $4 $2"; shift 4;;
esac
FITSfile="$1"
starmask="$2"
pngFile="$3"
[ -f "$FITSfile" ] || die "second arg must be FITs input file, and it must exist"
#[ ! -f "$pngFile" ] || die "third arg must be PNG output file"
#./convertFromFits.r123 "$FITSfile" "$starmask" 255 $PARAMS "$pngFile"
#exit

# WARNING: we used to assume galaxy names matched "[a-zA-Z0-9_-]+", but it seemed to restrictive, so now we're allowing
# galaxy names to contain any characters at all except dots, which are only allowed before the extension (.fits, .png, etc.)
NAME_REGEXP='[^.]+'     # was "[a-zA-Z0-9_-]+"

if echo $starmask | grep -i none >/dev/null; then
    "${exe_dir}/batchConvertFromFits.r120" "$fitsDir" "$NAME_REGEXP" .fits .png $PARAMS NONE NONE NONE $pngDir
else
    "${exe_dir}/batchConvertFromFits.r120" "$fitsDir" "$NAME_REGEXP" .fits .png $PARAMS "$pngDir" "$starmask" 255 $pngDir
fi
