#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#

die() { echo "$@" >&2; exit 1
}
# Variables you may need to change
#case `uname -a | awk '{print $3}'` in
#3.10.0*) export MCRROOT=/pkg/matlab/current; R=121;
#    ;;
#2.6.32*) export MCRROOT=$HOME/bin/ArcFinder/MCRlib/v717; R=120;
#    #unless it's odin
#    case `hostname` in
#    odin*) export MCRROOT=/pkg/matlab/current; R=121;
#	;;
#    esac
#    ;;
#esac
export MCRROOT=/pkg/matlab/R2017a;
#R="test";
R="123";
exe_dir=$SPARCFIRE_HOME/scripts/ArcPNGcreator

[ "$HOME" = "" ] && die "HOME not set!"

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
PARAMS="0.75 2 0.25" # on the command line this would be "0.25 0.75 2"
case "$1" in
-p) PARAMS="$3 $4 $2"; shift 4;;
esac
starmask="$1"
fitsDir="$2"
pngDir="$3"
[ -d "$fitsDir" ] || die "second arg must be FITs input directory, and it must exist"
[ -d "$pngDir" ] || die "third arg must be PNG output directory, and it must exist"
if echo $starmask | grep -i none >/dev/null; then
    "${exe_dir}/batchConvertFromFits.r$R" "$fitsDir" "[a-zA-Z0-9_]+" .fits .png $PARAMS NONE NONE NONE $pngDir
else
    "${exe_dir}/batchConvertFromFits.r$R" "$fitsDir" "[a-zA-Z0-9_]+" .fits .png $PARAMS "$pngDir" "$starmask" 255 $pngDir
fi
