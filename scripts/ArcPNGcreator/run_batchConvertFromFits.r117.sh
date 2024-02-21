#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
# Variables you may need to change
MCRROOT=$HOME/bin/ArcFinder/MCRlib
exe_dir=$HOME/bin/ArcPNGcreator

if [ "x$1" = "x" ]; then
    echo Usage:
    echo "$0 args"
else
    LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
    MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
    XAPPLRESDIR=${MCRROOT}/X11/app-defaults
    export LD_LIBRARY_PATH;
    export XAPPLRESDIR;
    #echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
    trap "/bin/rm -rf tmp.$$" 0 1 2 3 15
    mkdir -p tmp.$$
    cd tmp.$$ # to mimic being in hash directory until we know the name of the hash directory
    galName=`echo $1 | sed 's,[^0-9], ,g' | newlines | grep '^[0-9][0-9]*$' | head -1`
    hashDir=`echo $galName | sed 's/.*\(..\)$/\1/'`
    [ -d "$hashDir" ] || die "no hash dir <$hashDir>"
    outDir=${10}
    cd ../$outDir
    mkdir -p $hashDir
    N=1
    hostname >> $hashDir/$galName.matlab.log 2>&1
    # Recommended values from Darren for the SDSS run: "0.75 2 0.25"
    # batchConvertFromFits FITS_IN_DIR '\d+' r.fits .png 0.75 2 0.25 MASK_IN_DIR r_starmask.png 255 OUT_DIR
    until "${exe_dir}/batchConvertFromFits.r117b" $1 $2 $3 $4 $5 $6 $7 $8/$hashDir/ $9 $10 $hashDir >> $hashDir/$galName.matlab.log 2>&1
    do
	N=`expr $N - 1`
	if [ $N -eq 0 ]; then break; fi; 
    done
    [ -f "$hashDir/$galName.png" ] || die "no output galaxy"
    chgrp hayesgrp . $hashDir $hashDir/$galName.matlab.log $hashDir/$galName.png
    chmod g+rX . $hashDir $hashDir/$galName.matlab.log $hashDir/$galName.png
    mv $hashDir/$galName.matlab.log $hashDir/$galName`basename $3 .fits`.matlab.log
    mv $hashDir/$galName.png $hashDir/$galName`basename $3 .fits`.png
    #mv $hashDir/batchConvertFromFits.log $hashDir/$galName`basename $3 .fits`.log
fi
