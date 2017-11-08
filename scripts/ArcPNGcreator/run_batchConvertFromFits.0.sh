#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`pwd`
#exe_dir="/home/drdavis/ArcFinder"
#echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  #echo Setting up environment variables
  MCRROOT="$1"
  #echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
	MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export LD_LIBRARY_PATH;
  export XAPPLRESDIR;
  #echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  mkdir -p tmp
  cd tmp # to mimic being in hash directory until we know the name of the hash directory
  shift # get rid of MCRROOT
  galName=`echo $1 | sed 's,[^0-9], ,g' | newlines | grep '^[0-9][0-9]*$' | head -1`
  hashDir=`echo $galName | sed 's/.*\(..\)$/\1/'`
  outDir=${10}
  cd ../$outDir
  mkdir -p $hashDir
  N=1
  hostname >> $hashDir/$galName.matlab.log 2>&1
  # batchConvertFromFITS FITS_IN_DIR '\d+' r.fits .png 0.9 2 MASK_IN_DIR r_starmask.png 255 OUT_DIR
  until "${exe_dir}/batchConvertFromFits" $1 $2 $3 $4 $5 $6 $7/$hashDir/ $8 $9 $hashDir >> $hashDir/$galName.matlab.log 2>&1
  do
      N=`expr $N - 1`
      if [ $N -eq 0 ]; then break; fi; 
  done
  chgrp hayesgrp . $hashDir $hashDir/$galName.matlab.log $hashDir/$galName.png
  chmod g+rX . $hashDir $hashDir/$galName.matlab.log $hashDir/$galName.png
  mv $hashDir/$galName.matlab.log $hashDir/$galName`basename $3 .fits`.matlab.log
  mv $hashDir/$galName.png $hashDir/$galName`basename $3 .fits`.png
  #mv $hashDir/batchConvertFromFits.log $hashDir/$galName`basename $3 .fits`.log
fi
exit
