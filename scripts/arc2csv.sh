#!/bin/bash -x
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
die() { echo "ERROR: $@" >&2; exit 1
}
#[ $# = 1 ] || die "need exactly one argument, filename.fits"
USAGE="$0 {matout directory} {csvName}"
[ $# = 2 ] || die "$USAGE"

exe_name=$0
exe_dir="$SPARCFIRE_HOME/scripts/ArcServer"
#MCRROOT=/pkg/matlab/current
#MCRROOT=$HOME/global-data/MCRlib/v717
#MCRROOT=/pkg/matlab/R2014a
MCRROOT=/pkg/matlab/R2017a
#MCRROOT=/pkg/matlab/7.14_r2012a # What Araceli had
#MCRROOT=/pkg/matlib/7.11_r2010b

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
export MCRROOT;
#echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
N=1
until ${exe_dir}/writeGxyParamsToCsv.r125 "$@"
do
  N=`expr $N - 1`
  if [ $N -eq 0 ]; then break; fi; 
done
