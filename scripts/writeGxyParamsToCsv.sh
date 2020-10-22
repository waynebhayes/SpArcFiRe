#!/bin/bash
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
#HOME=/home/sparcfire   # Make this the base of your SpArcFiRe bin/ distribution, or leave commented to use default HOME.
unset DISPLAY # ensure no X windows pop up
die() { echo "$@" >&2; exit 1
}
exe_name=$0
exe_dir=`dirname "$0"`

# SETTING THE MCRROOT DEFINES WHICH VERSION OF THE MATLAB RUNTIME LIBRARY TO USE.
# USE THE FOLLOWING MCRROOT IF YOU DOWNLOADED THE DISTRIBUTION,
# WHICH USES THE 2013 MATLAB RUNTIME LIBRARY DISTRIBUTED WITH SPARCFIRE
#MCRROOT=$HOME/global-data/MCRlib/v717

#The MCRROOTs below are for running on the UCI Computer Science servers
MCRROOT=/pkg/matlab/R2017a
#MCRROOT=/pkg/matlab/current
#MCRROOT=/pkg/matlab/R2014a
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

USAGE="USAGE:
`basename $0` <matout dirname> <cssv filename>"

exec $HOME/bin/ArcServer/writeGxyParamsToCsv.r125 "$@"
