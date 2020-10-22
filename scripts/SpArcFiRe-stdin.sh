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
`basename $0` [options]
OPTIONS:
`grep ').*="$2";' $0 | egrep -v 'grep' | sed -e 's/).*/ VALUE/'`
Run with no options to see syntax of VALUE's (enter a blank line to exit)
Probably the most interesting options are:
    -unsharpMaskAmt: defaults to 10; increase it if arms are not being found
    -resizeDims: size of all output images; defaults to '[256 256]'.  Don't
		 increase it too much, since runtime scales as O(n^4), where
		 n is the side length of the output image. [default n=256].
		 (Runtime is O(N^2) where N is the total number of pixels in
		 the image.)

Each input line must consist of 5 columns:

    input_dir imageBaseName imageSuffix starMaskSuffix output_dir

starMaskSuffix should be 'NONE' to specify no starmask file.
Suffixes must include the '.', eg. the suffix for 'image.png' is '.png'.

Example input line if the input file is a PNG called 'ic2537_r.png' in the
current directory, and the starmask image is called 'ic2537_r_starmask.png'

    . ic2537 _r.png _r_starmask.png output_directory

A blank input line (or EOF on a pipe, or 'exit' or 'quit') causes the server to exit."

while [ "$1" != "" ]; do
    if [ "$#" -lt 2 ]; then
	echo "$USAGE"
	die "ERROR: need 2 arguments to specify an option." >&2
    fi
    case "$1" in
           -generateOrientationFieldPdf)              GENERATEORIENTATIONFIELDPDF="$2"; export GENERATEORIENTATIONFIELDPDF; shift 2;;
                  -groupOutputByProcess)                     GROUPOUTPUTBYPROCESS="$2"; export GROUPOUTPUTBYPROCESS; shift 2;;
               -groupOutputByInputImage)                  GROUPOUTPUTBYINPUTIMAGE="$2"; export GROUPOUTPUTBYINPUTIMAGE; shift 2;;
            -writeSettingsForEveryImage)               WRITESETTINGSFOREVERYIMAGE="$2"; export WRITESETTINGSFOREVERYIMAGE; shift 2;;
                        -writeBulgeMask)                           WRITEBULGEMASK="$2"; export WRITEBULGEMASK; shift 2;;
                                -useMex)                                   USEMEX="$2"; export USEMEX; shift 2;;
             -numOrientationFieldLevels)                NUMORIENTATIONFIELDLEVELS="$2"; export NUMORIENTATIONFIELDLEVELS; shift 2;;
                              -mirrorLR)                                 MIRRORLR="$2"; export MIRRORLR; shift 2;;
                        -useSubpixelCtr)                           USESUBPIXELCTR="$2"; export USESUBPIXELCTR; shift 2;;
                        -clusSizeCutoff)                           CLUSSIZECUTOFF="$2"; export CLUSSIZECUTOFF; shift 2;;
                       -minMinorAxisLen)                          MINMINORAXISLEN="$2"; export MINMINORAXISLEN; shift 2;;
                   -fitUsingNonUsmIVals)                      FITUSINGNONUSMIVALS="$2"; export FITUSINGNONUSMIVALS; shift 2;;
                     -allowArcBeyond2pi)                        ALLOWARCBEYOND2PI="$2"; export ALLOWARCBEYOND2PI; shift 2;;
               -useImageStandardization)                  USEIMAGESTANDARDIZATION="$2"; export USEIMAGESTANDARDIZATION; shift 2;;
                            -resizeDims)                               RESIZEDIMS="$2"; export RESIZEDIMS; shift 2;;
                            -medFiltRad)                               MEDFILTRAD="$2"; export MEDFILTRAD; shift 2;;
                     -useGalfitResidual)                        USEGALFITRESIDUAL="$2"; export USEGALFITRESIDUAL; shift 2;;
                    -generateFitQuality)                       GENERATEFITQUALITY="$2"; export GENERATEFITQUALITY; shift 2;;
                        -unsharpMaskAmt)                           UNSHARPMASKAMT="$2"; export UNSHARPMASKAMT; shift 2;;
                      -unsharpMaskSigma)                         UNSHARPMASKSIGMA="$2"; export UNSHARPMASKSIGMA; shift 2;;
                   -useDeProjectStretch)                      USEDEPROJECTSTRETCH="$2"; export USEDEPROJECTSTRETCH; shift 2;;
                           -fixToCenter)                              FIXTOCENTER="$2"; export FIXTOCENTER; shift 2;;
                 -useTwoStageCtrFinding)                    USETWOSTAGECTRFINDING="$2"; export USETWOSTAGECTRFINDING; shift 2;;
                          -lookForBulge)                             LOOKFORBULGE="$2"; export LOOKFORBULGE; shift 2;;
              -ctrDriftThresForStarMask)                 CTRDRIFTTHRESFORSTARMASK="$2"; export CTRDRIFTTHRESFORSTARMASK; shift 2;;
                         -barCandCutoff)                            BARCANDCUTOFF="$2"; export BARCANDCUTOFF; shift 2;;
                          -barDetCutoff)                             BARDETCUTOFF="$2"; export BARDETCUTOFF; shift 2;;
                                -nhSize)                                   NHSIZE="$2"; export NHSIZE; shift 2;;
                             -stopThres)                                STOPTHRES="$2"; export STOPTHRES; shift 2;;
                     -mergeChkMinClusSz)                        MERGECHKMINCLUSSZ="$2"; export MERGECHKMINCLUSSZ; shift 2;;
                   -balClusWtsInMerging)                      BALCLUSWTSINMERGING="$2"; export BALCLUSWTSINMERGING; shift 2;;
                         -errRatioThres)                            ERRRATIOTHRES="$2"; export ERRRATIOTHRES; shift 2;;
               -failWhenNoStarmaskFound)                  FAILWHENNOSTARMASKFOUND="$2"; export FAILWHENNOSTARMASKFOUND; shift 2;;
         -deleteClusterContainingCenter)            DELETECLUSTERCONTAININGCENTER="$2"; export DELETECLUSTERCONTAININGCENTER; shift 2;;
-ignoreJaggedBoundaryPixelsDuringMerges)   IGNOREJAGGEDBOUNDARYPIXELSDURINGMERGES="$2"; export IGNOREJAGGEDBOUNDARYPIXELSDURINGMERGES; shift 2;;
                       -recomputeCenter)                          RECOMPUTECENTER="$2"; export RECOMPUTECENTER; shift 2;;
                 -imageGuidingThreshold)                    IMAGEGUIDINGTHRESHOLD="$2"; export IMAGEGUIDINGTHRESHOLD; shift 2;;
                                      *) echo "$USAGE" >&2; die "unknown option '$1'";
    esac
done
echo "Enter 'help' for syntax; blank line or 'exit' or 'quit' exits. Wait for MATLAB to start..."
# only argument is a settings file, NONE for none.
set -x
tee SpArcFiRe-stdin.stdin.txt | exec $SPARCFIRE_HOME/scripts/ArcServer/run_findClusterArcsServer.sh $MCRROOT NONE
