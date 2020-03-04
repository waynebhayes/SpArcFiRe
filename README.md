# SpArcFiRe
SPiral ARC FInder and REporter - find spiral arcs in galaxies.... or maybe hurricanes too?
The "documentation" directory contains the original paper and some other (possibly irrelevant) info.
The documentaton below is up-to-date as of March 2020.
If you have any questions, please email whayes@uci.edu, or jeenglis@uci.edu


## How to Use SpArcFiRe ##

Usage

>~/SpArcFiRe/scripts/SpArcFiRe [-convert-FITS [-p UPPER LOWER N]] [-compute-starmask &lt;true | false>] [-ignore-starmask] [-elps_dir DIR] [-web] {inDir} {tmpDir} {outDir} [other optional options (see below)]

To use recommended arguments

> ~/bin/wschallo/SpArcFiRe-run.sh

This is the equivalent of

>~/ SpArcFiRe/scripts/SpArcFiRe -convert-FITS ~/SDSS/G.in ~/SDSS/G.tmp ~/SDSS/G.out -generateFitQuality 0 -writeBulgeMask 1


## Optional Arguments ##

**-convert-FITS [-p UPPER LOWER N]**

When converting from FITS files to PNG there is the option to cut out certain quantiles of brightnesses to increase the contrast of the image and remove background noise.  UPPER and LOWER are the values of the upper and lower quartile cutoffs, they should be a number in the range [0,1].  N is the number of times to apply the arcsine transform to the image.  If this flag is not set then it is assumed that the images have already been converted and are in {inDir} (see below), once they have been converted once there is no need to convert them again unless the -p options are being changed.

**-compute-starmask &lt;true | false>**

In order to avoid detecting noise or stars as part of the galaxy an alpha mask is created of what is best estimated to be the only galaxy.  This should always be set to true unless it has been computed on a previous run or if -ignore-starmask is set (see below).

**-ignore-starmask**

Enabling this will tell SpArcFiRe not to attempt to block out the noise or stars in the background of images.  This could lead to some false detections of arcs so it is best to use the star masks unless the images have low noise.

**-elps_dir DIR**

In order to view a galaxy straight on SpArcFiRe has to warp the original image of the galaxy into one where it appears circular.  This is done by fitting an ellipse around the galaxy and warping it into a circle, which is done automatically by SpArcFiRe unless this option is set.  If this option is set then a file following the format GALAXY_NAME_epls.txt is expected for each galaxy in DIR.


## Required Arguments ##


**inDir**

SpArcFire can handle image types of FITS, JPG, and PNG, it expects all galaxy images to be one of these types and in this directory.  A directory located at ~/SDSS/G.in has been created for this purpose.

**tmpDir**

All image types are converted into PNG files for processing, this directory stores these files along with the star masks.  A directory located at ~/SDSS/G.tmp has been created for this purpose.

**outDir**

All SpArcFire output will be saved here.  By default, each galaxy will have its own directory where images showing the steps SpArcFire takes to detect arm segments is visualized.  The data of each arm segment is located at outDir/galaxy_arc.csv and the data of each galaxy is located at outDir/galaxy.csv.  A directory located at ~/SDSS/G.out has been created for this purpose.


## Other Optional Arguments ##

**-generateOrientationFieldPdf &lt;0 | 1> {Default is 1}**

Controls whether an image showing the orientation field is saved to the output directory, the file is called GALAXY_NAME-O_orientation-field.png.

**-groupOutputByProcess &lt;0 | 1> {Default is 0}**

If set to 1 the galaxy output directories will be put into a directory called “proc” in outDir so that they can all be easily moved somewhere else when another batch of galaxies is run, if set to 0 all directories will be in outDir.

**-groupOutputByInputImage &lt;0 | 1>  {Default is 1}**

If set to 1 each galaxy will get its own output directory in outDir, if set to 0 then all of these files will be put directly into the outDir.

**-writeSettingsForEveryImage &lt;0 | 1> {Default is 1}**

If set to 1 the settings described here will be outputted into a file called GALAXY_NAME-S_settings.txt, if set to 0 this file won’t be created.

**-writeBulgeMask &lt;0 | 1> {Default is 0}**

If set to 1 then an image will be saved in the galaxy output directory representing where the estimated the core of the galaxy was, the image is called GALAXY_NAME-D1_bulgeMask.png.

**-numOrientationFieldLevels VALUE {Default is 3}**

Controls the amount of different orientation angles that are checked for in the galaxy.  In general, higher this value the more arcs are that are found, it is recommended that this value is left at 3 since anything higher has a good chance of detecting noise as arm segments.

**-mirrorLR &lt;0 | 1> {Default is 0}**

If set to 1 then the galaxy image will be flipped along the horizontal axis before being processed.

**-useSubpixelCtr &lt;0 | 1> {Default is 1}**

Setting this to 0 will tell SpArcFiRe to do all calculations along pixel blocks, this will change the output data slightly but shouldn’t change any conclusions drawn from it.

**-clusSizeCutoff VALUE {Default is 150}**

The value given is the minimum length in pixels for an arm segment to be considered an actual segment.

**-minMinorAxisLen VALUE {Default is 0}**

The minimum required length in pixels of the minor axis of a galaxy, if a galaxy does not meet this requirement then no arcs will be found on the galaxy.

**-allowArcBeyond2pi &lt;0 | 1> {Default is 1}**

Some arm segments might circle around a galaxy completely, by default SpArcFiRe allows these arms to exist but if set to 0 then it will cut them off arm at 2pi (or one complete revolution).

**-useImageStandardization &lt;0 | 1> {Default is 1}**

If set to 0 the intermediate image steps are not saved as files, this does not affect arc output.

**-resizeDims WIDTH HEIGHT {Default is 250 250}**

Controls the size of the output images, the default is 256x256.  These should be small images because the runtime will increase with O(N^2).

**-medFiltRad VALUE {Default is 1}**

Controls how much the image is blurred radially before being processed.  Some amount is recommended as an image that is too fine grain will produce many small arm segments where there should only be one.

**-useGalfitResidual &lt;0 | 1> {Default is 0}**

If set to 1 then it attempts to use the GalFit generated galaxy for processing.

**-unsharpMaskAmt VALUE {Default is 6}**

Controls how much the contrast of the image is increased prior to being processed by, some amount is recommended as it helps separate arcs from noise or dust.

**-unsharpMaskSigma VALUE {Default is 25}**

Controls the spread of how the image to be processed is unsharpened, the higher the value the higher the spread and vice versa.

**-useDeProjectStrech &lt;0 | 1> {Default is 1}**

If set to 1 then when cropping the image for processing the image will also be rotated and stretched to make the galaxy fill up as much of the image as possible.

**-fixToCenter &lt;0 | 1> {Default is 0}**

If set to 1 then all arm segments fit around the center of the image.  This should not have much of an effect on the output data since the galaxies are cropped so the center of the galaxy is in the center of the image.

**-useTwoStageCtrFinding &lt;0 | 1> {Default is 1}**

SpARcFiRE runs over the the input image in two stages to find clusters, if set to 0 then it will only go over the image once.  This could lead to some arcs not being found so it is recommended that this is left enabled.

**-lookForBulge &lt;0 | 1> {Default is 1}**

If set to 0 then no attempt is made to find the bulge at the center of a galaxy, this can effect arc and bar detection.

**-stopThres VALUE [0-1] {Default is 0.15}**

The value given is the minimum brightness needed for an arm segment to be considered an arc.

**-mergeChkMinClusSz VALUE {Default is 25}**

Controls the minimum value in pixels that a detected arm segment must be for it to be considered when trying to detect arm segments to merge.

**-balClusWtsInMerging &lt;0 | 1> {Default is 1}**

If set to 0 then clusters that are being analyzed for merging will weighted based on their size, this leads to smaller clusters being overpowered by larger ones, potentially missing entire arm segments so it is recommended that this is left enabled, although it can be used to simplify the galaxy output.

**-errRatioThres  VALUE {Default is 2.5}**

Controls the amount of error (or difference) that is allowed between clusters when merging, the higher the value the more merging that will occur and vise versa.

**-failWhenNoStarmaskFound &lt;0 | 1> {Default is 0}**

If set to 1 and no star mask file is found then then SpArcFiRe will fail instead of the default behavior which is just to create the star mask.

**-deleteClusterContainingCenter &lt;0 | 1> {Default is 1}**

If a cluster is found to contain the center of the galaxy then it will not be counted as an arm segment.

**-paApproxLevel &lt;0 | 1> {Default is 0}**

If set to 1 then an approximate pitch angle for each arc is used instead of the slower -deterministic method used by default.


## Features Under Maintenance ##

-web

-sleepSecondsAfterImageWrite VALUE {Default is 0}

-ignoreJaggedBoundaryPixelsDuringMerges &lt;0 | 1> {Default is 1}

-useMex &lt;0 | 1> {Default is 1}

-fitUsingNonUsmIVals &lt;0 | 1> {Default is 1}

-nhSize VALUE {Default is 1}

-barCandCutoff VALUE {Default is 7}

-barDetCutoff VALUE {Default is 2}

-ctrDriftThresForStarMask VALUE {Default is 2.5}

-generateFitQuality &lt;0 | 1> {Default is 1}
