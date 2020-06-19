function stgs = getDefaultSettings()

% whether to put output from each process in its own folder
stgs.groupOutputByProcess = false;
% whether to put output from each input image in its own folder
stgs.groupOutputByInputImage = true;
% whether to write settings information for each image
stgs.writeSettingsForEveryImage = false;
% whether to provide an orientation field PDF; will require several extra
% seconds to plot.
stgs.generateOrientationFieldPdf = false;

% How long to wait after writing an image, for contexts where each process
% should not use more than one CPU core (though it is not clear that Matlab
% actually uses non-blocking I/O).
stgs.sleepSecondsAfterImageWrite = 0.0;

% whether to use the C++ implementation where available (currently for the
% clustering only)
stgs.useMex = true;

% Number of image scales to use when computing the orientation field.
% The dimensions of the preprocessed image (see resizeDims) must be 
% divisible by 2^(numOrientationFieldLevels-1).
stgs.numOrientationFieldLevels = 3;

% Whether to flip the input image left-to-right (e.g., for chirality bias
% testing).
stgs.mirrorLR = false;

% Whether to skip fitting an image if a starmask suffix is provided but no
% starmask is found. If false under these conditions, arc-finding will 
% proceed as if no star mask suffix was provided, but different 
% images could be treated differently within the same image set (i.e., 
% where the code would use a starmask, the actual use of the starmask
% depends on whether the starmask file is present). 
stgs.failWhenNoStarmaskFound = true;

% DEPRECATED - should always be true
% Whether to track the center more precisely during preprocessing.
stgs.useSubpixelCtr = true;

% Whether to remove a cluster, after the main clustering and bar-cluster
% extraction but before secondary merging, if the cluster contains the
% center pixel.
stgs.deleteClusterContainingCenter = true;

% minimum size needed for a cluster to be included in the output 
stgs.clusSizeCutoff = 150;
% minimum minor axis length (in pixels) of the disk for trying to produce
% output for the galaxy/image
stgs.minMinorAxisLen = 0;

% whether to use the image intensities as the fit weights (instead of the
% intensities after the unsharp mask)
stgs.fitUsingNonUsmIVals = true;
% whether to try to fit a spiral with theta-span >2*pi, if there appears 
% to be more than one theta-revolution in the arc. This can capture a wider
% range of arcs, but more things can go wrong. Otherwise, for all arcs, the 
% attempted fit will be for a log spiral of theta-span <= 2*pi
stgs.allowArcBeyond2pi = true;

stgs.zeroThetaStart = true;

% In order for a cluster's theta-endpoints to be considered well defined,
% there must be a path of empty pixels that goes from the image center to
% a point on the image boundary, leaving at least the given number of 
% pixels between the pixel path and all cluster pixels. 
stgs.maxHalfGapFillForUndefinedBounds = 3;

% Removes jagged cluster edges to make the merge score slightly more
% accurate, but requires more least-squares fits and usually does not have
% a large effect on the final output.
stgs.ignoreJaggedBoundaryPixelsDuringMerges = true;

% recalculates the center of the galaxy and saves it to the elps file,
% regardless of other settings. If using another elps file, SpArcFiRe
% will still perform the ellipse fit but will disregard everything 
% but the center.
stgs.recomputeCenter = false;

% DEPRECATED - should always be 0
% How much to use an approximation (vs least-squares fitting) for the 
% pitch angle
% 0: no approximation; use the least-squares fit only
% 1: use the approximation as the initial state for least-squares fitting
% 2: use the approximation as the final value. Fit error is still
% least-squares error, so cluster error could be overestimated in
% non-uniform ways.
stgs.paApproxLevel = 0;

% ---------------------
% IMAGE STANDARDIZATION
% ---------------------

% Whether to actually perform the image standardization. If false, there
% will be no recentering, stretching, or resizing, and ellipse-fit 
% parameters will not be available, but an unsharp mask and/or median 
% filtering will still be applied (if the settings specify to do so).
stgs.useImageStandardization = true;
% The dimensions to resize the image to after preprocessing.
stgs.resizeDims = [256 256];
% Radius of the median filter to use (before preprocessing, including
% resizing), so 0 means no median filter, 1 means a 3x3 filter and so on
stgs.medFiltRad = 1;
% strength of unsharp mask during preprocessing; 0 means no USM
stgs.unsharpMaskAmt = 6;
% Gaussian scale during unsharp masking.
stgs.unsharpMaskSigma = 25;
% stgs.temp_useWideUsm = false;
% whether to stretch the image so that the fitted Gaussian ellipse would 
% become circular, in an attempt to remove inclination effects 
% (this also requires image rotation)
stgs.useDeProjectStretch = true;
% whether to assume the galaxy has been exactly centered already
stgs.fixToCenter = false;
% use two rounds to get a more precise center in some cases. No effect if
% fixToCenter = true
stgs.useTwoStageCtrFinding = true;
% Also report ellipse statistics for galaxy bulge. Requires more Gaussian 
% fit iterations, and thus more time (especially for higher-resolution
% input images).
stgs.lookForBulge = true;
% threshold at which the distance between the fitted and image center is
% high enough that interference from a star is a likely culprit.
stgs.ctrDriftThresForStarMask = 1 / 0.4;

% -------------
% BAR DETECTION
% -------------

% Minimum scores that the bar must meet in order to be available as an
% alternative to a log-spiral.
stgs.barCandCutoff = 7;
stgs.barDetCutoff = 2;

% -------------------------
% SIMILARITY MATRIX AND HAC
% -------------------------
% maximum horizontal, vertical, or diagonal distance for a pixel similarity
% value to be nonzero.  Due to time constraints, a value of 1 is 
% recommended if used for HAC.
stgs.nhSize = 1;
% termination condition for clustering; clustering stops when all 
% similarity values are below this threshold
stgs.stopThres = 0.15;

% ---------------
% CLUSTER MERGING
% ---------------

% Minimum size (in pixels) both clusters must reach in order to trigger a 
% fit-based merge check.
stgs.mergeChkMinClusSz = 25; %TODO: make a function of the image size
% Whether to make both clusters' weights sum to the same value when fitting
% arcs during fit-based merge checking. It is recommended to keep this set
% to true.
stgs.balClusWtsInMerging = true;
% Maximum ratio of merged to unmerged fit error. In order for a merge to
% proceed, this criterion must be met for both clusters involved in the
% potential merge.
stgs.errRatioThres = 2.5;

stgs.imageGuidingThreshold = -1;

end
