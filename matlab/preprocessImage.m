
function [img, imgNoUsm, gxyParams, fitParams, exactCtrR, exactCtrC, fromOrigImg, masked] = ...
    preprocessImage(img, stgs, starMask, prevFitParams, outputPath, gxyName)
% Transforms the given image into a standardized/enhanced form for
% clustering and arc-fitting.
% INPUTS:
%   img: the image to use for preprocessing
%   stgs: structure containing algorithm settings (see settings.m)
%   starMask: (optional) star mask where the lowest nonzero value indicates
%       an object and increasingly larger values indicate regions
%       progressively more aggressive in excluding stars (at the risk of 
%       excluding part of the galaxy)
%   prevFitParams: if present, uses the normalization (e.g., de-projection)
%       parameters previously calculated from another image (useful for
%       applying the same de-projection to multiple wavebands)
% OUTPUTS:
%   img: the preprocessed image
%   imgNoUsm: the preprocessed image, without applying the unsharp mask
%   gxyParams: structure containing some information about the ellipse fit

if nargin < 3
    starMask = [];
end

if nargin < 4
    prevFitParams = [];
end

if nargin < 5
    outputPath = '';
    gxyName = '';
end

masked = [];

plotFlag = false;

useImageStandardization = stgs.useImageStandardization;
resizeDims = stgs.resizeDims;
useDeProjectStretch = stgs.useDeProjectStretch;
fixToCenter = stgs.fixToCenter;
medFiltRad = stgs.medFiltRad;
ctrDriftThres = stgs.ctrDriftThresForStarMask;
recomputeCenter = stgs.recomputeCenter;

if ~isempty(resizeDims) && (resizeDims(1) ~= resizeDims(2))
    warning('resizeDims should be square (was %s)', mat2str(resizeDims))
end

[nRows, nCols] = size(img);
imgCtrR = nRows/2;
imgCtrC = nCols/2;
if stgs.useSubpixelCtr
    imgCtrR = imgCtrR + 0.5;
    imgCtrC = imgCtrC + 0.5;
else
    fprintf('WARNING: not using subpixel center')
end

if ~useImageStandardization
    if medFiltRad > 0
        img = medfilt2(img, [1 1] * (2*medFiltRad+1));
    end
    imgNoUsm = img;
	img = unsharpMask(img, stgs);
    img(img < 0) = 0; img(img > 1) = 1;
    gxyParams.warnings = {};
    gxyParams.starMaskUsed = 'not_applicable';
    gxyParams.noiseMaskUsed = 'not_applicable';
    gxyParams.numElpsRefits = [];
    gxyParams.diskAxisRatio = [];
    gxyParams.diskMinAxsLen = [];
    gxyParams.diskMajAxsLen = [];
    gxyParams.diskMajAxsAngleRadians = [];
    gxyParams.iptCtrXY = [imgCtrC imgCtrR];
    gxyParams.iptSz = [nRows nCols];
    gxyParams.muDist = [];
    gxyParams.muDistProp = [];
    gxyParams.wtdLik = [];
    gxyParams.likOfCtr = [];
    gxyParams.brtUnifScore = [];
    gxyParams.gaussLogLik = [];
    gxyParams.contourBrtRatio = [];
    gxyParams.badBulgeFitFlag = [];
    gxyParams.bulgeAxisRatio = [];
    gxyParams.bulgeMajAxsLen = [];
    gxyParams.bulgeMajAxsAngle = [];
    gxyParams.bulgeAvgBrt = [];
    gxyParams.bulgeDiskBrtRatio = [];
    fitParams.rotAngle = 0;
    fitParams.muFit = [(size(img, 2) / 2) + 0.5, (size(img, 1) / 2) + 0.5];
    fitParams.cropRad = min(ceil([size(img, 1), size(img, 2)] / 2));
    gxyParams.fitParams = fitParams;
    exactCtrR = imgCtrR;
    exactCtrC = imgCtrC;
    fromOrigImg = ones(size(img));
    masked = [];
    return
end

gxyParams.warnings = {};

if ~isempty(starMask)
    smVals = sort(unique(starMask(:)), 'ascend');
    if smVals(1) ~= 0
        % every pixel is considered part of an object ?!?
        gxyParams.warnings = [gxyParams.warnings
             'preprocessing:noZeroValuesInStarMask'];
    else
        smVals = smVals(2:end);
    end
    if any(~ismember(smVals, 255 * [1 2 3] / 3))
        error('unrecognized star mask encoding')
    end
    gxyRgnVal = 255;
    extRgnVal = 255 / 3;
    otherObjVal = (2 * 255) / 3;
end

% run_id = sprintf('%15.0f', int64(now()* 1e12));
% imwrite(img, [run_id '_deproject_00.png']);

if medFiltRad > 0
    img = medfilt2(img, [1 1] * (2*medFiltRad+1));
end

% imwrite(img, [run_id '_deproject_01.png']);

likCutoff = 10^-9;
if isempty(prevFitParams) | recomputeCenter
%%%perform the fit even if there is a prevFitParams to get the center of the galaxy%%%
    if fixToCenter
        ctrX = size(img, 2) / 2;
        ctrY = size(img, 1) / 2;
        if stgs.useSubpixelCtr
            ctrX = ctrX + 0.5;
            ctrY = ctrY + 0.5;
        end
        muFix = [ctrX, ctrY];
    else
        muFix = [];
    end

    starMaskLevel = 0;
    noiseMaskLevel = 0;
    masked = [];
    doneFitting = false;
    imgBeforeMasking = img;
    while ~doneFitting
        [muFit, covarFit, nonConvFlag, likFinal, gxyParams, bulgeMask] = ...
            iterFitGauss(img, stgs, muFix, gxyParams);
        ctrDrift = sqrt((muFit(1) - imgCtrC)^2 + (muFit(2) - imgCtrR)^2);
        [t1, t2, t3, curMajAxsLen] = findCovarElpsAxes(covarFit, likCutoff, size(img));
        if ctrDrift > ctrDriftThres
            if starMaskLevel == 0 % no star mask applied yet
                if isempty(starMask)
                    gxyParams.warnings = [gxyParams.warnings
                        'preprocessing:largeCenterDriftAndNoStarMaskAvailable'];
                    doneFitting = true;
                else
                    masked = (starMask == otherObjVal);
                    img(masked) = 0;
                    starMaskLevel = 1;
                end
            elseif starMaskLevel == 1 % already tried conservative mask
                masked = masked | (starMask == extRgnVal);
                img(masked) = 0;
                starMaskLevel = 2;
            elseif starMaskLevel == 2 % already tried aggressive mask
                masked = masked | (starMask ~= gxyRgnVal);
                img(masked) = 0;
                starMaskLevel = 3;
            elseif starMaskLevel == 3 % already tried last-resort mask
                doneFitting = true;
                starMaskLevel = 4;
            else
                error('internal error: unrecognized star mask level');
            end
        elseif curMajAxsLen > max(size(img))
            if noiseMaskLevel == 0 % no noise mask applied yet
                if isempty(starMask)
                    gxyParams.warnings = [gxyParams.warnings
                        'preprocessing:excessiveMajAxsLenAndNoNoiseMaskAvailable'];
                    doneFitting = true;
                else
                    masked = (starMask ~= gxyRgnVal) & (starMask ~= extRgnVal);
                    img(masked) = 0;
                    noiseMaskLevel = 1;
                end
            elseif noiseMaskLevel == 1 % already tried conservative mask
                masked = (starMask ~= gxyRgnVal);
                img(masked) = 0;
                noiseMaskLevel = 2;
            elseif noiseMaskLevel == 2 % already tried aggressive mask
                doneFitting = true;
                noiseMaskLevel = 3;
            else
                error('internal error: unrecognized noise mask level');
            end
        else
            doneFitting = true;
        end
    end
    img = imgBeforeMasking;

    [majAxsVec, majAxsAngle, axisRatio, majAxsLen] = ...
        findCovarElpsAxes(covarFit, likCutoff, size(img));
    %     minAxsLen = majAxsLen * axisRatio;
    semiMajAxsLen = majAxsLen / 2;

    starMaskLevelNames = {'none', 'conservative', 'aggressive', 'aggressive-exclusive', 'FAIL'};
    noiseMaskLevelNames = {'none', 'conservative-exclusive', 'aggressive-exclusive', 'FAIL'};
    if isempty(starMask)
        gxyParams.starMaskUsed = 'unavailable';
        gxyParams.noiseMaskUsed = 'unavailable';
    else
        gxyParams.starMaskUsed = starMaskLevelNames{starMaskLevel+1};
        gxyParams.noiseMaskUsed = noiseMaskLevelNames{noiseMaskLevel+1};
    end

    gxyParams.diskAxisRatio = axisRatio;
    gxyParams.diskMinAxsLen = majAxsLen * axisRatio;
    gxyParams.diskMajAxsLen = majAxsLen;
    gxyParams.diskMajAxsAngleRadians = majAxsAngle;
    gxyParams.iptCtrXY = muFit;
    gxyParams.iptSz = [nRows nCols];
    gxyParams.muDist = sqrt((muFit(1) - imgCtrC)^2 + (muFit(2) - imgCtrR)^2);
    gxyParams.muDistProp = gxyParams.muDist / min([nRows nCols]/2);
    gxyParams.wtdLik = sum(img(:) .* likFinal(:));
    gxyParams.likOfCtr = likFinal(round(imgCtrR), round(imgCtrC));
    gxyParams.brtUnifScore = getBrtUnifScore(img, likFinal, likCutoff, 10);
    gxyParams.gaussLogLik = sum(log(likFinal(:))); % remove this?

    % try to determine whether preprocessing went wrong by stretching in
    % between two bright sources
    % likInner = likFinal >= 10^-3.5;
    % if nnz(likInner) == 0
    %     likInner = likFinal >= 10^-4;
    % end
    innerLikRadiusRatio = 0.1;
    % in some cases with bad preprocessing, the ellipse is stretched between
    % two bright sources.  We may be able to measure this by comparing the
    % brightness in an inner elliptical contour with the brightness in an outer
    % elliptical contour.
    % the likelihood values already give us elliptical contours, so we use them
    % to find the cutoffs for the inner and outer regions
    innerLikContourPosn = (innerLikRadiusRatio * semiMajAxsLen) * majAxsVec;
    innerLikContourPosn = round([size(img, 1) - (muFit(2) + innerLikContourPosn(2)),...
        muFit(1) + innerLikContourPosn(1)]);
    innerLikContourVal = likFinal(innerLikContourPosn(1), innerLikContourPosn(2));
    likInner = likFinal >= innerLikContourVal;
    likOuter = likFinal >= 10^-9 & ~likInner;
    brtOuter = sort(img(likOuter), 'descend');
    contourBrtRatio = mean(img(likInner)) / mean(brtOuter(1:nnz(likInner)));

    if getenv('WRITEBULGEMASK')
        imwrite(bulgeMask, [outputPath '-D1_bulgeMask.png']);
    end

    % figure; imshow(img .* likOuter);

    % brtInner = sort(img(likInner), 'descend'); brtInner = brtInner(1:100);
    % brtOuter = brtOuter(1:100);
    % contourBrtRatio = mean(brtInner) / mean(brtOuter)

    gxyParams.contourBrtRatio = contourBrtRatio;

    % if abs(majAxsAngle) <= pi/4
    %     cropRad = abs(fzero(@(x)(mvnpdf([x, 0], [0, 0], covarFit) - likCutoff), 100))
    %     cropRad = majAxsLen * cos(majAxsAngle)
    % else
    %     cropRad = abs(fzero(@(x)(mvnpdf([0, x], [0, 0], covarFit) - likCutoff), 100))
    %     cropRad = abs(majAxsLen * sin(majAxsAngle))
    % end
    % cropRad = ceil(cropRad);

else %using elps file
%if ~isempty(prevFitParams)
    majAxsLen = prevFitParams.diskMajAxsLen;
    semiMajAxsLen = majAxsLen / 2;
    majAxsAngle = prevFitParams.diskMajAxsAngle;
    axisRatio = prevFitParams.diskAxisRatio;
    if ~recomputeCenter
        muFit = prevFitParams.muFit;
        gxyParams.iptCtrXY = muFit;
    end
    covarFit = prevFitParams.covarFit;
    contourBrtRatio = prevFitParams.contourBrtRatio;
    
    gxyParams.diskAxisRatio = axisRatio;
    gxyParams.diskMinAxsLen = majAxsLen * axisRatio;
    gxyParams.diskMajAxsLen = majAxsLen;
    gxyParams.diskMajAxsAngleRadians = majAxsAngle;
    %gxyParams.iptCtrXY = muFit;
    gxyParams.iptSz = [nRows nCols];
    gxyParams.muDist = sqrt((muFit(1) - imgCtrC)^2 + (muFit(2) - imgCtrR)^2);
    gxyParams.muDistProp = gxyParams.muDist / min([nRows nCols]/2);
    gxyParams.wtdLik = prevFitParams.wtdLik;
    gxyParams.likOfCtr = prevFitParams.likOfCtr; %not sure if this should be removed if we always find center
    gxyParams.brtUnifScore = prevFitParams.brtUnifScore;
    gxyParams.gaussLogLik = prevFitParams.gaussLogLik;
    gxyParams.contourBrtRatio = prevFitParams.contourBrtRatio;
    %WARNING: the prevFitParams might not have the following if it is old. The 2016 run will not have these
    gxyParams.diskMinAxsLen = prevFitParams.diskMinAxsLen;
    gxyParams.badBulgeFitFlag = prevFitParams.badBulgeFitFlag;
    gxyParams.bulgeAxisRatio = prevFitParams.bulgeAxisRatio;
    gxyParams.bulgeMajAxsLen = prevFitParams.bulgeMajAxsLen;
    gxyParams.bulgeMajAxsAngle = prevFitParams.bulgeMajAxsAngle;
    gxyParams.bulgeAvgBrt = prevFitParams.bulgeAvgBrt;
    gxyParams.bulgeDiskBrtRatio = prevFitParams.bulgeDiskBrtRatio;
    gxyParams.numElpsRefits = prevFitParams.numElpsRefits;
    gxyParams.bulgeMajAxsLen = prevFitParams.bulgeMajAxsLen;
    gxyParams.diskMajAxsAngleRadians = prevFitParams.diskMajAxsAngleRadians;
    gxyParams.starMaskUsed = prevFitParams.starMaskUsed;
    gxyParams.noiseMaskUsed = prevFitParams.noiseMaskUsed;
end
fitParams.diskMajAxsLen = majAxsLen;
fitParams.diskMajAxsAngle = majAxsAngle;
fitParams.diskAxisRatio = axisRatio;
fitParams.muFit = muFit;
fitParams.covarFit = covarFit;
fitParams.gaussLogLik = gxyParams.gaussLogLik;
fitParams.wtdLik = gxyParams.wtdLik;
fitParams.likOfCtr = gxyParams.likOfCtr;
fitParams.brtUnifScore = gxyParams.brtUnifScore;
fitParams.contourBrtRatio = contourBrtRatio;
%added the 2016 run might not have these in it's _elps files
fitParams.diskMinAxsLen = gxyParams.diskMinAxsLen;
fitParams.badBulgeFitFlag = gxyParams.badBulgeFitFlag;
fitParams.bulgeAxisRatio = gxyParams.bulgeAxisRatio;
fitParams.bulgeMajAxsLen = gxyParams.bulgeMajAxsLen;
fitParams.bulgeMajAxsAngle = gxyParams.bulgeMajAxsAngle;
fitParams.diskMajAxsAngleRadians = gxyParams.diskMajAxsAngleRadians;
fitParams.starMaskUsed = gxyParams.starMaskUsed;
fitParams.noiseMaskUsed = gxyParams.starMaskUsed;
fitParams.bulgeAvgBrt = gxyParams.bulgeAvgBrt;
fitParams.bulgeDiskBrtRatio = gxyParams.bulgeDiskBrtRatio;
fitParams.numElpsRefits = gxyParams.numElpsRefits;

clear imgCtrR
clear imgCtrC

[rInputPad, cInputPad] = detectConstantImagePadding(img);
if any(rInputPad) || any(cInputPad)
    fprintf('padding detected in input image (rows: %d, %d; cols: %d, %d)\n',...
        rInputPad(1), rInputPad(2), cInputPad(1), cInputPad(2));
    fprintf(['input padding will be replaced with replicated values to '...
        'avoid unsharp-mask artifacts.\n']);
    orig_size = size(img);
    img = img(1+rInputPad(1) : end-rInputPad(2), ...
        1+cInputPad(1) : end-cInputPad(2));
    img = padarray(img, [rInputPad(1) cInputPad(1)], 'replicate', 'pre');
    img = padarray(img, [rInputPad(2) cInputPad(2)], 'replicate', 'post');
    assert(all(size(img) == orig_size))
end


% split into preprocessGalfitImage.m
if getenv('USEGALFITRESIDUAL')
    img = preprocessGalfitImage(img,outputPath,gxyParams);
end

fromOrigImg = ones(size(img));
imgForUsm = img;

if useDeProjectStretch
    cropRad = ceil(semiMajAxsLen);
% otherwise...
% If we aren't doing a de-project stretch, we need to crop the image to
% inscribe an oblique ellipse.
% We use the parametric equations x(t), y(t), that describe the x- and y-
% coordinates of an ellipse for t in [0, 2*pi].  Derivatives of these are 
% used to find the maximum x- and y- distances from the center of the 
% ellipse, which we in turn use to find the "radius" (side half-length) 
% of the smallest square that contains all points in the ellipse 
% (likelihood contour of the fitted 2D Gaussian)
elseif abs(majAxsAngle) <= pi/4
    semiMinAxsLen = semiMajAxsLen * axisRatio;
    t = acos((semiMajAxsLen * cos(majAxsAngle)) / ...
        sqrt(semiMajAxsLen^2 * cos(majAxsAngle)^2 + ...
        semiMinAxsLen^2 * sin(majAxsAngle)^2));
    if mod(majAxsAngle, 2*pi) < pi
        % use the opposite angle, which still gives the same orientation
        % for the major axis of the ellipse, but restores the property of
        % pointing toward the positive y-axis
        t = -t;
    end
    xMax = abs(semiMajAxsLen * cos(t) * cos(majAxsAngle) - ...
        semiMinAxsLen * sin(t) * sin(majAxsAngle));
    cropRad = ceil(xMax);
else
    semiMinAxsLen = semiMajAxsLen * axisRatio;
	t = acos((semiMajAxsLen * sin(majAxsAngle)) / ...
        sqrt(semiMinAxsLen^2 * cos(majAxsAngle)^2 + ...
        semiMajAxsLen^2 * sin(majAxsAngle)^2));
	if (mod(majAxsAngle, 2*pi) > pi/2) && (mod(majAxsAngle, 2*pi) < 3*pi/2)
        % same type of adjustment, except ensuring that the angle points
        % toward the positive x-axis
        t = -t;
    end
	yMax = abs(semiMajAxsLen * cos(t) * sin(majAxsAngle) + ...
        semiMinAxsLen * sin(t) * cos(majAxsAngle));
	cropRad = ceil(yMax);
end

fitLik = calcLikelihoods([nRows, nCols], muFit, covarFit);

if ~useDeProjectStretch && plotFlag
    figure;
    hold on
    imagesc(flipud(img)); axis image; colormap gray
    % contour(flipud(fitLik), '-r')
    contour(flipud(fitLik), [10^-9, 10^-9], '-g', 'LineWidth', 2)
%     rectangle('Position', [muFit(1) - cropRad, muFit(2) - cropRad, ...
%         2*cropRad+1, 2*cropRad+1], 'EdgeColor', 'b', 'LineWidth', 2)
    axis off
    hold off
end

exactCtrC = muFit(1);
exactCtrR = size(img, 1) - muFit(2) + 1;
ctrC = round(exactCtrC);
ctrR = round(exactCtrR);
ctrAdjAmt = 0.5;

%%%%%%%%%%%%%

% ----- 1/31/20 - Matthew
% Code to print autocrop coordinates to file for galfit use
% Copying and pasting the crop code here shouldn't affect the later crop
% but I'll use different variables just in case.

if stgs.useSubpixelCtr
    xStart = (round(exactCtrR)-cropRad);
    xEnd = (round(exactCtrR)+cropRad);
    yStart = (round(exactCtrC)-cropRad);
    yEnd = (round(exactCtrC)+cropRad);
else
    xStart = (ctrR-cropRad);
    xEnd = (ctrR+cropRad);
    yStart = (ctrC-cropRad);
    yEnd = (ctrC+cropRad);
end

% check if the crop box would go out of range; if so, add more padding

cropRem = min([xStart - 1, size(img, 1) - xEnd,...
    yStart - 1, size(img, 2) - yEnd]);

if cropRem < 0

    xStart = xStart + (-cropRem);
    yStart = yStart + (-cropRem);

end

crop_coord_file = fopen([outputPath '_crop_coord.txt'],'at');
% fprintf(crop_coord_file, ['Path to output (for debugging):' outputPath '\n']);
fprintf(crop_coord_file, '%d\n', xStart, xEnd, yStart, yEnd);
fprintf(crop_coord_file, '\n');
fclose(crop_coord_file);

disp('Wrote autocrop coordinates to file.')

clear xStart
clear xEnd
clear yStart
clear yEnd

%%%%%%%%%%%%%%


% fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%         mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)
if useDeProjectStretch
    % pad image matrix so that fitted mean is at the center
    
    if stgs.useSubpixelCtr
        % rPad = round(exactCtrR - imgCtrR);
        rPad = (round(exactCtrR) - 1) - (nRows - round(exactCtrR));
    else
        rPad = (ctrR - 1) - (nRows - ctrR);
    end
%     fprintf('nRows = %2.4f ctrR = %2.4f rPad = %2.4f\n', nRows, ctrR, rPad);
    if rPad < 0
        rPadDir = 'pre';
        exactCtrR = exactCtrR + abs(rPad);
    else
        rPadDir = 'post';
    end
    img = padarray(img, [abs(rPad) 0], 0, rPadDir);
    fromOrigImg = padarray(fromOrigImg, [abs(rPad) 0], 0, rPadDir);
    imgForUsm = padarray(imgForUsm, [abs(rPad) 0], 'replicate', rPadDir);
%     imwrite(imgForUsm, [run_id '_deproject_03.png']);
    
    if stgs.useSubpixelCtr
        % cPad = round(exactCtrC - imgCtrC);
        cPad = (round(exactCtrC) - 1) - (nCols - round(exactCtrC));
    else
        cPad = (ctrC - 1) - (nCols - ctrC);
    end
%     fprintf('nCols = %2.4f ctrC = %2.4f cPad = %2.4f\n', nCols, ctrC, cPad);
    if cPad < 0
        cPadDir = 'pre';
        exactCtrC = exactCtrC + abs(cPad);
    else
        cPadDir = 'post';
    end
    img = padarray(img, [0 abs(cPad)], 0, cPadDir);
    fromOrigImg = padarray(fromOrigImg, [0 abs(cPad)], 0, cPadDir);
    imgForUsm = padarray(imgForUsm, [0 abs(cPad)], 'replicate', cPadDir);
%     imwrite(imgForUsm, [run_id '_deproject_04.png']);
    
%     fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%         mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)
else
    % pad image so that cropped region doesn't go out of bounds
    cPadAmt = max([0, 1 - (ctrC - cropRad), (ctrC + cropRad) - nCols]);
    rPadAmt = max([0, 1 - (ctrR - cropRad), (ctrR + cropRad) - nRows]);
    img = padarray(img, [rPadAmt, cPadAmt], 0);
    fromOrigImg = padarray(fromOrigImg, [rPadAmt, cPadAmt], 0);
    imgForUsm = padarray(imgForUsm, [rPadAmt, cPadAmt], 'replicate');
%     imwrite(imgForUsm, [run_id '_deproject_05.png']);
    ctrC = ctrC + cPadAmt; ctrR = ctrR + rPadAmt;
    exactCtrC = exactCtrC + cPadAmt; exactCtrR = exactCtrR + rPadAmt;
    
%     fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%         mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)
end

if useDeProjectStretch
    assert(all(size(img) == size(imgForUsm)))
    % rotate the image so that the major axis of the ellipse defined by the
    % fitted covariance matrix is aligned with the y-axis
    rotAngle = (pi/2 - majAxsAngle) * (180/pi);
    exactCtrX = exactCtrC - (size(img,2)/2 + ctrAdjAmt);
    exactCtrY = (size(img, 1) - exactCtrR + 1) - (size(img,1)/2 + ctrAdjAmt);
%     fprintf('ctrR = %2.2f exactCtrY = %2.2f ctrC = %2.2f exactCtrX = %2.2f\n',...
%         size(img,1)/2, exactCtrY, size(img,2)/2, exactCtrX)
    rotAngleDeg = rotAngle * (pi/180);
    rotMtx = [cos(rotAngleDeg), -sin(rotAngleDeg); sin(rotAngleDeg), cos(rotAngleDeg)];
    xy_rotated = rotMtx * [exactCtrX; exactCtrY]; 
    exactCtrX = xy_rotated(1);
    exactCtrY = xy_rotated(2);
%     fprintf('ctrR = %2.2f exactCtrY = %2.2f ctrC = %2.2f exactCtrX = %2.2f\n',...
%         size(img,1)/2, exactCtrY, size(img,2)/2, exactCtrX)
    exactCtrR = (size(img, 1) - (exactCtrY + (size(img,1)/2 + ctrAdjAmt)) + 1);
    exactCtrC = exactCtrX + (size(img,2)/2 + ctrAdjAmt);
    clear exactCtrX exactCtrY rotAngleDeg xy_rotated
    % TODO: determine if image size should be from original image or
    % rotated image (or if it matters since we currently adjust later)
%     fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%         mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)
    preRotateImgSz = size(img);
    img = imrotate(img, rotAngle, 'bilinear');
%     imgDbg = img;
%     rszAmt = [1 1] .* (mod(size(img), 2) == 0); % Octave
%     img = imresize(img, size(img) + rszAmt); % Octave
    fromOrigImg = imrotate(fromOrigImg, rotAngle, 'bilinear');
%     exactCtrR = exactCtrR * (size(img,1)/preRotateImgSz(1));
    exactCtrR = ((exactCtrR - ctrAdjAmt) * (size(img,1)/preRotateImgSz(1))) + ctrAdjAmt;
%     exactCtrC = exactCtrC * (size(img,2)/preRotateImgSz(2));
    exactCtrC = ((exactCtrC - ctrAdjAmt) * (size(img,2)/preRotateImgSz(2))) + ctrAdjAmt;
    clear preRotateImgSz
%     fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%         mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)

    % see written notes for calculation of usmImgPadAmt
    rotAngleForPadAmt = rotAngle * (pi/180);
    if rotAngleForPadAmt > pi/2
        rotAngleForPadAmt = rotAngleForPadAmt - pi/2;
    end
    usmImgPadAmt = ceil(max(size(imgForUsm)) * ...
        cos(rotAngleForPadAmt) * sin(rotAngleForPadAmt));
    imgForUsm = padarray(imgForUsm, [1 1] * usmImgPadAmt, 'replicate');
%     imwrite(imgForUsm, [run_id '_deproject_06.png']);
    imgForUsm = imrotate(imgForUsm, rotAngle, 'bilinear');
%     imwrite(imgForUsm, [run_id '_deproject_07.png']);
    usmImgCropAmt = round(size(imgForUsm) - size(img));
    % image rotate dimensions should be odd on both dimensions, so the
    % differences of the dimensions should be even
    % assert(all(mod(usmImgCropAmt, 2) == 0)); % This line causes csv writing errors. (tested 10/18/2017, not working)
    usmImgCropAmt = usmImgCropAmt / 2;
    imgForUsm = imgForUsm(round(1+usmImgCropAmt(1)):round(end-usmImgCropAmt(1)), ...
        round(1+usmImgCropAmt(2)):round(end-usmImgCropAmt(2)));
%     imwrite(imgForUsm, [run_id '_deproject_08.png']);
    assert(all(size(imgForUsm) == size(img)));
    
    % stretch the x-direction of the image so that the ellipse defined by
    % the fitted covariance matrix is isotropic
    oldNCols = size(img, 2);
    img = imresize(img, [size(img, 1) round(size(img, 2) / axisRatio)]);
    fromOrigImg = imresize(fromOrigImg, ...
        [size(fromOrigImg, 1) round(size(fromOrigImg, 2) / axisRatio)]);
    imgForUsm = imresize(imgForUsm, ...
        [size(imgForUsm, 1) round(size(imgForUsm, 2) / axisRatio)]);
%     imwrite(imgForUsm, [run_id '_deproject_09.png']);
%     figure; imshow(img); title('stretched image');
    nCols = size(img, 2);
%     exactCtrC = exactCtrC * (nCols/oldNCols);
    exactCtrC = ((exactCtrC - ctrAdjAmt) * (nCols/oldNCols)) + ctrAdjAmt;
    clear nCols oldNCols
%     fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%         mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)

    ctrR = round(size(img, 1) / 2);
    ctrC = round(size(img, 2) / 2);
end

assert(all(size(img) == size(imgForUsm)));
% fprintf('preprocessImage imgSz = %s\n', mat2str(size(img)))
% crop the image
if stgs.useSubpixelCtr
    rStart = (round(exactCtrR)-cropRad);
    rEnd = (round(exactCtrR)+cropRad);
    cStart = (round(exactCtrC)-cropRad);
    cEnd = (round(exactCtrC)+cropRad);
else
    rStart = (ctrR-cropRad);
    rEnd = (ctrR+cropRad);
    cStart = (ctrC-cropRad);
    cEnd = (ctrC+cropRad);
end
% check if the crop box would go out of range; if so, add more padding
cropRem = min([rStart - 1, size(img, 1) - rEnd,...
    cStart - 1, size(img, 2) - cEnd]);
if cropRem < 0
    img = padarray(img, -cropRem * [1 1], 0);
    fromOrigImg = padarray(fromOrigImg, -cropRem * [1 1], 0);
    imgForUsm = padarray(imgForUsm, -cropRem * [1 1], 'replicate');
%     imwrite(imgForUsm, [run_id '_deproject_10.png']);
    ctrR = ctrR + (-cropRem);
    ctrC = ctrC + (-cropRem);
    exactCtrR = exactCtrR + (-cropRem);
    exactCtrC = exactCtrC + (-cropRem);
    rStart = rStart + (-cropRem);
    cStart = cStart + (-cropRem);
end

% if stgs.useSubpixelCtr
%     % If the current image size is even, then there isn't really a pixel at
%     % the center and so we don't want to crop equally in both directions 
%     % from the "center" pixel.
%     % We still keep at least cropRad pixels in both directions;
%     % in the case of an odd size we have an "extra" pixel in the center.
%     if mod(size(img, 1), 2) == 0
%         % since we number pixels from 1, the "center" pixel will be in the
%         % first half of the image
%         rStart = rStart + 1;
%     end
%     if mod(size(img, 2), 2) == 0
%         cStart = cStart + 1;
%     end
%     % TODO: investigate if better option when one dimension is odd and the
%     % other is even
% end
img = img(rStart:rEnd, cStart:cEnd);
fromOrigImg = fromOrigImg(rStart:rEnd, cStart:cEnd);
imgForUsm = imgForUsm(rStart:rEnd, cStart:cEnd);
% imwrite(imgForUsm, [run_id '_deproject_11.png']);
exactCtrR = exactCtrR - (rStart-1);
exactCtrC = exactCtrC - (cStart-1);
% if useDeProjectStretch
%     ctrRforUsm = round(size(imgForUsm, 1) / 2);
%     ctrCforUsm = round(size(imgForUsm, 2) / 2);
% else
%     ctrRforUsm = ctrR;
%     ctrCforUsm = ctrC;
% end
% imgForUsm = imgForUsm((ctrRforUsm-cropRad):(ctrRforUsm+cropRad), ...
%     (ctrCforUsm-cropRad):(ctrCforUsm+cropRad));

% figure; subplot(1, 2, 1); imshow(img); subplot(1, 2, 2); imshow(imgForUsm);

if ~isempty(resizeDims)
    oldImSz = size(img);
    if (resizeDims(1) == resizeDims(2)) && (size(img, 1) ~= size(img, 2))
        gxyParams.warnings = [gxyParams.warnings
             ['preprocessing:preResizeDimsNonSquare:', mat2str(size(img))]];
    end
    img = imresize(img, resizeDims, 'bilinear');
    fromOrigImg = imresize(fromOrigImg, resizeDims, 'bilinear');
    imgForUsm = imresize(imgForUsm, resizeDims, 'bilinear');
%     imwrite(imgForUsm, [run_id '_deproject_12.png']);
%     exactCtrR = exactCtrR * (size(img,1)/oldImSz(1));
    exactCtrR = ((exactCtrR - ctrAdjAmt) * (size(img,1)/oldImSz(1))) + ctrAdjAmt;
%     exactCtrC = exactCtrC * (size(img,2)/oldImSz(2));
    exactCtrC = ((exactCtrC - ctrAdjAmt) * (size(img,2)/oldImSz(2))) + ctrAdjAmt;
    clear oldImSz
end

% if medFiltRad > 0
%     img = medfilt2(img, [1 1] * (2*medFiltRad+1));
%     imgForUsm = medfilt2(imgForUsm, [1 1] * (2*medFiltRad+1)); 
% end

imgNoUsm = img;
% imgDbg = imgForUsm;
if stgs.unsharpMaskAmt > 0
    %Octave chkPts = imerode(fromOrigImg > 0.99, strel('square', 3)) > 0;
    % if we use fromOrigImg==1, there are floating-point comparison issues
    chkPts = (fromOrigImg > .99);
    chkPts([1 end], :) = false; chkPts(:, [1 end]) = false;
    chkPts = imerode(chkPts, true(25, 25));
    assert(sum(chkPts(:) > 0) > 0);
    % make sure imgForUsm is the same as img in non-padding areas
%     assert(all(abs(imgForUsm(chkPts) - img(chkPts)) < 10^-4));
%     imgDbg = imgForUsm;
%     [img, tmp, imgDbg] = unsharpMask(imgForUsm, stgs);
    img = unsharpMask(imgForUsm, stgs);
%     imwrite(img, [run_id '_deproject_13.png']);
%     imgDbg = img;
    img(img < 0) = 0; img(img > 1) = 1;
    img = img .* fromOrigImg;
end
% img = real(log(img+1));
% img(img < 0) = 0;
% if useUnsharpMask
%     img = unsharpMask(img);
%     img(img < 0) = 0; img(img > 1) = 1;
% end

if isempty(prevFitParams) && nonConvFlag
    gxyParams.fit_state = 'preprocessing failure (ellipse fit did not converge within iteration limit)';
end

fitParams.cropRad = cropRad;
if useDeProjectStretch
    fitParams.rotAngle = rotAngle * (pi/180);
else
    fitParams.rotAngle = 0;
end

% remove negative values (from interpolation artifacts?)
img(img < 0) = 0;
imgNoUsm(imgNoUsm < 0) = 0;

% fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%     mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)

% exactCtrR = exactCtrR - ctrAdjAmt;
% exactCtrC = exactCtrC - ctrAdjAmt;

% fprintf('imgSz = %s, ctrR = %2.2f exactCtrR = %2.2f ctrC = %2.2f exactCtrC = %2.2f\n',...
%     mat2str(size(img)), size(img,1)/2, exactCtrR, size(img,2)/2, exactCtrC)

% img = imgDbg;

% if lookForBar
%     % find the bar length/angle for the standardized image
%     endR = barInfo.iptHalfLength * sin(barInfo.iptAngle);
%     endC = barInfo.iptHalfLength * cos(barInfo.iptAngle);
%     endR = endR * (resizeDims(1)/(2*cropRad+1));
%     endC = endC * (resizeDims(2)/(2*cropRad+1));
% 
%     barInfo.stdzCtrR = resizeDims(1)/2;
%     barInfo.stdzCtrC = resizeDims(2)/2;
%     barInfo.stdzAngle = barInfo.iptAngle + (rotAngle * (pi/180));
%     barInfo.stdzHalfLength = sqrt(endR.^2 + endC.^2);
% end

gxyParams.fitParams = fitParams;

end

function sc = getBrtUnifScore(img, lik, likCutoff, nBins)
    incl = lik >= likCutoff;
    brt = img(incl);
    dists = lik(incl);
    nPts = length(dists);
%     figure; scatter(dists, brt);
    [sortedDists, sI] = sort(dists, 'descend');
    sortedBrt = brt(sI);
    qBins = floor([0:1:(nPts-1)] * (nBins/nPts));
    
%     figure; hist(qBins);
%     qBrt = accumarray(qBins'+1, sortedBrt, [], @mean);
%     figure; plot(1:nBins, qBrt, '-s');
%     chgVals = -diff(qBrt) ./ qBrt(1:end-1);
%     isGoodChg = chgVals >= 0;
%     nGoodChg = sum(isGoodChg);
%     badChgAmt = sum(chgVals .* ~isGoodChg);

    qDiffs = accumarray(qBins'+1, sortedBrt, [], @(x)(diff(quantile(x, [0.1 0.9]))));
%     figure; plot(1:nBins, qDiffs, '-s');
    sc = max(qDiffs);
end

function [rInputPad, cInputPad] = detectConstantImagePadding(img)
    is_constant_row = all(img == repmat(img(:, 1), [1 size(img, 2)]), 2);
    is_constant_col = all(img == repmat(img(1, :), [size(img, 1) 1]), 1);
    rPrePadEnd = find(~is_constant_row, 1, 'first');
    rPostPadStart = find(~is_constant_row, 1, 'last');
    cPrePadEnd = find(~is_constant_col, 1, 'first');
    cPostPadStart = find(~is_constant_col, 1, 'last');
    if isempty(rPrePadEnd) || isempty(rPostPadStart)
        error('image values are all the same');
    end
    % all cols are constant iff all rows are constant, so we don't have to 
    % test them
    rInputPad = [rPrePadEnd - 1, size(img, 1) - rPostPadStart];
    cInputPad = [cPrePadEnd - 1, size(img, 2) - cPostPadStart];
end
