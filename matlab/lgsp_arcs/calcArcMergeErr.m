function [mergeErrRatio, anyFailed2rev] = ...
    calcArcMergeErr(clus1, params1, bounds1, err1, ...
    clus2, params2, bounds2, err2, ...
    ctrR, ctrC, barInfo, stgs, useTipWtg, thVals1, thVals2, plotFlag)

error(nargchk(12, 16, nargin));

% fprintf('calcArcMergeErr: ctrR=%2.2f, ctrC=%2.2f\n', ctrR, ctrC)

% TODO: improve handling of precomputed fit information (params, bounds,
%   thVals, etc, maybe move them into a struct since they're handled
%   together anyway)

if ~isempty(params1) && isempty(thVals1)
    error('if params1 specified, thVals1 must also be specified');
end
if ~isempty(params2) && isempty(thVals2)
    error('if params2 specified, thVals2 must also be specified');
end

if nargin < 13 || isempty(useTipWtg)
    useTipWtg = false;
end

if nargin < 16 || isempty(plotFlag)
    plotFlag = false;
end

global failed2revDuringMerge;

anyFailed2rev = false;
imgSz = size(clus1);

% if useTipWtg
%     wts1 = nonzeros(clus1);
%     wts2 = nonzeros(clus2);
%     [wts1, wts2] = calcTipWtsDbl(nonzeros(clus1), 
% end

if useTipWtg && (nargin < 15 || isempty(params1) || isempty(params2) || ...
        isempty(thVals1) || isempty(thVals2))
    error(['if tip weighting used, parameters and theta-values'...
        ' must be specified for each cluster']);
end

ignoreJaggedBoundaryPixelsDuringMerges = stgs.ignoreJaggedBoundaryPixelsDuringMerges;
outlierSd = 3;

balClusWtsInMerging = stgs.balClusWtsInMerging;

if useTipWtg
    [wts1, wts2, isOverlap] = ...
        calcTipWtsDbl(nonzeros(clus1), thVals1, params1, bounds1,...
        nonzeros(clus2), thVals2, params2, bounds2, imgSz);
%     if isOverlap && plotFlag
%         displayLgspOverlay(clus1 + clus2, [params1; params2], ctrR, ctrC, [bounds1; bounds2]);
%         title('before refit');
%     end
    clus1(clus1 > 0) = wts1;
    clus2(clus2 > 0) = wts2;
    if isOverlap % tip-weighting changed, so we need to refit
        [params1, bounds1, err1, used2Rev1, failed2rev1] = ...
            fitLogSpiral(clus1, ctrR, ctrC, stgs, 0, params1, bounds1);
        [params2, bounds2, err2, used2Rev2, failed2rev2] = ...
            fitLogSpiral(clus2, ctrR, ctrC, stgs, 0, params2, bounds2);
        anyFailed2rev = anyFailed2rev | failed2rev1 | failed2rev2;
%         if plotFlag
%             displayLgspOverlay(clus1 + clus2, [params1; params2], ctrR, ctrC, [bounds1; bounds2]);
%             title('after refit');
%         end
    end
end

% params1 = [];
% params2 = [];
if isempty(params1)
    [params1, bounds1, err1, used2Rev1, failed2rev1, tmp1, thVals1] = ...
        fitLogSpiral(clus1, ctrR, ctrC, stgs);
    anyFailed2rev = anyFailed2rev | failed2rev1;
end
if isempty(params2)
    [params2, bounds2, err2, used2Rev2, failed2rev2, tmp1, thVals2] = ...
        fitLogSpiral(clus2, ctrR, ctrC, stgs);
    anyFailed2rev = anyFailed2rev | failed2rev2;
end
if ignoreJaggedBoundaryPixelsDuringMerges
%     plotFlag = true;
    
    errs1 = getLgspFitErrs(stgs, clus1, ctrR, ctrC, params1, bounds1, thVals1);
    nerr1 = (errs1 - mean(errs1)) / std(errs1);
    isInlier1 = nerr1 <= outlierSd;
%     err1 = sum(errs1(isInlier1));
    
    errs2 = getLgspFitErrs(stgs, clus2, ctrR, ctrC, params2, bounds2, thVals2);
    nerr2 = (errs2 - mean(errs2)) / std(errs2);
    isInlier2 = nerr2 <= outlierSd;
%     err2 = sum(errs2(isInlier2));
    
    clus1(clus1 > 0) = clus1(clus1 > 0) .* isInlier1;
    clus2(clus2 > 0) = clus2(clus2 > 0) .* isInlier2;
    
    [params1, bounds1, err1, used2Rev1, failed2rev1, tmp, thVals1] = ...
        fitLogSpiral(clus1, ctrR, ctrC, stgs);
    [params2, bounds2, err2, used2Rev2, failed2rev2, tmp, thVals2] = ...
        fitLogSpiral(clus2, ctrR, ctrC, stgs);
    anyFailed2rev = anyFailed2rev | failed2rev1 | failed2rev2;
end

wtSum1 = sum(clus1(:));
wtSum2 = sum(clus2(:));
assert(wtSum1 ~= 0)
assert(wtSum2 ~= 0)
if balClusWtsInMerging
    wtSum = wtSum1 + wtSum1;
    clus1 = clus1 * (wtSum/wtSum1);
    clus2 = clus2 * (wtSum/wtSum2);
    
    wtSum1 = sum(clus1(:));
    wtSum2 = sum(clus2(:));
    assert(abs(wtSum1 - wtSum2) < 10^-6);  
    
    [params1, bounds1, err1, used2Rev1, failed2rev1, tmp, thVals1] = ...
        fitLogSpiral(clus1, ctrR, ctrC, stgs);
    [params2, bounds2, err2, used2Rev2, failed2rev2, tmp, thVals2] = ...
        fitLogSpiral(clus2, ctrR, ctrC, stgs);
    anyFailed2rev = anyFailed2rev | failed2rev1 | failed2rev2;
end
clusCombined = clus1 + clus2;

[fitParams1, arcBounds1, errFrom1, u2r1, failed2rev1, hasBadBounds1, thvFrom1, ptErrs1] = ...
    fitLogSpiral(clusCombined, ctrR, ctrC, stgs, [], params1);

[fitParams2, arcBounds2, errFrom2, u2r2, failed2rev2, hasBadBounds2, thvFrom2, ptErrs2] = ...
    fitLogSpiral(clusCombined, ctrR, ctrC, stgs, [], params2);

anyFailed2rev = anyFailed2rev | failed2rev1 | failed2rev2;

% if errFrom1 ~= errFrom2
%     fprintf('[errFrom1, errFrom2] = [%2.4f, %2.4f]\n', errFrom1, errFrom2);
% end

if errFrom1 <= errFrom2
    sglArcErr = errFrom1;
    params = fitParams1;
    bounds = arcBounds1;
    thv = thvFrom1;
    used2Rev = u2r1;
    ptErrs = ptErrs1;
else
    sglArcErr = errFrom2;
    params = fitParams2;
    bounds = arcBounds2;
    thv = thvFrom2;
    used2Rev = u2r2;
    ptErrs = ptErrs2;
end
% figure; plot(sort(thv));
% figure; hold all; plot(sort(thVals1)); plot(sort(thVals2)); plot(sort(thv)); legend('thVals1', 'thVals2', 'thv');

% [rC, cC, brtC] = find(clusCombined);
% [thC, rhoC] = cart2pol(cC - ctrC, -(rC - ctrR));
[r1, c1, brt1] = find(clus1);
% [th1, rho1] = cart2pol(c1 - ctrC, -(r1 - ctrR));
% fprintf('err1 = %16.16f\n', err1);
% % TODO: determine whether 2rev used via the variable used2Rev from
% %   fitLogSpiral.  This will be easier once all the fit information (params,
% %   bounds, err, etc) is moved into a struct
% cerr1 = min(sum(lgspErr(@logSpiralFxn, params1, thVals1, rho1, brt1).^2),...
%     sum(lgspErr(@logSpiralFxn2Rev, params1, thVals1, rho1, brt1).^2));
% fprintf('calculated err1 = %16.16f\n', cerr1);
% assert(abs(err1 - cerr1) < 10^-6);
[r2, c2, brt2] = find(clus2);
% [th2, rho2] = cart2pol(c2 - ctrC, -(r2 - ctrR));
% fprintf('err2 = %16.16f\n', err2);
% cerr2 = min(sum(lgspErr(@logSpiralFxn, params2, thVals2, rho2, brt2).^2),...
%     sum(lgspErr(@logSpiralFxn2Rev, params2, thVals2, rho2, brt2).^2));
% fprintf('calculated err2 = %16.16f\n', cerr2)
% assert(abs(err2 - cerr2) < 10^-6);

% TODO: determine whether 2rev used via the variable used2Rev from
%   fitLogSpiral.  This will be easier once all the fit information (params,
%   bounds, err, etc) is moved into a struct
% Get the errors for the combined-fit arc with respect to the individual
% clusters
from1 = (clus1 > 0); from1 = from1(clusCombined > 0);
% figure; hold all; plot(sort(thVals1)); plot(sort(thv(from1))); title('th vals clus 1');
% figure; hold all; plot(sort(thVals2)); plot(sort(thv(~from1))); title(sprintf('th vals clus 2 diff = %2.4f', min(thVals2) - min(thv(~from1))));
% if used2Rev
%     serr = ...
%         lgspErr(@logSpiralFxn2Rev, params, thv, rhoC, brtC).^2;
% else
%     serr = ...
%         lgspErr(@logSpiralFxn, params, thv, rhoC, brtC).^2;
% end
serr = ptErrs;
% th1Img = clus1; th1Img(th1Img > 0) = thVals1;
% th2Img = clus2; th2Img(th2Img > 0) = thVals2;
% thvImg = clusCombined; thvImg(thvImg > 0) = thv; 
% errImg = clusCombined; errImg(errImg > 0) = serr; 
% figure; subplot(2, 2, 1); imagesc(th1Img); axis image; subplot(2, 2, 2); imagesc(th2Img); axis image; subplot(2, 2, 3); imagesc(thvImg); axis image; subplot(2, 2, 4); imagesc(errImg); axis image; title(sprintf('serr, 2rev = %d', used2Rev));
assert(abs(sum(serr) - sglArcErr) / numel(serr) < 10^-8); % make sure we're correctly replicating the error
serr1 = sum(serr(from1));
serr2 = sum(serr(~from1));
% serr1 = min(sum(lgspErr(@logSpiralFxn, params, thv(from1), rho1, brt1).^2),...
%     sum(lgspErr(@logSpiralFxn2Rev, params, thv(from1), rho1, brt1).^2));
% sum(lgspErr(@logSpiralFxn, params, thv(from1), rho1, brt1).^2) / sum(brt1)
% sum(lgspErr(@logSpiralFxn2Rev, params, thv(from1), rho1, brt1).^2) / sum(brt1)
% serr2 = min(sum(lgspErr(@logSpiralFxn, params, thv(~from1), rho2, brt2).^2),...
%     sum(lgspErr(@logSpiralFxn2Rev, params, thv(~from1), rho2, brt2).^2));
% sum(lgspErr(@logSpiralFxn, params, thv(~from1), rho2, brt2).^2) / sum(brt2)
% sum(lgspErr(@logSpiralFxn2Rev, params, thv(~from1), rho2, brt2).^2) / sum(brt2)

barBetter1 = false; barBetter2 = false; barBetterM = false;
if ~isempty(barInfo) && barInfo.barDetected
    % if the bar-fit error is lower, use that error instead
    bErr1 = calcBarFitErr(clus1, barInfo.stdzCtrR, ...
        barInfo.stdzCtrC, barInfo.stdzAngle, barInfo.stdzHalfLength);
    if bErr1 < err1
        barBetter1 = true;
        err1 = bErr1;
    end
    bErr2 = calcBarFitErr(clus2, barInfo.stdzCtrR, ...
        barInfo.stdzCtrC, barInfo.stdzAngle, barInfo.stdzHalfLength);
    if bErr2 < err2
        barBetter2 = true;
        err2 = bErr2;
    end
    errS = calcBarFitErr(clusCombined, barInfo.stdzCtrR, ...
        barInfo.stdzCtrC, barInfo.stdzAngle, barInfo.stdzHalfLength);
    if errS < sglArcErr
        barBetterM = true;
        sglArcErr = errS;
    end
    
    serr1 = min(serr1, bErr1);
    serr2 = min(serr2, bErr2);
end
mse1 = (err1 / wtSum1);
mse2 = (err2 / wtSum2);
origErrSum = (mse1 + mse2);
sglArcMse = sglArcErr / (wtSum1 + wtSum2);

mergeErrRatio = sglArcMse / origErrSum;
if stgs.balClusWtsInMerging
    if barBetter1 || barBetter2 || barBetterM
        mergeErrRatio = max((bErr1/sum(brt1))/mse1, ...
            (bErr2/sum(brt2))/mse2);
    else
        mergeErrRatio = max((serr1/sum(brt1))/mse1, ...
            (serr2/sum(brt2))/mse2);
    end
end

if plotFlag
    % some old experimental measures that aren't being used anymore
    % so unless we want detailed info for plotting, don't bother to
    % calculate these
    errRatio = sglArcErr / (err1 + err2);
    errRatioToMin = sglArcMse / min((err1 / wtSum1), (err2 / wtSum2));

    altMeanErrRatio = ((sglArcMse / mse1) + (sglArcMse / mse2)) / 2;
    altMeanErrRatio = sqrt((sglArcMse / mse1) * (sglArcMse / mse2));

    [r1, c1, brt1] = find(clus1);
    r1 = r1 - ctrR; c1 = c1 - ctrC;
    [theta1, rho1] = cart2pol(c1, -r1);
    theta1 = mod(theta1, 2*pi);

    [r2, c2, brt2] = find(clus2);
    r2 = r2 - ctrR; c2 = c2 - ctrC;
    [theta2, rho2] = cart2pol(c2, -r2);
    theta2 = mod(theta2, 2*pi);

    err12 = sum( (brt2 .* (rho2 - logSpiralFxn(theta2, params1))).^2 ) / wtSum2;
    err21 = sum( (brt1 .* (rho1 - logSpiralFxn(theta1, params2))).^2 ) / wtSum1;

    minCrossErr = min(err12, err21) / origErrSum;
    maxCrossErr = max(err12, err21) / origErrSum;
    meanCrossErr = mean([err12, err21]) / origErrSum;

    paChg1 = abs(params(2) - params1(2));
    paChg2 = abs(params(2) - params2(2));
    maxPaChg = max(paChg1, paChg2);
    minPaChg = min(paChg1, paChg2);
    meanPaChg = (paChg1 + paChg2) / 2;
end

% paramDist = sqrt(sum((params1 - params2).^2 .* [1/2 5 1/150]));

% % extnProp = 0.25;
% % len1 = abs(bounds1(2) - bounds1(1)); extnAmt1 = extnProp * len1;
% % len2 = abs(bounds2(2) - bounds2(1)); extnAmt2 = extnProp * len2;
% extnAmt1 = pi/8; extnAmt2 = pi/8;
% exBnds1 = bounds1 + [-extnAmt1 extnAmt1];
% exBnds2 = bounds2 + [-extnAmt2 extnAmt2];
% ctrX = ctrC;
% ctrY = imgSz(1) - ctrR + 1;
% arcMtx1 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params1, bounds1);
% arcExMtx1 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params1, exBnds1);
% arcMtx2 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params2, bounds2);
% arcExMtx2 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params2, exBnds2);
% 
% arcMatch12 = bwdist(arcMtx2); 
% arcMatch12 = arcMatch12((arcExMtx1 - arcMtx1) > 0);
% % arcMatch12 = arcMatch12(arcExMtx1 > 0);
% % arcMatch12 = min(arcMatch12);
% arcMatch12 = sort(arcMatch12, 'ascend'); arcMatch12 = mean(arcMatch12(1:min(10, end)));
% 
% arcMatch21 = bwdist(arcMtx1); 
% arcMatch21 = arcMatch21((arcExMtx2 - arcMtx2) > 0);
% % arcMatch21 = arcMatch21(arcExMtx2 > 0);
% % arcMatch21 = min(arcMatch21);
% arcMatch21 = sort(arcMatch21, 'ascend'); arcMatch21 = mean(arcMatch21(1:min(10, end)));
% 
% arcDist = max([arcMatch12, arcMatch21]);
% if isempty(arcDist)
%     arcDist = inf;
% end

% if sum(clus1(:) > 0) > 150 || sum(clus2(:) > 0) > 150
%     plotFlag = true;
% end
if plotFlag
%     overlay = displayLgspOverlay(clusCombined, params, ctrR, ctrC, bounds);
%     figure('Position', [0 0 800 800]); imshow(overlay);
%     title(sprintf(...
%         ['errRatio = %2.4f, meanErrRatio = %2.4f\n' ...
%         'err12 = %2.4f, err21 = %2.4f, mean = %2.4f\n' ...
%         'sglArcErr = %2.4f, indErrs = {%2.4f, %2.4f}'], ...
%         errRatio, meanErrRatio, err12, err21, (err12 + err21) / 2, ...
%         sglArcErr, err1, err2));
    
    clusPlotImg = zeros([size(clus1) 3]);
    clusPlotImg(:, :, 1) = clus1;
    clusPlotImg(:, :, 3) = clus2;
    clusPlotImg(:, :, :) = clusPlotImg(:, :, :) + repmat(edge(clus1 + clus2 > 0) * 0.5, [1 1 3]);
%     overlay = displayLgspOverlay(clusPlotImg, [params1; params2; params], ctrR, ctrC, [exBnds1; exBnds2; bounds]);
    overlay = displayLgspOverlay(clusPlotImg, [params1; params2; params], ctrR, ctrC, [bounds1; bounds2; bounds]);  
    figure('Position', [0 0 800 800]); imshow(overlay);
    titleStr = sprintf('mergeErrRatio = %2.4f (%2.4f, %2.4f, %2.4f, %2.4f)',...
        mergeErrRatio, serr1/sum(brt1), mse1, serr2/sum(brt2), mse2);
    if ~isempty(barInfo) && barInfo.barDetected
        titleStr = [titleStr sprintf('\nbar errors: %2.4f, %2.4f, %2.4f',...
            bErr1/sum(brt1), bErr2/sum(brt2), errS/(sum(brt1)+sum(brt2)) )];
    end
    title(titleStr);
%     saveas(gca, ['D:\matlab-output-local\temp\newMergeErrCalc\' mat2str(round(rem(now, 1) * 10^14)) '.png']);
%     close gcf
end

failed2revDuringMerge = anyFailed2rev;

end