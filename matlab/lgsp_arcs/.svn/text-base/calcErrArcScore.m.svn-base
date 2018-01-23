function [propInCone, wtdPropInCone, propConeFill, ptConeItxSc, errConeArea, xFitErr, propConeFillNear, tipUsed] = ...
    calcErrArcScore(clusMtx1, clusMtx2, fitRes1, fitRes1L, fitRes1H, ...
    fitRes2, ctrR, ctrC, plotFlag)
% experimental; not currently used

% th2Vals useful for 2rev log-spirals, otherwise should not be given
% tipUsed: -1 if weighted to low-theta tip, 1 if weighted to high-theta
%   tip, 0 if no tip-weighting

% minThCutoff = pi/4;
% maxThCutoff = pi/2;
minAlenCutoff = size(clusMtx2, 1) / 10;
% maxAlenCutoff = size(clusMtx2, 1) / 4;
maxAlenCutoff = inf;

if nargin < 9
    plotFlag = false;
end

if ~isempty(fitRes1L) && ~isempty(fitRes1H)
    useTipWtg = true;
else
    useTipWtg = false;
    assert(isempty(fitRes1L) && isempty(fitRes1H));
end

if isempty(fitRes1) || isempty(fitRes2)
    % this should only be used when calling from the command line
    warning('recalculating fit parameters; assuming 2rev allowed');
    [fitParams1, arcBounds1, sumSqErr, used2Rev, failed2rev, hasBadBounds, fitTh, ptErrs, errArcs] = ...
        fitLogSpiral(clusMtx1, ctrR, ctrC, true, [], [], [], [], plotFlag);
    fitRes1 = fitResultsToStruct(fitParams1, arcBounds, sumSqErr, errArcs, fitTh);
    
    [fitParamsL, arcBoundsL, sumSqErr, used2Rev, failed2rev, hasBadBounds, fitTh, ptErrs, errArcs] = ...
        fitLogSpiral(clusMtx1, ctrR, ctrC, true, -1, fitParams1, arcBounds1, [], plotFlag);
    fitRes1L = fitResultsToStruct(fitParamsL, arcBounds, sumSqErr, errArcs, fitTh);
    
    [fitParamsH, arcBoundsH, sumSqErr, used2Rev, failed2rev, hasBadBounds, fitTh, ptErrs, errArcs] = ...
        fitLogSpiral(clusMtx1, ctrR, ctrC, true, 1, fitParams1, arcBounds1, [], plotFlag);
    fitRes1H = fitResultsToStruct(fitParamsH, arcBounds, sumSqErr, errArcs, fitTh);
    
    [fitParams2, arcBounds2, sumSqErr, used2Rev, failed2rev, hasBadBounds, fitTh, ptErrs, errArcs] = ...
        fitLogSpiral(clusMtx2, ctrR, ctrC, true, [], [], [], [], plotFlag);
    fitRes2 = fitResultsToStruct(fitParams2, arcBounds, sumSqErr, errArcs, fitTh);
end

% [rows, cols, brt1] = find(clusMtx1);
% [th1, rh1] = cart2pol(cols - ctrC, -(rows - ctrR));
par1 = fitRes1.params;
bnds1 = fitRes1.bounds;
% errArcs = fitRes1.errArcs;

[rows, cols, brt2] = find(clusMtx2);
[th2, rh2] = cart2pol(cols - ctrC, -(rows - ctrR));
par2 = fitRes2.params;
bnds2 = fitRes2.bounds;
fitTh2 = fitRes2.fitTh;

% displayLgspPlot([errArcs; par1], [bnds2; bnds2; bnds2] - par1(1) + par2(1) - 4*pi, clusMtx2, ctrR, ctrC);
% displayLgspPlot([errArcs; par1], [bnds2; bnds2; bnds2] - par1(1) + par2(1) - 2*pi, clusMtx2, ctrR, ctrC);
% displayLgspPlot([errArcs; par1], [bnds2; bnds2; bnds2] - par1(1) + par2(1), clusMtx2, ctrR, ctrC);
% displayLgspPlot([errArcs; par1], [bnds2; bnds2; bnds2] - par1(1) + par2(1) + 2*pi, clusMtx2, ctrR, ctrC);
% displayLgspPlot([errArcs; par1], [bnds2; bnds2; bnds2] - par1(1) + par2(1) + 4*pi, clusMtx2, ctrR, ctrC);

% % bounds of the second arc in terms of the offset for the first, up to a
% % multiple of 2*pi (we will find the correct multiple below).
% % note: we can only use these bounds with the parameters for the first
% % spiral; shifting the theta values of the second arc invalidates the
% % parameters of the second arc (changing the theta-offset doesn't just
% % rotate the arc)
% bnds2Trans = bnds2 - par1(1) + par2(1);

% Get the second cluster's theta-values into a contiguous range.

%     % begin TEMP
%     th2 = mod(th2, 2*pi);
%     th2srt = sort(th2, 'ascend');
%     gaps = diff(th2srt);
%     zGap = th2srt(1) + 2*pi - th2srt(end);
%     [maxGapVal, maxGapInd] = max(gaps);
%     if maxGapVal > zGap 
%         % gap doesn't go through zero, we need to move the discontinuity to
%         % where the gap is
%         gapMid = mean(th2srt([maxGapInd maxGapInd + 1]));
%         th2 = mod(th2 - gapMid, 2*pi) + gapMid;
%     end
%     [min(th2) max(th2)]
%     % end TEMP

if isempty(fitTh2)
    th2 = mod(th2, 2*pi);
    th2srt = sort(th2, 'ascend');
    gaps = diff(th2srt);
    zGap = th2srt(1) + 2*pi - th2srt(end);
    [maxGapVal, maxGapInd] = max(gaps);
    if maxGapVal > zGap 
        % gap doesn't go through zero, we need to move the discontinuity to
        % where the gap is
        gapMid = mean(th2srt([maxGapInd maxGapInd + 1]));
        th2 = mod(th2 - gapMid, 2*pi) + gapMid;
    end
else
    th2 = fitTh2;
end

% % now align the theta-values with the bounds, then adjust for theta-offset
% th2Trans =  th2 + (min(bnds2Trans) - min(th2)) + par1(1);
% % Find the right multiple by arc-fitting
% % rotAdjs = [-4*pi, -2*pi, 0, 2*pi, 4*pi];
% rotAdjs = [-2*pi, 0, 2*pi];
% rotAdjErrs = arrayfun(@(x)(sum(lgspErr(@logSpiralFxn2Rev, par1, th2Trans + x, rh2, brt2).^2)), rotAdjs);
% [minErr, minIdx] = min(rotAdjErrs);
% rotAdjByFitErr = rotAdjs(minIdx);
% 
% thDists = arrayfun(@(x)(min(abs([bnds1(1) - (bnds2Trans + x), bnds1(2) - (bnds2Trans + x)]))), rotAdjs);
% [minDist, minIdx] = min(thDists);
% rotAdjByThDist = rotAdjs(minIdx);
% 
% if rotAdjByFitErr ~= rotAdjByThDist
% %     figure; imshow(0.5 * edge(clusMtx1) + clusMtx2);
%     warning(['fit-based rotational adjustment (%2.4f) different than '...
%         'theta-distance rotational adjustment (%2.4f)'], ...
%         rotAdjByFitErr, rotAdjByThDist);
% end
% 
% rotAdjFinal = rotAdjByThDist;
% bnds2Trans = bnds2Trans + rotAdjFinal;
% th2Trans = th2Trans + rotAdjFinal;
% % fprintf('rotation adjustment: %d * 2*pi\n', rotAdjs(minIdx) / (2*pi));

[bnds2Trans, rotAdj, minDist, th2Trans] = translateThVals(par1, bnds1, par2, bnds2, th2);

bndsCommon = intervalIntersect(bnds1, bnds2Trans);
overlapAlen = 0;
if ~isnan(bndsCommon)
%     warning('clusters have overlapping theta-ranges, extending nearby-point cutoff');
    overlapAlen = calcLgspArcLengths(par2, bndsCommon + par1(1) - par2(1));
%     check theta-rotations if this is re-enabled
%     if ( abs(diff(bndsCommon)) - min(abs([diff(bnds1) diff(bnds2Trans)])) ) > -10e-6
%         fprintf('one interval contained entirely within another.  Not using nearby points (tip weighting).');
%         useTipWtg = false;
%     end

%     if plotFlag
%         displayLgspOverlay(clusMtx2 + 0.5 * edge(clusMtx1), par2, ctrR, ctrC,...
%             bndsCommon + par1(1) - par2(1) -rotAdj);
%         title(sprintf('calculating overlapAlen (rotAdj = %2.4f)', rotAdj));
%     end
end

th2forAlen =  th2Trans + (min(bnds2) - min(th2Trans)) + par2(1);
bnds2TransNear = sort(bnds2Trans, 'ascend');
withinCutoff = true(size(th2Trans));
alen2 = calcLgspArcLengths(par2, bnds2);
numPts2 = length(th2Trans);
alenCutoff = alen2 / 2;
alenCutoff = max(alenCutoff, minAlenCutoff);
alenCutoff = min(alenCutoff, maxAlenCutoff);
alenCutoff = alenCutoff + overlapAlen;
if ~useTipWtg
    errArcs = fitRes1.errArcs;
    tipUsed = 0;
elseif mean(bnds2Trans) < mean(bnds1)
    errArcs = fitRes1L.errArcs;
    ptAlens = calcLgspArcLengths(repmat(par2, numPts2, 1),...
        [repmat(max(th2forAlen), numPts2, 1), th2forAlen] - par2(1));

    assert(max(ptAlens) - alen2 < 5) % should be nearly equal
%     if abs(diff(bnds2TransNear)) > nearThCutoff
%         bnds2TransNear(1) = bnds2TransNear(2) - nearThCutoff;
%         withinCutoff = (th2Trans >= bnds2TransNear(1) + par1(1));
%     end
    withinCutoff = ptAlens <= alenCutoff;
    bnds2TransNear(1) = min(th2Trans(withinCutoff)) - par1(1);
    tipUsed = -1;
else
    errArcs = fitRes1H.errArcs;
    ptAlens = calcLgspArcLengths(repmat(par2, numPts2, 1),...
        [repmat(min(th2forAlen), numPts2, 1), th2forAlen] - par2(1));
%     if (max(ptAlens) - alen2 > 10e-6)
%         figure; imshow(clusMtx2);
%         figure; scatter(th2, th2);
%         figure; scatter(th2forAlen, th2forAlen)
%     end
    assert(max(ptAlens) - alen2 < 5) % should be nearly equal
%     if abs(diff(bnds2TransNear)) > nearThCutoff
%         bnds2TransNear(2) = bnds2TransNear(1) + nearThCutoff;
%         withinCutoff = (th2Trans <= bnds2TransNear(2) + par1(1));
%     end
    withinCutoff = ptAlens <= alenCutoff;
    bnds2TransNear(2) = max(th2Trans(withinCutoff)) - par1(1);
    tipUsed = 1;
end
% withinCutoffImg = -ones(size(clusMtx2));
% withinCutoffImg(clusMtx2 > 0) = withinCutoff;
% figure; imagesc(withinCutoffImg); axis image

errConeArea = lgspArea(errArcs(1, :), errArcs(2, :), bnds2Trans + par1(1));
assert(~isnan(errConeArea))
assert(errConeArea < inf)
nInClus = length(th2Trans);
gtBnd1 = rh2 > logSpiralFxn2Rev(th2Trans, errArcs(1, :));
gtBnd2 = rh2 > logSpiralFxn2Rev(th2Trans, errArcs(2, :));
isInCone = xor(gtBnd1, gtBnd2);
nInCone = sum(isInCone);

xFitErr = sum(lgspErr(@logSpiralFxn2Rev, par1, th2Trans, rh2, brt2).^2) / nInClus;
if nInCone > 0
    propConeFill = nInCone / errConeArea;
else
    propConeFill = 0;
end
assert(~isnan(propConeFill))
propInCone = nInCone / nInClus;
wtdPropInCone = sum(brt2(isInCone)) / sum(brt2);

nPtsNear = sum(withinCutoff);
nInConeNear = sum(isInCone & withinCutoff);
errConeAreaNear = lgspArea(errArcs(1, :), errArcs(2, :), bnds2TransNear + par1(1));
if nInConeNear > 0
    propConeFillNear = nInConeNear / errConeAreaNear;
else
    propConeFillNear = 0;
end
assert(~isnan(propConeFillNear))
propInConeNear = nInConeNear / nPtsNear;

ptConeItxSc = nInConeNear / (errConeAreaNear + (nPtsNear - nInConeNear));

gtCtr = rh2 > logSpiralFxn2Rev(th2Trans, par1);
isInCone1 = xor(gtBnd1, gtCtr) & withinCutoff;
nInCone1 = sum(isInCone1);
coneArea1 = lgspArea(errArcs(1, :), par1, bnds2TransNear + par1(1));
isInCone2 = xor(gtBnd2, gtCtr) & withinCutoff;
nInCone2 = sum(isInCone2);
coneArea2 = lgspArea(par1, errArcs(2, :), bnds2TransNear + par1(1));

ptConeItxSc1 = nInCone1 / (coneArea1 + (nPtsNear - nInCone1));
ptConeItxSc2 = nInCone2 / (coneArea2 + (nPtsNear - nInCone2));

ptConeItxScBest = max([ptConeItxSc, ptConeItxSc1, ptConeItxSc2]);

% icClusSimls = full(icSimls(clusMtx1 > 0, clusMtx2 > 0));

% displayLgspPlot([errArcs], [bnds2Trans; bnds2Trans], 0.5 * clusMtx1 + clusMtx2, ctrR, ctrC, [], [], {}, [], 2);
if plotFlag
    clusMtx2Disp = clusMtx2;
    clusMtx2Disp(clusMtx2Disp > 0) = (withinCutoff * 0.8 + 0.2) .* brt2;
    displayLgspPlot([errArcs; par1], [bnds2Trans; bnds2Trans; bnds2Trans],...
        0.5 * edge(clusMtx1) + clusMtx2Disp, ctrR, ctrC, [], [], {}, {':', ':', '-'}, 2);
    title(sprintf(['cross fit MSE = %2.4f\n'...
        'cone area = %2.4f, proportion with points = %2.4f\n'...
        '# points in cluster = %d, # in cone = %d, proportion in cone = %2.4f, wtd prop in cone = %2.4f\n'...
        '# in cone (near) = %d, prop in cone (near) = %2.4f, cone area (near) = %2.4f (%2.4f + %2.4f), proportion with points (near) = %2.4f\n'...
        'intersection score (full) = %2.4f, intersection score (best) = %2.4f, overlapAlen = %2.4f'],...
        xFitErr, errConeArea, propConeFill, nInClus, nInCone, propInCone, wtdPropInCone,...
        nInConeNear, propInConeNear, errConeAreaNear, coneArea1, coneArea2, propConeFillNear, ptConeItxSc, ptConeItxScBest, overlapAlen));
%         max(icClusSimls(:)), mean(icClusSimls(:)), min(icClusSimls(:))));
end

ptConeItxSc = ptConeItxScBest;

% figure; polar(th2(isInCone), rh2(isInCone), 'c.'); hold on; polar(th2(~isInCone), rh2(~isInCone), 'm.'); hold off


% 
% bnds2Trans - par2(1)
% fprintf('^^^\n');
% % deal with non-uniqueness of theta-values.  The theta-values could be
% % different by a multiple of 2*pi.  Find the right multiple by arc-fitting
% % 
% rotAdjs = [-2*pi, 0, 2*pi];
% % thErrFxn = @(theta)(lgspErr(@logSpiralFxn2pi, par1, theta, rh2, brt2));
% rotAdjErrs = arrayfun(@(x)(sum(lgspErr(@logSpiralFxn2Rev, par1, th2 + x, rh2, brt2).^2)), rotAdjs)
% [minErr, minIdx] = min(rotAdjErrs);
% bnds2Adj = bnds2Trans + rotAdjs(minIdx);
% 
% displayLgspPlot([errArcs; par1], [bnds2Adj; bnds2Adj; bnds2Adj], clusMtx2, ctrR, ctrC);

% if par1(1) > par2(1)
%     thShift = sign(par1(2)) * 2*pi;
% else
%     thShift = 0;
% end
% thShift
% 
% start1 = mod(par1(1) + bnds1(1), 2*pi)
% start2 = mod(par2(1) + bnds2(1), 2*pi)
% 
% if par1(2) > 0 && (start1(1) > start2(1))
%     thShift = 2*pi;
% elseif par1(2) < 0 && (start1(1) < start2(1))
%     thShift = -2*pi;
% else
%     thShift = 0;
% end
% thShift
% 
% displayLgspPlot(errArcs, [bnds2; bnds2] - par1(1) + par2(1) + thShift, clusMtx2, ctrR, ctrC);

end