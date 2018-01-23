function [clusMtxsMerged, nMerges] = mergeClusters(clusMtxs, ctrR, ctrC,...
    stgs, plotFlag)

allowArcBeyond2pi = stgs.allowArcBeyond2pi;

% this function has not been updated to use stgs.errRatioThres, etc, since
% this function is (was) experimental

if nargin < 5 || isempty(plotFlag)
    plotFlag = false;
end

% plotFlag = true;
methodNum = 12;

% tipWtSigma = pi/8;
useTipWtg = false;
inCtrRad = round(size(clusMtxs, 1) / 16);
% outCtrRad = round(size(clusMtxs, 1) / 8);
outCtrRad = inf;
maxClusThGap = pi/4; % maximum allowed theta-gap between clusters
maxClusThOverlap = pi/4; % maximum allowed theta-overlap between clusters

switch methodNum
    case 1
        clusDistFxn = @calcClusDistHybrid;
        stopThres = 1;
    case 2
        clusDistFxn = @maxCrossFitClusDist;
        stopThres = 50;
    case 3
        clusDistFxn = @meanCrossFitClusDist;
        stopThres = 50;
    case 4
        clusDistFxn = @worstPropConeFillClusDist;
        stopThres = 0.9;
    case 5
        clusDistFxn = @worstPropConeFillClusDist;
        stopThres = 0.8;
    case 6
        clusDistFxn = @errRatioAndConeFillClusDist;
        stopThres = 1;
    case 7
        clusDistFxn = @errRatioThresConeFillClusDist;
        stopThres = 0.9;
        meanErrRatioThres = 1;
    case 8
        clusDistFxn = @worstPropConeFillClusDist;
        stopThres = 0.9;
        useTipWtg = true;
    case 9
        baseClusDistFxn = @worstPropConeFillClusDist;
        clusDistFxn = @errRatioThresConeFillClusDist;
        stopThres = 0.9;
        meanErrRatioThres = 1;
        useTipWtg = true;
    case 10
        baseClusDistFxn = @worstItxScClusDist;
        clusDistFxn = @errRatioThresConeFillClusDist;
        stopThres = 0.9;
        meanErrRatioThres = 1;
        useTipWtg = true;
    case 11
        clusDistFxn = @worstItxScClusDist;
        stopThres = 0.8;
        useTipWtg = true;
    case 12
        clusDistFxn = @coneFillThresErrRatioClusDist;
        stopThres = 1;
        useTipWtg = true;
        maxConeFillDist = 0.9;
    otherwise
        error('unrecognized merge method number');
end

clusMtxsIn = clusMtxs;

tMerge = tic;

nClus = size(clusMtxs, 3);
imgSz = size(clusMtxs); imgSz = imgSz(1:2);

% TODO: this is recomputed and stored in allFitRes - change the earlier
% similarity function calcClusDistHybrid to use this instead
[lgspParams, lgspBounds, fitErrs] = ...
    fitLogSpiralsToClusters(clusMtxs, ctrR, ctrC, stgs);

emptyFitRes = struct('params', [], 'bounds', [], 'err', [], ...
    'errArcs', [], 'fitTh', [], 'fitWts', []);
allFitRes = repmat(emptyFitRes, nClus, 1);
allFitResL = repmat(emptyFitRes, nClus, 1);
allFitResH = repmat(emptyFitRes, nClus, 1);
for ii=1:1:nClus
    updateClusterFit(ii);
end

clusDists = inf * ones(nClus, nClus);

dists = zeros(size(clusMtxs(:, :, 1)));
for ii=1:1:nClus
    for jj=ii+1:1:nClus
        clusDists(ii, jj) = clusDistFxn(ii, jj);
    end
end

keepMerging = true;
deletedClusts = false(1, nClus);
nMerges = 0;
while keepMerging
    [vals, rows] = min(clusDists, [], 1);
    [val, cIdx1] = min(vals);
    cIdx2 = rows(cIdx1);
    
    if val <= stopThres
%         figure
%         subplot(1, 3, 1); imshow(clusMtxs(:, :, cIdx1));     
%         subplot(1, 3, 2); imshow(clusMtxs(:, :, cIdx2));
%         subplot(1, 3, 3); imshow(clusMtxs(:, :, cIdx1) + clusMtxs(:, :, cIdx2));
%         title(sprintf('val = %2.4f, sizes = %2.4f, %2.4f', ...
%             val, nnz(clusMtxs(:, :, cIdx1)), nnz(clusMtxs(:, :, cIdx2))));
        
        if plotFlag
            fprintf('merging with val = %2.4f\n', val);
            mergeImg = zeros([size(clusMtxs(:, :, cIdx1)) 3]);
            mergeImg(:, :, 2) = clusMtxs(:, :, cIdx1);
            mergeImg(:, :, 3) = clusMtxs(:, :, cIdx2);
            figure; imshow(mergeImg); title(sprintf('merge value %2.4f', val));
            if useTipWtg
                calcErrArcScore(clusMtxs(:, :, cIdx2), clusMtxs(:, :, cIdx1), allFitRes(cIdx2), allFitResL(cIdx2), allFitResH(cIdx2), allFitRes(cIdx1), ctrR, ctrC, true);
                calcErrArcScore(clusMtxs(:, :, cIdx1), clusMtxs(:, :, cIdx2), allFitRes(cIdx1), allFitResL(cIdx1), allFitResH(cIdx1), allFitRes(cIdx2), ctrR, ctrC, true);
            else
                calcErrArcScore(clusMtxs(:, :, cIdx2), clusMtxs(:, :, cIdx1), allFitRes(cIdx2), [], [], allFitRes(cIdx1), ctrR, ctrC, true);
                calcErrArcScore(clusMtxs(:, :, cIdx1), clusMtxs(:, :, cIdx2), allFitRes(cIdx1), [], [], allFitRes(cIdx2), ctrR, ctrC, true);
            end
        end

        clusMtxs(:, :, cIdx1) = clusMtxs(:, :, cIdx1) + clusMtxs(:, :, cIdx2);
        clusMtxs(:, :, cIdx2) = 0;
        deletedClusts(cIdx2) = true;
        clusDists(cIdx2, :) = inf;
        clusDists(:, cIdx2) = inf;
        nMerges = nMerges + 1;
        
%        	[fitParams, arcBounds, sumSqErr, fitTh, used2Rev, errArcs, fitWts] = ...
%             fitLogSpiral(clusMtxs(:, :, cIdx1), ctrR, ctrC, allowArcBeyond2pi);
%         allFitRes(cIdx1).params = fitParams;
%         allFitRes(cIdx1).bounds = arcBounds;
%         allFitRes(cIdx1).err = sumSqErr;
%         allFitRes(cIdx1).errArcs = errArcs;
%         allFitRes(cIdx1).fitTh = fitTh;
        updateClusterFit(cIdx1);
%         allFitRes(cIdx2) 
%         emptyFitRes
        allFitRes(cIdx2) = emptyFitRes;

        for idx=1:1:cIdx1-1
            if ~deletedClusts(idx)
                clusDists(idx, cIdx1) = clusDistFxn(idx, cIdx1);
            end
        end
        for idx=cIdx1+1:1:nClus
            if ~deletedClusts(idx)
                clusDists(cIdx1, idx) = clusDistFxn(cIdx1, idx);
            end
        end
    else
        if plotFlag
            mergeImg = zeros([size(clusMtxs(:, :, cIdx1)) 3]);
            mergeImg(:, :, 2) = clusMtxs(:, :, cIdx1);
            mergeImg(:, :, 3) = clusMtxs(:, :, cIdx2);
            figure; imshow(mergeImg); title(sprintf('first pair not merged (merge value %2.4f)', val));
        end
        keepMerging = false;
    end
%     showClustersFromMtxs(clusMtxs, imgSz);
%     pause
%     close all
end
clusMtxsMerged = clusMtxs(:, :, ~deletedClusts);

if plotFlag
    showClustersFromMtxs(clusMtxsIn, imgSz); 
    title(sprintf('%d input clusters', size(clusMtxsIn, 3)));
    showClustersFromMtxs(clusMtxsMerged, imgSz); 
    title(sprintf('%d clusters after merge', size(clusMtxsMerged, 3)));
end
    
fprintf('Time for arc-merging: \n');
toc(tMerge);

    function updateClusterFit(idx)
%         fprintf('updating cluster fit for %d\n', idx);
        toMin = -1;
        toMax = 1;
        
        [fitParams, arcBounds, err, used2Rev, failed2rev, hasBadBounds, fitTh, ptErrs, errArcs, fitWts] = ...
            fitLogSpiral(clusMtxs(:, :, idx), ctrR, ctrC, stgs);    
        allFitRes(idx) = struct('params', fitParams, 'bounds', arcBounds, ...
            'err', err, 'errArcs', errArcs, 'fitTh', fitTh, 'fitWts', fitWts);

        [fitParamsL, arcBoundsL, errL, used2Rev, failed2rev, hasBadBounds, fitThL, ptErrs, errArcsL, fitWtsL] = ...
            fitLogSpiral(clusMtxs(:, :, idx), ctrR, ctrC, stgs, ...
            toMin, fitParams, arcBounds);    
        allFitResL(idx) = struct('params', fitParamsL, 'bounds', arcBoundsL, ...
            'err', errL, 'errArcs', errArcsL, 'fitTh', fitThL, 'fitWts', fitWtsL);

        [fitParamsH, arcBoundsH, errH, used2Rev, failed2rev, hasBadBounds, fitThH, ptErrs, errArcsH, fitWtsH] = ...
            fitLogSpiral(clusMtxs(:, :, idx), ctrR, ctrC, stgs, ...
            toMax, fitParams, arcBounds);    
        allFitResH(idx) = struct('params', fitParamsH, 'bounds', arcBoundsH, ...
            'err', errH, 'errArcs', errArcsH, 'fitTh', fitThH, 'fitWts', fitWtsH);
    end

    function [dist, tip1, tip2] = worstItxScClusDist(idx1, idx2)
        clusMtx1 = clusMtxs(:, :, idx1);
        clusMtx2 = clusMtxs(:, :, idx2);
        
        if clusMtxIsThruCtr(clusMtx1, ctrR, ctrC, inCtrRad, outCtrRad) || ...
                clusMtxIsThruCtr(clusMtx2, ctrR, ctrC, inCtrRad, outCtrRad)
            dist = inf;
            tip1 = [];
            tip2 = [];
            return;
        end  
        
        fitRes1 = allFitRes(idx1);
        fitRes2 = allFitRes(idx2);
        if useTipWtg
            fitRes1L = allFitResL(idx1);
            fitRes1H = allFitResH(idx1);
            fitRes2L = allFitResL(idx2);
            fitRes2H = allFitResH(idx2);
        else
            fitRes1L = []; fitRes1H = []; fitRes2L = []; fitRes2H = [];
        end

        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc1, errConeArea, xFitErr, propConeFill1Near, tip1] = ...
            calcErrArcScore(clusMtx1, clusMtx2, fitRes1, fitRes1L, fitRes1H, fitRes2, ctrR, ctrC);
        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc2, errConeArea, xFitErr, propConeFill2Near, tip2] = ...
            calcErrArcScore(clusMtx2, clusMtx1, fitRes2, fitRes2L, fitRes2H, fitRes1, ctrR, ctrC);

        dist = max([1 - ptConeItxSc1, 1 - ptConeItxSc2]);
    end

    function [dist, tip1, tip2] = meanItxScClusDist(idx1, idx2)
        clusMtx1 = clusMtxs(:, :, idx1);
        clusMtx2 = clusMtxs(:, :, idx2);
        
        if clusMtxIsThruCtr(clusMtx1, ctrR, ctrC, inCtrRad, outCtrRad) || ...
                clusMtxIsThruCtr(clusMtx2, ctrR, ctrC, inCtrRad, outCtrRad)
            dist = inf;
            tip1 = [];
            tip2 = [];
            return;
        end  
        
        fitRes1 = allFitRes(idx1);
        fitRes2 = allFitRes(idx2);
        if useTipWtg
            fitRes1L = allFitResL(idx1);
            fitRes1H = allFitResH(idx1);
            fitRes2L = allFitResL(idx2);
            fitRes2H = allFitResH(idx2);
        else
            fitRes1L = []; fitRes1H = []; fitRes2L = []; fitRes2H = [];
        end

        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc1, errConeArea, xFitErr, propConeFill1Near, tip1] = ...
            calcErrArcScore(clusMtx1, clusMtx2, fitRes1, fitRes1L, fitRes1H, fitRes2, ctrR, ctrC);
        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc2, errConeArea, xFitErr, propConeFill2Near, tip2] = ...
            calcErrArcScore(clusMtx2, clusMtx1, fitRes2, fitRes2L, fitRes2H, fitRes1, ctrR, ctrC);

        dist = sqrt(prod([1 - ptConeItxSc1, 1 - ptConeItxSc2]));
    end

    function dist = coneFillThresErrRatioClusDist(idx1, idx2)
        if clusMtxIsThruCtr(clusMtxs(:, :, idx1), ctrR, ctrC, inCtrRad, outCtrRad) || ...
                clusMtxIsThruCtr(clusMtxs(:, :, idx2), ctrR, ctrC, inCtrRad, outCtrRad)
            dist = inf;
            return;
        end  
        
        [coneFillDist, tip1, tip2] = worstPropConeFillClusDist(idx1, idx2);
        
        cc = bwconncomp(clusMtxs(:, :, idx1) + clusMtxs(:, :, idx2) > 0);
        if cc.NumObjects == 1
            coneFillDist = 0;
        end
        
        if coneFillDist > maxConeFillDist
            dist = inf;
            return;
        end
        
        if useTipWtg
%             meanErrRatioNonWtd = getMeanErrRatio(idx1, idx2, false, tip1, tip2);
%             meanErrRatioWtd = getMeanErrRatio(idx1, idx2, true, tip1, tip2);
%             dist = min(meanErrRatioNonWtd, meanErrRatioWtd);
            dist = getMeanErrRatio(idx1, idx2, true, tip1, tip2);
        else
            dist = getMeanErrRatio(idx1, idx2, false, tip1, tip2);
        end
 
        if cc.NumObjects == 1 && (dist <= stopThres)
            dist = (stopThres/2) * dist;
        else
            dist = (stopThres/2) * dist + (stopThres/2);
        end
    end

    function dist = errRatioThresConeFillClusDist(idx1, idx2)
        [baseDist, tip1, tip2] = baseClusDistFxn(idx1, idx2);
        
        if baseDist > stopThres
            dist = inf; 
            return;
        end
        
%         if useTipWtg
%             mtx1 = zeros(imgSz);
%             mtx2 = zeros(imgSz);
%             if tip1 < 0
%                 mtx1(clusMtxs(:, :, idx1) > 0) = allFitResL(idx1).fitWts;
%                 fitRes1 = allFitResL(idx1);
%             elseif tip1 > 0
%                 mtx1(clusMtxs(:, :, idx1) > 0) = allFitResH(idx1).fitWts;
%                 fitRes1 = allFitResH(idx1);
%             else
%                 assert(false);
%             end
%             if tip2 < 0
%                 mtx2(clusMtxs(:, :, idx2) > 0) = allFitResL(idx2).fitWts;
%                 fitRes2 = allFitResL(idx2);
%             elseif tip2 > 0
%                 mtx2(clusMtxs(:, :, idx2) > 0) = allFitResH(idx2).fitWts;
%                 fitRes2 = allFitResH(idx2);
%             else
%                 assert(false);
%             end
%         else
%             mtx1 = clusMtxs(:, :, idx1);
%             mtx2 = clusMtxs(:, :, idx2);
%             fitRes1 = allFitRes(idx1);
%             fitRes2 = allFitRes(idx2);
%         end
% %         mtx1 = clusMtxs(:, :, idx1);
% %         mtx2 = clusMtxs(:, :, idx2);
 
        if useTipWtg
%             meanErrRatioNonWtd = getMeanErrRatio(idx1, idx2, false, tip1, tip2);
%             meanErrRatioWtd = getMeanErrRatio(idx1, idx2, true, tip1, tip2);
%             meanErrRatio = min(meanErrRatioNonWtd, meanErrRatioWtd);
            meanErrRatio = getMeanErrRatio(idx1, idx2, true, tip1, tip2);
        else
            meanErrRatio = getMeanErrRatio(idx1, idx2, false, tip1, tip2);
        end
        
        meanErrRatio
        if meanErrRatio > meanErrRatioThres
            dist = inf;
        else
            dist = baseDist;
        end
    end

    function mergeErrRatio = getMeanErrRatio(idx1, idx2, useTipWtg, tip1, tip2)
        if useTipWtg
            thVals1 = allFitRes(idx1).fitTh;
            thVals2 = allFitRes(idx2).fitTh;
            if tip1 < 0
                fitRes1 = allFitResL(idx1);
            elseif tip1 > 0
                fitRes1 = allFitResH(idx1);
            else % tip-weighting overridden 
%                 assert(false);
                fitRes1 = allFitRes(idx1);
                thVals1 = [];
                useTipWtg = false;
            end
            if tip2 < 0
                fitRes2 = allFitResL(idx2);
            elseif tip2 > 0
                fitRes2 = allFitResH(idx2);
            else
%                 assert(false);
                fitRes2 = allFitRes(idx2);
                thVals2 = [];
                useTipWtg = false;
            end
        else
            fitRes1 = allFitRes(idx1);
            fitRes2 = allFitRes(idx2);
            thVals1 = [];
            thVals2 = [];
        end
        mtx1 = clusMtxs(:, :, idx1);
        mtx2 = clusMtxs(:, :, idx2);
        
%         idx2
%         allFitResH(idx2).err
%         fprintf('stored errors: \n');
%         fitRes1.err
%         fitRes2.err
        mergeErrRatio = ...
            calcArcMergeErr(mtx1, ...
            fitRes1.params, fitRes1.bounds, fitRes1.err, ...
            mtx2, ...
            fitRes2.params, fitRes2.bounds, fitRes2.err, ...
            ctrR, ctrC, stgs, [], useTipWtg, thVals1, thVals2, plotFlag);
    end

    function dist = errRatioAndConeFillClusDist(idx1, idx2)        
        mergeErrRatio = ...
            calcArcMergeErr(clusMtxs(:, :, idx1), ...
            lgspParams(idx1, :), lgspBounds(idx1, :), fitErrs(idx1, :), ...
            clusMtxs(:, :, idx2), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), fitErrs(idx2, :), ...
            ctrR, ctrC, [], stgs);
        
        dist = mergeErrRatio * worstPropConeFillClusDist(idx1, idx2);
    end

    function [dist, tip1, tip2] = worstPropConeFillClusDist(idx1, idx2)
        clusMtx1 = clusMtxs(:, :, idx1);
        clusMtx2 = clusMtxs(:, :, idx2);
        
        if clusMtxIsThruCtr(clusMtx1, ctrR, ctrC, inCtrRad, outCtrRad) || ...
                clusMtxIsThruCtr(clusMtx2, ctrR, ctrC, inCtrRad, outCtrRad)
            dist = inf;
            tip1 = [];
            tip2 = [];
            return;
        end  
        
        fitRes1 = allFitRes(idx1);
        fitRes2 = allFitRes(idx2);
        if useTipWtg
            fitRes1L = allFitResL(idx1);
            fitRes1H = allFitResH(idx1);
            fitRes2L = allFitResL(idx2);
            fitRes2H = allFitResH(idx2);
        else
            fitRes1L = []; fitRes1H = []; fitRes2L = []; fitRes2H = [];
        end

        [bnds2Trans, rotAdj, clusThGap] = translateThVals(fitRes1.params, fitRes1.bounds, fitRes2.params, fitRes2.bounds, fitRes2.fitTh);
        if clusThGap > maxClusThGap
%             fprintf('rejecting merge, theta-gap too large (was %2.4f, max %2.4f)\n',...
%                 clusThGap, maxClusThGap);
            dist = inf;
            tip1 = [];
            tip2 = [];
            return
        end
        
%         thOverlap = abs(diff(intervalIntersect(fitRes1.bounds, bnds2Trans)));
        thOverlap = intervalIntersect(...
            repmat(fitRes1.bounds, 3, 1) + repmat([-2*pi;0;2*pi], 1, 2), ...
            repmat(bnds2Trans, 3, 1));
        thOverlap = max(abs(diff(thOverlap, [], 2)));
        if ~isnan(thOverlap) && (thOverlap > maxClusThOverlap)
            dist = inf;
            tip1 = [];
            tip2 = [];
            return
        end
        
        [propInCone, wtdPropInCone, propConeFill1, ptConeItxSc, errConeArea, xFitErr, propConeFill1Near, tip1] = ...
            calcErrArcScore(clusMtx1, clusMtx2, fitRes1, fitRes1L, fitRes1H, fitRes2, ctrR, ctrC);
        [propInCone, wtdPropInCone, propConeFill2, ptConeItxSc, errConeArea, xFitErr, propConeFill2Near, tip2] = ...
            calcErrArcScore(clusMtx2, clusMtx1, fitRes2, fitRes2L, fitRes2H, fitRes1, ctrR, ctrC);

%         dist = max([1 - propConeFill1, 1 - propConeFill2]);
        dist = max([1 - propConeFill1Near, 1 - propConeFill2Near]);
    end

    function dist = meanCrossFitClusDist(idx1, idx2)
        clusMtx1 = clusMtxs(:, :, idx1);
        clusMtx2 = clusMtxs(:, :, idx2);
        fitRes1 = allFitRes(idx1);
        fitRes2 = allFitRes(idx2);

        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc, errConeArea, xFitErr1] = ...
            calcErrArcScore(clusMtx1, clusMtx2, fitRes1, [], [], fitRes2, ctrR, ctrC);
        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc, errConeArea, xFitErr2] = ...
            calcErrArcScore(clusMtx2, clusMtx1, fitRes2, [], [], fitRes1, ctrR, ctrC);

        dist = mean([xFitErr1, xFitErr2]);
    end

    function dist = maxCrossFitClusDist(idx1, idx2)
        clusMtx1 = clusMtxs(:, :, idx1);
        clusMtx2 = clusMtxs(:, :, idx2);
        fitRes1 = allFitRes(idx1);
        fitRes2 = allFitRes(idx2);

        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc, errConeArea, xFitErr1] = ...
            calcErrArcScore(clusMtx1, clusMtx2, fitRes1, [], [], fitRes2, ctrR, ctrC);
        [propInCone, wtdPropInCone, propConeFill, ptConeItxSc, errConeArea, xFitErr2] = ...
            calcErrArcScore(clusMtx2, clusMtx1, fitRes2, [], [], fitRes1, ctrR, ctrC);

        dist = max([xFitErr1, xFitErr2]);

    end

    % distance by meanErrRatio if within pixel and arc distances, inf
    % otherwise
    function dist = calcClusDistHybrid(idx1, idx2)
%         calcClosestArcTips(lgspParams(idx1, :), lgspBounds(idx1, :), ...
%             lgspParams(idx2, :), lgspBounds(idx2, :));
        
        maxPixDist = mean([size(clusMtxs, 1), size(clusMtxs, 2)]) / 20; 
        maxArcDist = mean([size(clusMtxs, 1), size(clusMtxs, 2)]) / 10; 
        
        pixDists = bwdist(clusMtxs(:, :, idx1) > 0);
        pixDists = pixDists(clusMtxs(:, :, idx2) > 0);
        
        arcDist = calcArcDist(lgspParams(idx1, :), lgspBounds(idx1, :), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), imgSz, ctrR, ctrC);
        
        if isempty(pixDists) || min(pixDists) > maxPixDist || arcDist > maxArcDist
            dist = inf;
            return;
        end
        
        mergeErrRatio = ...
            calcArcMergeErr(clusMtxs(:, :, idx1), ...
            lgspParams(idx1, :), lgspBounds(idx1, :), fitErrs(idx1, :), ...
            clusMtxs(:, :, idx2), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), fitErrs(idx2, :), ...
            ctrR, ctrC, [], stgs);
        
        dist = mergeErrRatio;
%         if arcDist <= maxArcDist
%             dist = meanErrRatio;
%         else
%             dist = inf;
%         end
    end

%     % distance by meanErrRatio if within pixel and arc distances, inf
%     % otherwise
%     function dist = calcClusDistHybrid(clusMtx1, clusMtx2, params1, params2,...
%             bounds1, bounds2, errs1, errs2, allowArcBeyond2pi, ctrR, ctrC, imgSz)
%         
%         maxPixDist = mean([size(clusMtx1, 1), size(clusMtx1, 2)]) / 20; 
%         maxArcDist = mean([size(clusMtx1, 1), size(clusMtx1, 2)]) / 10; 
%         
%         
%         pixDists = bwdist(clusMtx1 > 0);
%         pixDists = pixDists(clusMtx2 > 0);
%         
%         arcDist = calcArcDist(params1, bounds1, ...
%             params2, bounds2, imgSz, ctrR, ctrC);
%         
%         if isempty(pixDists) || min(pixDists) > maxPixDist || arcDist > maxArcDist
%             dist = inf;
%             return;
%         end
%         
%         [errRatio, meanErrRatio] = ...
%             calcArcMergeErr(clusMtx1, params1, bounds1, errs1, ...
%             clusMtx2, params2, bounds2, errs2,...
%             ctrR, ctrC, allowArcBeyond2pi);
%         
%         dist = meanErrRatio;
% %         if arcDist <= maxArcDist
% %             dist = meanErrRatio;
% %         else
% %             dist = inf;
% %         end
%     end

    % distance by meanErrRatio if within pixel distance, inf otherwise
    function dist = calcClusDistByArcFit(idx1, idx2)
        dists = bwdist(clusMtxs(:, :, idx1) > 0);
        dists = dists(clusMtxs(:, :, idx2) > 0);
        
        if ~isempty(dists) && min(dists) <= maxPixDist
            mergeErrRatio = ...
                calcArcMergeErr(clusMtxs(:, :, idx1), [], [], [], ...
                clusMtxs(:, :, idx2), [], [], [], ctrR, ctrC,...
                [], stgs);
            dist = mergeErrRatio;
        else
            dist = inf;
        end
    end

    function dist = calcClusDistByArcExt(idx1, idx2)
%         [errRatio, meanErrRatio] = ...
%             calcArcMergeErr(clusMtxs(:, :, idx1), [], [], [], ...
%             clusMtxs(:, :, idx2), [], [], [], ctrR, ctrC,...
%             allowArcBeyond2pi);
        arcDist = calcArcDist(lgspParams(idx1, :), lgspBounds(idx1, :), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), imgSz, ctrR, ctrC);
        dist = arcDist;
        if meanErrRatio > 1
            dist = inf;
        end
    end

end