function [clusMtxs, nMerges] = mergeClustersByFit(clusMtxs, ctrR, ctrC,...
    barInfo, stgs)
% Merges clusters that are sufficiently close and retain a good log-spiral
%   fit. Unlike the main HAC step, this uses arc fits but not pixel 
%   similarities.
% INPUTS:
%   clusMtxs: a 3D array where nonzero elements of a page are members of a
%       particular cluster
%   ctrR: the row index of the center of the log-spirals
%   ctrC: the column index of the center of the log-spirals
%   barInfo: 
%   stgs: 

if size(clusMtxs, 3) == 0
    nMerges = 0;
    return;
end

tMerge = tic;

stopThres = stgs.errRatioThres;

imgSz = size(clusMtxs); imgSz = imgSz(1:2);
nClus = size(clusMtxs, 3);
% maxPixDist = 1;
maxPixDist = mean([size(clusMtxs, 1), size(clusMtxs, 2)]) / 20;
% maxArcDist = 10;
maxArcDist = mean([size(clusMtxs, 1), size(clusMtxs, 2)]) / 10;
% stopThres = 10;
% calcClusDist = @calcClusDistByArcExt;
% calcClusDist = @calcClusDistHybrid; 

[lgspParams, lgspBounds, fitErrs, used2rev, failed2rev, hasBadBounds, thVals] = ...
    fitLogSpiralsToClusters(clusMtxs, ctrR, ctrC, stgs);

clusDists = inf * ones(nClus, nClus);

dists = zeros(size(clusMtxs(:, :, 1)));
for ii=1:1:nClus
    for jj=ii+1:1:nClus
        clusDists(ii, jj) = calcClusDist(ii, jj);
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
%         title(sprintf('val = %2.4f, sizes = %2.4f, %2.4f', ...
%             val, nnz(clusMtxs(:, :, cIdx1)), nnz(clusMtxs(:, :, cIdx2))));
%         subplot(1, 3, 2); imshow(clusMtxs(:, :, cIdx2));
%         subplot(1, 3, 3); imshow(clusMtxs(:, :, cIdx1) + clusMtxs(:, :, cIdx2));
        
        clusMtxs(:, :, cIdx1) = clusMtxs(:, :, cIdx1) + clusMtxs(:, :, cIdx2);
        clusMtxs(:, :, cIdx2) = 0;
        deletedClusts(cIdx2) = true;
        clusDists(cIdx2, :) = inf;
        clusDists(:, cIdx2) = inf;
        nMerges = nMerges + 1;
        
        [paramsNew, boundsNew, errNew, tmp1, tmp2, tmp3, thNew] = ...
            fitLogSpiral(clusMtxs(:, :, cIdx1), ctrR, ctrC, stgs);
        lgspParams(cIdx1, :) = paramsNew;
        lgspBounds(cIdx1, :) = boundsNew;
        fitErrs(cIdx1, :) = errNew;
        thVals{cIdx1} = thNew;
        
        for idx=1:1:cIdx1-1
            if ~deletedClusts(idx)
                clusDists(idx, cIdx1) = calcClusDist(idx, cIdx1);
            end
        end
        for idx=cIdx1+1:1:nClus
            if ~deletedClusts(idx)
                clusDists(cIdx1, idx) = calcClusDist(cIdx1, idx);
            end
        end
    else
        keepMerging = false;
    end
end
clusMtxs = clusMtxs(:, :, ~deletedClusts);

fprintf('Time for arc-merging: \n');
toc(tMerge);

    % distance by meanErrRatio if within arc distance, inf otherwise
    function dist = calcClusDist(idx1, idx2)
        pixDists = bwdist(clusMtxs(:, :, idx1) > 0);
        pixDists = pixDists(clusMtxs(:, :, idx2) > 0);
        
        if isempty(pixDists) || min(pixDists) > maxPixDist
            dist = inf;
            return;
        end
        
        mergeErrRatio = ...
            calcArcMergeErr(clusMtxs(:, :, idx1), ...
            lgspParams(idx1, :), lgspBounds(idx1, :), fitErrs(idx1, :), ...
            clusMtxs(:, :, idx2), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), fitErrs(idx2, :), ...
            ctrR, ctrC, barInfo, stgs, [], thVals{idx1}, thVals{idx2});
        
        dist = mergeErrRatio;
    end

    % distance by meanErrRatio if within pixel and arc distances, inf
    % otherwise
    function dist = calcClusDistHybrid(idx1, idx2)
        pixDists = bwdist(clusMtxs(:, :, idx1) > 0);
        pixDists = pixDists(clusMtxs(:, :, idx2) > 0);
        
        arcDist = calcArcDist(lgspParams(idx1, :), lgspBounds(idx1, :), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), imgSz, ctrR, ctrC, stgs);
        
        if isempty(pixDists) || min(pixDists) > maxPixDist || arcDist > maxArcDist
            dist = inf;
            return;
        end
        
        mergeErrRatio = ...
            calcArcMergeErr(clusMtxs(:, :, idx1), ...
            lgspParams(idx1, :), lgspBounds(idx1, :), fitErrs(idx1, :), ...
            clusMtxs(:, :, idx2), ...
            lgspParams(idx2, :), lgspBounds(idx2, :), fitErrs(idx2, :), ...
            ctrR, ctrC, barInfo, stgs, [], thVals{idx1}, thVals{idx2});
        
        dist = mergeErrRatio;
%         if arcDist <= maxArcDist
%             dist = meanErrRatio;
%         else
%             dist = inf;
%         end
    end

    % distance by meanErrRatio if within pixel distance, inf otherwise
    function dist = calcClusDistByArcFit(idx1, idx2)
        dists = bwdist(clusMtxs(:, :, idx1) > 0);
        dists = dists(clusMtxs(:, :, idx2) > 0);
        
        if ~isempty(dists) && min(dists) <= maxPixDist
            mergeErrRatio = ...
                calcArcMergeErr(clusMtxs(:, :, idx1), [], [], [], ...
                clusMtxs(:, :, idx2), [], [], [], ctrR, ctrC,...
                barInfo, stgs);
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
            lgspParams(idx2, :), lgspBounds(idx2, :), imgSz, ctrR, ctrC, stgs);
        dist = arcDist;
        if meanErrRatio > errRatioThres
            dist = inf;
        end
    end

end