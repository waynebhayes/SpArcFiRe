function [clusMtx, convPts, clusFrom, clusTo] = ...
    meanShiftCluster(mtx, bandw, wrapCols)
% Performs mean-shift clustering on the input matrix, as follows: for each
%   point in the matrix, we trace a trajectory where we move to the mean of
%   the neighboring points, weighted by the values at those points.  When
%   we reach a rest position, the value of the rest position's index in the
%   output matrix is incremented by the value of the start position in the
%   input matrix.  Thus, each point votes for the mode it reaches during
%   mean-shift iteration, with vote weights determined by the value of the
%   point.
% INPUTS:
%   mtx: matrix for which we are performing mean-shift clustering
%   bandw: radius of points to consider when performing the next mean-shift
%     step.  Larger bandwidths tend to put more points in the same cluster.
%   wrapCols: whether to consider the first and last columns to be adjacent
%     to each other (default false)
% OUTPUTS:
%   clusMtx: mean-shift cluster result.  Points are nonzero if they were a
%     mode in the input matrix, with values equal to the sum of values of
%     points that were attracted to the mode.
%   convPts, clusFrom, clusTo...

% inclThres = quantile(mtx(:), .1);
inclThres = 0;
maxMoves = 100;
moveMinDist = 10^-1;
maxMoveCutoffs = 0;

% roundMeans = false;

[allR, allC] = find(mtx > inclThres);
allPts = [allR, allC];
numPts = length(allPts);
convPts = zeros(numPts, 2);

clusMtx = zeros(size(mtx));
clusFrom = zeros(size(mtx));

for sPt = 1:1:numPts
    prevMean = [Inf, Inf];
    curMean = allPts(sPt, :);
    numMoves = 0;
    while sqrt(sum((prevMean - curMean) .^ 2)) > moveMinDist
        if numMoves > maxMoves
            maxMoveCutoffs = maxMoveCutoffs + 1;
            break;
        end
        
        prevMean = curMean;
        curMean = msStep(mtx, bandw, curMean, wrapCols);
%         if roundMeans
%             curMean = round(curMean);
%         end
        numMoves = numMoves + 1;
    end
    % cluster matrix element at stationary point gets the votes of the
    % original matrix at the position of the starting point
    clusMtx(round(curMean(1)), round(curMean(2))) ...
        = clusMtx(round(curMean(1)), round(curMean(2))) + ...
        mtx(allPts(sPt, 1), allPts(sPt, 2));
    
    convPts(sPt, :) = curMean;
    clusFrom(allPts(sPt, 1), allPts(sPt, 2)) = sub2ind(size(mtx), round(curMean(1)), round(curMean(2)));
end

clusTo = zeros(size(mtx));
clusIdxs = find(clusMtx > 0);
numClus = size(clusIdxs, 1);
for ii=1:1:numClus
    clusTo(clusIdxs(ii)) = ii;
    clusFrom(clusFrom == clusIdxs(ii)) = ii;
end
clusFrom = clusFrom / numClus;
clusTo = clusTo / numClus;

maxMoveCutoffs

end