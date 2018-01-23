function [clusImg, clusColors] = showClustersFromMtxs(mtxs, imgSz, infoStr, colorThres, showBoundaries)
% Show clusters from cluster matrices instead of original structure array

if nargin < 2 || isempty(imgSz)
    imgSz = size(mtxs);
    imgSz = imgSz(1:2);
end

if nargin < 3
    infoStr = [];
end

if nargin < 4 || isempty(colorThres)
    colorThres = 0;
end

if nargin < 5 || isempty(showBoundaries)
    showBoundaries = false;
end

nClus = size(mtxs, 3);
clusts = repmat(struct('mergedPts', []), 1, nClus);
for clus=1:1:nClus
    clusts(clus).mergedPts = find(mtxs(:, :, clus));
end

if nargout < 1
    showClusters(clusts, imgSz, colorThres, infoStr);
else
    [clusImg, clusColors] = showClusters(clusts, imgSz, colorThres, infoStr);
end

if showBoundaries
    bdryBrt = 0.5;
    bdryImg = displayClusterOverlay(zeros(imgSz), mtxs);
    bdryImg = repmat(sum(bdryImg, 3) > 0, [1 1 3]);
    clusImg(bdryImg) = bdryBrt;
end
end