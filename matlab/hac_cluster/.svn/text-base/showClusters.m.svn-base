function [clusImg, clusterColors] = showClusters(clusters, imgSz, colorThres, infoStr)
% Produces an image showing cluster memberships of image pixels
% INPUTS:
%   clusters: HAC dendogram produced by genHacTree
%   imgSz: size of the image from which the HAC clusters were derived
%   colorThres: minimum size for clusters to get their own color
%   infoStr: plot title, if the image is to be displayed to the screen
% OUTPUTS:
%   clusImg: RGB image of cluster memberships; if not requested, cluster
%       memberships will instead be displayed to the screen in a new plot
%       window

if nargin < 4 || isempty(infoStr)
    infoStr = '';
end

% colorThres = 150;

nClusts = length(clusters);
% fprintf('%d clusters\n', nClusts);
clusterSizes = getClusterSizes(clusters);
clusterHasColor = clusterSizes > colorThres;
% colorClusts = clusters(getClusterSizes(clusters) > 1);
nColorClusts = sum(clusterHasColor);
% fprintf('%d colored clusters\n', nColorClusts);
% clusterSizes(clusterSizes > colorThres)
clusterColors = 0.2 * ones(nClusts, 3);
tempFig = figure; % suppress figure-generation by colormap function
clusterColors(clusterHasColor, :) = colormap(hsv(nColorClusts));
close(tempFig);

% clusImg = zeros([imgSz 3]);
r = zeros(imgSz);
g = zeros(imgSz);
b = zeros(imgSz);
% clusBoundaries = zeros(imgSz);
for clus=1:1:nClusts
%     clusterColors(clus, :)
    clusPts = clusters(clus).mergedPts;
%     [r, c] = ind2sub(imgSz, clusters(clus).mergedPts);
%     clusImg(r, c, :)
%     clusterColors(clus, :)
    r(clusPts) = clusterColors(clus, 1);
    g(clusPts) = clusterColors(clus, 2);
    b(clusPts) = clusterColors(clus, 3);
    
%     clusImg = zeros(imgSz); clusImg(clusPts) = 1;
%     clusBoundaries = clusBoundaries | ...
%         (clusImg & ~imerode(clusImg, strel('square', 3)));
end

clusImg = zeros([imgSz 3]);
clusImg(:, :, 1) = r;
clusImg(:, :, 2) = g;
clusImg(:, :, 3) = b;
% clusImg(repmat(clusBoundaries, [1, 1, 3]) > 0) = 0.5;

if nargout < 1
%     figure; imshow(clusImg); title(infoStr)
    figure; imagesc(clusImg);  axis image; axis off % TEMP - plot for paper
%     figure; hist(getClusterSizes(clusters)); title(infoStr)
end

end