function [clusMtxs, nMerges] = mergeClustersByGuide(clusMtxs, outputPath, threshold, guideImageFile)
% Merges clusters that are sufficiently close and retain a good log-spiral
%   fit. Unlike the main HAC step, this uses arc fits but not pixel 
%   similarities.
% INPUTS:
%   clusMtxs: a 3D array where nonzero elements of a page are members of a
%       particular cluster
%   outputPath: filename where SpArcFiRe galaxy info should be
%   threshold: a scalar that indicates to merge clusters after their overlap
%       passes the threshold
%   guideImageFile: path to image to be used as guide

fprintf('Using %s as guide image.\n',guideImageFile);
guideImage = imread(guideImageFile);
if size(clusMtxs, 3) == 0
    nMerges = 0;
    return;
end

tMerge = tic;

% Get guide clusters
%color = outputPath(end-1:end);
%meanPicPath = strcat(outputPath, '-H_clusMask-merged.png');
%meanPicPath = regexprep(meanPicPath, color, '_mean');
%guideClusMask = imread(meanPicPath);
%guideClusMask_columns = reshape(guideClusMask, [], 3);
guideClusMask_columns = reshape(guideImage, [], 3);
[unique_guide_clusters, m, n] = unique(guideClusMask_columns, 'rows');
guideClusters = zeros(256, 256, size(unique_guide_clusters,1)-1);
for g = 2:size(unique_guide_clusters,1)
    guide_cluster = n == g;
    guide_cluster = reshape(guide_cluster, 256, 256);
    guideClusters(:,:,g-1) = guide_cluster;
end

imgSz = size(clusMtxs); imgSz = imgSz(1:2);
nClus = size(clusMtxs, 3);

guideImgSz = size(guideClusters);
nGuideClus = size(guideClusters, 3);
guideImgSz = guideImgSz(1:2);

blob_list = zeros(guideImgSz(1), guideImgSz(2), nGuideClus);
deletedClusts = false(1, nGuideClus);
guideBinary = logical(guideClusters);
imgBinary = logical(clusMtxs);

nMerges = 0;
for i = 1:nClus
    maxClus = -1;
    maxIntersect = threshold;
    for r = 1:nGuideClus
        x = guideBinary(:,:,r) & imgBinary(:,:,i);
        total_intersect = nnz(x);
        total_img = nnz( imgBinary(:,:,i) );
        fracIntersect = double(total_intersect) / double(total_img);
        if fracIntersect > maxIntersect
            maxClus = r;
            maxIntersect = fracIntersect;
        end
    end
    if maxClus > 0
        blob_list(:,:,maxClus) = blob_list(:,:,maxClus) + clusMtxs(:,:,i);
        guideBinary(:,:,maxClus) = guideBinary(:,:,maxClus) + clusMtxs(:,:,i);
        nMerges = nMerges + 1;
    end
end

for i = 1:nGuideClus
    if nnz( blob_list(:,:,i) ) == 0
        deletedClusts(i) = true;
    end
end
clusMtxs = blob_list(:,:, ~deletedClusts);
fprintf("Time for arc-merging: \n");
toc(tMerge);

end
