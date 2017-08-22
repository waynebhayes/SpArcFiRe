function [mtxs] = hac2mtxs(clusters, img, clusSizeCutoff)
% Converts clustering results into a set of matrices, each containing the
% members of a single cluster.
% INPUTS:
%   clusters: clustering result from the HAC clusterer
%   img: the image from which the clusters were generated
%   clusSizeCutoff: minimum size needed for there to be an output matrix 
%       for the cluster (otherwise the cluster is discarded)
% OUTPUTS:
%   mtxs: 3D array where each page is a matrix where nonzero elements are
%       members of a single cluster.  Matrix positions correspond to pixel
%       positions in the image from which the clusters were generated.
%       Intensities of nonzero elements are the image intensities.

error(nargchk(nargin, 3, 3))

clusterSizes = getClusterSizes(clusters);
meetsCutoff = clusterSizes >= clusSizeCutoff;
clusters = clusters(meetsCutoff);
clusterSizes = clusterSizes(meetsCutoff);
[clusterSizes, order] = sort(clusterSizes, 'descend');
clusters = clusters(order);

nClusts = length(clusterSizes);
imgSz = size(img);
mtxs = zeros([imgSz nClusts]);

for clus=1:1:nClusts
    clusPts = clusters(clus).mergedPts;
    inClus = false(imgSz);
    inClus(clusPts) = true;
    mtxs(:, :, clus) = inClus .* img;
end

end