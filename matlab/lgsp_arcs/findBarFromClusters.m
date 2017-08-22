function [barInds, barAngles, barHalfLens, otherClusMtxs] = ...
    findBarFromClusters(clusMtxs, ctrR, ctrC)
% Obsolete, should NOT be used!
% Looks for bars in the given set of cluster matrices, using some
%   bar-detection heuristics
% INPUTS:
%   clusMtxs: a 3D array where nonzero elements of a page are members of a
%       particular cluster
%   ctrR: the row index of the center of the log-spirals
%   ctrC: the column index of the center of the log-spirals
% OUTPUTS:
%   barInds: the indices of cluster matrices where bars were detected
%   barAngles: the angles, in radians, of any bars indexed by barInds
%   barHalfLens: the half-lengths of any bars indexed by barInds
%   otherClusMtxs: the cluster matrices where bars where not found

error('this function is obsolete and should not be used.');

nClus = size(clusMtxs, 3);

barDetections = false(nClus, 1);
barAngles = zeros(nClus, 1);
barHalfLens = zeros(nClus, 1);

for clus = 1:1:nClus
    [curAngle, curLen] = ...
        fitBarUsingGaussian(clusMtxs(:, :, clus), ctrR, ctrC, false);
    if ~isempty(curAngle)
        barAngles(clus) = curAngle;
        barHalfLens(clus) = curLen;
        barDetections(clus) = true;
    end
end

barInds = find(barDetections);
barAngles = barAngles(barDetections);
barHalfLens = barHalfLens(barDetections);
otherClusMtxs = clusMtxs(:, :, ~barDetections);

if sum(barDetections) > 1
    warning('more than one bar was found')
end

end