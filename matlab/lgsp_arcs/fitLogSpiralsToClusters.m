function [params, thBounds, sumSqErr, used2rev, failed2rev, hasBadBounds, thVals] = ...
    fitLogSpiralsToClusters(clusMtxs, ctrR, ctrC, stgs)
% Fits a logarithmic spiral to each given cluster matrix
% INPUTS: 
%   clusMtxs: cluster matrices produced by hac2mtxs()
%   ctrR: the row position of the center of the log-spirals
%   ctrC: the column position of the center of the log-spirals
%   stgs: structure containing algorithm settings (see settings.m)
% OUTPUTS:
%   params: the parameters of the logarithmic spirals, one per row
%   thBounds: starting and ending theta-range of points to which
%       log-spirals were fitted, counterclockwise from the theta-offset
%       parameter.  Each row is a [start, end] pair.
%   sumSqErr: sum of squared errors for each fitted log-spiral

startParams = [0, 0, 10];

nClus = size(clusMtxs, 3);

params = zeros(nClus, 3);
thBounds = zeros(nClus, 2);
sumSqErr = zeros(nClus, 1);
used2rev = false(nClus, 1);
failed2rev = false(nClus, 1);
hasBadBounds = false(nClus, 1);
thVals = cell(nClus, 1);
for clus=1:1:nClus
    [params(clus, :), thBounds(clus, :), sumSqErr(clus), ...
            used2rev(clus), failed2rev(clus), hasBadBounds(clus), ...
            thVals{clus}] = ...
        fitLogSpiral(clusMtxs(:, :, clus), ctrR, ctrC, ...
            stgs, [], startParams);
end

end