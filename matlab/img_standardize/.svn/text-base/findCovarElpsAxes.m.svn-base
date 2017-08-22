function [majAxsVec, majAxsAngle, axisRatio, majAxsLength] = ...
    findCovarElpsAxes(covar, likCutoff, imgSz)
% Determines major and minor axis information from a given covariance
%   matrix
% INPUTS:
%   covar: the covariance matrix
%   likCutoff: the likelihood value to use as the ellipse contour; needed
%       if majAxsLength is specified
% OUTPUTS:
%   majAxsVec: unit vector pointing in the direction of the major axis
%   majAxsAngle: angle, in radians counterclockwise from the positive 
%       x-axis, specifying the direction of the major axis
%   axisRatio: number in [0, 1] specifying the length ratio of the minor 
%       to the major axis
%   majAxsLength: length of the major axis, using the elliptical contour
%       specified by the likelihood value likCutoff

if nargin < 1
    error('no covariance matrix given');
end

if (nargin < 2 || isempty(likCutoff)) && nargout >= 4
    error('need to specify likCutoff if majAxsLength needed')
end

[eigVec, eigVal] = eig(covar);
eigVal = diag(eigVal);
[maxEVal, mevIdx] = max(eigVal);
majAxsVec = eigVec(:, mevIdx);
% minAxsVec = eigVec(:, 3 - mevIdx) % 3-x converts from 2<->1

majAxsAngle = atan(majAxsVec(2)/majAxsVec(1));

% majAxsVar = majAxsVec' * covarFit * majAxsVec

axisRatio = sqrt(min(eigVal)) / sqrt(maxEVal);

if nargout >= 4
    startPt = [eps 5*max(imgSz)];
    majAxsLength = abs(fzero(@(x)(mvnpdf([majAxsVec(1) * x, majAxsVec(2) * x],...
        [0, 0], covar) - likCutoff), startPt));
    majAxsLength = 2 * majAxsLength;
end

end