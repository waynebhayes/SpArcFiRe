function [mu, covar] = calcWtdMeanCovar(img, muFix)
% Fits a Gaussian to the image pixels using the closed-form expressions.  
% Intensities are considered weights on data points in Cartesian x-y 
% coordinate space, with the origin at the bottom left.
% INPUTS:
%   img: image for which the weighted mean and covariance are to be
%       computed
%   muFix: optional; if set, fixes the mean to the given value and only
%       calculates the covariance from the data
% OUTPUTS:
%   mu: the mean of the data
%   covar: the covariance of the data (under the assumption of Gaussianity)

if nargin < 2
    muFix = [];
end

img = flipud(img); % y-values are reversed vs row values

[xVals, yVals] = meshgrid(1:1:size(img, 2), 1:1:size(img, 1));
sumWts = sum(img(:));
if isempty(muFix)
    xMean = sum(sum(img .* xVals)) / sumWts;
    yMean = sum(sum(img .* yVals)) / sumWts;
    mu = [xMean yMean];
else
    mu = muFix;
    xMean = mu(1);
    yMean = mu(2);
end

xDev = xVals - xMean;
yDev = yVals - yMean;
valsCov11 = xDev .^ 2;
valsCov22 = yDev .^ 2;
valsCov12 = xDev .* yDev;
cov11 = sum(sum(img .* valsCov11)) / sumWts;
cov22 = sum(sum(img .* valsCov22)) / sumWts;
cov12 = sum(sum(img .* valsCov12)) / sumWts;

covar = [cov11, cov12; cov12, cov22];

end