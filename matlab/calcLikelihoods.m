function lik = calcLikelihoods(imgSz, mean, covar)
% Calculates 2D Gaussian likelihoods at each pixels' point in Cartesian
% space
% INPUTS:
%   imgSz: size of the original image
%   mean: mean of the Gaussian
%   covar: covariance matrix of the Gaussian
% OUTPUTS:
%   lik: image matrix of likelihoods, where each pixel corresponds to the
%       likelihood of the corresponding Cartesian point, where the origin
%       is at the bottom left.

[xVals, yVals] = meshgrid(1:imgSz(2), 1:imgSz(1));

lik = mvnpdf([xVals(:) yVals(:)], mean, covar);

lik = reshape(lik, imgSz);

lik = flipud(lik); % y-values are reversed vs row values

end