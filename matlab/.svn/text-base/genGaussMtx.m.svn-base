function gaussMtx = genGaussMtx(sdev, msize)
% Generates a matrix of values of a Gaussian, with mean at the center
% INPUTS:
%   sdev: standard deviation of the Gaussian, in pixels
%   msize: optional parameter for the number of rows and columns of the
%       matrix
% OUTPUTS:
%   gaussMtx: the matrix of Gaussian values

% default size of the matrix: 
% 2 * radius of 2 standard deviations, plus center
if nargin < 2 || isempty(msize)
    msize = 2 * 2 * ceil(sdev) + 1;
end

% generate the Gaussian matrix
gaussMtx = zeros(msize, msize);
mu = floor(msize / 2) * [1;1] + [1;1];
sigma = sdev * eye(2);
for ii=1:1:msize
    for jj=1:1:msize
        gaussMtx(ii, jj) = mvnpdf([ii; jj], mu, sigma);
    end
end

end