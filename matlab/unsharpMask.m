function [img, adjustment] = unsharpMask(img, stgs)
% Performs an unsharp mask operation on the given image
% INPUTS: 
%   img: the image for which the unsharp mask is to be performed
%   stgs:
% OUTPUTS:
%   img: the result of the unsharp mask operation

amount = stgs.unsharpMaskAmt;
sigma = stgs.unsharpMaskSigma;
% temp_useWideUsm = stgs.temp_useWideUsm;

% if nargin < 3 || isempty(sigma)
%     sigma = 25;
% end

maxAdj = inf;
% maxAdj = 1;

% TODO: use a larger filter radius
% if temp_useWideUsm
%     filter_size = 2 * 2 * [sigma sigma] + 1;
% else
filter_size = 2 * [sigma sigma] + 1;
% end
gFilt = fspecial('gaussian', filter_size, sigma);

pImg = padarray(img, [sigma sigma], 'replicate');
bImg = conv2(pImg, gFilt, 'same');
bImg = bImg(sigma+1:end-sigma, sigma+1:end-sigma);

adjustment = amount * (img - bImg);
adjustment = min(maxAdj, abs(adjustment)) .* sign(adjustment);
img = img + adjustment;

% figure; imshow(img); impixelinfo
% figure; imagesc(adjustment); impixelinfo

% dnFilt = fspecial('gaussian', [5 5], 2);
% img = imfilter(img, dnFilt);

end