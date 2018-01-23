function arcImg = displayLgspOverlay(img, lsParams, ctrR, ctrC, thBounds)
% Displays or returns an input image combined with arcs from 
%   logarithmic spirals
% INPUTS:
%   img: the image to overlay
%   lsParams: parameters of the log-spirals, one per row
%   ctrR: the row position of the center of the log-spirals
%   ctrC: the column position of the center of the log-spirals
%   thBounds: starting and ending theta-range of points to which
%       log-spirals were fitted, counterclockwise from the theta-offset
%       parameter.  Each row is a [start, end] pair.
% OUTPUTS:
%   arcImg: the overlay image (displayed to screen if not requested as
%       function output)

nClus = size(lsParams, 1);
ctrX = ctrC;
ctrY = size(img, 1) - ctrR + 1;

if nargin < 5 || isempty(thBounds)
    thBounds = repmat([0, 2*pi], nClus, 1);
end

imgSz = size(img);
colorImg = (ndims(img) == 3);
if colorImg
    imgSz = imgSz(1:2);
end

arcsCW = zeros(imgSz);
arcsCCW = zeros(imgSz);
arcsN = zeros(imgSz);
for clus = 1:1:nClus
    if lsParams(clus, 2) > 0
        arcsCW = arcsCW + ...
        polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, ...
        lsParams(clus, :), thBounds(clus, :));
    elseif lsParams(clus, 2) < 0
        arcsCCW = arcsCCW + ...
        polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, ...
        lsParams(clus, :), thBounds(clus, :));
    else
        arcsN = arcsN + ...
        polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, ...
        lsParams(clus, :), thBounds(clus, :));
    end
end

% figure; imshow(arcs > 0);

arcOverlay = zeros([imgSz 3]);
arcOverlay(:, :, 1) = (arcsCW > 0) | (arcsN > 0);
arcOverlay(:, :, 3) = (arcsCCW > 0) | (arcsN > 0);
% arcOverlay(:, :, 3) = img;

if colorImg
    arcOverlay = sum(arcOverlay, 3) > 0;
end

% arcImg = arcOverlay;
arcPixels = sum(arcOverlay, 3) > 0;
bkgdImg = img;
if colorImg
    arcImg = bkgdImg;
    arcImg(repmat(arcPixels, [1 1 3])) = 1;
else
    bkgdImg(arcPixels) = 0;
    bkgdImg = repmat(bkgdImg, [1 1 3]);
    arcImg = arcOverlay + bkgdImg;
end

if nargout < 1
    figure;  imshow(arcImg);
end

end