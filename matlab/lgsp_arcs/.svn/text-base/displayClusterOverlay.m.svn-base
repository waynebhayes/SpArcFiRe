function clusOvlImg = displayClusterOverlay(img, clusMtxs, overlayColor)
% Displays or returns an image with cluster boundaries overlaid

if length(size(img)) ~= 2
    error('img must be a grayscale image');
end

if nargin < 3 || isempty(overlayColor)
    overlayColor = [0 1 1];
end

% allClus = sum(clusMtxs, 3) > 0;
% clusBounds = allClus & (imerode(allClus, strel('square', 3)) == 0);

clusBounds = false(size(img));
for cIdx=1:1:size(clusMtxs, 3)
    curClusMtx = clusMtxs(:, :, cIdx) > 0;
%     clusBounds = clusBounds | ...
%         curClusMtx & ~(imerode(curClusMtx, strel('square', 3)));
    clusBounds = clusBounds | ...
        imdilate(curClusMtx, strel('square', 3)) & ~curClusMtx;
end

img(clusBounds) = 0;
clusOvlImg = repmat(img, [1 1 3]);
clusOvlImg(:, :, 1) = clusOvlImg(:, :, 1) + clusBounds * overlayColor(1);
clusOvlImg(:, :, 2) = clusOvlImg(:, :, 2) + clusBounds * overlayColor(2);
clusOvlImg(:, :, 3) = clusOvlImg(:, :, 3) + clusBounds * overlayColor(3);

if nargout < 1
    figure; imshow(clusOvlImg);
end

end