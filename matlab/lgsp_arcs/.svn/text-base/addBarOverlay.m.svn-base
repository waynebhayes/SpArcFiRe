function [img, barOvlInds] = ...
    addBarOverlay(img, ctrR, ctrC, angles, halfLengths, allChannels)
% Modifies the given image to overlay the given bar(s), if any
% INPUTS:
%   img: the image to overlay
%   ctrR: the row index of the center of the log-spirals
%   ctrC: the column index of the center of the log-spirals
%   angles: the angles (in radians) of any bars, empty if no bar
%   halfLens: the half-lengths of any bars, empty if no bar
%   allChannels: whether the original (pre-arc-overlay) image is in color
%       if so, bar overlay will be white instead of green (optional,
%       default false)
% OUTPUTS:
%   img: the given image, modified to include overlays for any given bars
%   barOvlInds: the indices of the pixels where the bar overlay was added

if nargin < 6 || isempty(allChannels)
    allChannels = false;
end

if size(img, 3) == 1
    img = repmat(img, [1 1 3]);
end

if length(angles) > 1
    warning('more than one set of bar parameters');
end

overlay = zeros(size(img, 1), size(img, 2));
barOvlInds = cell(1, length(angles));
for bar = 1:1:length(angles)
    angle = angles(bar);
    halfLength = halfLengths(bar);
    rho = -halfLength:0.5:halfLength;
    th = angle * ones(1, length(rho));
    [x, y] = pol2cart(th, rho);
    rows = round(ctrR - y);
    cols = round(ctrC + x);
    out_of_range = (rows < 1) | (rows > size(img, 1)) |...
        (cols < 1) | (cols > size(img, 2));
    rows = rows(~out_of_range);
    cols = cols(~out_of_range);
    if any(out_of_range)
        warning('detected bar extends outside image range');
    end
    curOvlInds = sub2ind(size(img), rows, cols);
    overlay(curOvlInds) = 1;
    barOvlInds{bar} = curOvlInds;
end

overlay3 = repmat(overlay, [1 1 3]);
if allChannels
    img = img + overlay3;
else
    img(overlay3 > 0) = 0;
    img(:, :, 2) = img(:, :, 2) + overlay;
end

if length(barOvlInds) == 1
    barOvlInds = barOvlInds{1};
end

end