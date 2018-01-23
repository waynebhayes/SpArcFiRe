function [isInner, gapFail] = idInnerOuterSpiral(img, ctrR, ctrC, plotFlag)
% Determines which points are in the inner vs outer part of a spiral that 
%   may extend beyond one revolution (in theta-space).  There must be a
%   range of theta-values (even if a small one) where the two revolutions
%   do not overlap each other.
% INPUTS:
%   img:
%   ctrR:
%   ctrC:
%   plotFlag
% OUTPUTS:
%   isInner: Boolean vector indicating whether each nonzero point is
%       assigned to the inner or outer part of the spiral.  The order of
%       the vector corresponds to the order in which the nonzero points
%       would be returned by find().
%   gapFail: true if there is no single-revolution segment that can be used
%       to split into inner and outer regions (in which case, isInner is
%       undefined), false otherwise

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

if length(size(img)) ~= 2
    error('img must be a 2D array');
end

gapFail = false;

nRad = ceil(max(size(img)) / 2);
nTheta = 360;
% thBinVals = [0:nTheta-1]*(2*pi/nTheta);
thBinVals = [1:nTheta]*(2*pi/nTheta);

minAcceptableLength = 5 * ceil(nTheta/360);

% Look for known single-revolution theta-ranges.  Convert image coordinates
% to polar and look for columns (theta-slices) where there is exactly one
% transition from zero to nonzero, and exactly one transition from nonzero
% to zero.  Pad column edges with zero to ensure both these transitions are
% picked up with diff(), even at the edges.
imgp = fliplr(polartrans(img, nRad, nTheta, ctrC, size(img, 1) - ctrR, 'linear', 'full'));
imgp(isnan(imgp)) = 0;
imgp = imdilate(imgp, strel('square', 3));
imgp = padarray(imgp, [1 0]);
imgpD = diff(imgp > 0, 1, 1);
isSgl = (sum(imgpD == -1, 1) == 1) & (sum(imgpD == 1, 1) == 1);

% In addition to the nonzero regions in each theta-bin having exactly one
% start and end point, adjacent theta-bins must be contiguous. We can
% detect this by comparing maxima and minima of adjacent bins. 
% The two criteria are not necessarily the same because the inner or outer
% part could end at theta-bin i and the other could start at theta-bin i+1.
locs = repmat([1:size(imgp, 1)]', 1, size(imgp, 2));
imgpForMin = im2double(imgp > 0); imgpForMin(imgpForMin == 0) = inf;
minLocs = min(imgpForMin .* locs, [], 1); minLocs(minLocs == inf) = 0;
maxLocs = max((imgp > 0) .* locs, [], 1);
nbrMaxLocsL = circshift(maxLocs, [0 1]);
nbrMaxLocsR = circshift(maxLocs, [0 -1]);
nbrMinLocsL = circshift(minLocs, [0 1]);
nbrMinLocsR = circshift(minLocs, [0 -1]);
isSgl = isSgl ...
    & ( (minLocs <= nbrMaxLocsL) | (nbrMaxLocsL == 0) ) ...
    & ( (minLocs <= nbrMaxLocsR)  | (nbrMaxLocsR == 0) ) ...
    & ( (maxLocs >= nbrMinLocsL) | (nbrMinLocsL == 0) ) ...
    & ( (maxLocs >= nbrMinLocsR) | (nbrMinLocsR == 0) );

% Due to quantization and interpolation effects, the "border" theta-bins
% may not entirely be single-revolution, so we shrink the single-revolution
% regions on each side. A shrink amount of 2 is probably sufficient; but we
% apply a little more for safety (at the cost of requiring a slightly
% larger single-rev region).
shrinkAmt = 5;
% If the min acceptable length is reduced below shrinkAmt, we will need to 
% shift in increments of 1..shrinkAmt or use an erosion for the shrink.
isSgl = isSgl & circshift(isSgl, shrinkAmt * [1 1]) ...
    & circshift(isSgl, shrinkAmt * [-1 -1]);

if plotFlag
    figure; 
    subplot(2, 1, 1);
    imshow([flipud(imgp); repmat(isSgl, 5, 1)]);
    title('sgl-rev regions');
    subplot(2, 1, 2);
    imshow([flipud(imgp > 0); repmat(isSgl, 5, 1)]);
    title('sgl-rev regions');
    impixelinfo
end

isSglD = diff(isSgl);

stIdxs = find(isSglD == 1) + 1;
endIdxs = find(isSglD == -1);

hasWrap = false;
wrapStart = 0;
wrapEnd = 0;
if isempty(stIdxs) && isempty(endIdxs)
    if all(isSgl)
        % sgl-rev for the entire theta-range, so no endpoints picked up
        % (this could be a ring, but more likely it's a cluster in the
        % center)
        stIdxs = 1;
        endIdxs = length(isSgl);
    else
        warning('no sgl-rev regions in entire theta-range');
        isInner = true(size(find(img)));
        gapFail = true;
        return;
    end
end
if ~isempty(endIdxs) && (isempty(stIdxs) || (stIdxs(1) > endIdxs(1)) )
    % one of the sgl-rev regions wraps around; we'll deal with it later
    hasWrap = true;
    wrapEnd = endIdxs(1);
    endIdxs = endIdxs(2:end);
    assert(isempty(stIdxs) || isempty(endIdxs) || stIdxs(1) <= endIdxs(1));
end
if ~isempty(stIdxs) && (isempty(endIdxs) || (stIdxs(end) > endIdxs(end)) )
    hasWrap = true;
    wrapStart = stIdxs(end);
    stIdxs = stIdxs(1:end-1);
    assert(isempty(stIdxs) || isempty(endIdxs) || stIdxs(end) <= endIdxs(end));
end
assert(length(stIdxs) == length(endIdxs))
wrapLgth = 0;
if hasWrap
    if wrapStart == 0
        % last continuous sgl-rev is at the beginning, but doesn't 
        % actually wrap around
        stIdxs = [1 stIdxs];
        endIdxs = [wrapEnd endIdxs];
        hasWrap = false;
    elseif wrapEnd == 0
        % last continuous sgl-rev is at the end, but doesn't
        % actually wrap around
        stIdxs = [wrapStart stIdxs];
        endIdxs = [length(isSgl) endIdxs];
        hasWrap = false;
    else
        wrapStL = length(isSgl) - wrapStart + 1;
        wrapEndL = wrapEnd - 1 + 1;
        wrapLgth = wrapStL + wrapEndL;
    end
end

[rows, cols] = find(img);
ptIdxToImgIdx = sub2ind(size(img), rows, cols);
imgIdxToPtIdx(img > 0) = 1:1:length(rows);
[thVals, rVals] = cart2pol(cols - ctrC, -(rows - ctrR));
thVals = mod(thVals, 2*pi);

% use theta-values to split known sgl-rev range into rgn1 and rgn2
% after assigning all of the other points to one of these regions, these
% regions will become the inner and outer areas, but we don't know
% which is which yet
srLgths = endIdxs - stIdxs + 1;
[mxLgth, mlI] = max(srLgths);
if isempty(stIdxs) || (wrapLgth > mxLgth)
    maxLength = wrapLgth;
    wrapMidL = round(wrapLgth/2);
    thStart = thBinVals(wrapStart);
    thEnd = thBinVals(wrapEnd);
    inRgn = (thVals >= thStart) | (thVals < thEnd);
    if wrapStL >= wrapMidL 
        splitTh = thBinVals(wrapStart + wrapMidL - 1);
        inRgn1 = inRgn & (thVals >= thStart) & (thVals < splitTh);
        inRgn2 = inRgn & ~inRgn1;
    else
        splitTh = thBinVals(wrapMidL - wrapStL);
        inRgn1 = inRgn & (thVals >= splitTh) & (thVals < thEnd);
        inRgn2 = inRgn & ~inRgn1;
    end
else
    maxLength = mxLgth;
    thStart = thBinVals(stIdxs(mlI));
    thEnd = thBinVals(endIdxs(mlI));
    splitTh = thBinVals(round( (stIdxs(mlI) + endIdxs(mlI)) / 2 ));
    
    inRgn1 = (thVals >= thStart) & (thVals < splitTh);
    inRgn2 = (thVals >= splitTh) & (thVals < thEnd);
end
if maxLength < minAcceptableLength
%     figure; imshow(img); 
%     title('longest sgl-rev region length is below the minimum length');
    fprintf('Warning:idInnerOuterSpiral:longest sgl-rev region length (%d) is below the minimum length (%d)\n', ...
        maxLength, minAcceptableLength);
    isInner = true(size(find(img)));
    gapFail = true;
    return;
end

rgn1img = false(size(img)); rgn1img(ptIdxToImgIdx(inRgn1)) = true;
rgn2img = false(size(img)); rgn2img(ptIdxToImgIdx(inRgn2)) = true;
nonRgnImg = false(size(img)); nonRgnImg(ptIdxToImgIdx(~inRgn1 & ~inRgn2)) = true;

if plotFlag
    figure; subplot(2, 2, 1); 
    hold on
    imgRad = max(size(img))/2;
    imagesc(imgRad * [-1 1], imgRad * [1 -1], img); axis image; colormap gray

    hold on
    polar([0 eps] + splitTh, [0 max(size(img))/2]);
    polar([0 eps] + thStart, [0 max(size(img))/2], 'b:');
    polar([0 eps] + thEnd, [0 max(size(img))/2], 'b:');
    title(sprintf('thStart = %2.4f, splitTh = %2.4f, thEnd = %2.4f', ...
        thStart, splitTh, thEnd));
    hold off

    rgnAsgnImg = zeros([size(img) 3]);
    rgnAsgnImg(:, :, 1) = rgn1img; 
    rgnAsgnImg(:, :, 2) = rgn2img; 
    rgnAsgnImg(:, :, 3) = nonRgnImg;
    subplot(2, 2, 2); imshow(rgnAsgnImg); title('region assignment (initial)');
end

% If a connected component in the unassigned pixels is adjacent to one of
% the assigned regions, make it part of the region the component is
% adjacent to
r1dist = bwdist(rgn1img);
r2dist = bwdist(rgn2img);
cc = bwconncomp(nonRgnImg);
maxDiagDist = 1.5;
for ii=1:1:cc.NumObjects
    curInds = cc.PixelIdxList{ii};
    dist1 = min(r1dist(curInds));
    dist2 = min(r2dist(curInds));
    if dist1 < maxDiagDist % < sqrt(2)
        inRgn1(imgIdxToPtIdx(curInds)) = true;
    elseif dist2 < maxDiagDist
        inRgn2(imgIdxToPtIdx(curInds)) = true;
    end
end

rgn1img = false(size(img)); rgn1img(ptIdxToImgIdx(inRgn1)) = true;
rgn2img = false(size(img)); rgn2img(ptIdxToImgIdx(inRgn2)) = true;
nonRgnImg = false(size(img)); nonRgnImg(ptIdxToImgIdx(~inRgn1 & ~inRgn2)) = true;

if plotFlag
    rgnAsgnImg = zeros([size(img) 3]);
    rgnAsgnImg(:, :, 1) = rgn1img; 
    rgnAsgnImg(:, :, 2) = rgn2img; 
    rgnAsgnImg(:, :, 3) = nonRgnImg;
    subplot(2, 2, 3); imshow(rgnAsgnImg);
    title('region assignment (after handling contiguous regions)');
end

% Assign the remaining pixels according to closest distance to one of the
% two regions
r1dist = bwdist(rgn1img);
r2dist = bwdist(rgn2img);
cc = bwconncomp(nonRgnImg);
idx_lists = cc.PixelIdxList;
while ~isempty(idx_lists)
    dists1 = cellfun(@(x)(min(r1dist(x))), idx_lists);
    dists2 = cellfun(@(x)(min(r2dist(x))), idx_lists);
    [mV, minDistIdx] = min(min(dists1, dists2));
    curInds = idx_lists{minDistIdx};
    if dists1(minDistIdx) < dists2(minDistIdx)
        inRgn1(imgIdxToPtIdx(curInds)) = true;
        rgn1img(curInds) = true;
        r1dist = bwdist(rgn1img);
    else
        inRgn2(imgIdxToPtIdx(curInds)) = true;
        rgn2img(curInds) = true;
        r2dist = bwdist(rgn2img);
    end
    idx_lists(minDistIdx) = [];
end

assert(all(xor(inRgn1, inRgn2)))

% now one of the regions is the inner pixels and the other is the outer
% pixel; we simply use the average radius-value to figure out which is
% which
if mean(rVals(inRgn1)) < mean(rVals(inRgn2))
    isInner = inRgn1;
else
    isInner = inRgn2;
end

innerImg = false(size(img)); innerImg(ptIdxToImgIdx(isInner)) = true;
outerImg = false(size(img)); outerImg(ptIdxToImgIdx(~isInner)) = true;
if plotFlag
    rgnAsgnImg = zeros([size(img) 3]);
    rgnAsgnImg(:, :, 1) = outerImg;
    rgnAsgnImg(:, :, 2) = innerImg;
    rgnAsgnImg(:, :, 3) = innerImg + outerImg;
    subplot(2, 2, 4); imshow(rgnAsgnImg);
    title('region assignment (final)');
end

% gapFail
% resImg = img; resImg(img > 0) = 2*isInner - 1;
% figure; imagesc(resImg); axis image

end % idInnerOuterSpiral