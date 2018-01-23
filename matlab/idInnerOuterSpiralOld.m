function [isInner, pt2slice, numSlices] = ...
    idInnerOuterSpiralOld(img, ctrR, ctrC, plotFlag, inrOtrPref, forceSplit)
% Obsolete; do not use
% Heuristically attempts to determine which points are in the inner vs 
%   outer part of a spiral that may extend beyond one revolution (in
%   theta-space).  This assumes that there is a large-enough gap in
%   r-values between points with similar theta-values that are in different
%   revolutions of the spiral.  This function may not work if the arc
%   spans close to 2 full revolutions.
% INPUTS:
%   img:
%   ctrR:
%   ctrC: 
%   plotFlag
%   inrOtrPref: -1 to have as many inner points as possible, 1 to have as
%       many outer points as possible, 0 for no preference (0 recommended
%       in order to help avoid the theta-ranges of the inner and outer
%       segments getting close to 2*pi)
%   forceSplit: force an inner/outer split even if no second revolution is
%       detected
% OUTPUTS:
%   isInner: vector with entries corresponding to input r- and
%       theta-values; entries are true if the point is considered an inner
%       point and false otherwise.  If there is no detected second
%       revolution in the spiral, all points will be considered inner
%       points.
%   pt2slice: theta-slices each point was assigned to
%   numSlices: number of theta-slices that were used

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

if nargin < 5 || isempty(inrOtrPref)
    inrOtrPref = 0;
end

if nargin < 6 || isempty(forceSplit)
    forceSplit = false;
end

% TODO: redo/simplify this entire function by searching for a cutoff in
% polar space

% smaller slices require more function calls, result-merges 
% (because there are more such slices), but larger slices risk missed
% detections
numSlices = 16;
sliceAng = (2*pi) / numSlices;

[rows, cols] = find(img);
% inds also serves as ptIdxToImgIdx
inds = sub2ind(size(img), rows, cols);
[thVals, rVals] = cart2pol(cols - ctrC, -(rows - ctrR));
thVals = mod(thVals, 2*pi);
imgIdxToPtIdx = zeros(size(img));
pxIsNz = (img > 0);
imgIdxToPtIdx(pxIsNz) = 1:1:length(thVals);

imgCC = bwconncomp(pxIsNz);
isInner = false(size(inds));
pt2slice = zeros(size(inds));
if imgCC.NumObjects > 1
    warning('non-contiguous cluster');
    indInCC = cellfun(@(x)(ismember(inds, x)), imgCC.PixelIdxList, 'UniformOutput', false);
    ccSizes = cellfun(@nnz, indInCC);
    ccAvgR = cellfun(@(x)(mean(rVals(x))), indInCC);
    [maxCCSz, maxCCIdx] = max(ccSizes);
%     maxCCAvgR = ccAvgR(maxIdx);
%     % non-largest components favor inner if their average r-value is less
%     % than that of the largest component, and prefer outer otherwise
%     ccPrefs = sign(ccAvgR - maxCCAvgR);
%     % largest component prefers inner if most other components have higher
%     % average r-values, and prefers outer otherwise
%     ccPrefs(maxIdx) = -sign(sum(ccPrefs .* ccSizes));
    ccPrefs = zeros(size(ccSizes));
    
    ccLbls = labelmatrix(imgCC);
    for ccIdx=1:1:length(ccSizes)
        % TODO: for the non-largest connected components; we're really only
        % doing this for the pt2slice values; this calculation should be 
        % factored out
        [ccIsInner, ccPt2slice] = idInnerOuterSpiral(...
            img .* (ccLbls == ccIdx), ctrR, ctrC, false, ccPrefs(ccIdx));
        if ccIdx == maxCCIdx
            if all(ccIsInner == ccIsInner(1))
                % largest connected component isn't split, which may cause
                % all points in all connected components to have the same
                % inner/outer designation
                % so we redo, forcing a split
                [ccIsInner, ccPt2slice] = idInnerOuterSpiral(...
                    img .* (ccLbls == ccIdx), ctrR, ctrC, false, ...
                    ccPrefs(ccIdx), true);
            end
        end
        isInner(indInCC{ccIdx}) = ccIsInner;
        pt2slice(indInCC{ccIdx}) = ccPt2slice;
    end
    indInMaxCC = indInCC{maxCCIdx};
%     indMccInner = inds(indInMaxCC(isInner(indInMaxCC)));
%     indMccOuter = inds(indInMaxCC(~isInner(indInMaxCC)));
    indMccInner = inds(indInMaxCC); indMccInner = indMccInner(isInner(indInMaxCC));
    indMccOuter = inds(indInMaxCC); indMccOuter = indMccOuter(~isInner(indInMaxCC));
    mccInnerImg = zeros(size(img)); mccInnerImg(indMccInner) = 1;
    mccOuterImg = zeros(size(img)); mccOuterImg(indMccOuter) = 1;
%     mccImg = zeros(size(img)); mccImg(inds(indInMaxCC)) = 1;
%     figure; imshow(mccImg); title('mccImg');
%     figure; subplot(1, 2, 1); imshow(mccInnerImg); title('mcc inner'); subplot(1, 2, 2); imshow(mccOuterImg); title('mcc outer');
    
    innerDist = bwdist(mccInnerImg);
    outerDist = bwdist(mccOuterImg);
    for ccIdx=1:1:length(ccSizes)
        if ccIdx == maxCCIdx
            continue;
        end
        ccImgInds = inds(indInCC{ccIdx});
        isInner(indInCC{ccIdx}) = innerDist(ccImgInds) < outerDist(ccImgInds);
    end

    if plotFlag
        [xVals, yVals] = pol2cart(thVals, rVals);
        figure; hold on; axis equal;
        scatter(xVals(isInner), yVals(isInner), 'r.');
        scatter(xVals(~isInner), yVals(~isInner), 'b.');
        for slcNum = 1:1:numSlices
            polar((slcNum-1)*sliceAng + [0 eps], [0 max(rVals)], '-k');
        end
    end
    
    if plotFlag
        thAdj = ...
            removeThetaDiscontFor2rev(rVals, thVals, img, isInner, pt2slice, numSlices);
        thImg = zeros(size(img));
        thImg(inds) = thAdj;
        figure; imagesc(thImg); axis image
    end
    
    if plotFlag()
        plotIsInner();
    end
    return;
end

pt2slice = zeros(size(rVals));

% if there is only one connected component in a slice, we know the slice is
% part of only one revolution
% slcCCs = cell(numSlices, 1); % do we really need to store these, or just the counts?
slcCCcounts = zeros(numSlices, 1);
slcImgs = zeros([size(img) numSlices]);
nInSlice = zeros(numSlices, 1);
isSglRev = false(size(slcCCcounts));
emptySlc = false(size(slcCCcounts));
for slcNum = 1:1:numSlices
    slcStart = (slcNum-1)*sliceAng;
    slcEnd = slcNum*sliceAng;

    inSlice = (thVals >= slcStart) & (thVals < slcEnd);
    pt2slice(inSlice) = slcNum;
    
    if plotFlag
        fprintf('slice %d: [%2.4f, %2.4f)\n', slcNum, slcStart, slcEnd);
    end
    
    updateSliceInfo(slcNum);
end
% technically, empty slices are also sgl-rev, but this could mess up our
% propagation
% ~isSglRev' & ~emptySlc'

% thSorted = sort(thVals, 'ascend');
% maxThGap = max([diff(thSorted); thSorted(1) + 2*pi - thSorted(end)])
% if maxThGap > sliceAng
%     fprintf('sgl-rev arc\n');
%     if inrOtrPref <= 0
%         isInner = true(size(inds));
%     else
%         isInner = false(size(inds));
%     end
%     return;
% end

% diff(isSglRev([end 1:end 1]))'

if plotFlag
    plotSlices(); title('initial slices');
end

changeMade = true;
while changeMade
    changeMadeDbl = doDblSidedProp();
    changeMadeSgl = doSglSidedProp();
    changeMade = changeMadeSgl || changeMadeDbl;
end
fprintf('done with connected-region propagation\n');
% isSglRev'

if plotFlag
    plotSlices(); title('slices after propagation');
end

% use connected components info to determine whether each point is inner,
% outer, or could be either
inrOtrDes = zeros(size(rVals));
for slcIdx = 1:1:numSlices
    inSlice = (pt2slice == slcIdx);
    if isSglRev(slcIdx) || sum(inSlice) <= 1 % second condition unnecessary?
        inrOtrDes(inSlice) = 0;
    else
        slcR = rVals(inSlice);
        
        % split points by largest gap in r-values
        [srtR, rSrtInds] = sort(slcR, 'ascend');
        rGaps = diff(srtR);
        [srtGaps, gapSrtInds] = sort(rGaps, 'ascend'); % just get max?
        interGapSize = srtGaps(end); 
        gapSortPosn = gapSrtInds(end);
        
        slcInner = true(size(slcR));
        slcInner(rSrtInds(gapSortPosn+1:end)) = false;
        slcDes = zeros(size(slcR));
        slcDes(slcInner) = -1;
        slcDes(~slcInner) = 1;
        inrOtrDes(inSlice) = slcDes;
    end
end

% % how many more times wider the largest gap has to be,
% % vs the r-width of inner/outer regions, in order for the slice to be 
% % considered to have points from 2 revolutions of a spiral
% gapRatioThres = 1;  
% 
% inrOtrDes = zeros(size(rVals));
% 
% numSlicesWithInner = 0;
% numSlicesWithOuter = 0;
% sliceIsEither = false(1, numSlices);
% sliceWasSplit = false(1, numSlices);
% for slcNum = 1:1:numSlices
%     slcStart = (slcNum-1)*sliceAng;
%     slcEnd = slcNum*sliceAng;
%     if plotFlag
%         fprintf('slice %d: [%2.4f, %2.4f)\n', slcNum, slcStart, slcEnd);
%     end
% 
%     inSlice = (thVals >= slcStart) & (thVals < slcEnd);
%     pt2slice(inSlice) = slcNum;
%     slcDes = determineForSlice(inSlice);
%     inrOtrDes(inSlice) = slcDes;
%     if sum(slcDes == -1) > 0
%         numSlicesWithInner = numSlicesWithInner + 1;
%     end
%     if sum(slcDes == 1) > 0
%         numSlicesWithOuter = numSlicesWithOuter + 1;
%     end
% end

% % use neighbors to fix false splits
% % this won't always work when more than one contiguous slice is falsely
% % split, but most false splits occur at the tips of the spiral, where there
% % are not very many points in the slice
% chkSlices = find(sliceWasSplit)
% for ii=1:1:length(chkSlices)
%     slcIdx = chkSlices(ii);
%     nbrs = mod(slcIdx + [-1 1] - 1, numSlices) + 1;
%     if ~sliceIsEither(nbrs(1)) || (nInSlice(nbrs(1)) == 0)
%         nbrs(1) = 0;
%     end
%     if ~sliceIsEither(nbrs(2)) || (nInSlice(nbrs(2)) == 0)
%         nbrs(2) = 0;
%     end
%     nbrs = nbrs(nbrs > 0);
%     if isempty(nbrs)
%         continue;
%     end
%     nbrs
%     sum(ismember(pt2slice, [slcIdx nbrs]))
%     nbrImg = false(size(img));
%     nbrImg(inds(ismember(pt2slice, [slcIdx nbrs]))) = true;
%     figure; imshow(nbrImg);
%     ccs = bwconncomp(nbrImg);
%     ccs = ccs.PixelIdxList;
%     % if all pixels in the split slice share a connected component with
%     % non-split slices
% end

if plotFlag
    [xVals, yVals] = pol2cart(thVals, rVals);
    figure; hold on; axis equal;
    scatter(xVals(inrOtrDes == -1), yVals(inrOtrDes == -1), 'r.');
    scatter(xVals(inrOtrDes == 0), yVals(inrOtrDes == 0), 'g.');
    scatter(xVals(inrOtrDes == 1), yVals(inrOtrDes == 1), 'b.');
    for slcNum = 1:1:numSlices
        polar((slcNum-1)*sliceAng + [0 eps], [0 max(rVals)], '-k');
    end
end

% sliceIsEither = ismember(1:1:numSlices, pt2slice(inrOtrDes == 0))';
sliceIsEither = isSglRev;
numSlicesWithInner = sum(~isSglRev & (nInSlice > 0));
numSlicesWithOuter = numSlicesWithInner;

if sum(sliceIsEither) == numSlices
    % there's no second revolution of the spiral
    fprintf('all inner\n');
    if forceSplit
        fprintf('forcing a split anyway (due to input parameters)\n');
        isInner = createForcedSplit();
    elseif inrOtrPref <= 0
        isInner = true(size(rVals));
    else
        isInner = false(size(rVals));
    end
    if plotFlag()
        plotIsInner();
    end
    return;
end

% Set points to be either inner or outer.  Points that already have a
% designation keep that designation.  Points that could be either are
% split.
isInner = (inrOtrDes == -1);
% extraSlices = find(sliceIsEither);
% size(sliceIsEither)
% size(nInSlice)
extraSlices = find(sliceIsEither & (nInSlice > 0));
% don't include extra slices without any points
extraSlices = extraSlices(ismember(extraSlices, unique(pt2slice))); % necessary?
numExtraSlices = length(extraSlices);
hasExtraSlices = ~isempty(extraSlices);

if hasExtraSlices
    % rotate indices of extra slices so that the discontinuity is at the
    % array bounds
    [val, ind] = max(diff(extraSlices));
    if val > 1
        newPosns = mod([1:numExtraSlices] + ind - 1, numExtraSlices) + 1;
        extraSlices = extraSlices(newPosns);
    end

    % figure out which side of the undesignated slices should be designated
    % inner, by looking at which edge slice has points closer to the center
    innerFirst = mean(rVals(pt2slice == extraSlices(1))) < ...
            mean(rVals(pt2slice == extraSlices(end)));
    discontIdxs = find((diff(extraSlices)) ~= 1 & ...
        (abs(diff(extraSlices)) ~= (numSlices - 1)));
end

if sum(sliceIsEither) == sum(nInSlice > 0)
    fprintf('no splits needed\n');
    % no slices where points need to be split into inner and outer
    if forceSplit
        fprintf('forcing a split anyway (due to input parameters)\n');
        isInner = createForcedSplit();
    elseif inrOtrPref <= 0
        isInner(:) = true;
    else
        isInner(:) = false;
    end
elseif hasExtraSlices && sum(discontIdxs) == 0
    fprintf('one continuous gap\n');
    % one continuous gap of undesignated points - make all slices but
    % one inner slices
    if inrOtrPref == 0
        numExtraInner = ceil(numExtraSlices/2);
    elseif inrOtrPref < 0
        numExtraInner = numExtraSlices - 1;
    else
        numExtraInner = 1;
    end
    
    if innerFirst
        extraInnerSlices = extraSlices(1:numExtraInner);
    else
        extraInnerSlices = extraSlices(end-numExtraInner+1:end);
    end

%     numSlicesWithInner = numSlicesWithInner + numExtraInner;
%     numSlicesWithOuter = numSlicesWithOuter + (numExtraSlices - numExtraInner);
    isInner(ismember(pt2slice, extraInnerSlices)) = true;
elseif hasExtraSlices
    fprintf('multiple gaps\n');
    warning('multiple split regions');
    % multiple gaps of undesignated points - we'll have to find the
    % designations by looking at neighboring slices (propagating the
    % information from there)

    % this propagation prefers the designation of the slice preceding
    % the undesignated slice we started with, so start with the
    % designation that has fewer slices
    if ((numSlicesWithInner < numSlicesWithOuter) && ~innerFirst) ||...
        ((numSlicesWithInner > numSlicesWithOuter) && innerFirst)
        extraSlices = extraSlices(end:-1:1);
    end

    numSwitches = 0;
    for eslcIdx = 1:1:length(extraSlices)
        eslc = extraSlices(eslcIdx);
        beforeSlice = mod((eslc - 1) - 1, numSlices) + 1;
        afterSlice = mod((eslc + 1)  - 1, numSlices) + 1;
        isNearby = ismember(pt2slice, [beforeSlice, afterSlice]);
        sliceMeanR = mean(rVals(pt2slice == eslc));
        nearbyInnerR = rVals(isNearby & inrOtrDes == -1);
        nearbyOuterR = rVals(isNearby & inrOtrDes == 1);
        if isempty(nearbyInnerR)
            slcIsInner = 0;
        elseif isempty(nearbyOuterR)
            slcIsInner = 1;
        elseif abs(sliceMeanR - mean(nearbyInnerR)) < ...
                abs(sliceMeanR - mean(nearbyOuterR))
            slcIsInner = 1;
        else
            slcIsInner = 0;
        end

        if slcIsInner
            isInner(pt2slice == eslc) = true;
            inrOtrDes(pt2slice == eslc) = -1;
            numSlicesWithInner = numSlicesWithInner + 1;
        else
            inrOtrDes(pt2slice == eslc) = 1;
            numSlicesWithOuter = numSlicesWithOuter + 1;
        end

        if eslcIdx ~= 1 && (prevWasInner ~= slcIsInner)
            numSwitches = numSwitches + 1;        
        end
        prevWasInner = slcIsInner;     
    end
%     if numSwitches > 1
%         warning(['switched between inner/outer points more '...,
%                 'than once, designation may be inaccurate'])
%     end
end
% done assigning inner/outer to points in undesignated slices

% in a few cases, a few points may still be assigned incorrectly
imgInner = pxIsNz;
imgInner(pxIsNz) = isInner;
ccInner = bwconncomp(imgInner);
ccInnerSizes = cellfun(@numel, ccInner.PixelIdxList);
nccInner = length(ccInnerSizes);
[temp, mainIdxInner] = max(ccInnerSizes);
if nccInner > 1
    ccOKinner = false(1, nccInner);
    ccOKinner(mainIdxInner) = true;
    mainSlicesInner = unique(pt2slice(imgIdxToPtIdx(ccInner.PixelIdxList{mainIdxInner})));
%     for ii=1:1:nccInner
%         if ii == mainIdxInner
%             continue;
%         end
%     end
    curSlices = unique(pt2slice(imgIdxToPtIdx(ccInner.PixelIdxList{ii})));
    if sum(ismember(curSlices, mainSlicesInner)) == length(curSlices)
        ccOKinner(ii) = true;
    end
else
    ccOKinner = true;
end
imgOuter = pxIsNz - imgInner;
ccOuter = bwconncomp(imgOuter);
ccOuter.NumObjects
ccOuterSizes = cellfun(@numel, ccOuter.PixelIdxList);
nccOuter = length(ccOuterSizes);
[temp, mainIdxOuter] = max(ccOuterSizes);
if nccOuter > 1
    ccOKouter = false(1, nccOuter);
    ccOKouter(mainIdxOuter) = true;
    mainSlicesOuter = unique(pt2slice(imgIdxToPtIdx(ccOuter.PixelIdxList{mainIdxOuter})));
%     for ii=1:1:nccOuter
%         if ii == mainIdxOuter
%             continue;
%         end
%     end
    curSlices = unique(pt2slice(imgIdxToPtIdx(ccOuter.PixelIdxList{ii})));
    if sum(ismember(curSlices, mainSlicesOuter)) == length(curSlices)
        ccOKouter(ii) = true;
    end
else
    ccOKouter = true;
end
if sum(~ccOKinner) > 0 || sum(~ccOKouter) > 0
%     iiImg = img; iiImg(iiImg > 0) = (2*isInner - 1); figure; imagesc(iiImg); axis image; title('isInner before repair');
    warning('inner-outer designation non-contiguous, attempting repair.');
    rValsInner = rVals(imgIdxToPtIdx(ccInner.PixelIdxList{mainIdxInner}));
    rValsOuter = rVals(imgIdxToPtIdx(ccOuter.PixelIdxList{mainIdxOuter}));
    for ii=find(~ccOKinner)
        curInds = imgIdxToPtIdx(ccInner.PixelIdxList{ii});
        curMeanR = mean(rVals(curInds));
        if mean(abs(curMeanR - rValsInner)) > mean(abs(curMeanR - rValsOuter))
            assert(sum(isInner(curInds)) == length(curInds));
            isInner(curInds) = 0;
%             fprintf('reassigned %d points to inner\n', length(curInds));
%             chgImg = zeros(size(img)); chgImg(inds(curInds)) = 1; figure; imshow(chgImg); title('reassigned to outer');
        end
    end
    for ii=find(~ccOKouter)
        curInds = imgIdxToPtIdx(ccOuter.PixelIdxList{ii});
        curMeanR = mean(rVals(curInds));
        if mean(abs(curMeanR - rValsInner)) < mean(abs(curMeanR - rValsOuter))
            assert(sum(isInner(curInds)) == 0);
            isInner(curInds) = 1;
%             fprintf('reassigned %d points to outer\n', length(curInds));
%             chgImg = zeros(size(img)); chgImg(inds(curInds)) = 1; figure; imshow(chgImg); title('reassigned to inner');
        end
    end
    
    imgInner = false(size(img));
    imgInner(pxIsNz) = isInner;
    ccInner = bwconncomp(imgInner);
    nccInner = ccInner.NumObjects;
    ccInnerSizes = cellfun(@numel, ccInner.PixelIdxList);
    [temp, mainIdxInner] = max(ccInnerSizes);
    
    imgOuter = false(size(img));
    imgOuter(pxIsNz) = ~isInner;
    ccOuter = bwconncomp(imgOuter);
    nccOuter = ccOuter.NumObjects;
    ccOuterSizes = cellfun(@numel, ccOuter.PixelIdxList);
    [temp, mainIdxOuter] = max(ccOuterSizes);
    
    imgMccInner = false(size(img)); 
    imgMccInner(ccInner.PixelIdxList{mainIdxInner}) = true;
    distToInner = bwdist(imgMccInner);
    
    imgMccOuter = false(size(img));
    imgMccOuter(ccOuter.PixelIdxList{mainIdxOuter}) = true;
    distToOuter = bwdist(imgMccOuter);
    
    for ii=1:1:nccInner
        if ii == mainIdxInner
            continue;
        end
        curImgPts = ccInner.PixelIdxList{ii};
        if mean(distToInner(curImgPts)) > mean(distToOuter(curImgPts))
            isInner(imgIdxToPtIdx(curImgPts)) = false;
        end
    end
    for ii=1:1:nccOuter
        if ii == mainIdxOuter
            continue;
        end
        curImgPts = ccOuter.PixelIdxList{ii};
        if mean(distToOuter(curImgPts)) > mean(distToInner(curImgPts))
            isInner(imgIdxToPtIdx(curImgPts)) = true;
        end
    end
end

% innerLength = sliceAng * numSlicesWithInner;
% outerLength = sliceAng * numSlicesWithOuter;

if plotFlag
    [xVals, yVals] = pol2cart(thVals, rVals);
    figure; hold on; axis equal;
    scatter(xVals(isInner), yVals(isInner), 'r.');
    scatter(xVals(~isInner), yVals(~isInner), 'b.');
    for slcNum = 1:1:numSlices
        polar((slcNum-1)*sliceAng + [0 eps], [0 max(rVals)], '-k');
    end
%     title(sprintf('innerLength = %2.4f, outerLength = %2.4f',...
%         innerLength, outerLength));
end

function updateSliceInfo(slcIdx)
    inSlice = (pt2slice == slcIdx);
    
    nInSlice(slcIdx) = sum(inSlice);
    
    slcImg = (img > 0);
    slcImg(inds(~inSlice)) = 0;
%     figure; imshow(slcImg);
    slcImgs(:, :, slcIdx) = slcImg;
%     figure; imshow(slcImg); title(sprintf('%d', slcIdx));
    
    ccs = bwconncomp(slcImg);
%     figure; imagesc(labelmatrix(ccs)); axis image
%     slcCCs{slcIdx} = find(ismember(inds, ccs.PixelIdxList));
    slcCCcounts(slcIdx) = ccs.NumObjects;
    isSglRev(slcIdx) = (slcCCcounts(slcIdx) == 1);
    emptySlc(slcIdx) = (slcCCcounts(slcIdx) == 0);
    
    if plotFlag
        fprintf('\t%d connected components in slice\n', ccs.NumObjects);
    end
end

function changeMade = doSglSidedProp()
% propagate information from known-single-rev regions: if all connected
% components of a slice are connected to something from a *single*
% known-single-rev region, that slice is also a known-single-rev region
srBnd = find(~isSglRev & (circshift(isSglRev, [1 1]) | circshift(isSglRev, [-1 -1])));
fprintf('performing single-sided propagation...\n');
% srBnd'
changeMade = false;
while ~isempty(srBnd)
    isNxtBnd = false(size(isSglRev));

    for ii=1:1:length(srBnd)
        bIdx = srBnd(ii);
        if slcCCcounts(bIdx) == 0
            continue;
        end
        nbrs = mod(bIdx + [-1 1] - 1, numSlices) + 1;
        for jj = 1:1:length(nbrs)
            nbIdx = nbrs(jj);
            if ~isSglRev(nbIdx)
                continue;
            end
%             figure; imshow(sum(slcImgs(:, :, [bIdx nbIdx]), 3)); title(sprintf('bIdx = %d, nbIdx = %d', bIdx, nbIdx));
            retry = true;
            while retry
                retry = false;
                nbrCCs = bwconncomp(sum(slcImgs(:, :, [bIdx nbIdx]), 3));
                if nbrCCs.NumObjects == slcCCcounts(nbIdx)
                    % all of the boundary slice's connected components are
                    % connected to one of the neighbor's connected components
                    propNbrs = false(size(isSglRev));
                    propNbrs(nbrs) = true;
                    isNxtBnd = isNxtBnd | (propNbrs & ~isSglRev);
                    isSglRev(bIdx) = true;
                    fprintf('sgl-rev propagated to %d\n', bIdx);
                    changeMade = true;
                    break;
                else
                    % slice may still be sgl-rev if other side was falsely
                    % marked as split, and hasn't been repaired yet
                    nbImg = slcImgs(:, :, nbIdx);
                    notConnToSglRev = cellfun(@(x)(nnz(nbImg(x)) == 0), nbrCCs.PixelIdxList);
                    disconnPts = vertcat(nbrCCs.PixelIdxList{notConnToSglRev});
                    disconnPts = imgIdxToPtIdx(disconnPts);
                    if isempty(disconnPts)
%                         figure; imshow(img);
                        continue
                    end
%                     assert(~isempty(disconnPts));
                    assert(nnz(disconnPts == 0) == 0);
                    assert(length(nbrs) == 2);

                    % move the slice boundary so that the other (non-sgl-rev)
                    % slice gets the points we couldn't connect to the sgl-rev
                    % slice (move all within the boundary to avoid falsely
                    % propagating the sgl-rev region)
                    otherNbrIdx = nbrs(3-jj);
                    assert(ismember(abs(bIdx - otherNbrIdx), [1 numSlices-1]));
%                     assert(~isSglRev(otherNbrIdx));
                    if isSglRev(otherNbrIdx)
                        % can't propagate to other neighbor since it's also
                        % sgl-rev
                        continue;
                    end
    %                 if (bIdx == 1 && otherNbIdx == numSlices) || ...
    %                         (bIdx > 1 && (bIdx == otherNbrIdx+1))
                    if bIdx == mod(otherNbrIdx+1 - 1, numSlices) + 1
                        mvToOther = (pt2slice == bIdx) & (thVals <= max(thVals(disconnPts)));
                    elseif bIdx == mod(otherNbrIdx-1 - 1, numSlices) + 1
                        mvToOther = (pt2slice == bIdx) & (thVals >= min(thVals(disconnPts)));
                    else
                        assert(false);
                    end
                    
                    % the slice may have been given extra theta-values that
                    % wrap around 2*pi
                    if (bIdx == 1) && (otherNbrIdx == 2)
                        mvToOther = mvToOther & (thVals < pi);
                    elseif (bIdx == numSlices) && (otherNbrIdx == numSlices-1)
                        mvToOther = mvToOther & (thVals > pi);
                    end
                    
%                     assert(sum(mvToOther) > 0)
                    if sum(mvToOther == 0)
                        break;
                    end

                    if sum((pt2slice == bIdx) & ~mvToOther) > 0
    %                     plotSliceImg(); title(sprintf('slices before adjustment (bIdx = %d, otherNbrIdx = %d)', bIdx, otherNbrIdx));
                        pt2slice(mvToOther) = otherNbrIdx;
                        updateSliceInfo(otherNbrIdx);
                        updateSliceInfo(bIdx);
    %                     plotSliceImg(); title(sprintf('slices after adjustment (bIdx = %d, otherNbrIdx = %d)', bIdx, otherNbrIdx));
                        retry = true;

    %                     slcChgImg = false(size(img));
    %                     slcChgImg(img > 0) = mvToOther;
    %                     figure; imshow(slcChgImg); title('slice reassignments');
                    else
                        warning('single-sided propagation stopped to avoid emptying slice');
                    end
                    
%                     if sum(pt2slice == bIdx) == 0
%                         % TODO: handle this case if needed
%                         % will this happen whenever the slice is a true
%                         % split?  if so, detect that slice would be empty
%                         % before making the change, and skip this iteration
%                         error('empty slice');
%                     end
                end
            end
        end
    end

%     isSglRev'
    srBnd = find(isNxtBnd);
%     srBnd'
end
fprintf('done performing single-sided propagation\n');
end % doSglSidedProp

function changeMade = doDblSidedProp()
% check contiguous split regions to make sure they should really be split
% identify split regions
srDiff = diff(~isSglRev([end 1:end 1]) & ~emptySlc([end 1:end 1]));
srDiff = srDiff(2:end);
% srDiff'
startPts = find(srDiff == 1) + 1;
endPts = find(srDiff == -1);
if ~isempty(startPts) && startPts(1) > endPts(1)
    % one of the split regions wraps around
    endPts = endPts([2:end 1]);
end
% [startPts endPts]

changeMade = false;
for ii=1:1:length(startPts)
    sIdx = startPts(ii);
    eIdx = endPts(ii);
    fprintf('checking split region: %s\n', mat2str([sIdx eIdx]));
    if sIdx > eIdx
%         chkIdxs = mod([eIdx:1:eIdx+((sIdx+numSlices)-eIdx+1)] - 1, numSlices) + 1;
        chkIdxs = setdiff(1:1:numSlices, eIdx+1:1:sIdx-1);
    else
        chkIdxs = sIdx:1:eIdx;
    end
%     chkIdxs
    nbrIdxs = mod([sIdx eIdx] + [-1 1] - 1, numSlices) + 1;

    % %     % if neighboring slice isn't connected to split slice, don't count it
    % %     % as a neighbor
    %     if diff(nbrIdxs) == 3 || 
    %         fprintf('possible disconnected neighbor\n');
    %         splitImg = nbrImg + sum(slcImgs(:, :, chkIdxs), 3);
    %         splitCCs = bwconncomp(splitImg); splitCCs = splitCCs.NumObjects;
    %         nbrImg1 = slcImgs(:, :, nbrIdxs(1));
    %         nbr1ccs = bwconncomp(nbrImg1); nbr1ccs = nbr1ccs.NumObjects;
    %         nbrImg2 = slcImgs(:, :, nbrIdxs(2));
    %         nbr2ccs = bwconncomp(nbrImg2); nbr2ccs = nbr2ccs.NumObjects;
    %     end

    % if there are no points in the neighboring slice, don't count it as a
    % neighbor
    if slcCCcounts(nbrIdxs(1)) == 0
        nbrIdxs(1) = 0;
    end
    if slcCCcounts(nbrIdxs(2)) == 0
        nbrIdxs(2) = 0;
    end

    nbrIdxs = nbrIdxs(nbrIdxs > 0);
    if isempty(nbrIdxs)
        fprintf('\tkeeping split (no sgl-rev nbrs)\n');
        continue;
    elseif length(nbrIdxs) == 1
        % split region is actually joined iff everything in the split
        % region is connected to the single-rev region
        chkImg = sum(slcImgs(:, :, [chkIdxs nbrIdxs]), 3);
        chkCCs = bwconncomp(chkImg);
        chkCCs = chkCCs.NumObjects;
        if chkCCs == slcCCcounts(nbrIdxs)
            isSglRev(chkIdxs) = true;
            fprintf('\t making sgl-rev (connected to one-sided nbr)\n');
            changeMade = true;
        else
            fprintf('\t keeping split (not connected to one-sided nbr)\n');
        end
    else
        % split region is actually joined iff the split region connects
        % some from points from both neighbor slices
        nbrImg = sum(slcImgs(:, :, nbrIdxs), 3);
        if ismember(abs(diff(nbrIdxs)), [0 1 numSlices-1])
            warning('neighbor slices of split region are contiguous\n');
%             plotFlag = true;  % TODO: put a flag in the output instead
            % neighbors would already be connected; we need to fix this
            switch abs(diff(nbrIdxs))
            case 15
                excludeIdxs = inds(thVals > (2*pi - sliceAng/2) | ...
                thVals < sliceAng/2);
            case 0
                slcNum = nbrIdxs(1);
                excludeIdxs = inds(...
                    thVals >= ((slcNum-1)*sliceAng + sliceAng/4) & ...
                    thVals < (slcNum*sliceAng - sliceAng/4));   
            case 1
                excludeIdxs = inds(...
                    thVals >= ((slcNum-1)*sliceAng + sliceAng/4) & ...
                    thVals < (slcNum*sliceAng - sliceAng/4));
            otherwise
                assert false;
            end
            nbrImg(excludeIdxs) = false;
            figure; imshow(nbrImg); title('forced split of contiguous neighbors');
        end
        chkImg = nbrImg + sum(slcImgs(:, :, chkIdxs), 3);
        chkCCs = bwconncomp(chkImg);
        chkCCs = chkCCs.PixelIdxList;
        joined = false;
        for ccIdx = 1:1:length(chkCCs)
            ccSlices = unique(pt2slice(ismember(inds, chkCCs{ccIdx})));
            if sum(ismember(nbrIdxs, ccSlices)) == 2
                joined = true;
                break;
            end
        end
        if joined 
            isSglRev(chkIdxs) = true;
            fprintf('\t making sgl-rev (connects nbrs)\n');
            changeMade = true;
        else
            fprintf('\t keeping split (does not connect nbrs)\n');
        end
    end
end
end

function isInr = createForcedSplit()
    thv = mod(thVals, 2*pi); 
    ths = sort(thv, 'ascend');
    gaps = diff(ths);
    zGap = ths(1) + 2*pi - ths(end);
    % find the theta-value in the middle of the largest gap, and shift
    % it by pi to make an even (in theta-range) split of the 
    % cluster points
    [mgV, mgI] = max(gaps);
    if zGap > mgV
        % gap is through zero
        thBorder = mod((ths(1) + (-ths(end)) / 2) + pi, 2*pi);
    else
        thBorder = mod( ( (ths(mgI) + ths(mgI+1)) / 2) + pi, 2*pi);
    end
    thvr = mod(thv - thBorder, 2*pi);
    ltb = (thvr >= 0) & (thvr < pi);
    if mean(rVals(ltb)) < mean(rVals(~ltb))
        isInr = false(size(rVals));
        isInr(ltb) = true;
    else
        isInr = true(size(rVals));
        isInr(ltb) = false;
    end
    % TODO: instead of forcing isInner to conform to slice boundaries,
    % modify removeThetaDiscontFor2rev to remove this assumption
    thBorderSlice = floor(thBorder / sliceAng);
    inBorderSlice = (pt2slice == thBorderSlice);
    if nnz(inBorderSlice) > 0
        if sum(isInner(inBorderSlice)) > sum(~isInner(inBorderSlice))
            isInner(inBorderSlice) = true;
        else
            isInner(inBorderSlice) = false;
        end
    end
end

% for debugging
function plotSliceImg()
    slcImg = -ones(size(img));
    slcImg(img > 0) = pt2slice;
    figure; imagesc(slcImg); axis image
end

% for debugging
function plotSlices()
    [xVals, yVals] = pol2cart(thVals, rVals);
    figure; hold on; axis equal;
    
    scatter(xVals(isSglRev(pt2slice)), yVals(isSglRev(pt2slice)), 'g.');
    scatter(xVals(~isSglRev(pt2slice)), yVals(~isSglRev(pt2slice)), 'm.');
    for slcNum = 1:1:numSlices
        polar((slcNum-1)*sliceAng + [0 eps], [0 max(rVals)], '-k');
    end
end

function plotIsInner()
    iiImg = im2double(img > 0);
    iiImg(iiImg > 0) = (isInner*2 - 1);
    figure; imagesc(iiImg); axis image;
end

% function slcDes = determineForSlice(inSlice)
%     slcR = rVals(inSlice);
%     nInSlice = nnz(inSlice);
%     if nInSlice == 0
%         slcDes = [];
%         sliceIsEither(slcNum) = true;
%         return;
%     elseif nInSlice == 1
%         slcDes = zeros(size(slcR));
%         sliceIsEither(slcNum) = 1;
%         return;
%     end
%     
% %     slcImg = (img > 0);
% %     slcImg(inds(~inSlice)) = 0;
% %     figure; imshow(slcImg);
% 
%     % if there is a large gap in r-values, we consider the large
%     % r-values to be part of the inner spiral and small r-values to be
%     % part of the outer spiral
%     [srtR, rSrtInds] = sort(slcR, 'ascend');
%     rGaps = diff(srtR);
% 
%     [srtGaps, gapSrtInds] = sort(rGaps, 'ascend');
%     interGapSize = srtGaps(end); 
% 
%     % if there is an inner and outer spiral, the points will be divided
%     % according to their r-values, along the largest gap in these
%     % r-values.  We calculate this potential split here.
%     gapSortPosn = gapSrtInds(end);
%     isPossibleInner = true(size(slcR));
%     isPossibleInner(rSrtInds(gapSortPosn+1:end)) = false;
% 
%     % now we see if this split is big enough to consider it to be a
%     % split between inner and outer parts of the spiral
%     preGapR = slcR(isPossibleInner);
%     postGapR = slcR(~isPossibleInner);
% 
%     intraWidth = mean([ max(preGapR) - min(preGapR), ...
%         max(postGapR) - min(postGapR) ]);
% 
%     if (interGapSize / intraWidth) > gapRatioThres
%         % there is an inner and outer spiral; accept the split
%         slcDes(isPossibleInner) = -1;
%         slcDes(~isPossibleInner) = 1;
%         sliceWasSplit(slcNum) = true;
%     else
%         % no large gap; everything is part of the inner spiral or outer
%         % spiral
%         slcDes = zeros(size(isPossibleInner));
%         sliceIsEither(slcNum) = true;
%     end
%     
%     if plotFlag
%         fprintf('\tinter-gap, intra-width: %2.4f, %2.4f\n', ...
%             interGapSize, intraWidth);
%         fprintf('\tratio: %2.4f\n', interGapSize/intraWidth);
%         if sliceIsEither(slcNum)
%             fprintf('\tnot split\n');
%         else
%             fprintf('\tsplit\n');
%         end
%     end
% end % function determineForSlice

end