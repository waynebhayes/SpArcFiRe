function [clusMtxs1r, clusMtxs2r] = ...
    matchClusters(clusMtxs1, clusMtxs2, ctrR, ctrC, stgs, name, resPath)

if nargin >= 7 && ~isempty(resPath) && isempty(name)
    error('if resPath specified, name must also be specified');
end
if nargin >= 7 && ~isempty(resPath)
    resPath = [resPath filesep];
end
if nargin < 7
    resPath = [];
end

imgSz = size(clusMtxs1(:, :, 1));

nClus1 = size(clusMtxs1, 3);
nClus2 = size(clusMtxs2, 3);
otherClusts = (clusMtxs1 > 0);
closestClus = zeros(1, nClus2);
for ii=1:1:nClus2
    curClus = clusMtxs2(:, :, ii);
    overlaps = squeeze(sum(sum(repmat(curClus, [1 1 nClus1]) & otherClusts, 1), 2));
    [sV, sI] = sort(overlaps, 'descend');
%     figure;
%     for jj=1:1:3
%         subplot(1, 3, jj);
%         if sV(jj) > 0
%             ovlImg = makeOverlapImg(curClus, clusMtxs1(:, :, sI(jj)));
%             title(sprintf('overlap %d', sV(jj)));
%             imshow(ovlImg);
%         else
%             imshow(zeros(imgSz));
%         end
%     end
    if sV(1) > 0
        closestClus(ii) = sI(1);
    end
end

allClus1 = sum(clusMtxs1, 3) > 0;
allClus2 = sum(clusMtxs2, 3) > 0;
inBothClus = false(imgSz);
c1only = false(imgSz);
c2only = false(imgSz);
c1mismatch = false(imgSz);
c2mismatch = false(imgSz);
matchedCM2 = zeros(size(clusMtxs1));
for ii=1:1:nClus1
    curClus = clusMtxs1(:, :, ii);
    isMatch = (closestClus == ii);
    if sum(isMatch) > 0
        matchingClusts = sum(clusMtxs2(:, :, isMatch), 3);
    else
        matchingClusts = zeros(imgSz);
    end
    matchedCM2(:, :, ii) = matchingClusts;
%     figure; imshow(makeOverlapImg(curClus, matchingClusts));
    
    inBothClus = inBothClus | (curClus & matchingClusts);
    c1only = c1only | (curClus & ~matchingClusts);
    c2only = c2only | (~curClus & matchingClusts);
    c1mismatch = c1mismatch | (curClus & allClus2 & ~matchingClusts);
    c2mismatch = c2mismatch | (matchingClusts & allClus1 & ~curClus);
end
assert(sum(c1mismatch(:) == c2mismatch(:)) == numel(c1mismatch));
mismatch = c1mismatch;
noMatch = sum(clusMtxs2(:, :, (closestClus == 0)), 3) > 0;

mtchImg = im2double(repmat(inBothClus, [1 1 3]));
% mtchImg = mtchImg + repmat(0.5 * noMatch, [1 1 3]);
mtchImg(:, :, 1) = mtchImg(:, :, 1) + (c1only & ~mismatch);
mtchImg(:, :, 3) = mtchImg(:, :, 3) + (c2only & ~mismatch);
mtchImg(:, :, 1) = mtchImg(:, :, 1) + mismatch;
mtchImg(:, :, 2) = mtchImg(:, :, 2) + mismatch;
mtchImg(:, :, 2) = mtchImg(:, :, 2) + noMatch;

[clusMtxs1r, clusMtxs2r] = reassignMismatches(clusMtxs1, matchedCM2, ctrR, ctrC, stgs);

cImg1 = showClustersFromMtxs(clusMtxs1, imgSz);
cImg2 = showClustersFromMtxs(clusMtxs2, imgSz);
ccImg = showClustersFromMtxs(matchedCM2, imgSz, [], -1);
cImg1r = showClustersFromMtxs(clusMtxs1r, imgSz, [], -1);
cImg2r = showClustersFromMtxs(clusMtxs2r, imgSz, [], -1);
figure;
subplot(2, 3, 1); imshow(mtchImg); title('cluster match');
subplot(2, 3, 4); imshow(cImg2); title('band 2');
subplot(2, 3, 2); imshow(cImg1); title('band 1');
subplot(2, 3, 5); imshow(ccImg); title('band 2 (matched to band 1)');
subplot(2, 3, 3); imshow(cImg1r); title('band 1 (mismatches reassigned)');
subplot(2, 3, 6); imshow(cImg2r); title('band 2 (mismatches reassigned)');

if ~isempty(resPath)
    saveas(gcf, [resPath name '_clusMatch.png']);
end

showClustersFromMtxs(clusMtxs1r, imgSz, [], -1);
showClustersFromMtxs(clusMtxs2r, imgSz, [], -1);

    function ovlImg = makeOverlapImg(clus1, clus2)
        ovlImg = zeros([size(clus1) 3]);
        ovlImg(:, :, 1) = clus1;
        ovlImg(:, :, 3) = clus2;
    end
end