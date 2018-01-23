function [clusters, failed2revEncountered] = ...
    genHacTree(simlMtx, img, barInfo, ctrR, ctrC, stgs, tracePts, aviName, colorThres)
% Performs HAC clustering and returns the resulting dendogram
% INPUTS:
%   simlMtx: the similarity matrix to use for clustering
%   img: the image used to generate the similarity matrix
%   barInfo: 
%   ctrR:
%   ctrC: 
%   stgs: structure containing algorithm settings (see settings.m)
%   tracePts: optional parameter specifying points (similarity values) at
%       which the current clustering state should be saved as a frame in an
%       AVI video file (default: empty; no AVI file)
%   aviName: (optional) name of the AVI file, if such a file is to be saved
%   colorThres: minimum size for clusters to get their own color (if an AVI
%       is to be saved)
% OUTPUTS:
%   clusters: HAC dendogram, given as a structure array of root nodes for
%       each of the cluster trees.  Each node contains indices of the two
%       merged clusters, the similarity value between these clusters, and
%       all points from these clusters

stopThres = stgs.stopThres;
% useNewBarDet = stgs.useNewBarDet;
% allowArcBeyond2pi = stgs.allowArcBeyond2pi;
errRatioThres = stgs.errRatioThres;

mergeChkMinClusSz = stgs.mergeChkMinClusSz; 

error(nargchk(4, 7, nargin));

if isempty(barInfo)
    warning('barInfo is empty');
end

if nargin < 7
    tracePts = [];
%     tracePts = 1:-0.01:stopThres;
%     colorThres = 150;
end

if nargin < 8 || isempty(aviName)
    aviName = [mat2str(round(rem(now, 1) * 10^14)) '_trace.avi'];
end

failed2revEncountered = false;

% if ~useNewBarDet
%     barInfo = [];
% end

imgSz = size(img);
% ctrR = imgSz(1) / 2; ctrC = imgSz(2) / 2;
% if stgs.useSubpixelCtr
%     ctrR = ctrR + 0.5;
%     ctrC = ctrC + 0.5;
% end

% % delete self-similarities
simlMtx = simlMtx - diag(diag(simlMtx));

selPts = find(sum(simlMtx, 1));  % nzPts
fprintf('%d cluster points without nonzero similarities to other points\n', ...
    sum(img(selPts) == 0))
selPts = selPts(img(selPts) > 0);
simlMtx = simlMtx(selPts, selPts);

% numPts = size(simlMtx, 1);
numPts = size(simlMtx, 2);

% create a cluster for each point that has a nonzero similarity with at 
% least one other point
emptyClus = struct('cidx1', [], 'cidx2', [], 'mergeVal', [], 'mergedPts', []);
clusts = repmat(emptyClus, 1, numPts);
for pt=1:1:length(selPts)
    clusts(pt).mergedPts = selPts(pt);
end
clustsActv = true(1, numPts);

if ~isempty(tracePts)
    nextTracePt = tracePts(1);
    tracePtIdx = 1;
    aviTrace = avifile(aviName, 'fps', 1, 'compression', 'None', ...
        'quality', 100);
else
    nextTracePt = -inf;
end

% clusAsgns = 1:1:numPts;
% pqEmpty = PriorityQueue([], []);
% PQClusts = repmat(pqEmpty, numPts, 1);
% pqMainKeys = zeros(1, numPts);
% for ii=1:1:numPts
%     inds = find(simlMtx(:, ii));
%     PQClusts(ii) = PriorityQueue(full(simlMtx(inds, ii)), inds);
%     [key, val] = PQClusts(ii).peekMax();
%     pqMainKeys(ii) = key;
% end
% fprintf('building main priority queue...\n');
% PQMain = PriorityQueue(pqMainKeys, 1:1:numPts);
% clusLocs = PQMain.getStartLocs();

fprintf('performing clustering...\n'); tic
keepClustering = true;
arcMergeChecks = 0;
mergeStopCount = 0;

while keepClustering
    [vals, rows] = max(simlMtx, [], 1);
    [val, col] = max(vals);
    row = rows(col);
    cidx1 = row; cidx2 = col;
    
%     val = full(val);
%     [kMain, v] = PQMain.peekMax();
%     cidx1pq = v;
%     [k, v] = PQClusts(cidx1pq).takeMax();
%     cidx2pq = clusAsgns(v);
    
%     assert(kMain == k);
%     assert(val == k)
%     assert(round(abs(val - k) * 1000) == 0)
    
    clus1 = clusts(cidx1); clus2 = clusts(cidx2);
    if min(length(clus1.mergedPts), length(clus2.mergedPts)) > mergeChkMinClusSz 
        arcMergeChecks = arcMergeChecks + 1;
        mtxs = hac2mtxs([clus1 clus2], img, 0);
        [mergeErrRatio, anyFailed2rev] = ...
            calcArcMergeErr(mtxs(:, :, 1), [], [], [], ...
                mtxs(:, :, 2), [], [], [], ctrR, ctrC, barInfo, stgs);
        failed2revEncountered = failed2revEncountered | anyFailed2rev;
        if mergeErrRatio > errRatioThres
            mergeStopCount = mergeStopCount + 1;
            simlMtx(cidx1, cidx2) = 0;
            simlMtx(cidx2, cidx1) = 0;
            
%             % TODO: re-examine these once we take advantage of symmetry in
%             % the similarity matrix
%             fprintf('rejected merge between %d and %d\n', cidx1, cidx2);
%             keep1 = clusAsgns(PQClusts(cidx1).getVals()) ~= cidx2;
%             PQClusts(cidx1) = PriorityQueue.merge(PQClusts(cidx1), pqEmpty, keep1, []);
%             keep2 = clusAsgns(PQClusts(cidx2).getVals()) ~= cidx1;
%             PQClusts(cidx2) = PriorityQueue.merge(PQClusts(cidx2), pqEmpty, keep2, []);
%             if ~isEmpty(PQClusts(cidx1))
%                 PQMain.changeKeyAtLoc(PQClusts(cidx1).peekMax(), clusLocs(cidx1));
%             else
%                 PQMain.removeEltAtLoc(clusLocs(cidx1));
%             end
%             if ~isEmpty(PQClusts(cidx2))
%                 PQMain.changeKeyAtLoc(PQClusts(cidx2).peekMax(), clusLocs(cidx2));
%             else
%                 PQMain.removeEltAtLoc(clusLocs(cidx2));
%             end
            continue;
        end
    end
    
    if val >= stopThres      
        mergedClus = struct('cidx1', cidx1, 'cidx2', cidx2, ...
            'mergeVal', val, ...
            'mergedPts', [clusts(cidx1).mergedPts clusts(cidx2).mergedPts]);
        clusts(cidx1) = mergedClus;
        clusts(cidx2) = emptyClus;
        clustsActv(cidx2) = false;
        
        newSimls = max(simlMtx(:, cidx1), simlMtx(:, cidx2));
        % zero out merged similarity and self-similarity
        newSimls(cidx1) = 0; 
        
        [nzr1, nzc1, nzv1] = find(newSimls);
        [nzr2, nzc2, nzv2] = find(simlMtx(:, cidx1));
        % since all values nonnegative, the max vals in newSimls will never
        % have a zero element where at least one of the two original rows 
        % had a nonzero element, so indices nzr2 are a subset of indices
        % nzr1
        hasOldVal = ismember(nzr1, nzr2);
        chgVals = nzv1; chgVals(hasOldVal) = chgVals(hasOldVal) - nzv2;
        chgMtx = sparse([nzr1; repmat(cidx1, size(nzr1))], ...
            [repmat(cidx1, size(nzr1)); nzr1], repmat(chgVals, 2, 1), ...
            numPts, numPts);
        simlMtx = simlMtx + chgMtx;
        
%         simlMtx(:, cidx1) = newSimls;
%         simlMtx(cidx1, :) = newSimls;
%         simlMtx(:, cidx2) = 0;
%         simlMtx(cidx2, :) = 0;

        [nzr, nzc, nzv] = find(simlMtx(:, cidx2));
        % diagonal elements are zero, so we don't need to worry about
        % subtracting them twice
        chgMtx = sparse([nzr; repmat(cidx2, size(nzr))], ...
            [repmat(cidx2, size(nzr)); nzr], repmat(-nzv, 2, 1),...
            numPts, numPts);     
        simlMtx = simlMtx + chgMtx;
        
%         clusAsgns(clusAsgns == cidx2) = cidx1;
%         keep1 = clusAsgns(PQClusts(cidx1).getVals()) ~= cidx1;
%         keep2 = clusAsgns(PQClusts(cidx2).getVals()) ~= cidx1;
%         [PQClusts(cidx1), mgdSize] = PriorityQueue.merge(PQClusts(cidx1), PQClusts(cidx2), keep1, keep2);
%         PQClusts(cidx2).makeEmpty();
%         assert(isnan(PQMain.dbgGetValAtLoc(clusLocs(cidx1))) || ...
%             clusAsgns(PQMain.dbgGetValAtLoc(clusLocs(cidx1))) == clusAsgns(cidx1));
%         if mgdSize > 0
%             % queues merge to empty if the merged cluster is isolated
%             PQMain.changeKeyAtLoc(PQClusts(cidx1).peekMax(), clusLocs(cidx1));
%         else
%             PQMain.removeEltAtLoc(clusLocs(cidx1));
%         end
%         PQMain.removeEltAtLoc(clusLocs(cidx2));      
        
        if val < nextTracePt
            nextTracePt
            infoStr = sprintf('threshold %2.4f', nextTracePt);
            clusImg = showClusters(clusts(clustsActv == 1), imgSz, colorThres);
            figure; imshow(clusImg); title(infoStr);
            aviTrace = addframe(aviTrace, gcf);
            close
            tracePtIdx = tracePtIdx + 1;
            if tracePtIdx <= length(tracePts)
                nextTracePt = tracePts(tracePtIdx);
            else
                nextTracePt = inf;
            end
        end
    else
        keepClustering = false;
    end
%     if nnz(simlMtx) == 0
%         keepClustering = false;
%     end
end
if ~isempty(tracePts)
    aviTrace = close(aviTrace);
end
toc; fprintf('...done clustering\n');
clusters = clusts(clustsActv);

fprintf('arcMergeChecks = %d\n', arcMergeChecks);
fprintf('mergeStopCount = %d\n', mergeStopCount);

%     function dbgChkPQs()
%         keysMain = PQMain.dbgGetKeys();
%         valsMain = PQMain.getVals();
%         matchBad = false(1, length(keysMain));
%         for ii=1:1:length(keysMain)
%             matchBad(ii) = (PQClusts(valsMain(ii)).peekMax() ~= keysMain(ii));
%         end
%         find(matchBad)
%         if sum(matchBad) > 0
%             fprintf('FAIL check\n');
%         else
%             fprintf('PASS check\n');
%         end
%     end
end