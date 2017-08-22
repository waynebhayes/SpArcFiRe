function [scTbl, gxyTbl, arcsTbl] = genClusScoreTbl(inFolder, gxyNames, distVals, ctrR, ctrC, stgs)
% Not part of main flow - calculates cluster (contrast) scores based on
% cluster mask images read from disk
% OUTPUTS:
%   scTbl: one row per galaxy in gxyNames where the images could be read,
%       with the first column containing the names and the subsequent 
%       columns containing arc-length-ordered cluster-contrast scores, one
%       column for each distance value (used to determine the range of
%       outside pixels)
%   gxyTbl: cluster scores suitable for joining with a galaxy-level csv
%   arcsTbl: sorted arc scores suitable for joining with an _arcs.csv

% fnames = dir(inFolder);
% names = {fnames.name};
% names = names(cellfun(@(x)(~isempty(strfind(x, 'clusMask'))), names));

if inFolder(end) ~= filesep
    inFolder = [inFolder filesep];
end
    
scTbl = cell(length(gxyNames), length(distVals)+1);
nMissing = 0;
nTotalArcs = 0;
for ii=1:1:length(gxyNames)
    if mod(ii, 100) == 0
        fprintf('%d ', ii);
        if mod(ii, 1000) == 0
            fprintf('\n');
        end
    end
    try
        itsyImg = imread([inFolder gxyNames{ii} '-B_autoCrop.png']);
        itsyImg = im2double(itsyImg);
        maskImg = imread([inFolder gxyNames{ii} '-G_clusMask-merged.png']);
        
    catch ME
        nMissing = nMissing + 1;
        continue;
    end
    imgSz = size(itsyImg);
    scTbl{ii, 1} = gxyNames{ii};
%     colors = unique(reshape(permute(maskImg, [3 1 2]), prod(imgSz), 3), 'rows');
    pxlValList = reshape(maskImg(:), prod(imgSz), 3);
    colors = unique(pxlValList, 'rows');
    colors = setdiff(colors, [0 0 0], 'rows');
    [inClus, clusIdxs] = ismember(pxlValList, colors, 'rows');
    clusIdxs = reshape(clusIdxs, imgSz);
%     figure; imagesc(clusIdxs); axis image
    nClus = size(colors, 1);
    clusMtxs = zeros([imgSz nClus]);
    for jj=1:1:nClus
        clusMtxs(:, :, jj) = itsyImg .* (clusIdxs == jj);
    end
    [params, bounds] = fitLogSpiralsToClusters(clusMtxs, ctrR, ctrC, stgs);
    alens = calcLgspArcLengths(params, bounds);
    [temp, alenOrder] = sort(alens, 'descend');
    for jj=1:1:length(distVals)
        scores = calcClusScores(itsyImg, clusMtxs, struct('warnings', {}), distVals(jj));
        scTbl{ii, jj+1} = scores(alenOrder);
    end
    nTotalArcs = nTotalArcs + nClus;
end
fprintf('\n');

% delete empty rows
incl = ~cellfun(@isempty, scTbl(:, 1));
scTbl = scTbl(incl, :);
inclGxyNames = gxyNames(incl);

gxyTbl = cell(size(scTbl, 1) + 1, size(scTbl, 2));
gxyTbl{1, 1} = 'gxyName';
arcsTbl = cell(nTotalArcs + 1, length(distVals) + 2);
arcsTbl(1, 1:2) = {'gxyName', 'alenRank'};

for ii=1:1:length(distVals)
    colName = sprintf('medianClusContrastScore_dist%d', distVals(ii));
    gxyTbl{1, ii+1} = colName;
    arcsTbl{1, ii+2} = colName;
end

gxyTbl(2:end, 1) = scTbl(:, 1);
gxyTbl(2:end, 2:end) = cellfun(@median, scTbl(:, 2:end), 'UniformOutput', false);

curInd = 2;
for ii=1:1:length(inclGxyNames)
    numArcs = length(scTbl{ii, 2});
    curGxyInds = curInd:(curInd+numArcs-1);
    arcsTbl(curGxyInds, 1) = inclGxyNames(ii);
    arcsTbl(curGxyInds, 2) = mat2cell(1:numArcs, 1, ones(1, numArcs));
    for jj=1:1:length(distVals)
        arcsTbl(curGxyInds, jj+2) = mat2cell(scTbl{ii, jj+1}, ones(1, numArcs), 1);
    end
    curInd = curInd + numArcs;
end

fprintf('%d images were skipped (not found)\n', nMissing);

end