function [clusMtxList] = precomputeClusters(gxyBandPaths, outDir, stgs)
% gxyBandPaths should be a list of fits images for the same galaxy and
% resolution and different bands

stgs
length(gxyBandPaths)

convStgs = struct('blq', 0.5, 'wlq', 0.999, 'asinhBeta', 2);

skipped = false(1, length(gxyBandPaths));
imgs = cell(1, length(gxyBandPaths));
fimgs = cell(1, length(gxyBandPaths));
for ii=1:1:length(gxyBandPaths)
    curPath = gxyBandPaths{ii};
    if ~exist(curPath, 'file')
        fprintf('%s not found.  Skipping.\n', curPath);
        skipped(ii) = true;
        continue;
    end
    if isempty(regexpi(curPath, '\.fits', 'end'))
        fprintf('%s not a fits file. Skipping.\n', curPath);
        skipped(ii) = true;
        continue;
    end
    fprintf('loading %s ...\n', curPath);
    fimgs{ii} = im2double(imread(curPath));
    imgs{ii} = convertFromFits(fimgs{ii}, true, convStgs);
end

fprintf('combining (for standardization only) ... \n');
imgs = imgs(~skipped);
sumImg = imgs{1};
imgDims = size(imgs{1});
for ii=2:1:length(imgs)
    if sum(imgDims == size(imgs{ii})) ~= 2
        error('image resolutions should be the same');
    end
    sumImg = sumImg + imgs{ii};
end
meanImg = sumImg ./ length(imgs);

fprintf('standardizing image ... \n');
[meanImgS, t2, t3, t4, stdzParams] = ...
    preprocessImage(meanImg, stgs);

% fprintf('standardizing on first image ...\n');
% [t1, t2, t3, t4, stdzParams] = ...
%     preprocessImage(imgs{1}, stgs);

% figure; imshow(meanImg); title('combined image for de-projection');

simgs = cell(size(imgs));
stgsNoUsm = stgs;
stgsNoUsm.unsharpMaskAmt = 0;
for ii=1:1:length(imgs);
    curImg = preprocessImage(fimgs{ii}, stgsNoUsm, stdzParams);
%     simgs{ii} = convertFromFits(curImg, 0.2, 0.999);
    simgs{ii} = convertFromFits(curImg, [], convStgs);
%     figure; imshow(simgs{ii}); title(sprintf('standardized image: %s', gxyBandPaths{ii}));
%     size(simgs{ii})
end

clusMtxList = cell(size(simgs));
gxyParamsList = cell(size(simgs));
barInfoList = cell(size(simgs));
barIndsList = cell(size(simgs));
for ii=1:1:length(simgs)
%     gxyParams.warnings = {};
    fprintf('performing clustering on image %d of %d ... \n', ii, length(simgs));
%     curClusts = findClusters(unsharpMask(simgs{ii}, stgs), [], stgs);
%     curClusMtxs = hac2mtxs(curClusts, simgs{ii}, stgs.clusSizeCutoff);
%     showClustersFromMtxs(curClusMtxs, size(simgs{ii}));
    gxyOutputCtrl = struct('writeImages', true, 'displayFigures', false, 'writeTxt', false);
    tokens = regexp(gxyBandPaths{ii}, '\\', 'split'); % TODO: add case for UNIX
    gxyName = tokens{end};
    [lgspParams, lgspBounds, errs, used2rev, hasBadBounds, curBarInfo, ...
            curClusMtxs, curGxyParams, imgAutoCrop, curBarInds] = ...
        findClusterArcs(imgs{ii}, stgs, gxyName, gxyOutputCtrl, outDir, stdzParams);
    
    curGxyParams.name = gxyName;
    curGxyParams.fit_state = 'OK';
    curGxyParams.cpu_time = -1;
    curGxyParams.wall_time = -1;
    
    curGxyParams = getGalaxyParams(lgspParams, lgspBounds, errs,...
        curBarInfo, curBarInds, imgAutoCrop, curClusMtxs, stgs, curGxyParams);
    
    clusMtxList{ii} = curClusMtxs;
    barInfoList{ii} = curBarInfo;
    barIndsList{ii} = curBarInds;
    gxyParamsList{ii} = curGxyParams;
end

% clusFromCombined = findClusters(meanImgS, [], stgs);
% clusMtxsFromCombined = hac2mtxs(clusFromCombined, meanImgS, stgs.clusSizeCutoff);

figure; 
nPlots = length(clusMtxList)+1;
for ii=1:1:nPlots-1
    curClusImg = showClustersFromMtxs(clusMtxList{ii}, stgs.resizeDims);
    subplot(1, nPlots, ii);
    imshow(curClusImg);
    title(sprintf('band %d', ii));
end
% combinedClusImg = showClustersFromMtxs(clusMtxsFromCombined, stgs.resizeDims);
% subplot(1, nPlots, nPlots)
% imshow(combinedClusImg); title('combined');

names = gxyBandPaths;
for ii=1:1:length(names)
    fslocs = strfind(names{ii}, filesep);
    if ~isempty(fslocs)
        curname = names{ii};
        names{ii} = curname(fslocs(end)+1:end);
    end
    fextlocs = strfind(names{ii}, '.');
    if ~isempty(fextlocs)
        curname = names{ii};
        names{ii} = curname(1:fextlocs(end)-1);
    end
end

name1 = names{1};
resName = name1
for ii=1:1:length(names)
    curName = names{ii};
    len = min(length(name1), length(curName));
    match = (name1(1:len) == curName(1:len));
    if sum(match) < length(match)
        firstDiff = find(~match, 1);
        resName = [resName '_' curName(firstDiff:end)]
    elseif length(curName) > length(name1)
        resName = [resName '_' curName(length(name1)+1:end)]
    end
end

% save('fimgs', 'fimgs');
sfimgs = cell(length(gxyBandPaths), 1);
for ii=1:1:length(sfimgs)
    sfimg = preprocessImage(fimgs{ii}, stgsNoUsm, stdzParams);
    sfimgs{ii} = sfimg;
end

% save([outDir filesep resName], 'clusMtxList', 'clusMtxsFromCombined', 'sfimgs');
save([outDir filesep resName], 'clusMtxList', 'gxyParamsList', 'barInfoList', 'barIndsList', 'sfimgs');

%     clusts{ii} = findClusters(curImg, stgs);
%     sepLocs = strfind(curPath, filesep);
%     if isempty(sepLocs)
%         fname = img

end