imFolderPath = 'D:\Galaxy-Sets\Nair-100K\14K\100kpcrotate\';

[imgPaths, imgNames] = loadImageFolder(imFolderPath, 'jpg');
imgPaths = sortrows(imgPaths);

LoadNairAbrahamCatalog

idPos = strmatch('fpCid', NairAbrahamCatalog(1, :), 'exact');
NairAbrahamCatalog(2:end, :) = ...
    sortrows(NairAbrahamCatalog(2:end, :), idPos);
% sorted to align image path indices with catalog indices, and thus sample
% indices

% criteria for sample inclusion
ttypePos = strmatch('TType', NairAbrahamCatalog(1, :), 'exact');
pairPos = strmatch('Pairs', NairAbrahamCatalog(1, :), 'exact');
axsRatioPos = strmatch('bOverA', NairAbrahamCatalog(1, :), 'exact');

sampleSize = inf;
fprintf('sampling...');
[sampleIdxs, eligIdxs] = getImgSampleIdxs(NairAbrahamCatalog, ...
    {ttypePos, @(type)(type >= 1 && type <= 7);...
    pairPos, @(pair)(pair == 0);...
    axsRatioPos, @(axsRatio)(axsRatio > 0.3)},...
    sampleSize);
fprintf('done sampling\n');
fprintf('number of images eligible for sample: %d\n', length(eligIdxs));
fprintf('sample size: %d\n', sampleSize);

NairAbrahamSamplePool = NairAbrahamCatalog([1 1+eligIdxs'], :);
NairAbrahamCatalogSample = NairAbrahamCatalog([1 1+sampleIdxs'], :);
sampledImgPaths = imgPaths(sampleIdxs);
sampledImgNames = imgNames(sampleIdxs);

clear imFolderPath imgPaths idPos ttypePos pairPos axsRatioPos

fprintf('saving sample...');
save('NairAbrahamSample', ...
    'sampledImgPaths', 'sampledImgNames', 'NairAbrahamCatalogSample');
writeTableToFile(NairAbrahamCatalogSample, 'NairAbrahamSample.tbl', '\t');
fprintf('done saving sample\n');