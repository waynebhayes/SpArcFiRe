function batchAutoCropFits(inDir, outDir, convStgs)

stgsWithStretch = getDefaultSettings();
stgsWithStretch.resizeDims = [];
stgsWithStretch.useElpsCutoff = 0;
stgsWithStretch.medFiltRad = 0;
stgsWithStretch.unsharpMaskAmt = 0;

stgsNoStretch = stgsWithStretch;
stgsNoStretch.useDeProjectStretch = false;

fnames = dir(inDir);
fnames = {fnames.name};
fnames = fnames(~cellfun(@isempty, regexp(fnames, '\.fits$')));

for ii=1:1:length(fnames)
    fprintf('processing: %s\n', fnames{ii});
    imgF = imread([inDir filesep fnames{ii}]);
%     figure; imagesc(imgF); axis image; colorbar
    img = convertFromFits(imgF, false, convStgs);
    [img, imgNoUsm, gxyParams, fitParams] = preprocessImage(img, stgsNoStretch);
    imgC = preprocessImage(imgF, stgsNoStretch, fitParams);
    imgCS = preprocessImage(imgF, stgsWithStretch, fitParams);
%     figure; imagesc(imgC); axis image; colorbar
    prev_dir = pwd();
    cd(outDir);
    fitswrite(fliplr(flipud(imgC)'), ['autoCrop_' fnames{ii}]);
    fitswrite(fliplr(flipud(imgCS)'), ['autoCrop_deproj_' fnames{ii}]);
    cd(prev_dir);
end