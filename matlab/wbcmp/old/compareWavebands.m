function [clusMasks, b1ImgS, b2ImgS] = compareWavebands(b1path, b2path, stgs, name)

if nargin < 4
    name = [];
end

fitUsingNonUsmIVals = stgs.fitUsingNonUsmIVals;
clusSizeCutoff = stgs.clusSizeCutoff;

if ~fitUsingNonUsmIVals
    error('waveband comparisons must use original intensity values');
end

wtScalingAmt = 1;

b1ImgF = im2double(imread(b1path));
b2ImgF = im2double(imread(b2path));

b1ImgN = convertFromFits(b1ImgF);
b2ImgN = convertFromFits(b2ImgF);
aImgN = (b1ImgN + b2ImgN) / 2;

[aImgS, t1, aImgSnoUsm, t2, fitParams] = ...
    preprocessImage(aImgN, stgs);

ctrR = size(aImgS, 1) / 2; ctrC = size(aImgS, 2) / 2;

clusters = findClusters(aImgS, [], stgs);

figure; imshow(aImgSnoUsm); title('standardized input to clustering');
if ~isempty(name)
    export_fig([name '_input.png'], gcf);
end

clusMtxs = hac2mtxs(clusters, aImgS, clusSizeCutoff);
showClustersFromMtxs(clusMtxs, size(aImgS));

[barInds, barAngles, barHalfLens, clusMtxs] = ...
    findBarFromClusters(clusMtxs, ctrR, ctrC);
showClustersFromMtxs(clusMtxs, size(aImgS));

if stgs.useExtendedMerging
    clusMtxs = mergeClusters(clusMtxs, ctrR, ctrC, stgs);
else
    clusMtxs = mergeClustersByFit(clusMtxs, ctrR, ctrC, stgs);
end

[barInds, barAngles, barHalfLens, clusMtxs] = ...
    findBarFromClusters(clusMtxs, ctrR, ctrC);
showClustersFromMtxs(clusMtxs, size(aImgS));
if ~isempty(name)
    export_fig([name '_clusters.png'], gcf);
end

clusMasks = clusMtxs > 0;
inClus = sum(clusMasks, 3) > 0;

fitsPprocStgs = stgs;
fitsPprocStgs.useElpsCutoff = false;
fitsPprocStgs.useUnsharpMask = false;
b1ImgS = preprocessImage(b1ImgF, fitsPprocStgs, fitParams);
b2ImgS = preprocessImage(b2ImgF, fitsPprocStgs, fitParams);

figure; imagesc(b1ImgS, quantile(b1ImgS(:), [0.2, 0.999])); axis image; title('waveband 1'); impixelinfo; colormap(gray(256));
if ~isempty(name)
    export_fig([name '_waveband1.png'], gcf);
end
figure; imagesc(b2ImgS, quantile(b2ImgS(:), [0.2, 0.999])); axis image; title('waveband 2'); impixelinfo; colormap(gray(256));
if ~isempty(name)
    export_fig([name '_waveband2.png'], gcf);
end

% curImg = inClus .* (b1ImgS.^2); figure; imagesc(curImg, quantile(curImg(:), [0.2, 0.999])); axis image
% curImg = inClus .* (b2ImgS.^2); figure; imagesc(curImg, quantile(curImg(:), [0.2, 0.999])); axis image
% curImg = inClus .* exp(b1ImgS); figure; imagesc(curImg, quantile(curImg(:), [0.2, 0.999])); axis image
% curImg = inClus .* exp(b1ImgS); figure; imagesc(curImg, quantile(curImg(:), [0.2, 0.999])); axis image
% b1n = b1ImgS / median(b1ImgS(inClus)); b2n = b2ImgS / median(b2ImgS(inClus));
% curImg = inClus .* max(b1n - b2n, 0); figure; imagesc(curImg); axis image
% curImg = inClus .* max(b2n - b1n, 0); figure; imagesc(curImg); axis image

for ii=1:1:size(clusMasks, 3)
    curMask = clusMasks(:, :, ii);
%     nClus = nnz(curMask);
    b1clus = curMask .* b1ImgS;
    b2clus = curMask .* b2ImgS;
%     b1clus = b1clus * nClus / sum(b1clus(:)) * wtScalingAmt;
%     b2clus = b2clus * nClus / sum(b2clus(:)) * wtScalingAmt;
    b1clus = b1clus / median(b1clus(curMask)) * wtScalingAmt;
    b2clus = b2clus / median(b2clus(curMask)) * wtScalingAmt;
%     figure; imagesc(b1clus); axis image;
%     figure; imagesc(b2clus); axis image;

    plotClusBrt = 1;
    
    fimg1 = b1clus;
    fimg2 = b2clus;
%     subplot(1, 2, 1); imagesc(fimg1); axis image
%     subplot(1, 2, 2); imagesc(fimg2); axis image
    [params1, bounds1] = fitLogSpiral(fimg1, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    [params2, bounds2] = fitLogSpiral(fimg2, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    fh = figure('Position', [0 0 (size(fimg1, 2) .* 8) (size(fimg1, 1) .* 1.5)]);
    title(sprintf('cluster %d\n%s\n%s', ii, mat2str(params1, 6), mat2str(params2, 6)));
    subplot(1, 3, 1); displayLgspPlot(params1, bounds1, scaleForDisp(fimg1), ctrR, ctrC, [], [], {'m'}, {}, 1, false); axis image
    subplot(1, 3, 2); displayLgspPlot(params2, bounds2, scaleForDisp(fimg2), ctrR, ctrC, [], [], {'c'}, {}, 1, false); axis image
    subplot(1, 3, 3); displayLgspPlot([params1; params2], [bounds1; bounds2], plotClusBrt * clusMasks(:, :, ii), ctrR, ctrC, [], [], {'m', 'c'}, {}, 1, false); 
    colormap(gray(256));
    axis off
    if ~isempty(name)
%         saveas(fh, [name sprintf('_%03d-A.png', ii)]);
%         close(fh)
        export_fig([name sprintf('_%03d-A.png', ii)], '-a2', fh);
        close(fh);
    end

    fimg1 = max(b1clus - b2clus, eps) .* curMask;
    fimg2 = max(b2clus - b1clus, eps) .* curMask;
%     subplot(1, 2, 1); imagesc(fimg1); axis image
%     subplot(1, 2, 2); imagesc(fimg2); axis image
    [params1, bounds1] = fitLogSpiral(fimg1, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    [params2, bounds2] = fitLogSpiral(fimg2, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    fh = figure('Position', [0 0 (size(fimg1, 2) .* 8) (size(fimg1, 1) .* 1.5)]);
    title(sprintf('cluster %d\n%s\n%s', ii, mat2str(params1, 6), mat2str(params2, 6)));
    subplot(1, 3, 1); displayLgspPlot(params1, bounds1, scaleForDisp(fimg1), ctrR, ctrC, [], [], {'m'}, {}, 1, false); axis image
    subplot(1, 3, 2); displayLgspPlot(params2, bounds2, scaleForDisp(fimg2), ctrR, ctrC, [], [], {'c'}, {}, 1, false); axis image
    subplot(1, 3, 3); displayLgspPlot([params1; params2], [bounds1; bounds2], plotClusBrt * clusMasks(:, :, ii), ctrR, ctrC, [], [], {'m', 'c'}, {}, 1, false); 
    colormap(gray(256));
    axis off
    if ~isempty(name)
%         saveas(fh, [name sprintf('_%03d-A.png', ii)]);
%         close(fh)
        export_fig([name sprintf('_%03d-B.png', ii)], '-a2', fh);
        close(fh);
    end
    
    fimg1 = double(curMask); fimg2 = double(curMask);
    fimg1(curMask) = b1clus(curMask)./(b1clus(curMask)+b2clus(curMask));
    fimg2(curMask) = b2clus(curMask)./(b1clus(curMask)+b2clus(curMask));
% %     fimg1(curMask) = fimg1(curMask) - min(fimg1(curMask)) + eps;
% %     fimg2(curMask) = fimg2(curMask) - min(fimg2(curMask)) + eps;
%     subplot(1, 2, 1); imagesc(fimg1); axis image
%     subplot(1, 2, 2); imagesc(fimg2); axis image
    [params1, bounds1] = fitLogSpiral(fimg1, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    [params2, bounds2] = fitLogSpiral(fimg2, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    fh = figure('Position', [0 0 (size(fimg1, 2) .* 8) (size(fimg1, 1) .* 1.5)]);
    title(sprintf('cluster %d\n%s\n%s', ii, mat2str(params1, 6), mat2str(params2, 6)));
    subplot(1, 3, 1); displayLgspPlot(params1, bounds1, scaleForDisp(fimg1), ctrR, ctrC, [], [], {'m'}, {}, 1, false); axis image
    subplot(1, 3, 2); displayLgspPlot(params2, bounds2, scaleForDisp(fimg2), ctrR, ctrC, [], [], {'c'}, {}, 1, false); axis image
    subplot(1, 3, 3); displayLgspPlot([params1; params2], [bounds1; bounds2], plotClusBrt * clusMasks(:, :, ii), ctrR, ctrC, [], [], {'m', 'c'}, {}, 1, false); 
    colormap(gray(256));
    axis off
    if ~isempty(name)
%         saveas(fh, [name sprintf('_%03d-B.png', ii)]);
%         close(fh)
        export_fig([name sprintf('_%03d-C.png', ii)], '-a2', fh);
        close(fh);
    end

    fimg1 = double(curMask); fimg2 = double(curMask);
    fimg1(curMask) = b1clus(curMask)./b2clus(curMask);
    fimg2(curMask) = b2clus(curMask)./b1clus(curMask);
%     subplot(1, 2, 1); imagesc(fimg1); axis image
%     subplot(1, 2, 2); imagesc(fimg2); axis image
    [params1, bounds1] = fitLogSpiral(fimg1, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    [params2, bounds2] = fitLogSpiral(fimg2, ctrR, ctrC, stgs, [], [], [], @wbCmpLgspErr);
    fh = figure('Position', [0 0 (size(fimg1, 2) .* 8) (size(fimg1, 1) .* 1.5)]);
    title(sprintf('cluster %d\n%s\n%s', ii, mat2str(params1, 6), mat2str(params2, 6)));
    subplot(1, 3, 1); displayLgspPlot(params1, bounds1, scaleForDisp(fimg1), ctrR, ctrC, [], [], {'m'}, {}, 1, false); axis image
    subplot(1, 3, 2); displayLgspPlot(params2, bounds2, scaleForDisp(fimg2), ctrR, ctrC, [], [], {'c'}, {}, 1, false); axis image
    subplot(1, 3, 3); displayLgspPlot([params1; params2], [bounds1; bounds2], plotClusBrt * clusMasks(:, :, ii), ctrR, ctrC, [], [], {'m', 'c'}, {}, 1, false); 
    colormap(gray(256));
    if ~isempty(name)
        export_fig([name sprintf('_%03d-D.png', ii)], '-a2', fh);
        close(fh);
    end
end
end