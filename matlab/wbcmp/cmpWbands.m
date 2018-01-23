function xcMtx = cmpWbands(clusMtxsList, sfimgs, name)

if nargin < 3
    name = [];
end

if isempty(name)
    figVis = 'on';
else
    figVis = 'off';
end

useBrt = false;
wtRadVsAzi = false;
wtDiffFromDomPa = false;
useLgspPixels = true;
selectAllClusters = true;

% figure; 
% subplot(1, 2, 1); imshow(convertFromFits(sfimgs{1}, false)); 
% subplot(1, 2, 2); imshow(convertFromFits(sfimgs{2}, false));

imgB1 = convertFromFits(sfimgs{1});
imgB2 = convertFromFits(sfimgs{2});
if ~isempty(name)
    imwrite(imgB1, [name '_001-band1.png']);
    imwrite(imgB2, [name '_002-band2.png']);
end

screen_size = get(0, 'ScreenSize');

nRho = 256;
nTheta = 720;

clus1Mtxs = clusMtxsList{1};
clus2Mtxs = clusMtxsList{2};

stgs = getDefaultSettings();
ctrR = size(clus1Mtxs, 1) / 2;
ctrC = size(clus1Mtxs, 2) / 2;
[params1, bounds1] = fitLogSpiralsToClusters(clus1Mtxs, ctrR, ctrC, stgs);
lengths1 = calcLgspArcLengths(params1, bounds1);
[params2, bounds2] = fitLogSpiralsToClusters(clus2Mtxs, ctrR, ctrC, stgs);
lengths2 = calcLgspArcLengths(params2, bounds2);
chir1 = sign(params1(:, 2))' * lengths1;
chir2 = sign(params2(:, 2))' * lengths2;
if sign(chir1) ~= sign(chir2)
    warning('detected chirality is different across the two wavebands!');
end
if sign(chir1) > 0
    fprintf('CCW winding direction detected; flipping image\n');
    clusMtxsList = {flipdim(clusMtxsList{1}, 2), flipdim(clusMtxsList{2}, 2)};
    sfimgs = {fliplr(sfimgs{1}), fliplr(sfimgs{2})};
    clus1Mtxs = clusMtxsList{1};
    clus2Mtxs = clusMtxsList{2};
end
% b1arcs = struct('params', params1, 'bounds', bounds1, 'lengths', lengths1);
% b2arcs = struct('params', params2, 'bounds', bounds2, 'lengths', lengths2);

c1Img = showClustersFromMtxs(clus1Mtxs, size(sfimgs{1}));
c2Img = showClustersFromMtxs(clus2Mtxs, size(sfimgs{2}));
ovlpImg = zeros([size(sfimgs{1}) 3]);
ovlpImg(:, :, 1) = (sum(clus1Mtxs, 3) > 0);
ovlpImg(:, :, 3) = (sum(clus2Mtxs, 3) > 0);
figure('Visible', figVis); 
subplot(1, 3, 1); imshow(c1Img); title('clusters from band 1');
subplot(1, 3, 2); imshow(c2Img); title('clusters from band 2');
subplot(1, 3, 3); imshow(ovlpImg); title('cluster overlap');
if ~isempty(name)
    saveas(gca, [name '_003-clusters.png']);
end

if selectAllClusters
    b1sel = true(size(clus1Mtxs, 3), 1);
    figure('Visible', figVis); imshow(sum(clus1Mtxs, 3) > 0);
else
    b1sel = false(size(clus1Mtxs, 3), 1);
    iptFig = figure;
    selClus = getClusSelFromUser(clus1Mtxs, b1sel, iptFig);
    while ~isempty(selClus)
        if selClus > 0
            b1sel(selClus) = ~b1sel(selClus);
        else
            b1sel = ~b1sel;
        end
        selClus = getClusSelFromUser(clus1Mtxs, b1sel, iptFig);
    end
end
if ~isempty(name)
    saveas(gca, [name '_004-b1sel.png']);
end

if selectAllClusters
    b2sel = true(size(clus2Mtxs, 3), 1);
    figure('Visible', figVis); imshow(sum(clus1Mtxs, 3) > 0);
else
    b2sel = false(size(clus2Mtxs, 3), 1);
    iptFig = figure;
    selClus = getClusSelFromUser(clus2Mtxs, b2sel, iptFig);
    while ~isempty(selClus)
        if selClus > 0
            b2sel(selClus) = ~b2sel(selClus);
        else
            b2sel = ~b2sel;
        end
        selClus = getClusSelFromUser(clus2Mtxs, b2sel, iptFig);
    end
end
if ~isempty(name)
    saveas(gca, [name '_005-b2sel.png']);
end

figure('Visible', figVis); imshow(flipdim(polartransRGB(ovlpImg, nRho, nTheta), 1)); 
title('selected clusters (polar space)');
xlabel('\theta'); ylabel('\rho');
if ~isempty(name)
    saveas(gca, [name '_-006-clusOverlapPolar.png']);
end

% figure; imshow(flipdim(polartransRGB(ovlpImg, nRho, nTheta, [], [], 'log'), 1)); 
% title('selected clusters (log-polar space)');
% xlabel('\theta'); ylabel('log(\rho)');

b1mask = sum(clus1Mtxs(:, :, b1sel), 3) > 0;
b2mask = sum(clus2Mtxs(:, :, b2sel), 3) > 0;
if useBrt
    b1Wts = (sfimgs{1} - min(min(sfimgs{1})));
    b2Wts = (sfimgs{2} - min(min(sfimgs{2})));
else
    b1Wts = ones(size(sfimgs{1}));
    b2Wts = ones(size(sfimgs{2}));
end
if wtRadVsAzi
    for ii=1:1:size(clus1Mtxs, 3)
        curClus = clus1Mtxs(:, :, ii) > 0;
        b1Wts(curClus) = b1Wts(curClus) * sin(abs(params1(ii, 2)));
    end
    for ii=1:1:size(clus2Mtxs, 3)
        curClus = clus2Mtxs(:, :, ii) > 0;
        b2Wts(curClus) = b2Wts(curClus) * sin(abs(params2(ii, 2)));
    end
end
if wtDiffFromDomPa
    domPa1 = (params1(:, 2)' * lengths1) / sum(lengths1)
    for ii=1:1:size(clus1Mtxs, 3)
        curClus = clus1Mtxs(:, :, ii) > 0;
        b1Wts(curClus) = b1Wts(curClus) * cos(abs(params1(ii, 2) - domPa1));
    end
    domPa2 = (params2(:, 2)' * lengths2) / sum(lengths2)
    for ii=1:1:size(clus2Mtxs, 3)
        curClus = clus2Mtxs(:, :, ii) > 0;
        b2Wts(curClus) = b2Wts(curClus) * cos(abs(params2(ii, 2) - domPa2));
    end
end
if useLgspPixels
    ovl1 = displayLgspOverlay(zeros(size(b1Wts)), params1, ctrR, ctrC, bounds1);
    b1Wts = b1Wts .* (sum(ovl1, 3) > 0);
    ovl2 = displayLgspOverlay(zeros(size(b2Wts)), params2, ctrR, ctrC, bounds2);
    b2Wts = b2Wts .* (sum(ovl2, 3) > 0);
end
b1clus = b1Wts .* b1mask;
b2clus = b2Wts .* b2mask;
figure('Visible', figVis); subplot(1, 2, 1); imagesc(b1clus); axis image; colormap jet; title('b1clus');
subplot(1, 2, 2); imagesc(b2clus); axis image; colormap jet; title('b2clus');
if ~isempty(name)
    saveas(gca, [name '_-014-clusWts.png']);
end

b1pol = flipud(polartransRGB(b1clus, nRho, nTheta, [], [], 'linear', 'valid'));
b2pol = flipud(polartransRGB(b2clus, nRho, nTheta, [], [], 'linear', 'valid'));

% figure; imagesc(b1pol, quantile(nonzeros(b1pol(:)), [0.05 0.95]));
% % figure; imagesc(i1pol);
% axis image; colormap gray; xlabel('\theta'); ylabel('\rho');
% figure; imagesc(b2pol, quantile(nonzeros(b2pol(:)), [0.05 0.95]));
% % figure; imagesc(i2pol);
% axis image; colormap gray; xlabel('\theta'); ylabel('\rho');

b1logpol = flipud(polartransRGB(b1clus, nRho, nTheta, [], [], 'log', 'valid'));
b2logpol = flipud(polartransRGB(b2clus, nRho, nTheta, [], [], 'log', 'valid'));

xcMtx = analyzeXCorr(b1pol, b2pol, nTheta, name);
% analyzeXCorr(b1logpol, b2logpol, nTheta, name);
% xcMtx = analyzeXCorr(b1pol > 0, b2pol > 0, nTheta, name);
% analyzeXCorr(b1logpol > 0, b2logpol > 0);

nR = size(clusMtxsList{1}, 1);
radParams = [0 0 nR/8; 0 0 nR/4; 0 0 (3*nR)/8; 0 0 nR/2];
radBounds = repmat([0 2*pi], 4, 1);
ctrR = size(clusMtxsList{1}, 1)/2;
ctrC = size(clusMtxsList{1}, 2)/2;
ol1 = displayLgspOverlay(sum(clus1Mtxs(:, :, b1sel), 3), radParams, ctrR, ctrC, radBounds);
ol2 = displayLgspOverlay(sum(clus2Mtxs(:, :, b2sel), 3), radParams, ctrR, ctrC, radBounds);
figure('Visible', figVis); subplot(1, 2, 1); imshow(ol1); subplot(1, 2, 2); imshow(ol2);

ovlpImg = zeros([size(b1mask) 3]);
ovlpImg(:, :, 1) = b1mask;
ovlpImg(:, :, 3) = b2mask;
displayLgspOverlay(ovlpImg, radParams, ctrR, ctrC, radBounds);
title('selected clusters');
if ~isempty(name)
    saveas(gca, [name '_-013-selclus.png']);
end
% set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
% set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
% saveas(gca, '009_selclus.png');

    function selClus = getClusSelFromUser(clusMtxs, selected, iptFig)
        selImg = zeros(size(sfimgs{1}));
        selImg = selImg + (sum(clusMtxs(:, :, selected), 3) > 0);
        selImg = selImg + 0.5 * (sum(clusMtxs(:, :, ~selected), 3) > 0);
        selImg = displayClusterOverlay(selImg, clusMtxs);
        if gcf ~= iptFig
            figure(iptFig);
        end
        imshow(selImg); 
        title(sprintf(...
            ['Select/deselect clusters by clicking on one of their pixels. \n'...
            'Press enter when done.']));
        selCoords = round(ginput(1));
        if isempty(selCoords)
            selClus = [];
            return
        end
        clusMatchIdx = ...
            find(squeeze(clusMtxs(selCoords(2), selCoords(1), :)) > 0);
        assert(length(clusMatchIdx) <= 1)
        if isempty(clusMatchIdx)
            selClus = 0;
        else
            selClus = clusMatchIdx;
        end
    end

end