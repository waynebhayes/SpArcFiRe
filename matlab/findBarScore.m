function [barScore, barInfo, bestRadExt] = ...
    findBarScore(img, fitParams, ctrR, ctrC, candRad, candTh, candScore, stgs, plotFlag)

if nargin < 9 || isempty(plotFlag)
    plotFlag = false;
end

barCandCutoff = stgs.barCandCutoff;
barDetCutoff = stgs.barDetCutoff;
resizeDims = stgs.resizeDims;

dGran = 5;
% thGran = pi/90;
thGran = pi/(2*180);

% radSchAmt = ceil(stgs.unsharpMaskAmt);
radSchAmt = ceil(stgs.unsharpMaskSigma);

rotAngle = fitParams.rotAngle;

stgsNoDeproj = stgs;
stgsNoDeproj.useDeProjectStretch = false;
pprImg = preprocessImage(img, stgsNoDeproj, [], fitParams);
% figure; imshow(pprImg);

[imgR, imgC] = meshgrid(1:1:size(pprImg, 2), size(pprImg, 1):-1:1);
testRegion = sqrt((imgR-ctrR).^2+(imgC-ctrC).^2) <= candRad;
testImg = pprImg .* testRegion;

if plotFlag
    barImg = repmat(testImg, [1 1 3]);
    barImg(:, :, 3) = barImg(:, :, 3) + (~testRegion * 0.5);
    figure; imshow(barImg)
end

xOff = (size(pprImg, 2)/2 + 0.5) - ctrC;
yOff = (size(pprImg, 1)/2 + 0.5) - ctrR;
[ht, dVals, thVals] = houghTransform(testImg, [], dGran, thGran, xOff, yOff);

% if plotFlag
%     figure; imagesc(ht.^2); axis image; 
%     title('htSq'); xlabel('th'); ylabel('d'); colormap(hot(256)); 
%     axis off
% end

% adjustment by rotAngle since standardardized image rotated compared to
% input image
% adjustment by pi/2 due to Hough transform parameterization
thDev = abs(thVals - mod(candTh - rotAngle + pi/2, pi));
thDev(thDev > pi/2) = pi - thDev(thDev > pi/2);
[temp, barAngIdx] = min(thDev);
[temp, barOppIdx] = max(thDev);

% if plotFlag
%     figure; plot(thDev);
% end

perRadBarScores = zeros(1, radSchAmt + 1);
thMax = max(ht.^2, [], 1);
curBarScore = thMax(barAngIdx) / thMax(barOppIdx);
perRadBarScores(1) = curBarScore;

radExtVals = 0:1:radSchAmt;
for radExtIdx=2:1:length(perRadBarScores)
    radExt = radExtVals(radExtIdx);
    prevTestRegion = testRegion;
    prevHt = ht;
    testRegion = sqrt((imgR-ctrR).^2+(imgC-ctrC).^2) <= (candRad + radExt);
    testImg = pprImg .* (testRegion & ~prevTestRegion);
    ht = houghTransform(testImg, [], dGran, thGran, xOff, yOff) + prevHt;
    thMax = max(ht.^2, [], 1);
    curBarScore = thMax(barAngIdx) / thMax(barOppIdx);
    perRadBarScores(radExtIdx) = curBarScore;
end

if plotFlag
    figure; plot(radExtVals, perRadBarScores)
    xlabel('radius extension (pixels)');
    ylabel('bar score');
end

[mV, mI] = max(perRadBarScores);
barScore = mV;
bestRadExt = radExtVals(mI);

% barInfo.stdzCtrR = resizeDims(1)/2;
% barInfo.stdzCtrC = resizeDims(2)/2;
barInfo.stdzCtrR = ctrR;
barInfo.stdzCtrC = ctrC;
barInfo.stdzAngle = candTh;
barInfo.stdzHalfLength = candRad + bestRadExt;
barInfo.score = barScore;
barInfo.candScore = candScore;

mu = fitParams.muFit;
barInfo.iptCtrR = size(img, 1) - mu(2) + 1;
barInfo.iptCtrC = mu(1);

barInfo.iptAngle = mod(barInfo.stdzAngle - fitParams.rotAngle, 2*pi);

cropRad = fitParams.cropRad;
endR = candRad * sin(barInfo.iptAngle);
endC = candRad * cos(barInfo.iptAngle);
endR = endR * ((2*cropRad+1)/resizeDims(1));
endC = endC * ((2*cropRad+1)/resizeDims(2));
barInfo.iptHalfLength = sqrt(endR.^2 + endC.^2);

barInfo.barDetected = (candScore > barCandCutoff) & (barScore > barDetCutoff);

fitParams.muFit
if plotFlag
    bImg = preprocessImage(img, stgs, fitParams);
    bImg = addBarOverlay(bImg, ctrR, ctrC, barInfo.stdzAngle, barInfo.stdzHalfLength);
    figure; imshow(bImg);
    
%     displayLgspPlot([], [], img, barInfo.iptCtrR, barInfo.iptCtrC, barInfo.iptAngle, barInfo.iptHalfLength);
    [bImg, barOvlInds] = addBarOverlay(img, barInfo.iptCtrR, barInfo.iptCtrC, barInfo.iptAngle, barInfo.iptHalfLength);
    figure; imshow(bImg);
    
    [oppImg, oppOvlInds] = addBarOverlay(img, barInfo.iptCtrR, barInfo.iptCtrC, mod(barInfo.iptAngle + pi/2, pi), barInfo.iptHalfLength);
    figure; hold all
    plot(barOvlInds, img(barOvlInds)); 
    plot(oppOvlInds, img(oppOvlInds));
    title('brightnesses at bar overlay indices');
    legend('bar-line', 'perpindicular');
    axv = axis();
    axv(3) = 0; axv(4) = 1;
    axis(axv);
end

% change "bar half-length" to "bar radius"

end