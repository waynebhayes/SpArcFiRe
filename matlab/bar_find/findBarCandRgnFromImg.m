function [bestRad, bestTh, candScore, scorePerRad] = ...
    findBarCandRgnFromImg(img, ctrR, ctrC, stgs, plotFlag)
% Experimental function; don't use
% INPUTS:
%   iptImg: the image originally used as input, without any preprocessing
%   img:
%   ctrR: row-value for the galaxy center, in the standardized image
%   ctrC: column-value for the galaxy center, in the standardized image
%   fitParams: fit parameters used when standardizing the image
%   plotFlag: optional; whether or not to display additional figures (e.g.,
%       for debugging)

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

% minimum allowed minMaxRatio, to avoid giving a high score to low radius
% values that are almost completely filled with high-orientation pixels
minMaxRatioThres = 1.5; 
minMaxRatioThres = 0;
minRad = 5;

stgsNoDeproj = stgs; stgs.useDeProjectStretch = false;
[imgp, imgc] = preprocessImage(img, stgsNoDeproj);

% imgp = oriFieldimgpgths(ofld);

% changing these granularities may change the scale of the bar score, and
% thus require changing the bar-detection cutoff
dGran = 5;
thGran = pi/90;
thGran = pi/36;
[imgC, imgR] = meshgrid(1:1:size(imgp, 2), 1:1:size(imgp, 1));

radVals = minRad:1:floor(min(size(imgp))/2);
scores = zeros(size(radVals));
perRadBestTh = zeros(size(radVals));
minMaxRatios = zeros(size(radVals));
inRgnPrev = false(size(imgp));
htPrev = 0;
for ii=1:1:length(radVals)
%     [score, th, mmRatio] = calcBarScoreAtRadius(radVals(ii));
    
    curRad = radVals(ii);

    inRgn = sqrt((imgR-ctrR).^2+(imgC-ctrC).^2) <= curRad;
    testImg = imgp .* (inRgn & ~inRgnPrev);
    
%     testImg = (testImg - 0.5) * 2;
    
    [ht, dVals, thVals] = ...
        houghTransform(testImg, [], dGran, thGran);
    ht = ht+htPrev;
    nTh = length(thVals);
    
    thMax = max(ht, [], 1);  
    [mV, mI] = max(thMax);

%     scores(ii) = thMax(mI) / curRad;
%     scores(ii) = thMax(mI);
%     scores(ii) = thMax(mI) / sum(sum(imgp .* inRgn))^0.5;
    scores(ii) = thMax(mI) / (curRad^(0.5));
    scores(ii) = thMax(mI) / (curRad^(0.85));
%     scores(ii) = thMax(mI);
    perRadBestTh(ii) = thVals(mI);
    minMaxRatios(ii) = thMax(mI) / min(thMax);
    
    inRgnPrev = inRgn;
    htPrev = ht;
end
if plotFlag
    figure; 
    subplot(1, 2, 1);
    axs = plotyy(radVals, scores, radVals, perRadBestTh * (180/pi)); 
    xlabel('radius'); 
    set(get(axs(1),'Ylabel'),'String','bar score') 
    set(get(axs(2),'Ylabel'),'String','theta') 
    
    subplot(1, 2, 2);
    axs = plotyy(radVals, scores, radVals, minMaxRatios); 
    xlabel('radius');
    set(get(axs(1),'Ylabel'),'String','bar score') 
    set(get(axs(2),'Ylabel'),'String','min-max ratio')  
end

gvecRad = 5;
gvec = normpdf(-gvecRad:gvecRad, 0, 2); gvec = gvec / sum(gvec);
scsm = conv(padarray(scores, [0 2*gvecRad], 'replicate'), gvec); 
extAmt = (length(scsm) - length(scores))/2;
scsm = scsm(1+extAmt:end-extAmt);
scChgs = [-inf diff(scsm) inf];
scLocMaxIdxs = find((scChgs(1:end-1) > 0 & scChgs(2:end) < 0))
if plotFlag
    figure; hold all; plot(scores); plot(scsm); title('scores vs smoothed scores');
end
if plotFlag
    nMax = length(scLocMaxIdxs);
    angles = zeros(1, nMax);
    barRadii = zeros(1, nMax);
    for ii=1:1:nMax
        curIdx = scLocMaxIdxs(ii);
        angles(ii) = mod(perRadBestTh(curIdx) + pi/2, pi);
        barRadii(ii) = radVals(curIdx);
    end
    ovlImg = addBarOverlay(imgp, ctrR, ctrC, angles, barRadii);
    figure; imshow(ovlImg); title('bar candidates');
end

[bestSc, bestIdx] = max(scores .* (minMaxRatios >= minMaxRatioThres));
bestRad = radVals(bestIdx);
% due to (d, theta) parameterization of Hough transform, Hough transform
% angle is at a right angle to the bar angle
bestTh = mod(perRadBestTh(bestIdx) + pi/2, pi);

candScore = bestSc;
scorePerRad = scores;

if plotFlag
%     displayLgspPlot([], [], imgp, ctrR, ctrC, bestTh, bestRad);
%     title(sprintf('best bar (score = %2.4f, rad = %d, th = %2.4f)', ...
%         bestSc, bestRad, bestTh));
    
    [bImg, barOvlInds] = addBarOverlay(imgp, ctrR, ctrC, bestTh, bestRad);
    figure; imshow(bImg);
    title(sprintf('best bar (score = %2.4f, rad = %d, th = %2.4f)', ...
        bestSc, bestRad, bestTh));
    
    figure; hold all
    [barOvlR, barOvlC] = ind2sub(size(imgp), barOvlInds);
    [barDists, sI] = sort(sqrt((barOvlR - ctrR).^2 + (barOvlC - ctrC).^2), 'ascend');
    barIntensities = imgc(barOvlInds(sI));
    scatter(barDists, barIntensities, 'b.'); 
    
    hlfWay = floor(length(barDists)/2);
    distForFit = barDists(hlfWay:end);
    itsyForFit = barIntensities(hlfWay:end);
    coeff = polyfit(distForFit, itsyForFit, 1);
    itsyPred = (coeff(1) * distForFit + coeff(2));
    ferr = itsyForFit - itsyPred;
    ssr = sum(ferr.^2);
    sst = (length(distForFit)-1) * var(itsyForFit);
    rsq = 1 - ssr/sst;
    mse = ssr / length(distForFit);
    figure; hold all
    scatter(distForFit, itsyForFit);
    plot(distForFit, itsyPred);
    title(sprintf('MSE = %2.4f, rsq = %2.4f', mse, rsq));
end


% function [score, th, minMaxRatio] = calcBarScoreAtRadius(rad)
%     inRgn = sqrt((imgR-ctrR).^2+(imgC-ctrC).^2) <= rad;
%     testImg = imgp .* inRgn;
% %     testImg = im2double((imgp > 0.15) & inRgn);
%     
%     [ht, dVals, thVals] = ...
%         houghTransform(testImg, [], dGran, thGran);
%     
%     nTh = length(thVals);
%     if mod(nTh, 2) ~= 0
%         error('(internal): number of theta-bins needs to be even');
%     end
%     
% %     thMax = max(ht.^2, [], 1);  
%     thMax = max(ht, [], 1);  
%     
%     [mV, mI] = max(thMax);
% 
%     % index of theta-value with 90-degree offset
% %     oppI = mod(mI + (nTh/2) - 1, nTh) + 1;  
%     
%     score = thMax(mI) / rad;
%     score = thMax(mI) / sqrt(rad);
%     score = thMax(mI);
%     th = thVals(mI);
%     
%     minMaxRatio = thMax(mI) / min(thMax);
% 
% end % calcBarScoreAtRadius

end % findBarFromOriField