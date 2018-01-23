function [bestRad, bestTh, candScore] = ...
    findBarCandRgn(ofld, ctrR, ctrC, plotFlag)
% Finds the length and angle of the "best" bar in the image, and attempts 
% to quantify the linear strength of this bar. One best bar will be found 
% regardless of whether a bar actually exists.
% INPUTS:
%   ofld: orientation field (computed on the fully standardized image)
%   ctrR: row-value for the galaxy center, in the standardized image
%   ctrC: column-value for the galaxy center, in the standardized image
%   plotFlag: optional; whether or not to display additional figures (e.g.,
%       for debugging)
% OUTPUTS:
%   bestRad: radius of the best bar
%   bestTh:  angle of the best bar
%   candScore: one measure of the linear strength of the bar; an attempted
%       measure of whether a bar actually exists

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

% minimum allowed minMaxRatio, to avoid giving a high score to low radius
% values that are almost completely filled with high-orientation pixels
minMaxRatioThres = 1.5; 
minRad = 5;

ofldStren = oriFieldStrengths(ofld);

% changing these granularities may change the scale of the bar score, and
% thus require changing the bar-detection cutoff
dGran = 5;
% thGran = pi/90;
thGran = pi/(2*180);
[imgC, imgR] = meshgrid(1:1:size(ofldStren, 2), 1:1:size(ofldStren, 1));

radVals = minRad:1:floor(min(size(ofldStren))/2);
scores = zeros(size(radVals));
perRadBestTh = zeros(size(radVals));
minMaxRatios = zeros(size(radVals));
inRgnPrev = false(size(ofldStren));
htPrev = 0;
xOff = (size(ofldStren, 2)/2 + 0.5) - ctrC;
yOff = (size(ofldStren, 1)/2 + 0.5) - ctrR;
for ii=1:1:length(radVals)
%     [score, th, mmRatio] = calcBarScoreAtRadius(radVals(ii));
    
    curRad = radVals(ii);

    inRgn = sqrt((imgR-ctrR).^2+(imgC-ctrC).^2) <= curRad;
    testImg = ofldStren .* (inRgn & ~inRgnPrev);
    
    [ht, dVals, thVals] = ...
        houghTransform(testImg, [], dGran, thGran, xOff, yOff);
    ht = ht+htPrev;
    nTh = length(thVals);
    
    thMax = max(ht, [], 1);  
    [mV, mI] = max(thMax);

    scores(ii) = thMax(mI) / curRad;
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
    
%     figure;
%     hold all
%     plot(radVals, minMaxRatios);
%     xlabel('radius');
%     ylabel('max-min ratio');
%     plot(radVals, 1.25 * ones(size(radVals)));
%     plot(radVals, 1.5 * ones(size(radVals)));
end

[bestSc, bestIdx] = max(scores .* (minMaxRatios >= minMaxRatioThres));
bestRad = radVals(bestIdx);
% due to (d, theta) parameterization of Hough transform, Hough transform
% angle is at a right angle to the bar angle
bestTh = mod(perRadBestTh(bestIdx) + pi/2, pi);

candScore = bestSc;

if plotFlag
    displayLgspPlot([], [], ofldStren, ctrR, ctrC, bestTh, bestRad);
    title(sprintf('best bar (score = %2.4f, rad = %d, th = %2.4f)', ...
        bestSc, bestRad, bestTh));
end

function [score, th, minMaxRatio] = calcBarScoreAtRadius(rad)
    inRgn = sqrt((imgR-ctrR).^2+(imgC-ctrC).^2) <= rad;
    testImg = ofldStren .* inRgn;
    
    [ht, dVals, thVals] = ...
        houghTransform(testImg, [], dGran, thGran);
    
    nTh = length(thVals);
    if mod(nTh, 2) ~= 0
        error('(internal): number of theta-bins needs to be even');
    end
    
%     thMax = max(ht.^2, [], 1);  
    thMax = max(ht, [], 1);  
    
    [mV, mI] = max(thMax);

    % index of theta-value with 90-degree offset
%     oppI = mod(mI + (nTh/2) - 1, nTh) + 1;  
    
    score = thMax(mI) / rad;
    th = thVals(mI);
    
    minMaxRatio = thMax(mI) / min(thMax);

end % calcBarScoreAtRadius

end % findBarFromOriField