function [angle, halfLength] = ...
    fitBarUsingGaussian(img, ctrR, ctrC, plotFlag)
% Obsolete, should NOT be used!
% Attempts to heuristically determine whether the given cluster describes a
%   bar, by looking at distance from the center and the "lopsidedness" of
%   pixels along the major and minor axes of the fitted Gaussian.
% INPUTS:
%   img: the cluster matrx to check for a bar
%   ctrR: the row index of the center of the log-spirals
%   ctrC: the column index of the center of the log-spirals
%   plotFlag: whether to plot the bar-fit (optional, default false)
% OUTPUTS:
%   angle: the angle, in radians, of any detected bar, or empty if no bar
%   halfLength: the half-length of any detected bar, or empty if no bar

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

likCutoff = 10^-4;
maxMuDist = 10;
maxLop = 2;
ctrDist = round(size(img, 1) / 16);

imgCtr = size(img) / 2;
[mu, covar] = calcWtdMeanCovar(img);
% [mu, covar] = iterFitGauss(img, false);
muDist = sqrt(sum([mu - imgCtr].^2));

isCtr = clusMtxIsThruCtr(img, ctrR, ctrC, ctrDist);

if ~plotFlag && (muDist > maxMuDist || ~isCtr)
    angle = [];
    halfLength = [];
    return;
end

% if plotFlag
%     likAvg = sum(sum(calcLikelihoods(size(img), mu, covar) .* img)) / ...
%         sum(sum(img))
% end

[majAxsVec, angle, axisRatio] = findCovarElpsAxes(covar, likCutoff);

[r, c, brt] = find(img);
xVals = c - ctrC;
yVals = ctrR - r; % = (imgSz - r + 1) - (imgSz - ctrR + 1);

majAxProjns = [xVals yVals] * majAxsVec; % TODO: reverse majAxsVec?
% halfLength = max(abs(majAxProjns));
halfLength = min(max(majAxProjns(majAxProjns > 0)), ...
    max(-majAxProjns(majAxProjns < 0)));

sumGTz = sum(majAxProjns(majAxProjns > 0));
sumLTz = sum(abs(majAxProjns(majAxProjns < 0)));
lopsMajAxs = max(sumGTz, sumLTz) / min(sumGTz, sumLTz);
minAxsVec = (majAxsVec' * [0, -1; 1, 0])';
minAxProjns = [xVals yVals] * minAxsVec;
sumGTz = sum(minAxProjns(minAxProjns > 0));
sumLTz = sum(abs(minAxProjns(minAxProjns < 0)));
lopsMinAxs = max(sumGTz, sumLTz) / min(sumGTz, sumLTz);
% lopsidedness = max(lopsMajAxs, lopsMinAxs);
lopsidedness = (lopsMajAxs + lopsMinAxs) / 2;

if plotFlag
    figure; 
    imagesc([1 - imgCtr(2), size(img, 2) - imgCtr(2)], ...
        [1 - imgCtr(1), size(img, 1) - imgCtr(1)], img);
    axis image; colormap gray
    xLim = halfLength * cos(angle);
    yLim = halfLength * sin(angle);
    line([xLim -xLim], [-yLim, yLim], 'Color', 'g');
    title(sprintf(...
        'muDist = %2.4f, lopsidedness = %2.4f (%2.4f, %2.4f), axisRatio = %2.4f', ...
        muDist, lopsidedness, lopsMajAxs, lopsMinAxs, axisRatio));
end

if ~isCtr || (muDist > maxMuDist) || (lopsidedness > maxLop)
    angle = [];
    halfLength = [];
    return;
end

% gtz = zeros(size(img));
% gtz(sub2ind(size(img), r(majAxProjns > 0), c(majAxProjns > 0))) = brt(majAxProjns > 0);
% figure; imagesc(gtz); axis image
% ltz = zeros(size(img));
% ltz(sub2ind(size(img), r(majAxProjns < 0), c(majAxProjns < 0))) = brt(majAxProjns < 0);
% figure; imagesc(ltz); axis image
% 
% gtz = zeros(size(img));
% gtz(sub2ind(size(img), r(minAxProjns > 0), c(minAxProjns > 0))) = brt(minAxProjns > 0);
% figure; imagesc(gtz); axis image
% ltz = zeros(size(img));
% ltz(sub2ind(size(img), r(minAxProjns < 0), c(minAxProjns < 0))) = brt(minAxProjns < 0);
% figure; imagesc(ltz); axis image
% 
% majAxsVec
% minAxsVec



% imgSz = size(img);
% 
% [r, c, brt] = find(img);
% brtSum = sum(brt);
% xVals = c - ctrC;
% yVals = ctrR - r; % = (imgSz - r + 1) - (imgSz - ctrR + 1);
% 
% % if there is a bar, its center should be at the center of the galaxy
% mu = [ctrC (imgSz(1) - ctrR + 1)];
% 
% xDev = xVals - mu(1);
% yDev = yVals - mu(2);
% valsCov11 = xDev .^ 2;
% valsCov22 = yDev .^ 2;
% valsCov12 = xDev .* yDev;
% cov11 = sum(sum(brt .* valsCov11)) / brtSum;
% cov22 = sum(sum(brt .* valsCov22)) / brtSum;
% cov12 = sum(sum(brt .* valsCov12)) / brtSum;
% 
% covar = [cov11, cov12; cov12, cov22]
% 
% lik = calcLikelihoods(imgSz, mu, covar);
% 
% if plotFlag
%     figure; hold on
%     imagesc(flipud(img)); axis image; colormap gray; impixelinfo
%     contour(flipud(lik), '-r')
%     contour(flipud(lik), [10^-8, 10^-9], '-g')
%     hold off
% end

end