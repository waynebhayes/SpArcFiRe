function [barSc, barInfo] = ...
    houghFindBar(img, mu, covar, gxyMajAxsLen, stgs, plotFlag)
% Analyzes the given bar-candidate region, returning information about the
%   most likely bar, along with a score for how likely the galaxy 
%   actually has a bar.
% NOTE: this was used in an earlier bar-detection attempt, but is not used
%   in the current flow
% INPUTS:
%   img: grayscale galaxy image
%   mu: mean for the fitted Gaussian for the bar-candidate region
%   covar: covariance matrix for the fitted Gaussian for the bar-candidate
%       region
%   gxyMajAxsLen: major axis length, in pixels, of the ellipse fitted to 
%       the edge of the galaxy disk
%   stgs: 
%   plotFlag:
% OUTPUTS:
%   barSc: the score for how likely the galaxy actually has a bar with the
%       reported parameters
%   barInfo: the parameters for the bar, reported for ALL galaxies,
%       regardless of whether there is a bar present.  Use the barSc
%       parameter for a measure of whether the bar actually exists.

if nargin < 7 || isempty(plotFlag)
    plotFlag = false;
end

barDetContourLvl = stgs.barDetContourLvl;
barLenContourLvl = stgs.barLenContourLvl;

% barLikCutoff = 10^-6;
% barLikCutoff = 10^-9;
dGran = 5;
thGran = pi/90;
% dGran = 2

[majAxsVec, majAxsAngle, axisRatio, majAxsLen] = ...
        findCovarElpsAxes(covar, barDetContourLvl, size(img));
semiMajAxsLen = majAxsLen/2;
clear majAxsLen

gxySemiMajAxsLen = gxyMajAxsLen / 2;
clear gxyMajAxsLen;
% in preparation for deleting the redundant parameter gxySemiMajAxsLen
assert(abs(semiMajAxsLen - gxySemiMajAxsLen) < 10^-6);

lik = calcLikelihoods(size(img), mu, covar);

[imgR, imgC] = meshgrid(1:1:size(img, 2), size(img, 1):-1:1);
testRegion = sqrt((imgR-mu(1)).^2+(imgC-mu(2)).^2) <= semiMajAxsLen;

if plotFlag
    figure; hold on
    imagesc(flipud(0.5 * (testRegion .* img) + 0.5 * img)); axis image; colormap gray
    contour(flipud(lik), barDetContourLvl, '-g')
end

stgsUsmAdj = stgs;
stgsUsmAdj.unsharpMaskSigma = ...
    ceil(stgs.unsharpMaskSigma * (gxySemiMajAxsLen/mean(stgs.resizeDims)));
imgUsm = unsharpMask(img, stgsUsmAdj);
imgH = imgUsm .* (testRegion);
imgH(imgH > 1) = 1; imgH(imgH < 0) = 0;
[ht, dVals, thVals] = houghTransform(imgH, [], dGran, thGran);
% if plotFlag
%     displayTopNLines(imgH, ht, dVals, thVals);
% end

majAxsAngleH = mod(majAxsAngle + pi/2, pi);

if plotFlag
    figure; imagesc(ht.^2); axis image; 
    title('htSq'); xlabel('th'); ylabel('d'); colormap(hot(256)); 
    axis off
end

thDev = abs(thVals - majAxsAngleH);
thDev(thDev > pi/2) = pi - thDev(thDev > pi/2);
[temp, majAxsIdx] = min(thDev);
[temp, minAxsIdx] = max(thDev);

% squaring the HT values emphasizes concentrated points; avoids issue of
% all angles otherwise having equal votes
thMax = max(ht.^2, [], 1);  
axisVoteRatio = thMax(majAxsIdx) / thMax(minAxsIdx);

[tmp1, bdMajAxsAngle, tmp2, bdMajAxsLen] = ...
        findCovarElpsAxes(covar, barLenContourLvl, size(img));
barInfo.iptCtrR = mu(2);
barInfo.iptCtrC = mu(1);
barInfo.iptAngle = bdMajAxsAngle;
barInfo.iptHalfLength = bdMajAxsLen/2;
barSc = axisVoteRatio;

% figure; hold on
% imagesc(flipud(img)); axis image; colormap gray
% contour(flipud(lik), [barDetContourLvl, barLenContourLvl], 'g');

% figure; 
% imshow(img); hold on
% scatter(mu(1), size(img, 2) - mu(2), 'go')
rho = -semiMajAxsLen:0.5:semiMajAxsLen;
th = majAxsAngle * ones(1, length(rho));
[x, y] = pol2cart(th, rho);
rows = round((size(img, 1) - mu(2)) - y);
cols = round(mu(1) + x);
barLoc = zeros(size(img));
barLoc(sub2ind(size(img), rows, cols)) = 1;
ol = repmat(img, [1 1 3]);
ol(:, :, 2) = ol(:, :, 2) + barLoc;
if plotFlag
    figure; imshow(ol);
end

% figure; imshow(img); hold on; contour(lik, [10^-5, 10^-6, 10^-7], '-g')
% figure; imshow(img); hold on; contour(lik, [10^-5], '-g')
% figure; imshow(img); hold on; contour(lik, [10^-6], '-g')
% figure; imshow(img); hold on; contour(lik, [10^-7], '-g')

end