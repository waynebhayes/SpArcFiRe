function displayCz2CmpPlotsForCriterion(tightFrac, medFrac, looseFrac, paMeas, meetsCriterion, criterionName, outDir)

if nargin < 7 || isempty(outDir)
    writeOutput = false;
else
    writeOutput = true;
    outDir = [outDir filesep];
    plotNum = 1;
    screen_size = get(0, 'ScreenSize');
end

tightFrac = tightFrac(meetsCriterion);
medFrac = medFrac(meetsCriterion);
looseFrac = looseFrac(meetsCriterion);
paMeas = paMeas(meetsCriterion);

fprintf('comparisons for %s:\n', criterionName);

nCmpGxy = length(tightFrac);

[tmp, domTightVote] = max([tightFrac medFrac looseFrac], [], 2);

fprintf('GZ2 sum of vote fractions for tight arms: %2.4f\n', sum(tightFrac));
fprintf('GZ2 sum of vote fractions for medium arms: %2.4f\n', sum(medFrac));
fprintf('GZ2 sum of vote fractions for loose arms: %2.4f\n', sum(looseFrac));

gz2Tight = (domTightVote == 1);
gz2Med = (domTightVote == 2);
gz2Loose = (domTightVote == 3);
nTight = sum(gz2Tight);
nMed = sum(gz2Med);
nLoose = sum(gz2Loose);
fprintf('GZ2 number of galaxies with majority vote tight: %d of %d (%2.4f%%)\n', nTight, nCmpGxy, 100 * nTight/nCmpGxy);
fprintf('GZ2 number of galaxies with majority vote medium: %d of %d (%2.4f%%)\n', nMed, nCmpGxy, 100 * nMed/nCmpGxy);
fprintf('GZ2 number of galaxies with majority vote loose: %d of %d (%2.4f%%)\n', nLoose, nCmpGxy, 100 * nLoose/nCmpGxy);

paAlenwGz2Tight = paMeas(gz2Tight);
paAlenwGz2Med = paMeas(gz2Med);
paAlenwGz2Loose = paMeas(gz2Loose);

% figure; hist(paAlenwGz2Tight); title('pitch angle (alenw-dco) distribution for GZ2 tight arms');
% figure; hist(paAlenwGz2Med); title('pitch angle (alenw-dco) distribution for GZ2 medium arms');
% figure; hist(paAlenwGz2Loose); title('pitch angle (alenw-dco) distribution for GZ2 loose arms');

binGran = 2;
paBins = -90:binGran:90;
paBinIdxs = paMeas / binGran - (paBins(1) / binGran);
paBinIdxs = floor(paBinIdxs) + 1;

tCnts = accumarray(paBinIdxs, gz2Tight, [length(paBins) 1]);
mCnts = accumarray(paBinIdxs, gz2Med, [length(paBins) 1]);
lCnts = accumarray(paBinIdxs, gz2Loose, [length(paBins) 1]);
% % tCnts = histc(paAlenwGz2Tight, paBins);
% % mCnts = histc(paAlenwGz2Med, paBins);
% % lCnts = histc(paAlenwGz2Loose, paBins);
% figure; bar(paBins, [tCnts mCnts lCnts]); 
% title('measured pitch angle distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
% figure; bar(paBins, [tCnts/nTight mCnts/nMed lCnts/nLoose]); 
% title('measured pitch angle normalized distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
%close all
if writeOutput
    figure; hold all; 
    plot(paBins, tCnts); plot(paBins, mCnts); plot(paBins, lCnts);
    title('measured pitch angle distributions');
    legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'paDist.png']);
    plotNum = plotNum + 1;
    close gcf
end
if writeOutput
    figure; hold all; 
    plot(paBins, tCnts/nTight); plot(paBins, mCnts/nMed); plot(paBins, lCnts/nLoose);
    title('measured pitch angle normalized distributions');
    legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'paNormDist.png']);
    plotNum = plotNum + 1;
    close gcf
end

tWtdCnts = accumarray(paBinIdxs, tightFrac, [length(paBins) 1]);
mWtdCnts = accumarray(paBinIdxs, medFrac, [length(paBins) 1]);
lWtdCnts = accumarray(paBinIdxs, looseFrac, [length(paBins) 1]);
% figure; bar(paBins, [tWtdCnts mWtdCnts lWtdCnts]); 
% title('measured pitch angle distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
% figure; bar(paBins, [tWtdCnts/sum(tightFrac) mWtdCnts/sum(medFrac) lWtdCnts/sum(looseFrac)]); 
% title('measured pitch angle normalized distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
if writeOutput
    figure; hold all; 
    plot(paBins, tWtdCnts); plot(paBins, mWtdCnts); plot(paBins, lWtdCnts);
    title('measured pitch angle weighted distributions');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'paWtdDist.png']);
    plotNum = plotNum + 1;
    close gcf
end
if writeOutput
    figure; hold all; 
    plot(paBins, tWtdCnts/sum(tightFrac)); plot(paBins, mWtdCnts/sum(medFrac)); plot(paBins, lWtdCnts/sum(looseFrac));
    title('measured pitch angle normalized weighted distributions');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'paNmWtdDist.png']);
    plotNum = plotNum + 1;
    close gcf
end

gz2TightnessSc = 1 * tightFrac + 2 * medFrac + 3 * looseFrac;

if writeOutput
    figure; scatter(paMeas, gz2TightnessSc);
    xlabel('calculated pitch angle');
    ylabel('GZ2 vote-based tightness score');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'paTscScatter.png']);
    plotNum = plotNum + 1;
    close gcf
end


paBins = 0:binGran:90;
paBinIdxs = abs(paMeas) / binGran - (paBins(1) / binGran);
paBinIdxs = floor(paBinIdxs) + 1;
% paBinIdxsPrev = paBinIdxs;
% paBinIdxs = ones(size(paMeas));
% for ii=2:1:length(paBins)
%     paBinIdxs(abs(paMeas) >= paBins(ii)) = ii;
% end
% paBinIdxsTst = paBinIdxs;
% sum(paBinIdxs == paBinIdxsTst)
% paBinIdxs = paBinIdxsPrev;
% save('dbg', 'paBins', 'paMeas', 'paBinIdxs', 'paBinIdxsTst');
tCnts = accumarray(paBinIdxs, gz2Tight, [length(paBins) 1]);
mCnts = accumarray(paBinIdxs, gz2Med, [length(paBins) 1]);
lCnts = accumarray(paBinIdxs, gz2Loose, [length(paBins) 1]);
% figure; bar(paBins, [tCnts mCnts lCnts]); 
% title('measured pitch angle distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
% figure; bar(paBins, [tCnts/nTight mCnts/nMed lCnts/nLoose]); 
% title('measured pitch angle normalized distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
%close all
if writeOutput
    figure; hold all; 
    plot(paBins, tCnts); plot(paBins, mCnts); plot(paBins, lCnts);
    title('measured pitch angle distributions');
    legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'absPaDist.png']);
    plotNum = plotNum + 1;
    close gcf
end
if writeOutput
    figure; hold all; 
    plot(paBins, tCnts/nTight); plot(paBins, mCnts/nMed); plot(paBins, lCnts/nLoose);
    title('measured pitch angle normalized distributions');
    legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'absNmPaDist.png']);
    plotNum = plotNum + 1;
    close gcf
end

sumCnts = tCnts + mCnts + lCnts;
tProp = tCnts ./ sumCnts;
mProp = mCnts ./ sumCnts;
lProp = lCnts ./ sumCnts;

% error bars determined by multinomial standard deviations
% errorbar(paBins, tProp, tProp.*(1-tProp));
% errorbar(paBins, mProp, mProp.*(1-mProp));
% errorbar(paBins, lProp, lProp.*(1-lProp));
if writeOutput
    figure; 
%     set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    hold all;
    plot(paBins, tProp); plot(paBins, mProp); plot(paBins, lProp); 
    title('majority vote proportions by pitch angle');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'majVtPropByPa.png']);
    plotNum = plotNum + 1;
    close gcf
end

tWtdCnts = accumarray(paBinIdxs, tightFrac, [length(paBins) 1]);
mWtdCnts = accumarray(paBinIdxs, medFrac, [length(paBins) 1]);
lWtdCnts = accumarray(paBinIdxs, looseFrac, [length(paBins) 1]);
if writeOutput
    figure; hold all; 
    plot(paBins, tWtdCnts); plot(paBins, mWtdCnts); plot(paBins, lWtdCnts);
    title('measured pitch angle weighted distributions');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'absPaWtdDist.png']);
    plotNum = plotNum + 1;
    close gcf
end
tWtdNormCnts = tWtdCnts/sum(tightFrac); mWtdNormCnts = mWtdCnts/sum(medFrac); lWtdNormCnts = lWtdCnts/sum(looseFrac);
if writeOutput
    figure; hold all; 
    plot(paBins, tWtdNormCnts); plot(paBins, mWtdNormCnts); plot(paBins, lWtdNormCnts);
    title('measured pitch angle normalized weighted distributions');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'absPaNmWtdDist.png']);
    plotNum = plotNum + 1;
    close gcf
end

wtSums = tWtdCnts + mWtdCnts + lWtdCnts;
tProp = tWtdCnts ./ wtSums;
mProp = mWtdCnts ./ wtSums;
lProp = lWtdCnts ./ wtSums;
if writeOutput
    figure; hold all;
    plot(paBins, tProp); plot(paBins, mProp); plot(paBins, lProp);
    % plot(paBins, tProp); plot(paBins, mProp); plot(paBins, lProp);
    title('GZ2 weighted vote proportions by pitch angle');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'wtVtPropByPa.png']);
    plotNum = plotNum + 1;
    close gcf
end

wtNormSums = tWtdNormCnts + mWtdNormCnts + lWtdNormCnts;
tNormProp = tWtdNormCnts ./ wtNormSums;
mNormProp = mWtdNormCnts ./ wtNormSums;
lNormProp = lWtdNormCnts ./ wtNormSums;
if writeOutput
    figure; hold all;
    plot(paBins, tNormProp); plot(paBins, mNormProp); plot(paBins, lNormProp);
    title('GZ2 weighted normalized vote proportions by pitch angle');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'wtdNmPropByPa.png']);
    plotNum = plotNum + 1;
    close gcf
end

vts = sort([tightFrac medFrac looseFrac], 2, 'descend');
maxVts = vts(:, 1); midVts = vts(:, 2); minVts = vts(:, 3);
hsted = 0:0.05:1;
maxCts = histc(maxVts, hsted);
midCts = histc(midVts, hsted);
minCts = histc(minVts, hsted);
if writeOutput
    figure; hold all
    plot(hsted, maxCts); plot(hsted, midCts); plot(hsted, minCts);
    xlabel('vote proportions');
    ylabel('count');
    legend('highest vote fraction', '2nd highest vote fraction', 'lowest vote fraction');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'maxMidMinDist.png']);
    plotNum = plotNum + 1;
    close gcf
end


binCnts = accumarray(paBinIdxs, ones(size(paBinIdxs)), [length(paBins) 1]);
if writeOutput
    figure; bar(paBins + (binGran/2), binCnts); 
    title('distribution of measured pitch angles');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'distMeasPa.png']);
    plotNum = plotNum + 1;
    close gcf
end

meanMaxVts = accumarray(paBinIdxs, maxVts, [length(paBins) 1]) ./ binCnts;
meanMidVts = accumarray(paBinIdxs, midVts, [length(paBins) 1]) ./ binCnts;
meanMinVts = accumarray(paBinIdxs, minVts, [length(paBins) 1]) ./ binCnts;
means = zeros(size(paBinIdxs)); means = meanMaxVts(paBinIdxs);
stdMaxVts = sqrt(accumarray(paBinIdxs, (maxVts - means).^2, [length(paBins) 1]) ./ (binCnts-1));
means = zeros(size(paBinIdxs)); means = meanMidVts(paBinIdxs);
stdMidVts = sqrt(accumarray(paBinIdxs, (midVts - means).^2, [length(paBins) 1]) ./ (binCnts-1));
means = zeros(size(paBinIdxs)); means = meanMinVts(paBinIdxs);
stdMinVts = sqrt(accumarray(paBinIdxs, (minVts - means).^2, [length(paBins) 1]) ./ (binCnts-1));
if writeOutput
    figure; hold all
    % plot(paBins, meanMaxVts); plot(paBins, meanMidVts); plot(paBins, meanMinVts);
    errorbar(paBins, meanMaxVts, stdMaxVts);
    errorbar(paBins, meanMidVts, stdMidVts);
    errorbar(paBins, meanMinVts, stdMinVts);
    xlabel('pitch angle');
    ylabel('mean/sd vote proportion');
    legend('largest vote proportion', 'middle vote proportion', 'smallest vote proportion');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'meanSdVotePropByPa.png']);
    plotNum = plotNum + 1;
    close gcf
end

medianMaxVts = accumarray(paBinIdxs, maxVts, [length(paBins) 1], @median);
medianMidVts = accumarray(paBinIdxs, midVts, [length(paBins) 1], @median);
medianMinVts = accumarray(paBinIdxs, minVts, [length(paBins) 1], @median);
lqMaxVts = accumarray(paBinIdxs, maxVts, [length(paBins) 1], @(x)(quantile(x, 0.25)));
lqMidVts = accumarray(paBinIdxs, midVts, [length(paBins) 1], @(x)(quantile(x, 0.25)));
lqMinVts = accumarray(paBinIdxs, minVts, [length(paBins) 1], @(x)(quantile(x, 0.25)));
uqMaxVts = accumarray(paBinIdxs, maxVts, [length(paBins) 1], @(x)(quantile(x, 0.75)));
uqMidVts = accumarray(paBinIdxs, midVts, [length(paBins) 1], @(x)(quantile(x, 0.75)));
uqMinVts = accumarray(paBinIdxs, minVts, [length(paBins) 1], @(x)(quantile(x, 0.75)));
if writeOutput
    figure; hold all
    errorbar(paBins, medianMaxVts, medianMaxVts - lqMaxVts, uqMaxVts - medianMaxVts);
    errorbar(paBins, medianMidVts, medianMidVts - lqMidVts, uqMidVts - medianMidVts);
    errorbar(paBins, medianMinVts, medianMinVts - lqMinVts, uqMinVts - medianMinVts);
    xlabel('pitch angle');
    ylabel('median/quartile vote proportion');
    legend('largest vote proportion', 'middle vote proportion', 'smallest vote proportion');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'qtleVotePropByPa.png']);
    plotNum = plotNum + 1;
    close gcf
end

tscGran = 0.05;
tscBins = 1:tscGran:3;
tscBinIdxs = gz2TightnessSc / tscGran - (tscBins(1) / tscGran);
tscBinIdxs = floor(tscBinIdxs) + 1;
tCnts = accumarray(tscBinIdxs, gz2Tight, [length(tscBins) 1]);
mCnts = accumarray(tscBinIdxs, gz2Med, [length(tscBins) 1]);
lCnts = accumarray(tscBinIdxs, gz2Loose, [length(tscBins) 1]);
% figure; bar(tscBins, [tCnts mCnts lCnts]); 
% title('GZ2 tightness score distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
% figure; bar(tscBins, [tCnts/nTight mCnts/nMed lCnts/nLoose]); 
% title('GZ2 tightness score normalized distributions');
% legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
% %close all
if writeOutput
    figure; hold all; 
    plot(tscBins, tCnts); plot(tscBins, mCnts); plot(tscBins, lCnts);
    title('GZ2 tightness score distributions');
    legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'tscDist.png']);
    plotNum = plotNum + 1;
    close gcf
end
    if writeOutput
    figure; hold all; 
    plot(tscBins, tCnts/nTight); plot(tscBins, mCnts/nMed); plot(tscBins, lCnts/nLoose);
    title('GZ2 tightness score normalized distributions');
    legend('GZ2 majority tight', 'GZ2 majority med', 'GZ2 majority loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'tscNmDist.png']);
    plotNum = plotNum + 1;
    close gcf
end

sumCnts = tCnts + mCnts + lCnts;
tProp = tCnts ./ sumCnts;
mProp = mCnts ./ sumCnts;
lProp = lCnts ./ sumCnts;
if writeOutput
    figure; hold all;
    plot(tscBins, tProp); plot(tscBins, mProp); plot(tscBins, lProp); 
    title('majority vote proportions by tightness score');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'majVtPropByTsc.png']);
    plotNum = plotNum + 1;
    close gcf
end

tWtdCnts = accumarray(tscBinIdxs, tightFrac, [length(tscBins) 1]);
mWtdCnts = accumarray(tscBinIdxs, medFrac, [length(tscBins) 1]);
lWtdCnts = accumarray(tscBinIdxs, looseFrac, [length(tscBins) 1]);
wtSums = tWtdCnts + mWtdCnts + lWtdCnts;
tProp = tWtdCnts ./ wtSums;
mProp = mWtdCnts ./ wtSums;
lProp = lWtdCnts ./ wtSums;
if writeOutput
    figure; hold all;
    plot(tscBins, tProp); plot(tscBins, mProp); plot(tscBins, lProp);
    title('GZ2 weighted vote proportions by tightness score');
    legend('GZ2 weighted tight', 'GZ2 weighted med', 'GZ2 weighted loose');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'wtdVtPropByTsc.png']);
    plotNum = plotNum + 1;
    close gcf
end

vtsn = vts ./ repmat(sum(vts, 2), 1, 3);
ent = ent3(vtsn(:, 1), vtsn(:, 2), vtsn(:, 3));
if writeOutput
    figure; hist(ent, 100); title('entropy distribution');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'entropyDist.png']);
    plotNum = plotNum + 1;
    close gcf
end
fprintf('min, max possible entropy: %2.4f, %2.4f\n', ent3(1, 0, 0), ent3(1/3, 1/3, 1/3));
entQlvls = [0 0.25 .5 .75 1];
fprintf('entropy quantiles: \n');
for ii=1:1:length(entQlvls)
    fprintf('\t %2.4f%%\t %2.4f\n', 100 * entQlvls(ii), quantile(ent, entQlvls(ii)));
end

[p1g, p2g] = meshgrid(0:0.001:1, 0:0.001:1);
p3g = 1 - p1g - p2g;
feas = p3g >= 0;
% p1 = p1g(feas); p2 = p2g(feas); p3 = p3g(feas);
p1g(~feas) = NaN; p2g(~feas) = NaN; p3g(~feas) = NaN;
% entd = ent3(p1, p2, p3);
if writeOutput
    figure; contour(p1g, p2g, ent3(p1g, p2g, p3g), quantile(ent, [0.25, 0.5, 0.75]));
    colormap jet
    colorbar
    title(['entropy contours, quantile levels: ' sprintf('%2.2f ', quantile(ent, [0.25, 0.5, 0.75]))]);;
%     legend('min', 'first quartile', 'median', 'third quantile', 'max');
    saveas(gcf, [outDir criterionName sprintf('_%d_', plotNum) 'entropyContours.png']);
    plotNum = plotNum + 1;
    close gcf
end

end % displayPlotsForCriterion