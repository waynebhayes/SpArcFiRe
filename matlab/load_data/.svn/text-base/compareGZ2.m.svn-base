function compareGZ2(resultCsvPath, nSkip, cmpTbl, destDir)

% if nargin < 2 || isempty(cmpTbl)
if ~isempty(resultCsvPath)
%     outDirPath = [outDirPath filesep];

    [autoTbl, nSkip] = readTableFromFile(resultCsvPath, ',');
    gzTbl = readTableFromFile('D:\Galaxy-Sets\Bamford GZ2\gz2master.csv', ',');
    [autoTblI, gzTblI] = intersectRowsByField(autoTbl, gzTbl, 1, 1, true);
    clear autoTbl
    clear gzTbl
    cmpTbl = joinTables(autoTblI, gzTblI, 1, 1);
    clear autoTblI
    clear gzTblI
    save cmpTbl
    return;
end

if nargin < 4
    destDir = [];
else
    destDir = [destDir filesep];
end

% pctTightVote = (tightFrac + medFrac + looseFrac);
% figure; hist(pctTightVote);

nCmpGxy = size(cmpTbl, 1) - 1;

fprintf('%d galaxies skipped due to output unavailability\n', nSkip);
fprintf('availability: %d of %d (%2.4f%%)\n', nCmpGxy, (nCmpGxy+nSkip), 100 * (nCmpGxy/(nCmpGxy+nSkip)));

fprintf('comparing %d galaxies\n', nCmpGxy);

nameIdx = strmatch('name', cmpTbl(1, :), 'exact');
assert(length(nameIdx) == 1)
names = cmpTbl(2:end, nameIdx);

tightIdx = strmatch('t10_arms_winding_a28_tight_weighted_fraction', cmpTbl(1, :), 'exact');
medIdx = strmatch('t10_arms_winding_a29_medium_weighted_fraction', cmpTbl(1, :), 'exact');
looseIdx = strmatch('t10_arms_winding_a30_loose_weighted_fraction', cmpTbl(1, :), 'exact');

tightFrac = str2double(cmpTbl(2:end, tightIdx));
medFrac = str2double(cmpTbl(2:end, medIdx)); 
looseFrac = str2double(cmpTbl(2:end, looseIdx)); 

paAlenwDcoIdx = strmatch('pa_alenWtd_avg_domChiralityOnly', cmpTbl(1, :), 'exact');
paAlenwDco = str2double(cmpTbl(2:end, paAlenwDcoIdx));

vts = sort([tightFrac medFrac looseFrac], 2, 'descend');
vtsn = vts ./ repmat(sum(vts, 2), 1, 3);
ent = ent3(vtsn(:, 1), vtsn(:, 2), vtsn(:, 3));

maxVt = max(vts, [], 2);

% displayCz2CmpPlotsForCriterion(tightFrac, medFrac, looseFrac, paAlenwDco, true(1, nCmpGxy), 'allAvail', outDirPath);
% displayCz2CmpPlotsForCriterion(tightFrac, medFrac, looseFrac, paAlenwDco, ent <= median(ent), 'best_50pct_ent', outDirPath);
% displayCz2CmpPlotsForCriterion(tightFrac, medFrac, looseFrac, paAlenwDco, ent <= quantile(ent, 0.25), 'best_25pct_ent', outDirPath);
% displayCz2CmpPlotsForCriterion(tightFrac, medFrac, looseFrac, paAlenwDco, maxVt >= 1/2, 'highestVtGeHalf', outDirPath);
% displayCz2CmpPlotsForCriterion(tightFrac, medFrac, looseFrac, paAlenwDco, maxVt >= 2/3, 'highestVtGeTwoThirds', outDirPath);

paBins = [0:2:40 50 90];
paBins = [0:2:40 45 50 90];
paBins = [0:2:40 90];

paMeasAll = abs(paAlenwDco);
% paBinIdxsAll = ones(size(paMeasAll));
% for ii=2:1:length(paBins)
%     paBinIdxsAll(abs(paMeasAll) >= paBins(ii)) = ii;
% end
% % tightFracAll = tightFrac;
% % medFracAll = medFrac;
% % looseFracAll = looseFrac;
% [maxTightVote, domTightVoteAll] = max([tightFrac medFrac looseFrac], [], 2);
% figure; hist(maxTightVote, 100); title('strongest tightness vote');
% gz2TightAll = (domTightVoteAll == 1);
% gz2MedAll = (domTightVoteAll == 2);
% gz2LooseAll = (domTightVoteAll == 3);
% tCntsAll = accumarray(paBinIdxsAll, gz2TightAll, [length(paBins) 1]);
% mCntsAll = accumarray(paBinIdxsAll, gz2MedAll, [length(paBins) 1]);
% lCntsAll = accumarray(paBinIdxsAll, gz2LooseAll, [length(paBins) 1]);
% sumCntsAll = tCntsAll + mCntsAll + lCntsAll;
% tPropAll = tCntsAll ./ sumCntsAll;
% mPropAll = mCntsAll ./ sumCntsAll;
% lPropAll = lCntsAll ./ sumCntsAll;
% % f1 = figure; hold all;
% % plot(paBins, tPropAll); plot(paBins, mPropAll); plot(paBins, lPropAll); 
% % axis([0 50 0 1])
% % set(f1, 'Position', [100 100 500 250]);
paBinIdxsAll = getPaBinIdxs(paMeasAll, paBins);
plotPaFrequencies(paBins, paBinIdxsAll)
if ~isempty(destDir)
    export_fig([destDir 'pitch_angle_frequencies_all.pdf']);
end
[tPropAll, mPropAll, lPropAll] = getTightnessProps(paBinIdxsAll, tightFrac, medFrac, looseFrac, false);
% fprintf('majority-vote-medium pitch angle mean: %2.4f, sd: %2.4f\n',...
%     mean(paMeasAll(gz2TightAll)), std(paMeasAll(gz2TightAll)));

isTop = ent <= quantile(ent, 0.25);
% paMeasTop = abs(paAlenwDco(isTop));
% paBinIdxsTop = ones(size(paMeasTop));
% for ii=2:1:length(paBins)
%     paBinIdxsTop(abs(paMeasTop) >= paBins(ii)) = ii;
% end
% % assert(all(paBinIdxsTop == paBinIdxsAll(isTop)))
% tightFracTop = tightFrac(isTop);
% medFracTop = medFrac(isTop);
% looseFracTop = looseFrac(isTop);
% [tmp, domTightVoteTop] = max([tightFracTop medFracTop looseFracTop], [], 2);
% gz2TightTop = (domTightVoteTop == 1);
% gz2MedTop = (domTightVoteTop == 2);
% gz2LooseTop = (domTightVoteTop == 3);
% tCntsTop = accumarray(paBinIdxsTop, gz2TightTop, [length(paBins) 1]);
% mCntsTop = accumarray(paBinIdxsTop, gz2MedTop, [length(paBins) 1]);
% lCntsTop = accumarray(paBinIdxsTop, gz2LooseTop, [length(paBins) 1]);
% sumCntsTop = tCntsTop + mCntsTop + lCntsTop;
% tPropTop = tCntsTop ./ sumCntsTop;
% mPropTop = mCntsTop ./ sumCntsTop;
% lPropTop = lCntsTop ./ sumCntsTop;
% f2 = figure; hold all;
% plot(paBins, tPropTop); plot(paBins, mPropTop); plot(paBins, lPropTop); 
% axis([0 50 0 1])
% set(f2, 'Position', [100 100 500 200]);
% fprintf('majority-vote-medium pitch angle mean: %2.4f, sd: %2.4f\n',...
%     mean(paMeasTop(gz2TightTop)), std(paMeasTop(gz2TightTop)));
paBinIdxsTop = getPaBinIdxs(paMeasAll(isTop), paBins);
plotPaFrequencies(paBins, paBinIdxsTop)
if ~isempty(destDir)
    export_fig([destDir 'pitch_angle_frequencies_top_quartile.pdf']);
end
[tPropTop, mPropTop, lPropTop] = getTightnessProps(...
    paBinIdxsTop, tightFrac(isTop), medFrac(isTop), looseFrac(isTop), false);

% f3 = figure; hold all;
% cordr = get(gca, 'ColorOrder');
% set(gca, 'ColorOrder', cordr(1:3, :))
% plot(paBins, tPropAll, '-'); plot(paBins, mPropAll, '-'); plot(paBins, lPropAll, '-'); 
% plot(paBins, tPropTop, ':'); plot(paBins, mPropTop, ':'); plot(paBins, lPropTop, ':'); 
% axis([0 50 0 1])
% xlabel('measured pitch angle'); ylabel('vote proportion');
% set(f3, 'Position', [100 100 340 200]);
% 
% f4 = figure; hold all;
% cordr = get(gca, 'ColorOrder');
% set(gca, 'ColorOrder', cordr(1:3, :))
% plot(paBins, tPropAll, ':', 'LineWidth', 1.5); plot(paBins, mPropAll, ':', 'LineWidth', 1.5); plot(paBins, lPropAll, ':', 'LineWidth', 1.5); 
% plot(paBins, tPropTop, '-', 'LineWidth', 1.5); plot(paBins, mPropTop, '-', 'LineWidth', 1.5); plot(paBins, lPropTop, '-', 'LineWidth', 1.5); 
% axis([0 50 0 1])
% xlabel('measured pitch angle'); ylabel('vote proportion');
% set(f4, 'Position', [100 100 640 200]);

f5 = figure; hold all;
cordr = get(gca, 'ColorOrder');
set(gca, 'ColorOrder', cordr(1:3, :))
plot(paBins, tPropAll, ':', 'LineWidth', 1.5); plot(paBins, mPropAll, ':', 'LineWidth', 1.5); plot(paBins, lPropAll, ':', 'LineWidth', 1.5); 
plot(paBins, tPropTop, '-', 'LineWidth', 1.5); plot(paBins, mPropTop, '-', 'LineWidth', 1.5); plot(paBins, lPropTop, '-', 'LineWidth', 1.5); 
axis([0 50 0 1])
hXlbl = xlabel('Measured Pitch Angle (Degrees)'); 
hYlbl = ylabel('GZ Majority Vote Proportion');
xlim(paBins([1 end-1]))
set(f5, 'Position', [100 100 650 350]);
set(f5, 'Color', [1 1 1]);
set(hXlbl, 'FontSize', 18);
set(hYlbl, 'FontSize', 18);
set(gca, 'FontSize', 16);
if ~isempty(destDir)
    export_fig([destDir 'GZ2_tightness_comparison.pdf']);
end

[tPropAllWtd, mPropAllWtd, lPropAllWtd] = getTightnessProps(paBinIdxsAll, tightFrac, medFrac, looseFrac, true);
[tPropTopWtd, mPropTopWtd, lPropTopWtd] = getTightnessProps(...
    paBinIdxsTop, tightFrac(isTop), medFrac(isTop), looseFrac(isTop), true);

f6 = figure; hold all;
cordr = get(gca, 'ColorOrder');
set(gca, 'ColorOrder', cordr(1:3, :))
plot(paBins, tPropAllWtd, ':', 'LineWidth', 1.5); plot(paBins, mPropAllWtd, ':', 'LineWidth', 1.5); plot(paBins, lPropAllWtd, ':', 'LineWidth', 1.5); 
plot(paBins, tPropTopWtd, '-', 'LineWidth', 1.5); plot(paBins, mPropTopWtd, '-', 'LineWidth', 1.5); plot(paBins, lPropTopWtd, '-', 'LineWidth', 1.5); 
axis([0 50 0 1])
hXlbl = xlabel('Measured Pitch Angle (Degrees)'); 
hYlbl = ylabel('GZ Vote Weight Proportion');
xlim(paBins([1 end-1]))
set(f6, 'Position', [100 100 650 350]);
set(f6, 'Color', [1 1 1]);
set(hXlbl, 'FontSize', 18);
set(hYlbl, 'FontSize', 18);
set(gca, 'FontSize', 16);
title('weighted');

% ----- bar detection -----

barCandScoreIdx = strmatch('bar_cand_score_orifld', cmpTbl(1, :), 'exact');
assert(length(barCandScoreIdx) == 1)
barScoreIdx = strmatch('bar_score_img', cmpTbl(1, :), 'exact');
assert(length(barScoreIdx) == 1)
barAvailIdx = strmatch('bar_candidate_available', cmpTbl(1, :), 'exact');
assert(length(barAvailIdx) == 1)
barUsedIdx = strmatch('bar_used', cmpTbl(1, :), 'exact');
assert(length(barUsedIdx) == 1)

cmpBarFracIdx = strmatch('t03_bar_a06_bar_weighted_fraction', cmpTbl(1, :), 'exact');
assert(length(cmpBarFracIdx) == 1)
% % cmpBarFlagIdx = strmatch('t03_bar_a06_bar_flag', cmpTbl(1, :), 'exact');
% % assert(length(cmpBarFlagIdx) == 1)

barCandScores = str2double(cmpTbl(2:end, barCandScoreIdx));
barScores = str2double(cmpTbl(2:end, barScoreIdx));
barFracs = str2double(cmpTbl(2:end, cmpBarFracIdx));
% % barFlags = str2double(cmpTbl(2:end, cmpBarFlagIdx));
bar_used_strs = cmpTbl(2:end, barUsedIdx);
bar_used = strcmpi('true', bar_used_strs);
bar_not_used = strcmpi('false', bar_used_strs);
assert(all(bar_used == ~bar_not_used));
bar_avail_strs = cmpTbl(2:end, barAvailIdx);
bar_available = strcmpi('true', bar_avail_strs);
bar_not_available = strcmpi('false', bar_avail_strs);
assert(all(bar_available == ~bar_not_available));

bar_prop_edges = [-inf 0:0.1:1 inf];
GZ2frac_hist_bar_used = histc(barFracs(bar_used), bar_prop_edges);
GZ2frac_hist_bar_not_used = histc(barFracs(~bar_used), bar_prop_edges);
GZ2frac_hist_all = histc(barFracs, bar_prop_edges);
figure; hold all; 
plot(bar_prop_edges, GZ2frac_hist_bar_used); 
plot(bar_prop_edges, GZ2frac_hist_bar_not_used); 
plot(bar_prop_edges, GZ2frac_hist_all);
figure; hold all; 
plot(bar_prop_edges, GZ2frac_hist_bar_used / sum(GZ2frac_hist_bar_used)); 
plot(bar_prop_edges, GZ2frac_hist_bar_not_used / sum(GZ2frac_hist_bar_not_used)); 
plot(bar_prop_edges, GZ2frac_hist_all / sum(GZ2frac_hist_all));

bar_prop_edges = [-inf 0:0.05:1 inf];
p_bar_used = zeros(length(bar_prop_edges)-1, 1);
std_bar_used = zeros(length(bar_prop_edges)-1, 1);
p_bar_available = zeros(length(bar_prop_edges)-1, 1);
for ii=1:length(p_bar_used)
    in_cur_bin = (barFracs <= bar_prop_edges(ii)) & (barFracs < bar_prop_edges(ii+1));
    p_bar_used(ii) = sum(bar_used(in_cur_bin)) / sum(in_cur_bin);
    cur_p = p_bar_used(ii);
    std_bar_used(ii) = sqrt((cur_p * (1 - cur_p))/sum(in_cur_bin));
    p_bar_available(ii) = sum(bar_available(in_cur_bin)) / sum(in_cur_bin);
end
figure('Position', [100 100 650 300]); hold all
plot(bar_prop_edges(2:end-1), p_bar_used(2:end), 'b-', 'LineWidth', 1.5);
plot(bar_prop_edges(2:end-1), p_bar_available(2:end), 'g--', 'LineWidth', 1.5);
hXlbl = xlabel('GZ2 Bar Weighted Fraction'); 
hYlbl = ylabel('Bar Detection Rate');
set(gcf, 'Color', [1 1 1]);
set(hXlbl, 'FontSize', 18);
set(hYlbl, 'FontSize', 18);
set(gca, 'FontSize', 16);
if ~isempty(destDir)
    export_fig([destDir 'GZ2_bar_comparison.pdf']);
end

bad_bar_detections = (barFracs == 0) & (bar_used);
fprintf('worst bad bar detections (0%% bar vote):\n')
names(bad_bar_detections)
strong_misses = (barFracs == 1) & ~bar_used;
strong_miss_names = names(strong_misses);
fprintf('some strong bar misses(100%% bar vote):\n');
strong_miss_names(1:10)
fprintf('some good bar detections:\n');
good_bar_detection_names = names((barFracs == 1) & bar_used);
good_bar_detection_names(1:min(10, length(good_bar_detection_names)))
fprintf('some good bar negatives\n');
good_bar_negative_names = names((barFracs == 0) & ~bar_used);
good_bar_negative_names(1:10)
% std_bar_used
% figure; errorbar(bar_prop_edges(2:end-1), p_bar_used(2:end),
% std_bar_used(2:end));

% ----- arm count -----

close all;

oneArmIdx = strmatch('t11_arms_number_a31_1_weighted_fraction', cmpTbl(1, :), 'exact');
twoArmsIdx = strmatch('t11_arms_number_a32_2_weighted_fraction', cmpTbl(1, :), 'exact');
threeArmsIdx = strmatch('t11_arms_number_a33_3_weighted_fraction', cmpTbl(1, :), 'exact');
fourArmsIdx = strmatch('t11_arms_number_a34_4_weighted_fraction', cmpTbl(1, :), 'exact');
moreThanFourArmsIdx = strmatch('t11_arms_number_a36_more_than_4_weighted_fraction', cmpTbl(1, :), 'exact');
cannotTellArmsIdx = strmatch('t11_arms_number_a37_cant_tell_weighted_fraction', cmpTbl(1, :), 'exact');

oneArmFrac = str2double(cmpTbl(2:end, oneArmIdx));
twoArmsFrac = str2double(cmpTbl(2:end, twoArmsIdx));
threeArmsFrac = str2double(cmpTbl(2:end, threeArmsIdx));
fourArmsFrac = str2double(cmpTbl(2:end, fourArmsIdx));
moreThanFourArmsFrac = str2double(cmpTbl(2:end, moreThanFourArmsIdx));
cannotTellArmsFrac = str2double(cmpTbl(2:end, cannotTellArmsIdx));

% armcount_names = {'# arms >= 102.4', '# arms >= 51.2'};
% armcount_names = {'# arms (largest length gap)', '# arms (arc-length function flattening)',...
%     '# arms >= 307.2', '# arms >= 294.4', '# arms >= 281.6', '# arms >= 268.8',...
%     '# arms >= 256.0', '# arms >= 243.2', '# arms >= 230.4', '# arms >= 217.6',...
%     '# arms >= 204.8', '# arms >= 192.0', '# arms >= 179.2', '# arms >= 166.4',...
%     '# arms >= 153.6', '# arms >= 140.8', '# arms >= 128.0', '# arms >= 115.2',...
%     '# arms >= 102.4', '# arms >= 89.6', '# arms >= 76.8', '# arms >= 64.0',...
%     '# arms >= 51.2', '# arms >= 38.4', '# arms >= 25.6', '# arms >= 12.8', '# arms >= 0.0'};
% armcount_names = {'# arms (largest length gap)', '# arms (arc-length function flattening)',...
%     '# arms >= 166.4', '# arms >= 153.6', '# arms >= 140.8'};
% armcount_names = {'# arms >= 166.4', '# arms >= 153.6', '# arms >= 140.8'};
% armcount_names = {'# arms >= 153.6', '# arms >= 89.6', '# arms >= 76.8'};
% armcount_names = {'# arms >= 89.6'};
% armcount_names = {'numArcsGE085', 'numArcsGE090', 'numArcsGE095'};
% armcount_names = {'numArcsGE085', 'numArcsGE090', 'numArcsGE095', 'numDcoArcsGE085', 'numDcoArcsGE090', 'numDcoArcsGE095'};
armcount_names = {'numArcsGE060', 'numDcoArcsGE060', 'numArcsGE075', 'numDcoArcsGE075', 'numArcsGE090', 'numDcoArcsGE090'};
armcount_names = {'numDcoArcsGE075', 'numDcoArcsGE090', 'numDcoArcsGE100'};
armcount_names = {'numDcoArcsGE050', 'numDcoArcsGE055', 'numDcoArcsGE060', 'numDcoArcsGE065', 'numDcoArcsGE070', 'numDcoArcsGE075', 'numDcoArcsGE090', 'numDcoArcsGE100'};
armFracs = horzcat(oneArmFrac, twoArmsFrac, threeArmsFrac, fourArmsFrac, moreThanFourArmsFrac, cannotTellArmsFrac);
armFracs = armFracs ./ repmat(sum(armFracs, 2), [1, size(armFracs, 2)]);
entropy = armFracs .* log(armFracs);
entropy(armFracs == 0) = 0;
entropy = -sum(entropy, 2);
isTopArmCount = entropy <= quantile(entropy, 0.25);
all_armcounts = cell(1, length(armcount_names)+1);
for ii=1:length(armcount_names)
    armcount_idx = strmatch(armcount_names{ii}, cmpTbl(1, :), 'exact');
    armcounts = str2double(cmpTbl(2:end, armcount_idx));
    all_armcounts{ii} = armcounts;
    plot_title = armcount_names{ii};
%     fprintf([plot_title '\n']);
    plotArmCountComparison(armcounts, oneArmFrac, twoArmsFrac,...
        threeArmsFrac, fourArmsFrac, moreThanFourArmsFrac, cannotTellArmsFrac, ...
        plot_title, destDir);
    plot_title = [plot_title '-HighestEntropyQuartile'];
%     fprintf([plot_title '\n']);
    plotArmCountComparison(armcounts(isTopArmCount), oneArmFrac(isTopArmCount),...
        twoArmsFrac(isTopArmCount), threeArmsFrac(isTopArmCount),...
        fourArmsFrac(isTopArmCount), moreThanFourArmsFrac(isTopArmCount),...
        cannotTellArmsFrac(isTopArmCount), plot_title, destDir);
end
[max_gz_armcount_vote, max_gz_armcount] = max(armFracs, [], 2);
plotArmCountComparison(max_gz_armcount-1, oneArmFrac, twoArmsFrac,...
        threeArmsFrac, fourArmsFrac, moreThanFourArmsFrac, cannotTellArmsFrac, ...
        'Strongest GZ armcount', [], true);
xlabel('GZ majority vote');
set(gca, 'xTickLabel', {'1', '2', '3', '4', '>4', '?'})
vote_prop_bins = 0:0.01:1;
% figure('Position', [100 100 500 250]); 
% bar(vote_prop_bins, histc(max_gz_armcount_vote, vote_prop_bins));
% xlim([0.01, 1.01])
% % title('Max GZ armcount vote');
% set(gcf, 'Color', [1 1 1]);
% hXlbl = xlabel('Maximum Vote Proportion');
% hYlbl = ylabel('Frequency');
% set(hXlbl, 'FontSize', 18);
% set(hYlbl, 'FontSize', 18);
% set(gca, 'FontSize', 16);

gz_armcount_for_dist = max_gz_armcount;
% set "don't know" to -1
gz_armcount_for_dist(gz_armcount_for_dist == 6) = -1;
all_armcounts{end} = gz_armcount_for_dist;
printArmCountDistributionTable([armcount_names 'GZ2 Max Vote'], all_armcounts, isTopArmCount, length(oneArmFrac))

[max_gz_armcount_vote, max_gz_armcount] = max(armFracs, [], 2);
plotArmCountComparison(max_gz_armcount(isTopArmCount)-1, oneArmFrac(isTopArmCount),...
    twoArmsFrac(isTopArmCount), threeArmsFrac(isTopArmCount), ...
    fourArmsFrac(isTopArmCount), moreThanFourArmsFrac(isTopArmCount), ...
    cannotTellArmsFrac(isTopArmCount), ...
    'Strongest GZ armcount (highest entropy quartile)', [], true);
xlabel('GZ majority vote');
set(gca, 'xTickLabel', {'1', '2', '3', '4', '>4', '?'})

% figure('Position', [100 100 500 300]); 
% bar(vote_prop_bins, histc(max_gz_armcount_vote(isTopArmCount), vote_prop_bins));
% xlim([0.01, 1.01])
% title('Max GZ armcount vote (higest entropy quartile)');
% set(gcf, 'Color', [1 1 1]);
% hXlbl = xlabel('Maximum Vote Proportion');
% hYlbl = ylabel('Frequency');
% set(hXlbl, 'FontSize', 18);
% set(hYlbl, 'FontSize', 18);
% set(gca, 'FontSize', 16);

figure('Position', [100 100 750 300]);
max(max_gz_armcount_vote(isTopArmCount))
max(max_gz_armcount_vote(~isTopArmCount))
hb = bar(vote_prop_bins, horzcat(...
    histc(max_gz_armcount_vote(isTopArmCount), vote_prop_bins),...
    histc(max_gz_armcount_vote(~isTopArmCount), vote_prop_bins)),...
    'stack');
colormap([0, 0.75, 0; 0, 0, 1]);
set(hb, 'EdgeColor', 'none');
xlim([0.01, 1.01])
% title('Max GZ armcount vote');
set(gcf, 'Color', [1 1 1]);
hXlbl = xlabel('Maximum Vote Proportion');
hYlbl = ylabel('Frequency');
set(hXlbl, 'FontSize', 18);
set(hYlbl, 'FontSize', 18);
set(gca, 'FontSize', 16);
if ~isempty(destDir)
    export_fig([destDir 'Max_GZ_armcount_vote_distribution.pdf']);
end

fprintf('GZ2 average weight for dominant category: \n');
for ii=min(max_gz_armcount):max(max_gz_armcount)
    fprintf('\t%d:\t%2.4f\n', ii, mean(max_gz_armcount_vote(max_gz_armcount == ii)));
end

max_gz_armcount_vote_topArmcount = max_gz_armcount_vote(isTopArmCount);
max_gz_armcount_topArmcount = max_gz_armcount(isTopArmCount);
fprintf('GZ2 average weight for dominant category (highest entropy quartile): \n');
for ii=min(max_gz_armcount):max(max_gz_armcount)
    fprintf('\t%d:\t%2.4f\n', ii, mean(max_gz_armcount_vote_topArmcount(max_gz_armcount_topArmcount == ii)));
end

% arms_gt_102_4_idx = strmatch('# arms >= 102.4', cmpTbl(1, :), 'exact');
% arms_gt_102_4 = str2double(cmpTbl(2:end, arms_gt_102_4_idx));
% plotArmCountComparison(arms_gt_102_4, oneArmFrac, twoArmsFrac,...
%     threeArmsFrac, fourArmsFrac, moreThanFourArmsFrac, cannotTellArmsFrac);
% 
% arms_gt_51_2_idx = strmatch('# arms >= 51.2', cmpTbl(1, :), 'exact');
% arms_gt_51_2 = str2double(cmpTbl(2:end, arms_gt_51_2_idx));
% plotArmCountComparison(arms_gt_51_2, oneArmFrac, twoArmsFrac,...
%     threeArmsFrac, fourArmsFrac, moreThanFourArmsFrac, cannotTellArmsFrac);

function paBinIdxs = getPaBinIdxs(paMeas, paBins)
    paBinIdxs = ones(size(paMeas));
    for binIdx=2:1:length(paBins)
        paBinIdxs(abs(paMeas) >= paBins(binIdx)) = binIdx;
    end

end

function plotPaFrequencies(paBins, paBinIdxs)
    fprintf('pitch angle bins: %s\n', mat2str(paBins));
    fprintf('pitch angle bin counts: \n');
    paBinCounts = zeros(1, length(paBins)-1);
    barLabels = cell(size(paBinCounts));
    for jj=1:length(paBins)-1
        binCount = sum(paBinIdxs == jj);
        paBinCounts(jj) = binCount;
        barLabel = sprintf('[%d,%d)', paBins(jj), paBins(jj+1));
        fprintf('%s:\t%d\n', barLabel, binCount);
        barLabels{jj} = barLabel;
    end
    figure('Position', [100 100 1000 300]); 
    hb = bar(paBinCounts, 'histc');
    set(hb, 'FaceColor', [0 1 1]);
    set(gca, 'XTick', 1:length(paBins));
    set(gca, 'XTickLabel', paBins);
    xlim([1 length(paBins)]);
    ylim([0 1.05 * max(paBinCounts)]);
    hXlbl = xlabel('Measured Pitch Angle (Degrees)');
    hYlbl = ylabel('Frequency');
    set(gcf, 'Color', [1 1 1]);
    set(hXlbl, 'FontSize', 18);
    set(hYlbl, 'FontSize', 18);
    set(gca, 'FontSize', 16);
end

function [tProp, mProp, lProp] = getTightnessProps(...
        paBinIdxs, tightFrac, medFrac, looseFrac, weighted)
if weighted
    gz2Tight = tightFrac;
    gz2Med = medFrac;
    gz2Loose = looseFrac;
else
    [maxTightVote, domTightVote] = max([tightFrac medFrac looseFrac], [], 2);
    % figure; hist(maxTightVote, 100); title('strongest tightness vote');
    gz2Tight = (domTightVote == 1);
    gz2Med = (domTightVote == 2);
    gz2Loose = (domTightVote == 3);
end
if weighted
    vote_type = 'weighted';
else
    vote_type = 'majority';
end
vote_sum = sum([gz2Tight' gz2Med' gz2Loose']);
fprintf('overall %s vote percentages: \n', vote_type);
fprintf('\ttight: %2.4f\n\tmedium: %2.4f\n\tloose: %2.4f\n',...
    sum(gz2Tight) / vote_sum, ...
    sum(gz2Med) / vote_sum, ...
    sum(gz2Loose) / vote_sum);
tCnts = accumarray(paBinIdxs, gz2Tight, [length(paBins) 1]);
mCnts = accumarray(paBinIdxs, gz2Med, [length(paBins) 1]);
lCnts = accumarray(paBinIdxs, gz2Loose, [length(paBins) 1]);
sumCnts = tCnts + mCnts + lCnts;
tProp = tCnts ./ sumCnts;
mProp = mCnts ./ sumCnts;
lProp = lCnts ./ sumCnts;
end

function printArmCountDistributionTable(armcount_method_names, all_armcounts, isTopQuartile, n_inst)
    assert(length(armcount_method_names) == length(all_armcounts))
    if ~isempty(destDir)
        outFile = fopen([destDir 'armcount_categorization_distributions.tex'], 'wt');
    else
        outFile = 1;
    end
    fprintf('number of arm-count instances: %d\n', n_inst)
    fprintf(outFile, '{\\newcolumntype{C}{>{\\centering\\arraybackslash}X}');
    fprintf(outFile, '\\begin{tabularx}{\\textwidth}{|c|c|C C C C C C C|}\n');
    fprintf(outFile, '\\hline \n');
%     fprintf('armcount criterion, top quartile only?, cannot tell, 0 arms, 1 arm, 2 arms, 3 arms, 4 arms, more than 4 arms, not included\n');
    fprintf(outFile, ['\\textbf{Measurement} & \\textbf{Top?} & '...
        '\\textbf{0 Arms} & \\textbf{1 Arm} &'...
        '\\textbf{2 Arms} & \\textbf{3 Arms} & \\textbf{4 Arms} & \\textbf{$>$4 Arms} & \\textbf{?\\phantom{''}Arms} \\\\ \\hline \\hline\n']);
    for jj=1:length(armcount_method_names)
        cur_armcounts = all_armcounts{jj};
        cur_armcounts(cur_armcounts > 5) = 5;
%         category_counts = zeros(1, 8);
        category_counts = zeros(1, 7);
        for category_idx = -1:1:5
            category_counts(category_idx + 2) = sum(cur_armcounts == category_idx);
        end
        category_counts = category_counts([2:end 1]);
%         category_counts(end) = n_inst - length(cur_armcounts);
        name_for_table = strrep(armcount_method_names{jj}, 'numDcoArcsGE', 'DWD, $\ge$');
        if ~isempty(strfind(name_for_table, 'DWD, '))
            name_for_table = [name_for_table 'px'];
        end
        fprintf(outFile, '\\multirow{2}{*}{\\textbf{%s}} & \\textbf{N} & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f \\\\ \n',...
            name_for_table, 100 * category_counts / length(cur_armcounts));
        
        cur_armcounts = cur_armcounts(isTopQuartile);
%         category_counts = zeros(1, 8);
        category_counts = zeros(1, 7);
        for category_idx = -1:1:5
            category_counts(category_idx + 2) = sum(cur_armcounts == category_idx);
        end
        category_counts = category_counts([2:end 1]);
%         category_counts(end) = n_inst - length(cur_armcounts);
        fprintf(outFile, '\\cline{2-9} & \\textbf{Y} & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f & %2.2f \\\\ \\hline\n',...
            100 * category_counts / length(cur_armcounts));
    end
    fprintf(outFile, '\\end{tabularx}}\n');
    if ~isempty(destDir)
        fclose(outFile);
    end
end

end % function compareGZ2(cmpTbl)
