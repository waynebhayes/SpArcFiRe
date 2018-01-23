function compareChiralityGZ2(resultCsvPath, cmpCsvPath)

% if nargin < 1 || isempty(resultCsvPath)
%     resultCsvPath = 'D:\matlab-output-local\[paper] 2011-11-07 GZ2 secondRun\GZ2_full_secondRun.csv';
% end
% 
% if nargin < 2 || isempty(cmpCsvPath)
%     cmpCsvPath = 'D:\Galaxy-Sets\GalaxyZoo1_DR_table2.csv';
% end

[autoTbl, nSkip, skippedNames] = readTableFromFile(resultCsvPath, ',');
gzTbl = readTableFromFile(cmpCsvPath, ',');
skippedNames = [{'name'}; skippedNames];
[gzSkipTbl, tmp] = intersectRowsByField(gzTbl, skippedNames, 1, 1, true);
[autoTblI, gzTblI] = intersectRowsByField(autoTbl, gzTbl, 1, 1, true);
clear gzTbl
if size(autoTbl, 1) ~= size(autoTblI, 1)
    warning('some auto-generated entries not present in comparison file');
    pause
end
cmpTbl = joinTables(autoTblI, gzTblI, 1, 1);

nameIdx = strmatch('name', cmpTbl(1, :), 'exact');
assert(length(nameIdx) == 1)

cbrIdx = strmatch('contourBrtRatio', cmpTbl(1, :), 'exact');
assert(length(cbrIdx) == 1)
t2aIdx = strmatch('top2_chirality_agreement', cmpTbl(1, :), 'exact');
assert(length(t2aIdx) == 1)

pcwIdx = strmatch('P_CW', cmpTbl(1, :), 'exact');
assert(length(pcwIdx) == 1)
pacwIdx = strmatch('P_ACW', cmpTbl(1, :), 'exact');
assert(length(pacwIdx) == 1)

cMajIdx = strmatch('chirality_maj', cmpTbl(1, :), 'exact');
assert(length(cMajIdx) == 1)
cAlenwIdx = strmatch('chirality_alenWtd', cmpTbl(1, :), 'exact');
assert(length(cAlenwIdx) == 1)
cWpsIdx = strmatch('chirality_wtdPangSum', cmpTbl(1, :), 'exact');
assert(length(cWpsIdx) == 1)
cLongestIdx = strmatch('chirality_longestArc', cmpTbl(1, :), 'exact');
assert(length(cLongestIdx) == 1)

% cbr = str2double(cmpTbl(2:end, cbrIdx));
% cbrOK = (cbr > 1);
longestAgree = cellfun(@(x)(strcmp(x(2:end-1), 'agree')), cmpTbl(2:end, t2aIdx));

names = cmpTbl(2:end, nameIdx);
gzv = str2double(cmpTbl(2:end, [pcwIdx pacwIdx]));
chirGz = gzv(:, 1) > gzv(:, 2);
chirMaj = cmpTbl(2:end, cMajIdx);
chirAlenw = cmpTbl(2:end, cAlenwIdx);
chirWps = cmpTbl(2:end, cWpsIdx);
chirLongest = cmpTbl(2:end, cLongestIdx);
gzMaxChirVote = max(gzv, [], 2);

chirMatchMaj = (strcmp('S-wise', chirMaj) & chirGz) | (strcmp('Z-wise', chirMaj) & ~chirGz);
chirMatchAlenw = (strcmp('S-wise', chirAlenw) & chirGz) | (strcmp('Z-wise', chirAlenw) & ~chirGz);
chirMatchWps = (strcmp('S-wise', chirWps) & chirGz) | (strcmp('Z-wise', chirWps) & ~chirGz);
chirMatchLongest = (strcmp('S-wise', chirLongest) & chirGz) | (strcmp('Z-wise', chirLongest) & ~chirGz);

print_agreement_vs_discernability(gzMaxChirVote, chirMatchAlenw, 'alenWtd DCO, all categories');
gzChirOnlyVote = gzMaxChirVote ./ sum(gzv, 2);
gzChirOnlyVote(isnan(gzChirOnlyVote)) = 0;
print_agreement_vs_discernability(gzChirOnlyVote, chirMatchAlenw, 'alenWtd DCO, winding direction only');

% figure; hist(gzChirAgree, 100);

agrEdges = [0:0.01:1];
agrCnts = histc(gzMaxChirVote, agrEdges);
figure;
plot(agrEdges, cumsum(agrCnts) / numel(gzMaxChirVote));
title('max chirality vote distribution for galaxies in sample');
xlabel('proportion designating dominant chirality'); ylabel('proportion included');

thresVals = 0:0.01:1;

gz2MaxVtDist = getThresDistn(thresVals, gzMaxChirVote, chirMatchAlenw);
gz2MaxVtDistT2Aonly = getThresDistn(thresVals, gzMaxChirVote(longestAgree), chirMatchAlenw(longestAgree));
figure; hold all
plot(thresVals, gz2MaxVtDist);
plot(thresVals, gz2MaxVtDistT2Aonly);
xlabel('GZ2 min proportion designating dominant chirality');
ylabel('GZ2 match rate (alenwDco)');
legend('all', 'top 2 agree', 'Location', 'Best');
% axs = axis; axs(3:4) = [.8 1]; axis(axs);
axis([0 1 .75 1]);

gz2MaxVtDist = getThresDistn(thresVals, gzMaxChirVote, chirMatchLongest);
gz2MaxVtDistT2Aonly = getThresDistn(thresVals, gzMaxChirVote(longestAgree), chirMatchLongest(longestAgree));
figure; hold all
plot(thresVals, gz2MaxVtDist);
plot(thresVals, gz2MaxVtDistT2Aonly);
xlabel('GZ2 min proportion designating dominant chirality');
ylabel('GZ2 match rate (chirality of longest)');
legend('all', 'top 2 agree', 'Location', 'Best');
% axs = axis; axs(3:4) = [.8 1]; axis(axs);
axis([0 1 .75 1]);

gz2MaxVtDist = getThresDistn(thresVals, gzMaxChirVote, chirMatchMaj);
gz2MaxVtDistT2Aonly = getThresDistn(thresVals, gzMaxChirVote(longestAgree), chirMatchMaj(longestAgree));
figure; hold all
plot(thresVals, gz2MaxVtDist);
plot(thresVals, gz2MaxVtDistT2Aonly);
xlabel('GZ2 min proportion designating dominant chirality');
ylabel('GZ2 match rate (majority vote)');
legend('all', 'top 2 agree', 'Location', 'Best');
% axs = axis; axs(3:4) = [.8 1]; axis(axs);
axis([0 1 .75 1]);

gz2MaxVtDist = getThresDistn(thresVals, gzMaxChirVote, chirMatchWps);
gz2MaxVtDistT2Aonly = getThresDistn(thresVals, gzMaxChirVote(longestAgree), chirMatchWps(longestAgree));
figure; hold all
plot(thresVals, gz2MaxVtDist);
plot(thresVals, gz2MaxVtDistT2Aonly);
xlabel('GZ2 min proportion designating dominant chirality');
ylabel('GZ2 match rate wtd pitch angle sum)');
legend('all', 'top 2 agree', 'Location', 'Best');
% axs = axis; axs(3:4) = [.8 1]; axis(axs);
axis([0 1 .75 1]);

% proportion of people voting for the dominant chirality:
gzChirProp = max(gzv ./ repmat(sum(gzv, 2), 1, 2), [], 2); 
gzChirProp(isnan(gzChirProp)) = 0;
gzMcvGtH = (gzMaxChirVote > 0.5);

thresVals = 0.5:0.01:1;
gz2PropDist = getThresDistn(thresVals, gzChirProp(gzMcvGtH), chirMatchAlenw(gzMcvGtH));
gz2PropDistT2Aonly = getThresDistn(thresVals, gzChirProp(gzMcvGtH & longestAgree), chirMatchAlenw(gzMcvGtH & longestAgree));
figure; hold all
plot(thresVals, gz2PropDist);
plot(thresVals, gz2PropDistT2Aonly);
title('(at least one chirality gets >= 50% of the votes)');
xlabel('GZ2 min proportion of agreeing chirality votes');
ylabel('GZ2 match rate (alenwDco)');
axis([0.5 1 .75 1]);
legend('all', 'top 2 agree', 'Location', 'Best');

thresVals = 0.5:0.01:1;
gz2PropDist = getThresDistn(thresVals, gzChirProp(gzMcvGtH), chirMatchLongest(gzMcvGtH));
gz2PropDistT2Aonly = getThresDistn(thresVals, gzChirProp(gzMcvGtH & longestAgree), chirMatchLongest(gzMcvGtH & longestAgree));
figure; hold all
plot(thresVals, gz2PropDist);
plot(thresVals, gz2PropDistT2Aonly);
title('(at least one chirality gets >= 50% of the votes)');
xlabel('GZ2 min proportion of agreeing chirality votes');
ylabel('GZ2 match rate (chirality of longest)');
axis([0.5 1 .75 1]);
legend('all', 'top 2 agree', 'Location', 'Best');

thresVals = 0.5:0.01:1;
gz2PropDist = getThresDistn(thresVals, gzChirProp(gzMcvGtH), chirMatchMaj(gzMcvGtH));
gz2PropDistT2Aonly = getThresDistn(thresVals, gzChirProp(gzMcvGtH & longestAgree), chirMatchMaj(gzMcvGtH & longestAgree));
figure; hold all
plot(thresVals, gz2PropDist);
plot(thresVals, gz2PropDistT2Aonly);
title('(at least one chirality gets >= 50% of the votes)');
xlabel('GZ2 min proportion of agreeing chirality votes');
ylabel('GZ2 match rate (majority vote)');
axis([0.5 1 .75 1]);
legend('all', 'top 2 agree', 'Location', 'Best');

thresVals = 0.5:0.01:1;
gz2PropDist = getThresDistn(thresVals, gzChirProp(gzMcvGtH), chirMatchWps(gzMcvGtH));
gz2PropDistT2Aonly = getThresDistn(thresVals, gzChirProp(gzMcvGtH & longestAgree), chirMatchWps(gzMcvGtH & longestAgree));
figure; hold all
plot(thresVals, gz2PropDist);
plot(thresVals, gz2PropDistT2Aonly);
title('(at least one chirality gets >= 50% of the votes)');
xlabel('GZ2 min proportion of agreeing chirality votes');
ylabel('GZ2 match rate (weighted pitch angle sum)');
axis([0.5 1 .75 1]);
legend('all', 'top 2 agree', 'Location', 'Best');

% figure; hist(gzChirProp(gzMcvGtH), 100);
% figure; hist(gzChirProp(gzMcvGtH & longestAgree), 100);

nAvail = size(autoTblI, 1) - 1;
nImgs = nAvail + nSkip;
fprintf('\n\n');
fprintf('Availability: \n');
fprintf('All : %2.6f (%d of %d)\n', nAvail/nImgs, nAvail, nImgs);
fprintf('longest 2 agree: %2.6f (%d of %d)\n', sum(longestAgree)/nImgs, sum(longestAgree), nImgs);
fprintf('one chirality designation gets at least 50%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.5)/nImgs, sum(gzMaxChirVote > 0.5), nImgs);
fprintf('one chirality designation gets at least 67%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.67)/nImgs, sum(gzMaxChirVote > 0.67), nImgs);
fprintf('one chirality designation gets at least 75%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.75)/nImgs, sum(gzMaxChirVote > 0.75), nImgs);
fprintf('one chirality designation gets at least 90%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.9)/nImgs, sum(gzMaxChirVote > 0.9), nImgs);
fprintf('one chirality designation gets at least 50%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.95)/nImgs, sum(gzMaxChirVote > 0.95), nImgs);
fprintf('\n');

fprintf('Inclusion rates (images with output only): \n');
fprintf('longest 2 agree: %2.6f (%d of %d)\n', sum(longestAgree)/nAvail, sum(longestAgree), nAvail);
fprintf('one chirality designation gets at least 50%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.5)/nAvail, sum(gzMaxChirVote > 0.5), nAvail);
fprintf('one chirality designation gets at least 67%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.67)/nAvail, sum(gzMaxChirVote > 0.67), nAvail);
fprintf('one chirality designation gets at least 75%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.75)/nAvail, sum(gzMaxChirVote > 0.75), nAvail);
fprintf('one chirality designation gets at least 90%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.9)/nAvail, sum(gzMaxChirVote > 0.9), nAvail);
fprintf('one chirality designation gets at least 95%% of the votes: %2.6f (%d of %d)\n', ...
    sum(gzMaxChirVote > 0.95)/nAvail, sum(gzMaxChirVote > 0.95), nAvail);
fprintf('\n');

fprintf('Inclusion rates (images with output, one chirality gets >= 50%% of the votes): \n');
fprintf('longest 2 agree: %2.6f (%d of %d)\n', sum(longestAgree & gzMcvGtH)/sum(gzMcvGtH), sum(longestAgree & gzMcvGtH), sum(gzMcvGtH));
fprintf('at least 50%% of chirality votes agree: %2.6f (%d of %d)\n', ...
    sum(gzMcvGtH & gzChirProp > 0.5)/sum(gzMcvGtH), sum(gzMcvGtH & gzChirProp > 0.5), sum(gzMcvGtH));
fprintf('at least 67%% of chirality votes agree: %2.6f (%d of %d)\n', ...
    sum(gzMcvGtH & gzChirProp > 0.67)/sum(gzMcvGtH), sum(gzMcvGtH & gzChirProp > 0.67), sum(gzMcvGtH));
fprintf('at least 75%% of chirality votes agree: %2.6f (%d of %d)\n', ...
    sum(gzMcvGtH & gzChirProp > 0.75)/sum(gzMcvGtH), sum(gzMcvGtH & gzChirProp > 0.75), sum(gzMcvGtH));
fprintf('at least 90%% of chirality votes agree: %2.6f (%d of %d)\n', ...
    sum(gzMcvGtH & gzChirProp > 0.9)/sum(gzMcvGtH), sum(gzMcvGtH & gzChirProp > 0.9), sum(gzMcvGtH));
fprintf('at least 95%% of chirality votes agree: %2.6f (%d of %d)\n', ...
    sum(gzMcvGtH & gzChirProp > 0.95)/sum(gzMcvGtH), sum(gzMcvGtH & gzChirProp > 0.95), sum(gzMcvGtH));
fprintf('\n');

fprintf('Agreement Rates\t\t\t\t (All, longestAgree)\n');
fprintf('Majority-vote:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchMaj) / nAvail, ...
    sum(chirMatchMaj(longestAgree)) / sum(longestAgree));
fprintf('Arc-length-weighted:\t\t %2.6f, %2.6f\n', sum(chirMatchAlenw) / nAvail, ...
    sum(chirMatchAlenw(longestAgree)) / sum(longestAgree));
fprintf('Weighted-pitch-angle-sum:\t %2.6f, %2.6f\n', sum(chirMatchWps) / nAvail, ...
    sum(chirMatchWps(longestAgree)) / sum(longestAgree));
fprintf('Longest-arc:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchLongest) / nAvail, ...
    sum(chirMatchLongest(longestAgree)) / sum(longestAgree));
fprintf('\n');

fprintf('Agreement Rates, counting non-output as misclassification: (All, longestAgree)\n');
fprintf('Majority-vote:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchMaj) / nImgs, ...
    sum(chirMatchMaj(longestAgree)) / (sum(longestAgree) + nSkip));
fprintf('Arc-length-weighted:\t\t %2.6f, %2.6f\n', sum(chirMatchAlenw) / nImgs, ...
    sum(chirMatchAlenw(longestAgree)) / (sum(longestAgree) + nSkip));
fprintf('Weighted-pitch-angle-sum:\t %2.6f, %2.6f\n', sum(chirMatchWps) / nImgs, ...
    sum(chirMatchWps(longestAgree)) / (sum(longestAgree) + nSkip));
fprintf('Longest-arc:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchLongest) / nImgs, ...
    sum(chirMatchLongest(longestAgree)) / (sum(longestAgree) + nSkip));
fprintf('\n');

skippedGzv = str2double(gzSkipTbl(2:end, [strmatch('P_CW', gzSkipTbl(1, :), 'exact') strmatch('P_ACW', gzSkipTbl(1, :), 'exact')]));
skippedGzMaxChirVote = max(skippedGzv, [], 2);
% hstedg = 0:0.01:1;
figure;  
subplot(1, 2, 1); hist(skippedGzMaxChirVote, 100)
subplot(1, 2, 2); hist(gzMaxChirVote, 100)

fprintf('mean maxChiralityVote for images without output available: %2.6f\n', mean(skippedGzMaxChirVote));
fprintf('mean maxChiralityVote for images with output available: %2.6f\n', mean(gzMaxChirVote));

skippedGzChirProp = max(skippedGzv ./ repmat(sum(skippedGzv, 2), 1, 2), [], 2); 
skippedGzChirProp(isnan(skippedGzChirProp)) = 0;

figure;
subplot(1, 2, 1); hist(skippedGzChirProp, 100)
subplot(1, 2, 2); hist(gzChirProp, 100)

fprintf('proportion of agreeing chirality votes, images without output available: %2.6f\n', ...
    mean(skippedGzChirProp));
fprintf('proportion of agreeing chirality votes, images with output available: %2.6f\n', ...
    mean(gzChirProp));
fprintf('\n\n');

dscnVals = [0 0.6 0.8 0.9 0.95 1];
dscnVals = sort([0:0.1:1 .95], 'ascend');
dscnVals = sort([0 0.25 0.5:0.1:1 .95], 'ascend');
fh1 = figure; ah1 = gca; hold all
fh2 = figure; ah2 = gca; hold all
fprintf('\\begin{table*}[t]\n');
fprintf('\\begin{center}\n');
fprintf('\\setlength{\\tabcolsep}{5pt}\n');
fprintf('\\begin{tabular}{|l|%s|%s|}\n', repmat('c ', 1, length(dscnVals)), repmat('c ', 1, length(dscnVals)));
fprintf('\\hline\n');
fprintf('\\textbf{Min Discernibility}');
% fprintf(',%2.1f', [dscnVals(:)'; dscnVals(:)']); fprintf('\n');
fprintf(' & %d', 100 * repmat(dscnVals, 1, 2)); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Require Longest Agree}');
% dispLongestAgree = repmat({'F,', 'T,'}, 1, length(dscnVals));
% fprintf('%s', [dispLongestAgree{:}]); fprintf('\n');
fprintf('%s', [repmat(' & F', 1, length(dscnVals)), repmat(' & T', 1, length(dscnVals))]); fprintf(' \\\\ \\hline \\hline\n');
availDist = repmat(gzMaxChirVote, 1, 2*length(dscnVals)) >= repmat(dscnVals, size(gzMaxChirVote, 1), 2);
availDist = availDist & [true(size(gzMaxChirVote, 1), length(dscnVals)) repmat(longestAgree, 1, length(dscnVals))];
availDist = sum(availDist, 1) / (size(gzMaxChirVote, 1)+nSkip);
fprintf('\\textbf{Inclusion Rate}');
fprintf(' & %2.1f', 100 * availDist); fprintf(' \\\\ \\hline\n');
plot(ah1, dscnVals, availDist(1:length(dscnVals)), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
plot(ah2, dscnVals, availDist(length(dscnVals)+1:end), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', ':');
avgDscn1 = getThresAvgs(dscnVals, gzMaxChirVote, gzMaxChirVote);
avgDscn2 = getThresAvgs(dscnVals, gzMaxChirVote(longestAgree), gzMaxChirVote(longestAgree));
fprintf('\\textbf{Mean Discernibility}');
fprintf(' & %2.1f', 100* [avgDscn1(:)' avgDscn2(:)']); fprintf(' \\\\ \\hline \\hline\n');
fprintf('\\textbf{Majority Vote}');
distnMaj1 = getThresDistn(dscnVals, gzMaxChirVote, chirMatchMaj);
distnMaj2 = getThresDistn(dscnVals, gzMaxChirVote(longestAgree), chirMatchMaj(longestAgree));
% fprintf(',%2.1f', 100 * [distn1(:)'; distn2(:)']); fprintf('\n');
fprintf(' & %2.1f', 100 * [distnMaj1(:)' distnMaj2(:)']); fprintf(' \\\\ \\hline\n');
plot(ah1, dscnVals, distnMaj1, 'LineWidth', 1.5);
plot(ah2, dscnVals, distnMaj2, 'LineWidth', 1.5);
fprintf('\\textbf{Longest Arc}');
distnLongest1 = getThresDistn(dscnVals, gzMaxChirVote, chirMatchLongest);
distnLongest2 = getThresDistn(dscnVals, gzMaxChirVote(longestAgree), chirMatchLongest(longestAgree));
% fprintf(',%2.1f', 100 * [distn1(:)'; distn2(:)']); fprintf('\n');
fprintf(' & %2.1f', 100 * [distnLongest1(:)' distnLongest2(:)']); fprintf(' \\\\ \\hline\n');
plot(ah1, dscnVals, distnLongest1, 'LineWidth', 1.5);
plot(ah2, dscnVals, distnLongest2, 'LineWidth', 1.5);
fprintf('\\textbf{Length-weighted Vote}');
distnAlenw1 = getThresDistn(dscnVals, gzMaxChirVote, chirMatchAlenw);
distnAlenw2 = getThresDistn(dscnVals, gzMaxChirVote(longestAgree), chirMatchAlenw(longestAgree));
% fprintf(',%2.1f', 100 * [distn1(:)'; distn2(:)']); fprintf('\n');
fprintf(' & %2.1f', 100 * [distnAlenw1(:)' distnAlenw2(:)']); fprintf(' \\\\ \\hline\n');
plot(ah1, dscnVals, distnAlenw1, 'LineWidth', 1.5);
plot(ah2, dscnVals, distnAlenw2, 'LineWidth', 1.5);
fprintf('\\end{tabular}\n');
fprintf('\\end{center}\n');
fprintf('\\caption{}\n');
fprintf('\\label{tbl:ChiralityAgreementAll}\n');
fprintf('\\end{table*}\n');

fprintf('\n\n');
fprintf('not requiring agreement of two longest arcs:\n');

% only print the part where longest arcs aren't required to agree
fprintf('\\begin{table*}\n');
fprintf('\\begin{center}\n');
fprintf('\\setlength{\\tabcolsep}{5pt}\n');
fprintf('\\begin{tabular}{|l|%s|}\n', repmat('c ', 1, length(dscnVals)));
fprintf('\\hline\n');
fprintf('\\textbf{Min Discernibility}');
fprintf(' & %d', 100 * dscnVals); fprintf(' \\\\ \\hline \\hline\n');
fprintf('\\textbf{Inclusion Rate}');
fprintf(' & %2.1f', 100 * availDist(1:length(dscnVals))); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Mean Discernibility}');
fprintf(' & %2.1f', 100 * avgDscn1); fprintf(' \\\\ \\hline \\hline\n');
fprintf('\\textbf{Majority Vote}');
fprintf(' & %2.1f', 100 * distnMaj1); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Longest Arc}');
fprintf(' & %2.1f', 100 * distnLongest1); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Length-weighted Vote}');
fprintf(' & %2.1f', 100 * distnAlenw1); fprintf(' \\\\ \\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{center}\n');
fprintf('\\caption{}\n');
fprintf('\\label{tbl:ChiralityAgreementNoRequireLongestAgree}\n');
fprintf('\\end{table*}\n');

fprintf('\n\n');
fprintf('requiring agreement of two longest arcs:\n');

% only print the part where longest arcs are required to agree
fprintf('\\begin{table*}\n');
fprintf('\\begin{center}\n');
fprintf('\\setlength{\\tabcolsep}{5pt}\n');
fprintf('\\begin{tabular}{|l|%s|}\n', repmat('c ', 1, length(dscnVals)));
fprintf('\\hline\n');
fprintf('\\textbf{Min Discernibility}');
fprintf(' & %d', 100 * dscnVals); fprintf(' \\\\ \\hline \\hline\n');
fprintf('\\textbf{Inclusion Rate}');
fprintf(' & %2.1f', 100 * availDist(length(dscnVals)+1:end)); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Mean Discernibility}');
fprintf(' & %2.1f', 100 * avgDscn2); fprintf(' \\\\ \\hline \\hline\n');
fprintf('\\textbf{Majority Vote}');
fprintf(' & %2.1f', 100 * distnMaj2); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Longest Arc}');
fprintf(' & %2.1f', 100 * distnLongest2); fprintf(' \\\\ \\hline\n');
fprintf('\\textbf{Length-weighted Vote}');
fprintf(' & %2.1f', 100 * distnAlenw2); fprintf(' \\\\ \\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{center}\n');
fprintf('\\caption{}\n');
fprintf('\\label{tbl:ChiralityAgreementRequireLongestAgree}\n');
fprintf('\\end{table*}\n');

fprintf('\n\n');

isHighGZ2ConfidenceDisagreement = ~chirMatchAlenw & (gzMaxChirVote == 1);
highGZ2ConfidenceDisagreementInstances = names(isHighGZ2ConfidenceDisagreement & longestAgree);
fprintf('\nDisagreements with high GZ2 confidence, longest agree: \n');
for jj=1:length(highGZ2ConfidenceDisagreementInstances)
    fprintf('\t%s\n', highGZ2ConfidenceDisagreementInstances{jj});
end
fprintf('number of galaxies with high GZ2 confidence and longest arcs agreeing: %d\n', sum((gzMaxChirVote == 1) & longestAgree));
fprintf('\n');


figure(fh1); 
% hTitle = title({'Winding Direction Agreement:', 'Galaxy Zoo'});
hXlbl = xlabel('min discernibility rate');
hYlbl = ylabel('inclusion/agreement rate');
hLgd = legend('inclusion rate', 'majority vote', 'longest arc', 'length-weighted vote', 'Location', 'SouthWest');
set(gcf, 'Color', [1 1 1]);
set(gca, 'YTick', [0:0.1:1]);
psn = get(gcf, 'Position'); psn(3:4) = [650 300]; set(gcf, 'Position', psn);
% set(hTitle, 'FontSize', 18)
set(hXlbl, 'FontSize', 18)
set(hYlbl, 'FontSize', 18);
set(hLgd, 'FontSize', 18);
set(gca, 'FontSize', 12);
% set(hTitle, 'FontName', 'Arial');
set(hXlbl, 'FontName', 'Arial');
set(hYlbl, 'FontName', 'Arial');
set(hLgd, 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
set(gcf, 'PaperPositionMode', 'auto');


figure(fh2);
% hTitle = title({'Winding Direction Agreement', 'Galaxy Zoo, 2-longest-agree only'});
hXlbl = xlabel('min discernibility rate');
hYlbl = ylabel('inclusion/agreement rate');
hLgd = legend('inclusion rate', 'majority vote', 'longest arc', 'length-weighted vote', 'Location', 'SouthWest');
set(gcf, 'Color', [1 1 1]);
set(gca, 'YTick', [0:0.1:1]);
psn = get(gcf, 'Position'); psn(3:4) = [650 300]; set(gcf, 'Position', psn);
% set(hTitle, 'FontSize', 18)
set(hXlbl, 'FontSize', 18)
set(hYlbl, 'FontSize', 18);
set(hLgd, 'FontSize', 18);
set(gca, 'FontSize', 12);
% set(hTitle, 'FontName', 'Arial');
set(hXlbl, 'FontName', 'Arial');
set(hYlbl, 'FontName', 'Arial');
set(hLgd, 'FontName', 'Arial');
set(gca, 'FontName', 'Arial');
set(gcf, 'PaperPositionMode', 'auto');

chirBipGz = -chirGz * 2 + 1;

chirBipMaj = NaN * zeros(size(chirMaj));
chirBipMaj(strcmp('S-wise', chirMaj)) = -1;
chirBipMaj(strcmp('Z-wise', chirMaj)) = 1;
chirBipMaj(strcmp('EQ', chirMaj)) = 0;
print_chir_bias_stats(chirBipGz, chirBipMaj, 'majority-vote');

chirBipAlenw = NaN * zeros(size(chirAlenw));
chirBipAlenw(strcmp('S-wise', chirAlenw)) = -1;
chirBipAlenw(strcmp('Z-wise', chirAlenw)) = 1;
chirBipAlenw(strcmp('EQ', chirAlenw)) = 0;
print_chir_bias_stats(chirBipGz, chirBipAlenw, 'arc-length-weighted vote');

print_chir_bias_stats(chirBipGz(longestAgree), chirBipMaj(longestAgree), ...
    'majority-vote, SUBSET where longest two arcs agree');

print_chir_bias_stats(chirBipGz(longestAgree), chirBipAlenw(longestAgree), ...
    'arc-length-weighted vote, SUBSET where longest two arcs agree');

function distn = getThresDistn(thresVals, thresMeas, isMatch)
    thresMeas(isnan(thresMeas)) = -1;
    distn = zeros(size(thresVals));
    for ii=1:1:length(thresVals)
        withinThres = thresMeas >= thresVals(ii);
        distn(ii) = sum(isMatch(withinThres)) / sum(withinThres);
    end
end

function avgs = getThresAvgs(thresVals, thresMeas, vals)
    assert(all(~isnan(thresMeas)))
    avgs = zeros(size(thresVals));
    for ii=1:1:length(thresVals)
        withinThres = thresMeas >= thresVals(ii);
        avgs(ii) = mean(vals(withinThres));
    end
end

function print_agreement_vs_discernability(dscn, isChirMatch, name)
    [dscn, sI] = sort(dscn, 'descend');
    isChirMatch = isChirMatch(sI);
    avg_dscn = cumsum(dscn) ./ (1:length(dscn))';
    avg_agree = cumsum(isChirMatch) ./ (1:length(isChirMatch))';
    inclusion_rate = (1:length(dscn))' / length(dscn);
    agree_gt_dscn = avg_agree > avg_dscn;
    start_idx = find(agree_gt_dscn, 1, 'first');
    end_idx = find(~agree_gt_dscn, 1, 'last');
    fprintf('%s: \n', name);
    fprintf('first point where agreement > discernability: agreement = %f avg_dscn = %f min_dscn = %f incl = %d (%f%%)\n',...
        avg_agree(start_idx), avg_dscn(start_idx), dscn(start_idx), start_idx, 100 * start_idx / length(dscn));
    fprintf('last point where agreement <= discernability: agreement = %f avg_dscn = %f min_dscn = %f incl = %d (%f%%)\n',...
        avg_agree(end_idx), avg_dscn(end_idx), dscn(end_idx), end_idx, 100 * end_idx / length(dscn));
    
    crossover_x = end_idx / length(dscn);
    figure; hold all
    line(crossover_x * [1 1], [0, 1], 'Color', 'k');
    plot(inclusion_rate(end:-1:1), dscn(end:-1:1), 'm:', 'LineWidth', 1.5);
    plot(inclusion_rate(end:-1:1), avg_dscn(end:-1:1), 'm-', 'LineWidth', 1.5);
    plot(inclusion_rate(end:-1:1), avg_agree(end:-1:1), 'b-', 'LineWidth', 1.5);
    xlabel('inclusion rate');
    figure('Position', [100 100 650 350]); 
    hold all
    line(crossover_x * [1 1], [0, 1], 'Color', 'k');
%     plot(inclusion_rate(end:-1:1), dscn(end:-1:1), 'g:', 'LineWidth', 1.5);
    plot(inclusion_rate(end:-1:1), avg_dscn(end:-1:1), 'g-', 'LineWidth', 1.5);
    plot(inclusion_rate(end:-1:1), avg_agree(end:-1:1), 'b-', 'LineWidth', 1.5); 
    hXlbl = xlabel('Inclusion Rate');
    hYlbl = ylabel('Agreement and Discernibility');
    ylim([0.8 1.005]);
    set(gcf, 'Color', [1 1 1]);
    set(hXlbl, 'FontSize', 18);
    set(hYlbl, 'FontSize', 18);
    set(gca, 'FontSize', 16);
    
%     figure; hold all
%     plot(dscn(end:-1:1), inclusion_rate(end:-1:1));
%     plot(dscn(end:-1:1), avg_dscn(end:-1:1));
%     plot(dscn(end:-1:1), avg_agree(end:-1:1));
end

% agrBinGran = 0.01;
% agrBins = 0:agrBinGran:1;
% agrBinIdxs = gzMaxChirVote / agrBinGran - (agrBins(1) / agrBinGran);
% agrBinIdxs = floor(agrBinIdxs) + 1;
% 
% gz2MaxVtDist = accumarray(agrBinIdxs, ones(size(gzMaxChirVote)), [length(agrBins) 1]);
% gz2MatchDist = accumarray(agrBinIdxs, chirMatchAlenw, [length(agrBins) 1]);
% figure; plot(agrBins, cumsum(gz2MatchDist) ./ cumsum(gz2MaxVtDist));
% xlabel('GZ2 max proportion designating dominant chirality'); ylabel('GZ2 match rate');
% figure; hold all
% plot([0 1], [0 1]);
% matchPropDist = cumsum(gz2MatchDist(end:-1:1)) ./ cumsum(gz2MaxVtDist(end:-1:1));
% plot(agrBins(end:-1:1), matchPropDist);
% xlabel('GZ2 min proportion designating dominant chirality'); ylabel('GZ2 match rate');
% axis([0 1 min(matchPropDist) 1])
% 
% propBinIdxs = gzChirProp / agrBinGran - (agrBins(1) / agrBinGran);
% propBinIdxs = floor(propBinIdxs) + 1;
% 
% gzPropDist = accumarray(propBinIdxs(gzMaxGtHalf), ones(1, sum(gzMaxGtHalf)), [length(agrBins) 1]);
% gzMatchGthDist = accumarray(propBinIdxs(gzMaxGtHalf), chirMatchAlenw(gzMaxGtHalf), [length(agrBins) 1]);
% matchGthDist = cumsum(gzMatchGthDist(end:-1:1) ./ gzPropDist(end:-1:1));
% plot(agrBins(end:-1:1), matchGthDist);
% title('GZ2 agreement for galaxies where >= half of votes in dominant-chirality category');
% xlabel('GZ2 proportion dominant chirality vs both chirality'); ylabel('GZ2 match rate');
% 
% figure; plot(matchGthDist);

end
