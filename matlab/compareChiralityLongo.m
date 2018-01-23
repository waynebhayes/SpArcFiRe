% set resultCsvPath and cmpCsvPath before running

% resultCsvPath = 'D:\matlab-output-local\[paper] 2011-11-07 GZ2 secondRun\GZ2_full_secondRun.csv';
% cmpCsvPath = 'D:\Galaxy-Sets\Longo\LongoClassifications.csv';

[autoTbl, nSkip, skippedNames] = readTableFromFile(resultCsvPath, ',');
hTbl = readTableFromFile(cmpCsvPath, ',');
[skhTbl, tmp] = intersectRowsByField(hTbl, [{'name'}; skippedNames], 1, 1, true);
[autoTblI, hTblI] = intersectRowsByField(autoTbl, hTbl, 1, 1, true);
cmpTbl = joinTables(autoTblI, hTblI, 1, 1);
nLongo = size(hTbl, 1)-1;
nLongoSkip = size(skhTbl, 1) - 1;
nNonOverlap = size(hTbl, 1) - size(hTblI, 1) - nLongoSkip;
nCmp = size(cmpTbl, 1) - 1;
nNoOutput = length(skippedNames);

longoChirIdx = strmatch('HAND', cmpTbl(1, :), 'exact');
longoChir = cmpTbl(2:end, longoChirIdx);

cMajIdx = strmatch('chirality_maj', cmpTbl(1, :), 'exact');
cAlenwIdx = strmatch('chirality_alenWtd', cmpTbl(1, :), 'exact');
cWpsIdx = strmatch('chirality_wtdPangSum', cmpTbl(1, :), 'exact');
cLongestIdx = strmatch('chirality_longestArc', cmpTbl(1, :), 'exact');

chirMaj = cmpTbl(2:end, cMajIdx);
chirAlenw = cmpTbl(2:end, cAlenwIdx);
chirWps = cmpTbl(2:end, cWpsIdx);
chirLongest = cmpTbl(2:end, cLongestIdx);

chirMatchMaj = (strcmp('S-wise', chirMaj) & strcmp('R', longoChir)) | (strcmp('Z-wise', chirMaj) & strcmp('L', longoChir));
chirMatchAlenw = (strcmp('S-wise', chirAlenw) & strcmp('R', longoChir)) | (strcmp('Z-wise', chirAlenw) & strcmp('L', longoChir));
chirMatchWps = (strcmp('S-wise', chirWps) & strcmp('R', longoChir)) | (strcmp('Z-wise', chirWps) & strcmp('L', longoChir));
chirMatchLongest = (strcmp('S-wise', chirLongest) & strcmp('R', longoChir)) | (strcmp('Z-wise', chirLongest) & strcmp('L', longoChir));

t2aIdx = strmatch('top2_chirality_agreement', cmpTbl(1, :), 'exact');
longestAgree = cellfun(@(x)(strcmp(x(2:end-1), 'agree')), cmpTbl(2:end, t2aIdx));

fprintf('%d GZ2 sample galaxies without output\n', nNoOutput);
fprintf('\nOf the galaxies in the Longo table...\n');
fprintf('\t%02.4f%% (%d of %d) are included for comparison\n', 100 * nCmp/nLongo, nCmp, nLongo);
fprintf('\t%02.4f%% (%d of %d) do not overlap with the galaxies we have from GZ2\n', 100 * nNonOverlap/nLongo, nNonOverlap, nLongo);
fprintf('\t%02.4f%% (%d of %d) did not have output from our arc-finding\n', 100 * nLongoSkip/nLongo, nLongoSkip, nLongo);
fprintf('\n');

fprintf('Agreement Rates\t\t\t\t (All, longestAgree)\n');
fprintf('Majority-vote:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchMaj) / nCmp, ...
    sum(chirMatchMaj(longestAgree)) / sum(longestAgree));
fprintf('Arc-length-weighted:\t\t %2.6f, %2.6f\n', sum(chirMatchAlenw) / nCmp, ...
    sum(chirMatchAlenw(longestAgree)) / sum(longestAgree));
fprintf('Weighted-pitch-angle-sum:\t %2.6f, %2.6f\n', sum(chirMatchWps) / nCmp, ...
    sum(chirMatchWps(longestAgree)) / sum(longestAgree));
fprintf('Longest-arc:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchLongest) / nCmp, ...
    sum(chirMatchLongest(longestAgree)) / sum(longestAgree));
fprintf('\n');

fprintf('Agreement Rates, counting non-output as misclassification: (All, longestAgree)\n');
fprintf('Majority-vote:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchMaj) / (nCmp + nLongoSkip), ...
    sum(chirMatchMaj(longestAgree)) / (sum(longestAgree) + nLongoSkip));
fprintf('Arc-length-weighted:\t\t %2.6f, %2.6f\n', sum(chirMatchAlenw) / (nCmp + nLongoSkip), ...
    sum(chirMatchAlenw(longestAgree)) / (sum(longestAgree) + nLongoSkip));
fprintf('Weighted-pitch-angle-sum:\t %2.6f, %2.6f\n', sum(chirMatchWps) / (nCmp + nLongoSkip), ...
    sum(chirMatchWps(longestAgree)) / (sum(longestAgree) + nLongoSkip));
fprintf('Longest-arc:\t\t\t\t %2.6f, %2.6f\n', sum(chirMatchLongest) / (nCmp + nLongoSkip), ...
    sum(chirMatchLongest(longestAgree)) / (sum(longestAgree) + nLongoSkip));
fprintf('\n');

fprintf('%2.6f%% of galaxies have agreement between the longest arcs\n', 100*sum(longestAgree)/numel(longestAgree));

fprintf('\\begin{table}[htb]\n');
fprintf('\\begin{tabular}{|l|c|c|}\n');
fprintf('\\hline\n');
fprintf('\\textbf{} & \\multicolumn{1}{l|}{\\textbf{All}} & \\multicolumn{1}{l|}{\\textbf{Longest Agree}} \\\\ \\hline\n');
fprintf('\\textbf{Majority Vote} & %2.1f & %2.1f \\\\ \\hline\n', 100 * sum(chirMatchMaj) / nCmp, ...
    100 * sum(chirMatchMaj(longestAgree)) / sum(longestAgree));
fprintf('\\textbf{Longest Arc} & %2.1f & %2.1f \\\\ \\hline\n', 100 * sum(chirMatchLongest) / nCmp, ...
    100 * sum(chirMatchLongest(longestAgree)) / sum(longestAgree));
fprintf('\\textbf{Length-weighted Vote} & %2.1f & %2.1f \\\\ \\hline\n', 100 * sum(chirMatchAlenw) / nCmp, ...
    100 * sum(chirMatchAlenw(longestAgree)) / sum(longestAgree));
fprintf('\\end{tabular}\n');
fprintf('\\label{tbl:LongoChirality}\n');
fprintf('\\end{table}\n');
fprintf('\n');

chirBipLongo = zeros(nCmp, 1); 
chirBipLongo(strcmp('R', longoChir)) = -1;
chirBipLongo(strcmp('L', longoChir)) = 1;
assert(all(abs(chirBipLongo) == 1));

chirBipMaj = NaN * zeros(size(chirMaj));
chirBipMaj(strcmp('S-wise', chirMaj)) = -1;
chirBipMaj(strcmp('Z-wise', chirMaj)) = 1;
chirBipMaj(strcmp('EQ', chirMaj)) = 0;
print_chir_bias_stats(chirBipLongo, chirBipMaj, 'majority-vote');

chirBipAlenw = NaN * zeros(size(chirAlenw));
chirBipAlenw(strcmp('S-wise', chirAlenw)) = -1;
chirBipAlenw(strcmp('Z-wise', chirAlenw)) = 1;
chirBipAlenw(strcmp('EQ', chirAlenw)) = 0;
print_chir_bias_stats(chirBipLongo, chirBipAlenw, 'arc-length-weighted vote');

print_chir_bias_stats(chirBipLongo(longestAgree), chirBipMaj(longestAgree), ...
    'majority-vote, SUBSET where longest two arcs agree');

print_chir_bias_stats(chirBipLongo(longestAgree), chirBipAlenw(longestAgree), ...
    'arc-length-weighted vote, SUBSET where longest two arcs agree');
