% set resultCsvPath and cmpCsvPath before running
% cmpCsvPath = 'D:\Galaxy-Sets\GalaxyZoo1_DR_table2.csv'

[autoTbl, nSkip] = readTableFromFile(resultCsvPath, ',');
gzTbl = readTableFromFile(cmpCsvPath, ',');
[autoTblI, gzTblI] = intersectRowsByField(autoTbl, gzTbl, 1, 1, true);
clear gzTbl
if size(autoTbl, 1) ~= size(autoTblI, 1)
    warning('some auto-generated entries not present in comparison file');
    pause
end
cmpTbl = joinTables(autoTblI, gzTblI, 1, 1);

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

cbr = str2double(cmpTbl(2:end, cbrIdx));
cbrOK = (cbr > 1);
longestAgree = cellfun(@(x)(strcmp(x(2:end-1), 'agree')), cmpTbl(2:end, t2aIdx));

gzv = str2double(cmpTbl(2:end, [pcwIdx pacwIdx]));
chirGz = gzv(:, 2) > gzv(:, 1);
chirMaj = cmpTbl(2:end, cMajIdx);
chirAlenw = cmpTbl(2:end, cAlenwIdx);
chirWps = cmpTbl(2:end, cWpsIdx);
chirLongest = cmpTbl(2:end, cLongestIdx);

chirMatchMaj = (strcmp('S-wise', chirMaj) & chirGz) | (strcmp('Z-wise', chirMaj) & ~chirGz);
chirMatchAlenw = (strcmp('S-wise', chirAlenw) & chirGz) | (strcmp('Z-wise', chirAlenw) & ~chirGz);
chirMatchWps = (strcmp('S-wise', chirWps) & chirGz) | (strcmp('Z-wise', chirWps) & ~chirGz);
chirMatchLongest = (strcmp('S-wise', chirLongest) & chirGz) | (strcmp('Z-wise', chirLongest) & ~chirGz);

nAvail = size(autoTblI, 1) - 1;
nImgs = nAvail + nSkip;
fprintf('\n\n');
fprintf('Availability: \n');
fprintf('All : %2.6f (%d of %d)\n', nAvail/nImgs, nAvail, nImgs);
nAvailCbr = sum(cbrOK);
fprintf('ContourBrtRatio > 1 : %2.6f (%d of %d)\n', nAvailCbr/nImgs, nAvailCbr, nImgs);
nAvailCbrAg = sum(cbrOK & longestAgree);
fprintf('CBR > 1 and longest 2 agree : %2.6f (%d of %d)\n', nAvailCbrAg/nImgs, nAvailCbrAg, nImgs);
fprintf('\n');

fprintf('Agreement Rates (All, CBR > 1, CBR > 1 & longestAgree)\n');
fprintf('Majority-vote: %2.6f, %2.6f, %2.6f\n', sum(chirMatchMaj) / nAvail, ...
    sum(chirMatchMaj(cbrOK)) / nAvailCbr, sum(chirMatchMaj(cbrOK & longestAgree)) / nAvailCbrAg);
fprintf('Arc-length-weighted: %2.6f, %2.6f, %2.6f\n', sum(chirMatchAlenw) / nAvail, ...
    sum(chirMatchAlenw(cbrOK)) / nAvailCbr, sum(chirMatchAlenw(cbrOK & longestAgree)) / nAvailCbrAg);
fprintf('Weighted-pitch-angle-sum: %2.6f, %2.6f, %2.6f\n', sum(chirMatchWps) / nAvail, ...
    sum(chirMatchWps(cbrOK)) / nAvailCbr, sum(chirMatchWps(cbrOK & longestAgree)) / nAvailCbrAg);
fprintf('Longest-arc: %2.6f, %2.6f, %2.6f\n', sum(chirMatchLongest) / nAvail, ...
    sum(chirMatchLongest(cbrOK)) / nAvailCbr, sum(chirMatchLongest(cbrOK & longestAgree)) / nAvailCbrAg);

