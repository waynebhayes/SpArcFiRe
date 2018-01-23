
resultDir = 'D:\matlab-output-local\2011-11-07 GZ2 secondRun\';
outDir = 'D:\matlab-output-local\GZ2sample\';
resultCsvPath = [resultDir 'GZ2_full_secondRun.csv'];
cmpCsvPath = 'D:\Galaxy-Sets\GalaxyZoo1_DR_table2.csv';

% [autoTbl, nSkip, skippedNames] = readTableFromFile(resultCsvPath, ',');
% gzTbl = readTableFromFile(cmpCsvPath, ',');
% skippedNames = [{'name'}; skippedNames];
% [gzSkipTbl, tmp] = intersectRowsByField(gzTbl, skippedNames, 1, 1, true);
% [autoTblI, gzTblI] = intersectRowsByField(autoTbl, gzTbl, 1, 1, true);
% clear gzTbl
% if size(autoTbl, 1) ~= size(autoTblI, 1)
%     warning('some auto-generated entries not present in comparison file');
%     pause
% end
% cmpTbl = joinTables(autoTblI, gzTblI, 1, 1);

names = cmpTbl(2:end, 1);
pcwIdx = strmatch('P_CW', cmpTbl(1, :), 'exact');
pacwIdx = strmatch('P_ACW', cmpTbl(1, :), 'exact');
gzv = str2double(cmpTbl(2:end, [pcwIdx pacwIdx]));
gzMaxChirVote = max(gzv, [], 2);

splNames = randsample(names(gzMaxChirVote >= 0.9), 100);

subDirNames = cell(1, 6);
for ii=1:1:6
    fnames = dir([resultDir int2str(ii) filesep]);
    fnames = {fnames.name};
    fnames = regexp(fnames, '\d{5,}', 'match');
    fnames = fnames(~cellfun(@isempty, fnames));
    fnames = cellfun(@(x)(x{1}), fnames, 'UniformOutput', false);
    fnames = unique(fnames);
    subDirNames{ii} = fnames;
end

iptExt = '-A_input.png';
stdExt = '-B_autoCrop.png';
ovlExt = '-H_logSpiralArcs-merged.png';
for ii=1:1:length(splNames)
    curDir = '';
    for jj=1:1:6
%         tstDir = [resultDir int2str(jj) filesep];
%         tstPath = [tstDir names{jj} iptExt];
        if ismember(names{ii}, subDirNames{jj})
            curDir =  [resultDir int2str(jj) filesep];
        end
    end
    iptImg = imread([curDir names{ii} iptExt]);
    imwrite(iptImg, [outDir names{ii} '_01_input.png']);
    stdImg = imread([curDir names{ii} stdExt]);
    imwrite(stdImg, [outDir names{ii} '_02_standardized.png']);
    ovlImg = imread([curDir names{ii} ovlExt]);
    imwrite(ovlImg, [outDir names{ii} '_03_arc-overlay.png']);
end