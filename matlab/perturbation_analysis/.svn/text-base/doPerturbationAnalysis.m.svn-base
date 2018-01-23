function doPerturbationAnalysis(dirPath, baselineName, paramName, paramVals)

if dirPath(end) ~= filesep
    dirPath = [dirPath filesep];
end

outDirName = 'perturbation_analysis_out';
outDirPath = [dirPath outDirName filesep];
if ~exist(outDirPath, 'dir')
    mkdir(outDirPath);
end

fnames = dir(dirPath);
fnames = {fnames.name};
fnames = fnames(~cellfun(@isempty, regexp(fnames, 'csv$', 'match')));

bslIdx = strmatch(baselineName, fnames);
if isempty(bslIdx)
    error('no files matching the given baseline name');
elseif ~isscalar(bslIdx)
    error('multiple files matching the given baseline name');
end

nTbl = length(fnames);
gxyTbls = cell(size(fnames));
skipCounts = zeros(size(fnames));
allSkippedNames = {};
for ii=1:1:nTbl
    [curTbl, nSkip, skippedNames, wasSkipped] = ...
        readTableFromFile([dirPath fnames{ii}], ',');
    if ii == 1
        skippedInAny = wasSkipped;
    elseif any(size(wasSkipped) ~= size(skippedInAny))
        error('all tables should have the same number of rows');
    else
        skippedInAny = (skippedInAny | wasSkipped);
    end
    allSkippedNames = [allSkippedNames skippedNames'];
    
    skipCounts(ii) = nSkip;
    assert(~isempty(findstr('name', curTbl{1})), ...
        'first column of each csv should be the object name');
    curTbl(2:end, :) = sortrows(curTbl(2:end, :), 1);
    gxyTbls{ii} = curTbl;
end
fprintf('read %d files\n', nTbl);

% instances skipped in some runs need to be removed from all runs
unqSkippedNames = unique(allSkippedNames);
for ii=1:1:length(gxyTbls)
    curTbl = gxyTbls{ii};
    gxyTbls{ii} = curTbl([1; find(~ismember(curTbl(2:end, 1), unqSkippedNames))+1], :);
end

gxyNames = gxyTbls{1};
gxyNames = gxyNames(2:end, 1);
for ii=2:1:nTbl
    curNames = gxyTbls{ii};
    curNames = curNames(2:end, 1);
    assert(all(cellfun(@strcmp, gxyNames, curNames)));
end

allSkippedNames = sort(allSkippedNames);
try
    outFileForSkipped = fopen([outDirPath 'skippedCases.csv'], 'wt');
    fprintf(outFileForSkipped, 'name, skipCount\n');
    curName = '';
    curCount = 0;
    for ii=1:1:length(allSkippedNames)
        if ~isempty(strmatch(curName, allSkippedNames{ii}, 'exact'))
            curCount = curCount + 1;
        else
            if ii > 1
                fprintf(outFileForSkipped, '%s, %d\n', curName, curCount);
            end
            curName = allSkippedNames{ii};
            curCount = 1;
        end
    end
    if ii > 1
        fprintf(outFileForSkipped, '%s, %d\n', curName, curCount);
    end
    fclose(outFileForSkipped);
catch MErr
    MErr.stack
    warning('doPerturbationAnalysis:failedToWriteSkippedCasesFile', ...
        'could not write the skipped-cases file (%s)', MErr.message);
    fclose(outFileForSkipped);
end

nInst = size(gxyTbls{1}, 1) - 1; 

paAlenwDcoVals = getAllTableColumnsByName(...
    gxyTbls, 'pa_alenWtd_avg_domChiralityOnly', nTbl, nInst);
paAlenwDcoVals = str2double(paAlenwDcoVals);
analyzeDiffs('paAlenwDco', paAlenwDcoVals, bslIdx, paramName, paramVals, gxyNames, outDirPath);

vals = getAllTableColumnsByName(gxyTbls, 'pa_longest', nTbl, nInst);
vals = str2double(vals);
analyzeDiffs('paLongest', vals, bslIdx, paramName, paramVals, gxyNames, outDirPath);

end

function col = getTableColumnByName(tbl, name)
    colIdx = strmatch(name, tbl(1, :), 'exact');
    if isempty(colIdx)
        error('could not find a column with name "%s"', name);
    elseif ~isscalar(colIdx)
        error('multiple columns matching the name "%s"', name);
    end
    col = tbl(2:end, colIdx);
end

function vals = getAllTableColumnsByName(tbls, name, nTbl, nInst)
    vals = cell(nInst, nTbl);
    for ii=1:1:nTbl
        vals(:, ii) = getTableColumnByName(tbls{ii}, name);
    end
end

function analyzeDiffs(quantityName, vals, bslIdx, paramName, paramVals, gxyNames, outDirPath)
    bslVals = vals(:, bslIdx);
    diffFromBsl = vals - repmat(bslVals, 1, size(vals, 2));
%     figure; boxplot(diffFromBsl)
    
    qVals = quantile(diffFromBsl, [0.25 0.5 0.75]);
    figure; hold all
    for ii=1:1:size(qVals, 1)
        plot(qVals(ii, :));
    end
    xlabel(paramName)
    set(gca, 'XTick', 1:size(qVals, 2))
    set(gca, 'XTickLabel', paramVals)
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [100 100 1100 750]);
    ylabel(['difference in ' quantityName]);
    legend('1^{st} quartile', 'median', '3^{rd} quartile', 'Location', 'North');
    export_fig(gcf, [outDirPath paramName '_diff_in_' quantityName '.png']);
    
    qVals = quantile(abs(diffFromBsl), [0.25 0.5 0.75]);
    figure; hold all
    for ii=1:1:size(qVals, 1)
        plot(qVals(ii, :));
    end
    xlabel(paramName)
    set(gca, 'XTick', 1:size(qVals, 2))
    set(gca, 'XTickLabel', paramVals)
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [100 100 1100 750]);
    ylabel(['absolute difference in ' quantityName]);
    legend('1^{st} quartile', 'median', '3^{rd} quartile', 'Location', 'North');
    export_fig(gcf, [outDirPath paramName '_abs-diff_in_' quantityName '.png']);
    
    perGxySds = std(vals, [], 2);
    [perGxySds, sI] = sort(perGxySds, 'ascend');
    sdSortedGxyNames = gxyNames(sI);
    try
        outFile = fopen([outDirPath paramName '_perGxyStability_' quantityName '.csv'], 'wt');
        fprintf(outFile, 'name, stdev%s\n', quantityName);
        for ii=1:1:length(perGxySds)
            fprintf(outFile, '%s, %2.4f\n', sdSortedGxyNames{ii}, perGxySds(ii));
        end
        fclose(outFile);
    catch MErr
        warning('doPerturbationAnalysis:failedToWritePerGxyStabilityFile', ...
            'could not write a per-galaxy output-stability file (%s)', MErr.message);
        fclose(outFile)
    end   
end