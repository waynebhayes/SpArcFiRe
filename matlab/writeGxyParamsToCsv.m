function writeGxyParamsToCsv(resPath, imgSetName, outDir, aggregate_matfiles)

if nargin < 3 || isempty(outDir)
    outDir = resPath;
end

if nargin < 4 || isempty(aggregate_matfiles)
    aggregate_matfiles = false;
end

if ~isdir(resPath)
    error('%s is not a directory\n', resPath);
end

if outDir(end) ~= filesep
    outDir(end+1) = filesep;
end
if resPath(end) ~= filesep
    resPath(end+1) = filesep;
end

% fileNames = ls(resPath);
% fileNames = strtrim(mat2cell(fileNames, ...
%     ones(1, size(fileNames, 1)), size(fileNames, 2)));
fileNames = dir(resPath);
fileNames = {fileNames.name};

fileNames = fileNames(~cellfun(@isempty, ...
    strfind(fileNames, 'gxyParamsList-part-')))

next_file_idx = 1;
sampleGxyParams = struct();
while(~isfield(sampleGxyParams, 'length_thresholds') && (next_file_idx <= length(fileNames)))
    curParamsList = load([resPath fileNames{next_file_idx}]);
    curParamsList = curParamsList.gxyParamsList;
    next_gxy_idx = 1;
    while(~isfield(sampleGxyParams, 'length_thresholds') && (next_gxy_idx <= length(curParamsList)))
        sampleGxyParams = curParamsList{next_gxy_idx};
        next_gxy_idx = next_gxy_idx + 1;
    end
    next_file_idx = next_file_idx + 1;
end
if next_file_idx == 1
    error('no gxyParamsList files found')
end
if ~isfield(sampleGxyParams, 'length_thresholds')
    error('None of the given galaxies have any detected arcs')
end

gxyParamsNames = {'badBulgeFitFlag', 'bulgeAxisRatio', 'bulgeMajAxsLen',...
    'bulgeMajAxsAngle', 'bulgeAvgBrt', 'bulgeDiskBrtRatio', 'numElpsRefits',...
    'diskAxisRatio', 'diskMinAxsLen', 'diskMajAxsLen', 'diskMajAxsAngleRadians',...
    'inputCenterR', 'inputCenterC', 'iptSz', 'muDist', 'muDistProp', 'covarFit', ...
    'cropRad', 'wtdLik', 'likOfCtr',... % Added cropRad - Matthew P. 1/7/21
    'brtUnifScore', 'gaussLogLik', 'contourBrtRatio', 'standardizedCenterR',...
    'standardizedCenterC', 'hasDeletedCtrClus', ...
    'failed2revDuringMergeCheck', 'failed2revDuringSecondaryMerging',...
    'failed2revInOutput', ...
    'fitQualityF1','fitQualityPCC',...
    'cpu_time', 'wall_time',...
    'bar_candidate_available', 'bar_used',...
    'alenAt25pct', 'alenAt50pct', 'alenAt75pct', 'rankAt25pct', ...
    'rankAt50pct', 'rankAt75pct', 'avgArcLength', 'minArcLength', ...
    'lowerQuartileArcLength', 'medianArcLength', 'upperQuartileArcLength',... 
    'maxArcLength', 'totalArcLength', 'totalNumArcs', 'chirality_votes_maj',...
    'chirality_votes_alenWtd', 'alenWtdPangSum', 'top2_chirality_agreement',...
    'pa_longest', 'pa_avg', 'pa_avg_abs', 'pa_avg_domChiralityOnly', ...
    'pa_alenWtd_avg', 'paErr_alenWtd_stdev', ...
    'pa_alenWtd_avg_abs', 'paErr_alenWtd_stdev_abs', ...
    'pa_alenWtd_avg_domChiralityOnly', 'paErr_alenWtd_stdev_domChiralityOnly',...
    'pa_totBrtWtd', 'pa_avgBrtWtd', 'pa_alenWtd_median', ...
    'pa_alenWtd_lowQuartile', 'pa_alenWtd_highQuartile', ...
    'pa_alenWtd_lowDecile', 'pa_alenWtd_highDecile', ...
    'pa_alenWtd_median_domChiralityOnly', 'pa_alenWtd_lowQuartile_domChiralityOnly',...
    'pa_alenWtd_highQuartile_domChiralityOnly', 'pa_alenWtd_lowDecile_domChiralityOnly',...
    'pa_alenWtd_highDecile_domChiralityOnly', 'sorted_agreeing_pangs', ...
    'sorted_agreeing_arclengths', ...
    'numArcs_largest_length_gap', 'numDcoArcs_largest_length_gap',...
    'numArcs_arclength_function_flattening', 'numDcoArcs_arclength_function_flattening'};
gxyParamsNamesExcluded = setdiff(fields(sampleGxyParams), gxyParamsNames);
for ii=1:length(gxyParamsNamesExcluded)
    fprintf('%s\n', gxyParamsNamesExcluded{ii});
end
fields(sampleGxyParams)
armCountThresholds = sampleGxyParams.length_thresholds;

csvFile = fopen([outDir imgSetName '.csv'], 'wt');

fprintf(csvFile, 'name, fit_state, warnings, star_mask_used, noise_mask_used, chirality_maj, chirality_alenWtd, chirality_wtdPangSum, chirality_longestArc, ');
fprintf(csvFile, '%s, ', gxyParamsNames{:});
fprintf(csvFile, ['bar_score_img, bar_cand_score_orifld, bar_angle_input_img, bar_half_length_input_img, '...
    'bar_angle_standardized_img, bar_half_length_standardized_img, ']);
for thresIdx=1:1:length(armCountThresholds)
%     nxt_header_item = strrep(...
%         sprintf('numArcsGE%.1f, ', armCountThresholds(thresIdx)), '.', 'pt');
%     fprintf(csvFile, nxt_header_item);
    fprintf(csvFile, 'numArcsGE%03d, ', round(armCountThresholds(thresIdx)));
end
for thresIdx=1:1:length(armCountThresholds)
%     nxt_header_item = strrep(...
%         sprintf('numDcoArcsGE%.1f, ', armCountThresholds(thresIdx)), '.', 'pt');
%     fprintf(csvFile, nxt_header_item);
    fprintf(csvFile, 'numDcoArcsGE%03d, ', round(armCountThresholds(thresIdx)));
end
fprintf(csvFile, '\n');

arcCsvFile = fopen([outDir imgSetName '_arcs.csv'], 'wt');

fprintf(arcCsvFile, ['gxyName,alenRank,math_initial_theta,pitch_angle,math_initial_radius,'...
    'relative_theta_start,relative_theta_end,r_start,r_end,used2rev,failed2rev,hasUndefinedEndpoints,'...
    'arc_length,num_pixels,err_per_len,err_per_pixel,'...
    'mean_intensity,brt_ratio_score,anovaF_score,clusOutputColor\n']);

% for ii=1:1:length(fileNames)
%     paramsListPart = load([resPath fileNames{ii}], 'gxyParamsList');
%     paramsListPart = paramsListPart.gxyParamsList;
%     paramsListPart = paramsListPart(cellfun(@(x)(~isempty(x)), paramsListPart));
%     gxyParamsList = [gxyParamsList; paramsListPart];
% end

if aggregate_matfiles
    gxyParamsList = {};
end

imgNum = 0;
for fIdx=1:1:length(fileNames)
    curGxyParamsList = load([resPath fileNames{fIdx}], 'gxyParamsList');
    curGxyParamsList = curGxyParamsList.gxyParamsList;
    curGxyParamsList = curGxyParamsList(cellfun(@(x)(~isempty(x)), curGxyParamsList));
    for imgIdx = 1:1:length(curGxyParamsList)
        imgNum = imgNum + 1;
        if mod(imgNum, 100) == 0
            fprintf('processing image %d\n', imgNum);
        end
        curParams = curGxyParamsList{imgIdx};
        writeCsvEntry(curParams, csvFile, gxyParamsNames);
        writeArcsCsvEntry(curParams, arcCsvFile)
    end
    if aggregate_matfiles
        gxyParamsList = [gxyParamsList; paramsListPart];
    end
end

fclose(csvFile);
fclose(arcCsvFile);

fprintf('done writing CSV files.\n');

if aggregate_matfiles
    save([resPath 'gxyParamsList-combined'], 'gxyParamsList');
end

end

function writeCsvEntry(curParams, csvFile, gxyParamsNames)
    fprintf(csvFile, [curParams.name ', ']);
    
    if strmatch('OK', curParams.fit_state)
        fprintf(csvFile, 'OK, ');
    else
%         fprintf(csvFile, [curParams.fit_state ', ' repmat(', ', 1, ...
%             (length(fields(curParams)) + length(armCountThresholds) - 3)) '\n']);
        fprintf(csvFile, [strrep(strrep(curParams.fit_state, sprintf('\n'), ' '), ',', ' ') '\n']);
        return;
    end
    
    curWarnings = curParams.warnings;
    if isempty(curWarnings)
        fprintf(csvFile, '(no warnings),');
    else
        for ii=1:1:length(curWarnings)
            fprintf(csvFile, [curWarnings{ii} '|']);
        end
        fprintf(csvFile, ',');
    end
    
    fprintf(csvFile, '%s ,', curParams.starMaskUsed);
    fprintf(csvFile, '%s ,', curParams.noiseMaskUsed);
    
    chVotes = curParams.chirality_votes_maj;
    if chVotes(1) > chVotes(2)
        fprintf(csvFile, 'S-wise, ');
    elseif chVotes(1) < chVotes(2)
        fprintf(csvFile, 'Z-wise, ');
    else
        fprintf(csvFile, 'EQ, ');
    end
    chVotes = curParams.chirality_votes_alenWtd;
    if chVotes(1) > chVotes(2)
        fprintf(csvFile, 'S-wise, ');
    elseif chVotes(1) < chVotes(2)
        fprintf(csvFile, 'Z-wise, ');
    else
        fprintf(csvFile, 'EQ, ');
    end
    if curParams.alenWtdPangSum > 0
        fprintf(csvFile, 'S-wise, ');
    elseif curParams.alenWtdPangSum < 0
        fprintf(csvFile, 'Z-wise, ');
    else
        fprintf(csvFile, 'EQ, ');
    end
    if curParams.chirality_longest_arc > 0
        fprintf(csvFile, 'S-wise, ');
    elseif curParams.chirality_longest_arc < 0
        fprintf(csvFile, 'Z-wise, ');
    else
        fprintf(csvFile, 'EQ, ');
    end
    for fldIdx = 1:1:length(gxyParamsNames)
        fprintf(csvFile, '%s, ', mat2str(getfield(curParams, gxyParamsNames{fldIdx}), 10));
    end
    barInfo = curParams.bar_info;
    fprintf(csvFile, '%2.4f, ', barInfo.score);
    fprintf(csvFile, '%2.4f, ', barInfo.candScore);
%     if barInfo.barDetected
        fprintf(csvFile, '%2.4f, %2.4f, %2.4f, %2.4f, ', ...
            barInfo.iptAngle * (180/pi), barInfo.iptHalfLength,...
            barInfo.stdzAngle * (180/pi), barInfo.stdzHalfLength);
%     else
%         fprintf(csvFile, 'noBarDetected, noBarDetected, ');
%     end
   
    armCounts = curParams.thres_arm_counts;
%     armCounts = armCounts(:, 2);
    fprintf(csvFile, '%d, ', armCounts);
    armCountsDco = curParams.thres_arm_counts_dco;
    fprintf(csvFile, '%d, ', armCountsDco);
    fprintf(csvFile, '\n');
end

function writeArcsCsvEntry(curGxyParams, arcCsvFile)
    curName = curGxyParams.name;
    if ~isfield(curGxyParams, 'arc_params') % image rejection or preprocessing failure
        return;
    end
    pitchAngleIndex = 2;
    curArcParams = curGxyParams.arc_params;
    curArcParams(:, pitchAngleIndex) = ...
        atan(curArcParams(:, pitchAngleIndex) * (pi/180)) * (180/pi);
    curClusColors = curGxyParams.arcClusterColors;
    for arcIdx=1:1:size(curArcParams, 1)
        fprintf(arcCsvFile, '%s,%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%s\n',...
            curName, arcIdx, curArcParams(arcIdx, :), mat2str(curClusColors(arcIdx, :), 4));
    end
end
