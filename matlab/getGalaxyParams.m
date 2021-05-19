function gxyParams = ...
    getGalaxyParams(lgspParams, lgspBounds, lgspErrs, used2rev, failed2rev, hasBadBounds, ...
        barUsed, barInfo, barIdxs, img, clusMtxs, stgs, gxyParams)
% Calculates some estimates of galactic structure (e.g., chirality,
%   arm-count, pitch angle) from log-spiral and bar parameters, and
%   consolidates this information (including the raw log-spiral parameters)
%   into a struct
% INPUTS: 
%   lgspParams: parameters of the log-spirals, one per row
%   lgspBounds: starting and ending theta-range of points to which
%       log-spirals were fitted, counterclockwise from the theta-offset
%       parameter.  Each row is a [start, end] pair.
%   lgspErrs: 
%   barInfo: 
%   barIdxs: 
%   img: the (preprocessed) image to which the log-spiral arcs were fit
%   gxyParams: struct with any fields that should be included in the
%       output struct
% OUTPUT:
%   gxyParams: the struct with estimated galaxy parameter information

error(nargchk(13, 13, nargin));

len_chk_thres = min(stgs.resizeDims) / 2;
err_per_len_warn_thres = 500;

% gxyParams = struct('arc_params', {}, 'chirality_votes_maj', {}, ...
%     'chirality_votes_alenWtd', {}, 'pitch_angle_avg', {}, ...
%     'pitch_angle_alenWtd', {}, 'thresholded_arm_counts', {});
% gxyParams = struct([]);

% fprintf(outfile, 'number of arcs: %d\n\n', size(lgspParams, 1));

isBar = false(size(clusMtxs, 3), 1); isBar(barIdxs) = true;
lgspClusMtxs = clusMtxs(:, :, ~isBar);
lgspParams = lgspParams(~isBar, :);
lgspBounds = lgspBounds(~isBar, :);
lgspErrs = lgspErrs(~isBar);
used2rev = used2rev(~isBar);
hasBadBounds = hasBadBounds(~isBar);
assert(length(used2rev) == length(failed2rev));
assert(length(failed2rev) == length(hasBadBounds));
assert(length(hasBadBounds) == length(lgspErrs));

imgSz = size(img);
minImgDim = min(imgSz);

gxyParams.failed2revInOutput = any(failed2rev);

iptCtrXY = gxyParams.iptCtrXY;
iptImgSz = gxyParams.iptSz;
gxyParams.inputCenterR = iptImgSz(1) - iptCtrXY(2) + 1;
gxyParams.inputCenterC = iptCtrXY(1);

% 5/19 - Sahel's fix for csv's not being generated when
% image standardization is turned off. This in conjunction
% with the fixes to reProject now allow us to output csv's 
% with information in them... mostly arm level. 
%
% There are still bugs to hunt in the output of the regular galaxy csv's.

if isfield(gxyParams.fitParams, 'covarFit')
    gxyParams.covarFit = gxyParams.fitParams.covarFit;
else
    gxyParams.covarFit = 0
end

arcLengths = calcLgspArcLengths(lgspParams, lgspBounds);

lgspParams(:, 2) = lgspParams(:, 2) * 180/pi;
pixCounts = squeeze(sum(sum(lgspClusMtxs > 0, 1), 2));
avgBrt = squeeze(sum(sum(lgspClusMtxs, 1), 2)) ./ pixCounts;
% calcClusScores(unsharpMask(img), lgspClusMtxs)
[clusBrtRatioScores, clusAnovaFScores, gxyParams] = calcClusScores(img, lgspClusMtxs, gxyParams);
% clusScoresUsm = calcClusScores(unsharpMask(img, stgs), lgspClusMtxs);

[arcLengths, alenOrder] = sort(arcLengths, 'descend');
lgspErrs = lgspErrs(alenOrder);
pixCounts = pixCounts(alenOrder);
avgBrt = avgBrt(alenOrder);
lgspParams = lgspParams(alenOrder, :);
lgspBounds = lgspBounds(alenOrder, :);
hasBadBounds = hasBadBounds(alenOrder);
clusBrtRatioScores = clusBrtRatioScores(alenOrder);
clusAnovaFScores = clusAnovaFScores(alenOrder);
used2rev = used2rev(alenOrder);
failed2rev = failed2rev(alenOrder);

lgspRRanges = -ones(size(lgspBounds));
for ii=1:1:size(lgspRRanges, 1)
    lgspRRanges(ii, :) = ...
        logSpiralFxn2Rev(lgspBounds(ii, :)' + lgspParams(ii, 1), ...
        lgspParams(ii, :) .* [1 (pi/180) 1]);
end
reordRow = abs(lgspRRanges(:, 1)) > abs(lgspRRanges(:, 2));
lgspRRanges(reordRow, :) = lgspRRanges(reordRow, [2 1]);

% if this is reordered, update pitch angle index in writeGxyParamsToCsv
gxyParams.arc_params = [lgspParams, lgspBounds, lgspRRanges, ...
    used2rev, failed2rev, hasBadBounds, arcLengths, pixCounts, ...
    lgspErrs ./ arcLengths, lgspErrs ./ pixCounts, avgBrt,...
    clusBrtRatioScores, clusAnovaFScores];

above_chk_thres = (arcLengths > len_chk_thres);
maxErrPerLen = max(lgspErrs(above_chk_thres) ./ arcLengths(above_chk_thres));
if maxErrPerLen > err_per_len_warn_thres 
    newWarn = sprintf('lgspFitting:highErrPerLenForLongArc:%2.2f', maxErrPerLen);
    gxyParams.warnings = [gxyParams.warnings newWarn];
end

gxyParams.bar_candidate_available = barInfo.barDetected;
gxyParams.bar_used = barUsed;
gxyParams.bar_info = barInfo;

gxyParams.avg_brt_ratio_score = mean(clusBrtRatioScores);
gxyParams.avg_anovaF_score = mean(clusAnovaFScores);
gxyParams.median_brt_ratio_score = median(clusBrtRatioScores);
gxyParams.median_anovaF_score = median(clusAnovaFScores);

[lengths, ranks] = getGenomeScores(arcLengths, [0.25 0.50 0.75]);
gxyParams.alenAt25pct = lengths(1);
gxyParams.alenAt50pct = lengths(2);
gxyParams.alenAt75pct = lengths(3);
gxyParams.rankAt25pct = ranks(1);
gxyParams.rankAt50pct = ranks(2);
gxyParams.rankAt75pct = ranks(3);
gxyParams.avgArcLength = mean(arcLengths);
gxyParams.minArcLength = min(arcLengths);
gxyParams.lowerQuartileArcLength = quantile(arcLengths, 0.25);
gxyParams.medianArcLength = quantile(arcLengths, 0.5);
gxyParams.upperQuartileArcLength = quantile(arcLengths, 0.75);
gxyParams.maxArcLength = max(arcLengths);
gxyParams.totalArcLength = sum(arcLengths);
gxyParams.totalNumArcs = numel(arcLengths);

pitchAngles = lgspParams(:, 2);
pitchAngles = (180/pi) * atan(pitchAngles * (pi/180));

% chirality by majority vote
cwVotes = sum(pitchAngles > 0);
ccwVotes = sum(pitchAngles < 0);
if cwVotes > ccwVotes
    dctnVote = 1;
elseif ccwVotes > cwVotes
    dctnVote = -1;
else
    dctnVote = 0;
end
gxyParams.chirality_votes_maj = [cwVotes ccwVotes];

% chirality by arc-length-weighted vote
if sum(arcLengths < 0) > 0
    error('internal error: negative arc lengths!');
end
cwAlenVotes = sum((pitchAngles > 0) .* arcLengths);
ccwAlenVotes = sum((pitchAngles < 0) .* arcLengths);
if cwAlenVotes > ccwAlenVotes
    alenDctnVote = 1;
elseif ccwAlenVotes > cwAlenVotes
    alenDctnVote = -1;
else
    alenDctnVote = 0;
end
gxyParams.chirality_votes_alenWtd = [cwAlenVotes ccwAlenVotes];

% weighted sum of pitch angles, as another way to determine chirality
gxyParams.alenWtdPangSum = sum(pitchAngles .* arcLengths);

if ~isempty(pitchAngles)
    gxyParams.chirality_longest_arc = sign(pitchAngles(1));
else
    gxyParams.chirality_longest_arc = 0;
end

% are the top 2 longest arms fairly long, and disagree in pitch angle?
if length(arcLengths) < 2
    gxyParams.top2_chirality_agreement = '<2 arcs';
elseif arcLengths(1) < (minImgDim / 4)
    gxyParams.top2_chirality_agreement = 'all-short';
elseif arcLengths(2) < (minImgDim / 4)
    gxyParams.top2_chirality_agreement = 'one-long';
elseif sign(pitchAngles(1)) == sign(pitchAngles(2))
    gxyParams.top2_chirality_agreement = 'agree';
else
    gxyParams.top2_chirality_agreement = 'disagree';
end

if ~isempty(pitchAngles)
    gxyParams.pa_longest = pitchAngles(1);
else
    gxyParams.pa_longest = [];
end

% average pitch angle
% majVpitchAngles = pitchAngles(sign(pitchAngles) == dctnVote);
isDomChirality = (sign(pitchAngles) == alenDctnVote);
alenVpitchAngles = pitchAngles(isDomChirality);
% meanPangAll = mean(abs(pitchAngles));
% meanPangMajDctn = mean(abs(majVpitchAngles));
% meanPangAlenWtdDctn = mean(abs(alenVpitchAngles));
% gxyParams.pitch_angle_avg = ...
%     [meanPangAll meanPangMajDctn meanPangAlenWtdDctn];
gxyParams.pa_avg = mean(pitchAngles);
gxyParams.pa_avg_abs = mean(abs(pitchAngles));
gxyParams.pa_avg_domChiralityOnly = mean(alenVpitchAngles);

% arc-length-weighted pitch angle
% majVarcLengths = arcLengths(sign(pitchAngles) == dctnVote);
alenVarcLengths = arcLengths(isDomChirality);
% wMeanPangAll = sum(abs(pitchAngles) .* arcLengths) / sum(arcLengths);
% wMeanPangMajDctn = ...
%     sum(abs(majVpitchAngles) .* majVarcLengths) / sum(majVarcLengths);
% wMeanPangAlenWtdDctn = ...
%     sum(abs(alenVpitchAngles) .* alenVarcLengths) / sum(alenVarcLengths);
% gxyParams.pitch_angle_alenWtd = ...
%     [wMeanPangAll, wMeanPangMajDctn, wMeanPangAlenWtdDctn];
gxyParams.pa_alenWtd_avg = sum(pitchAngles .* arcLengths) / sum(arcLengths);
gxyParams.paErr_alenWtd_stdev = sqrt(...
    sum(arcLengths .* (pitchAngles - gxyParams.pa_alenWtd_avg).^2) / ...
    sum(arcLengths));
gxyParams.pa_alenWtd_avg_abs = sum(abs(pitchAngles) .* arcLengths) / sum(arcLengths);
gxyParams.paErr_alenWtd_stdev_abs = sqrt(...
    sum(arcLengths .* (abs(pitchAngles) - gxyParams.pa_alenWtd_avg_abs).^2) / ...
    sum(arcLengths));
gxyParams.pa_alenWtd_avg_domChiralityOnly = sum(alenVpitchAngles .* alenVarcLengths) / sum(alenVarcLengths);
gxyParams.paErr_alenWtd_stdev_domChiralityOnly = sqrt(...
    sum(alenVarcLengths .* (alenVpitchAngles - gxyParams.pa_alenWtd_avg_domChiralityOnly).^2) / ...
    sum(alenVarcLengths));

totBrtDomChirality = avgBrt .* pixCounts;
totBrtDomChirality = totBrtDomChirality(isDomChirality);
gxyParams.pa_totBrtWtd = sum(alenVpitchAngles .* totBrtDomChirality) / sum(totBrtDomChirality);
gxyParams.pa_avgBrtWtd = sum(alenVpitchAngles .* avgBrt(isDomChirality)) / sum(avgBrt(isDomChirality));

paqs = getPitchAngleQuantiles(...
    lgspParams .* repmat([1 (pi/180) 1], size(lgspParams, 1), 1), ...
    lgspBounds, ...
    [0.1 0.25 0.5 0.75 0.9]) * (180/pi);
if paqs(3) < 0
    % lower quartile/decile should still have smaller pitch angle
    paqs = paqs(end:-1:1);
end
paqs = (180/pi) * atan(paqs * (pi/180));
gxyParams.pa_alenWtd_median = paqs(3);
gxyParams.pa_alenWtd_lowQuartile = paqs(2);
gxyParams.pa_alenWtd_highQuartile = paqs(4);
gxyParams.pa_alenWtd_lowDecile = paqs(1);
gxyParams.pa_alenWtd_highDecile = paqs(5);

paqsDco = getPitchAngleQuantiles(...
    lgspParams(isDomChirality, :) .* repmat([1 (pi/180) 1], sum(isDomChirality), 1), ...
    lgspBounds(isDomChirality, :), ...
    [0.1 0.25 0.5 0.75 0.9]) * (180/pi);
if paqsDco(3) < 0
    % lower quartile/decile should still have smaller pitch angle
    paqsDco = paqsDco(end:-1:1);
end
paqsDco = (180/pi) * atan(paqsDco * (pi/180));
gxyParams.pa_alenWtd_median_domChiralityOnly = paqsDco(3);
gxyParams.pa_alenWtd_lowQuartile_domChiralityOnly = paqsDco(2);
gxyParams.pa_alenWtd_highQuartile_domChiralityOnly = paqsDco(4);
gxyParams.pa_alenWtd_lowDecile_domChiralityOnly = paqsDco(1);
gxyParams.pa_alenWtd_highDecile_domChiralityOnly = paqsDco(5);

% [alenVpitchAngles alenVarcLengths]
% we added a sort above, do we still need to sort here?
lenSort = sortrows([alenVpitchAngles alenVarcLengths], -2);
gxyParams.sorted_agreeing_pangs = lenSort(:, 1)';
gxyParams.sorted_agreeing_arclengths = lenSort(:, 2)';

% try to determine arm count by differences in arc lengths
% largest gap
gxyParams.numArcs_largest_length_gap = getAlenGapInd(arcLengths);
gxyParams.numDcoArcs_largest_length_gap = getAlenGapInd(alenVarcLengths);
% flattening out
gxyParams.numArcs_arclength_function_flattening = getAlenFlattenInd(arcLengths);
gxyParams.numDcoArcs_arclength_function_flattening = getAlenFlattenInd(alenVarcLengths);

% thresholds = minImgDim * [1.2:-0.05:0]';
% thresholds = [300:-15:0]';
thresholds = unique([0:20:300 10 50:5:100 350:50:600])';
gxyParams.length_thresholds = thresholds;
gxyParams.thres_arm_counts = ...
    getCumulativeCounts(thresholds, arcLengths);
gxyParams.thres_arm_counts_dco = ...
    getCumulativeCounts(thresholds, alenVarcLengths);
gxyParams.thres_pang_avgs = ...
    getCumulativeAverages(thresholds, arcLengths, pitchAngles);
gxyParams.thres_abs_pang_avgs = ...
    getCumulativeAverages(thresholds, arcLengths, abs(pitchAngles));
gxyParams.thres_pang_avgs_domChirality = ...
    getCumulativeAverages(thresholds, alenVarcLengths, alenVpitchAngles);

[clusImg, clusColors] = showClustersFromMtxs(clusMtxs, imgSz);
clusColors = clusColors(~isBar, :);
gxyParams.arcClusterColors = clusColors(alenOrder, :);

% gxyParams

    function counts = getCumulativeCounts(thresholds, measurements)
        counts = zeros(length(thresholds), 1);
        for thres=1:1:length(thresholds)
            counts(thres) = sum(measurements >= thresholds(thres));
        end
    end

    function avgs = getCumulativeAverages(thresholds, measurements, values)
        avgs = zeros(length(thresholds), 1);
        for thres=1:1:length(thresholds)
            avgs(thres) = mean(values(measurements >= thresholds(thres)));
            if isnan(avgs(thres))
                avgs(thres) = 0;
            end
        end
    end

    function mgap_ind = getAlenGapInd(arcLengths)
        [mgap_val, mgap_ind] = max(diff(arcLengths(end:-1:1)));
        mgap_ind = length(arcLengths) - mgap_ind; % arcLengths list was reversed
        if isempty(mgap_ind)
            mgap_ind = 0;
        end
    end

    function flatten_ind = getAlenFlattenInd(arcLengths)
        d2alens = diff(diff(arcLengths));
        if ~isempty(d2alens)
            if d2alens(1) < 0
                d2alens(1:find(d2alens >= 0, 1, 'first')-1) = 0;
            end
            flatten_ind = find(d2alens < 0, 1, 'first') - 1;
            if isempty(flatten_ind)
                flatten_ind = 0;
            end
        else
            flatten_ind = 0;
        end
    end
end