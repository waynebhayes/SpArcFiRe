function writeGalaxyParams(gxyParams, gxyName, outputDir)
% Writes the information in the gxyParams structure to a human-readable
%   text file

if nargin < 2
    gxyName = [];
end

if nargin < 3 || isempty(outputDir)
    outputDir = ['.' filesep];
elseif ~isdir(outputDir)
    error('%s is not a directory', outputDir);
end

if outputDir(end) ~= filesep
    outputDir(end+1) = filesep;
end

if ~isempty(gxyName)
    outfile = fopen([outputDir gxyName '-params.txt'], 'wt');
else
    outfile = 1;
end

cwStr = 'S-wise';
ccwStr = 'Z-wise';
eqStr = 'Equal';

fprintf(outfile, 'galaxy name: %s\n\n', gxyName);

fprintf(outfile, 'Gaussian-fit parameters: \n');
fprintf(outfile, 'axis ratio: %2.4f\n', gxyParams.diskAxisRatio);
fprintf(outfile, 'minor axis length: %2.4f\n', gxyParams.diskMinAxsLen);
fprintf(outfile, 'major axis length: %2.4f\n', gxyParams.diskMajAxsLen);
fprintf(outfile, 'major axis angle: %2.4f\n', gxyParams.diskMajAxsAngleRadians);
fprintf(outfile, '\n');

fprintf(outfile, 'mean brightness-ratio score: %2.4f\n', gxyParams.avg_brt_ratio_score);
fprintf(outfile, 'mean anovaF score: %2.4f\n', gxyParams.avg_anovaF_score);
fprintf(outfile, 'median brightness-ratio score: %2.4f\n', gxyParams.median_brt_ratio_score);
fprintf(outfile, 'median anovaF score: %2.4f\n', gxyParams.median_anovaF_score);
fprintf(outfile, '\n');

fprintf(outfile, 'arc lengths at 25%%, 50%%, and 75%% of total length: %s\n',...
    mat2str([gxyParams.alenAt25pct gxyParams.alenAt50pct gxyParams.alenAt75pct], 4));
fprintf(outfile, 'arc-length ranks at 25%%, 50%%, and 75%% of total length: %s\n',...
    mat2str([gxyParams.rankAt25pct gxyParams.rankAt50pct gxyParams.rankAt75pct], 4));
fprintf(outfile, 'mean arc length: %2.4f\n', gxyParams.avgArcLength);
fprintf(outfile, 'total arc length: %2.4f\n', gxyParams.totalArcLength);
fprintf(outfile, '\n');

fprintf(outfile, 'log-spiral arc parameters:\n');
fprintf(outfile, ['theta offset\tpitch angle\t\tinitial radius'...
    '\t\ttheta range\t\tradial range\t\tarc length\t\t# pixels\terr/len\t\terr/pixels\t\tmean intensity\t\tbrtRatioScore\t\tanovaFScore\n']);
fprintf(outfile, '%2.4f\t\t\t%8.4f\t\t%9.4f\t\t[%2.4f, %2.4f]\t[%2.4f,%2.4f]\t%8.4f\t\t%4d\t\t%8.4f\t%8.4f\t\t%6.4f\t\t\t\t%2.4f\t\t\t\t%2.4f\n',...
    gxyParams.arc_params');
fprintf(outfile, 'end of arc parameters\n\n');

barInfo = gxyParams.bar_info;
fprintf(outfile, 'bar score: %2.4f ', barInfo.score);
if barInfo.barDetected
    fprintf(outfile, '(bar detected)\n');
    fprintf(outfile, '\tcenter (original image): R = %2.4f, C = %2.4f\n',...
        barInfo.iptCtrR, barInfo.iptCtrC);
    fprintf(outfile, '\tangle (original image): %2.4f\n', ...
        barInfo.iptAngle * (180/pi) );
    fprintf(outfile, '\thalf-length (original image): %2.4f\n', ...
        barInfo.iptHalfLength);
    fprintf(outfile, '\tangle (standardized image): %2.4f\n', ...
        barInfo.stdzAngle * (180/pi) );
    fprintf(outfile, '\thalf-length (standardized image): %2.4f\n', ...
        barInfo.stdzHalfLength);
    fprintf(outfile, '\n');
else
    fprintf(outfile, '(no bar detected)\n\n');
end

cwVotes = gxyParams.chirality_votes_maj(1);
ccwVotes = gxyParams.chirality_votes_maj(2);
if cwVotes > ccwVotes
    dctnVoteStr = cwStr;
elseif ccwVotes > cwVotes
    dctnVoteStr = ccwStr;
else
    dctnVoteStr = eqStr;
end
numDctnVotes = cwVotes + ccwVotes;
fprintf(outfile, 'chirality (majority vote): %s\n', dctnVoteStr);
fprintf(outfile, 'S-wise votes: %d of %d (%2.4f%%)\n',...
    cwVotes, numDctnVotes, cwVotes/numDctnVotes * 100);
fprintf(outfile, 'Z-wise votes: %d of %d (%2.4f%%)\n',...
    ccwVotes, numDctnVotes, ccwVotes/numDctnVotes * 100);
fprintf(outfile, '\n');

cwAlenVotes = gxyParams.chirality_votes_alenWtd(1);
ccwAlenVotes = gxyParams.chirality_votes_alenWtd(2);
if cwAlenVotes > ccwAlenVotes
    alenDctnVoteStr = cwStr;
elseif ccwAlenVotes > cwAlenVotes
    alenDctnVoteStr = ccwStr;
else
    alenDctnVoteStr = eqStr;
end
numAlenDctnVotes = cwAlenVotes + ccwAlenVotes;
fprintf(outfile, 'chirality (arc-length-weighted vote): %s\n', alenDctnVoteStr);
fprintf(outfile, 'S-wise votes: %2.4f of %2.4f (%2.4f%%)\n',...
    cwAlenVotes, numAlenDctnVotes, cwAlenVotes/numAlenDctnVotes * 100);
fprintf(outfile, 'Z-wise votes: %2.4f of %2.4f (%2.4f%%)\n',...
    ccwAlenVotes, numAlenDctnVotes, ccwAlenVotes/numAlenDctnVotes * 100);
fprintf(outfile, '\n');

if gxyParams.alenWtdPangSum > 0
    alenWtdSumDctnStr = cwStr;
elseif gxyParams.alenWtdPangSum < 0
    alenWtdSumDctnStr = ccwStr;
else
    alenWtdSumDctnStr = eqStr;
end
fprintf(outfile, 'chirality (weighted-pitch-angle-sum): %s (%2.4f)\n\n',...
    alenWtdSumDctnStr, gxyParams.alenWtdPangSum);

fprintf(outfile, 'chirality agreement of 2 longest arcs?: %s\n\n', ...
    gxyParams.top2_chirality_agreement);

% fprintf(outfile, 'average pitch angle (all arcs): %2.4f\n', ...
%     gxyParams.pitch_angle_avg(1));
% fprintf(outfile, 'majority-vote chirality only: %2.4f\n', ...
%     gxyParams.pitch_angle_avg(2));
% fprintf(outfile, 'arc-length-weighted vote chirality only: %2.4f\n',...
%     gxyParams.pitch_angle_avg(3));
% fprintf(outfile, '\n');
% 
% fprintf(outfile, 'average pitch angle (arc-length-weighted): %2.4f\n',...
%     gxyParams.pitch_angle_alenWtd(1));
% fprintf(outfile, 'majority-vote chirality only: %2.4f\n',...
%     gxyParams.pitch_angle_alenWtd(2));
% fprintf(outfile, 'arc-length-weighted chirality only: %2.4f\n',...
%     gxyParams.pitch_angle_alenWtd(3));
% fprintf(outfile, '\n');

fprintf(outfile, 'pitch angle (longest arc): %2.4f\n', ...
    gxyParams.pa_longest);
fprintf(outfile, 'pitch angle (mean): %2.4f\n', ...
    gxyParams.pa_avg);
fprintf(outfile, 'pitch angle (mean-of-abs): %2.4f\n', ...
    gxyParams.pa_avg_abs);
fprintf(outfile, 'pitch angle (mean, detected chirality only): %2.4f\n',...
    gxyParams.pa_avg_domChiralityOnly);
fprintf(outfile, 'pitch angle (arc-length-weighted mean): %2.4f\n', ...
    gxyParams.pa_alenWtd_avg);
fprintf(outfile, 'pitch angle (arc-length-weighted mean-of-abs): %2.4f\n', ...
    gxyParams.pa_alenWtd_avg);
fprintf(outfile, 'pitch angle (weighted mean, detected chirality only): %2.4f\n',...
    gxyParams.pa_alenWtd_avg_domChiralityOnly);
fprintf(outfile, 'pitch angle (wtd by total brt, detected chirality only): %2.4f\n',...
    gxyParams.pa_totBrtWtd);
fprintf(outfile, 'pitch angle (wtd by avg brt, detected chirality only): %2.4f\n',...
    gxyParams.pa_avgBrtWtd);

fprintf(outfile, 'pitch angles, arc-length-weighted chirality only, sorted by length: \t%s\n\n',...
    mat2str(gxyParams.sorted_agreeing_pangs, 4));

fprintf(outfile, 'arm count by largest arc-length gap: %d\n', ...
    gxyParams.numArcs_largest_length_gap);
fprintf(outfile, 'arm count by arc-length function flattening: %d\n\n', ...
    gxyParams.numArcs_arclength_function_flattening);

% [gxyParams.length_thresholds, gxyParams.thres_arm_counts,...
%     gxyParams.thres_pang_avgs, gxyParams.thres_abs_pang_avgs]

fprintf(outfile, 'arm measures by length threshold: \n');
fprintf(outfile, 'threshold\tcount\t\tangle\t\tabs-angle\t\tdom-chirality-angle\n');
fprintf(outfile, '%6.2f\t\t%2d\t\t\t%8.4f\t\t\t%8.4f\t\t\t%8.4f\n', ...
    [gxyParams.length_thresholds, gxyParams.thres_arm_counts,...
    gxyParams.thres_pang_avgs, gxyParams.thres_abs_pang_avgs,...
    gxyParams.thres_pang_avgs_domChirality]');

fclose(outfile);

end