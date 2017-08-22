function [lengths, ranks] = getGenomeScores(arcLengths, totalLenRatios)
% Calculates information needed for measures of arc fragmentation
% INPUTS:
%   arcLengths: lengths of all measured log-spiral arcs for the image, in
%       any order
%   totalLenRatio: where, along the total length of all arcs, the measures
%       should be taken. These are expressed as a proportion of the total
%       length, so all values in this array should be in [0, 1].
% OUTPUTS:
%   lengths: lengths of the arcs at the given proportions of the total
%       length
%   ranks: arc-length ranks of the arcs at the given proportion of the
%       total length

if any(totalLenRatios < 0) || any(totalLenRatios > 1)
    error('totalLenRatios values must be in [0, 1]');
end

isBadArc = (arcLengths == 0) | isnan(arcLengths);
arcLengths = arcLengths(~isBadArc);
if numel(arcLengths) == 0
    lengths = -inf * ones(size(totalLenRatios));
    ranks = -inf * ones(size(totalLenRatios));
    return;
end

arcLengths = sort(arcLengths, 'descend');
totalLength = sum(arcLengths);
meanLength = mean(arcLengths);
alenCSum = cumsum(arcLengths);
tgtIdxs = arrayfun(@(x)(find(alenCSum >= (totalLength * x), 1, 'first')), totalLenRatios);
ranks = tgtIdxs;
tgtLenVals = arcLengths(tgtIdxs);
lengths = tgtLenVals;
% scs = tgtLenVals ./ meanLength;

% midLen = arcLengths(midLenIdx)
% midLen / meanLength

end