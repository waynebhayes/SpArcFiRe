function paqVals = getPitchAngleQuantiles(lgspParams, lgspBounds, qLvls, nSamples)

error(nargchk(3, 4, nargin))

if any((qLvls < 0) | (qLvls > 1))
    error('qLvls must be in the range [0, 1]');
end

arcLengths = calcLgspArcLengths(lgspParams, lgspBounds);

isBadArc = (arcLengths == 0) | isnan(arcLengths);
arcLengths = arcLengths(~isBadArc);
lgspParams = lgspParams(~isBadArc, :);
lgspBounds = lgspBounds(~isBadArc, :);
if numel(arcLengths) == 0
    paqVals = NaN * ones(size(qLvls));
    return;
end
if (nargin < 4)
    nSamples = ceil(sum(arcLengths));
end

arcLengths = sort(arcLengths, 'descend');
totalLength = sum(arcLengths);
alenCSum = cumsum(arcLengths);
% alenCSum
% totalLength
% nSamples
% linspace(0, 1, nSamples)
% arrayfun(@(x)(find(alenCSum >= (totalLength * x), 1, 'first')), linspace(0, 1, nSamples), 'UniformOutput', false)
splIdxs = arrayfun(@(x)(find(alenCSum >= (totalLength * x), 1, 'first')), linspace(0, 1, nSamples));
paSamples = lgspParams(splIdxs, 2);
paqVals = quantile(paSamples, qLvls);

end