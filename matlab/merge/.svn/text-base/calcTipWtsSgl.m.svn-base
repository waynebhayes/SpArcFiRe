function wts = calcTipWtsSgl(brt, theta, lgspParams, lgspBounds, dctn, imgSz, fullWtThSpan)

error(nargchk(6, 7, nargin));

minTipWtSigma = imgSz(1) / 25;
maxTipWtSigma = imgSz(1) / 4;

if nargin < 7 || isempty(fullWtThSpan)
    fullWtThSpan = 0;
end

assert(~isnan(fullWtThSpan))

if ~isscalar(fullWtThSpan)
    error('fullWtThSpan should be a scalar');
end

if isempty(dctn) || dctn == 0
    wts = brt;
    return;
elseif dctn < 0
%     tipEnd = min(theta);
    tipEnd = min(lgspBounds);
    fullWtAlen = calcLgspArcLengths(lgspParams, [tipEnd tipEnd + fullWtThSpan]);
else
%     tipEnd = max(theta);
    tipEnd = max(lgspBounds);
    fullWtAlen = calcLgspArcLengths(lgspParams, [tipEnd tipEnd - fullWtThSpan]);
end
assert(isscalar(fullWtAlen));

thetaForAlen = theta - min(theta) + min(lgspBounds);
ptAlens = calcLgspArcLengths(...
    repmat(lgspParams, length(theta), 1),...
    [repmat(tipEnd, length(theta), 1), thetaForAlen]);
[alenSrt, alenSrtIdxs] = sort(ptAlens);
[temp, restoreIdxs] = sort(alenSrtIdxs);
[maxGap, maxGapLoc] = max(diff(alenSrt));
closedGap = false;
while maxGap > 1
    % don't count arc-length that crosses gaps in the cluster (so that we
    %  fit to more than a fragment at the end of the cluster)
    closedGap = true;
    alenSrt(maxGapLoc+1:end) = alenSrt(maxGapLoc+1:end) - maxGap + eps;
    [maxGap, maxGapLoc] = max(diff(alenSrt));
end
if closedGap
    ptAlens = alenSrt(restoreIdxs);
end

ptAlens = max(0, ptAlens - fullWtAlen);
tipWtSigma = max(ptAlens)/4;
tipWtSigma = max(tipWtSigma, minTipWtSigma);
tipWtSigma = min(tipWtSigma, maxTipWtSigma);

wts = normpdf(ptAlens, 0, tipWtSigma);
wts = wts / max(wts(:));
wts = brt .* wts + eps;
% add eps so all points still detected as cluster members by find()

assert(sum(isnan(wts)) == 0);

end