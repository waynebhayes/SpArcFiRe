function [wts1, wts2, isOverlap] = ...
    calcTipWtsDbl(brt1, theta1, lgspParams1, bnds1, brt2, theta2, lgspParams2, bnds2, imgSz)

error(nargchk(9, 9, nargin));

maxOverlap = pi/8;

bnds2Trans = translateThVals(lgspParams1, bnds1, lgspParams2, bnds2);

thOverlaps = intervalIntersect(repmat(bnds1, 3, 1) + [-2*pi,-2*pi;0,0;2*pi,2*pi], repmat(bnds2Trans, 3, 1));
[mV, mI] = max(abs(diff(thOverlaps, [], 2)));
thOverlap = thOverlaps(mI, :);
% if ( abs(diff(thOverlap)) - min(abs([diff(bnds1) diff(bnds2Trans)])) ) > -10e-6
if ~isnan(thOverlap(1)) && ...
        ( abs(diff(thOverlap)) > maxOverlap || ...
        (abs(diff(thOverlap)) > (min(abs([diff(bnds1) diff(bnds2Trans)])) / 2)) )  % TODO: confirm these criteria
%     fprintf('one interval contained largely within another.  Not performing tip-weighting.\n');
    wts1 = brt1;
    wts2 = brt2;
    isOverlap = true;
    return;
end

toMinTip = -1;
toMaxTip = 1;
if mean(bnds2Trans) < mean(bnds1)
    dctn1 = toMinTip;
    dctn2 = toMaxTip;
else
    dctn1 = toMaxTip;
    dctn2 = toMinTip;
end
if isempty(thOverlap) || (sum(isnan(thOverlap)) > 0)
    overlapAmt = 0;
else
    overlapAmt = abs(diff(thOverlap));
end
wts1 = calcTipWtsSgl(brt1, theta1, lgspParams1, bnds1, dctn1, imgSz, overlapAmt);
wts2 = calcTipWtsSgl(brt2, theta2, lgspParams2, bnds2, dctn2, imgSz, overlapAmt);

isOverlap = (overlapAmt > 0);

end