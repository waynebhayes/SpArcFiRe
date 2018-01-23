function [bnds2Trans, rotAdj, minDist, th2Trans] = ...
    translateThVals(lgspParams1, bnds1, lgspParams2, bnds2, th2)
% Determines bounds (and optionally theta-values) of a second log-spiral 
% arc in terms of the parameters (theta-offset) of the first
% NOTE: we can only use these bounds with the parameters for the first
% spiral; shifting the theta values of the second arc invalidates the
% parameters of the second arc (changing the theta-offset doesn't just
% rotate the arc)

error(nargchk(4, 5, nargin));

if nargout >= 3 && nargin < 5
    error('must specify th2 if th2Trans is desired');
end

thOff1 = lgspParams1(1);
thOff2 = lgspParams2(1);

% bounds of the second arc in terms of the offset for the first, up to a
% multiple of 2*pi.
bnds2Trans = bnds2 - thOff1 + thOff2;

% deal with non-uniqueness of theta-values (find the correct multiple 
% of 2pi, and add it to the second cluster's bounds). 
rotAdjs = [-2*pi, 0, 2*pi];
thDists = arrayfun(@(x)(min(abs([bnds1(1) - (bnds2Trans + x), bnds1(2) - (bnds2Trans + x)]))), rotAdjs);
[minDist, minIdx] = min(thDists);
rotAdj = rotAdjs(minIdx);
bnds2Trans = bnds2Trans + rotAdj;

if nargout >= 3
    % align the theta-values with the bounds, then adjust for theta-offset
    th2Trans =  th2 + (min(bnds2Trans) - min(th2)) + thOff1;
end
    
end