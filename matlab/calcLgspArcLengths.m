function lengths = calcLgspArcLengths(lgspParams, lgspBounds)
% Calculates arc-lengths for the given set of log-spirals
% INPUTS:
%   lgspParams: parameters of the log-spirals, one per row
%   lgspBounds: starting and ending theta-range of points to which
%       log-spirals were fitted, counterclockwise from the theta-offset
%       parameter.  Each row is a [start, end] pair.
% OUTPUTS:
%   lengths: log-spiral arc lengths

% TODO: handle case where a is exactly zero

a = lgspParams(:, 2); ir = lgspParams(:, 3);
thStart = lgspBounds(:, 1); thEnd = lgspBounds(:, 2);

rStart = ir .* exp(-a .* thStart);
rEnd = ir .* exp(-a .* thEnd);

lengths = (sqrt(1 + a.^2) ./ a) .* (rStart - rEnd);
aIsZero = (a == 0);
lengths(aIsZero) = ir(aIsZero) .* (thEnd(aIsZero) - thStart(aIsZero));
lengths = abs(lengths);

end