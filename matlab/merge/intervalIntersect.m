function ivals = intervalIntersect(ivals1, ivals2)
% intersection of two intervals, each interval considered inclusive

if size(ivals1, 2) ~= 2 || size(ivals2, 2) ~= 2
    error('each interval set must have two columns');
end

if size(ivals1, 1) ~= size(ivals2, 1)
    error('interval sets must have the same number of rows');
end

nIvals = size(ivals1, 1);

ivals = zeros(nIvals, 2);
ivals(:, 1) = max([min(ivals1, [], 2), min(ivals2, [], 2)], [], 2);
ivals(:, 2) = min([max(ivals1, [], 2), max(ivals2, [], 2)], [], 2);

% ival = [0 0];
% ival(1) = max(min(ival1), min(ival2));
% ival(2) = min(max(ival1), max(ival2));

% if ival(1) > ival(2)
%     ival = [];
% end

emptyIval = ivals(:, 1) > ivals(:, 2);
ivals(emptyIval, 1) = NaN;
ivals(emptyIval, 2) = NaN;


end