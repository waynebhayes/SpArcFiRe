function [nbrInds, posns] = nbrIndsBR(ind, nRows, nEl)
% Determines the indexes of points below and to the right of the given
%   index (including diagonals)

nbrInds = ind + [1, nRows + [-1 0 1]];
inRange = true(1, 4) & ((ind <= nEl - nRows) | [true false(1, 3)]) & ...
    ((mod(ind, nRows) ~= 1) | [true false true true]) & ...
    ((mod(ind, nRows) ~= 0) | [false true true false]);
nbrInds = nbrInds(inRange);

if nargout > 1
    posns = find(inRange);
end

end