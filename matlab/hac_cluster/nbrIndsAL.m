function [nbrInds, posns] = nbrIndsAL(ind, nRows)
% Determines the indexes of points above and to the left of the given index
%   (including diagonals)

nbrInds = ind + [-nRows + [-1 0 1], -1];
inRange = true(1, 4) & ((ind > nRows) | [false(1, 3) true]) & ...
    ((mod(ind, nRows) ~= 1) | [false true true false]) & ...
    ((mod(ind, nRows) ~= 0) | [true true false true]);
nbrInds = nbrInds(inRange);

if nargout > 1
    posns = find(inRange);
end

end