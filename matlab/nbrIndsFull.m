function [inds, posns] = nbrIndsFull(ind, nRows, nEl)
% Determines the indices of points adjacent to the given index (including
%   diagonals)

[indsAL, posnsAL] = nbrIndsAL(ind, nRows);
[indsBR, posnsBR] = nbrIndsBR(ind, nRows, nEl);
inds = [indsAL indsBR];
posns = [posnsAL (4 + posnsBR)];

end