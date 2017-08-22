function locMax = nonmaxSuppress(mtx, thresMult, nhSize, wrapCols)
% Removes points from the given array that aren't local maxima
% INPUTS:
%   mtx: the matrix for which nonmax-suppression is to be performed
%   thresMult: multiple of the mean to use for thresholding local maxima;
%     no points with values below this threshold will be considered local
%     maxima
%   nhSize: neighboorhood radius to consider when looking for local maxima;
%     in order for a point (r, c) to be considered a local maximum, it
%     cannot be smaller than any point within [r(+/-)nhSize, c(+/-)nhSize],
%     inclusive.
%   wrapCols: whether to consider the first and last columns to be adjacent
%     to each other (default false)
% OUTPUTS:
%   locMax: matrix with the same dimensions as the input matrix, with
%     values equal to their original value if they are local maxima within
%     the specified neighborhood, and zero otherwise.

if length(size(mtx)) ~= 2 
    error('mtx must be 2D')
end

if nargin < 4 || isempty(wrapCols)
    wrapCols = false;
end

thres = thresMult * mean(mtx(:));

[mRows, mCols] = size(mtx);
maxs = zeros(mRows, mCols);
if wrapCols
    % consider wraparound columns to be neighbors
    mtxp = padarray(mtx, [0 nhSize], 'circular');
    % rows do not wrap around
    mtxp = padarray(mtxp, [nhSize 0], -inf);
else
    mtxp = padarray(mtx, nhSize * [1 1], -inf);
end
figure; imagesc(mtxp); axis image

for r = nhSize+1 : 1 : mRows + nhSize
    for c = nhSize+1 : 1 : mCols + nhSize
        nbrs = mtxp(r-nhSize:r+nhSize, c-nhSize:c+nhSize);
        % an element doesn't need to be bigger than itself:
        nbrs(1+nhSize,1+nhSize) = -Inf; 
        maxs(r - nhSize, c - nhSize) = max(max(nbrs));
    end
end

locMax = mtx .* (mtx >= maxs) .* (maxs > thres);

end