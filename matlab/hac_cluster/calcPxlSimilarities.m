function [fromInds, toInds, simlVals, numElts] = ...
    calcPxlSimilarities(img, oriMtx, stgs, simlCutoff)
% Calculates pixel similarities using dot products of neighbors'
% orientation field vectors (similarities are zero if pixels are not within
% an 8-neighborhood of each other). The outputs can be used directly, or
% used to construct a similarity matrix.
% INPUTS:
%   img: image used to generate the orientation field
%   oriMtx: orientation field generated from the image
%   stgs: structure containing algorithm settings (see getDefaultSettings.m)
%   simlCutoff: minimum value for a similarity score to be nonzero in the
%       similarity matrix (this is often related to, but not necessarily
%       the same as, the HAC stopping threshold in the 'stgs' struct)
% OUTPUTS:
%   fromInds, toInds, simlVals: nonzero pixel similarities, such that the
%       similarity from pixel fromInds(ii) to pixel toInds(ii) is
%       simlVals(ii)
%   numElts: total number of similarities (generally a much larger value
%       than the length of the other outputs, since only neighbors have
%       nonzero similarities)

if nargin < 3 || isempty(simlCutoff)
    simlCutoff = -inf;
end

nhSize = stgs.nhSize;
nhSize = round(nhSize);

[nRows, nCols] = size(img);

numVecs = size(oriMtx, 1) * size(oriMtx, 2);
ofldStren = oriFieldStrengths(oriMtx);

ofldVec = reshape(oriMtx, numVecs, 2);
nzIdxs = find(reshape(ofldStren, numVecs, 1));
numNzVecs = length(nzIdxs);
fprintf('%d nonzero orientation field vectors\n', numNzVecs);
fprintf('%d nonzero image pixels\n', nnz(img));

maxNelSm = 9 * numNzVecs;
rIdxs = zeros(maxNelSm, 1);
cIdxs = zeros(maxNelSm, 1);
smVals = zeros(maxNelSm, 1);

nxtIdx = 1;
% simls = ofldVec' * speye(numVecs) * ofldVec;
for ii=1:1:numNzVecs
    curInd = nzIdxs(ii);
    
    [r, c] = ind2sub([nRows nCols], curInd);
	cn = c + nhSize * [-1 0 1 -1 0];
	rn = r + nhSize * [-1 -1 -1 0 0];
	inRange = (cn >= 1) & (cn <= nCols) & (rn >= 1);
    nbrInds = sub2ind([nRows nCols], rn(inRange), cn(inRange));
        
    nbrInds = nbrInds(:);
    nbrs = ofldVec(nbrInds, :);
    nbrSimls = nbrs * ofldVec(curInd, :)';
    
    nbrSimls(nbrSimls < simlCutoff) = 0;

%   simls(curInd, nbrInds) = nbrSimls;
% 	simls(nbrInds(1:end-1), curInd) = nbrSimls(1:end-1);
    
    % assign similarities from current pixel to neighbors
	nNbrInd = numel(nbrInds);
    asgnIdxs = nxtIdx:(nxtIdx+nNbrInd-1);
	rIdxs(asgnIdxs) = curInd;
	cIdxs(asgnIdxs) = nbrInds;
    smVals(asgnIdxs) = nbrSimls;
    nxtIdx = nxtIdx + nNbrInd;
    
    % assign similarities from neighbors to current pixel
    nNbrInd = numel(nbrInds)-1;
    asgnIdxs = nxtIdx:(nxtIdx+nNbrInd-1);
    rIdxs(asgnIdxs) = nbrInds(1:end-1);
    cIdxs(asgnIdxs) = curInd;
    smVals(asgnIdxs) = nbrSimls(1:end-1);
    nxtIdx = nxtIdx + nNbrInd;
end
lastIdx = find(rIdxs, 1, 'last');

fromInds = rIdxs(1:lastIdx);
toInds = cIdxs(1:lastIdx);
simlVals = smVals(1:lastIdx);
numElts = numVecs;

% save('clusTst', 'fromInds', 'toInds', 'simlVals', 'numElts');

end