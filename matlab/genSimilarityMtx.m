function simls = genSimilarityMtx(img, oriMtx, stgs, simlCutoff)
% Generates a sparse pixel-to-pixel similarity matrix using an orientation
%   field derived from the input image.
% INPUTS:
%   img: image used to generate the orientation field
%   oriMtx: orientation field generated from the image
%   stgs: structure containing algorithm settings (see getDefaultSettings.m)
%   simlCutoff: minimum value for a similarity score to be nonzero in the
%       similarity matrix (this is often related to, but not necessarily
%       the same as, the HAC stopping threshold in the 'stgs' struct)
% OUTPUTS:
%   simls: the sparse pixel-to-pixel similarity matrix

if nargin < 3 || isempty(simlCutoff)
    simlCutoff = -inf;
end

[fromInds, toInds, simlVals, numElts] = ...
    calcPxlSimilarities(img, oriMtx, stgs, simlCutoff);

simls = sparse(fromInds, toInds, simlVals, numElts, numElts);

% simlVec = sort(nonzeros(simls));
% figure; hist(simlVec, 100);

end