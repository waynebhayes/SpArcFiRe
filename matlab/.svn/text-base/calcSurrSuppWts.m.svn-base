function voteWts = calcSurrSuppWts(edgeImg, sdev, strength)
% Calculates the voting strength of each edge point according to the
% "surround suppression" voting scheme described in 
%  T. Pridmore, Y. Kong, and X. Zhang, “An improved Hough transform voting 
%   scheme utilizing surround suppression,” Pattern Recognition Letters, 
%   vol. 30, no. 13, Oct. 2009, pp. 1241-1252.
% INPUTS:
%   edgeImg: grayscale (or binary) image of edge intensities
%   sdev: size of surround suppression filter; larger values mean larger
%    neighborhoods will be considered when considering the texture
%    complexity of an edge point's region
%   strength: surround suppression strength - higher values mean that edges
%    in high-texture-complexity regions will have stronger down-weighting
% OUTPUTS:
%   voteWts: voting strengths for each edge pixel, with the same indexing
%    as the edge intensity image

if length(size(edgeImg)) ~= 2
    error('edgeImg should be a matrix')
end

if nargin < 2 || isempty(sdev)
    sdev = 1;
end

if nargin < 3 || isempty(strength)
    strength = 1;
end

% DoG radius 4*sdev - DoG radius sdev, with the truncation on the latter 
%   adjusted so that it has the same size as the former
dog4 = genGaussMtx(4 * sdev);
dog1 = genGaussMtx(1 * sdev, size(dog4, 1));
wtFxn = dog4 - dog1;
wtFxn = wtFxn .* (wtFxn > 0);
wtFxn = wtFxn / norm(wtFxn, 1);

issVals = conv2(edgeImg, wtFxn, 'same');
% issVals = imfilter(edgeStrs, wtFxn, 'conv', 'replicate');
figure; imagesc(flipud(issVals)); title('issVals');
sslopeVals = atan(issVals ./ edgeImg);
sslopeVals(edgeImg == 0) = atan(inf);

voteWts = cos(sslopeVals) .^ strength;

end 