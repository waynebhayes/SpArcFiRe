function fimgs = genOriFilteredImgs(imgMtx, stgs, hilbt)
% Computes the filter responses to the 9 orientation filters described in
%  "Inferring Galaxy Morphology Through Texture Analysis" (K. Au 2006).
% INPUTS:
%   imgMtx: image to filter (as a matrix of doubles)
%   stgs: structure containing algorithm settings (see settings.m)
%   hilbt: whether to use the Hilbert transform of the filter instead
% OUTPUTS:
%   fimgs: filtered images, where fimgs(:, :, k) is the image response to
%    the orientation filter at (k*pi)/9

if nargin < 3 || isempty(hilbt)
    hilbt = false;
end

fimgs = zeros([size(imgMtx) 9]);

for angIdx=0:1:8
    angle = (angIdx * pi) / 9;
    oriFilt = genOriFilterFxn(angle, [], hilbt);
    fimgs(:, :, angIdx + 1) = conv2(imgMtx, oriFilt, 'same');
end

end