function [oriVecs, strengths, oris, ofldLvls] = genOriField(imgMtx, stgs, lightDark)
% Generates an orientation field as described in the PhD thesis
% "Inferring Galaxy Morphology Through Texture Analysis" (K. Au 2006).
% All steps (3-level orientation field extraction, merging, and de-noising)
% are performed.
% INPUTS:
%   imgMtx: the m-by-n image mxn for which an orientation field is to be
%       generated
%   stgs: structure containing algorithm settings (see settings.m)
%   lightDark: 1 to obtain responses only to bright linear regions, -1 to
%       obtain responses only to dark linear regions, 0 for both
% OUTPUTS:
%   oriVecs: the final m-by-n-by-2 orientation field matrix
%   strengths: m-by-n matrix of orientation field strengths (vector
%       magnitudes)
%   oris: m-by-n matrix of orientation field directions, in radians, from 
%       0 to pi since orientation field vectors have an orientation but not
%       a direction 

if nargin < 3 || isempty(lightDark)
    lightDark = 1;
end

numOrientationFieldLevels = stgs.numOrientationFieldLevels;

ofldLvls = cell(1, numOrientationFieldLevels);
for ii=1:1:numOrientationFieldLevels
    ofldLvls{ii} = genSglLvlOriFld(imresize(imgMtx, 1/(2^(numOrientationFieldLevels-ii))), stgs, lightDark);
end
oriVecs = ofldLvls{1};
for ii=2:1:numOrientationFieldLevels
    oriVecs = mergeOriFields(oriVecs, ofldLvls{ii});
end

oriVecs = deNoiseOrientationField(oriVecs);
strengths = oriFieldStrengths(oriVecs);
oris = oriFieldDirections(oriVecs);

% figure; imagesc(strengths); axis image; colorbar; title(strengths);

end