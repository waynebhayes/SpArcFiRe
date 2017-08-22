function [orivecs, energy] = genSglLvlOriFld(img, stgs, lightDark)
% Generates one level of an orientation field as described in the PhD 
% thesis "Inferring Galaxy Morphology Through Texture Analysis" 
% (K. Au 2006).
% INPUTS:
%   imgMtx: the MxN image for which an orientation field is to be generated
%   stgs: structure containing algorithm settings (see settings.m)
%   lightDark: 1 to obtain responses only to bright linear regions, -1 to
%       obtain responses only to dark linear regions, 0 for both
% OUTPUTS:
%   oriVecs: the MxNx2 orientation field matrix

% if nargin < 3 || isempty(darkPhaseWt)
%     darkPhaseWt = 1/2;
% end

if nargin < 3 || isempty(lightDark)
    lightDark = 1;
end

fimgs = genOriFilteredImgs(img, stgs);
if lightDark > 0
    fimgs(fimgs < 0) = 0;
elseif lightDark < 0
    fimgs(fimgs > 0) = 0;
end

energy = fimgs .^ 2;
for angIdx=0:1:8;
    energy(:, :, angIdx+1) = ...
        energy(:, :, angIdx+1) * exp(1i * 2 * (angIdx * pi)/9);
end
energy = sum(energy, 3);

strengths = abs(energy); 
% figure; imagesc(strengths); axis image

% strengths = stgs.propWtEnergyWithImg * (strengths .* img) + (1 - stgs.propWtEnergyWithImg) * strengths;
% strengths(energy < 0) = stgs.darkPhaseWt * strengths(energy < 0);

% strengths = abs(energy);
directions = angle(energy) / 2;
% directions = mod(directions, pi);

orivecs = zeros([size(energy) 2]);
orivecs(:, :, 1) = strengths .* cos(directions);
orivecs(:, :, 2) = strengths .* sin(directions);

end