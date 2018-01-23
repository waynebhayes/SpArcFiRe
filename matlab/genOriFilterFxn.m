function filter = genOriFilterFxn(theta, radius, hilbt)
% Generates an orientation field filter matrix described in the PhD thesis
%  "Inferring Galaxy Morphology Through Texture Analysis" (K. Au 2006).
%  The filter is a 1D LoG filter extended in 2D along an angle theta, such
%  that the filter response is strongest for that angle.
% INPUTS: 
%   theta: angle of the filter orientation
%   radius: optional parameter specifying where to truncate the filter 
%    values; matrix size will be 2*radius+1
%   hilbt: whether to use the Hilbert transform of the filter instead

if nargin < 1 || isempty(theta)
    error('theta not specified')
end

if nargin < 2 || isempty(radius)
    radius = 5;
%     radius = 7;
end

if nargin < 3 || isempty(hilbt)
    hilbt = false;
end

% if nargin < 3 || isempty(lightDark)
%     lightDark = 1;
% end

% if nargin < 3 || isempty(sigma)
%     sigma = 2.5;
% end

if radius <= 0
    warning('invalid size of %d, resetting to 5', radius)
    radius = 5;
end

% [cVals, rVals] = meshgrid([-radius:1:radius], [-radius:1:radius]);
msize = 2*ceil(radius)+1;
% range = pi;
range = pi * ((2*radius)/(2*radius + 1)); % sample [-pi, pi], in pixel middles
% range = (3*pi)/4
[cVals, rVals] = meshgrid(linspace(-range, range, msize), linspace(-range, range, msize));
% [cVals, rVals] = meshgrid(linspace(-pi/2, pi/2, msize), linspace(-pi/2, pi/2, msize));
% cVals;
dVals = rVals * cos(theta) + cVals * sin(theta);
% figure; imagesc(dSquared); axis image; title('dSquared');
dSquared = dVals .^2;
filter = (2 / sqrt(3)) * (pi^(-1/4)) * (1 - dSquared) .* exp(-dSquared/2);
if hilbt
    filterH = imag(-cmhf(dVals)) / sqrt(3);
    filterH = filterH .* (sum(abs(filter(:))) / sum(abs(filterH(:))));
    filter = filterH;
end
% % filter = mhf(dVals);
% fprintf('filter change!\n');
% normConst = mhf(0) / real(cmhf(0));
% normConst = 1;
% if lightDark > 0
%     filter = real(cmhf(dVals)) * normConst;
% else
%     filter = imag(cmhf(dVals)) * normConst;
% end
% % filter = real(cmhf(dVals)) * normConst;
% % sum(sum(abs(filter - filterNew))) /numel(filter)
% % filter = imag(cmhf(dVals)) * normConst;
% filter = filter - (sum(filter(:)) / numel(filter)); %%%%%%
% if matchFilt
%     sigma = 1.5;
%     gwind = (1 / sqrt(2*pi*(sigma^2))) * exp((-1/(2*(sigma^2))) * (cVals .^ 2 + rVals .^ 2));
% else
%     gwind = (1 / sqrt(2)) * exp((-1/2) * (cVals .^ 2 + rVals .^ 2));
% end
% sigma = pi/2;
% sigma = 1.315975053682;
% sigma = 1.171378597280;
% sigma = 1.171388779811;
% sigma = 1.171932801054609100;
sigma = range / 2;
% sigma = pi/2
% gwind = (1 / sqrt(2*pi*(sigma^2))) * exp((-1/(2*(sigma^2))) * (cVals .^ 2 + rVals .^ 2));
% gwind = (1 / sqrt(2*pi*(sigma^2))) * exp((-1/(2*(sigma^2))) * (cVals .^ 2));
wVals = rVals * cos(theta + pi/2) + cVals * sin(theta + pi/2);
% imagesc(wVals); axis image; figure;
gwind = (1 / sqrt(2*pi*(sigma^2))) * exp((-1/(2*(sigma^2))) * (wVals .^ 2));
% imagesc(gwind); axis image; figure;
% figure; imagesc(gwind); axis image; title('gwind');
% [cVals, rVals] = meshgrid(linspace(-pi/2, pi/2, msize), linspace(-pi/2, pi/2, msize));
% gwind = (1 / sqrt(2)) * exp((-1/(2*(2)^2)) * (cVals .^ 2 + rVals .^ 2)); % TEMP

% gwind = (1 / sqrt(2)) * exp(-0.5*(rVals.^2+cVals.^2));
% gwind = (1 / sqrt(2*pi*(sigma^2))) * exp((-1/(2*(sigma^2))) * (rVals .^ 2 + cVals .^2));

filter = filter .* gwind;  
% filter = filter - (sum(filter(:)) / numel(filter));
% filter = filter - sum(filter(:));
% filter = filter - sum(filter(:).^2);
% filter = filter - (sum(filter(:).^2) / numel(filter));
filter = filter / sqrt(sum(filter(:).^2));

%figure; imagesc(filter); axis image;
% filter = filter .* genGaussMtx(sigma, 2*radius+1);
% filter = filter / sum(sum(filter));

end