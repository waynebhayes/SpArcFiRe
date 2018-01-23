function [x, y] = powerlawSpiralXY(params, thBounds)
% params: thOff, k, ir

if nargin < 2 || isempty(thBounds)
    thBounds = [0, 2*pi];
end

lb = thBounds(1); ub = thBounds(2);

theta = [lb:pi/360:ub];
thOff = params(1); k = params(2); ir = params(3);

ftheta = ir .* theta.^(-k);
theta = theta + thOff;

[x, y] = pol2cart(theta, ftheta);

end