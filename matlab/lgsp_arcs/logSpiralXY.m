function [x, y] = logSpiralXY(params, thBounds)
% params: thOff, a, ir

if nargin < 2 || isempty(thBounds)
    thBounds = [0, 2*pi];
end

lb = thBounds(1); ub = thBounds(2);
thOff = params(1); a = params(2); ir = params(3);

% if a > 0
%     theta = [lb:pi/360:ub];
% else
%     theta = [-ub:pi/360:lb];  % :-pi/360:  ???
% end
theta = [lb:pi/360:ub];

ftheta = ir .* exp(-a .* theta);
theta = theta + thOff;

[x, y] = pol2cart(theta, ftheta);


end