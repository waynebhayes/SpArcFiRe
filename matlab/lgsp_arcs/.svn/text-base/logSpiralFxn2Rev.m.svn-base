function [rho, thetaAdj] = logSpiralFxn2Rev(theta, params)

% Can calculate log-spiral function for multiple parameter sets
% Theta is a column vector
% Params has one parameter set per row
% Result matrix rows correspond to theta-values, columns correspond to
% parameter values

nParamSets = size(params, 1);
nPts = length(theta);
thOff = params(:, 1); pitchAngle = params(:, 2); initRadius = params(:, 3);

% thetaAdj = mod(theta - thOff, 2*pi);
thetaAdj = repmat(theta, 1, nParamSets) - repmat(thOff', nPts, 1);
% figure; scatter(thetaAdj, thetaAdj); title('thetaAdj');
% thetaAdj = mod(thetaAdj - thStart, (thEnd - thStart)) + thStart;

rho = repmat(initRadius', nPts, 1) .* exp(repmat(-pitchAngle', nPts, 1) .* thetaAdj);

end