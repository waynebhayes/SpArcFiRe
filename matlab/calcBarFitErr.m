function err = calcBarFitErr(img, bCtrR, bCtrC, bAngle, bHalfLen)
% Calculates the sum of squared Euclidean distances from each pixel to the
% (closest point on the) given bar, weighted by pixel brightness

[r, c, brt] = find(img);
xVals = c - bCtrC;
yVals = bCtrR - r;

ctrDists = sqrt(xVals.^2 + yVals.^2);
angles = atan2(yVals, xVals);

% Measure distance along an orthogonal coordinate system with one axis
% parallel to the bar; distance along the bar-parallel axis is zero until
% going beyond the length of the bar
perpDists = ctrDists .* sin(angles - bAngle);
alongDists = max(0, ctrDists .* abs(cos(angles - bAngle)) - bHalfLen);

barDists = perpDists.^2 + alongDists.^2; % no sqrt since we want squared dist

err = sum(sum(brt .* barDists));

end