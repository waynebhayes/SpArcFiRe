function [angle, halfLength, sumSqErr] = ...
    fitBarUsingLine(img, ctrR, ctrC, allowArcBeyond2pi, plotFlag)
% not used in current workflow
% halfLength is not a free parameter

if nargin < 5 || isempty(plotFlag)
    plotFlag = false;
end

[r, c, brt] = find(img);
xVals = c - ctrC;
yVals = ctrR - r;
wtdXmean = sum(xVals .* brt) / sum(brt);
wtdYmean = sum(yVals .* brt) / sum(brt);
muDist = sqrt(wtdXmean^2 + wtdYmean^2);
[thVals, rVals] = cart2pol(xVals, yVals);
[angle, sumSqErr] = ...
    lsqnonlin(...
    @(angle)(sqrt(brt) .* (angleDevFromLine(thVals, angle) .* rVals)),...
    [0], [], [], ...
    optimset('MaxFunEvals', 1000, 'Display', 'off'));
halfLength = max(abs(xVals .* cos(angle) + yVals .* sin(angle)));

if plotFlag
    [t1, t2, sumSqErrLgsp] = fitLogSpiral(img, ctrR, ctrC, ...
        allowArcBeyond2pi);
    
    figure
    imagesc([1 - ctrC, size(img, 2) - ctrC], ...
        [1 - ctrR, size(img, 1) - ctrR], img);
    axis image; colormap gray
    xLim = halfLength * cos(angle);
    yLim = halfLength * sin(angle);
    line([xLim -xLim], [-yLim, yLim], 'Color', 'g');
    title(sprintf(...
        ['muDist = %2.4f, sumSqErr = %2.4f, angle = %2.4f, halfLength = %2.4f\n'...
        'MSE err ratio vs log-spiral: %2.4f (%2.4f vs %2.4f)'], ...
        muDist, sumSqErr, angle, halfLength, sqrt(sumSqErr/sum(brt))/sqrt(sumSqErrLgsp/sum(brt)), ...
        sqrt(sumSqErr/sum(brt)), sqrt(sumSqErrLgsp/sum(brt))));  
end

    function dists = angleDevFromLine(thVals, lineAng)
        % rotate the angles so that line deviations are deviations from a
        % vertical line
        lineAng = mod(lineAng, 2*pi);
        thVals = mod(thVals - lineAng + pi/2, 2*pi);
        
        dists = inf * ones(size(thVals));
        upperQuadrants = thVals < pi;
        dists(upperQuadrants) = abs(thVals(upperQuadrants) - pi/2);
        dists(~upperQuadrants) = abs(thVals(~upperQuadrants) - (3*pi)/2);
%         figure; hist(dists);
    end

end