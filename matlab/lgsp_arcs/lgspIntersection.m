function [itxTh, itxR] = lgspIntersection(paramsA, paramsB, bnds, useStrict)
% Finds all intersection points for the given log-spirals, for theta-values
% in the given range.
% INPUTS:
%   paramsA: parameters [theta-offset, pitch-angle, initial-radius] for
%       first arc
%   paramsB: parameters for second arc
%   bnds: theta-interval to look for intersections (must be finite, or
%       there could be an infinite number of intersections to enumerate!)
%   useStrict: true if the intersection points must be within the given
%       theta-range for both arcs, false if the intersection points must be
%       within the given theta range for at least one arc
% OUTPUTS:
%   itxTh: theta-values of the intersections, if any
%   itxR: r-values of the intersections, if any

if nargin < 4 || isempty(useStrict)
    useStrict = true;
end

bnds = sort(bnds, 'ascend');

% ensure that params1 has the smaller pitch angle
% if params1(2) < params2(2)
% % if params1(2) > params2(2)
%     temp = params1;
%     params1 = params2;
%     params2 = temp;
% end

% thOff1 = params1(1); pa1 = params1(2); ir1 = params1(3);
% thOff2 = params2(1); pa2 = params2(2); ir2 = params2(3);

arcApa = paramsA(2); arcBpa = paramsB(2);
arcAir = paramsA(3); arcBir = paramsB(3);

if (arcApa == arcBpa) || (arcAir == 0) || (arcBir == 0)
    itxR = [];
    itxTh = []; 
    return;
end

itxThAB = itxOnArc(paramsA, paramsB);
if ~isempty(itxThAB)
    itxRhAB = logSpiralFxn2Rev(itxThAB, paramsB);
else
    itxRhAB = [];
end
itxThBA = itxOnArc(paramsB, paramsA);
if ~isempty(itxThBA)
    itxRhBA = logSpiralFxn2Rev(itxThBA, paramsA);
else
    itxRhBA = [];
end

% itxR1 = logSpiralFxn2Rev(itxTh, params1);
% itxR2 = logSpiralFxn2Rev(itxTh, params2);

% itxR1
% itxR2
% round((itxR1 - itxR2) * 1000)
% assert(sum(round((itxR1 - itxR2) * 1000) == 0) == numel(itxR1), ...
%     'r-values for intersection theta-values were different');

itxRhABrnd = round(itxRhAB * 1024) / 1024;
itxRhBArnd = round(itxRhBA * 1024) / 1024;
if useStrict
    % easier to use r-values since they're unambiguous and there can be at
    % most one arc-intersection at a given r-value (when there are two
    % arcs)
    [temp, ia, ib] = intersect(itxRhABrnd, itxRhBArnd);
else
    [temp, ia, ib] = union(itxRhABrnd, itxRhBArnd);
end
itxTh = [itxThAB(ia); itxThBA(ib)];
itxR = [itxRhAB(ia); itxRhBA(ib)];

% itxTh = itxThBA;
% itxR = itxRhBA;

% itxTh = [itxThAB; itxThBA];
% itxR = [itxRhAB; itxRhBA];
% [itxTh itxR]
% union


% itxTh = itxThAB;
% itxR = itxRhAB;

% xR = exp(...
%     ((pa1 * pa2) * (thOff2 - thOff1) - pa2 * log(ir1) + pa1 * log(ir2))...
%     / (pa1 - pa2) )
% 
% xTh = ((pa1 * thOff1) - (pa2 * thOff2) + log(ir1/ir2)) / (pa1 - pa2)


% theta = [-2*pi:pi/360:6*pi]';
% f1theta = logSpiralFxn2Rev(theta, params1);
% f2theta = logSpiralFxn2Rev(theta, params2);
% figure; hold on; plot(theta, f1theta); plot(theta, f2theta); hold off
% ir1
% ir2
% log(ir1)
% log(ir2)

% imgSz = [256 256];
% bg = zeros(imgSz);
% plotParams = [paramsA; paramsB];
% plotBounds = repmat([bnds(1) - paramsA(1), bnds(2) - paramsB(1)], 2);
% displayLgspPlot(plotParams, plotBounds, bg, imgSz(1) / 2, imgSz(2) / 2, [], [], [], 1)

% hold on
% [itxX, itxY] = pol2cart(itxThAB, itxRhAB);
% scatter(itxX, itxY, 'ms');
% [itxX, itxY] = pol2cart(itxThBA, itxRhBA);
% scatter(itxX, itxY, 'co');
% [itxX, itxY] = pol2cart(itxTh, itxR);
% scatter(itxX, itxY, 'y*');
% hold off
% th = [-2*pi:pi/360:2*pi]';
% figure; hold on; plot(th, log(logSpiralFxn2Rev(th, params1))); plot(th, log(logSpiralFxn2Rev(th, params2))); hold off

    function itxTh = itxOnArc(params1, params2)
    % absolute theta-values where the first arc, within the theta-range
    % specified by the parameter "bnds," intersects with the second arc
    % (anywhere on the second arc)
    
    itxTh = [];

    curShift = 0;
    curItxTh = itxForShift(params1, params2, curShift);
    while (curItxTh > bnds(1)) && (curItxTh < bnds(2))
        itxTh = [itxTh curItxTh];
        curShift = curShift - 2*pi;
        curItxTh = itxForShift(params1, params2, curShift);
    end
    curShift = 2*pi;
    curItxTh = itxForShift(params1, params2, curShift);
    while (curItxTh > bnds(1)) && (curItxTh < bnds(2))
        itxTh = [itxTh curItxTh];
        curShift = curShift + 2*pi;
        curItxTh = itxForShift(params1, params2, curShift);
    end 
    
    itxTh = itxTh';
    
    end

    function itxTh = itxForShift(params1, params2, thOff1Shift)
        thOff1 = params1(1); pa1 = params1(2); ir1 = params1(3);
        thOff2 = params2(1); pa2 = params2(2); ir2 = params2(3);
        itxTh = ...
            ((pa1 * (thOff1 + thOff1Shift)) - (pa2 * thOff2) + log(ir1/ir2))...
            / (pa1 - pa2);
    end

end