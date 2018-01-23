function area = lgspArea(params1, params2, bnds)
% Calculates the area between two log-spiral arcs
% INPUTS:
%   paramsA: parameters [theta-offset, pitch-angle, initial-radius] for
%       first arc
%   paramsB: parameters for second arc
%   bnds: theta-range for which the area is to be computed ("absolute"
%       theta-values; not relative to the theta-offsets of the arcs)

bnds = sort(bnds, 'ascend');

itxTh = lgspIntersection(params1, params2, bnds);
if ~isempty(itxTh)
%     warning('error arcs intersect within interval');
    segBnds = sort([bnds'; itxTh], 'ascend');
    area = 0;
    for segIdx = 1:1:(length(segBnds)-1)
        area = area + abs(...
            areaUnderSegment(params1, segBnds(segIdx), segBnds(segIdx+1))...
            - areaUnderSegment(params2, segBnds(segIdx), segBnds(segIdx+1)));
    end
else
    area = abs(areaUnderSegment(params1, bnds(1), bnds(2)) - ...
        areaUnderSegment(params2, bnds(1), bnds(2)));
end
% check our result against Matlab's numerical integration
% area2 = (1/2) * quad(@(x)(abs(logSpiralFxn2Rev(x', params1).^2 - logSpiralFxn2Rev(x', params2).^2)), bnds(1), bnds(2))

    function segArea = areaUnderSegment(params, ivlStart, ivlEnd)
        thOff = params(1); pa = params(2); ir = params(3);
        if pa == 0
            segArea = ((ivlEnd - ivlStart) * ir^2) / 2;
            return;
        end
        segArea = abs(...
            (((exp(2*pa*(thOff-ivlStart)) - exp(2*pa*(thOff-ivlEnd)))...
            * ir^2) / (4*pa)) );
    end

%     function ivlArea = areaInInterval(ivlStart, ivlEnd)
%         term1 = (((exp(2*pa1*(thOff1-ivlStart)) - exp(2*pa1*(thOff1-ivlEnd))) ...
%             * ir1^2) / pa1);
%         term2 = (((exp(2*pa2*(thOff2-ivlEnd)) - exp(2*pa2*(thOff2-ivlStart))) ...
%             * ir2^2) / pa2);
%         ivlArea = abs((term1 + term2)/4);
%     end

end