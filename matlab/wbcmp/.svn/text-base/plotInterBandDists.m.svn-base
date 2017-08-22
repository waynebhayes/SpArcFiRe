function plotInterBandDists(clusMtxs1, clusMtxs2, imgSz, ctrR, ctrC, stgs, clustersMatch, name, resPath)

if nargin < 7 || isempty(clustersMatch)
    clustersMatch = false;
end

if nargin >= 9 && ~isempty(resPath) && isempty(name)
    error('if resPath specified, name must also be specified');
end
if nargin >= 9 && ~isempty(resPath)
    resPath = [resPath filesep];
end

if nargin < 8
    name = '';
end
if nargin < 9
    resPath = [];
end

if clustersMatch && (size(clusMtxs1, 3) ~= size(clusMtxs2, 3))
    error('clustersMatch = true, but cluster sets have different sizes');
end

[params1, bounds1] = fitLogSpiralsToClusters(clusMtxs1, ctrR, ctrC, stgs);
[params2, bounds2] = fitLogSpiralsToClusters(clusMtxs2, ctrR, ctrC, stgs);
lengths1 = calcLgspArcLengths(params1, bounds1);
lengths2 = calcLgspArcLengths(params2, bounds2);
chir1 = sign(params1(:, 2))' * lengths1;
chir2 = sign(params2(:, 2))' * lengths2;
if sign(chir1) ~= sign(chir2)
    warning('detected chirality is different across the two wavebands!');
end
if sign(chir1) > 0
    fprintf('CCW winding direction detected; flipping image\n');
    clusMtxs1 = flipdim(clusMtxs1, 2);
    clusMtxs2 = flipdim(clusMtxs2, 2);
    [params1, bounds1] = fitLogSpiralsToClusters(clusMtxs1, ctrR, ctrC, stgs);
    [params2, bounds2] = fitLogSpiralsToClusters(clusMtxs2, ctrR, ctrC, stgs);
end

lgsp(bounds1(1)+params1(1),params1,bounds1)
lgsp(bounds1(2)+params1(1),params1,bounds1)
lgsp(bounds2(1)+params2(1),params2,bounds2)
lgsp(bounds2(2)+params2(1),params2,bounds2)

cImg1 = showClustersFromMtxs(clusMtxs1, imgSz);
cImg2 = showClustersFromMtxs(clusMtxs2, imgSz);
ovlpImg = zeros([imgSz 3]);
ovlpImg(:, :, 1) = (sum(clusMtxs1, 3) > 0);
ovlpImg(:, :, 3) = (sum(clusMtxs2, 3) > 0);
figure;
subplot(1, 3, 1); imshow(cImg1); title('waveband 1');
subplot(1, 3, 2); imshow(cImg2); title('waveband 2');
subplot(1, 3, 3); imshow(ovlpImg); title('cluster overlap');
if ~isempty(resPath)
    saveas(gca, [resPath name '_001_clusters.png']);
end

cbndImg1 = displayClusterOverlay(zeros(imgSz), clusMtxs1);
cbndImg2 = displayClusterOverlay(zeros(imgSz), clusMtxs2);
% figure; imshow(cbndImg1(:, :, [2 1 3]) + cbndImg2);
figure; imshow(cbndImg1(:, :, [2 1 1]) + cbndImg2(:, :, [1 1 2]));
title([name ' : cluster boundaries'], 'Interpreter', 'none');

colordef black

[avrPlot, thDistPlot, rDistPlot] = ...
    plotOneDirection(params1, bounds1, params2, bounds2, 'waveband 1 to waveband 2');
if ~isempty(resPath)
    figure(avrPlot); set(gcf, 'InvertHardcopy', 'off'); 
    saveas(gcf, [resPath name '_002_arcs1-vs-arcs2.png']);
    figure(thDistPlot); set(gcf, 'InvertHardcopy', 'off'); 
    saveas(gcf, [resPath name '_003_thDist-1vs2.png']);
    figure(rDistPlot); set(gcf, 'InvertHardcopy', 'off'); 
    saveas(gcf, [resPath name '_004_rDist-1vs2.png']);
end
[avrPlot, thDistPlot, rDistPlot] = ...
    plotOneDirection(params2, bounds2, params1, bounds1, 'waveband 2 to waveband 1');
if ~isempty(resPath)
    figure(avrPlot); set(gcf, 'InvertHardcopy', 'off'); 
    saveas(gcf, [resPath name '_005_arcs2-vs-arcs1.png']);
    figure(thDistPlot); set(gcf, 'InvertHardcopy', 'off'); 
    saveas(gcf, [resPath name '_006_thDist-2vs1.png']);
    figure(rDistPlot); set(gcf, 'InvertHardcopy', 'off'); 
    saveas(gcf, [resPath name '_007_rDist-2vs1.png']);
end

colordef white

    function [avrPlot, thDistPlot, rDistPlot] = ...
            plotOneDirection(params1, bounds1, params2, bounds2, dctnName)
        nArcs1 = size(params1, 1);
        nArcs2 = size(params2, 1);

        lineColors = colormap(hsv(256));
        lineColors = lineColors(round(linspace(1, 256, nArcs1)), :);

        displayLgspPlot([params2; params1], [bounds2; bounds1], zeros(256), ...
            ctrR, ctrC, [], [], [ones(nArcs2, 3); lineColors], [], 1);
        title([name ' : distance-measured arcs vs reference arcs (' dctnName ')'],...
            'Interpreter', 'none')
        avrPlot = gcf;

        maxR = sqrt(sum((imgSz/2).^2));
        thDistPlot = figure; hold all
        line([0 maxR], [0 0], 'Color', 'w');
        axis([0 maxR -180 180])
        xlabel('r');
        ylabel('th-dist');
        title([name ' : r vs th-dist to nearest (' dctnName ')'], 'Interpreter', 'none');
        rDistPlot = figure; hold all
        line([0 maxR], [0 0], 'Color', 'w');
        axis([0 maxR -maxR maxR])
        xlabel('r');
        ylabel('r-dist')
        title([name ' : r vs r-dist to nearest (' dctnName ')'], 'Interpreter', 'none');
        for arcIdx = 1:1:nArcs1
            thGran = pi/360;
            thStart = bounds1(arcIdx, 1) + params1(arcIdx, 1);
            thEnd = bounds1(arcIdx, 2) + params1(arcIdx, 1);
            % find start to nearest granularity increment, make sure within bounds
            thStart = thGran * ceil(thStart/thGran);
            thEnd = thGran * floor(thEnd/thGran);

            thVals = thStart:thGran:thEnd;
            rVals = NaN * ones(size(thVals));
            thDistVals = NaN * ones(size(thVals));
            rDistVals = NaN * ones(size(thVals));
            for thIdx = 1:length(thVals)
                curTh = thVals(thIdx);
                curR = lgsp(curTh, params1(arcIdx, :), bounds1(arcIdx, :));
                assert(~isnan(curR));
                rVals(thIdx) = curR;
                otherTh = NaN * zeros(1, nArcs2);
                if clustersMatch
                    otherTh(arcIdx) = ...
                        invLgsp(curR, params2(arcIdx, :), bounds2(arcIdx, :));
                else
                    for ii=1:1:nArcs2
                        otherTh(ii) = ...
                            invLgsp(curR, params2(ii, :), bounds2(ii, :));
                    end
                end

                thDists = (mod(otherTh, 2*pi) - mod(curTh, 2*pi));
                thDctns = sign(thDists);
                thDists = abs(thDists);
                thruZero = (thDists > pi);
                thDists(thruZero) = 2*pi - thDists(thruZero);
                thDctns(thruZero) = -thDctns(thruZero);
                [mVal, mIdx] = min(thDists);
                thDistVals(thIdx) = mVal * thDctns(mIdx);

                otherR = NaN * zeros(1, nArcs2);
                if clustersMatch
                    otherR(arcIdx) = ...
                        lgsp(curTh, params2(arcIdx, :), bounds2(arcIdx, :), true);
                else
                    for ii=1:1:nArcs2
                        otherR(ii) = ...
                            lgsp(curTh, params2(ii, :), bounds2(ii, :), true);
                    end
                end
                rDists = (curR - otherR);
                [mVal, mIdx] = min(abs(rDists));
                rDistVals(thIdx) = rDists(mIdx);
            end

            figure(thDistPlot);
            plot(rVals, thDistVals * (180/pi), 'Color', lineColors(arcIdx, :));

            figure(rDistPlot);
            plot(rVals, rDistVals, 'Color', lineColors(arcIdx, :));
        end

        h=figure;
        copyobj(get(thDistPlot,'children'),h);
        axis([0 maxR -15 15])
        thDistPlot = h;

        h=figure;
        copyobj(get(rDistPlot,'children'),h);
        axis([0 maxR -maxR/8 maxR/8])
        rDistPlot = h;
    end

    function r = lgsp(thv, params, bounds, matchRev)
        diffTol = 10^-12;
        if nargin < 3
            bounds = [-inf inf];
        end
        if nargin < 4 || isempty(matchRev)
            matchRev = false;
        end
        if matchRev
%             fprintf('thv: %2.4f\t bounds+thOff: %s\n', thv, mat2str(bounds + params(1)));
%             revDiff = fix((bounds(1)+params(1))/(2*pi)) - fix(thv/(2*pi));
%             thv = thv + (2*pi) * revDiff;
            while thv < bounds(1)+params(1)
                thv = thv + 2*pi;
            end
            while thv > bounds(2)+params(1)
                thv = thv - 2*pi;
            end
%             fprintf('revDiff: %2.2f\t new thv: %2.4f\n', revDiff, thv);
        end
        thOff = params(1);
        a = params(2);
        ir = params(3);
        if ((thv - thOff) < bounds(1) - diffTol) ||...
                ((thv - thOff) > bounds(2) + diffTol)
            r = NaN;
        return;
        end
        r = ir * exp(-a * (thv - thOff));
    end

    function thv = invLgsp(rv, params, bounds)
        diffTol = 10^-12;
        if nargin < 3
            bounds = [-inf inf];
        end
        thOff = params(1);
        a = params(2);
        ir = params(3);
        if (a == 0) && (rv ~= ir)
            % params give a circle with radius different than given r-value
            thv = NaN;
            return;
        end
        thv = (-log(rv/ir))/a + thOff;
        if thv < (bounds(1) + thOff - diffTol) ||...
                thv > (bounds(2) + thOff + diffTol)
            % outside the bounds of the other log spiral
            thv = NaN;
            return
        end
%         % theta value was the difference from the theta-offset, now we make
%         % it the theta-value from normal polar coordinates
%         thv = thv + thOff;
    end

end