function analyzeGaussFit(xcMtx, lags, nTheta, name)

if nargin < 4
    name = [];
end

if isempty(name)
    figVis = 'on';
else
    figVis = 'off';
end

radTickPcts = [0 0.25 0.5 0.75 1];
lagTickPcts = linspace(0, 1, 7);
screen_size = get(0, 'ScreenSize');

minXcorr = min(xcMtx(:));
if minXcorr < 0
    warning('cross-correlation matrix has negative values (smallest is %e)', minXcorr);
    xcMtx(xcMtx <= 0) = 0;
end
xcMtx(isnan(xcMtx)) = 0; % zero out areas coming from places with no data values

nLagVals = size(xcMtx, 1);
nRadVals = size(xcMtx, 2);

degPerPixel = (360/nTheta);
maxLagDeg = ((nLagVals-1)/2) * degPerPixel;

xcRadNorm = xcMtx;
for radIdx = 1:1:nRadVals;
    curMin = min(xcRadNorm(:, radIdx));
    if curMin > 0
        xcRadNorm(:, radIdx) = xcRadNorm(:, radIdx) - curMin;
    end
end
xcRadNorm = xcRadNorm ./ repmat(max(xcRadNorm, [], 1), size(xcRadNorm, 1), 1);

makeXcContourPlot(xcMtx, 10, 'cross-correlation values')
if ~isempty(name)
    saveas(gca, [name '_007-xcorrContours.png']);
end
makeXcContourPlot(xcRadNorm, 5, 'xcorr (independently normalized radii)')
if ~isempty(name)
    saveas(gca, [name '_008-xcorrContoursColNormalized.png']);
end
makeXcColorPlot(xcMtx, 'cross-correlation values')
makeXcColorPlot(xcRadNorm, 'xcorr (independently normalized radii)')

% offsets = inf * ones(1, nRadVals);
% startCoeffs = [0 -1 -1 0];
% fxnVals = NaN * ones(size(xcMtx));
% fitOptions = optimset('MaxFunEvals', 5000, 'MaxIter', 2000);
% for ii=1:1:nRadVals
%     params = lsqcurvefit(@xcParabola, startCoeffs, lags', xcMtx(:, ii),...
%         [], [], fitOptions);
%     offsets(ii) = params(1);
%     fxnVals(:, ii) = xcParabola(params, lags');
% end
% figure; plot(offsets);
% figure; imagesc(fxnVals); axis image; colormap gray

lagMtx = repmat(lags', 1, nRadVals);

xcWts = (xcMtx) .^ 2;
means = sum(xcWts .* lagMtx, 1) ./ sum(xcWts, 1);
sdevs = sqrt(sum(xcWts .* (lagMtx - repmat(means, nLagVals, 1)).^2) ./ sum(xcWts, 1));
% plotMax(xcMtx, means - min(lags) + 1, sdevs, [], '');
plotMax(xcMtx, means - min(lags) + 1, sdevs, [], '', 'b.'); colormap hot
title('cross-correlation overlay, wts = xcorr^2');
if ~isempty(name)
    saveas(gca, [name '_009-xcorrOvlWtPow2.png']);
end
plotMax(xcRadNorm, means - min(lags) + 1, sdevs, [], '', 'b.'); colormap hot
title('column-normalized cross-correlation overlay, wts = xcorr^2');
if ~isempty(name)
    saveas(gca, [name '_010-xcorrOvlColNormWtPow2.png']);
end

xcWts = (xcMtx) .^ 10;
means = sum(xcWts .* lagMtx, 1) ./ sum(xcWts, 1);
sdevs = sqrt(sum(xcWts .* (lagMtx - repmat(means, nLagVals, 1)).^2) ./ sum(xcWts, 1));
% plotMax(xcMtx, means - min(lags) + 1, sdevs, [], '');
plotMax(xcMtx, means - min(lags) + 1, sdevs, [], '', 'b.'); colormap hot
title('cross-correlation overlay, wts = xcorr^{10}');
if ~isempty(name)
    saveas(gca, [name '_011-xcorrOvlWtPow10.png']);
end
plotMax(xcRadNorm, means - min(lags) + 1, sdevs, [], '', 'b.'); colormap hot
title('column-normalized cross-correlation overlay, wts = xcorr^{10}');
if ~isempty(name)
    saveas(gca, [name '_012-xcorrOvlColNormWtPow10.png']);
end

figure('Visible', figVis); plot(sdevs)

    function makeXcContourPlot(xcVals, nContours, titleStr)
        if nargin < 2
            titleStr = '';
        end
        figure('Visible', figVis); hold on; 
        % xcVals(:, lags == 0) = 0;
        imagesc(xcVals); axis image; colormap gray
        line([0 nRadVals], find(lags==0) * [1 1])
        contour(xcVals, nContours, '-r');
        title(titleStr);
        setXcPlotAxes()
    end

    function makeXcColorPlot(xcVals, titleStr)
        if nargin < 2
            titleStr = '';
        end
        figure('Visible', figVis); hold on; 
        imagesc(xcVals); axis image;
        line([0 nRadVals], find(lags==0) * [1 1], 'Color','k')
%         contour(xcVals, 10, '-k');
        title(titleStr);
        setXcPlotAxes()
    end

    function setXcPlotAxes()
        xlabel('rIdx');
        ylabel('lag');
        set(gca, 'XTick', radTickPcts * nRadVals); 
        set(gca,'YDir','normal')
        set(gca, 'YTick', lagTickPcts * nLagVals);
        set(gca, 'YTickLabel', linspace(-maxLagDeg, maxLagDeg, length(lagTickPcts)));
        set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    end

    function plotMax(xcVals, idx, sd, sdThres, maxName, markerType)
        if nargin < 2
            idx = [];
        end
        if nargin < 3
            sd = [];
        end
        if nargin < 4 || isempty(sdThres)
            sdThres = inf;
        end
        if nargin < 5
            maxName = [];
        end
        if nargin < 6 || isempty(markerType)
            markerType = 'r.';
        end
        if isempty(sd) && ~isempty(sdThres)
            error('if sdThres specified, sd must also be specified');
        end
        hfig = figure('Visible', figVis); hold on; 
        imagesc(xcVals); axis image; colormap gray
        line([0 nRadVals], find(lags==0) * [1 1], 'Color', 'b')
        xlabel('radius-index');
        ylabel('lag');
        set(gca,'YDir','normal')
        set(gca, 'XTick', radTickPcts * nRadVals)
        set(gca, 'YTick', lagTickPcts * nLagVals);
%         lbls = get(gca, 'YTickLabel');
%         set(gca, 'YTickLabel', lags(str2double(mat2cell(lbls, ones(size(lbls, 1), 1), size(lbls, 2)))));
        set(gca, 'YTickLabel', linspace(-maxLagDeg, maxLagDeg, length(lagTickPcts)));
        radVals = 1:nRadVals;
        incl = (sd <= sdThres);
        if ~isempty(idx) && isempty(sd)
            scatter(radVals(incl), idx(incl), markerType, 'SizeData', 5^2);
        elseif ~isempty(idx) && ~isempty(sd)
            errorbar(radVals(incl), idx(incl), sd(incl), markerType);
        end
        if ~isempty(maxName)
            title(maxName);
        end
        set(hfig, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    end

    function y = xcParabola(params, xdata)
        assert(length(params) == 4);
        xoff = params(1);
        coeffs = params(2:4);
        y = coeffs(1) * (xdata - xoff).^2 + coeffs(2) * (xdata - xoff) + coeffs(3);
    end

end