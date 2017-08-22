function findXcorrMax(xCorrMtx, lags)

nLagVals = size(xCorrMtx, 1);
nRadVals = size(xCorrMtx, 2);

radTickPcts = [0 0.25 0.5 0.75 1];

minXcorr = min(xCorrMtx(:));
if minXcorr < 0
    warning('cross-correlation matrix has negative values (smallest is %e)', minXcorr);
    xCorrMtx(xCorrMtx <= 0) = 0;
end
xCorrMtx(isnan(xCorrMtx)) = 0; % zero out areas coming from places with no data values

% figure; hold all;
% for radIdx = 50:25:200
%     plot(lags, xCorrMtx(:, radIdx));
% end


[maxVal, maxIdx] = max(xCorrMtx, [], 1);
plotMax(maxIdx, [], 'arrayMax');
isMax = (xCorrMtx == repmat(maxVal, nLagVals, 1));
allMaxImg = xCorrMtx; allMaxImg(isMax) = 0; 
allMaxImg = repmat(allMaxImg, [1 1 3]); 
allMaxImg((lags==0), :, 1:2) = 0;
allMaxImg(:, :, 1) = allMaxImg(:, :, 1) + isMax;
figure; imshow(allMaxImg); title('allMax'); axis on
set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
% saveas(gca, '005_allMax.png');

xcnlz = xCorrMtx;
xcnlz(lags == 0, :) = 0;
[val, idx] = max(xcnlz, [], 1);
plotMax(idx, [], 'arrayMax (excluding lag0)');

maxXCorr = max(xCorrMtx, [], 1);
brtFrac = 0.75;
isPeak = (xCorrMtx >= (brtFrac * repmat(maxXCorr, nLagVals, 1)));
%         figure; imshow(isPeak); title('isPeak');
% plotMainPeakMeanSd(isPeak, false, ...
%     'mean/sd, peak by xcorr-fraction');
plotMainPeakMeanSd(isPeak, true, ...
    'mean/sd, peak by xcorr-fraction (normalized)');
% saveas(gca, '007_meanSdPeakByXcorrFraction.png');

%-- maybe try this again later
% [mxcV, mxcIdx] = max(xCorrMtx, [], 1);
% idxs = repmat([1:nLagVals]', 1, nRadVals);
% isPeak = (idxs >= repmat(mxcIdx - 5, nLagVals, 1)) & ...
%     (idxs <= repmat(mxcIdx + 5, nLagVals, 1));
% % plotMainPeakMeanSd(isPeak, false, 'mean/sd, distance cutoff from max');
% plotMainPeakMeanSd(isPeak, true, 'mean/sd, distance cutoff from max (normalized)');

% plotLagDctnProportions(lags, xCorrMtx, 'lag direction proportions: raw xcorr scores')
% plotLagDctnProportions(lags, xCorrMtx .* repmat(normpdf(lags, 0, 2)', 1, nRadVals),...
%     'lag direction proportions: gaussWtd')

% [msS, msN] = meanShift(1);
% [msS, msN] = meanShift(5);
displayMeanShift(2);
displayMeanShift(5);
% displayMeanShift(10);
% [msS, msN] = meanShift(25);
% [mV, mI] = max(msS, [], 1);
% plotMax(mI, 'msS bw=25');
% [mV, mI] = max((msS > 0) .* xCorrMtx, [], 1);
% plotMax(mI, 'ms peaks bw=25');
% [msS, msN] = meanShift(50);
% [msS, msN] = meanShift(256);
% save('ms', 'msS', 'msN', 'basins');

    function plotLagDctnProportions(lags, wtdXcorr, plotTitle)
        if size(lags, 1) < size(lags, 2)
            lags = lags';
        end
        lagMtx = repmat(lags, 1, nRadVals);
        negLagSums = sum(wtdXcorr .* (lagMtx < 0), 1);
        posLagSums = sum(wtdXcorr .* (lagMtx > 0), 1);
        dctns = sign(posLagSums - negLagSums);
        figure; stem(dctns .* ((max(posLagSums, negLagSums) ./ min(posLagSums, negLagSums)) - 1));
        title(plotTitle);
    end

    function displayMeanShift(bwidth)
        [msS, msN, basins] = meanShift(bwidth);
        [mV, mI] = max(msS, [], 1);
        plotMax(mI, [], sprintf('msS bw=%d', bwidth));
        [mV, mI] = max((msS > 0) .* xCorrMtx, [], 1);
        plotMax(mI, [], sprintf('ms peaks bw=%d', bwidth));
        mspImg = repmat(xCorrMtx .* (msN <= 0), [1 1 3]);
        mspImg((lags==0), :, 1:2) = 0;
        mspImg(:, :, 1) = mspImg(:, :, 1) + (msN > 0);
        figure; imshow(mspImg); title(sprintf('ms all peaks bw=%d', bwidth));
        axis on; set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
%         saveas(gca, sprintf('006_msAllPeaks_bw%d.png', bwidth));
        
%         [mV, mI] = max(msS, [], 1);
%         inMainPeak = (basins == repmat(mI, size(msS, 1), 1));
%         plotMainPeakMeanSd(inMainPeak, sprintf('peak-basin means bw=%d', bwidth));
        
        figure; imshow(maxBasins(basins, msS)); axis on; axis tight
        set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
        title('max-basins (by sum)');
        
        figure; imshow(maxBasins(basins, msN)); axis on; axis tight
        set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
        title('max-basins (by count)');
        
%         plotMainPeakMeanSd(maxBasins(basins, msS), false, ...
%             sprintf('peak-basin means bw=%d', bwidth));
        plotMainPeakMeanSd(maxBasins(basins, msS), true, ...
            sprintf('peak-basin means bw=%d (normalized)', bwidth));
%         saveas(gca, sprintf('008_meanSdPeakByMeanShiftBasin_bw%d.png', bwidth));
        
        rPlotVals = [177 178 191 240];
        for ii=1:1:length(rPlotVals)
            figure; plotyy(1:nLagVals, xCorrMtx(:, rPlotVals(ii)),...
                1:nLagVals, msS(:, rPlotVals(ii))); 
            title(sprintf('cross section rIdx %d', rPlotVals(ii)));
        end
        
%         figure; plotyy(1:nLagVals, xCorrMtx(:, 178),...
%             1:nLagVals, msS(:, 178));
%         figure; plotyy(1:nLagVals, xCorrMtx(:, 191),...
%             1:nLagVals, msS(:, 191));
%         figure; plotyy(1:nLagVals, xCorrMtx(:, 240),...
%             1:nLagVals, msS(:, 240));
%         figure; plotyy(1:nLagVals, basins(:, 240),...
%             1:nLagVals, msS(:, 240));
    end

    function plotMainPeakMeanSd(inMainPeak, normalize, plotTitle)
        xcvPeak = xCorrMtx .* inMainPeak;
        if normalize
%             figure; imagesc(xcvPeak); title(['beforeNorm: ' plotTitle]); axis image; set(gca,'YDir','normal'); impixelinfo
%             mins = zeros(1, nRadVals);
%             maxs = zeros(1, nRadVals);
            for ii=1:1:nRadVals
                curXcvPeak = xcvPeak(:, ii);
                isPos = curXcvPeak > 0;
                if sum(isPos) > 0
%                     mins(ii) = min(nonzeros(curXcvPeak));
%                     maxs(ii) = max(nonzeros(curXcvPeak));
                    xcvPeak(isPos, ii) = xcvPeak(isPos, ii) - min(xcvPeak(isPos, ii));
                    xcvPeak(isPos, ii) = xcvPeak(isPos, ii) ./ max(xcvPeak(isPos, ii));
                end
            end
%             figure; plot(mins); title(sprintf('%s', mat2str([min(mins) max(mins)])));
%             figure; plot(maxs); title(sprintf('%s', mat2str([min(maxs) max(maxs)])));
        end
        idxs = [1:nLagVals]';
        idxs = repmat(idxs, 1, nRadVals);
%         figure; imagesc(xcvPeak); title(['xcvPeak: ' plotTitle]); axis image; set(gca,'YDir','normal')
        means = sum(xcvPeak .* idxs, 1) ./ sum(xcvPeak, 1);
        sdevs = sqrt(sum(xcvPeak .* (idxs - repmat(means, nLagVals, 1)).^2) ./ sum(xcvPeak, 1));
        figure; imshow(inMainPeak); title(['in peak: ' plotTitle]); set(gca,'YDir','normal')
        plotMax(means, sdevs, plotTitle);
    end

    function plotMax(idx, sd, maxName)
        if nargin < 1
            idx = [];
        end
        if nargin < 2
            sd = [];
        end
        if nargin < 3
            maxName = [];
        end
        figure; hold on; 
        imagesc(xCorrMtx); axis image; colormap gray
        line([0 nRadVals], find(lags==0) * [1 1])
        xlabel('rIdx');
        ylabel('lag');
        set(gca,'YDir','normal')
        set(gca, 'XTick', radTickPcts * nRadVals)
        lbls = get(gca, 'YTickLabel');
        set(gca, 'YTickLabel', lags(str2double(mat2cell(lbls, ones(size(lbls, 1), 1), size(lbls, 2)))));
        if ~isempty(idx) && isempty(sd)
            scatter(1:nRadVals, idx, 'r.', 'SizeData', 5^2);
        elseif ~isempty(idx) && ~isempty(sd)
            errorbar(1:nRadVals, idx, sd, 'r.');
        end
        if ~isempty(maxName)
            title(maxName);
        end
    end

    function [msS, msN, basins] = meanShift(bwidth)
        bwidth = round(bwidth);
        msS = zeros(size(xCorrMtx));
        msN = zeros(size(xCorrMtx));
        basins = zeros(size(xCorrMtx));
        nIdx = nLagVals;
%         fromTo = zeros(0, 2);
        for radIdx=1:1:size(msS, 2)
%             fprintf('radIdx = %d\n', radIdx);
            curXc = xCorrMtx(:, radIdx);
            examined = false(size(curXc));
            while(nnz(~examined) > 0)
%                 fprintf('\tnnex = %d\n', nnz(~examined));
                nexmI = find(~examined);
                [tmp, idx] = min(curXc(~examined));
                stIdx = nexmI(idx);
                curIdx = stIdx;
                prevIdx = [];
                wasPrev = false(size(examined));
                while(~wasPrev(curIdx))
                    prevIdx = curIdx;
                    wasPrev(prevIdx) = true;
                    idxs = [max(curIdx-bwidth, 1):1:min(curIdx+bwidth, nIdx)]';
%                     if radIdx == 63
%                         fprintf('\//\//\///n');
%                         curIdx
% %                         idxs
%                     end
                    if sum(curXc(idxs)) == 0
                        curMean = stIdx;
%                         if radIdx == 63
%                             curMean
%                         end
                    else
                        xcScores = curXc(idxs);
                        xcScores = xcScores - min(xcScores);
                        xcScores = xcScores / max(xcScores);
%                         wts = curXc(idxs).^2;
%                         wts = exp(curXc(idxs));
%                         wts = exp(curXc(idxs).^2);
%                         wts = curXc(idxs).^2 .* exp((idxs - curIdx).^2);
                        wts = xcScores.^2 .* exp((idxs - curIdx).^2);
                        wts = xcScores.^2;
%                         [mV, mI] = max(curXc(idxs));
%                         wts = zeros(size(idxs)); wts(mI) = 1;
                        curMean = sum(idxs .* wts) ./ sum(wts);
%                         if radIdx == 63
%                             wts
%                             curMean
%                         end
                    end
                    shiftAmt = (curMean - curIdx);
                    % limit max step size to 1 in order to avoid
                    % oscillations
                    if abs(shiftAmt) >= 0.5
                        curIdx = curIdx + sign(shiftAmt);
                    end
%                     curIdx = round(curMean);
%                     if radIdx == 63
%                         curIdx
%                         fprintf('---\n');
%                     end
%                     curIdx = round(curMean);
%                     if isnan(curIdx)
%                         curXc(idxs)
%                         curXc(idxs).^2
%                         (idxs - curIdx).^2
%                         exp((idxs - curIdx).^2)
%                         idxs
%                         wts
%                         sum(idxs .* wts)
%                         curIdx
%                     end
%                      if radIdx == 63
%                          fromTo = [fromTo; prevIdx curIdx];
%                      end
                end
%                 if radIdx == 63
%                 	fromTo = [fromTo; 0 0];
%                 end
                wasPrev(stIdx) = true;
                if curIdx ~= prevIdx
%                     find(wasPrev)
                    warning('mean-shift cycle detected at radIdx %d, starting lagIdx %d', ...
                        radIdx, stIdx);
                    % make the convergence (peak) point the visited point
                    % with maximum value, rather than an arbitrary point in
                    % the cycle
                    [mV, mI] = max(curXc .* wasPrev);
                    curIdx = mI;
                end
%                 inCurBasin = false(size(curXc)); % TODO: better name
%                 if curIdx <= stIdx
%                     inCurBasin(curIdx:stIdx) = true;
%                 else
%                     inCurBasin(stIdx:curIdx) = true;
%                 end
                msS(curIdx, radIdx) = msS(curIdx, radIdx) + sum(curXc(wasPrev & ~examined));
                msN(curIdx, radIdx) = msN(curIdx, radIdx) + sum(nnz(wasPrev & ~examined));
%                 oldBasinRgn = basins(wasPrev, radIdx);
%                 prevSetBsns = oldBasinRgn > 0;
%                 if nnz(prevSetBsns) > 0
%                     stIdx
%                     radIdx
%                     find(wasPrev)
%                     find(wasPrev & basins(:, radIdx) > 0)
%                     oldBasinRgn(prevSetBsns)
%                     oldBasinRgn
%                 end
%                 if radIdx == 63
%                     fprintf('setting basins tracker\n');
%                     find(wasPrev)
%                     basins(wasPrev, radIdx)
%                     curIdx
%                 end
                basins(wasPrev, radIdx) = curIdx;
%                 newBasinRgn = basins(wasPrev, radIdx);
%                 if nnz(prevSetBsns) > 0
%                     newBasinRgn(prevSetBsns)
%                     newBasinRgn
%                     fromTo
%                 end
%                 assert(nnz(prevSetBsns) == 0 || sum(oldBasinRgn(prevSetBsns) ~= newBasinRgn(prevSetBsns)) == 0);
%                 sum(wasPrev)
                examined = examined | wasPrev;
            end
        end
%         save('fromTo', 'fromTo');
        figure; imagesc(msS); axis image;  
        set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
        title(sprintf('mean-shift sums (bandwidth %d)', bwidth));
        
        figure; imagesc(msN); axis image; 
        set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
        title(sprintf('mean-shift #pts (bandwidth %d)', bwidth));
        
        figure; imagesc(basins); colormap gray; axis image; title('basins');
        set(gca, 'XTick', radTickPcts * nRadVals); set(gca,'YDir','normal')
        
%         [tmp, mxcIdxs] = max(xCorrMtx, [], 1);
%         figure; imshow(basins == repmat(mxcIdxs, nLagVals, 1));

%         [tmp, mxcIdxs] = max(msS, [], 1);
%         figure; imshow(basins == repmat(mxcIdxs, nLagVals, 1)); axis on; axis tight
%         set(gca, 'XTick', radTickPcts * nRadVals)
%         title('points reaching max mode');
        
        
%         xcImg = repmat(xCorrMtx .* (msN == 0), [1 1 3]);
%         xcImg(:, :, 1) = xcImg(:, :, 1) + (msN / max(max(msN)));
%         figure; imshow(xcImg);
%         
%         xcImg = repmat(xCorrMtx .* (msS == 0), [1 1 3]);
%         xcImg(:, :, 1) = xcImg(:, :, 1) + (msS / max(max(msS)));
%         figure; imshow(xcImg);
    end

    function isMaxBasin = maxBasins(basins, peakStrengths)
        isMaxBasin = false(size(basins));
        [maxVals, maxIdxs] = max(peakStrengths, [], 1);
        peakBorders = diff([zeros(1, nRadVals); peakStrengths > 0; zeros(1, nRadVals)], 1, 1);
        for ii=1:1:nRadVals
            peakStarts = find(peakBorders(:, ii) == 1);
            peakEnds = find(peakBorders(:, ii) == -1) - 1;
            for jj=1:1:length(peakStarts)
                if peakStarts(jj) <= maxIdxs(ii) && maxIdxs(ii) <= peakEnds(jj)
                    isMaxBasin(:, ii) = ...
                        (peakStarts(jj) <= basins(:, ii) & ...
                        basins(:, ii) <= peakEnds(jj));
                end
            end
        end
    end
end