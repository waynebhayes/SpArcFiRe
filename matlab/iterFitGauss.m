function [muFinal, covarFinal, nonConvFlag, likFinal, gxyParams, bulgeMask] = ...
    iterFitGauss(img, stgs, muFix, gxyParams, plotFlag)
% Iteratively fits a Gaussian to the image pixels (re-weighting according
% to an estimate of the probability that the pixel value actually belongs 
% to the Gaussian).  Intensities are considered weights on data points in 
% Cartesian x-y coordinate space, with the origin at the bottom left.
% INPUTS:
%   img: image for which the Gaussian parameters are to be computed
%   stgs: 
%   muFix: optional; if set, fixes the mean to the given value and only
%       calculates the covariance from the data
%   gxyParams: structure used for collecting information about the output
%   plotFlag: optional; whether to plot some tracing information about the
%       fit (default false)
% OUTPUTS:
%   muFinal: the mean of the fitted Gaussian
%   covarFinal: the covariance matrix of the fitted Gaussian
%   nonConvFlag: true if the fit failed to converge within the allowed
%       number of iterations
%   likFinal:
%   barSc: the score for how likely the galaxy actually has a bar with the
%       reported parameters. 
%   barInfo: the parameters for the bar, reported for ALL galaxies,
%       regardless of whether there is a bar present.  Use the barSc
%       parameter for a measure of whether the bar actually exists.
%   gxyParams:

if nargin < 3
    muFix = [];
end

% fprintf('fixing mu to image center!!!\n');
% muFix = size(img)/2;

if nargin < 5 || isempty(plotFlag)
    plotFlag = false;
end
% plotFlag = true

lookForCtr = stgs.useTwoStageCtrFinding && isempty(muFix);
if lookForCtr
    origPlotFlag = plotFlag;
    plotFlag = false;
    convTol = 0.02;
end

lookForBulge = stgs.lookForBulge && ~lookForCtr;
bulgeFound = false;

% lookForBar = ~lookForCtr && (nargout > 4);

% barDetContourLvl = stgs.barDetContourLvl;

likCutoff = 10^-9;

% starMaskThres = 5
% starBrtThres = inf;

maxRefitsDisk = 25;
if lookForCtr
    maxRefits = 100;
%     maxRefits = 50;
elseif lookForBulge
    maxRefits = 250;
else
    maxRefits = maxRefitsDisk;
end

nonConvFlag = false;

% if lookForBar
% %     minAr = 1;
% %     muMinAr = [];
% %     covarMinAr = [];
% %     likMinAr = [];
% %     maxBarSc = -inf;
% %     maxBarInfo = -inf;
%     eligBarCand = false(1, maxRefits);
%     muTrace = cell(1, maxRefits);
%     covarTrace = cell(1, maxRefits);
% end

% [imgUsm, adjAmt] = unsharpMask(img);
% img(adjAmt > starMaskThres) = 0;
% img(img > starBrtThres) = 0;
% figure; imshow(img);

% img = img .* imextendedmin(img, 0.6);

% initConvCount = 2;
% convCount = initConvCount;

pixWts = (1 / numel(img)) .* img;
% pixWts = calcLikelihoods(size(img), size(img) / 2, 16 * size(img) * eye(2)) .* img;

[mu, covar] = calcWtdMeanCovar(pixWts, muFix);
% mcPrev = [mu covar(1) covar(2) covar(4)]';

% wtTrace = [];
% figure
if plotFlag
    sumWtsTrace = NaN * zeros(1, maxRefits);
    diffTrace = NaN * zeros(1, maxRefits);
    diff2Trace = NaN * zeros(1, maxRefits);
    areaTrace = NaN * zeros(1, maxRefits);
    iterConts = figure;
end
% dbgBarScTrace = NaN * zeros(1, maxRefits);
% dbgMinArIter = -1;
% dbgMaxBarScIter = -1;
if plotFlag
    axisRatioTrace = NaN * zeros(1, maxRefits);
end
converged = false;
for refit = 1:1:maxRefits % change this condition
%     if converged && lookForBar
%         egv = eig(covar);
%         if max(egv) < 1
%             break;
%         end
%     end

    lik = calcLikelihoods(size(img), mu, covar + eps * eye(2));

    titleStr = '';
    
    pixWts = img .* lik;
    
    sumWtsCur = sum(pixWts(:));
    if plotFlag
        sumWtsTrace(refit) = sumWtsCur;
        areaTrace(refit) = sum(lik(:) > 10^-9);
    end
    if plotFlag
        d = eig(covar); 
        curAxisRatio = sqrt(min(abs(d))) / sqrt(max(abs(d)));
        axisRatioTrace(refit) = curAxisRatio;
        muTrace{refit} = mu;
        covarTrace{refit} = covar;
    end
%     if lookForBar && converged
%         [t1, t2, t3, maLenBarCand] = findCovarElpsAxes(covar, barDetContourLvl);
%         
%         if maLenBarCand < (minAxsLenFinal/2)
%             eligBarCand(refit) = true;
% %             [tmpSc, tmpInfo] = ...
% %                 houghFindBar(img, mu, covar, majAxsLenFinal, getDefaultSettings(), barDetContourLvl);
% %             dbgBarScTrace(refit) = tmpSc;
% %             if tmpSc > maxBarSc
% %                 maxBarSc = tmpSc;
% %                 maxBarInfo = tmpInfo;
% %                 dbgMaxBarScIter = refit;
% %             end
%         end
%     end

%     if lookForBar && converged && curAxisRatio < minAr
%         [t1, t2, t3, maLenBarCand] = findCovarElpsAxes(covar, barDetContourLvl);
%         
%         if maLenBarCand < (minAxsLenFinal/2)
%             minAr = curAxisRatio;
% %             muMinAr = mu;
% %             covarMinAr = covar;
% %             likMinAr = lik;
%         end
%     end
    if refit > 1
        sumWtsDiff = sumWtsCur - sumWtsPrev;
%         titleStr = [titleStr sprintf('diff = %2.4f ', sumWtsDiff)];
        diffTrace(refit) = sumWtsDiff;
    end
    if refit > 2
        sumWtsDiff2 = sumWtsCur - 2*sumWtsPrev + sumWtsPrev2;
%         titleStr = [titleStr sprintf('diff2 = %2.4f', sumWtsDiff2)];
        diff2Trace(refit) = sumWtsDiff2;
    end
    
%     if plotFlag
%         set(0,'CurrentFigure', iterConts);
% %         subplot(ceil(maxRefits/5), 5, refit)
%         figure; hold on
%         hold on
%         imagesc(flipud(img)); axis image; colormap gray
%         contour(flipud(lik), 3, '-r', 'LineWidth', 2)
% %         contour(flipud(lik), [10^-8, 10^-9], '-g')
%         contour(flipud(lik), [10^-9], '-g', 'LineWidth', 2)
% %         title(titleStr)
% %         title([sprintf('refit %d ', refit) titleStr])
%         axis off
%         hold off
%         set(gcf, 'color', 'w');
%         export_fig(sprintf('refit%02d.pdf', refit), '-m0.5', '-a1');
%     end
    
    if refit > 3 && ~converged && ~lookForCtr ...
            && convergenceCriteriaMet(sumWtsDiff2Prev, sumWtsDiff2)
        % converged
        muFinal = mu;
        covarFinal = covar;
        likFinal = lik;
%         logLik = sum(log(likFinal(:)));
%         gaussMatch = sum(sum( (lik / sum(lik(:))) .* img ));
%         gaussMatch2 = sum(sum( (lik / sum(lik(:))) .* (img / sum(img(:))) ));
        refitConv = refit;
%         fprintf('refitConv=%2.4f\n', refitConv);
        converged = true;
        
%         [t1, t2, axisRatio, majAxsLenFinal] = findCovarElpsAxes(covarFinal, likCutoff, size(img));
%         minAxsLenFinal = majAxsLenFinal * axisRatio;
        
        if ~plotFlag && ~lookForBulge
            break
        end
    elseif ~converged && ~lookForCtr && refit == maxRefitsDisk
        warning('failed to converge after %d re-fits', maxRefits);
        warn_msg = sprintf(...
            'preprocessing:ellpseFinding:failedToConvergeAfter%dRefits',...
            maxRefits);
        gxyParams.warnings = [gxyParams.warnings warn_msg];
        
        nonConvFlag = true;
        muFinal = mu;
        covarFinal = covar;
        likFinal = lik;
        refitConv = refit;
%         logLik = sum(log(likFinal(:)));
%         gaussMatch = sum(sum( (lik / sum(lik(:))) .* img ));
%         gaussMatch2 = sum(sum( (lik / sum(lik(:))) .* (img / sum(img(:))) ));
    end
    
    if lookForBulge && converged && (refit - refitConv) > 1 && ~bulgeFound &&...
            (sumWtsDiff2 - sumWtsDiff2Prev) >= 0 && sumWtsDiff2 < 0 && sumWtsDiff2 >= -0.000001
            
        covarBulge = covar + eps * eye(2);
        likBulge = lik;
        [bulgeMajAxsVec, bulgeMajAxsAngle, bulgeAxisRatio, bulgeMajAxsLength] = findCovarElpsAxes(...
            covarBulge, likCutoff, size(img));
        if bulgeAxisRatio > 0.75
            gxyParams.badBulgeFitFlag = false;
            gxyParams.bulgeAxisRatio = bulgeAxisRatio;
            gxyParams.bulgeMajAxsLen = bulgeMajAxsLength;
            gxyParams.bulgeMajAxsAngle = bulgeMajAxsAngle;
            bulgeFound = true;
            if ~lookForCtr
                break
            end
        else
            clear covarBulge
            clear likBulge
        end
    end
    
    if refit > 1
        sumWtsPrev2 = sumWtsPrev;
    end
    if refit > 2
        sumWtsDiff2Prev = sumWtsDiff2;
    end
    sumWtsPrev = sumWtsCur;
    
    if lookForCtr
        muPrev = mu;
    end
    
    [mu, covar] = calcWtdMeanCovar(pixWts, muFix);
    
    if lookForCtr
        muChg = sqrt(sum((mu - muPrev).^2));
        if (muChg < convTol) || (refit == maxRefits)
            [muFinal, covarFinal, nonConvFlag, likFinal, gxyParams, bulgeMask] = ...
                iterFitGauss(img, stgs, mu, gxyParams, origPlotFlag);
            return;
        end
    end
end

if lookForBulge 
    if ~bulgeFound
        covarBulge = covar + eps * eye(2);
        likBulge = lik;
        [bulgeMajAxsVec, bulgeMajAxsAngle, bulgeAxisRatio, bulgeMajAxsLength] = findCovarElpsAxes(...
            covarBulge, likCutoff, size(img));
        gxyParams.badBulgeFitFlag = true;
        gxyParams.bulgeAxisRatio = bulgeAxisRatio;
        gxyParams.bulgeMajAxsLen = bulgeMajAxsLength;
        gxyParams.bulgeMajAxsAngle = bulgeMajAxsAngle;
    end
    bulgeRgn = likBulge >= likCutoff;
    bulgeMask = zeros(size(likBulge, 1), size(likBulge, 2));
    bulgeMask(bulgeRgn) = 1;
    diskRgn = (likFinal >= likCutoff) & ~bulgeRgn;
    diskAvgBrt = sum(img(diskRgn)) / sum(diskRgn(:));
    gxyParams.bulgeAvgBrt = sum(img(bulgeRgn)) / sum(bulgeRgn(:));
    gxyParams.bulgeDiskBrtRatio = gxyParams.bulgeAvgBrt / diskAvgBrt;
else
    bulgeMask = zeros(size(img));
end

if ~lookForBulge && ~isfield(gxyParams, 'badBulgeFitFlag')
    gxyParams.badBulgeFitFlag = [];
    gxyParams.bulgeAxisRatio = [];
    gxyParams.bulgeMajAxsLen = [];
    gxyParams.bulgeMajAxsAngle = [];
    gxyParams.bulgeAvgBrt = [];
    gxyParams.bulgeDiskBrtRatio = [];
end

if ~lookForCtr
    gxyParams.numElpsRefits = refitConv;
end

if plotFlag
    figure; hold on
    imagesc(flipud(img)); axis image; colormap gray
%     contour(flipud(likFinal), '-r')
%     contour(flipud(likFinal), [10^-3.5, 10^-4, 10^-8, 10^-9], '-g')
    contour(flipud(likFinal), [10^-9, 10^-12, 10^-15], '-g', 'LineWidth', 2)
%     scatter(muFinal(1), muFinal(2), 'r*');
    if nonConvFlag
        title(sprintf('cutoff iteration (not converged)'))
    else
        title(sprintf('iteration %d', refitConv))
    end
    set(gcf, 'color', 'w');
    axis off
    hold off
%     export_fig(sprintf('final_iter%02d.pdf', refitConv), '-m0.5', '-a1');
%     figure; plot(wtTrace); title('weight trace');
    figure; plot(1:1:maxRefits, sumWtsTrace); xlabel('iter'); ylabel('sumWts');
    figure; plot(1:1:maxRefits, diffTrace); xlabel('iter'); ylabel('diff');
    figure; plot(1:1:maxRefits, diff2Trace); xlabel('iter'); ylabel('diff2');
    figure; plot(1:1:maxRefits, areaTrace); xlabel('iter'); ylabel('area');
    figure; plot(1:1:maxRefits, axisRatioTrace); xlabel('iter'); ylabel('axis ratio');
%     figure; plotyy(1:1:maxRefits, axisRatioTrace, 1:1:maxRefits, eligBarCand); xlabel('iter'); ylabel('bar-candidate eligibility');
    if ~nonConvFlag
        figure; plot(1:1:maxRefits, areaTrace/areaTrace(refitConv)); xlabel('iter'); ylabel('area / (area at convergence)');
    end
    
%     if lookForBulge
%         figure; hold on
%         imagesc(flipud(img)); axis image; colormap gray
%         contour(flipud(likFinal), [10^-9], '-g', 'LineWidth', 1)
%         contour(flipud(likBulge), [10^-9], '-b', 'LineWidth', 1)
%         export_fig(gcf, ['D:\matlab-output-local\temp\' gxyParams.name '_BulgeFit.png']);
%         figure; plot(1:1:maxRefits, diff2Trace); xlabel('iter'); ylabel('diff2');
%         ylim([-0.00001 0])
% %         export_fig(gcf, ['D:\matlab-output-local\temp\' gxyParams.name '_measTrace.png']);
%     end
end

% if lookForBar
%     % look for local minima in the eligible part of the axis ratio trace.
%     % First point considered a local minimum if the next point is greater;
%     % last point consdiered a local mimimum if the previous point is less
%     arChgs = [-inf diff(axisRatioTrace) inf];
%     arLocMins = find((arChgs(1:end-1) < 0 & arChgs(2:end) > 0) & eligBarCand);
%     if numel(arLocMins) == 0
% %         warning(['no axis-ratio local minima in eligible region, using '...
% %             'smallest for bar detection']);
%         tmp = axisRatioTrace; tmp(~eligBarCand) = inf; [v, idx] = min(tmp);
%         arLocMins = idx;
%     end
%     nLocMin = length(arLocMins);
%     lmBarScs = zeros(1, nLocMin);
%     lmBarInfo = cell(1, nLocMin);
%     for ii=1:1:nLocMin
%         lmIdx = arLocMins(ii);
%         [curSc, curInfo] = houghFindBar(img, muTrace{lmIdx}, covarTrace{lmIdx},...
%             majAxsLenFinal, stgs, plotFlag);
%         lmBarScs(ii) = curSc;
%         lmBarInfo{ii} = curInfo;
%     end
%     lmBarScs
%     lmBarInfo
%     [mV, mI] = max(lmBarScs);
%     barSc = lmBarScs(mI);
%     barInfo = lmBarInfo{mI};
% end

% figure; hold on
% imagesc(img); axis image; colormap gray
% scatter(muFinal(1), muFinal(2), 'b.');

% figure; hold on
% imagesc(flipud(img)); axis image; colormap gray
% muTrace = muTrace(cellfun(@(x)(~isempty(x)), muTrace));
% muChgs = NaN * zeros(1, length(muTrace));
% for ii=2:1:length(muTrace)
%     muCur = muTrace{ii};
%     muPrev = muTrace{ii-1};
%     muChgs(ii) = sqrt(sum((muCur - muPrev).^2));
% end
% stIdx = find(muChgs < 0.02, 1, 'first');
% scatter(cellfun(@(x)(x(1)), muTrace), cellfun(@(x)(x(2)), muTrace), 'g.');
% scatter(muFinal(1), muFinal(2), 'r.');
% scatter(cellfun(@(x)(x(1)), muTrace(stIdx)), cellfun(@(x)(x(2)), muTrace(stIdx)), 'b.');
% figure; plot(muChgs);

% figure; hold on
% imgUsmTemp = unsharpMask(img, getDefaultSettings()); imgUsmTemp(imgUsmTemp < 0) = 0; imgUsmTemp(imgUsmTemp > 1) = 1;
% imagesc(flipud(imgUsmTemp)); axis image; colormap gray
% contour(flipud(likMinAr), [10^-6, 10^-9], '-g')

% barSc = maxBarSc;
% barInfo = maxBarInfo;
% fprintf('dbgMaxBarScIter = %d\n', dbgMaxBarScIter);

% figure; plot(1:1:maxRefits, axisRatioTrace); xlabel('iter'); ylabel('axis ratio');
% figure; plot(1:1:maxRefits, dbgBarScTrace); xlabel('iter'); ylabel('barSc');
% figure; plotyy(1:1:maxRefits, axisRatioTrace, 1:1:maxRefits, dbgBarScTrace);
% dbgMinArIter

% save('axisRatioTrace', 'axisRatioTrace');

% if lookForBar
% %     barScE = houghFindBarExperimental(img, muMinAr, covarMinAr, majAxsLenFinal, getDefaultSettings())
%     [barSc, barInfo] = ...
%         houghFindBar(img, muMinAr, covarMinAr, majAxsLenFinal, getDefaultSettings(), barDetContourLvl, plotFlag);    
% end

end

function tf = convergenceCriteriaMet(sumWtsDiff2Prev, sumWtsDiff2)
    tf = (sumWtsDiff2 - sumWtsDiff2Prev) >= 0 && sumWtsDiff2 < 0 && sumWtsDiff2 >= -0.005;
end
