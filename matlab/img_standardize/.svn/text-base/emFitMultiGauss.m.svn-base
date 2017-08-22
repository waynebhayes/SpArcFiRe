function [bestMean, bestCovar, ptInSameCC, starSuppressImg] = emFitMultiGauss(img, nGauss, plotFlag)
% experimental star-resistance code; not currently in main flow, but an
% improvement on this could be used in the future

if nargin < 3 || isempty(plotFlag)
    plotFlag = false;
end

maxIter = 25;
convTol = 10^-6;

brtCCthres = 0.2;
brtCCs = bwconncomp(img > brtCCthres);

inclNoiseCmpt = true;

if inclNoiseCmpt
    nCmpt = nGauss+1;
else
    nCmpt = nGauss;
end

imgSz = size(img);
nPts = numel(img);

startVar = eye(2) * imgSz(1) * 2;

% noiseDens = ones(size(img)) / numel(img);

% minMixCoeff = 0.05; minMixCoeff = 0;
% maxNoiseCoeff = 0.01;

% means = [noiseMean; (imgSz / 2) + (imgSz / 4); (imgSz / 2) - (imgSz / 4)]; 
% means = means(2:end, :);
means = getStartMeans();
covars = repmat(startVar, [1 1 nGauss]);
if inclNoiseCmpt
    noiseMean = mean(img(:));
    noiseSd = std(img(:));
else
    noiseMean = [];
    noiseSd = [];
end
mixCoeffs = ones(1, nGauss) / nGauss;
dVals = getAllDens(means, covars, mixCoeffs, noiseMean, noiseSd);
sum(sum(log(sum(dVals(:, :, 1:end-1), 3))))
loglik = sum(sum(log(sum(dVals, 3))));
if plotFlag
    showDens(dVals); title('start');
end

for iter=1:1:maxIter
    mprobs = getMembProbs(dVals);  % E-step
    [means, covars, mixCoeffs, noiseMean, noiseSd] = ...
        getAllGaussCmpts(means, covars, mprobs); % M-step
    dVals= getAllDens(means, covars, mixCoeffs, noiseMean, noiseSd);
    loglikPrev = loglik;
    loglik = sum(sum(log(sum(dVals, 3))));
    fprintf('loglikelihood change: %2.4f\n', loglik - loglikPrev);
    if abs(loglik - loglikPrev) < convTol
        fprintf('converged on iteration %d\n', iter);
        break;
    end
%     if mod(iter, 1) == 0
%         showMprobs(mprobs);
%         showDens(dVals); title(sprintf('iteration %d', iter));
%     end
end

distsToCtr = sum( ((means - repmat(imgSz/2, nGauss, 1)).^2), 2);
[minDist, minIdx] = min(distsToCtr);
bestMean = means(minIdx, :);
bestCovar = covars(:, :, minIdx);

if plotFlag
    figure; imshow(img); hold on; 
    contour(dVals(:, :, minIdx) / mixCoeffs(minIdx), [10^-4, 10^-5, 10^-6], 'r-');
    % figure; imshow(img); hold on; contour(calcLikelihoods(imgSz, bestMean, bestCovar), [10^-4, 10^-5, 10^-6], 'r-');
end

starSuppressImg = img .* (mprobs(:, :, minIdx) + mprobs(:, :, end));
if plotFlag
    figure; imshow(starSuppressImg); colormap gray;
    showMprobs(mprobs); 
    showDens(dVals); title('final');
end

nParam = nGauss * 6 - 1;
if inclNoiseCmpt
    nParam = nParam + 3; % mean, sd, and mixing coefficient
end

bic = -2 * loglik + nParam * log(nPts);

fprintf('lik = %2.8e\n', exp(loglik));
fprintf('loglik = %2.4f\n', loglik);

fprintf('BIC = %2.4f\n', bic);

indsOfMeans = sub2ind(imgSz, imgSz(2) - round(means(:, 2)), round(means(:, 1)));
ccMembs = cellfun(@(x)(ismember(indsOfMeans, x)), brtCCs.PixelIdxList, 'UniformOutput', false);
if max(sum(cell2mat(ccMembs), 1)) > 1 % more than one point in same CC
    ptInSameCC = true;
else
    ptInSameCC = false;
end

% ----- Helper Functions -----

function means = getStartMeans()
    if brtCCs.NumObjects >= nGauss
        nPix = cellfun(@numel,brtCCs.PixelIdxList);
        [temp, srtIdxs] = sort(nPix, 'descend');
        pixRC = cellfun(@ccCtrs, brtCCs.PixelIdxList(srtIdxs(1:nGauss)), 'UniformOutput', false);
        pixRC = cell2mat(pixRC');
        pixXY = [pixRC(:, 2) imgSz(2) - pixRC(:, 1)];
        means = pixXY;
%         means = cell2mat(pixRC') + repmat(imgSz/2, nGauss, 1)
%         means(:, 1) = imgSz(2) - means(:, 1);
%         error('temp');
    else
        warning('more Gaussians than bright regions, arranging initial means in circle');
        r = ones(nGauss, 1) * (imgSz(1) / 4);
        th = (2*pi/nGauss) * [0:1:nGauss-1];
        [x, y] = pol2cart(th', r);
        means = [x y] + repmat(imgSz/2, nGauss, 1);
    end
end

function rc = ccCtrs(idxList)
    [rows, cols] = ind2sub(imgSz, idxList);
    rc = mean([rows, cols], 1);
end

function dVals = getAllDens(means, covars, mixCoeffs, noiseMean, noiseSd)
    dVals = zeros([imgSz nCmpt]);
    assert(sum(isnan(means(:))) == 0);
    assert(sum(isnan(covars(:))) == 0);
    for ii=1:1:nGauss
        dVals(:, :, ii) = mixCoeffs(ii) * ...
            calcLikelihoods(imgSz, means(ii, :), covars(:, :, ii));
    end
    if inclNoiseCmpt
        dVals(:, :, nGauss+1) = normpdf(img, noiseMean, noiseSd);
    end
end

function [means, covars, mixCoeffs, noiseMean, noiseSd] = getAllGaussCmpts(means, covars, mprobs)
%     fprintf('VVVVV\n');
    mixCoeffs = squeeze(sum(sum(mprobs, 1), 2));
    mixCoeffs = mixCoeffs / sum(mixCoeffs);
%     mixCoeffs = max(mixCoeffs, minMixCoeff);
%     mixCoeffs = mixCoeffs / sum(mixCoeffs);
%     if inclNoiseCmpt
% %         mixCoeffs(end) = min(mixCoeffs(end), sum(mixCoeffs(1:end-1)) * maxNoiseCoeff);
% %         mixCoeffs = mixCoeffs / sum(mixCoeffs);
%     end
%     fprintf('^^^^^\n');
    for ii=1:1:nGauss
        [newMean, newCov] = calcWtdMeanCovar(img .* mprobs(:, :, ii));
        means(ii, :) = newMean;
        covars(:, :, ii) = newCov;
    end
    if inclNoiseCmpt
        wts = mprobs(:, :, end);
        sumWts = sum(wts(:));
        noiseMean = sum(sum(img .* wts)) / sumWts;
        noiseSd = sqrt((nPts/(nPts-1)) * (sum(sum(wts .* (img - noiseMean).^2)) / sumWts));
    else
        noiseMean = [];
        noiseSd = [];
    end
end

function mprobs = getMembProbs(dVals)
    mprobs = dVals ./ repmat(sum(dVals, 3), [1 1 nCmpt]);
end

function showDens(dVals)
    figure; imshow(img); axis image; hold on
%     imshow(img); axis image
%     contour(sum(dVals, 3), 'b');
    for ii=1:1:nGauss
        contour(dVals(:, :, ii), 'c');
    end
    hold off
end

function showMprobs(mprobs)
    figure;
    subplot(1, nCmpt+1, 1); imagesc(img); axis image
    for ii=1:1:nCmpt
        subplot(1, nCmpt+1, ii+1); imagesc(mprobs(:, :, ii) .* img); axis image
    end
    figure; imagesc(mprobs(:, :, end)); axis image; title('noise probabilities');
end

end