function [meanRatioScores, anovaFScores, gxyParams] = calcClusScores(img, clusMtxs, gxyParams, outsideClusDist)

defaultDilAmt = 5;

if nargin < 4
%     outsideClusDist = 5;
    outsideClusDist = [];
end
%     clusSizes = squeeze(sum(sum(clusMtxs > 0, 1), 2));
%     clusWidths = clusSizes ./ arcLengths;
%     outsideDists = ceil(clusWidths / 2)
% elseif ~isempty(arcLengths)
%     error('outsideClusDist and arcLengths should not both be specified');
% else
%     outsideDists = outsideClusDist * ones(1, size(clusMtxs, 3));
% end

nClus = size(clusMtxs, 3);
allClus = sum(clusMtxs, 3);

meanRatioScores = zeros(nClus, 1);
anovaFScores = zeros(nClus, 1);
% dilAmt = outsideClusDist;
for cIdx = 1:1:nClus
%     dilAmt = outsideDists(cIdx);
    curClus = clusMtxs(:, :, cIdx) > 0;
    if sum(curClus(:)) == 0
        meanRatioScores(cIdx) = -inf;
        anovaFScores(cIdx) = -inf;
        continue;
    end
    if isempty(outsideClusDist)
        inClus = curClus > 0;
        nBorder = sum(sum(imdilate(inClus, strel('square', 3)) & ~allClus));
        if nBorder == 0
            newWarn = sprintf(...
                'getGxyParams:noBorderPixelsOutsideACluster:arc%dof%d',...
                cIdx, nClus);
            gxyParams.warnings = [gxyParams.warnings newWarn];
            meanRatioScores(cIdx) = NaN;
            anovaFScores(cIdx) = NaN;
            return;
        end
        dilAmt = ceil(sum(inClus(:))/nBorder);
    else
        dilAmt = defaultDilAmt;
    end
    extGrad = imdilate(curClus, strel('square', 2*dilAmt+1)) & ~allClus; % - curClus for true external gradient
    extGrad = extGrad > 0;
    insideItsy = img(curClus);
    outsideItsy = img(extGrad);
    
    meanRatioScores(cIdx) = mean(insideItsy) / mean(outsideItsy);
    
    vals = [insideItsy; outsideItsy];
    group = [zeros(size(insideItsy)); ones(size(outsideItsy))];
    [p,table] = anova1(vals, group, 'off');
    anovaFScores(cIdx) = table{2, 5};
end

end