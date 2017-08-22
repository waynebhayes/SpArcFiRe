function ssImg = emStarSuppress(img)
% experimental star-resistance code; not currently in main flow, but an
% improvement on this could be used in the future

maxGauss = 8;

[meanPrev, covPrev, hasSameCC, ssImgPrev] = emFitMultiGauss(img, 1);
nGauss = 2;
hasSameCC = false;
while nGauss < maxGauss && ~hasSameCC
    [meanCur, covCur, hasSameCC, ssImgCur] = emFitMultiGauss(img, nGauss);
    if hasSameCC
        meanCur = meanPrev;
        covCur = covPrev;
        ssImgCur = ssImgPrev;
    else
        meanPrev = meanCur;
        covPrev = covCur;
        ssImgPrev = ssImgCur;
    end
    nGauss = nGauss + 1;
end
ssImg = ssImgCur;

end