function arcDist = calcArcDist(params1, bounds1, params2, bounds2, imgSz, ctrR, ctrC, stgs)
% This function is obsolete; arc-distance is no longer used

if ~stgs.useArcDist
    arcDist = -inf;
    return;
end

% extnProp = 0.25;
% len1 = abs(bounds1(2) - bounds1(1)); extnAmt1 = extnProp * len1;
% len2 = abs(bounds2(2) - bounds2(1)); extnAmt2 = extnProp * len2;
extnAmt1 = pi/8; extnAmt2 = pi/8;
exBnds1 = bounds1 + [-extnAmt1 extnAmt1];
exBnds2 = bounds2 + [-extnAmt2 extnAmt2];
ctrX = ctrC;
ctrY = imgSz(1) - ctrR + 1;
arcMtx1 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params1, bounds1);
arcExMtx1 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params1, exBnds1);
arcMtx2 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params2, bounds2);
arcExMtx2 = polXyMtx(imgSz, ctrX, ctrY, @logSpiralXY, params2, exBnds2);

arcMatch12 = bwdist(arcMtx2); 
arcMatch12 = arcMatch12((arcExMtx1 - arcMtx1) > 0);
% arcMatch12 = arcMatch12(arcExMtx1 > 0);
% arcMatch12 = min(arcMatch12);
arcMatch12 = sort(arcMatch12, 'ascend'); arcMatch12 = mean(arcMatch12(1:min(10, end)));

arcMatch21 = bwdist(arcMtx1); 
arcMatch21 = arcMatch21((arcExMtx2 - arcMtx2) > 0);
% arcMatch21 = arcMatch21(arcExMtx2 > 0);
% arcMatch21 = min(arcMatch21);
arcMatch21 = sort(arcMatch21, 'ascend'); arcMatch21 = mean(arcMatch21(1:min(10, end)));

arcDist = max([arcMatch12, arcMatch21]);
if isempty(arcDist)
    arcDist = inf;
end

end