function errs = getLgspFitErrs(stgs, clusMtx, ctrR, ctrC, lgspParams, lgspBounds, fitTh)

allowArcBeyond2pi = stgs.allowArcBeyond2pi;

if allowArcBeyond2pi
    lgspFxn = @logSpiralFxn2Rev;
else
    lgspFxn = @logSpiralFxn;
end

[rows, cols, brt] = find(clusMtx);
[theta, rho] = cart2pol(cols - ctrC, -(rows - ctrR));

fitTh = fitTh - min(fitTh) + lgspParams(1) + lgspBounds(1);

errs = lgspErr(lgspFxn, lgspParams, fitTh, rho, ones(size(brt)));
errs = errs.^2;

end