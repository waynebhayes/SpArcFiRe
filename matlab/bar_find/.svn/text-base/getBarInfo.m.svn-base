function barInfo = getBarInfo(img, stgs, ctrR, ctrC, plotFlag)
% for testing bar detection on its own; not used in main flow

[imgp, tmp1, tmp2, gxyParams, fitParams] = ...
    preprocessImage(img, stgs);

fprintf('calculating orientation field...\n'); tic
ofld = genOriField(imgp, stgs);
% fix for NaNs from zero-valued image pixels
% TODO: fix this in orientation field code (fixed already?)
ofld(isnan(ofld)) = 0;
toc; fprintf('...done calculating orientation field\n')

% ctrR = size(imgp, 1) / 2; ctrC = size(imgp, 2) / 2;
% if stgs.useSubpixelCtr
%     ctrR = ctrR + 0.5;
%     ctrC = ctrC + 0.5;
% end
[candRad, candTh, candScore] = findBarCandRgn(ofld, ctrR, ctrC, plotFlag);
[barScore, barInfo] = findBarScore(...
    img, fitParams, ctrR, ctrC, candRad, candTh, candScore, stgs, plotFlag);

barInfo

end