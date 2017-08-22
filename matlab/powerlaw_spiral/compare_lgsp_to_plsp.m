function [plspParams, plspBounds, lgspParams, lgspBounds] = ...
    compare_lgsp_to_plsp(img, stgs, ctrR, ctrC)

[lgspParams, lgspBounds, lgspErr] = fitLogSpiral(img, ctrR, ctrC, stgs);
[plspParams, plspBounds, plspErr] = fitPowerlawSpiral(img, ctrR, ctrC, stgs);

display_arc_plot(lgspParams, lgspBounds, @logSpiralXY, img, ctrR, ctrC);
title(sprintf('log spiral, thOff=%2.4f, a=%2.4f, ir=%2.4f\nerror = %2.4f',...
    lgspParams(1), lgspParams(2), lgspParams(3), lgspErr))
display_arc_plot(plspParams, plspBounds, @powerlawSpiralXY, img, ctrR, ctrC);
title(sprintf('powerlaw spiral, thOff=%2.4f, k=%2.4f, ir=%2.4f\nerror = %2.4f',...
    plspParams(1), plspParams(2), plspParams(3), plspErr))

end