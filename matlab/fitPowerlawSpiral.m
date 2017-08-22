function [plspParams, plspBounds, sumSqErr, fitTh, used2Rev] = ...
    fitPowerlawSpiral(img, ctrR, ctrC, stgs)

[rows, cols, brt] = find(img);
[tmp, rho] = cart2pol(cols - ctrC, -(rows - ctrR));

[lgspParams, lgspBounds, lgspErr, used2Rev, tmp1, tmp2, fitTh] = fitLogSpiral(...
    img, ctrR, ctrC, stgs);

lsqOpts = optimset('MaxFunEvals', 10000, ...
    'TolFun', 1e-4, 'TolX', 1e-4, 'Display', 'off');

% if lgspParams(2) > 0
%     fitThStart = min(fitTh);
%     fitTh = -fitTh;
%     fitTh = fitTh - min(fitTh) + fitThStart;
% end

% if lgspParams(1) < min(fitTh)
%     lb = [-inf, -inf, -inf];
%     ub = [min(fitTh), inf, inf];
% elseif lgspParams(1) > max(fitTh)
%     lb = [max(fitTh), -inf, -inf];
%     ub = [inf, inf, inf];
% else
%     error('log-spiral theta-offset parameter within range of theta-values')
% end
% [plspParams, sumSqErr] = lsqnonlin(@plspErr, [lgspParams(1) 1 1], lb, ub, lsqOpts);

[min(fitTh) max(fitTh)]

% thImg = zeros(size(img));
% thImg(img > 0) = fitTh;
% figure; imagesc(thImg); axis image; title('thImg'); impixelinfo

[fitParams1, sumSqErr1] = lsqnonlin(@plspErr, [lgspParams(1) sign(lgspParams(2)) lgspParams(3)],...
    [-inf, -inf, 1], [min(fitTh), inf, inf], lsqOpts);
[fitParams2, sumSqErr2] = lsqnonlin(@plspErr, [lgspParams(1) sign(lgspParams(2)) lgspParams(3)],...
    [max(fitTh), -inf, 1], [inf, inf, inf], lsqOpts);
% fitParams1
% fitParams2
% sumSqErr1
% sumSqErr2
if sumSqErr1 < sumSqErr2
    fitParams = fitParams1;
    sumSqErr = sumSqErr1;
else
    fitParams = fitParams2;
    sumSqErr = sumSqErr2;
end

% [fitParams, sumSqErr] = lsqnonlin(@plspErr, [lgspParams(1) sign(lgspParams(2)) 1],...
%     [], [], lsqOpts);

plspParams = fitParams;
plspBounds = [min(fitTh) max(fitTh)] - plspParams(1);

% display_arc_plot(fitParams1, [min(fitTh) max(fitTh)] - fitParams1(1), @powerlawSpiralXY, img, ctrR, ctrC);
% display_arc_plot(fitParams2, [min(fitTh) max(fitTh)] - fitParams2(1), @powerlawSpiralXY, img, ctrR, ctrC);

function err = plspErr(params)
    thOff = params(1);
    k = params(2);
    ir = params(3);
    ftheta = ir * ((fitTh - thOff).^(-k));
    err = sqrt(brt) .* (rho - ftheta);

%     if ~isreal(err)
%         max(abs(imag(err)))
%         params
%     end
    err = real(err);
end

end