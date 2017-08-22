function [ir, err] = bestIr(pa)
% minimal-error initial radius for a given pitch angle

lfv = lgspFxn(theta, [thOff pa 1]);
ir = sum(rho .* lfv .* brt) / sum(brt .* (lfv .^ 2));
err = sum(lgspErr(lgspFxn, [thOff pa ir], theta, rho, brt).^2);
    
end