function thruCtr = clusMtxIsThruCtr(clusMtx, ctrR, ctrC, inRad, outRad)

if nargin < 5 || isempty(outRad)
    outRad = inf;
end

[rows, cols] = find(clusMtx);
[theta, rho] = cart2pol(cols - ctrC, -(rows - ctrR));
thruCtr = (mean(rho) < inRad);



% [r, c] = find(clusMtx);
% dists = sqrt(sum(([r c] - repmat([ctrR ctrC], length(r), 1)).^2, 2));
% 
% thruCtr = (min(dists) < inRad) && (max(dists) < outRad);



% ctrPixs = clusMtx(ctrR + [-rad:1:rad], ctrC + [-rad:1:rad]);
% thruCtr =  max(ctrPixs(:)) > 0;

end