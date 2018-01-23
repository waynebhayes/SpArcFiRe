function asgns = calcClusAsgns(img, lsParams, ctrR, ctrC, thBounds)
% experimental code; used by emRefineArcs

maxDist = 5;

nClus = size(lsParams, 1);
ctrX = ctrC;
ctrY = size(img, 1) - ctrR + 1;

asgns = zeros([size(img) nClus]);
for clus = 1:1:nClus
    asgns(:, :, clus) = ...
        polXyMtx(size(img), ctrX, ctrY, @logSpiralXY, ...
        lsParams(clus, :), thBounds(clus, :));
    asgns(:, :, clus) = maxDist - bwdist(asgns(:, :, clus));
end
asgns(asgns < 0) = 0;

asgns = asgns .* (asgns ./ repmat(sum(asgns, 3), [1, 1, nClus]));

% happens when a pixel is more than maxDist away from any cluster arc
% maybe this should be fixed more elegantly?
asgns(isnan(asgns)) = 0; 

end