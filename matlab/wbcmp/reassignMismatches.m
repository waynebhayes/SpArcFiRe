function [cm1r, cm2r] = reassignMismatches(clusMtxs1, matchedCM2, ctrR, ctrC, stgs)
% Changes pixel-to-cluster assignments to remove mismatched pixels.
% Mismatched pixels are in a cluster in both sets of cluster matrices, but
% the corresponding clusters are not matched with each other. 
% INPUTS:
%   clusMtxs1: first set of clusters
%   matchedCM2: second set of clusters, matched so that the cluster
%       clusMtx1(:, :, ii) matches the cluster matchedCM2(:, :, ii)

if size(clusMtxs1, 3) ~= size(matchedCM2, 3)
    error('matched cluster sets do not have the same number of clusters');
end

nClus = size(clusMtxs1, 3);
clusMemb2 = (matchedCM2 > 0);

clusSizes1 = squeeze(sum(sum(clusMtxs1 > 0, 1), 2));
params1 = NaN*ones(nClus, 3); bounds1 = NaN*ones(nClus, 2); errs1 = NaN*ones(nClus, 1);
nzClus1 = clusSizes1 > 0;
[p1, b1, e1] = fitLogSpiralsToClusters(clusMtxs1(:, :, nzClus1), ctrR, ctrC, stgs);
params1(nzClus1, :) = p1; bounds1(nzClus1, :) = b1; errs1(nzClus1) = e1;
mse1 = errs1 ./ squeeze(sum(sum(clusMtxs1, 1), 2));
clusSizes2 = squeeze(sum(sum(matchedCM2 > 0, 1), 2));
params2 = NaN*ones(nClus, 3); bounds2 = NaN*ones(nClus, 2); errs2 = NaN*ones(nClus, 1);
nzClus2 = clusSizes2 > 0;
[p2, b2, e2] = fitLogSpiralsToClusters(matchedCM2(:, :, nzClus2), ctrR, ctrC, stgs);
params2(nzClus2, :) = p2; bounds2(nzClus2, :) = b2; errs2(nzClus2) = e2;
mse2 = errs2 ./ squeeze(sum(sum(matchedCM2, 1), 2));

for ii=1:1:nClus
    curClus = clusMtxs1(:, :, ii);
    isMism = repmat((curClus > 0), [1 1 nClus]) & clusMemb2;
    for jj=[1:(ii-1) (ii+1):nClus]
        nMism = nnz(isMism(:, :, jj));
        if nMism > 0
            assert(nnz(clusMtxs1(:, :, ii)) > 0)
            assert(nnz(matchedCM2(:, :, jj)) > 0)
            % determine whether it's better to do the reassignment of 
            % mismatched pixels in the first cluster set, or the second
            
            isCurMism = isMism(:, :, jj);
            
            mvInS1Clus = clusMtxs1(:, :, jj) + (clusMtxs1(:, :, ii) .* isCurMism);
            [s1params, s1bounds, s1err] = fitLogSpiral(mvInS1Clus, ctrR, ctrC, stgs);
            mseMv1 = s1err / sum(mvInS1Clus(:));
            errRatioMv1 = mseMv1 / mse1(jj);
            
            mvInS2Clus = matchedCM2(:, :, ii) + (matchedCM2(:, :, jj) .* isCurMism);
            [s2params, s2bounds, s2err] = fitLogSpiral(mvInS2Clus, ctrR, ctrC, stgs);
            mseMv2 = s2err / sum(mvInS2Clus(:));
            errRatioMv2 = mseMv2 / mse2(ii);
            
            assert(~isnan(errRatioMv1) || ~isnan(errRatioMv2));
            if isnan(errRatioMv1)
                errRatioMv1 = inf;
            end
            if isnan(errRatioMv2)
                errRatioMv2 = inf;
            end
            
            % TODO: deal with special case where all the cluster's pixels
            % are reassigned
            if errRatioMv1 < errRatioMv2
                clusMtxs1(:, :, jj) = mvInS1Clus;
                clusMtxs1(:, :, ii) = clusMtxs1(:, :, ii) .* ~isCurMism;
                
                clusSizes1(jj) = clusSizes1(jj) + nMism;
                clusSizes1(ii) = clusSizes1(ii) - nMism;
                params1(jj, :) = s1params; bounds1(jj, :) = s1bounds; errs1(jj) = s1err;
                [newParams, newBounds, newErr] = fitLogSpiral(clusMtxs1(:, :, ii), ctrR, ctrC, stgs);
                params1(ii, :) = newParams; bounds1(ii, :) = newBounds; errs1(ii) = newErr;
                mse1(ii) = errs1(ii) / sum(sum(clusMtxs1(:, :, ii)));
            else
                matchedCM2(:, :, ii) = mvInS2Clus;
                matchedCM2(:, :, jj) = matchedCM2(:, :, jj) .* ~isCurMism;
                
                clusSizes2(ii) = clusSizes2(ii) + nMism;
                clusSizes2(jj) = clusSizes2(jj) - nMism;
                params2(ii, :) = s2params; bounds2(ii, :) = s2bounds; errs2(ii) = s2err;
                [newParams, newBounds, newErr] = fitLogSpiral(matchedCM2(:, :, jj), ctrR, ctrC, stgs);
                params2(jj, :) = newParams; bounds2(jj, :) = newBounds; errs2(jj) = newErr;
                mse2(jj) = errs2(jj) / sum(sum(matchedCM2(:, :, jj)));
            end
        end
    end
end

cm1r = clusMtxs1;
cm2r = matchedCM2;

end % reassignMismatches