function clusMtxs = asgnsToMtxs(asgnImg, brtImg, clusSizeCutoff)
% Converts cluster assignments (from the C++ implementation) to cluster
%   matrices

asgnIds = unique(asgnImg(:)); 
asgnIds = asgnIds(arrayfun(...
    @(x)((x > 0) & (sum(asgnImg(:) == x) >= clusSizeCutoff)), asgnIds));

nClus = length(asgnIds);
clusMtxs = zeros([size(brtImg) nClus]);
for ii=1:1:nClus
  clusMtxs(:, :, ii) = brtImg .* (asgnImg == asgnIds(ii));
end

end