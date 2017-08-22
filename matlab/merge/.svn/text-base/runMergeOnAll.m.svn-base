function runMergeOnAll(allClusMtxs, ctrR, ctrC, imgSz, allowArcBeyond2pi, outputDir, plotFlag)

if nargin < 6 || isempty(outputDir)
    outputDir = [];
elseif outputDir(end) ~= filesep
    outputDir = [outputDir filesep];
end

if nargin < 7 || isempty(plotFlag)
    plotFlag = false;
end

for ii=1:1:length(allClusMtxs)
    close all
    clusMtxs = allClusMtxs{ii};
    clusMtxsMerged = mergeClusters(clusMtxs, ctrR, ctrC, allowArcBeyond2pi, plotFlag);
    if ~isempty(outputDir)
        img = showClustersFromMtxs(clusMtxs, imgSz);
        imwrite(img, [outputDir sprintf('%03d', ii) '_A_pre-merge.png']);
        img = showClustersFromMtxs(clusMtxsMerged, imgSz);
        imwrite(img, [outputDir sprintf('%03d', ii) '_B_post-merge.png']);
    end
end

end