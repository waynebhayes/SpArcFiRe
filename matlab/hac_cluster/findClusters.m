function clusters = findClusters(img, barInfo, stgs, plotFlag)
% Driver function for finding and displaying HAC clusters derived from an
%   image.  Image preprocessing is not performed here.

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

useMex = stgs.useMex;
stopThres = stgs.stopThres;
colorThres = stgs.clusSizeCutoff;

tStart = tic;

fprintf('calculating orientation field...\n'); tic
ofld = genOriField(img, stgs);
% fix for NaNs from zero-valued image pixels
% TODO: fix this in orientation field code
ofld(isnan(ofld)) = 0;
toc; fprintf('...done calculating orientation field\n')

simls = genSimilarityMtx(img, ofld, stgs, stopThres);

clusters = genHacTree(simls, img, barInfo, stgs);


if plotFlag
    showClusters(clusters, size(img), colorThres);
end

fprintf('Time for all clustering steps: \n');
toc(tStart)

end