function run_all(inDir, outDir)

files = dir(inDir);
% pangs = NaN * ones(length(files), 2);
% stgs = getDefaultSettings();
for ii=1:1:length(files)
    files(ii).name
    if regexpi(files(ii).name, 'NGC\w*\.mat')
        files(ii).name
        load([inDir filesep files(ii).name])
%         ctrR = size(clusMtxList{1}, 1) / 2;
%         ctrC = size(clusMtxList{1}, 2) / 2;
%         [params1, bounds1] = fitLogSpiralsToClusters(clusMtxList{1}, ctrR, ctrC, stgs);
%         lengths1 = calcLgspArcLengths(params1, bounds1);
%         alenwPa1 = (params1(:, 2)' * lengths1) / sum(lengths1);
%         [params2, bounds2] = fitLogSpiralsToClusters(clusMtxList{2}, ctrR, ctrC, stgs);
%         lengths2 = calcLgspArcLengths(params2, bounds2);
%         alenwPa2 = (params2(:, 2)' * lengths2) / sum(lengths2);
%         pangs(ii, 1) = alenwPa1;
%         pangs(ii, 2) = alenwPa2;
        
        name = regexpi(files(ii).name, 'NGC(\w)*', 'match');
        name = name{1};
%         cmpWbands(clusMtxList, sfimgs, [outDir filesep name]);
        cmsz = size(clusMtxList{1});
        cmsz = cmsz(1:2);
        
        img1vis = convertFromFits(sfimgs{1});
        img2vis = convertFromFits(sfimgs{2});
        figure; 
        subplot(1, 2, 1); imshow(img1vis); title('band 1');
        subplot(1, 2, 2); imshow(img2vis); title('band 2');
        saveas(gca, [outDir filesep name '_000_images.png']);
        
        stgs = getDefaultSettings();
        
        fprintf('merging clusters in first band\n');
        clusMtxList{1} = mergeClustersByFit(clusMtxList{1}, cmsz(1)/2, cmsz(2)/2, [], stgs);
        
        [cmm1, cmm2] = matchClusters(clusMtxList{1}, clusMtxList{2}, cmsz(1)/2, cmsz(2)/2, stgs);
        hasBlueMtch = squeeze(sum(sum(cmm2, 1), 2)) > 0;
        cmm1 = cmm1(:, :, hasBlueMtch);
        cmm2 = cmm2(:, :, hasBlueMtch);
        showClustersfromMtxs(cmm1, [512 512]); % temp
        showClustersfromMtxs(cmm2, [512 512]); % temp
        cmm1 = wtClustMtxsByFitsBrt(cmm1, sfimgs{1});
        cmm2 = wtClustMtxsByFitsBrt(cmm2, sfimgs{2});
%         cmm1 = im2double(cmm1 > 0);
%         cmm2 = im2double(cmm2 > 0);
        plotInterBandDists(cmm1, cmm2, cmsz, cmsz(1)/2, cmsz(2)/2, getDefaultSettings(), true, name, outDir);

%         plotInterBandDists(clusMtxList{1}, clusMtxList{2}, cmsz, cmsz(1)/2, cmsz(2)/2, getDefaultSettings(), false, name, outDir);
        
        
%         matchClusters(clusMtxList{1}, clusMtxList{2}, cmsz(1)/2, cmsz(2)/2, stgs, name, outDir);
        close all
    end
end
% pangs = pangs(~isnan(pangs(:, 1)), :);
% pangs = pangs * (180/pi);
% figure; hold on
% scatter(1:size(pangs, 1), abs(pangs(:, 1)), 'r.');
% scatter(1:size(pangs, 1), abs(pangs(:, 2)), 'b.');

end