function batchComparePowerlawSpiral(in_dir, out_dir)

fnames = dir(in_dir);
fnames = {fnames.name};
valid = cellfun(@(x)(~isempty(strfind(x, 'png'))), fnames);
fnames = fnames(valid);

if ~ischar(out_dir) || ~isdir(out_dir)
    error([out_dir ' is not a directory'])
end

if out_dir(end) ~= filesep
    out_dir = [out_dir filesep];
end

stgs = getDefaultSettings();

allowArcBeyond2pi = stgs.allowArcBeyond2pi;

gxyOutputCtrl = struct('writeImages', true, ...
        'displayFigures', false, 'writeTxt', false);

for idx = 1:1:length(fnames)    
    gxyName = regexp(fnames{idx}, '\.', 'split');
    gxyName = gxyName{1};
    
    curImg = imread([in_dir filesep fnames{idx}]);
    curImg = im2double(curImg);
    if length(size(curImg)) > 2
        curImg = rgb2gray(curImg);
    end
    
    [lgspParams, lgspBounds, lgspErr, used2rev, hasBadBounds, barInfo, ...
            clusMtxs, gxyParams, imgAutoCrop, barInds] = ...
        findClusterArcs(curImg, stgs, gxyName, gxyOutputCtrl, out_dir);
    
    ctrR = 128; ctrC = 128;
    
    hasBar = false(1, size(clusMtxs, 3));
    hasBar(barInds) = true;
    lgspParams = lgspParams(~hasBar, :);
    lgspBounds = lgspBounds(~hasBar, :);
    lgspErr = lgspErr(~hasBar);
    clusMtxs = clusMtxs(:, :, ~hasBar);
    
    lgspAlens = calcLgspArcLengths(lgspParams, lgspBounds);
    [lgspAlens, alenOrder] = sort(lgspAlens, 'descend');
    lgspParams = lgspParams(alenOrder, :);
    lgspBounds = lgspBounds(alenOrder, :);
    lgspErr = lgspErr(alenOrder);
    clusMtxs = clusMtxs(:, :, alenOrder);
    
    nClus = size(clusMtxs, 3);
    curPlspParams = zeros(nClus, 3);
    curPlspBounds = zeros(nClus, 2);
    curPlspErr = zeros(nClus, 1);
    for ii=1:1:nClus
        [curPlspParams(ii, :), curPlspBounds(ii, :), curPlspErr(ii)] = ...
            fitPowerlawSpiral(clusMtxs(:, :, ii), ctrR, ctrC, stgs);
    end
    
%     [tmp1, tmp2, imgc] = preprocessImage(curImg, stgs);
    showClustersFromMtxs(clusMtxs, [256 256]);
    saveas(gcf, [out_dir gxyName '__clusMxs_alenOrder.png'])
    display_arc_plot(lgspParams, lgspBounds, @logSpiralXY, imgAutoCrop, ctrR, ctrC);
    saveas(gcf, [out_dir gxyName '__lgspArcs.png'])
    display_arc_plot(curPlspParams, curPlspBounds, @powerlawSpiralXY, imgAutoCrop, ctrR, ctrC);
    saveas(gcf, [out_dir gxyName '__plspArcs.png'])
    
    outfile = fopen([out_dir gxyName '_param_cmp.txt'], 'wt');
    fprintf(outfile, 'arcNum,lgspAlen,lgspThOff,lgspPa,lgspIr,lgspErr,plspThOff,plspK,plspIr,plspErr\n');
    for ii=1:1:length(lgspAlens)
        fprintf(outfile, ['%d,', repmat('%2.4f,', 1, 9), '\n'], ...
            [ii, lgspAlens(ii), lgspParams(ii, :), lgspErr(ii), curPlspParams(ii, :), curPlspErr(ii)]);
    end
    fclose(outfile);
    close all
end

end