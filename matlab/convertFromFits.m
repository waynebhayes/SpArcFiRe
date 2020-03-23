function imgs = convertFromFits(imgOrPath, starMaskPath, starMaskVal, qLvl, nReps, blackQlvl, outPath)
% paramEstMask: which pixels to use when determining brightness
% distribution for conversion parameter decisions

error(nargchk(6, 7, nargin))
% if nargin < 4 || isempty(qLvl)
%     qLvl = 0.9;
% end
% if nargin < 5 || isempty(nReps)
%     nReps = 2;
% end
% if nargin < 6 || isempty(blackQlvl)
%     blackQlvl = 0.05;
% end
if nargin < 7
    outPath = [];
end

if strcmpi(starMaskPath, 'NONE')
    starMaskPath = [];
end
if strcmpi(starMaskVal, 'NONE')
    starMaskVal = [];
end
if isempty(starMaskPath) ~= isempty(starMaskVal)
    error('if starMaskPath or starMaskVal specified, the other must be specified')
end
if ischar(starMaskVal)
    fprintf('starMaskVal given as string; interpreting as floating point number\n');
    starMaskValIn = starMaskVal;
    starMaskVal = str2double(starMaskValIn);
    if isnan(starMaskVal)
        error('cannot interpret %s as a number\n', starMaskValIn);
    end
end

paramEstMask = [];
if ~isempty(starMaskPath)
    try
        mask = imread(starMaskPath);
        paramEstMask = (mask == starMaskVal);
    catch ME
        err_str = sprintf('unable to read mask image %s: %s)',... 
            starMaskPath, ME.message);
        fprintf('WARNING: %s\n', err_str);
%         fprintf(errs_file, '%s\n', err_str);
%         bad_read = true;
    end
end
if ischar(qLvl)
    fprintf('qLvl given as string; interpreting as floating point number\n');
    qLvlIn = qLvl;
    qLvl = str2double(qLvl);
    if isnan(qLvl)
        error('cannot interpret %s as a number\n', qLvlIn);
    end
end
if ischar(nReps)
    fprintf('nReps given as string; interpreting as integer\n');
    nRepsIn = nReps;
    nReps = round(str2double(nReps));
    if isnan(nReps)
        error('cannot interpret %s as a number\n', nRepsIn)
    end
end
if ischar(blackQlvl)
    fprintf('blackQlvl given as string; interpreting as floating point number\n');
    blackQlvlIn = blackQlvl;
    blackQlvl = str2double(blackQlvl);
    if isnan(blackQlvl)
        error('cannot interpret %s as a number\n', blackQlvlIn)
    end
end

if qLvl < 0 || qLvl > 1
    error('qLvl must be between 0 and 1 (inclusive)')
end
if blackQlvl < 0 || blackQlvl > 1
    error('blackQlvl must be between 0 and 1 (inclusive)')
end

% convert image input to standard form (cell array of matrices)
iptWasCell = false;
if ischar(imgOrPath)
    curImg = imread(imgOrPath);
%     inImgs = {im2double(handleImgMask(curImg, imgOrPath))};
    inImgs = {im2double(zeroConstantRowsAndColumns(curImg))};
elseif isnumeric(imgOrPath)
%     inImgs = {handleImgMask(imgOrPath)};
    inImgs = {zeroConstantRowsAndColumns(imgOrPath)};
elseif iscell(imgOrPath)
    iptWasCell = true;
    inImgs = cell(size(imgOrPath));
    for ii=1:1:length(inImgs)
        if ischar(imgOrPath{ii})
            curImg = imread(imgOrPath{ii});
%             inImgs{ii} = im2double(handleImgMask(curImg, imgOrPath{ii}));
            inImgs{ii} = im2double(zeroConstantRowsAndColumns(curImg));
        else
%             inImgs{ii} = im2double(handleImgMask(imgOrPath{ii}));
            inImgs{ii} = im2double(zeroConstantRowsAndColumns(imgOrPath{ii}));
        end
    end
else
    error('unrecognized input type for imgOrPath');
end

nImgs = length(inImgs);
imgSz = size(inImgs{1});
if length(imgSz) ~= 2
    error('images must be 2D matrices');
end

meanImg = zeros(imgSz);
imgsRs = cell(size(inImgs));
for ii=1:1:nImgs
    curImg = inImgs{ii};
    if sum(size(curImg) == imgSz) ~= 2
        error('input images have different sizes');
    end
%     figure; imagesc(curImg); axis image; impixelinfo; title(sprintf('curImg %d', ii));
%     curImg = curImg - min(curImg(:));
%     figure; imagesc(curImg); axis image; impixelinfo
    imgsRs{ii} = curImg;
    meanImg = meanImg + curImg;
%     figure; imshow(meanImg); impixelinfo; title(sprintf('meanImg after %d', ii));
end
meanImg = meanImg / nImgs;
% figure; imshow(meanImg); impixelinfo; title('meanImg');
% figure; imagesc(meanImg); axis image; colormap gray; impixelinfo; title('meanImg')
% figure; imshow(isnan(meanImg)); title('isnan meanImg');
mImgStch = imgStretchFxn(meanImg);
% figure; imshow(mImgStch); impixelinfo; title('mImgStch');
% figure; imagesc(mImgStch); axis image; colormap gray; impixelinfo; title('mImgStch')
% figure; imshow(isnan(mImgStch)); title('isnan mImgStch');
imgs = cell(size(imgsRs));
for ii=1:1:nImgs
%     figure; subplot(1, 3, 1); hist(curImg(:));
    curImg = imgsRs{ii} .* mImgStch ./ meanImg;
%     subplot(1, 3, 2); hist(curImg(:));
    curImg(meanImg == 0) = 0;
%     curImg = curImg - min(curImg(:));
%     curImg = curImg / max(curImg(:));
%     curImg = handleImgMask(curImg);
%     subplot(1, 3, 3); hist(curImg(:));
    imgs{ii} = curImg;
    assert(all(~isnan(curImg(:))), 'all image pixels should be non-NaN')
end

if ~iptWasCell
    imgs = imgs{1};
end

if ~isempty(outPath)
    if iscell(imgs)
        error('outPath can only be used with one input image');
    else
        imwrite(imgs, outPath);
    end
end

    function img = imgStretchFxn(img)
%         img = img - min(img(:));
        img(img < 0) = 0;
        is_nz = (img(:) ~= 0);
        for i=1:1:nReps
            blackLevel = quantile(img(is_nz), blackQlvl);
            fprintf('blackLevel = %2.4f\n', blackLevel);
            img = max(img - blackLevel, 0);
            incl = is_nz;
            if ~isempty(paramEstMask)
                fprintf('asinhBeta (if no star masking were used) = %2.4f\n',...
                    quantile(img(incl), qLvl));
                incl = incl & paramEstMask(:);
            end
            if isempty(incl)
                error('no nonzero pixels within region used to estimate asinhBeta\n');
            end
            asinhBeta = quantile(img(incl), qLvl);
            fprintf('asinhBeta = %2.4f\n', asinhBeta);
            img = asinh(img/asinhBeta)*asinhBeta;
%             img = img - min(img(:));
            img = img / max(img(incl));
            img(img > 1) = 1;
%             img = handleImgMask(img);
        end
    end

end

function [img, bad_rows, bad_cols] = zeroConstantRowsAndColumns(img)

    bad_rows = all(img == repmat(img(:, 1), [1 size(img, 2)]), 2);
    bad_cols = all(img == repmat(img(1, :), [size(img, 1) 1]), 1);
    img(bad_rows, :) = 0;
    img(:, bad_cols) = 0;
    
end

% NOTE: no longer used
% try to deal with cases where the image has been masked or padded with
% (possibly nonzero) values
function img = handleImgMask(img, imgName)
    if nargin < 2 || isempty(imgName)
        imgName = '<unknown image name>';
    end
    assert(all(~isnan(img(:))), 'all image pixels should be non-NaN')
    img = double(img);
%     figure; imagesc(img); axis image; impixelinfo
    modeFreqThres = 0.05;
    imgMin = min(img(:));
    imgMode = mode(img(:));
    subtAmt = imgMode;
    if abs(imgMin - imgMode) > 1
%         fprintf('warning: %s: imgMin different than imgMode\n', imgName);
        modeFreq = sum(abs(img - imgMode) < 10^-4);
        if modeFreq < modeFreqThres
            fprintf('%s: mode frequency proportion (%2.4f) less than %2.4f, using min',...
                imgName, modeFreq, modeFreqThres);
            subtAmt = imgMin;
        end
    end
    img = img - subtAmt;
    img = img / max(img(:));
    if any(isnan(img(:)))
        fprintf(['WARNING: NaN pixels found; zeroing them. '...
            'This can happen if all input pixel values are the same within a band.\n'])
        img(isnan(img)) = 0;
    end
%     figure; imagesc(img); axis image; impixelinfo
end
