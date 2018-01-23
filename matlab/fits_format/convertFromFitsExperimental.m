function imgs = convertFromFitsExperimental(imgOrPath, blq, wlq, plotFlag)
% Not used in main flow

if nargin < 2 
    blq = [];
end

if nargin < 3
    wlq = [];
end

if nargin < 4 || isempty(plotFlag)
    plotFlag = false;
end

% plotFlag = true;
% wlq = [];

whiteLvlThres = 0.05;
% whiteLvlThres = 0.01;
% whiteLvlThres = 0.1;

% blq = 0.5;
% wlq = 0.999;
imgSubAmt = 0;
asinhBeta = 1;
asinhBeta = 10;
% asinhBeta = 2 for GZ2 (?)

% q = 1000;
% asinhAlpha = 0.1;

% q = 100;
% asinhAlpha = 2;

% q = 10;
% asinhAlpha = 100;

q = 100;
asinhAlpha = 10;

q = 100;
asinhAlpha = 2;

% q = 25;
% asinhAlpha = 2;

iptWasCell = false;
if ischar(imgOrPath)
    inImgs = {im2double(imread(imgOrPath))};
elseif isnumeric(imgOrPath)
    inImgs = {imgOrPath};
elseif iscell(imgOrPath)
    iptWasCell = true;
    inImgs = cell(size(imgOrPath));
    for ii=1:1:length(inImgs)
        if ischar(imgOrPath{ii})
            inImgs{ii} = im2double(imread(imgOrPath{ii}));
        else
            inImgs{ii} = imgOrPath{ii};
        end
    end
else
    error('unrecognized input type for imgOrPath');
end

% % cutoffPct = 0.01;
% blackLevelQuantile = 0.2;
% % blackLevelQuantile = 0.001;
% % blackLevelQuantile = 0.1;
% whiteLevelQuantile = 0.999;

stretchRangeMin = 0;
stretchRangeMax = 100;
% quant = quantile(img(:), [cutoffPct, 1-cutoffPct]);
% quant = quantile(img(:), [.2, 0.995]);

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
    
    blackLevelQuantile = blq;
    whiteLevelQuantile = wlq;
%     imgv = medfilt2(curImg);
%     imgv = nonzeros(imgv(:));
    imgv = nonzeros(curImg(:));
    if isempty(blq) || isempty(wlq)
        qLvls = 0:0.0001:1;
        qVals = quantile(imgv, qLvls);
        qDiff = diff(qVals);
        qRelDiff = diff(qVals) ./ qVals(1:(end-1));
        [mV, mI] = min(qDiff);
        if isempty(blq)
            blackLevelQuantile = qLvls(mI+1);
        end
        wlqIdx = find(qRelDiff(mI:end) > whiteLvlThres, 1, 'first') + (mI-1) + 1;
        if isempty(wlq)
            whiteLevelQuantile = qLvls(wlqIdx);
        end
    end
    
%     curImg = medFilt2(curImg, [5 5]); % TEMP
    imgMin = min(curImg(:));
    imgMax = max(curImg(:));
%     whiteLevelQuantile
%     blackLevelQuantile
    if isempty(whiteLevelQuantile)
        whiteLevelQuantile = 1;
    end
    if isempty(blackLevelQuantile)
        blackLevelQuantile = 0;
    end
    blackLevel = quantile(imgv, blackLevelQuantile);
    % blackLevel = imgMin;
    whiteLevel = quantile(imgv, whiteLevelQuantile);

%     fprintf('overriding white level!');
%     whiteLevel
%     whiteLevel = blackLevel + sinh(q)/(asinhAlpha*q)
    
    curImg(curImg < blackLevel) = blackLevel;
    curImg(curImg > whiteLevel) = whiteLevel;
    
    if plotFlag
%         cmap = zeros(256, 3); cmap(:, 2) = linspace(0, 1, 256);
        figure; imagesc(curImg); axis image; colormap gray; impixelinfo; colorbar
        title(sprintf(['truncated linear conversion from FITS (blackLevelQuantile = %2.4f, whiteLevelQuantile = %2.4f)\n'...
            'min = %2.4f, blackLevel = %2.4f, whiteLevel = %2.4f, max = %2.4f'], ...
            blackLevelQuantile, whiteLevelQuantile, imgMin, blackLevel, whiteLevel, imgMax));
    end

    curImg = curImg - blackLevel;
%     curImg = curImg * ((stretchRangeMax - stretchRangeMin) / max(curImg(:)));
    curImg = curImg / max(curImg(:));
%     curImg = curImg + stretchRangeMin;
    
    imgsRs{ii} = curImg;
    meanImg = meanImg + curImg;
end
meanImg = meanImg / nImgs;

mImgStch = imgStretchFxn(meanImg);
imgs = cell(size(imgsRs));
for ii=1:1:nImgs
    curImg = imgsRs{ii} .* mImgStch ./ meanImg;
%     curImg = inImgs{ii} .* mImgStch ./ meanImg;
    curImg(meanImg == 0) = 0;
    imgs{ii} = curImg;
end

if ~iptWasCell
    imgs = imgs{1};
end

function img = imgStretchFxn(img)
    % img = sqrt(img);
%     img = log(img + 1);
    % img = sqrt(img);
%     img = asinh(img) - asinh(stretchRangeMin);



%     img = asinh(q*(img - imgSubAmt)/asinhBeta)/q;
    img = asinh(q*asinhAlpha*(img - imgSubAmt))/q;
    img = img - min(img(:));
    img = img / max(img(:));
end


% if plotFlag
%     figure; imagesc(img); axis image; colormap gray; impixelinfo; colorbar
%     title(sprintf(['final image\n'...
%         '(max before rescaling = %2.4f, pre-stretch range = [%2.2f, %2.2f])'],...
%         whiteLevel, stretchRangeMin, stretchRangeMax));
% end

end