function imgs = convertFromFitsOld(imgOrPath, plotFlag, convStgs)
% Rescales FITS-image intensities to a standard range, stretching the 
% brightness values for human viewability (as well as accessibility for 
% standard image-processing approaches)
% INPUTS:
%   imgOrPath: the FITS image matrix, or a path to a FITS image
%   plotFlag: whether to display additional information about the process
%       (mainly for debugging)
%   blq: (optional) black-level quantile to use when rescaling
%   wlq: (optional) white-level quantile to use when rescaling

isRealScalar = @(x)(isreal(x) & isscalar(x));

if nargin < 2 || isempty(plotFlag)
    plotFlag = false;
end

% ip = inputParser;
% ip.addRequired('imgOrPath');
% ip.addOptional('plotFlag', false, @(x)(islogical(x) & isscalar(x)));
% ip.addParamValue('blq', [], isRealScalar);
% ip.addParamValue('wlq', [], isRealScalar);
% ip.addParamValue('blackLevel', [], isRealScalar);
% ip.addParamValue('whiteLevel', [], isRealScalar);
% ip.addParamValue('asinhBeta', [], isRealScalar);
% ip.addParamValue('stretchRangeMin', 0, isRealScalar);
% ip.addParamValue('stretchRangeMax', 100, isRealScalar);
% 
% ip.parse(imgOrPath, plotFlag, varargin{:});
% blq = ip.Results.blq;
% wlq = ip.Results.wlq;
% blackLevel = ip.Results.blackLevel;
% whiteLevel = ip.Results.whiteLevel;
% asinhBeta = ip.Results.asinhBeta;
% stretchRangeMin = ip.Results.stretchRangeMin;
% stretchRangeMax = ip.Results.stretchRangeMax;

blq = []; wlq = []; blackLevel = []; whiteLevel = []; asinhBeta = []; stretchRangeMin = 0; stretchRangeMax = 100;
if isfield(convStgs, 'blq')
    blq = convStgs.blq;
end
if isfield(convStgs, 'wlq')
    wlq = convStgs.wlq;
end
if isfield(convStgs, 'blackLevel')
    blackLevel = convStgs.blackLevel;
end
if isfield(convStgs, 'whiteLevel')
   whiteLevel = convStgs.whiteLevel;
end
if isfield(convStgs, 'asinhBeta')
   asinhBeta = convStgs.asinhBeta;
end
if isfield(convStgs, 'stretchRangeMin')
   stretchRangeMin = convStgs.stretchRangeMin;
end
if isfield(convStgs, 'stretchRangeMax')
   stretchRangeMax = convStgs.stretchRangeMax;
end

if xor(isempty(blq), isempty(wlq))
    error('black and while quantiles must be specified together or not at all');
end
if xor(isempty(blackLevel), isempty(whiteLevel))
    error('black and white levels must be specified together or not at all');
end
if ~isempty(blackLevel) && ~isempty(blq)
    error('brightness quantiles and levels cannot both be specified');
end

if ~isempty(blackLevel)
    brtLvlMethod = 'direct';
elseif ~isempty(blq)
    brtLvlMethod = 'quantile';
else
    warning('attemptying to determine black/white levels automatically');
    brtLvlMethod = 'auto';
end

% whiteLvlThres for auto-quantile-determination only
whiteLvlThres = 0.05;
% whiteLvlThres = 0.01;
% whiteLvlThres = 0.1;

imgSubAmt = 0;

if isempty(asinhBeta)
    warning('no value specified for asinhBeta');
    asinhBeta = 50;
    asinhBeta = 10;
    % asinhBeta = 2 for GZ2 (?)
    % asinhBeta = 1 for 
    asinhBeta = 2;
end

% stretchRangeMin = 0;
% stretchRangeMax = 100;

% convert image input to standard form (cell array of matrices)
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
    
    if strcmp(brtLvlMethod, 'quantile')
        blackLevelQuantile = blq;
        whiteLevelQuantile = wlq;
    else
        blackLevelQuantile = [];
        whiteLevelQuantile = [];
    end
    imgv = nonzeros(curImg(:));
    if strcmp(brtLvlMethod, 'auto')
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
    
    imgMin = min(curImg(:));
    imgMax = max(curImg(:));
%     if isempty(whiteLevelQuantile)
%         whiteLevelQuantile = 1;
%     end
%     if isempty(blackLevelQuantile)
%         blackLevelQuantile = 0;
%     end
    if strcmp(brtLvlMethod, 'quantile')
        blackLevel = quantile(imgv, blackLevelQuantile);
        whiteLevel = quantile(imgv, whiteLevelQuantile);
    end

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
    curImg = curImg * ((stretchRangeMax - stretchRangeMin) / max(curImg(:)));
    curImg = curImg + stretchRangeMin;
    
    imgsRs{ii} = curImg;
    meanImg = meanImg + curImg;
end
meanImg = meanImg / nImgs;

mImgStch = imgStretchFxn(meanImg, asinhBeta, imgSubAmt);
imgs = cell(size(imgsRs));
for ii=1:1:nImgs
    curImg = imgsRs{ii} .* mImgStch ./ meanImg;
    curImg(meanImg == 0) = 0;
    imgs{ii} = curImg;
end

if ~iptWasCell
    imgs = imgs{1};
end

function img = imgStretchFxn(img, asinhBeta, imgSubAmt)
    % img = sqrt(img);
%     img = log(img + 1);
    % img = sqrt(img);
%     img = asinh(img) - asinh(stretchRangeMin);
    img = asinh((img - imgSubAmt)/asinhBeta);
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