function batchConvertFromFits(inDirPath, nameRegexp, inSfxs, outSfx, qLvl, nReps, blackQlvl, maskDirPath, maskSuffix, maskVal, outDirPath)
% INPUTS:
%   inDirPath: directory containing the fits files to convert
%   nameRegexp: regular expression describing the name of each object.
%       input files with the same match to this regular expression will be
%       combined into a single composite image.
%   inSfxs: cell array of suffixes for each waveband. For each set of 
%       images that has the same match to nameRegexp, there should be a 
%       1:1 correspondence between the remainder of the file name and an
%       entry in this suffix list.
%   outSfx: suffix for output images (the earlier part of the file name is
%       the match to nameRegexp), including the file extension
%   qLvl: quantile level at which the asinh beta parameter will be set
%   nReps: number of times to repeat the image stretching
%   maskDirPath: directory containing the star masks (empty if star masks
%       should not be used or if the directory is the same as the input 
%       image directory)
%   maskSuffix: suffix for the star masks (empty if star masks should not
%       be used)
%   maskVal: the value in the mask image indicating the pixels from which
%       the brightness distribution should be determined (empty if star
%       masks should not be used
%   outDirPath: directory where the output files should be placed

error(nargchk(11, 11, nargin))

% if nargin < 7
%     maskDirPath = [];
% end
% 
% if nargin < 8
%     maskSuffix = [];
% end
% 
% if nargin < 9
%     maskVal = [];
% end

% hack so that an empty value can be passed when used on the command line
if strcmpi(maskDirPath, 'NONE') || isempty(maskDirPath)
    maskDirPath = [];
end
if strcmpi(maskSuffix, 'NONE') || isempty(maskSuffix)
    maskSuffix = [];
end
if strcmpi(maskVal, 'NONE') || isempty(maskVal)
    maskVal = [];
end

if isempty(maskSuffix) ~= isempty(maskVal)
    error('if maskSuffix or maskVal specified, the other must be specified')
end

if isempty(maskDirPath) && ~isempty(maskSuffix)
    maskDirPath = [inDirPath filesep];
end

if inDirPath(end) ~= filesep
    inDirPath = [inDirPath filesep];
end
if outDirPath(end) ~= filesep
    outDirPath = [outDirPath filesep];
end
    
logPath = [outDirPath 'batchConvertFromFits.log'];
fprintf('output will be written to %s\n', logPath)
errFilePath = [outDirPath 'batchConvertFromFits.err'];
fprintf('galaxy-level errors will also be written to %s\n', errFilePath);
errs_file = fopen(errFilePath, 'wt');
if errs_file < 0
    error('cannot write to %s\n', errFilePath);
end
diary(logPath)

if ischar(inSfxs)
    fprintf('inSfxs given as string; interpreting as single-element array\n');
    inSfxs = {inSfxs};
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
if ischar(maskVal)
    fprintf('maskVal given as string, intrepreting as number\n');
    maskValIn = maskVal;
    maskVal = str2double(maskVal);
    if isnan(maskVal)
        error('cannot interpret %s as a number\n', maskValIn)
    end
end

files = dir(inDirPath);
fnames = {files.name};

objNames = regexp(fnames, ['^' nameRegexp], 'match');
fileSkipped = cellfun(@isempty, objNames);

fprintf('%d matching files, %d non-matching files\n', ...
    sum(~fileSkipped), sum(fileSkipped));

if sum(fileSkipped) < 100
    fprintf('skipped files: \n');
    for ii=find(fileSkipped)
        fprintf('\t%s\n', fnames{ii});
    end
end

objNames = objNames(~fileSkipped);
objNames = cellfun(@(x)(x{1}), objNames, 'UniformOutput', false);
objNames = unique(objNames);

fprintf('%d unique objects\n', length(objNames));

% stgs = getDefaultSettings();
% stgs.resizeDims = [];
% stgs.unsharpMaskAmt = 0;
% stgsDproj = stgs;
% stgsNDproj = stgsDproj; stgsNDproj.useDeProjectStretch = false;

fprintf('qLvl is %2.4f\n', qLvl);
fprintf('nReps is %d\n', nReps);
fprintf('black-level quantile is %2.4f\n', blackQlvl);

if ~isempty(maskSuffix)
    fprintf('using galaxy mask, suffix = %s, maskVal = %d\n', maskSuffix, maskVal)
end

for objIdx = 1:1:length(objNames)
    curName = objNames{objIdx};
    fprintf('processing: %s\n', curName)
    
    fimgs = cell(1, length(inSfxs));
    bad_read = false;
    for ii=1:1:length(inSfxs)
        fpath = [inDirPath curName inSfxs{ii}];
        fprintf('\treading: %s\n', fpath)
        try
            fimgs{ii} = imread(fpath);
        catch ME
            err_str = sprintf('could not convert %s (unable to read %s: %s)',...
                curName, fpath, ME.message);
            fprintf('WARNING: %s\n', err_str);
            fprintf(errs_file, '%s\n', err_str);
            bad_read = true;
            break;
        end
        if ~isempty(maskSuffix)
            maskPath = [maskDirPath filesep curName maskSuffix];
%             try
%                 mask = imread(maskPath);
%             catch ME
%                 err_str = sprintf('could not convert %s (unable to read mask image %s: %s)',...
%                     curName, maskPath, ME.message);
%                 fprintf('WARNING: %s\n', err_str);
%                 fprintf(errs_file, '%s\n', err_str);
%                 bad_read = true;
%                 break;
%             end
%             mask = (mask == maskVal);
        else
            maskPath = [];
        end
        
    end
    
    if bad_read
        continue;
    end
    
    try
        imgs = convertFromFits(fimgs, maskPath, maskVal, qLvl, nReps, blackQlvl);
    catch ME
        err_str = sprintf('could not convert %s (%s)', curName, ME.message);
        fprintf('WARNING: %s\n', err_str);
        fprintf(errs_file, '%s\n', err_str);
        continue
    end
    
    if length(inSfxs) == 3
        cimg = zeros([size(imgs{1}) 3]);
        cimg(:, :, 1) = imgs{1}; cimg(:, :, 2) = imgs{2}; cimg(:, :, 3) = imgs{3};
    elseif length(inSfxs) == 1
        cimg = imgs{1};
    else
        fprintf('WARNING: inSfxs does not have length 1 or 3; converting to grayscale via averaging\n');
        cimg = zeros(size(imgs{1}));
        for ii=1:length(imgs)
            cimg = cimg + imgs{1};
        end
        cimg = cimg ./ length(imgs);
    end
    
    outPath = [outDirPath curName outSfx];
%     if exist(outPath, 'file')
%         error('%s already exists', outPath);
%     end
    imwrite(cimg, outPath);
%     imwrite(medfilt2(rgb2gray(cimg), [3 3]), outPath);
end

fclose(errs_file);
diary off

end
