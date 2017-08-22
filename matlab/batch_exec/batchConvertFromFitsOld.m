function batchConvertFromFitsOld(inDirPath, outDirPath, nameRegexp, inSfxs, outSfx, convStgs)
% INPUTS:
%   inDirPath: directory containing the fits files to convert
%   outDirPath: directory where the output files should be placed
%   nameRegexp: regular expression describing the name of each object.
%       input files with the same match to this regular expression will be
%       combined into a single composite image.
%   inSfxs: cell array of suffixes for each waveband. For each set of 
%       images that has the same match to nameRegexp, there should be a 
%       1:1 correspondence between the remainder of the file name and an
%       entry in this suffix list.
%   outSfx: suffix for output images (the earlier part of the file name is
%       the match to nameRegexp), including the file extension

% blq = 0.2;
% % wlq = 0.999;
% wlq = 1;

convStgs

inDirPath = [inDirPath filesep];
outDirPath = [outDirPath filesep];

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

for objIdx = 1:1:length(objNames)
    curName = objNames{objIdx};
    
    fimgs = cell(1, length(inSfxs));
    for ii=1:1:length(inSfxs)
        fpath = [inDirPath curName inSfxs{ii}];
        try
            fimgs{ii} = im2double(imread(fpath));
        catch ME
            error('unable to read %s', fpath);
        end
    end
    
%     for ii=1:1:length(fimgs)
%         fimgs{ii} = medfilt2(fimgs{ii}, [7 7]);
%     end

%     varargin{:}
    imgs = convertFromFits(fimgs, false, convStgs);
    
    if length(inSfxs) == 3
        cimg = zeros([size(imgs{1}) 3]);
        cimg(:, :, 1) = imgs{1}; cimg(:, :, 2) = imgs{2}; cimg(:, :, 3) = imgs{3};
    elseif length(inSfxs) == 1
        cimg = imgs{1};
    else
        error('inSfxs must have length 1 or 3');
    end
    
    outPath = [outDirPath curName outSfx];
%     if exist(outPath, 'file')
%         error('%s already exists', outPath);
%     end
    imwrite(cimg, outPath);
%     imwrite(medfilt2(rgb2gray(cimg), [3 3]), outPath);
end

end