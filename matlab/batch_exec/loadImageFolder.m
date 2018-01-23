function [imgPaths, imgNames, imgFileExts, imgMtxs] = ...
    loadImageFolder(folderPath, fileTypes)
% TODO: allow for other types of ls output (from non-Windows systems)

if ~isdir(folderPath)
    error('%s is not a folder', folderPath);
end

if folderPath(end) ~= filesep
    folderPath = [folderPath filesep];
end

prevDir = pwd();
cd(folderPath);

fileNames = ls(folderPath);
nFiles = size(fileNames, 1);

imgNames = cell(nFiles, 1);
types = cell(nFiles, 1);
imgMtxs = cell(nFiles, 1);
isImage = false(nFiles, 1);

for ii=1:1:nFiles
    curFileStr = strtrim(fileNames(ii, :));
    extPosn = strfind(curFileStr, '.');
    if isempty(extPosn)
        imgNames{ii} = curFileStr;
        types{ii} = [];
    else
        extPosn = extPosn(end);
        imgNames{ii} = curFileStr(1:extPosn-1);
        types{ii} = curFileStr(extPosn+1:end);
    end
    
    curIsImg = ~isempty(strmatch(types{ii}, fileTypes, 'exact'));
    if nargout >= 4 && curIsImg
        fprintf('Loading image: %s\n', curFileStr);
        img = im2double(imread(curFileStr));
        imgSz = size(img);
        if (length(imgSz) == 3 && imgSz(3) == 3) % color image, make gray
            img = rgb2gray(img);
        end
        imgMtxs{ii} = img;
        isImage(ii) = true;
    elseif curIsImg
        fprintf('Found image: %s\n', curFileStr);
        isImage(ii) = true;
    else
        fprintf('Skipping file: %s\n', curFileStr);
    end
end

imgNames = imgNames(isImage);
imgFileExts = types(isImage);
imgPaths = cellfun(@(name, extn)(strcat(folderPath, name, '.', extn)),...
    imgNames, imgFileExts, 'UniformOutput', false);
if nargout >= 4
    imgMtxs = imgMtxs(isImage);
end

cd(prevDir);

end