function compareFitsConversion(dirPath1, dirPath2, nameRegexp, sfx1, sfx2, outDir)
% Not used in main flow

dirPath1 = [dirPath1 filesep];
dirPath2 = [dirPath2 filesep];
outDir = [outDir filesep];

files1 = dir(dirPath1);
fnames1 = {files1.name};
objNames1 = regexp(fnames1, ['^' nameRegexp], 'match');
objNames1 = objNames1(~cellfun(@isempty, objNames1));
objNames1 = cellfun(@(x)(x{1}), objNames1, 'UniformOutput', false);
objNames1 = unique(objNames1);

files2 = dir(dirPath2);
fnames2 = {files2.name};
objNames2 = regexp(fnames2, ['^' nameRegexp], 'match');
objNames2 = objNames2(~cellfun(@isempty, objNames2));
objNames2 = cellfun(@(x)(x{1}), objNames2, 'UniformOutput', false);
objNames2 = unique(objNames2);

mtchNames = objNames1(ismember(objNames1, objNames2));

for ii=1:1:length(mtchNames)
    img1 = imread([dirPath1 mtchNames{ii} sfx1]);
    img2 = imread([dirPath2 mtchNames{ii} sfx2]);
    if size(img1, 3) == 3
        img1 = rgb2gray(img1);
    end
    if size(img2, 3) == 3
        img2 = rgb2gray(img2);
    end
    outName1 = [outDir mtchNames{ii} '_set1.png'];
    outName2 = [outDir mtchNames{ii} '_set2.png'];
    if exist(outName1, 'file')
        error('%s already exists', outName1);
    end
    if exist(outName2, 'file')
        error('%s already exists', outName2);
    end

    imwrite(img1, outName1);
    imwrite(img2, outName2);
end

end