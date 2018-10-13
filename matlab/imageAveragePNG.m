function imageAveragePNG(files,cleanUp)
%INPUTS
%   files - a matlab cell array of filenames
%       the first cell should contain the full path of the output file
%       the rest of the cells should be full filenames of .FITS files
if eq(length(files),0)
	errorStruct.message = 'FATAL ERROR: no files inputted';
	errorStruct.identifier = 'imageAveragePNG:noInputFiles';
	error(errorStruct);
end

if eq(length(files),1)
	errorStruct.message = 'FATAL ERROR: no images inputted';
	errorStruct.identifier = 'imageAveragePNG:noInputImages';
	error(errorStruct);
end

outpath = files{1};
images = files(2:length(files));

fprintf('making mean image %s \n', outpath);


if exist(outpath, 'file') == 2 || exist(outpath, 'file') == 7
	errorStruct.message = 'FATAL ERROR: output file exists; not overwritting';
	errorStruct.identifier = 'imageAveragePNG:outpathExists';
	error(errorStruct);
end


nImages = length(images);
firstImageRead = false;
%sum all the png images
for i = 1:length(images)
	try
		img = imread(images{i});
		if eq(firstImageRead,false)
			imgSize = size(img);
			sumImage = double(img);
			firstImageRead = true;
		else
			if ~eq(size(img), imgSize)
				errorStruct.message = ...
					sprintf('FATAL ERROR: %s and %s are different dimensions, cannot make mean',images{1},images{i});
				errorStruct.identifier = 'imageAveragePNG:wrongDimensions';
				error(errorStruct);
			end
			sumImage = sumImage + double(img);
		end
	catch ME
		if strcmp(ME.identifier,'imageAveragePNG:wrongDimensions')
			rethrow(ME);
		end
		warning('Error loading file %s; skipping file : original error message(%s)',images{i},ME.message);	
		nImages = nImages - 1;
	end
end

if eq(nImages,0)
	errorStruct.message = 'FATAL ERROR: no images successfully loaded';
	errorStruct.identifier = 'imageAverageFITS:noFilesLoaded';
	error(errorStruct);
end

fprintf('making mean from %d PNG images \n', nImages)

meanImage = uint8(sumImage / nImages);

fprintf('%s\n', outpath);
imwrite(meanImage, outpath);
