function imageAverageFITS(files,cleanUp)
%INPUTS
%   files - a matlab cell array of the files names that are to be averaged
%       the first cell should contain the full path for the output file
%       the rest of the cells should contain full path names to .FITS images.
%   cleanUp - used to clean up tempGZdir if error occurs
if eq(length(files),0)
	errorStruct.message = 'FATAL ERROR: no files inputted';
	errorStruct.identifier = 'imageAverageFITS:noInputFiles';
	error(errorStruct);
end

if eq(length(files),1)
	errorStruct.message = 'FATAL ERROR: no images inputted';
	errorStruct.identifier = 'imageAverageFITS:noInputImages';
	error(errorStruct);
end

outpath = files{1};
images = files(2:length(files));

fprintf('making mean image %s \n', outpath);


if exist(outpath, 'file') == 2 || exist(outpath, 'file') == 7
	errorStruct.message = 'FATAL ERROR: output file exists; not overwritting';
	errorStruct.identifier = 'imageAverageFITS:outpathExists';
	error(errorStruct);
end


nImages = length(images);
firstImageRead = false;
%sum all the fits images
for i = 1:length(images)
	try
		img = fitsread(images{i});
		if eq(firstImageRead,false)
			imgSize = size(img);
			sumImage = double(img);
			firstImageRead = true;
		else
			if ~eq(size(img), imgSize)
				errorStruct.message = ...
					sprintf('FATAL ERROR: %s and %s are different dimensions, cannot make mean',images{1},images{i});
				errorStruct.identifier = 'imageAverageFITS:wrongDimensions';
				error(errorStruct);
			end
			sumImage = sumImage + double(img);
		end
	catch ME
		if strcmp(ME.identifier,'imageAverageFITS:wrongDimensions')
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

fprintf('making mean from %d FITS images \n', nImages)

meanImage = (sumImage / nImages);

fprintf('%s\n', outpath);
fitswrite(meanImage, outpath);
