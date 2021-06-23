function [] = createMeans(varargin)
%reads from stdin

fprintf('--- Creating Means ---\n')

tempGZdir = '/tmp/SpArcFiRe-mean-unzipped-files'; %this is where .gz files will be unzipped to
%tempGZdir = '/home/barkleya/SpArcFiRe/SpArcFiRe-new-code/theta/temp';
cleanUp = onCleanup(@() removeTempFiles(tempGZdir));
if length(varargin) > 0 %use arguements

	%fprintf('ARGS\n');
	%check stuff
	checkOutpath(varargin{1});
	if size(varargin,2) < 2
		errorStruct.message = 'FATAL ERROR: NO IMAGES INPUTTED';
		errorStruct.identifier = 'createMeans:noInputImages';
		error(errorStruct);
	end
	fprintf('checking %s\n', varargin{1});

	[flags,varargin] = checkConsistentFileTypes(varargin,tempGZdir);
	if eq(flags(1),1)
		imageAveragePNG(varargin,cleanUp);
	elseif eq(flags(2),1)
		imageAverageFITS(varargin,cleanUp);
	else
		errorStruct.message = 'FATAL ERROR: ALL FILES MUST BE ALL PNGS OR ALL FITS';
		errorStruct.identifier = 'createMeans:inconsistentFileTypes';
		error(errorStruct);
	end

else %read from stdin
	
	while true 
		try
			fprintf('-----\n');
			line = input(' ','s');
			%fprintf('%s\n',line);
			%parse line
			files = strsplit(line);
			if eq(size(files),0)
				errorStruct.message = 'FATAL ERROR: NO FILES INPUTTED';
				errorStruct.identifier = 'createMeans:noInputFiles';
				error(errorStruct);
			end
			fprintf('making %s\n', files{1});
			checkOutpath(files{1});

			[flags,files] = checkConsistentFileTypes(files,tempGZdir);
			if eq(flags(1),1)
				imageAveragePNG(files,cleanUp);
			elseif eq(flags(2),1)
				imageAverageFITS(files,cleanUp);
			else
				errorStruct.message = 'FATAL ERROR: ALL FILES MUST BE ALL PNGS OR ALL FITS';
				errorStruct.identifier = 'createMeans:inconsistentFileTypes';
				error(errorStruct);
			end
		catch ME
			if (strcmp(ME.identifier, 'MATLAB:REGEXP:invalidInputs'))%exit if reach end of stdin
				break
			%{elseif startsWith(ME.identifier,'createMeans:') | ...
				startsWith(ME.identifier,'imageAveragePNG') | ...
				startsWith(ME.identifier,'imageAverageFITS')%change to warning
				if (strcmp(ME.identifier,'createMeans:noInputFiles'))
					warning(ME.message);
				end
				warning('Skipping %s ; orginal error message (%s)', files{1},ME.message);
				%}
			else
				if (strcmp(ME.identifier,'createMeans:noInputFiles'))
					warning(ME.message);
				end
				warning('Skipping %s ; orginal error message (%s)', files{1},ME.message);
			end
		end
	end
	exit;

end

end


function [] = checkOutpath (outpath)
	fprintf('checking outpath\n');
	if exist(outpath, 'file' ) == 2 | exist (outpath, 'file') == 7 
		errorStruct.message = 'FATAL ERROR: OUTPUT FILE EXISTS; NOT OVERWRITTING';
		errorStruct.identifier = 'createMeans:outpathExists';
		error(errorStruct);
	end
end

function [flags,files] = checkConsistentFileTypes(files,tempGZdir)
	fprintf('checking file types\n');
	flags = [true,true]; %all PNG, ALL FITS
    %disp(files);
	for i = 1:size(files,2)
    file = files{i};
        try
            if length(file) >=2 & eq(file(end-2:end),'.gz') & ~eq(i,1)
                if exist(tempGZdir) ~= 7
                    mkdir(tempGZdir);
                end
                newFile = gunzip(file,tempGZdir);
                file = newFile{1};
                files{i} = newFile{1};
                    
            end

            if length(file) <= 3 | ~eq(file(end-3:end),'.png')
                flags(1) = 0;
            end

            if length(file) <= 4 | ~eq(file(end-4:end),'.fits')
                flags(2) = 0;
            end
         catch ME
            warning('Skipping %s ; orginal error message(%s)',file,ME.message);
        end

	end

end

function [] = removeTempFiles(tempGZdir)
    fprintf('removing %s\n',tempGZdir);
    if exist(tempGZdir) == 7
        delete(sprintf('%s/*',tempGZdir));
        rmdir(tempGZdir);
    end
end
