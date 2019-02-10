function measureTheta(outdir, galaxyID, clusNames, guideClusMaskFile, clusMaskFiles)
%
%	For each radius, finds the difference of the theta of clusters in the colorband and the guide image clusters
%	if the two are associated with one another (if they overlap at any point).
%	INPUTS
%		outdir - directory where csv file will be outputted (must exist)
%		galaxyID - number of the galaxy
%		guideClusMaskFile - file of mean image clusters (should end in -H_clusMask-merged.png)
%		clusNames - cell array with the colorbands. The colorbands should be in an order that is associated
%           with the files.
%		clusMaskFiles - cell array of files that will have their theta measured
%	
%	OUTPUTS
%		creates a csv of the theta values. if the files already exists, then it is appended to said file
%		
%

nfiles = size(clusMaskFiles,2);

if eq(nfiles,0)
	errorStruct.message = 'FATAL ERROR: no files were inputted, cannot find theta';
	errorStruct.identifier = 'MATLAB:MeasureTheta:noFilesInputted';
	error(errorStruct);
end

if exist(outdir) ~= 7
	errorStruct.message = sprintf('FATAL ERROR: %s does not exist',outdir);
	errorStruct.identifier = 'MATLAB:MeasureTheta:OutdirDoesNotExist';
	error(errorStruct)
end

if size(clusNames,2) ~= nfiles
	errorStruct.message = 'FATAL ERROR: there must be equal amount of colorbands and files';
	errorStruct.identifier = 'MATLAB:MeasureTheta:incorrectNumberOfColorbands';
	error(errorStruct);
end

try
	guideClusMask = imread(guideClusMaskFile);
catch ME
	if startsWith(ME.identifier,'MATLAB:imagesci:imread')
		msg = ['Unable to load guideClusMask, unable to perform theta calculations'];
		causeException = MException('MATLAB:measureTheta:failToLoadGuideClusMask',msg);
		ME = addCause(ME,causeException);
	end
	rethrow(ME);
end
clusList = cell(1,nfiles);

for i = 1:nfiles
try
	clusMask = imread(clusMaskFiles{i});
	clusList{i} = clusMask;
catch ME
	warning(sprintf('error loading file, skipping; original error message (%s)',ME.message));
	nfiles = nfiles-1;
	clusNames{i} = [];
	end
end

if eq(nfiles,0)
	errorStruct.message = 'FATAL ERROR: no files were loaded, cannot find theta';
	errorStruct.identifier = 'MATLAB:MeasureTheta:noFilesLoaded';
	error(errorStruct);
end

clusList = clusList(~cellfun('isempty',clusList));
clusNames = clusNames(~cellfun('isempty',clusNames));

outputFilename = sprintf('%s/theta.csv',outdir);

if exist(outputFilename) == 2
	fileID = fopen(outputFilename, 'a');
else
	fileID = fopen(outputFilename, 'w');
	fprintf(fileID, 'galaxyID,radius,guideCluster,colorband,cluster,meanGuideTheta,meanColorTheta\n');
end

radiusMtxs = createCircleMtxs(size(clusList{1},2),size(clusList{1},1));
for r = 1:127
	radiusMatrix = radiusMtxs(:,:,r+1); %r+1 because the radiusMtxs start at a radius of 0
	

    % find guide clusters that touch r
    guideClusMask_columns = reshape(guideClusMask, [], 3);
    [unique_guide_clusters, m, n] = unique(guideClusMask_columns, 'rows');

    for g = 2:size(unique_guide_clusters,1)  % first element/color is just (0,0,0), no use for it
        % find guide clusters that touch r
        guide_cluster = n == g;
        guide_cluster = reshape(guide_cluster, 256, 256);
        guideOverlapAtRadius = radiusMatrix & guide_cluster;
        if nnz(guideOverlapAtRadius) == 0
            continue;
        end

        mean_guide_theta = getMeanTheta(guideOverlapAtRadius);

        numColors = length(clusList);
        for w = 1:numColors
            waveband = clusList{w};

            % iterate over clusters in color
            waveband_column = reshape(waveband, [], 3);
            [unique_color_clusters, m_cc, n_cc] = unique(waveband_column, 'rows');

            for c = 2:size(unique_color_clusters,1)
                color_cluster = n_cc == c;
                color_cluster = reshape(color_cluster,256,256);
                logicOverlap = color_cluster & guide_cluster;
                radiusOverlap = color_cluster & radiusMatrix;
                if nnz(logicOverlap) == 0 || nnz(radiusOverlap) == 0
                    continue;
                end
                mean_color_theta = getMeanTheta(radiusOverlap);

                fprintf(fileID, '%s,%d,%d,%c,%d,%6.5f,%6.5f\n',galaxyID, r, g-1, clusNames{w}, c-1, mean_guide_theta, mean_color_theta);
            end
        end
    end
end

fclose(fileID);
end

%this version finds the theta at every point and finds the mean of their thetas 
%while the version below it finds the mean position of all the nonzero elements 
%and then finds the theta of that position
%they begin to differ at the 3rd decimal place
%function mean_theta = getMeanTheta(cluster)
%    theta = 0;
%    [row, col, v] = find(cluster);
%    for i = 1:length(row)
%        y = 128 - row(i);
%        x = col(i) - 128;
%        theta = theta + atan2(y, x);
%    end
%    mean_theta = theta / length(row);
%end

function mean_theta = getMeanTheta(cluster)
    theta = 0;
    [row, col, v] = find(cluster);
	x = zeros(size(col));
	y = zeros(size(row));
    for i = 1:length(row)
        y(i) = 128 - row(i);
        x(i) = col(i) - 128;
    end
	mean_theta = atan2(mean(y),mean(x));
end



function[circleMtxs] = createCircleMtxs(sizeX,sizeY)
%createCircleMtxs makes matricies with circles of radius 0-129 pixels 
%used to compare against the colorband clusters
%INPUTS:
%	sizeX: number of cells along the x-axis of the picture
%	sizeY: number of cells along the y-axis of the picture
%Outputs:
%	circleMtxs: 3-D matrix where the nth matrix in the 3rd dimension 
%		contains a matrix with a circle of radius of n pixels

circleMtxs = zeros(sizeX,sizeY);
	for radius=0:sizeX/2 %assuming square
		circleMtx = circle(sizeX/2,sizeY/2,radius);
		circleMtxs(:,:,radius+1) = circleMtx;
	end
end

function[circle] = circle(x,y,radius)
%Creates a matrix with a circle with a specificied radius
%INPUTS:
%	x: x-position of the center of the picture
%	y: Y-poistion of the center of the picture
%	radius: radius of the circle
%OUPUTS:
%	circle: the matrix that contains a circle with a 
%		specified radius

circle = zeros(x*2,y*2);
for cellX=1:x*2
	for cellY=1:y*2
		dist = round(sqrt((cellX-x)^2 + (cellY-y)^2));
		if dist == radius
			circle(cellX,cellY) = 1;
		end
	end
end
end
