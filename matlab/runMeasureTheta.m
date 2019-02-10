function runMeasureTheta(indir)
%
%   Runs measureTheta.m on all the galaxies in a standard SpArcFiRe output directory
%   INPUTS
%       indir - directory where SpArcFiRe has output
%           (The .csv will be put in this directory)
%   NOTE
%       this program assumes the diretories follow the following format:
%           galaxyID_colorband
%       where galaxyID is any number of digits and colorband is at least
%       1 character that specifies the colorband
%           EX) 1237661361303650360_r would be a valid directory name
%               for galaxy 1237661361303650360 in the red colorband 
%       if your directories are not in the specified format, this 
%       program may need to be rewritten (see measureTheta.m to see
%       what this program has to pass to it)


files = dir(indir);

galaxyDirs = regexp({files.name},'^[\d]+_[\w]+[\/]?$','match');
galaxyDirs = fixCellArray(galaxyDirs);
if isempty(galaxyDirs)
    %warnStruct.messageID = 'no files were successfully loaded, check that the directories in indir follow the correct format';
    %warnStruct.state = 'on';
    %warnStruct.identifier = 'MATLAB:runMeasureTheta:noGalaxyDirectories';
    %warning(warnStruct);
    warning('no files were successfully loaded, check that the directories in indir follow the correct format');
end

galaxies = regexp((galaxyDirs),'^[\d]+','match');
galaxies = fixCellArray(galaxies);
galaxies = unique(galaxies);
if isempty(galaxies)
    warning('no galaxies were found in indir, check that the directories in indir follow the correct format');
end 

for g = 1:length(galaxies) %TODO add error catching
    try
        galaxyID = galaxies{g};
        disp(galaxyID);

        clusNames = regexp(galaxyDirs,sprintf('^%s.*',galaxyID),'match');
        clusNames = fixCellArray(clusNames);
        clusNames = extractAfter(clusNames,'_'); %this assumes there is no _ in the galaxyID
        clusNames(strcmp(clusNames,'mean')) = []; %TODO check if no mean
        disp(clusNames);

        guideClusMaskFile = sprintf('%s/%s_mean/%s_mean-H_clusMask-merged.png',indir,galaxyID,galaxyID);

        clusMaskFiles = cellfun(@(clusName) sprintf('%s/%s_%s/%s_%s-H_clusMask-merged.png',indir,galaxyID,clusName,galaxyID,clusName),clusNames,'UniformOutput',false);
        for i = 1:length(clusMaskFiles)
            disp(clusMaskFiles{i});
        end
        measureTheta(indir,galaxyID,clusNames,guideClusMaskFile,clusMaskFiles)
    catch ME
        warning(sprintf('skipping galaxy %s; orginal error message(%s)',galaxyID,ME.message));
    end
end

end

function[fixed] = fixCellArray(cellArray)
%converts ouput of regexp to a cell array of char row vectors with no empty cells
    fixed = cellArray(~cellfun('isempty',cellArray));
    fixed = cellfun(@(c) c{1},fixed,'UniformOutput',false);
end
