function genCsvFiles(dirPath)

    contents = dir(dirPath);
    directories = contents(arrayfun(@(x)(x.isdir && ~isempty(findstr('=', x.name))), contents));
    for ii = 1:1:length(directories)
        dirName = directories(ii).name;
        [dirPath filesep dirName]
        writeGxyParamsToCsv([dirPath filesep dirName], dirName);
        copyfile([dirPath filesep dirName filesep dirName '.csv'],...
            [dirPath filesep dirName '.csv']);
    end
    
end