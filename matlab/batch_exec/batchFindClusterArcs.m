function [gxyParamsList, allClusMtxs] = ...
    batchFindClusterArcs(iptImgDir, imgNames, imgSuffix, starMaskSuffix, procNum, numProcs, stgs, outputDir, params_filename, recoveryFileName, out_id)

usage_msg = ['parameters: input_img_dir img_names img_suffix ',...
    'star_mask_suffix proc_num num_procs stgs output_dir\n',...
    '\tinput_img_dir is a path to the directory where the input PNG or JPG images are to be found\n',...
    '\timg_names is a text-file path or cell array giving the image names\n',...
    '\timg_suffix is what should be appended to each image name ',...
    'to produce the image file name (within input_img_dir)\n',...
    '\tstar_mask_suffix is the suffix for star mask files, or "NONE"\n',...
    '\tproc_num is the number of this process, starting at 1\n',...
    '\tnum_procs is the total number of processes being run. If there is more ',...
    'than one process, each should be given the same img_names.\n',...
    '\tstgs is a settings struct or path to a settings file, or "NONE" to use defaults.\n',...
    '\toutput_dir is where the output files should be placed.\n',...
    'If this is being run in standalone mode, the corresponding script may ',...
    'need additional parameters.\n'];
% out_id is for findClusterArcsServer and should be specified iff the arc
% finding is used as a server and procNum/numProcs are not specified
% TODO: out_id no longer used -- remove it

if nargin < 8
    fprintf(usage_msg)
end

error(nargchk(8, 11, nargin));

if strcmpi(stgs, 'NONE') || isempty(stgs)
    fprintf('No fit settings specified. Using defaults.\n');
    stgs = getDefaultSettings();
elseif ischar(stgs)
    fprintf('stgs specified as a path string, attempting to load from disk.\n');
    load(stgs, 'stgs');
end
if ~isempty(procNum)
    stgs = processOverridesFromEnv(stgs);
end

if ~isdir(outputDir)
    error('%s is not a directory', outputDir);
end

if iptImgDir(end) ~= filesep
    iptImgDir = [iptImgDir filesep];
end

if outputDir(end) ~= filesep
    outputDir = [outputDir filesep];
end

if strcmpi(params_filename, 'NONE');
        stdzParamsList = {};
end

if ~strcmpi(params_filename, 'NONE') 
        paramsPath = [iptImgDir params_filename];
    if exist(paramsPath, 'file')
        elps_fit_file = fopen(paramsPath, 'rt');
        elps_fit_str = fgetl(elps_fit_file);
        fclose(elps_fit_file);
        stdzParamsList = {stringToDblStruct(elps_fit_str)};
    elseif ~exist(paramsPath, 'file')
        error('File "%s" does not exist.\n', paramsPath);
    end
end

if nargin < 9
    stdzParamsList = {};
end

if nargin < 10
    recoveryFileName = [];
end

% hack so that an empty value can be passed when used on the command line
if strcmpi(starMaskSuffix, 'NONE')
    starMaskSuffix = [];
end

if strcmpi(params_filename, 'NONE');
        params_filename = {};
end

validateattributes(imgNames, {'cell', 'char'}, {});

if ischar(procNum)
    procNumIn = procNum;
    procNum = str2double(procNum);
    if isnan(procNum)
        fprintf('could not interpret procNum as a number, checking environment variable\n');
        envVarName = procNumIn;
        procNum = getenv(envVarName);
        if isempty(procNum)
            error('no environment variable %s for procNum\n', envVarName);
        else
            fprintf('found environment variable for procNum: %s=%s',...
                envVarName, procNum);
            procNum = str2double(procNum);
            if isnan(procNum)
                error('environment-variable value for procNum is not interpretable as a number');
            end
        end
    end
end
if ischar(numProcs)
    numProcs = str2double(numProcs);
    if isnan(numProcs)
        error('could not interpret numProcs as a string');
    end
end

if ischar(imgNames)
    fprintf('imgNames specified as a path string, attempting to load from disk.\n');
    load_ok = true;
    try
        lres = load(imgNames);
        flds = fields(lres);
    catch ME
        load_ok = false;
        fprintf('could not read mat-file for image names (%s)\n', ME.message)
    end
    if load_ok
        ispath = cellfun(@(x)(~isempty(regexpi(x, 'Names$'))), flds);
        if sum(ispath) == 0
            error('no variables with suffix "Names" found in mat-file %s',...
                imgNames)
        end
        if sum(ispath) > 1
            error('multiple variables with suffix "Names" found in mat-file %s',...
                imgNames)
        end
        fldname = flds(ispath);
        imgNames = getfield(lres, fldname{1});
    else
        namesFilename = imgNames;
        fprintf('could not read %s as MAT-file, trying as text.\n', namesFilename);
        imgNames = {};
        namesFile = fopen(namesFilename, 'rt');
        nxtLine = fgetl(namesFile);
        linesRead = 0;
        while nxtLine ~= -1
            imgNames{end+1} = nxtLine;
            nxtLine = fgetl(namesFile);
            linesRead = linesRead + 1;
        end
        imgNames = imgNames(~cellfun(@isempty, imgNames));
        fclose(namesFile);
    end
end
% if ischar(imgPaths)
%     fprintf('imgPaths specified as a path string, attempting to load from disk.\n');
%     load(imgPaths, 'imgPaths');
% end
% if ischar(imgNames)
%     fprintf('imgNames specified as a path string, attempting to load from disk.\n');
%     load(imgNames, 'imgNames');
% end
if ischar(stdzParamsList)
    fprintf('stdzParamsList specified as a path string, attempting to load from disk.\n');
    load(stdzParamsList, 'stdzParamsList');
    % load([iptImgDir stdzParamsList], 'stdzParamsList'); % ADDED BY ARA
end

% if length(imgPaths) ~= length(imgNames)
%     error('image array and image name array must have same length');
% end

if nargout >= 2
    allClusMtxs = cell(length(imgNames), 1);
end

% allowArcBeyond2pi = stgs.allowArcBeyond2pi;
groupOutputByProcess = stgs.groupOutputByProcess;
groupOutputByInputImage = stgs.groupOutputByInputImage;
failWhenNoStarmaskFound = stgs.failWhenNoStarmaskFound;

rootOutputDir = outputDir;
if groupOutputByProcess
    outputDir = [outputDir sprintf('proc%04d', procNum) filesep];
    mkdir(outputDir);
end
paramsMatOutputDir = [outputDir 'matout' filesep];
if ~exist(paramsMatOutputDir, 'dir')
    mkdir(paramsMatOutputDir);
end

if ~isempty(procNum)
    if isempty(numProcs)
        error('if procNum is specified, numProcs must also be specified');
    end
    outfileFmt = '%04d_settings.txt';
    if ~isempty(recoveryFileName)
        outfileFmt = '%04d_settings_recovery.txt';
    end
    diary([rootOutputDir filesep sprintf(outfileFmt, procNum)]);
    fprintf('%s\n', datestr(clock));
    fprintf('procNum: %d\n', procNum);
    fprintf('numProcs: %d\n', numProcs);
    fprintf('input image directory: %s\n', iptImgDir);
    if isempty(starMaskSuffix)
        fprintf('not using star masks\n');
    else
        fprintf('star mask suffix: %s\n', starMaskSuffix)
    end
    stgs
    diary off
    
    imgNames = imgNames(procNum:numProcs:length(imgNames));
    if ~isempty(stdzParamsList)
        stdzParamsList = stdzParamsList(procNum:numProcs:length(stdzParamsList));
        size(stdzParamsList)
    end
else
    if (nargin < 11) || isempty(out_id)
        error('if procNum and numProcs not given, then out_id must be given');
    end
    if iscell(imgNames) && (length(imgNames) ~= 1)
        error('if procNum and numProcs not given, then imgNames must have exactly one element');
    end
end

gxyParamsList = cell(length(imgNames), 1);
startIdx = 1;
if ~isempty(recoveryFileName)
    load([outputDir filesep recoveryFileName]);
    startIdx = find(cellfun(@isempty, gxyParamsList), 1, 'first');
end
if isempty(stdzParamsList)
    stdzParamsList = cell(size(gxyParamsList));
end

imgSuffixLen = length(imgSuffix);
gxyOutputCtrl = struct('writeImages', true, ...
        'displayFigures', false, 'writeTxt', false);
for imgIdx = startIdx:1:length(imgNames)
    gxyName = imgNames{imgIdx};
    gxyPath = [iptImgDir gxyName imgSuffix];
    fprintf('processing %s (image %d of %d, %s)...\n', ...
        gxyName, imgIdx, length(imgNames), gxyPath);
    
    if ~isempty(starMaskSuffix) 
        starMaskPath = [iptImgDir gxyName starMaskSuffix];
        if exist(starMaskPath, 'file')
            fprintf('star mask found: %s\n', starMaskPath);
            starMask = imread(starMaskPath);
        elseif failWhenNoStarmaskFound
            error('star mask suffix specified, but no star mask found at %s\n',...
                starMaskPath)
        else
            fprintf(2, ['WARNING: %s: star mask suffix specified, but no star mask found. '...
                'Settings specify to continue anyway.\n'], gxyName);
            starMask = [];
        end
    else
        fprintf('no star mask suffix given, not looking for a star mask\n')
        starMask = [];
    end
    
    curImg = im2double(imread(gxyPath));
    imgSz = size(curImg);
    if (length(imgSz) == 3 && imgSz(3) == 3) % color image, make gray
       curImg = rgb2gray(curImg);
    end
    
    tStart = tic;
    cpuTStart = cputime();
    if groupOutputByInputImage
        imgOutputDir = [outputDir gxyName filesep];
        [success, message] = mkdir(imgOutputDir);
        if ~success
            fprintf(2, 'WARNING: %s: could not create output directory (%s)',...
                gxyName, message);
            imgOutputDir = outputDir;
        end
    else
        imgOutputDir = outputDir;
    end
    try
%         stdzParamsList
        [lgspParams, lgspBounds, errs, used2rev, failed2rev, hasBadBounds, barInfo, ...
                clusMtxs, gxyParams, imgAutoCrop, barInds, barUsed] = ...
            findClusterArcs(curImg, stgs, ...
            gxyName, gxyOutputCtrl, imgOutputDir, starMask, stdzParamsList{imgIdx});
%         isBar = false(size(clusMtxs, 3), 1); isBar(barInds) = true;
        if nargout >= 2
            allClusMtxs{imgIdx} = clusMtxs;
        end    
    catch ME
        gxyParams.name = gxyName;
%         if ~isempty(strmatch('input rejected', ME.message))
        if ~isfield('gxyParams', 'fit_state');
            gxyParams.fit_state = ME.message;
        else
            gxyParams.fit_state = sprintf('internal error (%s)',...
                ME.message);
        end
        gxyParams.error = ME;
        gxyParamsList{imgIdx} = gxyParams;
        fprintf(2, 'WARNING: %s: could not perform fitting (%s)\n', ...
            gxyName, ME.message);
        for ii=1:1:length(ME.stack)
            stack_item = ME.stack(ii);
            fprintf(2, '\tin function %s, line %d\n', ...
                stack_item.name, stack_item.line);
        end
        continue;
    end
    cpu_time = cputime() - cpuTStart;
    wall_time = toc(tStart);
    gxyParams.name = gxyName;
    gxyParams.fit_state = 'OK';
    gxyParams.cpu_time = cpu_time;
    gxyParams.wall_time = wall_time;
    % TODO: current code's method of splitting bar/non-bar entries is ugly
    % - fix this!
%     showClustersFromMtxs(clusMtxs, [size(clusMtxs(:, :, 1))]); % can this be removed?
    gxyParams = getGalaxyParams(lgspParams, lgspBounds, errs,...
        used2rev, failed2rev, hasBadBounds, barUsed, barInfo, barInds, imgAutoCrop, ...
        clusMtxs, stgs, gxyParams);
    if gxyOutputCtrl.writeTxt
        writeGalaxyParams(gxyParams, gxyName, imgOutputDir);
    end
%     close all
    gxyParamsList{imgIdx} = gxyParams;
    if mod(imgIdx, 100) == 0 && ~isempty(procNum)
        save([outputDir 'gxyParamsList-part-' sprintf('%04d', procNum) ...
            '-of-' sprintf('%04d', numProcs) '-checkpoint-to-'...
            sprintf('%04d', imgIdx)], 'gxyParamsList');
    end
    close all
end

if ~isempty(procNum)
    out_name = ['gxyParamsList-part-' sprintf('%04d', procNum) ...
        '-of-' sprintf('%04d', numProcs)];
else
%     out_name = ['gxyParamsList-part-' out_id];
    out_name = ['gxyParamsList-part-' imgNames{1}];
end
save([paramsMatOutputDir out_name], 'gxyParamsList');

end
