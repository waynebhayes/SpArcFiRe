function [result] = getGalfitFitQuality(img,clusReproj,outputPath,gxyParams)
% Performs and return a GALFIT model fitting on the input image.
% INPUTS:
%   img: the image to use for preprocessing
%   outputPath: the path to write and store the GALFIT .feedme/.fits input
%       and .fits output files
%   gxyParams: galaxy parameters used for generating the GALFIT initial
%       fit .feedme input file.
% OUTPUTS:
%   img: the preprocessed image
%   imgNoUsm: the preprocessed image, without applying the unsharp mask
%   gxyParams: structure containing some information about the ellipse fit 

    % Big edits coming in - Matthew 1/12/21
    % Fixing where GALFIT runs and retaining only the final FITS
    % for posterity and potential analysis/reporduction

    % Changing outputPath to drop everything into tmp
    
    OGoutputPath = outputPath;
    gal_name = split(string(outputPath), '/');
    gal_name = gal_name(end);

    [status, outputPath] = system('mktemp -d /tmp/galfit_junk_XXXXXX');
    outputPathGname = strtrim(string(outputPath)) + '/' + string(gal_name);

    disp('Checking outputPathGname - ' + outputPathGname)
    % outputPath = char(outputPath);
    fitswrite(img, [char(outputPathGname) '_galfit_input.fits']);
    %galfitTemplateFilename = ['/home/' getenv('USER') '/bin/GalfitTemplates/template.feedme']; %added this to call galfit using correct path
    galfitTemplateFilename = string([getenv('SPARCFIRE_HOME') '/GalfitTemplates/template.feedme']);
    disp('Reading GALFIT template file: ' + galfitTemplateFilename)
    galfitTemplate = fopen(galfitTemplateFilename, 'r');
    text = fread(galfitTemplate, '*char')';
    fclose(galfitTemplate);

    % Fill out template by replacing variables with their actual values
    muFit = gxyParams.iptCtrXY;
    nRows = gxyParams.iptSz(1);
    nCols = gxyParams.iptSz(2);
    text = strrep(text, '$input_name', outputPathGname + '_galfit_input.fits');
    text = strrep(text, '$output_name', [OGoutputPath '_galfit_output.fits']);
    text = strrep(text, '$x_center', num2str(muFit(1)));
    text = strrep(text, '$y_center', num2str(size(img, 1) - muFit(2) + 1));
    text = strrep(text, '$radius_bulge', num2str(gxyParams.bulgeMajAxsLen / 4));
    text = strrep(text, '$radius_disk', num2str(gxyParams.diskMajAxsLen / 4));
    text = strrep(text, '$sersic_index', num2str(1.0));
    text = strrep(text, '$axis_ratio', num2str(gxyParams.diskAxisRatio));
    text = strrep(text, '$position_angle', num2str(90 - rad2deg(gxyParams.diskMajAxsAngleRadians)));
    text = strrep(text, '$x_max', num2str(nRows));
    text = strrep(text, '$y_max', num2str(nCols));

    disp('Writing GALFIT template file: ' + outputPathGname + '.feedme')
    galfitInput = fopen(outputPathGname + '.feedme','wt');
    fwrite(galfitInput, text);
    fclose(galfitInput);

    % Run galfit
    %galfitCommand = ['/home/' getenv('USER') '/bin/galfit ' outputPathGname '.feedme'];
    % ^Will 9/30/19: This line of code was added a few months ago because of how matlab evaluates path's at compile time
    % see documentation/stuck_try_this.txt for more explanation
    galfitCommand = char(getenv('SPARCFIRE_HOME') + "/scripts/galfit " + outputPathGname + ".feedme");

    % Grabbing current directory before cd-ing to tmp
    former_dir = cd(char(strtrim(outputPath)));

    system(galfitCommand);

    %disp(former_dir)
    %disp(class(former_dir))
    cd(former_dir);

    % Retrieve input/model subtraction
    model = mat2gray(fitsread([OGoutputPath '_galfit_output.fits'], 'image', 2));
    residual = mat2gray(fitsread([OGoutputPath '_galfit_output.fits'], 'image', 3));
    binarizedClusters = imbinarize(rgb2gray(clusReproj), 0.1);
    binarizedModel = imbinarize(model);
    binarizedResidual = imbinarize(residual);
    relevantElements = binarizedModel.*binarizedResidual;
    
   % disp(mean2(model));
   % disp(mean2(residual));
   % disp(mean2(binarizedResidual));
   % residualMedian = medfilt2(residual, 'symmetric');
   % residual2Median = medfilt2(binarizedResidual, 'symmetric');
   % figure(1);
   % imshow(residual);
   % figure(2);
   % imshow(residualMedian);
   % figure(3);
   % imshow(residual2Median);
   %  figure(3);
   %  imshow(binarizedResidual);
   %  figure(4);
   %  imshow(relevantElements);
    
    % Calculate number of true positives.
    selectedElements = binarizedClusters;
    %relevantElements = binarizedModel.*binarizedResidual;
    truePositiveElements = times(relevantElements, selectedElements);
    truePositives = sum(sum(times(relevantElements, selectedElements)));
    disp(['     truePositives: ' int2str(truePositives)])

    % Calculate F1 confidence.
    precision = truePositives / sum(sum(selectedElements));
    recall = truePositives / sum(sum(relevantElements));
    disp(['     precision: ' num2str(precision)]);
    disp(['     recall: ' num2str(recall)]);
    f1confidence = (2 * precision * recall) / (precision + recall);
    disp(['     f1confidence: ' num2str(f1confidence)]);
    gxyParams.fitQualityF1 = f1confidence;

    pearsonCorrelationCoefficient = corr2(residual, selectedElements);
    gxyParams.fitQualityPCC = pearsonCorrelationCoefficient;
    result = gxyParams;

    % Writing images
    grouped =  cat(2,model,residual,relevantElements,selectedElements,truePositiveElements);
    %imwrite(model, [outputPathGname '-L1_model.png']);
    %imwrite(residual, [outputPathGname '-L2_residual.png']);
    %imwrite(relevantElements, [outputPathGname '-L3_maskedResidual.png']);
    %imwrite(selectedElements, [outputPathGname '-L4_clusMask.png']);
    %imwrite(truePositiveElements, [outputPathGname '-L5_maskedClustersResidual.png']);
    %imwrite(grouped, [outputPathGname '-L_fitQuality.png']);

    % Cleanup files - let's save it for the end shall we? Matthew 1/12/21
    %disp('Removing galfit.* files...')
    %delete([outputPathGname 'galfit.*']);
    %delete(['galfit.*']); 
    %delete([outputPathGname '.feedme']);
    %delete([outputPathGname '_galfit_input.fits']);
    %delete([outputPathGname '_galfit_output.fits']);
end

