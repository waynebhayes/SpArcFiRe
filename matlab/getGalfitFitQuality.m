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

    
    %  
    fitswrite(img, [outputPath '_galfit_input.fits']);
    galfitTemplate = fopen('~/sparcfire/galfit/template.feedme','r');
    text = fread(galfitTemplate);
    fclose(galfitTemplate);

    % Fill out template by replacing variables with their actual values
    muFit = gxyParams.iptCtrXY;
    nRows = gxyParams.iptSz(1);
    nCols = gxyParams.iptSz(2);
    text = strrep(text, '$input_name', [outputPath '_galfit_input.fits']);
    text = strrep(text, '$output_name', [outputPath '_galfit_output.fits']);
    text = strrep(text, '$x_center', num2str(muFit(1)));
    text = strrep(text, '$y_center', num2str(size(img, 1) - muFit(2) + 1));
    text = strrep(text, '$radius_bulge', num2str(gxyParams.bulgeMajAxsLen / 4));
    text = strrep(text, '$radius_disk', num2str(gxyParams.diskMajAxsLen / 4));
    text = strrep(text, '$sersic_index', num2str(1.0));
    text = strrep(text, '$axis_ratio', num2str(gxyParams.diskAxisRatio));
    text = strrep(text, '$position_angle', num2str(90 - rad2deg(gxyParams.diskMajAxsAngleRadians)));
    text = strrep(text, '$x_max', num2str(nRows));
    text = strrep(text, '$y_max', num2str(nCols));

    galfitInput = fopen([outputPath '.feedme'],'wt');
    fwrite(galfitInput, text);
    fclose(galfitInput);

    % Run galfit
    galfitCommand = ['!~/sparcfire/galfit/galfit ' outputPath '.feedme'];
    eval(galfitCommand)

    % Retrieve input/model subtraction
    model = mat2gray(fitsread([outputPath '_galfit_output.fits'], 'image', 2));
    residual = mat2gray(fitsread([outputPath '_galfit_output.fits'], 'image', 3));
    binarizedClusters = imbinarize(rgb2gray(clusReproj), 0.1);
    binarizedModel = imbinarize(model);
    binarizedResidual = imbinarize(residual);
    relevantElements = binarizedModel.*binarizedResidual;
%     figure(1);
%     imshow(binarizedClusters);
%     figure(2);
%     imshow(binarizedModel);
%     figure(3);
%     imshow(binarizedResidual);
%     figure(4);
%     imshow(relevantElements);
    
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
    
    gxyParams.fitQuality = f1confidence;
    result = gxyParams;



    % Writing images
    grouped =  cat(2,model,residual,relevantElements,selectedElements,truePositiveElements);
    imwrite(model, [outputPath '-L1_model.png']);
    imwrite(residual, [outputPath '-L2_residual.png']);
    imwrite(relevantElements, [outputPath '-L3_maskedResidual.png']);
    imwrite(selectedElements, [outputPath '-L4_clusMask.png']);
    imwrite(truePositiveElements, [outputPath '-L5_maskedClustersResidual.png']);
    imwrite(grouped, [outputPath '-L_fitQuality.png']);

    % Cleanup files
    %disp('Removing galfit.* files...')
    %delete([outputPath 'galfit.*']);
    %delete([outputPath '.feedme']);
    %delete([outputPath '_galfit_input.fits']);
    %delete([outputPath '_galfit_output.fits']);
end

