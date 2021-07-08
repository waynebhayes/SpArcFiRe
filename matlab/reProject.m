function [ repI ] = reProject(I, diskAxisRatio, diskMajorAxisLength, diskMajorAxisAngleRadians, inputCenter, iptSz)

% Documentation added by Matthew P. --- 5/15/21
% 
% This function calculates a reprojection of the image using galaxy parameters.
% If image standardization is turned off, the function returns the original
% image back so that SpArcFiRe may proceed. Note, this seems to return an error
% for parameter covarFit during runtime but it seems to output to the csv correctly.

% -------- INPUT --------
% inputCenter = [inputCenterR, inputCenterC], other arguments are same as
% them in galaxy.csv
% -----------------------
%
% -------- OUTPUT -------
% With image standardization:
%     reprojected image, repI
% Without:
%     original image, I (as repI)
% -----------------------

    fprintf('Calculating reprojection...\n');
   
    paddedI = padarray(I, [384 384]);
    ratio = diskAxisRatio;
    angle = pi/2-diskMajorAxisAngleRadians;
    scaleX = diskMajorAxisLength/size(I, 2);
    scaleY = diskMajorAxisLength/size(I, 1);
    
    % Added by Matthew P. --- 5/15/21
    % Checking if diskAxisRatio is empty (which happens when 
    % image standardization is turned off) to return an
    % unprojected image back to the calling function.
    if isempty(ratio)
        fprintf('Image standardization is turned off. No reprojection needed.\n')
        repI = I;
        return
    end
    % Note, this means there's a bug in the calling script
    % but this seems to solve the issue... for now.    
    
    m = affine2d([scaleY*ratio*cos(angle), scaleY*ratio*sin(angle), 0; 
        -scaleX*sin(angle), scaleX*cos(angle), 0; 
        0, 0, 1]);
    
    limits = imref2d([1024 1024], [-512 512], [-512 512]);
    %repI = imwarp(paddedI, m, 'UData', [-512 512], 'VData', [-512 512],...
    %                   'XData', [-512 512], 'YData', [-512 512]);
    repI = imwarp(paddedI, limits, m, 'OutputView', limits);
    
    leftBorder = round(512 - inputCenter(2));
    %leftBorder
    rightBorder = round(512 + iptSz(2) - inputCenter(2) - 1);
    %rightBorder
    topBorder = round(512 - inputCenter(1));
    %topBorder
    bottomBorder = round(512 + iptSz(1) - inputCenter(1) - 1);
    %bottomBorder
    repI = repI(topBorder:bottomBorder, leftBorder:rightBorder,:);    
    fprintf('...done calculating reprojection\n');
%     centerShiftXL = round(inputCenter(2)-size(repI,2)/2);
%     centerShiftXR = round(iptSz(2)-inputCenter(2)-size(repI,2)/2);
%     centerShiftYL = round(inputCenter(1)-size(repI,1)/2);
%     centerShiftYR = round(iptSz(1)-inputCenter(1)-size(repI,1)/2);
%     repI = cat(2, zeros(size(repI,1), centerShiftXL, size(repI,3)), repI,...
%                     zeros(size(repI,1), centerShiftXR, size(repI,3)));
%     repI = cat(1, zeros(centerShiftYL, size(repI,2), size(repI,3)), repI,...
%                     zeros(centerShiftYR, size(repI,2), size(repI,3)));
end

