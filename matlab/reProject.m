function [ repI ] = reProject(I, diskAxisRatio, diskMajorAxisLength, diskMajorAxisAngleRadians, inputCenter, iptSz)
%inputCenter = [inputCenterR, inputCenterC], other arguments are same as
%them in galaxy.csv
    fprintf('Calculating reprojection...\n');
    paddedI = padarray(I, [384 384]);
    ratio = diskAxisRatio;
    angle = pi/2-diskMajorAxisAngleRadians;
    scaleX = diskMajorAxisLength/size(I, 2);
    scaleY = diskMajorAxisLength/size(I, 1);
    
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

