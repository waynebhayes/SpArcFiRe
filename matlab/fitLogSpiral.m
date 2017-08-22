function [fitParams, arcBounds, sumSqErr, used2Rev, failed2rev, hasBadBounds, fitTh, ptErrs, errArcs, fitWts] = ...
    fitLogSpiral(img, ctrR, ctrC, stgs, tipWtDctn, initParams, initBounds, lgspErrFxn, plotFlag, fixedPa)
% Fits a 3-parameter logarithmic spiral model to the data points given via
% an input image.  
% INPUTS: 
%   img: image giving the data points as nonzero pixels
%   ctrR: row of the polar origin for the log spiral
%   ctrC: column of the polar origin for the log spiral
%   stgs: structure containing algorithm settings (see settings.m)
%   tipWtDctn: 
%   initParams: optional guess (starting parameters) for the log spiral
%   initBounds: (needed if tip weighting used)
%   plotFlag: whether to plot the fit results (optional, default false)
%   fixedPa: if set, gives the fit error for the given pitch angle
%       instead of finding the optimal one
% OUTPUTS:
%   fitParams: 3-element row vector of the fit parameters 
%       [theta-offset, pitch-angle, initial-radius]
%   arcBounds: start and ending angle of the data points, counterclockwise
%       from the theta-offset value
%   sumSqErr: sum of squared errors for the arc fit
%   fitTh: the theta-values used for fitting (avoid using for now (TODO:
%       get this working by getting rid of rotAmt, if possible))
%   used2Rev
%   ptErrs: 
%   errArcs
%   fitWts

% err = MSE w/o normalization by weight sum

allowArcBeyond2pi = stgs.allowArcBeyond2pi;
paApproxLevel = stgs.paApproxLevel;
zeroThetaStart = stgs.zeroThetaStart;
maxHalfGapFillForUndefinedBounds = stgs.maxHalfGapFillForUndefinedBounds;

imgSz = size(img);

error(nargchk(4, 10, nargin))

if nargin < 5 || isempty(tipWtDctn);
    tipWtDctn = [];
end

% figure; imshow(img); title(sprintf('input image (tipWtDctn = %s)', mat2str(tipWtDctn, 2)));

if ~isempty(tipWtDctn) && ~isscalar(tipWtDctn)
    error('tipWtDctn must be a scalar');
end

if nargin < 7 || isempty(initParams) || isempty(initBounds)
    if ~isempty(tipWtDctn)
        error(['must specify initial log-spiral parameters and bounds if using '...
            'tip-weighting']);
    end
end

if nargin < 6 || isempty(initParams)
    initParams = [0, 0, 10];
end
if length(initParams) ~= 3
    error('length of initParams must be 3 (thOff, pitchAngle, initRad)');
end

if nargin < 7
    initBounds = [];
end

if nargin < 8 || isempty(lgspErrFxn)
    lgspErrFxn = @lgspErr;
end

if nargin < 9 || isempty(plotFlag)
    plotFlag = false;
end

if nargin < 10
    fixedPa = [];
end

fitErrArcs = (nargout >= 9);

lsqOpts = optimset('MaxFunEvals', 1000*length(initParams), 'TolFun', 1e-4, 'TolX', 1e-4, 'Display', 'off',...
    'DerivativeCheck', 'on', 'Jacobian', 'off');
% lsqOpts = optimset(lsqOpts, 'Algorithm', 'levenberg-marquardt'); %Octave-test

[rows, cols, brt] = find(img);
[theta, rho] = cart2pol(cols - ctrC, -(rows - ctrR));
theta = mod(theta, 2*pi);
hasBadBounds = false; % true if no well-defined theta-bounds can be determined
failed2rev = false; % true if 2rev fit attempted but a sgl-rev region could not be found

if isempty(brt)
    warning('empty cluster');
    fitParams = [NaN NaN NaN];
    arcBounds = [NaN NaN];
    % the worst possible error, so that this cluster will never be selected
    % over a nonempty one
    sumSqErr = inf; 
    fitTh = [];
    used2Rev = [];
    ptErrs = [];
    errArcs = [];
    fitWts = [];
    return;
end

if allowArcBeyond2pi
    [lBounds, uBounds] = calculateBounds(theta);
    ctrInClus = (img(round(ctrR), round(ctrC)) ~= 0);
   
    if ~isempty(lBounds) 
        % theta-gap large enough to get bounds (on the theta-offset) so
        % we don't need more than one theta-revolution
        allowArcBeyond2pi = false;
    elseif ctrInClus || clusterHasNoEndpoints(img)
%         fprintf('Warning:fitLogSpiral:cluster includes the center; not using 2-revolution fit\n');
        allowArcBeyond2pi = false;
    else  
        [isInner, gapFail] = idInnerOuterSpiral(img, ctrR, ctrC, plotFlag);
        nInner = sum(isInner);
        failed2rev = gapFail;
        if gapFail || (nInner == 0) || (nInner == numel(isInner))
            allowArcBeyond2pi = false;
%             lBounds = [];
%             uBounds = [];
        else
            thAdj = removeThetaDiscontFor2rev(theta, img, isInner);
            % TODO: make theta and thAdj one variable after making sure that
            % bounds calculation doesn't need a different (unaltered) value
            theta = thAdj; 
        end
    end
end

% To avoid a discontinuity of the log-spiral curve in the region where we
% are fitting, we want to restrict the range of theta-offset values to
% exclude theta-values within the range of the data points.
% We thus need to find the range of polar coordinate theta-values occupied
% by the data points, which we do by looking for the largest gap between
% data points.  This is safe if the data points are generated from HAC
% clustering (and are thus close to other points) and the data points
% occupy a range of theta-values that is less than 2*pi

if allowArcBeyond2pi
    lgspFxn = @logSpiralFxn2Rev;
else
    lgspFxn = @logSpiralFxn;
end

if allowArcBeyond2pi && isempty(fixedPa)
    % fitting depends on the arc bounds, but we don't know what the arc
    % bounds are (i.e., whether to calculate the bounds from the inner or
    % outer points) until we know the chirality.  Since we only know the
    % chirality after the fit, we do two fits, one assuming CW and the
    % other assuming CCW, and take the one with the better error.
    
    % For the >2*pi case, we don't do mod 2*pi in the rotations because
    % theta-values can be outside the range [0, 2*pi], and discontinuities
    % in theta-values can severely impact the fitting.
    
    [lBoundsCW, uBoundsCW, rotAmtCW] = calculateBounds(theta(~isInner));
    theta = theta - rotAmtCW;
    lBoundsCW(2) = 0; uBoundsCW(2) = inf;
    thOffCW = (lBoundsCW(1) + uBoundsCW(1)) / 2;
    thOff = thOffCW;
    if fitErrArcs
        pixWts = calcTipWtsSgl(brt, theta - thOff, initParams, initBounds, tipWtDctn, imgSz); 
    end
    [fitParamsCW, sumSqErrCW] = lsqnonlin(@errFxnPa,...
        initParams(2), lBoundsCW(2), uBoundsCW(2), lsqOpts);
    ptErrsCW = errFxnPa(fitParamsCW).^2;
    fitParamsCW = [fitParamsCW bestIr(fitParamsCW, thOff)];  
    theta = theta + rotAmtCW;
    fitParamsCW = [thOff fitParamsCW];
    
    [lBoundsCCW, uBoundsCCW, rotAmtCCW] = calculateBounds(theta(isInner));
    theta = theta - rotAmtCCW;
    lBoundsCCW(2) = -inf; uBoundsCCW(2) = 0;
    thOffCCW = (lBoundsCCW(1) + uBoundsCCW(1)) / 2;
    thOff = thOffCCW;
    if fitErrArcs
        pixWts = calcTipWtsSgl(brt, theta - thOff, initParams, initBounds, tipWtDctn, imgSz); 
    end
    [fitParamsCCW, sumSqErrCCW] = lsqnonlin(@errFxnPa,...
    	initParams(2), lBoundsCCW(2), uBoundsCCW(2), lsqOpts);
    ptErrsCCW = errFxnPa(fitParamsCCW).^2;
    fitParamsCCW = [fitParamsCCW bestIr(fitParamsCCW, thOff)];
    theta = theta + rotAmtCCW;
    fitParamsCCW = [thOff fitParamsCCW];
    
    if sumSqErrCW < sumSqErrCCW
        fitParams = fitParamsCW;
        sumSqErr = sumSqErrCW;
        theta = theta - rotAmtCW;
        rotAmt = rotAmtCW;
        thOff = thOffCW;
        ptErrs = ptErrsCW;
    else
        fitParams = fitParamsCCW;
        sumSqErr = sumSqErrCCW;
        theta = theta - rotAmtCCW;
        rotAmt = rotAmtCCW;
        thOff = thOffCCW;
        ptErrs = ptErrsCCW;
    end
else
    [lBounds, uBounds, rotAmt, maxGapSize] = calculateBounds(theta);
    theta = mod(theta - rotAmt, 2*pi);
    if ~isempty(lBounds)
        thOff = (lBounds(1) + uBounds(1)) / 2;
    else
        % no good split (entire theta-range covered), and we're not
        % allowing more than one theta-revolution), but we have to set the
        % theta-offset to something
%         fprintf('Warning:fitLogSpiral:low maxGapSize (%2.4f)\n', maxGapSize);
        [minVal, minIdx] = min(rho);
        thOff = theta(minIdx);
    end
%     pixWts = calcPixWts(mod(theta - thOff, 2*pi), initParams);
    if fitErrArcs
        pixWts = calcTipWtsSgl(brt, mod(theta - thOff, 2*pi), initParams, initBounds, tipWtDctn, imgSz); 
    end
    if ~isempty(fixedPa)
        fitParams = fixedPa;
        sumSqErr = sum(errFxnPa(fixedPa).^2);
        ptErrs = errFxnPa(fitParams).^2;
    elseif isempty(lBounds)
%         fprintf('fitting without theta-gap\n');
%         minThOff = 0;
%         maxThOff = 2*pi;
%         initPa = lsqnonlin(@errFxnPa, initParams(2), [], [], lsqOpts);
%         initParams(1) = (minThOff + maxThOff)/2;
%         initParams(2) = initPa;
%         initParams(1:2)
%         [fitParams, sumSqErr] = lsqnonlin(@errFxnThoffPa, ...
%             initParams(1:2), [minThOff, -inf], [maxThOff, inf], lsqOpts);
%         ptErrs = errFxnThoffPa(fitParams).^2;
        hasBadBounds = true;
%         [maxVal, maxIdx] = max(rho);
%         initParams(1) = theta(maxIdx);
%         [fitParams, sumSqErr] = lsqnonlin(@errFxnLgspFull, initParams,...
%             [], [], lsqOpts);
%         ptErrs = errFxnLgspFull(fitParams).^2;
        bad_bounds_pa = 0;
        bad_bounds_thOff = 0;
        [ir, err] = bestIr(bad_bounds_pa, bad_bounds_thOff);
        ptErrs = err.^2;
        sumSqErr = sum(ptErrs);
        fitParams = [bad_bounds_thOff, bad_bounds_pa, ir];
    else
        if paApproxLevel > 0
            eigPa = calcPitchAngle(theta, rho, brt);
        end

        if paApproxLevel < 1
            [fitParams, sumSqErr] = lsqnonlin(@errFxnPa, ...
                initParams(2), [], [], lsqOpts);
        elseif paApproxLevel < 2
            [fitParams, sumSqErr] = lsqnonlin(@errFxnPa, ...
                eigPa, [], [], lsqOpts);
        else
            fitParams = eigPa;
        end
        ptErrs = errFxnPa(fitParams).^2;
        if paApproxLevel >= 2
            sumSqErr = sum(ptErrs);
        end
        
%         fprintf('fitted pa = %2.4f, atan(fitted pa) = %2.4f, eig pa = %2.4f, pa + eig pa = %2.4f, pa - eig pa = %2.4f\n',...
%             fitParams * (180/pi), atan(fitParams) * (180/pi), eigPa *
%             (180/pi), (fitParams + eigPa) * (180/pi), (fitParams - eigPa) * (180/pi));
    end
    if ~isempty(fixedPa) || ~isempty(lBounds)
        fitParams = [fitParams bestIr(fitParams, thOff)];
        fitParams = [thOff fitParams];
    end
    theta = mod(theta + rotAmt, 2*pi); % TODO: do we really need this mod operation?
end
assert(abs(sum(ptErrs) - sumSqErr) / numel(ptErrs) < 10^-4);

% adjust the theta-offset parameter to account for the rotation (which was
% only done so that the allowable bounds for the theta-offest could be
% expressed as a single interval)
if allowArcBeyond2pi
    % why do we need mod 2*pi?  Is it for the condition 
    % fitParams(1) > min(thAdj) ?
    fitParams(1) = mod(fitParams(1) + rotAmt, 2*pi);
else
    fitParams(1) = fitParams(1) + rotAmt;
end

% find the polar-coordinate theta-range occupied by the data points
if allowArcBeyond2pi
    % why is thAdj used to calculate the bounds?  We won't need this if we
    % get rid of rotAmt, since thAdj and theta differ by rotAmt
    if fitParams(1) > min(thAdj)
        arcBounds = [min(thAdj) max(thAdj)] - fitParams(1);
    else
        arcBounds = [min(thAdj) max(thAdj)] - (fitParams(1) + 2*pi);
    end
else
    arcSize = 2*pi - (ub - lb);
    arcStart = min(mod(angleDist(fitParams(1), [lb, ub]) + rotAmt, 2*pi));
    % note: lb, ub aren't really meaningful if the points extend beyond 2pi -
    % these bounds were for the optimizer
    arcBounds = [arcStart arcStart + arcSize];
end

if zeroThetaStart
    if ~allowArcBeyond2pi
        theta = mod(theta - fitParams(1), 2*pi) + fitParams(1);
    end
    thOffShift = min(arcBounds);
    newThOff = fitParams(1) + thOffShift;
    theta = theta + (newThOff - min(theta));
    if hasBadBounds
        [ir, err] = bestIr(0, newThOff);
    else
        [ir, err] = bestIr(fitParams(2), newThOff);
    end
    fitParams(1) = newThOff;
    fitParams(3) = ir;
    arcBounds = arcBounds - thOffShift;
    newSumSqErr = sum(err.^2);
    sqErrDiffPerPixel = abs(newSumSqErr - sumSqErr) / length(theta);
    if (sqErrDiffPerPixel > 1e-10)
        error('internal error: inconsistent fit when eliminating theta_start (diff = %g)\n', sqErrDiffPerPixel);
    end
    sumSqErr = newSumSqErr;
end

if nargout >= 7  % fitTh
    if allowArcBeyond2pi
%         fitTh = thAdj;
        fitTh = theta;
%         figure; hold all; plot(sort(theta)); plot(sort(thAdj)); title(sprintf('theta vs thAdj (diff = %2.4f, thOff = %2.4f, rotAmt = %2.4f)', min(theta) - min(thAdj), fitParams(1), rotAmt));
    else
%         thvImg = img; thvImg(thvImg > 0) = theta; figure; imagesc(thvImg); axis image; title('theta'); impixelinfo
%         thvImg = img; thvImg(thvImg > 0) = theta - thOff; figure; imagesc(thvImg); axis image; title('theta - thOff'); impixelinfo
        thOff = fitParams(1);
%         thOffA = mod(thOff, 2*pi);
        fitTh = mod((theta - thOff), 2*pi) + thOff;
%         thvImg = img; thvImg(thvImg > 0) = fitTh; figure; imagesc(thvImg); axis image; title('fitTh'); impixelinfo
    end
    % use theta (not thAdj) for both?
    
end

if nargout >= 9 % errArcs
    % errMult used to be specified in settings, before extended merging 
    % was removed. If error arcs are brought back, errMult should be
    % an input parameter.
    errMult = 2; 
    if allowArcBeyond2pi
        thForErrArcs = thAdj;
    else
        thForErrArcs = [];
    end
    wtImg = zeros(size(img));
    wtImg(img > 0) = pixWts;
    errArcs = fitErrArcs(wtImg, ctrR, ctrC, arcBounds, errMult, fitParams, ...
        allowArcBeyond2pi, thForErrArcs, plotFlag);
end

if nargout >= 10 % fitWts
%     fitWtImg = zeros(size(img));
%     fitWtImg(img > 0) = pixWts;
    fitWts = pixWts;
end

if plotFlag
%     figure; hold on; axis equal; scatter(cols, -rows, 'g.');
    figure; hold on; axis equal; scatter(cols - ctrC, ctrR - rows, 'g.');
    if allowArcBeyond2pi
        plThetaRange = arcBounds + fitParams(1);
        plTheta = [plThetaRange(1):pi/360:plThetaRange(2)]';
        polar(plTheta, logSpiralFxn2Rev(plTheta, fitParams));
    else
        plTheta = 0:pi/360:2*pi;
        polar(plTheta, logSpiralFxn(plTheta, fitParams));
    end
    polar(0, 0);
    thetaInBounds = [arcBounds(1):pi/360:arcBounds(2)] + fitParams(1);
    polar(thetaInBounds, 10 * ones(size(thetaInBounds)))
    polar(arcBounds(1) + [0 eps] + fitParams(1), [0 max(rho)])
    polar(arcBounds(2) + [0 eps] + fitParams(1), [0 max(rho)]) 
    title(sprintf(['[thOff, a, ir] = %s\n '...
        '[lb, ub] = %s, thOff = %2.4f, arcBounds = %s, TZ = %d'], ...
        mat2str(fitParams, 4), mat2str([lb, ub], 4), fitParams(1), ...
        mat2str(arcBounds, 4), gapIsThruZero));
    hold off
    
    displayLgspOverlay(img, fitParams, ctrR, ctrC, arcBounds);
end

used2Rev = allowArcBeyond2pi;

%     function wts = calcPixWts(theta, lgspParams)
%         if isempty(tipWtDctn) % || (tipWtDctn == 0)
%             wts = brt;
%         else
%             if tipWtDctn < 0
%                 tipEnd = min(theta);
%             else
%                 tipEnd = max(theta);
%             end
%             ptAlens = calcLgspArcLengths(...
%                 repmat(lgspParams, length(theta), 1),...
%                 [repmat(tipEnd, length(theta), 1), theta]);
%             
%             tipWtSigma = max(ptAlens)/4;
%             wts = normpdf(ptAlens, 0, tipWtSigma);
%             wts = wts / max(wts(:));
%             wts = brt .* wts + eps;
%             % add eps so all points still detected as cluster members by find()
%         end
%         if plotFlag
%             makeThetaPlot(theta, 'theta-values when constructing weights');
%         end
%     end

    function pa = calcPitchAngle(theta, rho, brt)
        if length(theta) > 2
            th_s = sort(theta);
            absdiff = abs(diff(th_s));
            zero_gap = 2*pi - (th_s(end) - th_s(1));
            if max(absdiff) > zero_gap
                theta = mod(theta, 2*pi);
                [mV, mI] = max(absdiff);
                theta = mod(theta - th_s(mI+1), 2*pi);
            end
        end
        
        logRho = log(rho + 1);
        wtdMean = [sum(logRho .* brt), sum(theta .* brt)] / sum(brt);
        ctrdRho = logRho - wtdMean(1);
        ctrdTh = theta - wtdMean(2);
        ctrdPts = [ctrdRho(:),ctrdTh(:)];
        wtdCov = (repmat(brt(:), 1, 2)' .* ctrdPts') * ctrdPts;
        assert(all(size(wtdCov) == [2,2]))
        [eVec, eVal] = eig(wtdCov);
        if eVal(1) > eVal(2)
            peVec = eVec(1, :);
        else
            peVec = eVec(2, :);
        end
        pa = atan(peVec(2)/peVec(1));
        
        err = sum(errFxnPa(pa).^2);
        errOpp = sum(errFxnPa(-pa).^2);
        if errOpp < err
            pa = -pa;
        end
    end

    function [err, jac] = errFxnPa(pa)
        [ir, err] = bestIr(pa, thOff);
        if nargout >= 2
            jac = jacPa(pa, ir);
            fprintf('ir = %2.4f\n', ir);
            fprintf('pa = %2.4f\n', pa);
            fprintf('jac = %s\n', mat2str(jac, 4));
            fprintf('%s - %s\n', mat2str(errFxnPa(pa + sqrt(eps)), 4), mat2str(errFxnPa(pa - sqrt(eps)), 4));
            fprintf('efd = %s\n', mat2str(errFxnPa(pa + sqrt(eps)) - errFxnPa(pa - sqrt(eps)) / (2*sqrt(eps)), 4));
        end
    end
    function err = errFxnThoffPa(params)
        [ir, err] = bestIr(params(2), params(1));
    end
    function err = errFxnLgspFull(params)
        err = lgspErrFxn(lgspFxn, params, theta, rho, brt);
    end
    function [ir, err] = bestIr(pa, thOff)
        lfv = lgspFxn(theta, [thOff pa 1]);
        ir = sum(rho .* lfv .* brt) / sum(brt .* (lfv .^ 2));
        err = lgspErrFxn(lgspFxn, [thOff pa ir], theta, rho, brt);
        
%         minIr = 1;
%         maxIr = 4 * size(img, 1);
%         bestIrOpts = optimset('TolFun', 1e-4, 'TolCon', 1e-4, 'TolX', 1e-4);
%         [ir, err] = fminbnd(@(x)(sum(lgspErrFxn(lgspFxn, [thOff pa x], theta, rho, brt).^2)), minIr, maxIr, bestIrOpts);
    end
    function jac = jacPa(pa, ir)
%         jac = 2 * sqrt(brt) .* exp(-pa.*(theta-thOff)).*(theta-thOff).*ir.*(rho-exp(-pa.*(theta-thOff)).* ir);
%         jac = sqrt(2 .* brt .* exp(-pa.*(theta-thOff)) .* (theta-thOff) .* ir .* (rho - exp(-pa .* (theta-thOff)) .* ir));
        jac = sqrt(brt) .* exp(-pa.*(theta-thOff)) .* (theta-thOff) .* ir;
    end

    function [lBounds, uBounds, rotAmt, maxGapSize] = calculateBounds(thForGap)
        if any(size(thForGap) == 0)
            figure; imshow(img);
        end
        thSortedForGap = sort(thForGap, 'ascend');
        
        gaps = diff(thSortedForGap);
        zGap = thSortedForGap(1) + 2*pi - thSortedForGap(end);
        gapIsThruZero = zGap > max(gaps);
%         maxGap = max([gaps; zGap]);

        % The optimization function lets us restrict theta-offset values by
        % specifying lower and upper bounds.  If the range of allowable 
        % values goes through zero, then this range gets split into two 
        % parts, which we can't express with a single pair of bounds.  In 
        % this case, we temporarily rotate the points to allievate this
        % problem, fit the log-spiral model to the set of rotated points, 
        % and then reverse the rotation on the fitted model.
        if ~gapIsThruZero
            rotAmt = 0;
            [maxGapSize, maxGapLoc] = max(gaps);
            lb = thSortedForGap(maxGapLoc);
            ub = thSortedForGap(maxGapLoc+1);
        else
            rotAmt = thSortedForGap(1);
            lb = mod(thSortedForGap(end) - rotAmt, 2*pi);
            ub = 2*pi;
            maxGapSize = ub - lb;
        end

        if maxGapSize > 0.1
            lBounds = [lb, -inf, -inf];
            uBounds = [ub, inf, inf];
        else
            lBounds = [];
            uBounds = [];
        end
%         assert(abs(maxGapSize - maxGap) < 1e-6)
    end

    function TF = clusterHasNoEndpoints(clusImg)
        % See if the cluster has actual spiral endpoints by seeing if it is
        % possible to "escape" from the center point to the image boundary,
        % considering non-cluster pixels as empty pixels.
        TF = false;
        inClus = (clusImg > 0);
        ctrPixR = round(ctrR);
        ctrPixC = round(ctrC);
        if inClus(ctrPixR, ctrPixC)
            % We don't need to find holes in this case
            TF = true;
        end
        se_size = 2 * maxHalfGapFillForUndefinedBounds + 1;
        isHoleOrClus = imfill(imclose(inClus, ones(se_size, se_size)), 'holes');
        if isHoleOrClus(ctrPixR, ctrPixC)
            TF = true;
        end
    end

    % for debugging
    function makeThetaPlot(thv, titleStr)
        if nargin < 2
            titleStr = '';
        end
    	thvImg = -10 * ones(size(img));
%         thvImg(sub2ind(size(img), rows + ctrR, cols + ctrC)) = thv;
        thvImg(img > 0) = thv;
        figure; imagesc(thvImg); axis image; title(titleStr);
    end

end