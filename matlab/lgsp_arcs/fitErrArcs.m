function errArcs = ...
    fitErrArcs(clusMtx, ctrR, ctrC, bndryTh, errMult, bestParams, allowArcBeyond2pi, theta, plotFlag)

% wtSt - standard deviation to use if using direction weighting (empty for
% no weighting)

error(nargchk(7, 9, nargin));

bestIrOpts = optimset('TolFun', 1e-4, 'TolCon', 1e-4, 'TolX', 1e-4);

if nargin < 9 || isempty(plotFlag)
    plotFlag = false;
end

if allowArcBeyond2pi && (nargin < 8 || isempty(theta))
    error('theta-values must be specified for 2rev fits');
end

if length(size(clusMtx)) ~= 2
    error('clusMtx must be a 2D matrix');
end

if ~isvector(bndryTh) || numel(bndryTh) ~= 2
    error('bndryTh must be 2-element vector');
end

if size(bndryTh, 1) == 1
    bndryTh = bndryTh';
end

if allowArcBeyond2pi
    lgspFxn = @logSpiralFxn2Rev;
else
    lgspFxn = @logSpiralFxn;
end

[r, c, brt] = find(clusMtx);
r = r - ctrR; c = c - ctrC;
[th, rho] = cart2pol(c, -r);
if ~allowArcBeyond2pi
    theta = mod(th, 2*pi);
end

ignorePts = (brt < 0.001);  % save time by ignoring points with very low weight
brt = brt(~ignorePts);
rho = rho(~ignorePts);
theta = theta(~ignorePts);

thOff = bestParams(1);

% totBrt = sum(brt)

% bestSse = sum( (sqrt(brt) .* (rho - lgspFxn(theta, bestParams))).^2 );
bestSse = sum(lgspErr(lgspFxn, bestParams, theta, rho, brt).^2);
% bestMseRecalc = bestMseRecalc / totBrt

% optFxn = @(prms)(lgspFxn(bndryTh(1), [thOff prms]));
% optFxn = @(prms)(diff(lgspFxn(bndryTh, [thOff prms])));
% optFxn = @(prms)(sum(lgspFxn(theta, [thOff prms])));
optFxn = @(prms)(prms(2));
% optFxn = @(prms)(prms(1));

minIr = 1;
maxIr = 2 * size(clusMtx, 1);
maxIr = 4 * size(clusMtx, 1);
% maxIr = size(clusMtx, 1);
% lBounds = [-inf -maxIr]; 
lBounds = [-inf minIr]; 
uBounds = [inf maxIr];

% opts = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'Display', 'notify', 'Algorithm', 'interior-point');

opts = optimset('MaxFunEvals', 1000, 'Display', 'notify', 'Algorithm', 'active-set');
opts = optimset('MaxFunEvals', 500, 'MaxIter', 500, 'TolFun', 1e-4, 'TolX', 1e-4, 'Display', 'notify', 'Algorithm', 'active-set');
% opts = optimset('MaxFunEvals', 500, 'MaxIter', 500, 'TolFun', 10e-4, 'TolX', 10e-4, 'Display', 'notify');

% bdry1par = fmincon(optFxn, bestParams(2:end), ...
%     [], [], [], [], lBounds, uBounds, ...
%     @errcon, opts); %, optimset('Algorithm', 'interior-point', 'AlwaysHonorConstraints', 'bounds'))
% 
% bdry2par = fmincon(@(prms)(-optFxn(prms)), bestParams(2:end), ...
%     [], [], [], [], lBounds, uBounds, ...
%     @errcon, opts);
% 
% bdry1par = [thOff bdry1par];
% bdry2par = [thOff bdry2par];

maxPa = pi;
minPa = -maxPa;
% maxPa = [];
% minPa = [];
[bdry1pang, fval1, outFlag1] = fmincon(@(x)(x), bestParams(2), [], [], [], [], minPa, maxPa, @errconPa, opts);
% fprintf('*******\n');
% fprintf('start = %2.4f\n', bestParams(2));
[bdry2pang, fval2, outFlag2] = fmincon(@(x)(-x), bestParams(2), [], [], [], [], minPa, maxPa, @errconPa, opts);
bdry1par = [thOff bdry1pang bestIr(bdry1pang)];
bdry2par = [thOff bdry2pang bestIr(bdry2pang)];
if outFlag1 < 0 || outFlag2 < 0 || abs(bdry1pang) > pi || abs(bdry2pang) > pi
    outFlag1 < 0
    outFlag2 < 0
    abs(bdry1pang) > pi
    abs(bdry2pang) > pi
    assert(false);
end


    function [c, ceq] = errcon(iterParams)
        if isnan(iterParams)
            error('NaN encountered');
        end
%         err = sum( (sqrt(brt) .* (rho - lgspFxn(theta, [thOff iterParams]))).^2 );
        err = sum(lgspErr(lgspFxn, [thOff iterParams], theta, rho, brt).^2);
%         err = err / totBrt
%         errMult * bestSse
        c = err - (errMult * bestSse);  % must have err <= (errMult*bestMse)
        ceq = 0; % no equality constraints
        
%         displayLgspOverlay(clusMtx, [thOff iterParams], ctrR, ctrC, [-2*pi 2*pi]);
%         title(sprintf('err = %2.4f, params = %s', err, mat2str(iterParams, 8)));
%         title(sprintf('err = %2.4f, vals = %s', err, mat2str(logSpiralFxn2Rev(bndryTh, [thOff iterParams]), 4)));
    end

    function [c, ceq] = errconPa(pa)
        [ir, err] = bestIr(pa);
%         displayLgspPlot([thOff pa ir], bndryTh', clusMtx, 128, 128);
%         fprintf('pa = %2.4f\n', pa);
        
        c = err - (errMult * bestSse);  % must have err <= (errMult*bestMse)
%         fprintf('c = %2.4f = (%2.4f - %2.4f)\n', c, err, (errMult * bestSse));
        ceq = 0; % no equality constraints
    end

    function [ir, err] = bestIr(pa)
        lfv = lgspFxn(theta, [thOff pa 1]);
        ir = sum(rho .* lfv .* brt) / sum(brt .* (lfv .^ 2));
        err = sum(lgspErr(lgspFxn, [thOff pa ir], theta, rho, brt).^2);
%         [ir, err] = fminbnd(@(x)(sum(lgspErr(lgspFxn, [thOff pa x], theta, rho, brt).^2)), minIr, maxIr, bestIrOpts)
    end

if plotFlag
    bndExt = [-pi/2 pi/2];
%     bndExt = [-2*pi, 2*pi];
    displayLgspPlot([bestParams; bdry1par; bdry2par], repmat(bndryTh' + bndExt, 3, 1),...
        clusMtx, 128, 128, [], [], {}, {'-', ':', ':'}, 2);
end

errArcs = [bdry1par; bdry2par];

% if ~isempty(wtSd)
%     brtWtdL = brt .* normpdf(theta, min(theta), wtSd);
%     brtWtdL = brtWtdL / max(brtWtdL(:));
%     clusMtxWtdL = zeros(size(clusMtx));
%     clusMtxWtdL(clusMtx > 0) = brtWtdL;
%     errArcsL = fitErrArcs(clusMtxWtdL, ctrR, ctrC, bndryTh, errMult, ...
%         bestParams, allowArcBeyond2pi, wtSd, theta, plotFlag);
% end

end