function [params, thBounds] = emRefineArcs(img, params, ctrR, ctrC, thBounds, aviName)
% experimental; not currently used

tStart = tic;

maxIter = 20;
nClus = size(params, 1);

error(nargchk(5, 6, nargin));

useTrace = true;
if nargin < 6
    aviName = [];
    useTrace = false;
end

if useTrace
	aviTrace = avifile([tempdir aviName], 'fps', 1, ...
        'compression', 'None', 'quality', 100, 'KeyFramePerSec', 1)
end

for iter = 1:1:maxIter
    tic; asgns = calcClusAsgns(img, params, ctrR, ctrC, thBounds); toc
    wts = asgns .* repmat(img, [1, 1, nClus]);
%     if useTrace
% %         figure('Position', [1, 1, 800, 800]); 
%         figure
%         imagesc(sum(wts, 3)); axis image; colormap gray
%         aviTrace = addframe(aviTrace, gcf);
%         close
%     end
    tic; [params, thBounds] = fitLogSpiralsToClusters(wts, ctrR, ctrC); toc
%     if useTrace
%         arcImg = displayLgspOverlay(img, params, ctrR, ctrC, thBounds);
% %         figure('Position', [1, 1, 800, 800]); 
%         figure
%         imagesc(arcImg); axis image; colormap gray
%         aviTrace = addframe(aviTrace, gcf); close
%     end
    if useTrace
        figure('Position', [1 1 1200 600]);
        subplot(1, 2, 1); imagesc(sum(wts, 3)); axis image; colormap gray
        arcImg = displayLgspOverlay(img, params, ctrR, ctrC, thBounds);
        subplot(1, 2, 2); imagesc(arcImg); axis image; colormap gray
        aviTrace = addframe(aviTrace, gcf); close
    end
end

if useTrace
    aviTrace = close(aviTrace);
end

displayLgspOverlay(img, params, ctrR, ctrC, thBounds);

fprintf('Time for all %d iterations:\n');
toc(tStart)

end