function xcMtx = analyzeXCorr(b1pol, b2pol, nTheta, name)

if nargin < 4
    name = [];
end

if isempty(name)
    figVis = 'on';
else
    figVis = 'off';
end

nPadReps = 20;
radTickPcts = [0 0.25 0.5 0.75 1];

if length(size(b1pol)) ~= 2 || length(size(b2pol)) ~= 2
    error('b1pol and b2pol must be 2D arrays');
end
 
if sum(size(b1pol) == size(b2pol)) ~= 2
    error('b1pol and b2pol must have the same dimensions');
end

% size(b1pol, 2) = one full revolution
maxLag = ceil(size(b1pol, 2)/2);
maxLag = ceil(size(b1pol, 2)/8);
maxLag = 30;
xcMtx = zeros(2*maxLag+1, size(b1pol, 1));
for radIdx=size(b1pol, 1):-1:1
    a1 = b1pol(radIdx, :); a2 = b2pol(radIdx, :);
    [xc, lags] = ...
        xcorr(repmat(a1, 1, 1+2*nPadReps), repmat(a2, 1, 1+2*nPadReps), ...
        maxLag, 'coeff');
    xcMtx(:, size(b1pol, 1)-radIdx+1) = xc;
end
[min(lags) max(lags)]

figure('Visible', figVis); imshow(xcMtx); set(gca,'YDir','normal')

% figure; hold on; 
% imagesc(xcMtx); axis image; colormap gray
% % contour(flipud(xcMtx), 3, '-g');
% lbls = get(gca, 'YTickLabel');
% set(gca, 'YTickLabel', lags(str2double(mat2cell(lbls, ones(size(lbls, 1), 1), size(lbls, 2)))));

% figure; hold on; 
% % xcMtx(:, lags == 0) = 0;
% imagesc(xcMtx); axis image; colormap gray
% line([0 size(xcMtx, 2)], find(lags==0) * [1 1])
% xlabel('rIdx');
% ylabel('lag');
% lbls = get(gca, 'YTickLabel');
% set(gca, 'YTickLabel', lags(str2double(mat2cell(lbls, ones(size(lbls, 1),
% 1), size(lbls, 2)))));

% findXcorrMax(xcMtx, lags);

analyzeGaussFit(xcMtx, lags, nTheta, name);

end