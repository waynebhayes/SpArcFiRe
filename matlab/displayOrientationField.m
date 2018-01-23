function fh = displayOrientationField(ofld, removeZeroElts, visible)

if nargin < 3 || isempty(visible)
    visible = true;
end
if visible
    visible_str = 'on';
else
    visible_str = 'off';
end

isNz = oriFieldStrengths(ofld) > 0;
if removeZeroElts
    [nzr, nzc, nzv] = find(isNz);
    nzi = sub2ind(size(isNz), nzr, nzc);
    ofld1 = ofld(:, :, 1);
    ofld2 = -ofld(:, :, 2);
    fh = figure('visible', visible_str); 
    quiver(nzc, nzr, ofld1(nzi), ofld2(nzi), 'MaxHeadSize', 0);
else
    fh = figure('visible', visible_str); 
    quiver(ofld(:, :, 1), -ofld(:, :, 2), 'MaxHeadSize', 0, 'LineWidth', 2); 
end

axis image; axis ij
axis([1 size(isNz, 2) 1 size(isNz, 1)]);

end