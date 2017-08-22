function img = broaden_arc_width(img, arcColor, newColor, tol)
% Makes plotted arcs thicker, and optionally changes their color, for 
% better visibility in presented materials.

error(nargchk(2, 4, nargin))

imgcls = class(img);
img = double(img);

if nargin < 3 || isempty(newColor)
    newColor = arcColor;
end

if nargin < 4 || isempty(tol)
    tol = 0;
end

if length(arcColor) ~= 3 || length(newColor) ~= 3
    error('arcColor and newColor must have length 3')
end

arcRegion = true(size(img, 1), size(img, 2));
for ii=1:1:3
    arcRegion = arcRegion & (abs(img(:, :, ii) - arcColor(ii)) <= tol);
end
if sum(find(arcRegion)) == 0
    fprintf('warning: no pixels match the arc color.\n')
end
arcRegion = imdilate(arcRegion, strel('square', 3));

for ii=1:1:3
    imgD = img(:, :, ii);
    imgD(arcRegion) = newColor(ii);
    img(:, :, ii) = imgD;
end

img = cast(img, imgcls);

end