function display_arc_plot(params, bounds, arcFxn, img, ctrR, ctrC, barAngles, barHalfLens, lineColors, lineStyles, lineWidth, makeFig)
% Displays an image overlayed with plots of parameterized arcs and bars
% INPUTS:
%   params: parameters of the log-spirals, one per row
%   bounds: starting and ending theta-range of points to which
%       log-spirals were fitted, counterclockwise from the theta-offset
%       parameter.  Each row is a [start, end] pair.
%   arcFxn: handle for a function of the form 
%       [x, y] = arcFxn(params, bounds)
%   img: the image to overlay
%   ctrR: the row position of the center of the log-spirals
%   ctrC: the column position of the center of the log-spirals
%   barAngles: the angles (in radians) of any bars, empty if no bar
%   barHalfLens: the half-lengths of any bars, empty if no bar
%   lineStyles:
%   makeFig: whether to put the plot in a new figure, rather than using the
%       current one default true

error(nargchk(5, 11, nargin));

nArcs = size(params, 1);

if nargin < 7 || isempty(barAngles)
    barHalfLens = [];
end

if nargin < 8 || isempty(barHalfLens)
    barHalfLens = [];
end

if nargin < 9
    lineColors = {};
elseif ~isempty(lineColors)
    if size(lineColors, 1) < nArcs
        lineColors = repmat(lineColors, ceil(nArcs/size(lineColors, 1)), 1);
        lineColors = lineColors(1:nArcs, :);
    end
end

if nargin < 10 || isempty(lineStyles)
    lineStyles = {'-'};
else
    validateattributes(lineStyles, {'cell'}, {});
    if length(lineStyles) < nArcs
        lineStyles = repmat(lineStyles(:), ceil(nArcs/length(lineStyles)), 1);
        lineStyles = lineStyles(1:nArcs);
    end
end

if nargin < 11 || isempty(lineWidth)
    lineWidth = 3;
    lineWidth = 2;
    lineWidth = 1;
end

if nargin < 12 || isempty(makeFig)
    makeFig = true;
end

if size(lineColors, 2) > 1
    if size(lineColors, 1) > 1
        error('lineColors should be a vector');
    end
    lineColors = lineColors';
end

if length(lineStyles) == 1
    lineStyles = repmat(lineStyles, 1, size(params, 1));
end

% ctrX = ctrC;
% ctrY = size(img, 1) - ctrR + 1;

if makeFig
    fh = figure; 
else
    fh = gcf;
end
figure(fh); axh = gca();
hold on
hold all
% if isnumeric(lineColors)
%     set(fh,'DefaultAxesColorOrder',lineColors)
% end
clims = [0 1];
% clims = [0 max(img(:))];
if makeFig
    if length(size(img)) > 2
        imshow(flipdim(img, 1), 'XData', [1 size(img, 2)] - ctrC, 'YData', [1 size(img, 1)] - ctrR);
    else
        imagesc([1 size(img, 2)] - ctrC, [1 size(img, 1)] - ctrR, flipud(img), clims);
    end
end
for arcIdx=1:1:size(params, 1);
    if ~isempty(lineColors) && isnumeric(lineColors{1})
        % line colors already set
        lineType = lineStyles{arcIdx};
    elseif ~isempty(lineColors)
        lineType = [lineColors{arcIdx} lineStyles{arcIdx}];
    elseif params(arcIdx, 2) > 0
        lineType = ['r' lineStyles{arcIdx}];
    else
        lineType = ['c' lineStyles{arcIdx}];
    end

    [x, y] = arcFxn(params(arcIdx, :), bounds(arcIdx, :));
%     if isnumeric(lineColors)
    if ~isempty(lineColors)
        plot(axh, x, y, lineType, 'Color', lineColors{arcIdx, :}, 'LineWidth', lineWidth);
    else
        plot(axh, x, y, lineType, 'LineWidth', lineWidth);
    end
    
%     plot([0 100*cos(params(arcIdx, 1))], [0 100*sin(params(arcIdx, 1))], lineType, 'LineWidth', lineWidth);
end
for bar = 1:1:length(barHalfLens)
    angle = barAngles(bar);
    halfLength = barHalfLens(bar);
    rho = -halfLength:0.5:halfLength;
    th = angle * ones(1, length(rho));
    [x, y] = pol2cart(th, rho);
    plot(x, y, 'g-', 'LineWidth', lineWidth);
end
axis image
% if makeFig
    colormap gray
% end
hold off

% figure;
% polar(th, r);
% axis image

% [x, y] = pol2cart(th, r);
% figure;
% plot(x, y)
% axis image

end