function pim = polartransRGB(im, nrad, ntheta, cx, cy, linlog, shape)

if nargin < 4 || isempty(cx)
    cx = size(im, 2) / 2 + 0.5;
end

if nargin < 5 || isempty(cy)
    cy = size(im, 1) / 2 + 0.5;
end

if nargin < 6 || isempty(linlog)
    linlog = 'linear';
end

if nargin < 7 || isempty(shape)
    shape = 'full';
end

if length(size(im)) == 2
    pim = polartrans(im, nrad, ntheta, cx, cy, linlog, shape);
    return
end

pim = zeros([nrad ntheta 3]);
pim(:, :, 1) = ...
    polartrans(im(:, :, 1), nrad, ntheta, cx, cy, linlog, shape);
pim(:, :, 2) = ...
    polartrans(im(:, :, 2), nrad, ntheta, cx, cy, linlog, shape);
pim(:, :, 3) = ...
    polartrans(im(:, :, 3), nrad, ntheta, cx, cy, linlog, shape);

end