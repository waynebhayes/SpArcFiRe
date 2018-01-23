function [ht, dVals, thVals] = houghTransform(img, oris, ...
    dGran, thGran, xOff, yOff, wtByIntensities, plotFlag)
% Computes the Hough transform of an image
% INPUTS:
%   img: grayscale image for which the Hough transform is to be computed
%   oris: optional parameter specifying the orientations at each point 
%       (e.g., from Canny edge detection).  If these orientations are 
%       given, points will give higher weights to lines with slopes similar
%       to the orientation at the point.
%   dGran: optional parameter specifying the spacing between d-values in
%       the Hough accumulator matrix.  Defaults to 2 and may be adjusted so
%       that the range of d-values covers the image exactly and uniformly.
%   thGran: optional parameter specifying the spacing between theta-values
%       in the Hough accumulator matrix.  Defaults to pi/360 and may be
%       adjusted so that the range of theta-values covers the interval 
%       [0, pi] exactly and uniformly.
%   xOff: optional parameter specifying the offset to be added to the
%       x-positions inferred from the pixel positions
%   yOff: optional parameter specifying the offset to be added to the
%       y-positions inferred from the pixel positions
%   wtByIntensities: optional parameter specifying whether points with
%       higher intensities give stronger votes.  Defaults to true.
%   plotFlag: 
% OUTPUTS:
%   ht: Hough transform output matrix.  Rows correspond to d-values, which
%       range from 0 to the Euclidean length of the image diagonal,
%       measured in pixels.  Columns correspond to theta-values, which
%       range from 0 to pi.  The value of any element corresponds to the
%       vote total for the presence of a line with the (d, theta) parameter
%       values corresponding to the element indices.  d-values are measured
%       from the bottom left of the image.
%   dVals: d-Values used in the Hough transform, so that dVals(i) is the
%       d-value of the i^th row in the Hough transform matrix. The
%       granularity may have been adjusted (vs the given granularity) to
%       cover the image exactly and uniformly.
%   thVals: theta-values used in the Hough transform, so that thVals(i) is
%       the theta-value of the i^th column in the Hough transform matrix.
%       The granularity may have been adjusted (vs the given granularity)
%       to cover the range of theta-values exactly and uniformly.

if nargin < 1
    error('img needs to be specified')
end

if nargin < 2 || isempty(oris)
    use_oris = false;
else
    use_oris = true;
end

if nargin < 3 || isempty(dGran)
    dGran = 2;
end

if nargin < 4 || isempty(thGran)
    thGran = pi/(2 * 180);
end

if nargin < 5 || isempty(xOff)
    xOff = 0;
end

if nargin < 6 || isempty(yOff)
    yOff = 0;
end

if nargin < 7 || isempty(wtByIntensities)
    wtByIntensities = true;
end

if nargin < 8 || isempty(plotFlag)
    plotFlag = false;
end

if length(size(img)) ~= 2
    error('img needs to be a m x n matrix')
end

if ~isa(img, 'double')
    img = im2double(img);
end

% image indexing to Cartesian indexing for the rows
img = flipud(img);

% increase granularity a bit if needed to make the bins fit exactly
dMax = ceil(sqrt(sum(size(img) .^2)));
dBins = ceil((2 * dMax) / dGran) + 1;
dGran = (2 * dMax) / (dBins - 1);
dVals = [0:1:dBins-1] * dGran - dMax;

thBins = ceil(pi/thGran);
thGran = pi / thBins;
thVals = [0:1:thBins-1] * thGran;

[nzY, nzX] = find(img);
nzNum = length(nzY);

ht = zeros(dBins, thBins);

if plotFlag
    fprintf('{dMax, dBins, thBins} = {%d, %d, %d}\n', dMax, dBins, thBins);
    fprintf('%d nonzero edge values\n', nzNum);
end

voteSizes = ones(size(img));
if wtByIntensities
    voteSizes = voteSizes .* img;
end

for eIdx=1:1:nzNum
    x = nzX(eIdx);
    y = nzY(eIdx);
    for thIdx=1:1:thBins
        th = thVals(thIdx);
        d = (x + xOff) * cos(th) + (y + yOff) * sin(th);
%         dIdx = uint16(round((d - dVals(1)) / dGran) + 1);
        dPos = (d - dVals(1)) / dGran;
        dIdxL = uint16(floor(dPos) + 1);
        dIdxH = uint16(ceil(dPos) + 1);
        lWt = 1 - (dPos - floor(dPos));
        hWt = 1 - (ceil(dPos) - dPos);
        if use_oris
            vote_amt = voteSizes(y, x) * abs(cos(th - oris(y, x)));
        else
            vote_amt = voteSizes(y, x);
        end
%         ht(dIdx, thIdx) = ht(dIdx, thIdx) + vote_amt;
        if (dIdxL < 0) || (dIdxL > size(ht, 1))
            hWt = hWt + lWt;
            lWt = 0;
        end
        if (dIdxH < 0) || (dIdxH > size(ht, 1))
            lWt = lWt + hWt;
            hWt = 0;
        end
        if lWt > 0
            ht(dIdxL, thIdx) = ht(dIdxL, thIdx) + lWt * vote_amt; 
        end
        if hWt > 0
            ht(dIdxH, thIdx) = ht(dIdxH, thIdx) + hWt * vote_amt;
        end
    end
end

end % ht