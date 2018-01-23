function [ht, thVals, aVals, irVals] = ...
    logSpiralHoughTransform(img, ctrX, ctrY, oris)

if nargin < 1
    error('img needs to be specified')
end

if nargin < 3
    error('x-center and y-center need to be specified')
end

if nargin < 4 || isempty(oris)
    use_oris = false;
else
    use_oris = true;
end

if length(size(img)) ~= 2
    error('img needs to be a m x n matrix')
end

if ~isa(img, 'double')
    img = im2double(img);
end

thGran = pi/16;

aMax = 10;
aGran = 0.5;

irMax = 100;
irGran = 2;

excludeRadius = 2;

% image indexing to Cartesian indexing for the rows
img = flipud(img);
ctrY = size(img, 2) - ctrY + 1;

% increase granularity a bit if needed to make the bins fit exactly
thBins = ceil((2*pi)/thGran);
thGran = (2*pi) / thBins;
thVals = [0:1:thBins-1] * thGran;

aBins = ceil((2 * aMax) / aGran) + 1;
aGran = (2 * aMax) / (aBins - 1);
aVals = [0:1:aBins-1] * aGran - aMax;

irBins = ceil((irMax)/irGran) + 1;
irGran = (irMax) / (irBins - 1);
irVals = [0:1:irBins-1] * irGran;

fprintf('{aMax, irMax} = {%d, %d}\n', aMax, irMax);
fprintf('{thBins, aBins, irBins} = {%d, %d, %d}\n', thBins, aBins, irBins);

[nzY, nzX] = find(img);
nzNum = length(nzY);
fprintf('%d nonzero edge values\n', nzNum);

ht = zeros(thBins, aBins, irBins);

voteSizes = img;

for pIdx=1:1:nzNum
    x = nzX(pIdx) - ctrX;
    y = nzY(pIdx) - ctrY;
    [ptTh, ptR] = cart2pol(x, y);
    ptTh = mod(ptTh, 2*pi);
    if ptR <= excludeRadius
        continue;
    end
    for thIdx=1:1:thBins
        thShift = thVals(thIdx);
        for aIdx=1:1:aBins 
            a = aVals(aIdx);
            ir = ptR * exp(a * (ptTh - thShift));
            % assign vote to index corresponding to closest value for 
            % initial radius
            irIdx = round(ir / irGran) + 1;
            % the initial radius may be beyond what we're considering
            if irIdx <= irBins 
                vote_amt = voteSizes(nzY(pIdx), nzX(pIdx));
                ht(thIdx, aIdx, irIdx) = ht(thIdx, aIdx, irIdx) + vote_amt;
            end
        end
    end
end

end % logSpiralHoughTransform