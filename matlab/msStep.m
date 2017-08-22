function newMean = msStep(mtx, bandw, oldMean, wrapCols, traceFlag)

if nargin < 3
    error(['first 3 arguments (matrix, bandwidth, old mean)' ...
        'must be specified']);
end

if nargin < 4 || isempty(wrapCols)
    wrapCols = false;
end

if nargin < 5 || isempty(traceFlag)
    traceFlag = false;
end

% get a bounding box around the bandwidth circle, to avoid doing
% distance calculations that we know we won't need
minRow = round(oldMean(1)-bandw)-0;
maxRow = round(oldMean(1)+bandw)+0;
minCol = round(oldMean(2)-bandw)-0;
maxCol = round(oldMean(2)+bandw)+0;
rStart = max(1,minRow);
cStart = max(1,minCol);
nearby = ...
    mtx(rStart:min(end,maxRow), cStart:min(end,maxCol));

% add wraparound columns if needed
mtxCols = size(mtx, 2);
colWrapDiff = 0;
if wrapCols && minCol < 1
    colWrapDiff = minCol - 1;
    nearby = [mtx(rStart:min(end,maxRow), end+colWrapDiff:end) nearby];
elseif wrapCols && maxCol > mtxCols
    colWrapDiff = maxCol - mtxCols;
    nearby = [nearby mtx(rStart:min(end,maxRow), 1:colWrapDiff)];
end

% nearby = mtx; rStart = 1; cStart = 1; % disable the speed trick to see if it is harming accuracy
[rPts, cPts] = find(nearby > 0);

% shift the window indices to be the same as the matrix indices
rPts = rPts + rStart - 1;
cPts = cPts + cStart - 1;
if colWrapDiff < 0
    % we added wrapped columns on the left, compensate for the change in
    % indexing
    cPts = cPts + colWrapDiff;
end
pts = [rPts, cPts];
wts = nearby(nearby > 0); % weight by the number of votes

% sqDists = sum((repmat(curMean',1,numPts)' - pts).^2, 2);
sqDists = distToAll(pts, round(oldMean)); % ......
distWts = (sqDists <= bandw^2) .* wts;
% distWts = (1./(sqDists + 1)) .* wts;
if traceFlag
    nearby'
    reshape(wts, size(nearby))'
    reshape(distWts, size(nearby))'
    reshape(rPts, size(nearby))'
    reshape(cPts, size(nearby))'
    ([sum(rPts .* distWts), sum(cPts .* distWts)] / sum(distWts)) - oldMean
    [sum((rPts - round(oldMean(1))) .* distWts), sum((cPts - round(oldMean(2))) .* distWts)] / sum(distWts)
end

rDiffs = rPts - round(oldMean(1));
cDiffs = cPts - round(oldMean(2)); 
% rDiffs = rPts - oldMean(1);
% cDiffs = cPts - oldMean(2); 
newMean = oldMean + ...
    ( [sum(rDiffs .* distWts), sum(cDiffs .* distWts)] / sum(distWts) );
if wrapCols
    % if a column index went out of bounds, wrap it around, keeping in mind
    % that means will eventually be rounded to the nearest index (integer)
    newMean(2) = mod(newMean(2)-1, mtxCols - 0.5) + 1;
end

end