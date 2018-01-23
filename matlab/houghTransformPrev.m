function [ht, dVals, thVals] = houghTransform(img, oris, ...
    dGran, thGran, wtByIntensities, ssEdges)
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
%   wtByIntensities: optional parameter specifying whether points with
%       higher intensities give stronger votes.  Defaults to true.
%   ssEdges: optional parameter specifying edges for surround-suppression
%       vote weighting (see calcSurrSuppWts for details).  If unset,
%       surround suppression will not be used.
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

if nargin < 5 || isempty(wtByIntensities)
    wtByIntensities = true;
end

if nargin < 6
    ssEdges = [];
end
if length(ssEdges) == 1 && isa(ssEdges, 'logical')
    warning('ssEdges is a logical, should be matrix or empty');
    ssEdges = [];
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
% thVals = [0:thNum-1] * (pi/thNum);
thVals = [0:1:thBins-1] * thGran;

fprintf('{dMax, dBins, thBins} = {%d, %d, %d}\n', dMax, dBins, thBins);

[nzY, nzX] = find(img);
nzNum = length(nzY);

% tPts = zeros(length(nzY) * thNum, 2); %temp
% tPts = [0 0];
ht = zeros(dBins, thBins);

fprintf('%d nonzero edge values\n', nzNum);

voteSizes = ones(size(img));
if wtByIntensities
    voteSizes = voteSizes .* img;
end
if ~isempty(ssEdges);
    voteSizes = voteSizes .* calcSurrSuppWts(ssEdges, 2, 2);
end
figure; imagesc(flipud(voteSizes)); title('vote sizes');

for eIdx=1:1:nzNum
    x = nzX(eIdx);
    y = nzY(eIdx);
    for thIdx=1:1:thBins
        th = thVals(thIdx);
        d = x * cos(th) + y * sin(th);
        dIdx = uint16(round((d - dVals(1)) / dGran) + 1);
        % tPts(thBins * (eIdx - 1) + thIdx, :) = [dVals(dIdx), th];
        if use_oris
            vote_amt = voteSizes(y, x) * abs(cos(th - oris(y, x)));
        else
            vote_amt = voteSizes(y, x);
        end
        ht(dIdx, thIdx) = ht(dIdx, thIdx) + vote_amt;
    end
end

% figure; scatter(tPts(:, 1), tPts(:, 2))

% FIXME (hack)
% thWts = repmat(abs(mod(thVals + pi/4, pi/2) - pi/4), [dBins, 1]);
% ht = ht .* thWts;
% ht = nonmaxSuppress(ht, 1, 5);

% figure; imagesc(ht); title('hough transform'); axis image; axis off;
% fprintf('\nunprocessed:\n');
% printSlopeScores(ht);
% 
% mthres_ht = ht .* ( ht > mean(ht(:)) );
% figure; imagesc(mthres_ht); 
% title('mean-thresholded SHT'); 
% axis image; axis off;
% fprintf('\nmean-thresholded:\n');
% printSlopeScores(mthres_ht);
% 
% qval = .9;
% medthres_ht = ht .* (ht > quantile(ht(:), qval));
% figure; imagesc(medthres_ht);
% title(sprintf('%1.4f quantile thresholded SHT', qval));
% axis image; axis off;
% fprintf('\n%1.4f quantile thresholded SHT:\n', qval)
% printSlopeScores(medthres_ht);
% 
% qval = .95;
% medthres_ht = ht .* (ht > quantile(ht(:), qval));
% figure; imagesc(medthres_ht);
% title(sprintf('%1.4f quantile thresholded SHT', qval));
% axis image; axis off;
% fprintf('\n%1.4f quantile thresholded SHT:\n', qval)
% printSlopeScores(medthres_ht);
% 
% qval = .99;
% medthres_ht = ht .* (ht > quantile(ht(:), qval));
% figure; imagesc(medthres_ht);
% title(sprintf('%1.4f quantile thresholded SHT', qval));
% axis image; axis off;
% fprintf('\n%1.4f quantile thresholded SHT:\n', qval)
% printSlopeScores(medthres_ht);
% 
% qval = .999;
% medthres_ht = ht .* (ht > quantile(ht(:), qval));
% figure; imagesc(medthres_ht);
% title(sprintf('%1.4f quantile thresholded SHT', qval));
% axis image; axis off;
% fprintf('\n%1.4f quantile thresholded SHT:\n', qval)
% printSlopeScores(medthres_ht);
% 
% sorted_ht_vals = sort(ht(:), 'descend');
% % sorted_ht_vals(1:10)
% 
% ht_top_20 = ht .* (ht > sorted_ht_vals(21));
% figure; imagesc(ht_top_20);
% title('top 20 points');
% axis image; axis off;
% fprintf('\ntop-20:\n')
% printSlopeScores(ht_top_20);
% 
% ht_top_100 = ht .* (ht > sorted_ht_vals(101));
% figure; imagesc(ht_top_100);
% title('top 100 points');
% axis image; axis off;
% fprintf('\ntop-100:\n')
% printSlopeScores(ht_top_100);
% 
% thres = numel(ht) * 0.0001 * mean(ht(:))
% ht_thres = ht .* (ht > thres);
% figure; imagesc(ht_thres);
% title('points with .01%%+ mean-votes')
% axis image; axis off;
% fprintf('\n .01%% threshold:\n');
% printSlopeScores(ht_thres);
% 
% % figure; plot(thVals)
% % figure; plot(mod(thVals, pi/2))
% % figure; plot(abs(mod(thVals, pi/2) - pi/4))
% % figure; plot(abs(mod(thVals + pi/4, pi/2) - pi/4))
% imagesc(log(repmat(abs(mod(thVals + pi/4, pi/2) - pi/4), [dBins, 1])));
% 
% % mean_nz = sum(ht(:)) / sum(ht(:) > 0);
% % sfactor = 2; fprintf('\nmean-threshold scale factor: %f\n', sfactor);
% % mthres_ht = ht .* (ht > (sfactor * mean_nz));
% % figure; imagesc(mthres_ht); 
% % title(sprintf('mean-thresholded SHT (%d X mean)', sfactor)); 
% % axis image; axis off;
% % printSlopeScores(mthres_ht);
% % 
% % mean_nz = sum(ht(:)) / sum(ht(:) > 0);
% % sfactor = 5; fprintf('\nmean-threshold scale factor: %f\n', sfactor);
% % mthres_ht = ht .* (ht > (sfactor * mean_nz));
% % figure; imagesc(mthres_ht); 
% % title(sprintf('mean-thresholded SHT (%d X mean)', sfactor)); 
% % axis image; axis off;
% % printSlopeScores(mthres_ht);
% 
% % nms_ht = nonmaxSuppress(ht, 2, 1);
% % figure; imagesc(nms_ht); 
% % title('nonmax-suppressed SHT'); axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT:\n');
% % printSlopeScores(nms_ht);
% % 
% % nms_ht = nonmaxSuppress(ht, 2, 2);
% % figure; imagesc(nms_ht); 
% % title('nonmax-suppressed SHT (neighborhood size 2)'); axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (neighborhood size 2):\n');
% % printSlopeScores(nms_ht);
% % 
% % nms_ht = nonmaxSuppress(ht, 2, 5);
% % figure; imagesc(nms_ht); 
% % title('nonmax-suppressed SHT (neighborhood size 5)'); axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (neighborhood size 5):\n');
% % printSlopeScores(nms_ht);
% % 
% % nms_ht = nonmaxSuppress(ht, 2, 10);
% % figure; imagesc(nms_ht); 
% % title('nonmax-suppressed SHT (neighborhood size 10)'); 
% % axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (neighborhood size 10):\n');
% % printSlopeScores(nms_ht);
% % 
% % nms_ht = nonmaxSuppress(ht, 2, 20);
% % figure; imagesc(nms_ht); 
% % title('nonmax-suppressed SHT (neighborhood size 20)'); 
% % axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (neighborhood size 20):\n');
% % printSlopeScores(nms_ht);
% 
% 
% 
% 
% % nms_htb = nonmaxSuppress(conv2(ht, genGaussMtx(1), 'same'), 2);
% % figure; imagesc(nms_htb); 
% % title('nonmax-suppressed SHT (r1 blur)'); axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (r1 blur):\n');
% % printSlopeScores(nms_htb);
% % 
% % nms_htb = nonmaxSuppress(conv2(ht, genGaussMtx(2), 'same'), 2);
% % figure; imagesc(nms_htb); 
% % title('nonmax-suppressed SHT (r2 blur)'); axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (r2 blur):\n');
% % printSlopeScores(nms_htb);
% % 
% % nms_htb = nonmaxSuppress(conv2(ht, genGaussMtx(3), 'same'), 2);
% % figure; imagesc(nms_htb); 
% % title('nonmax-suppressed SHT (r3 blur)'); axis image; axis off;
% % fprintf('\nnonmax-suppressed SHT (r3 blur):\n');
% % printSlopeScores(nms_htb);

end % ht