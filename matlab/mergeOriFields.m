function oriMerged = mergeOriFields(oriCoarse, oriFine)
% Merges two orientation fields at different resolutions, as described in 
% the PhD thesis "Inferring Galaxy Morphology Through Texture Analysis" 
% (K. Au 2006).
% The "fine" orientation field should have twice the resolution (in each
% dimension) as the "coarse" orientation field
% INPUTS:
%   oriCoarse: the orientation field generated from a lower-resolution
%       version of the image
%   oriFine: the orientation field generated from a higher-resolution
%       version of the image
% OUTPUTS:
%   oriMerged: result from merging the two orientation fields

% coarseStrengths = oriFieldStrengths(oriCoarse);
fineStrengths = oriFieldStrengths(oriFine);

[fineRows, fineCols] = size(fineStrengths);
if sum(mod([fineRows, fineCols], 2)) ~= 0
    error('matrix dimensions not divisible by 2: %s', ...
        mat2str([fineRows, fineCols]))
end
oriCoarseOld = oriCoarse;
oriCoarse = zeros(size(oriFine));
oriCoarse(:, :, 1) = imresize(oriCoarseOld(:, :, 1), 2);
oriCoarse(:, :, 2) = imresize(oriCoarseOld(:, :, 2), 2);
coarseStrengths = oriFieldStrengths(oriCoarse);

gains = fineStrengths ./ (coarseStrengths + fineStrengths);
gains(isnan(gains)) = 0;
oriMerged = orientationAdd(oriCoarse, repmat(gains, [1 1 2]) .* orientationSubtract(oriFine, oriCoarse));
% figure; imagesc(oriFieldStrengths(oriMerged)); axis image; title('merged');

% Adds a and b, interpreting the vectors as orientations. Unlike
% vector subtraction, orientations differing by 180 degrees are considered
% equivalent. a(i, j) and b(i, j) should be 2-element vectors.
function sum = orientationAdd(a, b)
    sum = orientationAddOrSubtract(a, b, true);
end

function sum = orientationSubtract(a, b)
    sum = orientationAddOrSubtract(a, b, false);
end

function res = orientationAddOrSubtract(a, b, add)
    neg_b2 = b(:, :, 2) < 0;
    neg_b2 = repmat(neg_b2, [1 1 2]);
    b(neg_b2) = -b(neg_b2);
    vec_sum = a + b;
    vec_diff = a - b;
    vec_sum_lengths = sqrt(sum(vec_sum.^2, 3));
    vec_diff_lengths = sqrt(sum(vec_diff.^2, 3));
    sumGreater = vec_sum_lengths > vec_diff_lengths;
    sumGreater = repmat(sumGreater, [1 1 2]);
    res = zeros(size(a));
    if add
        res(sumGreater) = vec_sum(sumGreater);
        res(~sumGreater) = vec_diff(~sumGreater);
    else
        res(~sumGreater) = vec_sum(~sumGreater);
        res(sumGreater) = vec_diff(sumGreater);
    end
end

end