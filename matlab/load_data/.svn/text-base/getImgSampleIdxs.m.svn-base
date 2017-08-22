function [sampleIdxs, eligIdxs] = ...
    getImgSampleIdxs(imgInfo, selCriteria, sampleSize)
% Produces a random sample of indices of images meeting the given criteria
% Sample use: 
%   sampleIdxs = getImgSample(tbl, {49, @(type)(type >= 1 && type <= 7); 
%       55, @(pairs)(pairs == 0)}, 15);
% INPUTS:
%   imgInfo: 2D cell array of strings where the first row contains 
%       attribute descriptors, the following rows correspond to images and 
%       columns correspond to attributes (as produced by readTableFromFile)
%   selCriteria: n-by-2 cell array describing criteria images must meet in
%       order to be eligible for inclusion in the sample.  The first column
%       of selCriteria is the column index of an attribute in the imgInfo 
%       table; the second column of selCriteria is a function of the
%       numerical value of that attribute that returns true if the image
%       remains eligible for inclusion after checking that attribute
%   sampleSize: the sample size.  If empty or larger than the number of 
%       images meeting the given criteria, all eligible images will be 
%       returned
% OUTPUTS
%   sampleIdxs: indices of the images in the sample, where indices are the
%       images' positions in the imgInfo table, starting at 1 and not
%       including the table's header line
%   eligIdxs: indices of the images that were eligible for the sample

    error(nargchk(1, 3, nargin));
    if nargin < 2
        selCriteria = {};
    end
    if nargin < 3
        sampleSize = [];
    end

    imgInfo = imgInfo(2:end, :);
    nImgs = size(imgInfo, 1);
    
    isEligible = true(nImgs, 1);
    for critIdx=1:1:size(selCriteria, 1)
        isEligible = isEligible .* arrayfun(selCriteria{critIdx, 2}, ...
            str2double(imgInfo(:, selCriteria{critIdx, 1})) );
    end
    eligIdxs = find(isEligible);
    
    if ~isempty(sampleSize)
        sampleSize = min(sampleSize, length(eligIdxs));
        sampleIdxs = randsample(eligIdxs, sampleSize);
    else
        sampleIdxs = eligIdxs;
    end
end