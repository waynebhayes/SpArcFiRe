function region = indexRegion(arr, minIdxs, maxIdxs, wrapDims, truncate)
% Extracts a region of the input array equivalent to
%   arr(minIdxs(1):maxIdxs(1), ... , minIdxs(end):maxIdxs(end))
% INPUTS:
%   arr: the array to index from
%   minIdxs: the minimum indices along each dimension
%   maxIdxs: the maximum indices along each dimension
%   wrapDims: whether to wrap indices if they are out of range, expressed
%       either as a vector (specifying this setting per dimension) or as a
%       scalar (specifying this setting for all dimensions)
%   truncate: whether to ignore out of range indices if wrapDims is not
%       set, expressed as either a vector or a scalar
% OUTPUTS:
%   region: the input array section specified by the given indices


if nargin < 4
    wrapDims = false;
end

if nargin < 5
    truncate = false;
end

ndim = length(size(arr));

if ~isvector(minIdxs) || ~isvector(maxIdxs)
    error('minIdxs and maxIdxs need to be vectors')
end

if ~isvector(wrapDims)
    error(['wrapDims must be a vector or scalar specifying whether ' ...
        'to wrap indices']);
end

if ~isa(wrapDims, 'logical')
    error('wrapDims needs to be of type "logical"')
end

if length(minIdxs) ~= ndim || length(maxIdxs) ~= ndim
    error(['minIdxs and maxIdxs need to have one element per array ' ...
        'dimension'])
end

if length(wrapDims) == 1
    wrapDims = wrapDims .* true(1, length(minIdxs));
end

if length(wrapDims) ~= ndim
    error(['wrapDims needs to either be a scalar, or have one element ' ...
        'per array dimension'])
end
    
if ~isscalar(truncate) || ~isa(truncate, 'logical')
    error('input "truncate" needs to be a scalar logical')
end

idxs = cell(1, ndim);
for dim=1:1:ndim
    dimIdxs = minIdxs(dim):maxIdxs(dim);
    dimLen = size(arr, dim);
    if wrapDims(dim)
        dimIdxs = mod(dimIdxs - 1, dimLen) + 1;
    elseif truncate
        dimIdxs = dimIdxs(dimIdxs > 0 & dimIdxs <= dimLen);
    elseif sum(dimIdxs < 1 | dimIdxs > dimLen)
        error(['out-of-range index for non-wrapped dimension. ' ...
            'to ignore and omit these indices, set truncate=true.'])
    end
    idxs{dim} = dimIdxs;
end

region = arr(idxs{:});

end