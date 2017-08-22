function [tbl1, tbl2, nDrop1, nDrop2] = intersectRowsByField(tbl1, tbl2, idx1, idx2, hasHeader)

if hasHeader
    firstRow1 = tbl1(1, :); tbl1 = tbl1(2:end, :);
    firstRow2 = tbl2(1, :); tbl2 = tbl2(2:end, :);
end

% TODO: put this all in one pass if this takes too long (check correctness
% if this change made)

keys1 = tbl1(:, idx1);
keys2 = tbl2(:, idx2);
% save('keys', 'keys1', 'keys2');

% if length(unique(keys1)) < length(keys1)
%     error('non-unique elements in keys from tbl1');
% elseif length(unique(keys2)) < length(keys2)
%     error('non-unique elements in keys from tbl2');
% end

% commonElts = intersect(keys1, keys2);

inIntersect1 = ismember(keys1, keys2);
inIntersect2 = ismember(keys2, keys1);

fprintf('table 1 entries not in intersection: \n')
keys1(~inIntersect1)
fprintf('table 2 entries not in intersection: \n')
keys2(~inIntersect2)
tbl1 = [firstRow1 ; tbl1(inIntersect1, :)];
tbl2 = [firstRow2 ; tbl2(inIntersect2, :)];
nDrop1 = sum(~inIntersect1);
nDrop2 = sum(~inIntersect2);

end