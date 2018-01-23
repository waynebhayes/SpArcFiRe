function agreeingRows = ...
    getAgreeingRows(tbl, fld1, fld2, matches, hasHeader)

error(nargchk(5, 5, nargin));

if hasHeader
    header = tbl(1, :);
    tbl = tbl(2:end, :);
end

isAgreeingRow = false(size(tbl, 1), 1);
for ii=1:1:size(matches, 1)
    isAgreeingRow = isAgreeingRow | ...
        ( strcmp(matches{ii, 1}, tbl(:, fld1)) & ...
        strcmp(matches{ii, 2}, tbl(:, fld2)) );
end

agreeingRows = tbl(isAgreeingRow, :);

if hasHeader
    agreeingRows = [header; agreeingRows];
end

end