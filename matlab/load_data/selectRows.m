function tbl = selectRows(tbl, fxns, fldIdxs, hasHeader)

if hasHeader
    header = tbl(1, :);
    tbl = tbl(2:end, :);
end

isSelected = false(size(tbl, 1), 1);
for ii=1:1:length(fxns)
    values = str2double(tbl(:, fldIdxs(ii)));
    isSelected = isSelected | arrayfun(fxns{ii}, values);
end

tbl = tbl(isSelected, :);

if hasHeader
    tbl = [header; tbl];
end

end