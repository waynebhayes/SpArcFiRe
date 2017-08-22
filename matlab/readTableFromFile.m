function [tbl, nSkip, col1ofSkipped, wasSkipped] = readTableFromFile(fname, fldDelim)
% Reads an ASCII table and returns a cell array with each entry
% corresponding to a field in the input table.  All rows in the table must
% have the same number of columns.

if nargin < 2 || isempty(fldDelim)
    fldDelim = ' \b\n\r\t';
end

fprintf('reading table: %s\n', fname);

fid = fopen(fname, 'rt');
lines = textscan(fid, '%[^\n]', 'BufSize', 10000);
fclose(fid);

lines = lines{1};
header = lines(1);
% lines = lines(2:end);

headerFlds = textscan(header{1}, '%s', 'Delimiter', fldDelim);
headerFlds = headerFlds{1};
nFlds = length(headerFlds);
fprintf('number of fields: %d\n', nFlds);

tbl = cell(length(lines), length(headerFlds));
nSkip = 0;
wasSkipped = false(length(lines), 1);
col1ofSkipped = cell(length(lines), 1);
for ii=1:1:length(lines)
    curFlds = textscan(lines{ii}, '%s', 'Delimiter', fldDelim);
    curFlds = curFlds{1};
    if nFlds ~= length(curFlds)
        curFlds
        fprintf([lines{ii} '\n'])
        warning(['line %d of file has wrong number of fields (expected %d, ',...
            'was %d). line skipped.'], ii, nFlds, length(curFlds))
        nSkip = nSkip + 1;
        wasSkipped(ii) = true;
        col1ofSkipped{ii} = curFlds{1};
        continue;
    end
    tbl(ii, :) = curFlds;
end
fprintf('%d lines skipped.\n', nSkip);
col1ofSkipped = col1ofSkipped(wasSkipped);

fprintf('found %d rows and %d columns\n', size(tbl));

tbl = tbl(~wasSkipped, :);

end