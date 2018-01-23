function tbl = joinTables(tbl1, tbl2, idx1, idx2)
% Joins two tables on a column that gives 1:1 correspondence between table
% rows.
% INPUTS:
%   tbl1, tbl2: the tables to join, as 2D cell arrays with the first row
%       containing titles
%   idx1, idx2: the columns to use for joining the tables.  These columns
%       must be equal except possibly for the column titles and element
%       ordering
% OUTPUTS:
%   tbl: the joined table

rowstart = 2;

tbl1(rowstart:end, :) = sortrows(tbl1(rowstart:end, :), idx1);
tbl2(rowstart:end, :) = sortrows(tbl2(rowstart:end, :), idx2);

[isEq, neqIdx] = ...
    allEqual(tbl1(rowstart:end, idx1), tbl2(rowstart:end, idx2));
if ~isEq
    if neqIdx < 0
        neqReason = 'row counts different';
    else
        neqReason = sprintf('different at sorted index %d', neqIdx);
    end
    error('1:1 correspondence not found for given columns (%s)', neqReason)
end

fprintf('joining tables on columns "%s" and "%s"\n', ...
    tbl1{1, idx1}, tbl2{1, idx2});

t1cols = true(1, size(tbl1, 2));
t1cols(idx1) = false;
t2cols = true(1, size(tbl2, 2));
t2cols(idx2) = false;

tbl = [tbl1(:, idx1), tbl1(:, t1cols), tbl2(:, t2cols)];

    function [isEq, neqIdx] = allEqual(clarr1, clarr2)
        isEq = true;
        neqIdx = -1;
        if length(clarr1) ~= length(clarr2)
            isEq = false;
            neqIdx = -1;
            return;
        end
        for ii=1:1:length(clarr1)
            if ~strcmp(clarr1{ii}, clarr2{ii})
                isEq = false;
                neqIdx = ii;
                return;
            end
        end
    end
end