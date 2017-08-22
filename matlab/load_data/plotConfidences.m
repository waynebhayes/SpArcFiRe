function plotConfidences(tbl, fld1, fld2, hasHeader)

if hasHeader
    tbl = tbl(2:end, :);
end

if isscalar(fld1)
    conf1 = cellfun(@str2num, tbl(:, fld1), 'UniformOutput', false);
    conf1 = cell2mat(conf1);
else
    conf1 = tbl(:, ismember(1:1:size(tbl, 1), fld1));
    conf1 = cellfun(@str2double, conf1);
end

% normalize votes/probabilities/confidences
conf1 = conf1 ./ repmat(sum(conf1, 2), 1, 2);
conf1 = max(conf1, [], 2);

if isscalar(fld2)
    conf2 = cellfun(@str2num, tbl(:, fld2), 'UniformOutput', false);
    conf2 = cell2mat(conf2);
else
    conf2 = tbl(:, ismember(1:1:size(tbl, 1), fld2));
    conf2 = cellfun(@str2double, conf2);
end
% normalize votes/probabilities/confidences
conf2 = conf2 ./ repmat(sum(conf2, 2), 1, 2);
conf2 = max(conf2, [], 2);

scatter(conf1, conf2);

end