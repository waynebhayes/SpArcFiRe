function compareArmCounts(tbl, autoCountInds, autoThresholds, naCountInds, inclCrit)
% inclCriteria: struct array of idx and criFxn

includeInstance = true(size(tbl, 1) - 1, 1);
for ii=1:1:length(inclCrit)
    includeInstance = includeInstance & arrayfun(inclCrit(ii).criFxn, ...
        cellfun(@str2double, tbl(2:end, inclCrit(ii).idx)));
end
fprintf('%d of %d instances rejected due to inclusion criteria\n',...
    sum(~includeInstance), numel(includeInstance));
includeInstance = [false; includeInstance]; % skip header

autoCounts = cellfun(@str2double, tbl(includeInstance, autoCountInds));
nThres = size(autoCounts, 2);

naCounts = cellfun(@str2double, tbl(includeInstance, naCountInds));
naCounts = (naCounts == 1);  % 1 for true, 999999 for false
maxCount = size(naCounts, 2)
naCountsVec = zeros(size(naCounts, 1), 1);
for count=1:1:maxCount
    naCountsVec(naCounts(:, count)) = count;
end
naCounts = naCountsVec;

% eliminate cases where we don't have a clear designation
autoCounts = autoCounts(naCounts > 0, :);
naCounts = naCounts(naCounts > 0);
nInst = length(naCounts);
isLtMax = (naCounts < maxCount);
% isLtMax = isLtMax & (autoCounts < maxCount)';

% sum(isLtMax, 1)

% if we say more than maxCount arms, reduce that to maxCount
autoCounts(autoCounts > maxCount) = maxCount;

countDiffs = autoCounts - repmat(naCounts, 1, nThres);

size(autoThresholds)
size(countDiffs)

figure
% plot(autoThresholds, sum(abs(countDiffs), 1) / nInst)
plot(autoThresholds, sum(countDiffs, 1) / nInst)
xlabel('threshold')
ylabel('average arm-count difference');

figure
% plot(autoThresholds, sum(abs(countDiffs(isLtMax, :)), 1) ./ sum(isLtMax, 1));
plot(autoThresholds, sum(countDiffs(isLtMax, :), 1) ./ sum(isLtMax, 1));
xlabel('threshold')
ylabel('average arm-count difference (non-max-count only)')

figure
plot(autoThresholds, sum(countDiffs == 0, 1) / nInst);
xlabel('threshold');
ylabel('exact-match proportion');

figure
plot(autoThresholds, sum(countDiffs(isLtMax, :) == 0, 1) ./ sum(isLtMax, 1));
xlabel('threshold');
ylabel('exact-match proportion (non-max-count only)');

end