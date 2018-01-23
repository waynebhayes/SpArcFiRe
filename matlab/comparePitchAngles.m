function comparePitchAngles(tbl, paIdxs, refPas, tblHasHeader, measNames)

error(nargchk(5, 5, nargin))

% check to make sure than refPas is a vector

if tblHasHeader
    tbl = tbl(2:end, :);
end

% measPas = str2num(tbl(:, paIdxs));
measPas = cellfun(@str2double, tbl(:, paIdxs));
nImgs = size(measPas, 1);
nMeas = size(measPas, 2);

measPas = abs(measPas); % reference pitch angles are all unsigned

errs = measPas - repmat(refPas, 1, nMeas);

markerTypes = {'+', 'o', '*', 'x', 's', 'd', 'p'};
% colors = {'r', 'g', 'b', 'c', 'm', 'k'};
colors = {'b', 'g', 'r', 'c', 'm', 'k'};

% perImgErr = figure;
figure; hold on
for meas=1:1:nMeas
    scatter(1:1:nImgs, abs(errs(:, meas))', markerTypes{meas});
end
legend(measNames{:}, 'Location', 'Best')
xlabel('image');
ylabel('absolute-error');

figure; hold on
for meas = 1:1:nMeas
    plot(1:1:nImgs, sort(abs(errs(:, meas)))', [colors{meas} '-']);
end
legend(measNames{:}, 'Location', 'NorthWest');
xlabel('error rank');
ylabel('error value');

figure; hold on
for meas = 1:1:nMeas
    plot(1:1:nImgs, cumsum(sort(abs(errs(:, meas))))', [colors{meas} '-']);
end
legend(measNames{:}, 'Location', 'NorthWest');
xlabel('error rank');
ylabel('cumulative error');

% figure; hold on
% for meas=1:1:nMeas
%     scatter(measPas(:, meas), refPas, markerTypes{meas});
% end

for meas = 1:1:nMeas
    curErrs = errs(:, meas);
    figure; hist(abs(curErrs))
    title(sprintf('absolute errors for %s (mean-abs: %.4f, mean: %.4f)',...
        measNames{meas}, mean(abs(curErrs)), mean(curErrs)));
end

end