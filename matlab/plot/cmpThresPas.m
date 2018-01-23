function cmpThresPas(gxyParamsList, paFieldName, truePas, thresholds)

nInst = length(gxyParamsList);
nThres = length(thresholds);

pas = cellfun(@(x)(getfield(x, paFieldName)), gxyParamsList, 'UniformOutput', false);
pas = reshape(cell2mat(pas), nThres, nInst)';


% validateattributes(truePas, {}, {'column'});
% 
% nInst = size(pas, 1);
% nThres = length(thresholds);
% 
% if nInst ~= length(truePas)
%     error(['number of pitch-angle vectors must match number of true '...
%         'pitch angles'])
% end
% 
% if size(pas, 2) ~= nThres
%     error(['number of thresholded pitch angle measurements must equal '...
%         'the number of thresholds'])
% end

errs = abs(abs(pas) - repmat(truePas, 1, nThres));
inclErrs = (pas ~= 0); % TODO: change this to some other indicator
% meanThresErrs = mean(errs, 1);
% medianThresErrs = median(errs, 1);
meanThresErrs = zeros(nThres, 1);
medianThresErrs = zeros(nThres, 1);
lqThresErrs = zeros(nThres, 1);
uqThresErrs = zeros(nThres, 1);

for meas=1:1:nThres
    curErrs = errs(inclErrs(:, meas), meas);
    meanThresErrs(meas) = mean(curErrs);
    medianThresErrs(meas) = median(curErrs);
    lqThresErrs(meas) = quantile(curErrs, 0.25);
    uqThresErrs(meas) = quantile(curErrs, 0.75);
end

figure; hold on
plot(thresholds, sum(inclErrs, 1) / nInst);
xlabel('threshold');
ylabel('measurement availability rate');

figure; hold on
plot(thresholds, meanThresErrs, 'r-');
plot(thresholds, medianThresErrs, 'b-');
plot(thresholds, lqThresErrs, 'b:');
plot(thresholds, uqThresErrs, 'b:');
title(regexprep(paFieldName, '_', '-'))
xlabel('threshold');
ylabel('mean absolute error in pitch angle');
legend('mean', 'median', '1st quartile', '3rd quartile');

figure; hold on
plot(thresholds, meanThresErrs);
plot(thresholds, medianThresErrs);
title(regexprep(paFieldName, '_', '-'))
xlabel('threshold');
ylabel('mean absolute error in pitch angle');
for ii=1:1:size(errs, 1)
    scatter(thresholds(inclErrs(ii, :)), errs(ii, inclErrs(ii, :)));
end

end