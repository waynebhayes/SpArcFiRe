% set arcsCsvPath before running, and optionally plotPath and/or outCsvPath

% arcsCsvPath = 'D:\matlab-output-local\2011-11-10 JunMa composite\combined\_junMa_arcs.csv';
% outCsvPath = 'D:\matlab-output-local\2011-11-10 JunMa composite\combined\_junMa_conf2.csv';

if ~exist('plotPath', 'var')
    plotPath = [];
else
    plotPath = [plotPath filesep];
end

arcsCsv = readTableFromFile(arcsCsvPath, ',');

nameIdx = strmatch('gxyName', arcsCsv(1, :), 'exact');
paIdx = strmatch('pitch_angle', arcsCsv(1, :), 'exact');
alenIdx = strmatch('arc_length', arcsCsv(1, :), 'exact');
pixCntIdx = strmatch('num_pixels', arcsCsv(1, :), 'exact');

names = arcsCsv(2:end, nameIdx);
pas = str2double(arcsCsv(2:end, paIdx));
alens = str2double(arcsCsv(2:end, alenIdx));
pixCnts = str2double(arcsCsv(2:end, pixCntIdx));

errs = atan(pixCnts./(2*(alens.^2))) * (180/pi);

nameList = unique(names);

paMeans = zeros(1, length(nameList));
paSds = zeros(1, length(nameList));
arcCounts = zeros(1, length(nameList));
arcErrMeans = zeros(1, length(nameList));
arcErrMeansAlenw = zeros(1, length(nameList));
paMeansDco = zeros(1, length(nameList));
paSdsDco = zeros(1, length(nameList));
arcCountsDco = zeros(1, length(nameList));
arcErrMeansDco = zeros(1, length(nameList));
arcErrMeansAlenwDco = zeros(1, length(nameList));
awps = zeros(1, length(nameList));
x = -180:0.001:180;
for ii=1:1:length(nameList)
    curIdxs = strmatch(nameList{ii}, names, 'exact');
    curPas = pas(curIdxs);
    curAlens = alens(curIdxs);
    curErrs = errs(curIdxs);
    assert(all(curErrs >= 0))
    if ~isempty(plotPath)
        figure; subplot(2, 1, 1); hold all;
    end
    mv = zeros(size(x));
    for jj=1:1:length(curIdxs)
        cv = curAlens(jj) * normpdf(x, curPas(jj), curErrs(jj));
        if ~isempty(plotPath)
            plot(x, cv);
        end
        mv = mv + cv;
    end
    mu = sum(x .* mv) / sum(mv);
    sigma = sqrt(sum(mv .* (x - mu).^2) / sum(mv));
    paMeans(ii) = mu;
    paSds(ii) = sigma;
    arcCounts(ii) = length(curIdxs);
    arcErrMeans(ii) = mean(curErrs);
    arcErrMeansAlenw(ii) = sum(curErrs .* curAlens) / sum(curAlens);
    if ~isempty(plotPath)
%         title(regexprep(nameList{ii}, '_', '-'));
        title([regexprep(nameList{ii}, '_', '-') sprintf(': \\mu = %2.4f, \\sigma = %2.4f', mu, sigma)]);
        subplot(2, 1, 2);
        plot(x, mv);
        saveas(gca, [plotPath nameList{ii} '_errMeas.png']);
        close all
    end
    
%     nameList{ii}
%     curPas
%     curAlens
%     curPas .* curAlens
%     sum(curPas .* curAlens)

    curAwps = sum(curPas .* curAlens);
    curDomChir = sign(curAwps);
    isDco = (sign(curPas) == curDomChir);
    dcoPas = curPas(isDco);
    dcoAlens = curAlens(isDco);
    dcoErrs = curErrs(isDco);
    assert(all(dcoErrs >= 0));
    
    mv = zeros(size(x));
    for jj=1:1:length(dcoPas)
        mv = mv + ( dcoAlens(jj) * normpdf(x, dcoPas(jj), dcoErrs(jj)) );
    end
    mu = sum(x .* mv) / sum(mv);
    sigma = sqrt(sum(mv .* (x - mu).^2) / sum(mv));
    paMeansDco(ii) = mu;
    paSdsDco(ii) = sigma;
    arcCountsDco(ii) = sum(isDco);
    awps(ii) = curAwps;
    arcErrMeansDco(ii) = mean(dcoErrs);
    arcErrMeansAlenwDco(ii) = sum(dcoErrs .* dcoAlens) / sum(dcoAlens);
end

if exist('outCsvPath', 'var')
    outCsvFile = fopen(outCsvPath, 'wt');
    fprintf(outCsvFile, ['gxyName,pitchAngleMean,pitchAngleSd,arcErrMean,arcErrMeanAlenw,numArcs,'...
        'pitchAngleMeanDco,pitchAngleSdDco,arcErrMeanDco,arcErrMeanAlenwDco,numArcsDco,alenWtdPangSum\n']);
    for ii=1:1:length(nameList)
        fprintf(outCsvFile, '%s,%f,%f,%f,%f,%d,%f,%f,%f,%f,%d,%f\n', ...
            nameList{ii}, paMeans(ii), paSds(ii), arcErrMeans(ii), arcErrMeansAlenw(ii), arcCounts(ii), ...
            paMeansDco(ii), paSdsDco(ii), arcErrMeansDco(ii), arcErrMeansAlenwDco(ii), arcCountsDco(ii), awps(ii));
    end
end

fclose(outCsvFile);
        
        