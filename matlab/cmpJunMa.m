function cmpJunMa(resultCsvPath, cmpCsvPath)

% resultCsvPath = 'D:\matlab-output-local\[paper] 2011-11-10 JunMa composite\combined\_junMa.csv';
% confCsvPath = 'D:\matlab-output-local\[paper] 2011-11-10 JunMa composite\combined\_junMa_conf.csv';
% cmpCsvPath = 'D:\My Dropbox\Research\Galaxy\Workspace\Main\matlab\load_data\JunMaMeas.csv';
% outDir = 'D:\matlab-output-local\temp';

[autoTbl, nSkip] = readTableFromFile(resultCsvPath, ',');
autoTbl(2:end, 1) = strrep(autoTbl(2:end, 1), '-red', '');
jmTbl = readTableFromFile(cmpCsvPath, ',');
[autoTblI, jmTblI] = intersectRowsByField(autoTbl, jmTbl, 1, 1, true);
cmpTbl = joinTables(autoTblI, jmTblI, 1, 1);
% confTbl = readTableFromFile(confCsvPath, ',');
% cmpTbl = joinTables(cmpTbl, confTbl, 1, 1);

paAlenwDcoIdx = strmatch('pa_alenWtd_avg_domChiralityOnly', cmpTbl(1, :), 'exact');
paAlenwDco = str2double(cmpTbl(2:end, paAlenwDcoIdx));

paLongestIdx = strmatch('pa_longest', cmpTbl(1, :), 'exact');
paLongest = str2double(cmpTbl(2:end, paLongestIdx));

paJunMaMeanIdx = strmatch('JunMaPaMean', cmpTbl(1, :), 'exact');
paJunMaMean = str2double(cmpTbl(2:end, paJunMaMeanIdx));

confJunMaIdx = strmatch('JunMaPaHalfWidth', cmpTbl(1, :), 'exact');
confJunMa = 2*str2double(cmpTbl(2:end, confJunMaIdx));

[paAlenwDco, paJunMaMean] = alignChiralityConventions(paAlenwDco, paJunMaMean);
[paLongest, paJunMaMean] = alignChiralityConventions(paLongest, paJunMaMean);

has2pa = confJunMa ~= 0;
paSpread = abs(confJunMa(has2pa));

distEdges = [0:1:40]';
distEdges = [0:0.01:40]';

fprintf('found %d galaxies for comparison\n', length(paJunMaMean));

% half-width gives confidence interval width, so need to double to get a
% comparable difference measure
figure; hist(paSpread, -39:2:39); title('Jun Ma within-galaxy pitch angle deviations');
%saveas(gcf, [outDir 'JunMaWithinGalaxyPitchAngleDeviations.png']);

cnts = histc(paSpread, distEdges);
figure; hold on
line([0 40], [1 1], 'LineStyle', '-', 'Color', 'r');
plot(distEdges, cumsum(cnts) / sum(has2pa));
xlabel('difference threshold'); ylabel('proportion included');
axis([0 40 0 1.1]);

figure; hold on;
plot([0 40], [0 40]);
scatter(abs(paAlenwDco), abs(paJunMaMean), 'b.'); axis([0 40 0 40]); axis square
legend('y = x', 'Location', 'Best');
xlabel('paAlenwDco'); ylabel('paJunMaMean');
%saveas(gcf, [outDir 'JunMaPaVsAlenwDcoPa.png']);

figure; hist(abs(paAlenwDco) - abs(paJunMaMean), -39:2:39); title('paAlenwDco - paJunMaMean');
%saveas(gcf, [outDir 'paDiffHistAlenwDco.png']);

cnts = histc(abs(abs(paAlenwDco) - abs(paJunMaMean)), distEdges);
figure; hold on
line([0 40], [1 1], 'LineStyle', '-', 'Color', 'r');
plot(distEdges, cumsum(cnts) / length(paJunMaMean));
xlabel('difference threshold'); ylabel('proportion included');
axis([0 40 0 1.1]);

figure; hold on
plot([0 40], [0 40]);
scatter(abs(paJunMaMean), abs(paLongest), 'b.'); axis([0 40 0 40]); axis square
xlabel('paJunMaMean'); ylabel('paLongest');

figure; hist(abs(paLongest) - paJunMaMean, -39:2:39); title('paLongest - paJunMaMean');

cnts = histc(abs(abs(paLongest) - paJunMaMean), distEdges);
figure; hold on
line([0 40], [1 1], 'LineStyle', '-', 'Color', 'r');
plot(distEdges, cumsum(cnts) / length(paJunMaMean));
xlabel('difference threshold'); ylabel('proportion included');
axis([0 40 0 1.1]);

figure; hold all
cntsJunMaInternal = histc(paSpread, distEdges);
cntsAlenwDco = histc(abs(abs(paAlenwDco) - abs(paJunMaMean)), distEdges);
cntsLongest = histc(abs(abs(paLongest) - abs(paJunMaMean)), distEdges);
line([0 40], [1 1], 'LineStyle', '-', 'Color', 'k');
h1 = plot(distEdges, cumsum(cntsJunMaInternal) / sum(has2pa));
h2 = plot(distEdges, cumsum(cntsAlenwDco) / length(paJunMaMean));
h3 = plot(distEdges, cumsum(cntsLongest) / length(paLongest));
legend([h1 h2 h3], 'JunMa internal difference', '|JunMaMean - paAlenwDco|', '|JunMaMean - paLongest|', 'Location', 'SouthEast');
axis([0 40 0 1.1]);
xlabel('difference threshold'); ylabel('proportion included');
%saveas(gcf, [outDir 'DistanceDistributions.png']);

cntsAlenwDco2paOnly = histc(abs(abs(paAlenwDco(has2pa)) - abs(paJunMaMean(has2pa))), distEdges);

figure; hold all
cordr = get(gca, 'ColorOrder');
set(gca, 'ColorOrder', cordr([3 1 2 4:end], :))
plot(distEdges, cumsum(cntsJunMaInternal) / sum(has2pa), 'LineWidth', 1.5);
plot(distEdges, cumsum(cntsAlenwDco) / length(paJunMaMean), 'LineWidth', 1.5);
plot(distEdges, cumsum(cntsAlenwDco2paOnly) / sum(has2pa), 'LineWidth', 1.5);
hXlbl = xlabel('difference threshold'); 
hYlbl = ylabel('proportion included');
hLgdLbl = legend('within-galaxy discrepancy', 'between-method discrepancy', 'between-method discrepancy (2 meas. available)', 'Location', 'SouthEast');
axis([0 20 0 1.01])
set(gcf, 'Position', [100 100 650 300]);
% set(gca, 'YTick', 0:0.1:1);
set(gca, 'YTick', 0:0.25:1);
set(gcf, 'Color', [1 1 1]);
set(gca, 'FontSize', 16);
set(hXlbl, 'FontSize', 18);
set(hYlbl, 'FontSize', 18);
set(hLgdLbl, 'FontSize', 12);

junMaPa1idx = strmatch('JunMaPitchAngle1', cmpTbl(1, :), 'exact');
junMaPa1 = str2double(cmpTbl(2:end, junMaPa1idx));
junMaPa2idx = strmatch('JunMaPitchAngle2', cmpTbl(1, :), 'exact');
junMaPa2 = str2double(cmpTbl(2:end, junMaPa2idx));

fprintf('number of galaxies (in Jun Ma catalog) with 2 pitch angle measurements: %d\n', sum(has2pa));

junMaMinPa = min([junMaPa1(:) junMaPa2(:)], [], 2);
junMaMaxPa = max([junMaPa1(:) junMaPa2(:)], [], 2);

inPaRangeAlenwDco = has2pa & (abs(paAlenwDco) > junMaMinPa) & (abs(paAlenwDco) < junMaMaxPa);

fprintf('number with measured pitch angle (alenWtdDco) inside range: %d of %d (%2.4f%%)\n',...
    sum(inPaRangeAlenwDco), sum(has2pa), 100*sum(inPaRangeAlenwDco)/sum(has2pa));

figure; hold all
scatter(1:sum(has2pa), junMaMinPa(has2pa), 'o');
scatter(1:sum(has2pa), junMaMaxPa(has2pa), 'o');
scatter(1:sum(has2pa), abs(paAlenwDco(has2pa)), 'x');
legend('JunMaMin', 'JunMaMax', 'paAlenwDco', 'Location', 'Best');
axis([0 45 0 40]);
title('pitch angle values (alenwDco and JunMa)');
%saveas(gcf, [outDir 'paValsAlenwDco.png']);

inPaRangeLongest = has2pa & (abs(paLongest) > junMaMinPa) & (abs(paLongest) < junMaMaxPa);

fprintf('number with measured pitch angle (longest) inside range: %d of %d (%2.4f%%)\n',...
    sum(inPaRangeLongest), sum(has2pa), 100*sum(inPaRangeLongest)/sum(has2pa));

figure; hold all
scatter(1:sum(has2pa), junMaMinPa(has2pa), 'o');
scatter(1:sum(has2pa), junMaMaxPa(has2pa), 'o');
scatter(1:sum(has2pa), abs(paLongest(has2pa)), 'x');
legend('JunMaMin', 'JunMaMax', 'paAlenwDco', 'Location', 'Best');
axis([0 45 0 40]);

ivlDistAlenwDco = min(abs(abs(paAlenwDco) - abs(junMaPa1)), abs(abs(paAlenwDco) - abs(junMaPa2)));
ivlDistAlenwDco(inPaRangeAlenwDco) = 0;
cnts = histc(ivlDistAlenwDco, distEdges);
figure; hold on
line([0 40], [1 1], 'LineStyle', '-', 'Color', 'r');
plot(distEdges, cumsum(cnts) / length(paAlenwDco));
title('distance from nearest Jun Ma measurement (0 if in interval), paAlenwDco');
xlabel('difference threshold'); ylabel('proportion included');
axis([0 40 0 1.1]);

[sV, sI] = sort(ivlDistAlenwDco, 'descend');
names = cmpTbl(2:end, 1);
fprintf('interval distances, worst to best:\n');
% fprintf('name, dist, junMaSpread\n');
for ii=1:1:length(sI)
    fprintf('%s,\t %2.4f\n', names{sI(ii)}, ivlDistAlenwDco(sI(ii)));
end

ivlDistLongest = min(abs(abs(paLongest) - abs(junMaPa1)), abs(abs(paLongest) - abs(junMaPa2)));
ivlDistLongest(inPaRangeLongest) = 0;
cnts = histc(ivlDistLongest, distEdges);
figure; hold on
line([0 40], [1 1], 'LineStyle', '-', 'Color', 'r');
plot(distEdges, cumsum(cnts) / length(paLongest));
title('distance from nearest Jun Ma measurement (0 if in interval), paLongest');
xlabel('difference threshold'); ylabel('proportion included');
axis([0 40 0 1.1]);

% paSdIdx = strmatch('pitchAngleSd', cmpTbl(1, :), 'exact');
paSdIdx = strmatch('paErr_alenWtd_stdev_domChiralityOnly', cmpTbl(1, :), 'exact');
paSd = str2double(cmpTbl(2:end, paSdIdx));

figure; scatter(paSd(has2pa), paSpread);
title('spread comparison (excluding galaxies with only one Jun Ma measure)');
xlabel('pitch angle sd'); ylabel('JunMa pitch angle spread');
% saveas(gcf, [outDir 'confMeasCmp.png']);

function [chirality1, chirality2] = alignChiralityConventions(chirality1, chirality2)
    if sum(sign(chirality1) == sign(chirality2)) < (length(chirality1)/2)
        chirality1 = -chirality1;
    end
end

end
