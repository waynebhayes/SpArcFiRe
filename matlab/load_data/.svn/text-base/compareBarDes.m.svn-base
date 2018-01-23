
resultCsvPath = 'D:\matlab-output-local\2011-10-03 2measBarScore scFaceOn\2measBarScore.csv';

autoTbl = readTableFromFile(resultCsvPath, ',');
gzTbl = readTableFromFile('D:\Galaxy-Sets\Bamford GZ2\gz2master.csv', ',');
[autoTblI, gzTblI] = intersectRowsByField(autoTbl, gzTbl, 1, 1, true);
clear autoTbl
clear gzTbl
cmpTbl = joinTables(autoTblI, gzTblI, 1, 1);

gzBarVProp = str2double(cmpTbl(2:end, 145));
barScores = str2double(cmpTbl(2:end, 34));
barCandScores = str2double(cmpTbl(2:end, 35));
desBar = (barScores >= 2) & (barCandScores > 7);
cbr = str2double(cmpTbl(2:end, 12));
cbrOK = (cbr > 1);
gzIsBar = gzBarVProp > 0.5;
figure; hold all
scatter(barScores(gzIsBar & cbrOK), barCandScores(gzIsBar & cbrOK), 'r.');
scatter(barScores(~gzIsBar & cbrOK), barCandScores(~gzIsBar & cbrOK), 'b.');