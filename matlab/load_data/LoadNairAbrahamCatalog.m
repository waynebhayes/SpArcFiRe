cattbl = readTableFromFile('D:\Galaxy-Sets\NairAbrahamMorphology.cat');
idtbl = readTableFromFile('D:\Galaxy-Sets\SDSSids.tbl');

idtbl(2:end, 3) = regexprep(idtbl(2:end, 3), 'sp-|\.id', '');

NairAbrahamCatalog = joinTables(cattbl, idtbl, 13, 3);

clear cattbl idtbl