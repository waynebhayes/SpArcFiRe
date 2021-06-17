#!/bin/sh
DIR="`dirname $0`"
awk "`cat $DIR/misc.awk`"'
    BEGIN{MAX_REL_ERR=0.01}
    FNR==1{
	F[ARGIND]=FILENAME;
	if(index(FILENAME,".tsv")){
	    ++numTSVs;
	    for(i=1;i<=NF;i++)varCols[$i][ARGIND]=i; # record column number of each variable name
	}
    }
    {for(i=1;i<=NF;i++)L[ARGIND][FNR][i]=$i}
    function CheckCol(l, varName, col1, col2) {
	if(!index(varName,"time") && L[1][l][col1]!=L[2][l][col2]) {
	    printErrMsg = 0;
	    isNum = 1*L[1][l][col1] + 1*L[2][l][col2];
	    if(isNum) {
		if(L[1][l][col1]==int(L[1][l][col1]) && L[2][l][col2]!=int(L[2][l][col2])) printErrMsg=1 # both are integers
		else if(1*L[2][l][col2]) { # correct one is non-zero, so we can compute a relative error
		    rel_err = (L[1][l][col1] - L[2][l][col1])/L[2][l][col2];
		    if(ABS(rel_err) > MAX_REL_ERR) printErrMsg=1
		}
		else printErrMsg=1 # correct one is zero but new one is not
	    }
	    if(printErrMsg) {
		++diff;
		Warn(sprintf("line %d columns %d,%d (%s): \"%s\" <-> \"%s\"",l,col1,col2, varName, L[1][l][col1], L[2][l][col2]))
	    }
	}
    }
    END{
	diff=0;
	ASSERT(ARGIND==2, "expecting exactly two filenames as input");
	# Check that both files have the same number of columns, especially on the header line
	if(length(L[1]) != length(L[2])){++diff; Warn(sprintf("Line count mismatch:\n\t%d lines in \"%s\"\n\t%d lines in \"%s\"", length(L[1]), F[1], length(L[2]), F[2]))}

	if(isarray(varCols)) {
	    ASSERT(numTSVs==ARGIND, "all files must be TSVs, or none");
	    delete headerMismatch;
	    PROCINFO["sorted_in"]="@ind_num_asc"; #traverse for loop based on integer VALUE (not INDEX) of elements
	    for(i in varCols) if(length(varCols[i])!=numTSVs){
		Warn(sprintf("header column name \"%s\" does not appear in all input files",i));
		headerMismatch[varCols[i]]=i; # record column number and variable name
		delete varCols[i];
	    }
	    if(!isarray(headerMismatch)) { # only makes sense to compare header column-by-column if they have same # of columns
		if(length(L[1][1]) != length(L[2][1]))
		Warn(sprintf("Column count mismatch:\n\t%d columns in \"%s\"\n\t%d columns in \"%s\"",
		    length(L[1][1]), F[1], length(L[2][1]), F[2]));
		for(i=0;i<=1;i++) for(j in L[i+1][1]){
		    other=1-i;
		    if(L[i+1][1][j] != L[other+1][1][j]) {
			++diff;
			Warn(sprintf("column %d of \"%s\" is named \"%s\", but the other file has name \"%s\" there",
			    j,F[i+1],L[i+1][1][j], L[other+1][1][j]))
		    }
		}
	    }
	}
	else
	    ASSERT(numTSVs==0, "all files must be TSVs, or none");

	if(isarray(headerMismatch) || diff) {
	    print "\n****************\n****************\n**************** ABOVE ERRORS MAY BE FATAL; SUPPRESSING FURTHER WARNINGS\n****************\n****************\n" > "/dev/fd/1";
	    exit(length(headerMismatch)+diff);
	}

	for(l=2;l<=length(L[1]);l++)
	    if(numTSVs)
		for(v in varCols) CheckCol(l, v, varCols[v][1], varCols[v][2]);
	    else
		for(i=1;i<=length(L[1][l]);i++) CheckCol(l, "", i, i)
	exit(diff);
    }' "$@"
exit $?
