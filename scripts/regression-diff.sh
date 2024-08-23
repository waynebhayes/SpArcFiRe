#!/bin/sh
die() { echo "$USAGE${NL}FATAL ERROR: $@" >&2; exit 1;}
TestGawk(){ echo hello | gawk '{++array[$1][$1]}' || die "Upgrade to a [g]awk that supports multi-dimensional arrays"; }
TestGawk

DIR="`dirname $0`"
awk "`cat $DIR/misc.awk`"'
    BEGIN{MAX_REL_ERR=0.02}
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
		else if(ABS(L[1][l][col1])>1e-12) # correct one is zero and new one is sig larger than machine EPS.
		    printErrMsg=1
	    }
	    if(printErrMsg) {
		++diff;
		if(isNum) Warn(sprintf("line %d columns %d,%d (%s): numerical values beyond tolerance: \"%s\" <-> \"%s\"",
			l,col1,col2, varName, L[1][l][col1], L[2][l][col2]))
		else Warn(sprintf("line %d columns %d,%d (%s): non-numerical values differ: \"%s\" <-> \"%s\"",
		    l,col1,col2, varName, L[1][l][col1], L[2][l][col2]))
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
	    for(v in varCols) if(length(varCols[v])!=numTSVs){
		Warn(sprintf("header column name \"%s\" does not appear in all input files",v));
		yes=no="";
		for(i=1;i<=ARGC;i++)if(i in varCols[v])yes=yes" "F[i];else no=no" "F[i];
		Warn(sprintf("It exists in (%s) but is missing in (%s)",yes,no));
		++headerMismatch[v]; # record column number (but not variable name)
		delete varCols[v];
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

	headErr = (isarray(headerMismatch)?length(headerMismatch):0);
	if(headErr || diff) {
	    Warn(sprintf("\n****************\n****************\n**************** ABOVE ERRORS (%d header, %d diff) MAY BE FATAL; SUPPRESSING FURTHER WARNINGS\n****************\n****************\n", headErr, diff));
	    exit(headErr + diff);
	}

	for(l=2;l<=length(L[1]);l++)
	    if(numTSVs)
		for(v in varCols) CheckCol(l, v, varCols[v][1], varCols[v][2]);
	    else
		for(i=1;i<=length(L[1][l]);i++) CheckCol(l, "", i, i)
	exit(diff);
    }' "$@"
exit $?
