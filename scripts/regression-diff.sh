#!/bin/sh
DIR="`dirname $0`"
awk "`cat $DIR/misc.awk`"'
    BEGIN{MAX_REL_ERR=1e-4}
    {for(i=1;i<=NF;i++)L[ARGIND][FNR][i]=$i}
    END{
	diff=0;
	ASSERT(ARGIND==2, "expecting exactly two filenames as input");
	a=ARGIND;
	for(l=1;l<=length(L[a]);l++) for(i=1;i<=length(L[a][l]);i++) {
	    colName=L[1][1][i];
	    if(!index(colName,"time") && L[1][l][i]!=L[a][l][i]) {
		printErrMsg = 0;
		isNum = 1*L[1][l][i] + 1*L[a][l][i];
		if(isNum) {
		    if(L[1][l][i]==int(L[1][l][i]) && L[a][l][i]!=int(L[a][l][i])) printErrMsg=1 # both are integers
		    else if(1*L[a][l][i]) { # correct one is non-zero, so we can compute a relative error
			rel_err = (L[1][l][i] - L[a][l][i])/L[a][l][i];
			if(ABS(rel_err) > MAX_REL_ERR) printErrMsg=1
		    }
		    else printErrMsg=1 # correct one is zero but new one is not
		}
		if(printErrMsg) {diff=1; printf "line %d column %d (%s)\n<%s\n>%s\n",l,i, L[1][1][i], L[1][l][i], L[2][l][i] >"/dev/stderr"}
	    }
	}
	exit(diff);
    }' "$@"
exit $?
