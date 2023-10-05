#!/bin/sh
TMP=/tmp/Cora-gawk.$$
trap "/bin/rm -f $TMP" 0 1 2 3 15

gawk 'NR==1{print};NR>1{print toupper($0)}' "$@" | sed 's///' > $TMP # ensure newlines and no carriage returns

gawk -F"	" '
    function ASSERT(cond, s) { if(cond) return; print s >"/dev/stderr"; exit(1)}
    function CommonChars(s,t,    i,j,sum) { # return the number of common characters in two strings
	sum=0
	for(i=1;i<=length(s);i++) for(j=1;j<=length(t);j++) sum += (substr(s,i,1)==substr(t,j,1))
	return sum
    }

    BEGIN{printf "good?    name\tdark\tBEST\tpScores\n"}
    NR==1 { # get the names of the columns
	#printf "reading Header of %d columns\n", NF
	for(col=1;col<=NF;col++) {
	    nf=split($col,a," ")
	    #printf "column %d is \"%s\" nf=%d with elements:",col,$col,nf
	    #for(i=1;i<=nf;i++)printf " %s",a[i]; print ""
	    if(a[2]=="label") {
		ASSERT(nf==2,"label name without 2-element array");
		++bandpairs[a[1]];
		labelCol[a[1]]=col;
	    } else if(index($col,"_p")) {
		ASSERT(nf==2,"_p name with nf="nf);
		pValCol[a[1]]=col
	    } else if(nf==1) colNum[a[1]]=col
	}
	#print "parsed"
	#for(pair in labelCol) printf "pair %s has labels in column %d\n",pair,labelCol[pair]
	#for(pair in pValCol)  printf "pair %s has pVals in column %d\n",pair,pValCol[pair]
    }

    NR>1{
	delete logPsum;
	for(pair in bandpairs)
	    logPsum[$labelCol[pair]] += log($pValCol[pair]);
	which=minP=0;
	scoreString="";
	for(label in logPsum){
	    scoreString = scoreString sprintf("\t%s %g",label,logPsum[label]);
	    if(logPsum[label]<minP){
		minP=logPsum[label];which=label
	    }
	}
	printf "%d %13s\t%s\t%s%s\n", !!CommonChars(which,$colNum["dark"]), $colNum["name"],
	    $colNum["dark"],which, scoreString
    }' $TMP
