#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [-h] [-v] [data file]
PURPOSE: compute holistic p-value of possibly correlated variates using Empirical Brown's Method (Poole 2016)
OPTIONS:
-h: print help
-v: verbose output
INPUT:
m+1 lines, with n+2 columns each.
The top line is a 'header' line; this line is completely ignored, and can be whatever you want
Each line represents a variable; the first entry is the 'name' (whatever you want, but must be non-empty string);
then the p-value; and then the n raw samples of that variable.
The Empirical Brown's Method (Poole 2016) computes the covariance of all the variables and spits out a holistic
p-value that is >= product of p-values (ie., less significant), accounting approximately for inter-dependencies.
NOTE: There is a compiled C version of this code, which runs about 50-100x faster, in my libwayne repo in tests/ebm.c.
This script will check to see if a compiled version named 'ebm' is in the PATH, verifies it's correctness if it
exists, and quietly runs it if it checks out."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" = skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME = "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# = 0 ] && tty --silent && warn "reading stdin from tty; press ^D to finish, or ^C to exit"

VERBOSE=0
ALLOW_EXE=true
while true; do
    case "$1" in
    -h) die "<printing help message only>" ;;
    -[Vv]*) VERBOSE=1; shift;;
    -noexe) ALLOW_EXE=false; shift;;
    -*) die "option '$1' not supported";;
    *) break;;
    esac
done

export VERBOSE
DIRNAME=`dirname "$0"`
EBM_EXE=`(/bin/which ebm || /usr/bin/which ebm) 2>/dev/null | head -1`
if $ALLOW_EXE && [ -x "$EBM_EXE" ] && VERBOSE=1 "$EBM_EXE" <"$DIRNAME/ebm.test.in" 2>/dev/null | cmp - "$DIRNAME/ebm.test.out" >/dev/null 2>&1; then
    if [ "$EBM_SH_TRYING_EXECUTABLE" = "" ]; then
	[ "$VERBOSE" -gt 1 ] && echo "exec'ing $EBM_EXE" >&2
	EBM_SH_TRYING_EXECUTABLE=true; export EBM_SH_TRYING_EXECUTABLE
	exec "$EBM_EXE" "$@"
    fi
fi

tail -n +2 "$@" > $TMPDIR/input

hawk 'function TransformData(n,     j,cdf) {
	for(j=1;j<=n;j++) {
	    cdf=StatHistECDF(NR, data[NR][j]);
	    ASSERT((NR in _statN) && _statN[NR], "TransformData: \""NR"\" has no samples");
	    if(cdf==0) cdf = 1.0/(_statN[NR]*_statN[NR]); # really small but not zero
	    data[NR][j] = -2*log(cdf);
	}
    }
	NR==1{n=NF-2}
	{ # every line including NR==1, because the header line was discarded using tail(1) above.
	    ASSERT(NF-2==n,"number of columns must be constant; first column had "n" but this one has "NF-2);
	    name[NR]=$1; pVal[NR]=$2;
	    for(i=1;i<=n;i++) {
		data[NR][i]=sample=$(i+2)
		StatAddSample(NR,sample);
	    }
	    #printf "DATA[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " %g",data[NR][i]>"/dev/stderr";print "">"/dev/stderr";
	    for(i=1;i<=n;i++){
		data[NR][i] = (data[NR][i] - StatMean(NR))/StatStdDev(NR);
		StatHistAddSample(NR, data[NR][i]);
	    }
	    #printf "NORM[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " %g",data[NR][i]>"/dev/stderr";print "">"/dev/stderr";
	    StatHistMakeCDF(NR);
	    #printf "ECDF[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " [%d][%d][%g]=%g",NR,i,data[NR][i],StatHistECDF(NR,data[NR][i])>"/dev/stderr";print "">"/dev/stderr";
	    TransformData(n); # note it always works on the current line (NR)
	    #printf "-LOG[%d]",NR>"/dev/stderr";for(i=1;i<=n;i++)printf " %g",data[NR][i]>"/dev/stderr";print "">"/dev/stderr";
	}
END{
    #printf "last line has pval[%d]=%g\n", NR,pVal[NR] > "/dev/stderr";
    m=NR
    ASSERT(m>0, "need at least one variable");
    # Compute covariances across all samples of all input variables.
    for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) for(k=1;k<=n;k++) CovarAddSample(i" "j, data[i][k], data[j][k]);
    # Now perform Empirical Browns Method
    df_fisher = Expected = 2.0*m;
    cov_sum = 0;
    for(i=1;i<=m;i++) { 
	#printf "row %4d", i > "/dev/stderr"
	for(j=i+1;j<=m;j++) {
	    covar=CovarCompute(i" "j);
	    #printf "\t[%d,%d] %.5f", i,j, covar > "/dev/stderr";
	    cov_sum += covar;
	}
	#print "" > "/dev/stderr"
    }
    #printf "cov_sum.sh %s %g\n", FILENAME, cov_sum > "/dev/stderr";
    Var = 4.0*m+2*cov_sum;
    c = Var/(2.0*Expected);
    df_brown = 2.0*Expected**2/Var;
    if(df_brown > df_fisher) {
	df_brown = df_fisher
	c = 1.0
    }
    #printf "c = %g\n", c > "/dev/stderr"
    pProd=1; log_pProd=0;
    x=0; # twice the sum of logs of p-values
    for(i=1;i<=m;i++)
	if(1*pVal[i]>0) {
	    x += -log(pVal[i]); log_pProd+=log(pVal[i]); pProd *= pVal[i]; #printf "x[%d]=%g\n",i,x
	}
	else if('$VERBOSE') Warn(sprintf("skipping pVal[%d]=%g",i,pVal[i]));
    x *= 2;
    #printf("x = %g\n", x) > "/dev/stderr";
    log_p_brown = logChi2_pair(int(df_brown+0.5), 1.0*x/c)
    p_brown = Exp(log_p_brown)

    ASSERT(p_brown >= pProd,
	"Oops, something wrong: p_brown should be < product(pVals), but p_brown ="p_brown" product ="pProd);
    if('$VERBOSE') fmt="p-value < %g = bitscore %g (product gives %g ; bitscore %g )\n"
    else fmt="%g %g %g %g\n";
    printf(fmt, p_brown, -log_p_brown/log(2), pProd, -log_pProd/log(2));
}' $TMPDIR/input
