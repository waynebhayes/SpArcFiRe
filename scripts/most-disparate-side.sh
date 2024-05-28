#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [-ebm] normed-ellipse-flux.tsv
PURPOSE: given a file (from Cora Schallock) with fluxes-per-pixel (normalized separately per-waveband inside the
    chosen ellipse of the bulge), compute the K-S test p-value, for each waveband pair, that they are different distributions.
    If given the -ebm argument, use Empirical Brown's Method to compute the combined p-value across all waveband pairs for
    each side."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
BIGTMP=`for i in /var/tmp/SF$$ /tmp/SF$$; do mkdir -p "$i" && (df $i | awk 'NR==1{for(av=1;av<=NF;av++)if(match($av,"[Aa]vail"))break;}NR>1{print $av,"'"$i"'"}'); done 2>/dev/null | sort -nr | awk 'NR==1{print $2}'`
[ "$MYTMP" ] || MYTMP="$BIGTMP"
TMPDIR=`mktemp -d $MYTMP/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR $BIGTMP; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# -gt 0 ] || die "expecting at least one argument"

echo TMPDIR is $TMPDIR >&2

EBM=0
EBM_NOEXE=""

case "$1" in
-h) die "help displayed";;
-ebm-noexe) EBM=1; EBM_NOEXE="-noexe"; shift;;
-ebm) EBM=1; shift;;
esac

for i
do
    echo --- "$i" --- >&2
    awkcel 'BEGIN{EBM='$EBM'; outFile="'$TMPDIR'/EBM.in"}
	{
	    for(i=1;i<NF-1;i++) for(j=i+1;j<NF;j++) {
		name=i j side;
		StatHistAddSample(name,$i-$j);
		StatAddSample(name,$i-$j)
		if(EBM) px[i j][NR]=($i-$j);
	    }
	}
	END{

	    sumLog[0]=sumLog[1]=0;
	    if(EBM) for(vote=0;vote<=1;vote++) {
		printf "XXX\tks_p_val" > outFile vote
		for(l=1;l<=NR;l++) printf "\tpx%d", l > outFile vote
		print "" > outFile vote
	    }
	    for(i=1;i<NF-1;i++) for(j=i+1;j<NF;j++) {
		n0=i j 0; n1=i j 1; # names
		M0=StatMean(n0); SD0=StatStdDev(n0); V0=StatVar(n0);
		M1=StatMean(n1); SD1=StatStdDev(n1); V1=StatVar(n1);
		vote=M0 < M1;
		gap=ABS(M0-M1);
		ks=KStest(n0, n1)/1.3;

		# See if the confidence intervals overlap, tuned to failure rate fR
		fR = 1e-2; # failure rate, eg 1e-2=99% confidence, 1e-3 = 99.9% confidence
		ci0 = StatConfidenceInterval(n0,1-fR);
		ci1 = StatConfidenceInterval(n1,1-fR);
		good = (gap > MAX(ci0,ci1));
		if(good) sumLog[vote]+=log10(ks);
		status = good ? "success" : "fail";
		if(EBM) {
		    printf "%s_%s\t%g", colName[i], colName[j], ks > outFile vote
		    for(l=1;l<=NR;l++) printf "\t%g", px[i j][l] > outFile vote
		    print "" > outFile vote
		}
		printf "%s-%s\t%d (gap %g, CI %g, %7s)\t%g\n", colName[i], colName[j], vote, gap, MAX(ci0, ci1), status, ks
	    }
	    winner=sumLog[1]<sumLog[0];
	    printf "WINNER side%d with log10(p)=%g, loser log10(p)=%g LOG10 DIFF %g\n", winner, sumLog[winner], sumLog[1-winner],
		sumLog[winner]-sumLog[1-winner]
	}' "$i"
    [ "$EBM" -ne 1 ] && continue
    #set -x
    for vote in 0 1; do
	echo -n "side $vote: "
	if [ `wc -l < $TMPDIR/EBM.in$vote` -gt 1 ]; then
	    ebm.sh -v $EBM_NOEXE < $TMPDIR/EBM.in$vote
	else
	    echo "no votes"
	fi
    done
done
