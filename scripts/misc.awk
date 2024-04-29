BEGIN{PI=M_PI=3.14159265358979324;BIGNUM=1*1e30; for(i=0;i<256;i++)ASCII[sprintf("%c",i)]=i}

function ASSERT(cond,str){if(!cond){s=sprintf("ASSERTION failure, line %d of input file %s: %s\nInput line was:\n<%s>\n", FNR,FILENAME,str,$0); print s >"/dev/stderr"; exit 1}}
function WARN(cond,str,verbose){if(!cond){s=sprintf("WARNING: line %d of input file %s: %s",FNR,FILENAME,str); if(verbose)s=s sprintf("\nInput line was:\n<%s>\n", $0); print s >"/dev/stderr"}}
function ABS(x){return x<0?-x:x}
function SIGN(x){return x==0?0:x/ABS(x)}
function MAX(x,y){return x>y?x:y}
function MIN(x,y){return x<y?x:y}

function nul(s){} # do nothing; useful for holding temp strings in code on the command line.

# The default srand() uses time-of-day, which only changes once per second. Not good enough for paraell runs.
function GetFancySeed(           seed,hostname,n,h,i){
    seed=systime()+PROCINFO["gid"]+PROCINFO["uid"]+PROCINFO["pgrpid"]+PROCINFO["ppid"]+PROCINFO["pid"];
    "hostname"|getline hostname; n=length(hostname); for(i=1;i<=n;i++) seed+=ASCII[substr(hostname,i,1)];
    return seed;
}
function strip(s) { gsub("  *"," ",s); sub("  *$","",s); return s}
function randsort(i1,v1,i2,v2) { return rand()-0.5 } # use this to have for loops go in random order
function Srand() { return srand(GetFancySeed());}
function RandInt(mean,radius) {return mean+2*radius*(rand()-0.5)}

function ftos(f){f=sprintf("%.3g",f); gsub("e[+]0","e+",f); return f} # remove leading 0s from exponent
function floor(x) {if(x>=0) return int(x); else return int(x)-1}
function ceil(x) {if(x==int(x)) return x; else return floor(x)+1}
function int2binary(i,l, _s){if(i<0)return "nan";_s="";while(i){_s=(i%2)""_s;i=int(i/2)};while(length(_s)<l)_s="0"_s;return _s}
function Fatal(msg){printf "FATAL ERROR: %s\n",msg >"/dev/stderr"; exit(1);}
function Warn(msg){printf "Warning: %s\n",msg >"/dev/stderr"}
function NormDotProd(u,v,    _dot,_dot1,_dot2,i){_dot=_dot1=_dot2=0;
    for(i in u){_dot+=u[i]*v[i];_dot1+=u[i]*u[i];_dot2+=v[i]*v[i]};
    return _dot/sqrt(_dot1*_dot2);
}
function inarray(element,array,      i){for(i=1;i<=length(array);i++)if(element==array[i])return 1; return 0}
function IsPrime(N,   i){if(N<2)return 0; for(i=2;i<=sqrt(N); i++)if(N/i==int(N/i))return 0; return 1}
function NSORT(a,ai,   i,NsortTc){delete sortTb;delete sortTc; for(i in a)sortTb[a[i]*(1+1e-7*rand())]=i;NsortTc=asorti(sortTb,sortTc);for(i=1;i<=NsortTc;i++)ai[i]=sortTb[sortTc[i]];return NsortTc}
#Bubble Sort: assumes 1-indexed arrays!
function bsort(array,outindices,    N,i,j,temp){
    delete outindices;delete iSortValue;
    N=0; for(i in array){N++;outindices[N]=i;iSortValue[N]=array[i]}
    for(i=2;i<=N;i++)
    {
	# Invariant: array from 1 .. i-1 is sorted, now bubble the next element into its place
	j=i;
	while(j>1 && iSortValue[j] < iSortValue[j-1]) {
	    temp=iSortValue[j-1];iSortValue[j-1]=iSortValue[j];iSortValue[j]=temp;
	    temp=outindices[j-1];outindices[j-1]=outindices[j];outindices[j]=temp;
	    j--;
	}
    }
    return N;
}
function Factor(n,  i,s){s="";i=int(sqrt(n)); while(i>1){while(n/i==int(n/i)){s=s" "i;n/=i} i--}; return s}
function asin(x) { return atan2(x, sqrt(1-x*x)) }
function acos(x) { return atan2(sqrt(1-x*x), x) }
function atan(x) { return atan2(x,1) }
function sind(x) { return sin(x/180*PI) }
function cosd(x) { return cos(x/180*PI) }
function tand(x) { return tan(x/180*PI) }
function asind(x) { return asin(x)/PI*180 }
function acosd(x) { return acos(x)/PI*180 }
function atand(x) { return atan(x)/PI*180 }

#Given x, compute ln(1+x), using the Taylor series if necessary
function AccurateLog1(x,    absX,n,term,sum){
    #return log(1+x); # fuck it
    absX=ABS(x);
    if(absX<1e-16) return x; # close to machine eps? it is just x
    if(absX>4e-6) return log(1+x); # built-in one is very good in this range
    ASSERT(x>=-0.5&&x<=1,"AccurateLog1("x") will not converge");
    if(x in _memAccLog1) return _memAccLog1[x];
    sum=0;
    n=1; term=x;
    delete _log1Terms;
    # This first loop is just to get the terms, not actually computing the true sum
    while(n==1 || ABS(term/sum)>1e-20){sum+=ABS(term); _log1Terms[n++]=term; term*=x/n}
    # Now sum the terms smallest-to-largest, keeping the two signs separate
    for(i=n;i>0;i--)if(_log1Terms[i]<0)_log1Terms[-1]+=_log1Terms[i]; else _log1Terms[0]+=_log1Terms[i]
    sum = _log1Terms[0] + _log1Terms[-1];
    sum -= sum*sum; # I am not sure why, but this gives a MUCH better approximation???
    if(n>_nMaxAccLog1){_nMaxAccLog1=n;
	#printf "AccurateLog: memoize log(1+%.16g)=%.16g, nMax %d\n",x,sum,_nMaxAccLog1 >"/dev/fd/2"
    }
    return (_memAccLog1[x]=sum);
}
# Assuming S=a+b, but we only have log(a) and log(b), we want to compute log(S)=log(a+b)=log(a(1+b/a))=log(a)+log(1+b/a)
# And then we can call AccurateLog1(b/a), except we do not have b/a explicitly, but we can get it with Exp(log_b-log_a)
function LogSumLogs(log_a,log_b,    truth, approx) {
    m=MIN(log_a,log_b)
    M=MAX(log_a,log_b)
    ASSERT(M>=m,"BUG: M is not greater than m in LogSumLogs");
    #return M+log(1+Exp(m-M))
    if(M-m > 37) return M; # m < M*machine_eps, so m will not change M.
    # fuck it
    approx = M+AccurateLog1(Exp(m-M))
    if(ABS(log_a)<700 && ABS(log_b) < 700){
	truth=log(exp(log_a)+exp(log_b));
	if(ABS((approx-truth)/truth)>1e-10)
	    printf "LogSumLogs badApprox: log_a %g log_b %g M %g m %g approx %g truth %g\n",log_a,log_b,M,m,approx,truth > "/dev/stderr"
    }
    return approx
}

# Given the logarithm log(p) of an arbitrarily large or small number p, return a string
# looking like printf("%.Ng", p) printing p with N significant digits (just like printf).
function logPrint(logp,digits,   l10_p,intLog,offset,mantissa,fmt) {
    if(!digits)digits=6; # same as the default for printf
    if(ABS(logp)<100) {
	fmt = sprintf("%%.%dg",digits);
	return sprintf(fmt, exp(logp));
    } 
    else fmt = sprintf("%%.%dge%%+d",digits);
    l10_p = logp/log(10);
    if(l10_p<0)offset = floor(l10_p);
    intLog=int(l10_p-offset)
    mantissa=10^(l10_p-intLog-offset)
    return sprintf(fmt, mantissa, intLog+offset);
}
function fact(k)    {if(k<=0)return 1; else return k*fact(k-1)}
function logFact(k) { if(k in _memLogFact) return _memLogFact[k];
    if(k<=0)return 0; else return (_memLogFact[k]=log(k)+logFact(k-1));
}
function log2NumSimpleGraphs(n) { # log2(number of non-isomorphic graphs on n nodes), ie. number of bits needed to count them
    # see https://oeis.org/A000088 for comparison for n=0,.. 19
    return choose(n,2)-logFact(n)/log(2) # https://cw.fel.cvut.cz/b211/_media/courses/b4m33pal/lectures/isom_notes.pptx
}
function fact2(k)    {if(k<=1)return 1; else return k*fact2(k-2)}
function logFact2(k) {if(k<=1)return 0; else return log(k)+logFact2(k-2)}
# see Reza expansion: (n k) = ((n-1) (k-1)) + ((n-1) k)
function choose(n,k,     r,i) {if(0<=k&&k<=n){r=1;for(i=1;i<=k;i++)r*=(n-(k-i))/i;} else {r=0; Warn("choose: ("n" choose "k") may not make sense; returning 0")}; return r}
function logChoose(n,k) {
    if(n<k) return log(0); else ASSERT(0<=k && k <=n,"invalid logChoose("n","k")");
    return logFact(n)-logFact(k)-logFact(n-k);
}
function logChooseClever(n,k,     r,i) {
    ASSERT(0<=k&&k<=n,"impossible parameters to logChoose "n" "k)
    if(n in _logChooseMemory && k in _logChooseMemory[n]) return _logChooseMemory[n][k];
    r=0;for(i=1;i<=k;i++)r+=log(n-(k-i))-log(i)
    _logChooseMemory[n][k]=r
    return r;
}
function HalfGamma(k)   {return     sqrt(PI) *   fact2(k-2)/sqrt(2)^(k-1)}
function logHalfGamma(k){return log(sqrt(PI))+logFact2(k-2)-(k-1)*log(sqrt(2))}
function Gamma(x)    {if(x==int(x)) return    fact(x-1);if(2*x==int(2*x))return    HalfGamma(2*x);else ASSERT(0,"Gamma only for integers and half-integers")}
function logGamma(x) {if(x==int(x)) return logFact(x-1);if(2*x==int(2*x))return logHalfGamma(2*x);else ASSERT(0,"Gamma only for integers and half-integers")}

function IncGamma(s,x){
    ASSERT(s==int(s)&&s>=1,"IncGamma(s="s",x="x"): s must be int>=1 for now");
    if(s==1)return Exp(-x)
    else return (s-1)*IncGamma(s-1,x) + x^(s-1)*Exp(-x)
}
function logIncGamma(s,x){
    ASSERT(s==int(s)&&s>=1,"logIncGamma: s must be int>=1 for now");
    if(s==1)return -x;
    else {
	ASSERT(x>0,"logIncGamma: x=" x " must be > 0");
	log_a = log(s-1)+logIncGamma(s-1,x)
	log_c = (s-1)*log(x)-x
	return LogSumLogs(log_a,log_c)
    }
}

# Since gawk cannot pass arrays as parameters, we usurp the global array _Chi2_bins[*][*]. The first index of this array
# is NAME; for a fixed name, the second index is the bins, which are assumed to be equally probable.
function Chi2_stat(name,   bin,X2,avg) { ASSERT(name in _Chi2_bins && isarray(_Chi2_bins[name]), "Chi2_Stat: _Chi2_bins["name"] must be an array of your bin counts");
    _Chi2_n[name]=0; for(bin in _Chi2_bins[name])_Chi2_n[name]+=_Chi2_bins[name][bin];
    avg=_Chi2_n[name]/length(_Chi2_bins[name]);
    X2=0; for(bin in _Chi2_bins[name]) X2+=(_Chi2_bins[name][bin]-avg)^2/avg;
    return X2;
}
function    Chi2_pair2(df,X2){ASSERT(df%2==0,"Chi2_pair2 df "df" must be even"); return    IncGamma(df/2,X2/2) /  Gamma(df/2)}
function logChi2_pair2(df,X2){ASSERT(df%2==0,"Chi2_pair2 df "df" must be even"); return logIncGamma(df/2,X2/2)-logGamma(df/2)}
function    Chi2_pair (df,X2){return df%2==0 ? Chi2_pair2(df,X2) : sqrt(Chi2_pair2(df-1, X2)*Chi2_pair2(df+1,X2))}
function logChi2_pair (df,X2){return df%2==0 ? logChi2_pair2(df,X2) : (logChi2_pair2(df-1, X2)+logChi2_pair2(df+1,X2))/2}
function Chi2_tail_raw(df, x){return  Chi2_pair(df, x)}
function    Chi2_tail(name)  {return  Chi2_pair(length(_Chi2_bins[name]),Chi2_stat(name))}
function logChi2_tail(name){return logChi2_pair(length(_Chi2_bins[name]),Chi2_stat(name))}

function NumBits(n,    b) {b=0;while(n>0){if(n%2==1)b++;n=int(n/2)}; return b}
function log2(n){return log(n)/log(2)}
function log10(n){return log(n)/log(10)}

# res1 is your variable, where the output set goes; it will be nuked and replaced with the set intersection of T1 and T2.
function SetIntersect(res1,T1,T2,
    g){delete res1;if(length(T1)<length(T2)){for(g in T1)if(g in T2)res1[g]=1}
						       else{for(g in T2)if(g in T1)res1[g]=1}}
# same as above but for set union, and res2 is the result.
function SetUnion(res2,T1,T2,
    g){delete res2;for(g in T1)res2[g]=1;for(g in T2)res2[g]=1}
# cumulative add set T to res3
function SetCumulativeUnion(res3,T, g){for(g in T)res3[g]=1}
function SetCopy(dest,src,   g){delete dest;for(g in src)dest[g]=1}

function Jaccard(T1,T2,   i,u){SetIntersect(i,T1,T2); SetUnion(u,T1,T2); return length(i)/length(u);}

# And now counting the info in an edge list. One way to view the info is simply the number of edges.
# Another is to view each nodes adjacency list as having log(n) bits for each of its neighbors.
# This then says that the amount of info is as follows: the end of each edge is listed twice (ie is in two neighbor lists),
# but only one is strictly needed. And each entry is log2(n) bits. So the total info is just log2(n)*numEdges.
# But this number is *way* bigger than the number of edges for all 2018 BioGRID networks, so clearly it is too high.
function netbits(n1,n2,     i,bits){if(n2==0)n2=n1; bits=0;for(i=0;i<MIN(n1,n2);i++)bits+=log(MAX(n1,n2)-i)/log(2); return bits}

function logb(b,x){return log(x)/log(b)}
function dtob(n,   s,sgn) {n=1*n;if(!n)return "0";s=sgn="";if(n<0){sgn="-";n=-n};while(n){s=sprintf("%d%s",(n%2),s); n=int(n/2)}; return sgn s}
function btod(n) {}

function LSPredict(n, x, y, xIn,      SUMx,SUMy,SUMxy,SUMxx,i,slope,y_intercept,x_intercept) {
    SUMx=SUMy=SUMxy=SUMxx=0;
    for(i=0;i<n;i++)
    {
	SUMx += x[i];
	SUMy += y[i];
	SUMxy += x[i]*y[i];
	SUMxx += x[i]*x[i];
    }
    if(n>0 && (SUMx*SUMx - n*SUMxx) != 0) {
	slope = ( SUMx*SUMy - n*SUMxy ) / ( SUMx*SUMx - n*SUMxx );
	y_intercept = ( SUMy - slope*SUMx ) / n;
	x_intercept = BIGNUM;
	if(slope != 0) x_intercept = -y_intercept / slope;
	return slope*xIn + y_intercept;
    }
}

# Queue functions: call QueueAlloc(name) to allocate a queue with name "name"; then Add(name) and Next(name) do the obvious.
function QueueAlloc(name) { _queueFirst[name]=1; _queueLast[name]=0; _queueVals[name][1]=1; delete _queueVals[name][1];}
function QueueDelloc(name) { delete _queueFirst[name]; delete _queueLast[name]; delete _queueVals[name] }
function QueueLength(name) {return _queueLast[name]-_queueFirst[name]+1;}
function QueueAdd(name, val) {_queueVals[name][++_queueLast[name]]=val;}
function QueueNext(name,	val) {
    val =  _queueVals[name][_queueFirst[name]  ];
    delete _queueVals[name][_queueFirst[name]++];
    return val;
}

# if quantiles is true (anything nonzero or nonempty string), remember everyting so we can retrieve quantiles later.
function StatReset(name, quantiles) {
    _statQuantiles[name]=quantiles;
    _statN[name] = _statSum[name] = _statSum2[name] = 0;
    _statMin[name]=BIGNUM;_statMax[name]=-BIGNUM;
    _statmin[name]=BIGNUM;
}

function StatHistReset(name) {
    delete _statHist[name];
    delete _statHistN[name];
    delete _statHistMin[name];
    delete _statHistCDF[name];
    delete _statHistCDFix[name];
}
function StatHistAddSample(name, x) {
    if(!(name in _statHistMin)) _statHistMin[name]=1*BIGNUM;
    if(1*x < _statHistMin[name]) _statHistMin[name]=1*x;
    ++_statHist[name][1*x];
    ++_statHistN[name];
}
function StatHistMakeCDF(name,    n,x,prevX,PMF,prevSort) {
    delete _statHistCDF[name];
    delete _statHistCDFix[name];
    prevX=-1*(BIGNUM); # very very negative number
    n=0; _statHistCDFix[name][0] = prevX;
    prevSort=PROCINFO["sorted_in"];
    PROCINFO["sorted_in"]="@ind_num_asc"; # traverse the array in numerical ascending order by index (ie., x)
    for(x in _statHist[name]) {
	_statHistCDFix[name][++n] = 1*x;
	ASSERT(1*x > 1*prevX, "oops, StatHistMakeCDF found non-incrementing x: "prevX" to "x);
	PMF = _statHist[name][1*x]/(_statHistN[name]);
	_statHistCDF[name][1*x] = _statHistCDF[name][1*prevX] + PMF;
	#printf "_statHistCDF[%s][%g]=%g\n", name, x, _statHistCDF[name][x] >"/dev/stderr";
	prevX = x;
    }
    PROCINFO["sorted_in"]=prevSort;
    # _statHistCDF[name][prevX] may be above 1 due to numerical error; give it some leeway here.
    ASSERT(_statHistCDF[name][prevX]<1+1e-6/_statHistN[name], "_statHistCDF["name"]["prevX"]-1="_statHistCDF[name][prevX]-1" which is too far above 1");
    _statHistCDF[name][prevX]=1;
    # remove the -infinity elements created above
    delete _statHistCDF[name][-1*BIGNUM]; # remove the array element that was created in the first loop above.
    delete _statHistCDFix[name][0];
}

# return the m with x closest to z with x<=z
function StatHistBinarySearch(name,z,    i,n,L,R,m,x) {
    n=length(_statHistCDFix[name]);
    for(i=1;i<=n;i++) ASSERT(i in _statHistCDFix[name], i" is not in F_ix out of "n);
    if(z < _statHistCDFix[name][1]) return 0;
    if(z >= _statHistCDFix[name][n]) return n;
    L=1; R=n;
    while(L < R) {
        m = int((L + R) / 2);
	ASSERT(m>0 && m<=n, "m "m" is out of bounds for n "n" L "L" R "R);
	ASSERT(m in _statHistCDFix[name],"oops, m is "m" out of n="n);
	x=1*_statHistCDFix[name][m];
	ASSERT(x==0|| (x in _statHistCDF[name]), "oops, x "x" is not in _statHistCDF["name"]");
	if(x < z) L = m + 1
        else if(x > z) R = m - 1
        else return m
    }
    # At this point, the value was not found, so return the m just below z
    m = int((L + R) / 2);
    if(m>0 && m<=n) {
	ASSERT(m in _statHistCDFix[name],"oops, m is "m" out of n="n);
	x=_statHistCDFix[name][m];
	ASSERT(x in _statHistCDF[name], "oops, x "x" is not in _statHistCDF["name"]");
	while(m>0 && _statHistCDFix[name][m] > z) --m;
    }
    #printf "FOUND x %g at m %d from n %d L %d R %d\n", x,m,n,L,R
    return m;
}

function StatHistInterpSearch(name,z,    frac,i,n,L,R,m,x) {
    n=length(_statHistCDFix[name]);
    for(i=1;i<=n;i++) ASSERT(i in _statHistCDFix[name], i" is not in F_ix out of "n);
    if(z < _statHistCDFix[name][1]) return 0;
    if(z >= _statHistCDFix[name][n]) return n;
    L=1; R=n; frac=0.5;
    while(L < R) {
	ASSERT(0<=frac && frac<=1,"frac "frac" out of bounds");
        m = int(frac*(L + R));
	ASSERT(m>0 && m<=n, "m "m" is out of bounds for n "n" L "L" R "R);
	ASSERT(m in _statHistCDFix[name],"oops, m is "m" out of n="n);
	x=1*_statHistCDFix[name][m];
	ASSERT(x==0|| (x in _statHistCDF[name]), "oops, x "x" is not in _statHistCDF["name"]");
	if(x < z) L = m + 1
        else if(x > z) R = m - 1
        else return m
	if(_statHistCDFix[name][R] == _statHistCDFix[name][L]) frac=0.5;
	else frac=(x-_statHistCDFix[name][L])/(_statHistCDFix[name][R]-_statHistCDFix[name][L]);
	if(frac<=0 || frac>=1) frac=0.5;
    }
    # At this point, the value was not found, so return the m just below z
    m = int((L + R) / 2);
    if(m>0 && m<=n) {
	ASSERT(m in _statHistCDFix[name],"oops, m is "m" out of n="n);
	x=_statHistCDFix[name][m];
	ASSERT(x in _statHistCDF[name], "oops, x "x" is not in _statHistCDF["name"]");
	while(m>0 && _statHistCDFix[name][m] > z) --m;
    }
    #printf "FOUND x %g at m %d from n %d L %d R %d\n", x,m,n,L,R
    return m;
}

# Return the value in [0,1] of the empirical CDF of [name]
function StatHistECDF(name,z,  n,x,prevX,frac,h1,h2,interp,prevSort,m) {
    z=1*z;
    ASSERT(name in _statHist, "StatHistECDF: no such histogram "name);
    if(!(name in _statHistCDF)) StatHistMakeCDF(name);
    if(z<=_statHistMin[name]) return 0;
    n=length(_statHistCDFix[name]);
    m=StatHistBinarySearch(name,z);
    #printf "z %g i %d x %g\n", z, m, _statHistCDF[name][_statHistCDFix[name][m]]
    # in the following, h1 and h2 are actually x values
    if(m<1) h1=-BIGNUM; else h1=_statHistCDFix[name][m];
    if(m>=n) h2=BIGNUM; else h2=_statHistCDFix[name][m+1];
    frac=(z-h1)/(h2-h1);
    # Now convert the x values to histogram values
    interp=_statHistCDF[name][h1]+frac*(_statHistCDF[name][h2] - _statHistCDF[name][h1]);
    return interp; ######### COMMENT OUT THIS LINE TO CHECK THIS VALUE AGAINST OLD CORRECT CODE BELOW
    prevSort=PROCINFO["sorted_in"];
    PROCINFO["sorted_in"]="@ind_num_asc";
    for(x in _statHistCDF[name]){
	if(1*x>z) {
	    if(m>0) ASSERT(h1==prevX,"m is "m" with x1 "h1" but new is "prevX);
	    if(m>1) ASSERT(h2==    x,"m is "m" with x2 "h2" but new is "x);
	    if(m>0 && m<=n) ASSERT(frac==(z-prevX)/(x-prevX), "frac disagreement "frac" vs "(z-prevX)/(x-prevX));
	    h1=_statHistCDF[name][prevX]; h2=_statHistCDF[name][x];
	    PROCINFO["sorted_in"]=prevSort;
	    ASSERT(interp == h1+frac*(h2-h1), "interp disagreement "interp" vs "h1+frac*(h2-h1));
	    return h1+frac*(h2-h1);
	}
	prevX=x;
    }
    PROCINFO["sorted_in"]=prevSort;
    return 1;
}

# Return the K-S (Kolmogorov-Smirnnov) statistic: the maximum distance between the empirical CDFs of name1 and name2
# FIXME: time is O((n1+n2)^2) since we loop through every value of the histogram, calling ECDF which ALSO does the SAME loop
# It can be done in time O(n1+n2) if we are a bit more clever.
function KSstat(name1,name2,   x,maxD,diff,sign) {
    ASSERT(isarray(_statHist[name1]) && isarray(_statHist[name2]), "KSstat needs stats with histograms");
    maxD=0;
    StatHistMakeCDF(name1);
    StatHistMakeCDF(name2);
    prevX=_statHistMin[name1];
    prevSort=PROCINFO["sorted_in"];
    PROCINFO["sorted_in"]="@ind_num_asc";
    for(x in _statHistCDF[name1]) {
	diff = _statHistCDF[name1][x] - StatHistECDF(name2, x);
	if(ABS(diff) > ABS(maxD)) {sign=1; maxD = diff}
    }
    prevX=_statHistMin[name2];
    for(x in _statHistCDF[name2]) {
	diff = _statHistCDF[name2][x] - StatHistECDF(name1, x);
	if(ABS(diff) > ABS(maxD)) {sign=-1; maxD = diff}
    }
    PROCINFO["sorted_in"]=prevSort;
    return sign*maxD;
}

function KSpvalue(ks_stat,n1,n2, C) {
    C=ks_stat/sqrt((n1+n2)/(n1*n2));
    return 2*Exp(-2*C*C);
}

function KStest(name1, name2) { return KSpvalue(KSstat(name1,name2), _statHistN[name1], _statHistN[name2]); }

function StatQuantile(name,q,   i,which,where,oldWhere,prevSort) {
    ASSERT(_statQuantiles[name], "StatQuantile called on name "name", but _statQuantiles[name] is <"_statQuantiles[name]">");
    ASSERT(0<= q && q<=1, "StatQuantile called with quantile q="q" which is not in [0,1]");
    where=0;
    which=int(q*_statN[name]+0.5);
    #print "StatQuantile called with q="q" on "_statN[name]" elements; which is set to "which
    prevSort=PROCINFO["sorted_in"];
    PROCINFO["sorted_in"]="@ind_num_asc"; # traverse history in numerical order of the indices.
    for(x in _statValue[name]){
	oldWhere=where;
	where += _statValue[name][x];
	if(oldWhere <= which && which <=where) return x;
    }
    PROCINFO["sorted_in"]=prevSort;
}

function StatMedian(name) { return StatQuantile(name,0.5);}
function StatLowerQuartile(name) { return StatQuantile(name,0.25);}
function StatUpperQuartile(name) { return StatQuantile(name,0.75);}

function StatAddSample(name, x) {
    if(1*_statN[name]==0 && !_statQuantiles[name])StatReset(name);
    _statN[name]++;
    _statSum[name]+=x;
    _statSum2[name]+=x*x;
    _statMin[name]=MIN(_statMin[name],x);
    if(x)_statmin[name]=MIN(_statmin[name],x);
    _statMax[name]=MAX(_statMax[name],x);
    if(_statQuantiles[name]) ++_statValue[name][x];
}
function StatAddWeightedSample(name, x, w) {
    if(1*_statN[name]==0)StatReset(name);
    _statN[name]+=w;
    _statNunWtd[name]++;
    _statSum[name]+=w*x;
    _statSum2[name]+=w*x*x;
    _statMin[name]=MIN(_statMin[name],x);
    if(x)_statmin[name]=MIN(_statmin[name],x);
    _statMax[name]=MAX(_statMax[name],x);
}
function StatMean(name) {
    return _statSum[name]/_statN[name];
}
function StatMin(name) {
    return _statMin[name];
}
function Statmin(name) {
    return _statmin[name];
}
function StatMax(name) {
    return _statMax[name];
}
function StatVar(name) {
    if(_statN[name]<2)return 0;
    return (_statSum2[name] - _statSum[name]*_statSum[name]/_statN[name]) / (_statN[name]-1);
}
function StatStdDev(name,     x) {
    x=StatVar(name);if(x<0)return 0;
    return sqrt(StatVar(name));
}
function StatN(name) {
    return _statN[name];
}
function Norm2(n,vec,      i,sum2)
{
    sum2 = 0;
    for(i=0;i<n;i++)sum2+=vec[i]*vec[i];
    return sqrt(sum2);
}
function NormalPtoZ(quantile,    q,z1,n,d)
{
    #printf "NormalPtoZ input is %g\n", quantile
    q = quantile > 0.5 ? (1 - quantile) : quantile;
    #printf "NormalPtoZ q is %g\n", q
    z1 = sqrt (-2.0 * log (q));
    n = (0.010328 * z1 + 0.802853) * z1 + 2.515517;
    d = ((0.001308 * z1 + 0.189269) * z1 + 1.43278) * z1 + 1.0;
    z1 -= n / d;
    return (quantile > 0.5 ? -z1 : z1);
}
function Exp(x){
    if(x < -745) return 5e-324
    else if(x > 707) return BIGNUM7;
    else return exp(x);
}
function NormalPhi(x,    arg)
{
    arg=-x*x/2;
    if(arg<-723) return 1e-314;
    return 0.39894228040143267794*Exp(arg)
}
function NormalDist(mu,sigma,x,  E){E=(mu-x)^2/(2*sigma*sigma);return Exp(-E)/(sigma*2.506628274631000502415765284811)}
function NormalZ2P(x,    b0,b1,b2,b3,b4,b5,t,paren)
{
    if(x<0) return 1-NormalZ2P(-x);
    b0 = 0.2316419; b1 = 0.319381530; b2 = -0.356563782; b3 = 1.781477937; b4 = -1.821255978; b5 = 1.330274429;
    t=1/(1+b0*x);
    paren = t*(b1+t*(b2+t*(b3+t*(b4+t*b5))));
    return NormalPhi(x)*paren;
}
function StatTDistPtoZ(quantile, freedom,        z1,z2,h,x,i)
{
    #printf "StatTDistPtoZ input is (%g,%d)\n",quantile,freedom;
    z1 = ABS(NormalPtoZ(quantile));
    z2 = z1 * z1;
    h[0] = 0.25 * z1 * (z2 + 1.0);
    h[1] = 0.010416667 * z1 * ((5.0 * z2 + 16.0) * z2 + 3.0);
    h[2] = 0.002604167 * z1 * (((3.0 * z2 + 19.0) * z2 + 17.0) * z2 - 15.0);
    h[3] = z1 * ((((79.0 * z2 + 776.0) * z2 + 1482.0) * z2 - 1920.0) * z2 - 945.0);
    h[3] *= 0.000010851;

    x = 0.0;
    for(i = 3; i >= 0; i--)
    x = (x + h[i]) / freedom;
    z1 += x;
    return (quantile > 0.5 ? -z1 : z1);
}
function StatTDistZtoP(t,n,     z,y,b,a)
{
    if(t<0) return 1-StatTDistZtoP(-t,n);
    z=1; t=t*t;y=t/n;b=1+y;
    if(n>=20&&t<n||n>200){
	if(y>1e-6)y=log(b);
	a=n-.5;b=48*a*a;y=a*y;
	y=(((((-.4*y-3.3)*y-24)*y-85.5)/(.8*y*y+100+b)+y+3)/b+1)*sqrt(y);
	return NormalPhi(-y);
    } else {
	a=y=sqrt(y);if(n==1)a=0;
	n-=2;
	while(n>1) {
	    a=a*(n-1)/(b*n)+y;
	    n=-2;
	}
	a=(n==0?a/sqrt(b):(atan2(y,1)+a/b)*0.63661977236);
	return ABS(z-a)/2;
    }
}
function StatConfidenceInterval(name,conf)
{
    return StatTDistPtoZ((1-conf)/2, _statN[name] - 1) * sqrt(StatVar(name) / _statN[name]);
}

function BinomialPMF(p,n,k) { return choose(n,k)*p^k*(1-p)^(n-k)}
function logBinomialPMF(p,n,k) { return logChoose(n,k)+log(p)*k+log(1-p)*(n-k)}
function BinomialCDF(p,n,k, i,sum) {sum=0;
    # Sum terms smallest to largest, along the shortest path to the endpoints since PMF is symmetric.
    if(k>=n/2) for(i=n;i>=k;i--) sum+=BinomialPMF(p,n,i);
    else       for(i=0;i<=k;i++) sum+=BinomialPMF(1-p,n,n-i);
    return sum
}
function logBinomialCDF(p,n,k, i,logSum) {
    # Sum terms smallest to largest, along the shortest path to the endpoints since PMF is symmetric.
    if(k>=n/2){logSum=logBinomialPMF(  p,n,k);for(i=n;i> k;i--) logSum=LogSumLogs(logSum, logBinomialPMF(  p,n,i))}
    else      {logSum=logBinomialPMF(1-p,n,n);for(i=1;i<=k;i++) logSum=LogSumLogs(logSum, logBinomialPMF(1-p,n,n-i))}
    return logSum
}
function Pearson2T(n,r){if(r==1)return BIGNUM; else return r*sqrt((n-2)/(1-r^2))}
# The Poisson1_CDF is 1-CDF, and sums terms smallest to largest; near CDF=1 (ie., 1-CDF=0) it is accurate well below eps_mach.
function PoissonCDF(l,k, sum, term, i){sum=term=1;for(i=1;i<=k;i++){term*=l/i;sum+=term}; return sum*Exp(-l)}
function PoissonPMF(l,k, r,i){if(l>723)return NormalDist(l,sqrt(l),k);r=Exp(-l);for(i=k;i>0;i--)r*=l/i;return r} 
function LogPoissonPMF(l,k, r,i){r=-l;for(i=k;i>0;i--)r+=log(l/i);return r} 
function Poisson1_CDF(l,k, i,sum,psum){psum=-1;sum=0;for(i=k;psum!=sum;i++){psum=sum;sum+=PoissonPMF(l,i)};
    if(sum==0 && k<l) return 1; # this means the numbers are so big the sum got zero but we got less than expected.
    else return sum
}
function LogPoisson1_CDF(l,k, i,sum,pmax,max){pmax=2;max=-BIGNUM;for(i=k;pmax!=max;i++){pmax=max;max=MAX(max,LogPoissonPMF(l,i))};
    if(max==1 && k<l) return 0; # this means the numbers are so big the sum got zero but we got less than expected.
    else return max/.894
}

# Return the logarithm of the area under the tail of the Normal(0,1) distribution from -infinity up to z
function logPhi(z){
    # The first is derived by L hopitals rule, and returns the log of the tail of the normal distribution at z.
    # At a p-value of 0.01, the error is about 3%, and the percent error drops by a factor of about 1.5 per factor of 10
    # decrease in p-value: so it is about 1.3% at p=1e-3; 0.7% at 1e-4; 0.43% at 1e-5; 0.3% at 1e-6; 0.22% at 1e-7, etc.
    if(z<-1.40532)return -(log(sqrt(2*PI))+z*z/2+log(-z));
    #The one below is a quadratic fit in the range sigma [-1.4,1] (ie., p-values in the range about 0.05 up to .8),
    # and can be off by a few tens of percent but for p-values so large... do we really care?
    else if(z<1.00268) return (-0.976263+((z+1.40532)/2.408)^(0.75)*(-.0686+0.9763))*log(10)
    # Otherwise this is an excellent approximation for large z
    else return -NormalZ2P(z)
}

function Log10Poisson1_CDF(l,k, i,sum,psum){return LogPoisson1_CDF(l,k, i,sum,psum)/2.302585092994046}

# Hypergeometric distribution: Given: total population N, K of which have desired property.
# What is the probability of exactly k successes in n draws, without replacement?
function HyperGeomPMF(k,n,K,N){return Exp(logHyperGeomPMF(k,n,K,N))}
function logHyperGeomPMF(k,n,K,N) {
    #print "logHyperGeomPMF",k,n,K,N
    return logChoose(K,k)+logChoose(N-K,n-k)-logChoose(N,n);
    #return logChoose(n,k)+logChoose(N-n,K-k)-logChoose(N,K);
}
function HyperGeomTail(k,n,K,N, sum,term,i,logTerm) {
    #print "HyperGeomTail",k,n,K,N
    # 166 1113 2141 3436
    ASSERT(k<=K && k<=n && K<=N && n<=N && n-k<=N-K && K-k<=N-n,"HyperGeomTail: impossible values "k"/"K","n"/"N)
    if(k==0 && K>0) return 1;
    if(k in _hyperGeomMem && n in _hyperGeomMem[k] && K in _hyperGeomMem[k][n] && N in _hyperGeomMem[k][n][K])
	return _hyperGeomMem[k][n][K][N]
    logTerm = logHyperGeomPMF(k,n,K,N)
    sum = term = Exp(logTerm)
    for(i=k+1; (i<=n && i<=K && sum < 0.999999) && (logTerm<723 || (sum && (term/sum > 1e-20))); i++) {
	#print i,logTerm,term,sum,sum-1
	logTerm = logHyperGeomPMF(i,n,K,N)
	term = Exp(logTerm)
	sum += term
    }
    #print "DONE!!!",sum
    if(sum==0)sum=1e-320
    _hyperGeomMem[k][n][K][N] = sum
    return sum
}

function logHyperGeomTail(k,n,K,N, logSum,logTerm,i) {
    #print "logHyperGeomTail",k,n,K,N
    ASSERT(k<=K && k<=n && K<=N && n<=N && n-k<=N-K && K-k<=N-n,"logHyperGeomTail: impossible values "k"/"K","n"/"N)
    if(k==0 && K>0) return 0;
    if(k in _logHyperGeomMem && n in _logHyperGeomMem[k] && K in _logHyperGeomMem[k][n] && N in _logHyperGeomMem[k][n][K])
	return _logHyperGeomMem[k][n][K][N]
    logTerm = logHyperGeomPMF(k,n,K,N)
    logSum = logTerm
    for(i=k+1; i<=n && i<=K && logTerm-logSum > -40; i++) {
	logTerm = logHyperGeomPMF(i,n,K,N)
	logSum  = LogSumLogs(logSum,logTerm)
    }
    _logHyperGeomMem[k][n][K][N] = logSum
    return logSum
}

function StatRV_Normal(){if(!_StatRV_which) { do { _StatRV_v1 = 2*rand()-1; _StatRV_v2 = 2*rand()-1; _StatRV_rsq = _StatRV_v1^2+_StatRV_v2^2; } while(_StatRV_rsq >= 1 || _StatRV_rsq == 0); _StatRV_fac=sqrt(-2*log(_StatRV_rsq)/_StatRV_rsq); _StatRV_next = _StatRV_v1*_StatRV_fac; _StatRV_which = 1; return _StatRV_v2*_StatRV_fac; } else { _StatRV_which = 0; return _StatRV_next; } }
function NormalRV(mu,sigma) { return mu+StdNormRV()*sigma }

# The Spearman correlation is just the Pearson correlation of the rank. It measures monotonicity, not linearity.
# Unfortunately it means we need to store every sample, and sort them into rank order when we want the coefficient.
function SpearmanAddSample(name,X,Y) {
    delete _SpComputeResult[name];
    _SpN=(_Spearman_N[name]++) # 1-indexed, not zero.
    _SpearmanSampleX[name][_SpN]=X;
    _SpearmanSampleY[name][_SpN]=Y;
}
function SpearmanCompute(name, i) {
    ASSERT(name in _Spearman_N, "SpearmanCompute: no such data "name);
    if(name in _SpComputeResult) return _SpComputeResult[name];
    ASSERT(length(_SpearmanSampleX[name])==length(_SpearmanSampleY[name]), "SpearmanCompute: input arrays are different lengths");
    # Too hard to do this in awk, just run external spearman program
    _SpCommand = "spearman"
    for(i=1;i<=_Spearman_N[name];i++) print _SpearmanSampleX[name][i],_SpearmanSampleY[name][i] |& _SpCommand;
    close(_SpCommand,"to");
    _SpCommand |& getline _SpComputeResult[name]
    close(_SpCommand,"from");
    n=split(_SpComputeResult[name],a);
    ASSERT(a[1]==_Spearman_N[name],"SpearmanCompute: first field returned by external command "_SpCommand" is not _Spearman_N["name"]="_Spearman_N[name]);
    _Spearman_rho[name]=a[2];
    _Spearman_p[name]=a[3];
    _Spearman_t[name]=a[4];
    delete a
    return _SpComputeResult[name];
}
function SpearmanPrint(name) { return SpearmanCompute(name) }

function CovarReset(name) {
    delete _Covar_sumX[name]
    delete _Covar_sumY[name]
    delete _Covar_sumXY[name]
    delete _Covar_N[name]
}
function CovarAddSample(name,X,Y) {
    _Covar_sumX[name]+=X
    _Covar_sumY[name]+=Y
    _Covar_sumXY[name]+=X*Y
    _Covar_N[name]++;
}

function CovarCompute(name){
    ASSERT(1*_Covar_N[name]>1, "CovarCompute requires N>=1 but it is "_Covar_N[name]);
    return (_Covar_sumXY[name]-_Covar_sumX[name]*_Covar_sumY[name]/_Covar_N[name])/(_Covar_N[name]-1);
}

function PearsonReset(name) {
    delete _Pearson_sumX[name]
    delete _Pearson_sumY[name]
    delete _Pearson_sumXY[name]
    delete _Pearson_sumX2[name]
    delete _Pearson_sumY2[name]
    delete _Pearson_N[name]
    delete _Pearson_rho[name]
    delete _Pearson_t[name]
    delete _Pearson_p[name]
}
function PearsonAddSample(name,X,Y) {
    _PearsonComputeValid[name]=0;
    _Pearson_sumX[name]+=X
    _Pearson_sumY[name]+=Y
    _Pearson_sumXY[name]+=X*Y
    _Pearson_sumX2[name]+=X*X
    _Pearson_sumY2[name]+=Y*Y
    _Pearson_N[name]++;
}

function PearsonCompute(name,     numer,DX,DY,denom,z,zse,F){
    if(!_Pearson_N[name])return 0;
    if(_PearsonComputeValid[name]) return 1;
    numer=_Pearson_sumXY[name]-_Pearson_sumX[name]*_Pearson_sumY[name]/_Pearson_N[name]
    DX=_Pearson_sumX2[name]-_Pearson_sumX[name]*_Pearson_sumX[name]/_Pearson_N[name]
    DY=_Pearson_sumY2[name]-_Pearson_sumY[name]*_Pearson_sumY[name]/_Pearson_N[name]
    #print DX,DY >"/dev/stderr"
    denom=sqrt(ABS(DX*DY)); # ABS since sometimes it is very slightly negative due to rounding errors
    _Pearson_rho[name]=0; if(denom)_Pearson_rho[name]=numer/denom;
    _Pearson_t[name]=Pearson2T(_Pearson_N[name],_Pearson_rho[name]);
    if(_Pearson_t[name]<0)_Pearson_t[name]=-_Pearson_t[name];
    # Fisher R-to-z
    z=0.5*log((1+_Pearson_rho[name])/(1-_Pearson_rho[name]))
    zse=1/sqrt(ABS(_Pearson_N[name]-3))
    _Pearson_p[name]=F=2*MIN(NormalDist(0,zse,z),NormalDist(0,zse,-z))
    # We seem to be at least 100x too small according to Fisher
    if(_Pearson_p[name]>1)_Pearson_p[name]=1-1/_Pearson_p[name]
    _PearsonComputeValid[name]=1;
    return 1
}

function PearsonPrint(name, logp){
    #if(!_Pearson_N[name]) return;
    PearsonCompute(name);
    TINY=1e-200; # using the fancy log algorithm if p-value is smaller than this
    logp = -logPhi(-_Pearson_t[name]); # working with the negative log is easier (so log is positive)
    if(logp < -log(TINY))
	return sprintf("%d\t%.4g\t%.4g\t%.4f", _Pearson_N[name], _Pearson_rho[name], _Pearson_p[name], _Pearson_t[name])
    else {
	#printf "t %g p %g log10p %g logp %g", _Pearson_t[name], _Pearson_p[name], logp/log(10), logp > "/dev/stderr"
	logp = (logp - 8.28931 - logp/65.1442)/0.992 # Empirical correction to get in line with Fisher for small p-values
	#printf " (logp corrected %g %g)\n", logp/log(10), logp > "/dev/stderr"
	return sprintf("%d\t%.4g\t%s\t%.4f (using log)", _Pearson_N[name], _Pearson_rho[name], logPrint(-logp,4), _Pearson_t[name]);
	#p=10^-logp; print "log-over-Fisher", p/F # Sanity check
    }
}

# Functions for computing the AUPR
function AUPR_add(name, value, thresh, truth){
    _AUPR_N[name]++;
    if(value > thresh){ # predicted
	if(truth)++_AUPR_TP[name];else ++_AUPR_FP[name]
    }else{ # not predicted
	if(truth)++_AUPR_FN[name];else ++_AUPR_TN[name]
    }
}
function AUPR_Prec(name,  TP, FP){TP=_AUPR_TP[name];FP=_AUPR_FP[name];if(TP+FP==0)return 1; else return TP/(TP+FP)}
function AUPR_Rec(name,     TP, FN){TP=_AUPR_TP[name];FN=_AUPR_FN[name];if(TP+FN==0)return 1; else return TP/(TP+FN)}
function AUPR_F1(name,         Prec,Rec){Prec=AUPR_Prec(name); Rec=AUPR_Rec(name);
    if(Prec+Rec==0)return 0; else return 2*Prec*Rec/(Prec+Rec)
}
function AUPR_TPR(name,  TP, FN){TP=_AUPR_TP[name];FN=_AUPR_FN[name]; return TP/(TP+FN)}
function AUPR_FPR(name,  FP, TN){FP=_AUPR_FP[name];TN=_AUPR_TN[name]; return FP/(FP+TN)}


# The method of least squares is a standard technique used to find
#  the equation of a straight line from a set of data. Equation for a
#  straight line is given by 
#	 y = mx + b
#  where m is the slope of the line and b is the y-intercept.
#
#  Given a set of n points {(x1,y1), x2,y2),...,xn,yn)}, let
#      SUMx = x1 + x2 + ... + xn
#      SUMy = y1 + y2 + ... + yn
#      SUMxy = x1*y1 + x2*y2 + ... + xn*yn
#      SUMxx = x1*x1 + x2*x2 + ... + xn*xn
#
#  The slope and y-intercept for the least-squares line can be 
#  calculated using the following equations:
#        slope (m) = ( SUMx*SUMy - n*SUMxy ) / ( SUMx*SUMx - n*SUMxx ) 
#  y-intercept (b) = ( SUMy - slope*SUMx ) / n
# AUTHOR: Dora Abdullah (Fortran version, 11/96)
# REVISED: RYL (converted to C, 12/11/96)
# Converted to awk: Wayne Hayes (2020-March-09)

function LeastSquaresReset(x,y) {_LS_SUMx=_LS_SUMy=_LS_SUMxy=_LS_SUMxx=_LS_n=0; delete _LS_x; delete _LS_y}
function LeastSquaresSample(x,y) {_LS_SUMx+=x; _LS_SUMy+=y; _LS_SUMxy+=x*y;_LS_SUMxx+=x*x; _LS_x[_LS_n]=x; _LS_y[_LS_n]=y; ++_LS_n}
function LeastSquaresSlope(){_LS_slope=(_LS_SUMx*_LS_SUMy - _LS_n*_LS_SUMxy )/( _LS_SUMx*_LS_SUMx - _LS_n*_LS_SUMxx)
    return _LS_slope
}
function LeastSquares_y_intercept() { return ( _LS_SUMy - LeastSquaresSlope()*_LS_SUMx ) / _LS_n}
function LeastSquaresMSR(  i) {
  SUMres = 0;
  SUMres2 = 0;
  slope = LeastSquaresSlope();
  y_intercept = LeastSquares_y_intercept();
  for (i=0; i<_LS_n; i++) {
    y_estimate = slope*_LS_x[i] + y_intercept;
    res = _LS_y[i] - y_estimate;
    SUMres += res;
    SUMres2 += res*res;
  }
  return (SUMres2)/_LS_n;
  variance = (SUMres2 - SUMres*SUMres/_LS_n)/(_LS_n-1);
}


################# GRAPH ROUTINES ##################

#Input: edgeList; a single node, u, to start the BFS; and an (optional) "searchNode" to stop at.
#Output: array dist[] contains shortest paths from u to all nodes reachable from u within maxDist; includes dist[u]=0.
#        Call with maxDist=n (size of network) to get the BFS distance to everybody
function BFS(edgeList,u,searchNode,dist,   V,Q,m,M,x,y) {
    ASSERT(isarray(edge), "BFS: edgeList must be binary symmetric 2D array");
    delete V; # visited
    delete Q; # queue
    delete dist; # distance from u
    dist[u]=0;
    m=M=0;
    Q[M++]=u; # the BFS queue runs from m [inclusive] to M-1, and we increment m as we dequeue elements
    while(M>m) {
	x = Q[m++];
	ASSERT(x in dist, x" in Q but not in distance array");
	if(!(x in V)) {
	    V[x]=1;
	    ASSERT(isarray(edge[x]), "edge["x"] is not an array");
	    for(y in edge[x]) if(!(y in V)) {
		if(y in dist) dist[y]=MIN(dist[y],dist[x]+1);
		else dist[y]=dist[x]+1;
		if(y==searchNode) return;
		Q[M++]=y;
	    }
	}
    }
}
function MakeEmptySet(S){delete S; S[0]=1; delete S[0]}
function InducedEdges(edge,T,D,       u,v,m) { # note you can skip passing in D
    MakeEmptySet(D);
    for(u in T) for(v in T) if((u in edge) && (v in edge[u])) { ++D[u]; ++D[v]; ++m; }
    for(u in T) { ASSERT(D[u]%2==0, "InducedEdges: D["u"]="D[u]); D[u]/=2; }
    ASSERT(m%2==0, "m is not even");
    return m/2;
}
function InducedWeightedEdges(edge,T,D,       u,v,m,all1) { # note you can skip passing in D
    MakeEmptySet(D); all1=1;
    for(u in T) for(v in T) if((u in edge) && (v in edge[u])) {
	if(edge[u][v] != 1) all1 = 0;
	D[u]+=edge[u][v]; D[v]+=edge[u][v]; m+=edge[u][v];
    }
    for(u in T) { if(all1) ASSERT(D[u]%2==0, "InducedEdges: D["u"]="D[u]); D[u]/=2; }
    if(all1) ASSERT(m%2==0, "m is not even");
    return m/2;
}

# Note: Possible sort orders are: "@unsorted",
# "@ind_str_asc",	"@ind_num_asc",	 "@val_type_asc",  "@val_str_asc",  "@val_num_asc",
# "@ind_str_desc",	"@ind_num_desc", "@val_type_desc", "@val_str_desc", "@val_num_desc",
# This implementation allows multiple elements with the same priority... and even multiple [p][element] duplicates
function PQpush(name, pri, element) { ++_PQ_[name][pri][element]; _PQ_size[name]++ }

function PQpop(name,    prevSort, element, p) {
    prevSort=PROCINFO["sorted_in"]; # remember sort order to restore it afterwards
    PROCINFO["sorted_in"]="@ind_num_desc";
    for(p in _PQ_[name]) {
        # Note that if multiple elements have the same priority, we will return them in SORTED order not INSERTION order
        for(element in _PQ_[name][p]) {
            if(--_PQ_[name][p][element]==0) {
                delete _PQ_[name][p][element];
                if(length(_PQ_[name][p])==0) delete _PQ_[name][p];
            }
            break; # exit at first iteration
        }
        break; # exit at first iteration
    }
    PROCINFO["sorted_in"]=prevSort; # restore sort order
    _PQ_size[name]--
    return element;
}

function PQlength(name) { return _PQ_size[name]; }
function PQalloc(name) { _PQ_size[name]=0;PQ_[name][0][0]=1; delete PQ_[name][0] }
function PQdelloc(name) { delete PQ_size[name]; delete PQ_[name] }
function PQfree(name) { PQdelloc(name); }
