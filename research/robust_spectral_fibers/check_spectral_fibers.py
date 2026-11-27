"""Exact synthetic checks for spectral dispersion; no curve/contractor code."""
from fractions import Fraction as F
from itertools import product, combinations
from collections import Counter
from math import factorial, isqrt
from random import Random


def variance_num(xs):
    n=len(xs)
    return n*sum(x*x for x in xs)-sum(xs)**2


def minvar(xs,r):
    """All cut patterns occur at an occupied residue (interval half-open)."""
    return F(min(variance_num([((x-a)%r) for x in xs]) for a in set(xs)),len(xs)**2*r*r)


def atan_bounds(d,n=24):
    s=sum((F((-1)**j,(2*j+1)*d**(2*j+1)) for j in range(n)),F(0))
    t=s+F((-1)**n,(2*n+1)*d**(2*n+1))
    return min(s,t),max(s,t)


def pi_bounds():
    a,b=atan_bounds(5); c,d=atan_bounds(239)
    return 16*a-4*d,16*b-4*c


Q=1<<48

def cos_table(r):
    plo,phi=pi_bounds(); p=(plo+phi)/2
    out=[]
    for j in range(r):
        if j==0:
            out.append((Q,Q));continue
        x=2*p*F(min(j,r-j),r)
        s=sum(((-1)**k*x**(2*k)/factorial(2*k) for k in range(21)),F(0))
        err=x**42/factorial(42)+(phi-plo)*F(min(j,r-j),r)
        lo=max(F(-1),s-err);hi=min(F(1),s+err)
        out.append(((lo*Q).numerator//(lo*Q).denominator,-((-hi*Q).numerator//(-hi*Q).denominator)))
    return out


def sq_bounds(hist,cos):
    r=len(hist);den=sum(hist)**2*Q
    assert den
    a=[sum(hist[x]*hist[(x+s)%r] for x in range(r)) for s in range(r)]
    lo=F(sum(a[s]*cos[s][0] for s in range(r)),den)
    hi=F(sum(a[s]*cos[s][1] for s in range(r)),den)
    return max(F(0),lo),min(F(1),hi)


def sqrt_upper(f,den=1<<24):
    """Smallest rational of the given denominator >= sqrt(f)."""
    assert f>=0
    n=isqrt((f.numerator*den*den)//f.denominator)
    if F(n*n,den*den)<f:n+=1
    return F(n,den)


def conv1(a,b):
    out=Counter()
    for x,c in a.items():
        for y,d in b.items():out[x+y]+=c*d
    return out


def diffs(a):return conv1(a,{-x:c for x,c in a.items()})


def moment(hist,degree):
    return F(sum(c*x**degree for x,c in hist.items()),sum(hist.values()))


def even_moments(hist,H):
    # Iterative powers avoid repeated large exponentiation.
    xs=list(hist); weights=[hist[x] for x in xs];sq=[x*x for x in xs]
    powers=[1]*len(xs); out=[]
    for j in range(H//2+1):
        out.append(sum(c*y for c,y in zip(weights,powers)))
        powers=[y*z for y,z in zip(powers,sq)]
    return out


def p_expect_bounds(hist,r,H):
    assert H%4==0
    plo,phi=pi_bounds();lo=F(0);hi=F(0)
    moments=even_moments(hist,H); mass=sum(hist.values())
    for j in range(1,H//2+1):
        k=2*j; c=F(moments[j],mass*r**k*factorial(k))
        a=c*(2*plo)**k;b=c*(2*phi)**k
        if j%2:lo+=a;hi+=b
        else:lo-=b;hi-=a
    return lo,hi


def variance_tests():
    rng=Random(45021);checks=0;cut_checks=0
    for r in (5,7,11,13,17,31):
        cos=cos_table(r)
        for _ in range(90):
            xs=[rng.randrange(r) for _ in range(rng.randrange(1,13))]
            v=minvar(xs,r)
            assert v<=F(1,4)
            # Theoretical V <= (1-|M|)/8. Verify with certified upper |M|.
            hist=[xs.count(a) for a in range(r)]
            _,ub=sq_bounds(hist,cos)
            # Avoid rounding issues in a point mass: use exact 8V<=1-|M|,
            # equivalently |M|^2 <= (1-8V)^2 when V<=1/8.
            assert 1-8*v>=0
            assert ub <= (1-8*v)**2
            checks+=1
            # Independently check all possible cuts, not just occupied ones.
            ref=F(min(variance_num([(x-a)%r for x in xs]) for a in range(r)),len(xs)**2*r*r)
            assert v==ref
            cut_checks+=r
    print(f'Circular variance: {checks} certified phase/variance bounds; {cut_checks} independent cut trials.')
    for h in (256,257,511,512,1024,4096,65536):
        # h >= 256 makes this conservative rational remainder guard tiny.
        assert F(5000*h,1<<(2*h//15))<F(1,4)
        m=max(h+1,2*h)
        # Upper bound on perturbation ratio after replacing (7/100)^h by 1.
        assert F(6*m*m,(100*m)**3)<F(1,4)
    print('Proof constants: seven exact conservative remainder and perturbation guards.')


def separation_test():
    # Abstract labels, not actual affine-mask output. All projective combinations
    # of (lambda,u) are tested, and arbitrary scaling preserves the cardinality
    # lower bound on width. This is a test of necessary conditions, not a kernel.
    r=1009;m=1024;h=256
    slopes=[1+i%(r-1) for i in range(m)]
    # Group values can all be nonidentity by choosing a reference offset that
    # is not in {0,-s,-2s}; block reference subtraction makes lambda=s*j.
    for s in slopes:
        ref=next(t for t in range(1,4) if all((t+s*j)%r for j in range(3)))
        assert len({(ref+s*j)%r for j in range(3)})==3
    min_modal=m;minwidth=10**9
    for a,c in [(1,c) for c in range(r)]+[(0,1)]:
        nonconstant=sum((a*s+c)%r!=0 for s in slopes)
        # Three distinct residues in every nonconstant block; a retained
        # 3/4-fraction of a three-point block must retain all three.
        E=F(2,3)*nonconstant
        width=2*(nonconstant-h//4)
        assert E>F(h,16) and width>r
        min_modal=min(min_modal,E);minwidth=min(minwidth,width)
    v=m*minvar([0,1,2],r)
    assert 0<v<F(h,60000)
    # Exact high-degree minorant on a wide-lift low-variance independent law.
    # Uniform labels 0,1,2; total histogram via ordinary convolution.
    hist=Counter({0:1})
    for _ in range(m):
        new=Counter()
        for x,c in hist.items():
            new[x]+=c;new[x+1]+=c;new[x+2]+=c
        hist=new
    # For the minorant it is enough to check order four rigorously here;
    # the high-order identity is separately tested on exact OA controls.
    d=diffs(hist);lo,hi=p_expect_bounds(d,r,4)
    assert lo>0
    print(f'Wide-lift separation: r={r}, m={m}, h={h}, {r+1} projective slopes; modal deficit >= {min_modal}, trimmed-width lower bound {minwidth} > {r}, variance {v} < {F(h,60000)}. Certified positive degree-four independent minorant.')


def high_moment_controls():
    # A modular parity array: each half is uniform on a sum-zero fiber of
    # F3, length h+1. The two simultaneous constraints leave exact h-wise
    # marginals but many high-variance directions. Integer histogram DP
    # proves the paired moments agree through h without enumerating 3^514 rows.
    r=3;h=256;n=h+1
    half=Counter({0:1})
    for _ in range(n):half=conv1(half,{-1:1,0:1,1:1})
    restricted=Counter({x:c for x,c in half.items() if x%r==0})
    assert sum(restricted.values())*r==sum(half.values())
    # u uses + in one half, - in other; symmetry of each half permits same histogram.
    assert restricted==Counter({-x:c for x,c in restricted.items()})
    ind=diffs(conv1(half,half));cond=diffs(conv1(restricted,restricted))
    I=even_moments(ind,h);C=even_moments(cond,h)
    for a,b in zip(I,C):assert a*C[0]==b*I[0]
    assert all(x%r==0 for x in cond)
    v=2*n*minvar([-1%r,0,1],r)
    assert v>F(h,60000)
    print(f'Exact two-constraint OA: m={2*n}, h={h}, {len(I)} paired even-moment identities through degree {h}; V={v} above the exclusion threshold. Abstract conditional array, not a full curve fiber.')


def approx_marginal_tests():
    rng=Random(45022);checks=0;marg=0
    for _ in range(32):
        m=5;r=17;h=4;delta=F(1,10**8)
        labels=[[rng.randrange(r) for _ in range(2)] for i in range(m)]
        words=list(product(range(2),repeat=m));target=rng.choice(words)
        # Mix a tiny atom into uniform; all marginal multiplicative errors <=delta.
        den=10**8*len(words)**2; base=den//len(words)
        # uniform component weight den-len(words), plus atom len(words)
        weights={x:base-1+len(words)*(x==target) for x in words}
        assert sum(weights.values())==den
        for l in range(h+1):
            for J in combinations(range(m),l):
                hits=Counter()
                for x,c in weights.items():hits[tuple(x[j] for j in J)]+=c
                for a in product(range(2),repeat=l):
                    assert abs(F(hits[a],den)*2**l-1)<=delta
                    marg+=1
        a=Counter();b=Counter()
        for x in words:
            y=sum(labels[i][x[i]] for i in range(m));a[y]+=1;b[y]+=weights[x]
        v=sum(F((s[0]-s[1])**2,4*r*r) for s in labels)
        aa=diffs(a);bb=diffs(b)
        for j in range(1,h//2+1):
            ma=moment(aa,2*j)/r**(2*j);mb=moment(bb,2*j)/r**(2*j)
            assert abs(ma-mb)<=6*delta*m**(2*j-1)*v
            checks+=1
    print(f'Approximate pair-law perturbations: {marg} exact marginal tests; {checks} even-moment error bounds.')


def conv2(a,b,r):
    out=Counter()
    for (g,s),c in a.items():
        for (h,t),d in b.items():out[((g+h)%r,(s+t)%r)]+=c*d
    return out


def projection(a,r):
    p=[0]*r
    for (g,s),c in a.items():p[g]+=c
    return p


def char_hist(a,r,x,y):
    out=[0]*r
    for (g,s),c in a.items():out[(x*g+y*s)%r]+=c
    return out


def corr(hist):
    r=len(hist);return [sum(hist[x]*hist[(x+j)%r] for x in range(r)) for j in range(r)]


def mixing_tests():
    rng=Random(45023);cases=targets=biaschecks=parseval=failchecks=pins=0;strong=0
    for r in (5,7,11,13):
        cos=cos_table(r)
        for trial in range(18):
            m=5
            lists=[]
            for i in range(m):
                # Distinct group values in each abstract indexed block.
                gs=rng.sample(range(r),rng.randrange(max(2,r//2),r+1))
                lists.append(Counter({(g,rng.randrange(r)):1 for g in gs}))
            def build(retain):
                out=Counter({(0,0):1})
                for i in retain:out=conv2(out,lists[i],r)
                return out
            A=build(range(m));mass=sum(A.values());p=projection(A,r)
            chi=F(r*sum(x*x for x in p),mass*mass)-1
            if chi>=F(1,4):continue
            # Exact Parseval in the cyclotomic group ring: two coefficient
            # arrays represent the same value at zeta_r iff their difference
            # is constant (multiple of 1+z+...+z^(r-1)).
            best=F(0)
            for b in range(1,r):
                lhs=[0]*r;rhs=[0]*r
                for a in range(r):
                    h=char_hist(A,r,a,b)
                    _,ub=sq_bounds(h,cos);best=max(best,ub)
                    for j,x in enumerate(corr(h)):lhs[j]+=x
                for g in range(r):
                    h=[0]*r
                    for s in range(r):h[(b*s)%r]+=A[g,s]
                    for j,x in enumerate(corr(h)):rhs[j]+=r*x
                assert len({x-y for x,y in zip(lhs,rhs)})==1
                assert rhs[0] <= r*sum(x*x for x in p)
                parseval+=1
            rho=sqrt_upper(best)
            if rho>=1:continue
            cases+=1
            B=conv2(A,A,r);C=conv2(A,B,r);D=conv2(A,C,r)
            for l,dist,T in ((3,C,B),(4,D,C)):
                eta=(1+chi)/(1-chi)*rho**(l-2)
                pp=projection(dist,r);tot=sum(pp)
                good=[F(p[g],mass)>=F(1,2*r) for g in range(r)]
                beta=F(good.count(False),r)
                bound=(1+chi)*beta+2*(1+chi)*(F(1,r)+eta)
                for t in range(r):
                    assert 1-chi<=F(r*pp[t],tot)<=1+chi
                    row=[dist[t,s] for s in range(r)]
                    assert max(F(n,pp[t]) for n in row)<=F(1,r)+eta
                    for b in range(1,r):
                        hh=[row[(j*pow(b,-1,r))%r] for j in range(r)]
                        _,ub=sq_bounds(hh,cos)
                        assert ub<=eta**2
                        biaschecks+=1
                    targets+=1
                    # Any actual reference score is a valid comparison value.
                    for ref in range(r):
                        fail=0
                        for (g,s),c in T.items():
                            v=(t-g)%r
                            if not good[v] or A[v,(ref-s)%r]==p[v]:fail+=c
                        assert F(fail,sum(T.values()))<=bound
                        failchecks+=1
                    if F(1,r)+eta<F(1,2):strong+=1
            # Heterogeneous/pinned copy bound: remove one slot from each of
            # three copies; arbitrary fixed choices contribute constant shifts.
            I=[build([i for i in range(m) if i!=j]) for j in range(3)]
            ps=[projection(x,r) for x in I];ms=[sum(z) for z in ps]
            e=max(F(r*sum(x*x for x in z),mass*mass)-1 for z,mass in zip(ps,ms))
            if e>=F(1,4):continue
            rho2=F(0)
            for x in I:
                for a in range(r):
                    for b in range(1,r):rho2=max(rho2,sqrt_upper(sq_bounds(char_hist(x,r,a,b),cos)[1]))
            if rho2>=1:continue
            C=conv2(conv2(I[0],I[1],r),I[2],r);eta2=(1+e)/(1-e)*rho2
            for t in range(r):
                row=[C[t,s] for s in range(r)];s=sum(row)
                assert (1-e)<=F(r*s,sum(C.values()))<=1+e
                assert F(max(row),s)<=F(1,r)+eta2
                for b in range(1,r):
                    hh=[row[(j*pow(b,-1,r))%r] for j in range(r)]
                    assert sq_bounds(hh,cos)[1]<=eta2**2
                    pins+=1
    assert cases>0 and strong>0 and pins>0
    # Deliberate affine-only negative: no independent score direction.
    r=7;cos=cos_table(r);A=Counter({(g,g):1 for g in range(r)})
    assert sqrt_upper(sq_bounds(char_hist(A,r,-1,1),cos)[1])==1
    print(f'Conditional spectral mixing: {cases} synthetic laws, {parseval} exact cyclotomic Parseval identities, {targets} target/arity pairs, {biaschecks} certified bias inequalities, {failchecks} exact translation-failure bounds, {strong} target/arity pairs with atom bound <1/2; {pins} heterogeneous pin-copy bias checks. Affine-only rho=1 negative control retained.')


if __name__=='__main__':
    variance_tests()
    separation_test()
    high_moment_controls()
    approx_marginal_tests()
    mixing_tests()
