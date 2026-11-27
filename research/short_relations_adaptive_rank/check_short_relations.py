"""Exact finite tests for short relation and conditional-mean certificates.
These enumerate toy tuples, not the asymptotic contractor. Python 3.10+.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import prod
from random import Random


def rank(rows, r, n):
    basis = {}
    for row in rows:
        v = [x % r for x in row]
        for j in sorted(basis):
            if v[j]:
                a = v[j]
                v = [(x-a*y) % r for x,y in zip(v,basis[j])]
        k = next((j for j,x in enumerate(v) if x), None)
        if k is not None:
            t = pow(v[k],-1,r)
            basis[k] = [x*t % r for x in v]
    return len(basis)


def short_fixture(r, blocks, cycle=None, add=None):
    """Count pinned mass and incidence vectors at every displacement."""
    sizes = list(map(len,blocks)); starts = [sum(sizes[:i]) for i in range(len(sizes))]
    vals = [a for B in blocks for a in B]; n = len(vals); m = len(blocks)
    ix = [range(j,j+k) for j,k in zip(starts,sizes)]
    T = list(product(*ix)); W = len(T)
    if cycle is None:
        value = lambda X: sum(vals[j] for j in X) % r
        sub = lambda y,x: (y-x)%r
        disps = list(range(r))
    else:
        def value(X):
            S = None
            for j in X: S = add(S,cycle[vals[j]])
            return S
        def sub(y,x):
            return add(y,None if x is None else (x[0],(-x[1])%13))
        disps = cycle
    tv = [(X,value(X)) for X in T]
    pins = [(sgn,j) for sgn in (-1,1) for j in range(n)]
    H = {g:{pin:[0,[0]*n] for pin in pins} for g in disps}
    d = Counter(g for _,g in tv)
    for X,x in tv:
        for Y,y in tv:
            t = sub(y,x)
            a = [0]*n
            for j in X: a[j] -= 1
            for j in Y: a[j] += 1
            for pin in [(-1,j) for j in X]+[(1,j) for j in Y]:
                h = H[t][pin]; h[0]+=1
                h[1]=[u+v for u,v in zip(h[1],a)]
    # Verify conditional mean rows annihilate true logs and block constants.
    rows = []; by_t = []
    for t in disps:
        tab=H[t]
        ref=next((v for v in tab.values() if v[0]%r), None)
        rr=[]
        if ref:
            c,v=ref
            for a,u in tab.values():
                row=[(c*x-a*y)%r for x,y in zip(u,v)]
                assert sum(x*y for x,y in zip(row,vals))%r==0
                assert all(sum(row[j:j+k])%r==0 for j,k in zip(starts,sizes))
                rr.append(row)
        by_t.append(rr); rows.extend(rr)
    missing=max(sum(H[t][pin][0]%r==0 for t in disps) for pin in pins)
    return dict(rows=rows,by_t=by_t,missing=missing,n=n,W=W,d=d,H=H,vals=vals,disps=disps)


def curve13():
    p=13; a=0; c=2
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        s=(3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p)
        X=(s*s-x-u)%p;return X,(s*(x-X)-y)%p
    pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-c)%p==0]
    P=pts[0]; G=[None]
    while True:
        Q=add(G[-1],P)
        if Q is None:break
        G.append(Q)
    assert len(G)==19
    return G,add


def mean_checks():
    rng=Random(3102); fixtures=0; dirs=0; minbad=1.; elliptic=0
    for r,m,k,attempts in [(19,2,5,20),(31,2,6,20),(43,3,4,20),(59,3,5,10)]:
        for it in range(attempts):
            B=[tuple(sorted(rng.sample(range(1,r),k))) for _ in range(m)]
            cy,add=curve13() if r==19 else (None,None)
            f=short_fixture(r,B,cy,add)
            if 6*m*f['missing']>=r or 2*len(f['d'])<=r:continue
            assert rank(f['rows'],r,f['n'])==f['n']-m-1
            fixtures+=1;elliptic+=r==19
            for _ in range(8):
                u=[rng.randrange(r) for _ in range(f['n'])]
                bad=sum(any(sum(a*b for a,b in zip(row,u))%r for row in rr) for rr in f['by_t'])
                if not bad:continue
                assert bad*12*m>=r # deliberately weaker than theorem at delta<=1/(6m)
                dirs+=1;minbad=min(minbad,bad/r)
    assert fixtures and elliptic and dirs
    print(f'Conditional means: {fixtures} dense-support fixtures ({elliptic} elliptic); {dirs} nonaffine directions; minimum bad-displacement fraction {minbad:.6f}.')


def robust_checks():
    # Full nonidentity blocks on C19: energy tiny, not the LZR large-input masks.
    r=19; B=[list(range(1,r))]*2; n=36; W=18**2
    T=list(product(range(18),range(18))); vals=[(x+y+2)%r for x,y in T]
    d=Counter(vals); chi=r*sum(x*x for x in d.values())/W**2-1
    assert Fraction(r*sum(x*x for x in d.values())-W*W,W*W)<=Fraction(1,4096)
    rng=Random(321); minrate=1.
    for _ in range(50):
        y=[rng.randrange(r) for _ in range(n)]
        if y==B[0]+B[1]:continue
        H=[Counter() for _ in range(r)]
        for (a,b),g in zip(T,vals): H[g][(y[a]+y[18+b]-g)%r]+=1
        def fail(s,e):
            den=sum(d[g]*d[(g+s)%r] for g in range(r))
            good=sum(v*H[(g+s)%r][(a+e)%r] for g in range(r) for a,v in H[g].items())
            return Fraction(den-good,den)
        best=max([fail(s,0) for s in range(r)]+[fail(g,(v-g)%r) for v,g in zip(y,B[0]+B[1])])
        assert best>=Fraction(1,8); minrate=min(minrate,float(best))
    print(f'Short-tuple gap: 50 false candidates, chi-square={chi:.8f}, minimum best violation fraction={minrate:.6f}.')


def moment_checks():
    rng=Random(433); checks=0
    for r in (5,7,11,19,31,101):
        for J in range(1,min(5,r-1)+1):
            for _ in range(30):
                R=rng.sample(range(1,r),J); vals=[rng.choice([0]+R)+r*rng.randrange(-3,4) for _ in range(rng.randrange(1,min(r*r,100)))];
                if not any(x%r for x in vals):continue
                ks=[k for k in range(2,2*(J+1)) if sum(x**k for x in vals)%(r*r)]
                assert ks
                k=ks[0]; curr=vals
                while len(curr)>1:
                    h=len(curr)//2
                    curr=curr[:h] if sum(x**k for x in curr[:h])%(r*r) else curr[h:]
                    assert sum(x**k for x in curr)%(r*r)
                assert curr[0]%r
                checks+=1
    assert (3**2+4**2)%25==0 and (3**4+4**4)%25!=0
    r=101; h=10; R=[a for a in range(1,r) if pow(a,h,r)==1]
    assert len(R)==h and h*r<r*r
    assert all(r*sum(a**j for a in R)%(r*r)==0 for j in range(1,h))
    print(f'Moment support: {checks} exact randomized fixtures with deterministic extraction; quadratic cancellation control and multiplicative-subgroup obstruction passed.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    mean_checks();robust_checks();moment_checks()
