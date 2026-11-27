"""Exact finite tests of a structural counterexample, not an asymptotic algorithm."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product, combinations
from math import prod
from random import Random
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def dictionary(r,m,extra=0,seed=17):
    base=1<<((r.bit_length()+m-1)//m)
    D=base
    assert D<r and D**m>r
    for kap in range(1,r):
        pairs=[[(j*D**i,1) for j in range(D+1)]+[(0,-1),(0,2),(0,-2)] for i in range(m)]
        A=[[(kap*c+d)%r for c,d in B] for B in pairs]
        if all(0 not in B and len(B)==len(set(B)) for B in A):break
    else: raise ValueError('no dictionary')
    rng=Random(seed)
    for i,B in enumerate(A):
        for a in rng.sample([a for a in range(1,r) if a not in set(B)],extra):
            pairs[i].append(((a-2)*pow(kap,-1,r)%r,2))
    return D,kap,[[((kap*c+d)%r,c,d) for c,d in B] for B in pairs]


def rank_mod(rows,p,n):
    piv={}
    for rr in rows:
        a=np.array(rr,dtype=np.int64)%p
        for j in sorted(piv):
            if a[j]:a=(a-int(a[j])*piv[j])%p
        ix=np.flatnonzero(a)
        if len(ix):
            j=int(ix[0]);piv[j]=a*pow(int(a[j]),-1,p)%p
    return len(piv)


def digits(c,D,m):
    out=[]
    for i in range(m):c,j=divmod(c,D);out.append(j)
    assert c==0
    return out


def rows_for(r,D,kap,B,elliptic=None):
    m=len(B);n=sum(len(F)-1 for F in B); offsets=[]; s=0
    for F in B:offsets.append(s);s+=len(F)-1
    u=[c%r for F in B for _,c,d in F[1:]]
    lam=[(a-1)%r for F in B for a,c,d in F[1:]]
    zero=(0,)*(2*m);rows=[]
    def change(X,updates):
        Y=list(X)
        for i,slot,j in updates:Y[2*i+slot]=j
        return tuple(Y)
    def inc(X):
        a=[0]*n
        for i,j in enumerate(X):
            if j:a[offsets[i//2]+j-1]+=1
        return a
    def value(X):return sum(B[i//2][j][0] for i,j in enumerate(X))%r
    def add(X,Y):
        row=[y-x for x,y in zip(inc(X),inc(Y))]
        assert (value(Y)-value(X))%r==1
        assert sum(a*c for a,c in zip(row,u))%r==0
        assert sum(a*c for a,c in zip(row,lam))%r==1
        if elliptic:
            G,law=elliptic
            def eval_t(T):
                R=None
                for i,j in enumerate(T):R=law(R,G[B[i//2][j][0]])
                return R
            assert eval_t(Y)==law(eval_t(X),G[1])
        rows.append(row)
    def contrast(X,Xp):
        # Both are padded doubled tuples; a reference in X provides a successor.
        i=next(j for j,v in enumerate(X) if v==0)
        Y=change(X,[(i//2,i%2,D+2)])
        assert value(X)==value(Xp)
        assert sum(a*b for a,b in zip(inc(X),u))%r==sum(a*b for a,b in zip(inc(Xp),u))%r
        add(X,Y);add(Xp,Y)
    for i,F in enumerate(B):
        neg,tw,neg2=D+1,D+2,D+3
        add(zero,change(zero,[(i,0,tw)]))
        add(change(zero,[(i,0,neg2)]),change(zero,[(i,0,neg)]))
        add(change(zero,[(i,0,neg),(i,1,tw)]),zero)
        for j in range(2,D+1):
            contrast(change(zero,[(i,0,j)]),change(zero,[(i,0,j-1),(i,1,1)]))
        if i:
            contrast(change(zero,[(i-1,0,D)]),change(zero,[(i,0,1)]))
        for j,(a,c,d) in enumerate(F[D+4:],D+4):
            assert d==2
            X=change(zero,[(i,0,j)])
            Xp=change(zero,[(k,0,v) for k,v in enumerate(digits(c,D,m))]+[(0,1,tw)])
            contrast(X,Xp)
    contrast(zero,change(zero,[(i,0,j) for i,j in enumerate(digits(r,D,m))]))
    assert rank_mod(rows,r,n)==n-1
    return rows,u,lam


def elliptic13():
    p=7
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;z,w=Q
        if x==z and (y+w)%p==0:return None
        a=((3*x*x)*pow(2*y,-1,p) if P==Q else (w-y)*pow(z-x,-1,p))%p
        X=(a*a-x-z)%p;return X,(a*(x-X)-y)%p
    pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-3)%p==0]
    P=pts[0];G=[None]
    for _ in range(12):G.append(add(G[-1],P))
    assert len(pts)==12 and set(G)=={None,*pts} and add(G[-1],P) is None
    return G,add


def shift_support(r,B,copies):
    states={(0,0)}
    for F in B:
        for _ in range(copies):
            states={((a+c)%r,b+d) for a,b in states for _,c,d in F}
    by=defaultdict(set)
    for c,d in states:by[c].add(d)
    return {(d-e)%r for U in by.values() for d in U for e in U}


def convolve_slots(r,slots):
    x=np.zeros(r,dtype=np.int64);x[0]=1
    for F in slots:
        y=np.zeros_like(x)
        for a in F:y+=np.roll(x,a)
        x=y
    return x


def pointwise_counts(r,B):
    F=[[a for a,c,d in block] for block in B]
    d=convolve_slots(r,F);W=prod(map(len,F));assert int(d.sum())==W and min(d)>0
    pin=[]
    for i in range(len(F)):
        # A target-side pin only translates the convolution of the complete
        # source half reflected, with every remaining target slot.
        q=convolve_slots(r,[[-a%r for a in A] for A in F]+[A for j,A in enumerate(F) if j!=i]);assert int(q.sum())==W*W//len(F[i]);pin.append((int(q.min()),int(q.max())))
    return (int(d.min()),int(d.max())),(min(a for a,b in pin),max(b for a,b in pin)),W


def lattice_tests():
    for r,m,extra in [(13,4,1),(23,3,2),(31,2,3),(101,3,3),(257,3,3),(1009,2,4),(4093,3,2)]:
        D,kap,B=dictionary(r,m,extra)
        A,u,lam=rows_for(r,D,kap,B,elliptic13() if r==13 else None)
        n=len(u);s=smith_normal_form(sp.Matrix(A),domain=sp.ZZ)
        assert [abs(int(s[i,i])) for i in range(n)]==[1]*(n-1)+[r]
        targets=shift_support(r,B,2);bound=min(r,16*m+1)
        assert len(targets)<=bound and all(min(q,r-q)<=8*m for q in targets)
        fc,pc,W=pointwise_counts(r,B)
        print(f'r={r}, m={m}, D={D}, block sizes {[len(F) for F in B]}: {len(A)} genuine rows, rational rank {n}, modular rank {n-1}, selected lattice index {r}; individual targets {len(targets)}/{r} (bound {bound}); short fibers {fc}; pin counts {pc}.')
        if r in (1009,4093):assert 0<pc[0]<=pc[1]<r
        if r==13: print('  Every selected row also checked with actual group operations on E/F7 of order 13.')


def partitions(v):
    D={};return tuple(D.setdefault(x,len(D)) for x in v)


def bell(h):
    vals=[1]
    from math import comb
    for n in range(1,h+1):vals.append(sum(comb(n-1,k)*vals[k] for k in range(n)))
    return vals[-1]


def expected_centered(pats,h,k,theta):
    ans=Q(0)
    for mask in range(1<<h):
        U=[j for j in range(h) if mask>>j&1]
        ell=sum(len({p[j] for j in U}) for p in pats)
        ans+=(-1)**(h-len(U))*theta**(ell+k*(h-len(U)))
    return ans


def weighted_moments():
    r=5;k=2;ambient=[(1,2,3),(1,3,4)];T=list(product(*ambient));w=[2,1,0,3,0];theta=Q(1,2);S=sum(x>0 for x in w);M=max(w);q=r*theta;mu=q**k*S/r
    seeds=list(product((0,1),repeat=6));expectations=0;patterns=0
    for g in range(r):
        valid=[(X,w[(g-sum(X))%r]) for X in T if w[(g-sum(X))%r]]
        counts=[]
        for seed in seeds:
            counts.append(sum(v*prod(seed[3*i+ambient[i].index(x)] for i,x in enumerate(X)) for X,v in valid))
        E=Q(sum(counts),len(counts));assert E==sum(v for X,v in valid)*theta**k
        for h in (2,4):
            direct=Q(sum((Q(c)-E)**h for c in counts),len(counts))
            raw=Q(sum(c**h for c in counts),len(counts));tab=Counter();num=Counter()
            for z in product(valid,repeat=h):
                X=[t[0] for t in z];pat=tuple(partitions(v) for v in zip(*X));tab[pat]+=prod(t[1] for t in z);num[pat]+=1
            via=sum(v*expected_centered(pat,h,k,theta) for pat,v in tab.items())
            rawvia=sum(v*theta**sum(1+max(p) for p in pat) for pat,v in tab.items())
            assert via==direct and rawvia==raw
            for pat,narr in num.items():
                cols=[[int(p[j]==a) for j in range(h)] for p in pat for a in range(1+max(p))];R=rank_mod(list(zip(*cols)),r,len(cols));ell=len(cols)
                assert narr<=r**(ell-R)*S**R
                assert all(1+max(p)<=R for p in pat)
                isolated=any(all(p.count(p[j])==1 for p in pat) for j in range(h))
                if isolated:assert expected_centered(pat,h,k,theta)==0
                else:assert 4*(h-R)>=h or 4*(k*h-ell-k*(h-R))>=h
                patterns+=1
            assert raw<=M**h*bell(h)**k*max(Q(1),mu)**h
            # K*=1 is sufficient here. The asymptotic proof uses large K*.
            assert direct<=2**h*M**h*bell(h)**k*mu**h
            expectations+=1
    print(f'Weighted moment checks: {expectations} exact centered/raw identities on {len(seeds)} seeds each; {patterns} partition rank/count checks.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    lattice_tests();weighted_moments()
