"""Exact toy checks for the two-coordinate dictionary and simultaneous criterion obstruction."""
from itertools import product
from random import Random
from math import prod
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def rank_mod(rows,p,n):
    piv={}
    for row in rows:
        a=np.asarray(row,dtype=np.int64).copy()%p
        for j in sorted(piv):
            if a[j]:a=(a-int(a[j])*piv[j])%p
        ix=np.flatnonzero(a)
        if len(ix):
            j=int(ix[0]);piv[j]=a*pow(int(a[j]),-1,p)%p
    return len(piv)


def digits(a,D,m):
    z=[]
    for _ in range(m):a,j=divmod(a,D);z.append(j)
    assert not a
    return z


def dictionary(r,m,extra,seed):
    D=1<<((r.bit_length()+m-1)//m)
    if 2*(D+1)+1+extra>=r:raise ValueError('too small')
    assert D**m>r
    forbidden={j*D**i%r for i in range(m) for j in range(-D,D+1)}
    zeta=next(z for z in range(1,r) if z not in forbidden)
    T=[[(0,j*D**i) for j in range(D+1)]+[(j*D**i,j*D**i+zeta) for j in range(D+1)] for i in range(m)]
    used={t%r for B in T for _,t in B}
    tau=next(t for t in range(1,r) if pow(t,-1,r) not in used and -pow(t,-1,r)%r not in used)
    B=[]
    for pairs in T:
        pairs=pairs+[(0,pow(tau,-1,r))]
        pts=[((1+tau*t)%r,c%r,t%r) for c,t in pairs]
        assert 0 not in [a for a,c,t in pts] and len({a for a,c,t in pts})==len(pts)
        B.append(pts)
    rng=Random(seed)
    for pts in B:
        used={a for a,c,t in pts}
        for a in rng.sample([x for x in range(1,r) if x not in used],extra):
            pts.append((a,rng.randrange(r),(a-1)*pow(tau,-1,r)%r))
    return D,zeta,tau,B


def selected_rows(r,m,D,zeta,tau,B):
    n=sum(len(F)-1 for F in B);starts=[];s=0
    for F in B:starts.append(s);s+=len(F)-1
    u=[c for F in B for a,c,t in F[1:]];lam=[(a-1)%r for F in B for a,c,t in F[1:]]
    zero=(0,)*(2*m);rows=[];off=D+1;special=2*(D+1)
    def change(X,updates):
        Y=list(X)
        for i,j,k in updates:Y[2*i+j]=k
        return tuple(Y)
    def inc(X):
        a=[0]*n
        for i,j in enumerate(X):
            if j:a[starts[i//2]+j-1]+=1
        return a
    def vals(X):
        return tuple(sum(B[i//2][j][k] for i,j in enumerate(X))%r for k in range(3))
    def add(X,Y):
        x,c,t=vals(X);y,d,v=vals(Y)
        assert (y-x)%r==1 and c==d and (v-t)%r==pow(tau,-1,r)
        row=[b-a for a,b in zip(inc(X),inc(Y))]
        assert sum(a*b for a,b in zip(row,lam))%r==1
        assert sum(a*b for a,b in zip(row,u))%r==0
        rows.append(row)
    def contrast(X,Xp):
        assert vals(X)==vals(Xp)
        i=next(i for i,v in enumerate(X) if v==0)
        Y=change(X,[(i//2,i%2,special)])
        add(X,Y);add(Xp,Y)
    def dtup(t):return change(zero,[(i,0,j) for i,j in enumerate(digits(t,D,m))])
    def ctup(c):return change(zero,[(i,0,off+j) for i,j in enumerate(digits(c,D,m))])
    def ctup_with_t(c,t):
        return change(ctup(c),[(i,1,j) for i,j in enumerate(digits((t-c-m*zeta)%r,D,m))])
    for i,F in enumerate(B):
        add(zero,change(zero,[(i,0,special)]))
        # d and c digits, including c offset at j=0.
        contrast(change(zero,[(i,0,off)]),dtup(zeta))
        contrast(change(zero,[(i,0,special)]),dtup(pow(tau,-1,r)))
        for j in range(2,D+1):
            contrast(change(zero,[(i,0,j)]),change(zero,[(i,0,j-1),(i,1,1)]))
            contrast(change(zero,[(i,0,off+j),(i,1,off)]),change(zero,[(i,0,off+j-1),(i,1,off+1)]))
        if i:
            contrast(change(zero,[(i-1,0,D)]),change(zero,[(i,0,1)]))
            contrast(change(zero,[(i-1,0,off+D)]),change(zero,[(i,0,off+1)]))
        for j,(a,c,t) in enumerate(F[special+1:],special+1):
            contrast(change(zero,[(i,0,j)]),ctup_with_t(c,t))
    contrast(zero,dtup(r));contrast(ctup(0),ctup(r))
    assert rank_mod(rows,r,n)==n-1
    H=smith_normal_form(sp.Matrix(rows),domain=sp.ZZ)
    assert [abs(int(H[i,i])) for i in range(n)]==[1]*(n-1)+[r]
    return rows,u,lam,ctup_with_t,vals


def min_range(r,blocks,tau):
    """Exhaust ALL nonzero (h,k) in the augmented solution plane.
    A cyclic range is the shortest integer-lift interval, allowing wrap cuts.
    With the zero reference present, its endpoint can be shifted to include lift 0.
    """
    minR=r*len(blocks);arg=None;perblock=[r]*len(blocks)
    k=np.arange(r,dtype=np.int64)[:,None]
    arrays=[(np.array([c for a,c,t in F],dtype=np.int64)[None,:],np.array([(a-1)%r for a,c,t in F],dtype=np.int64)[None,:]) for F in blocks]
    for h in range(r):
        scores=np.zeros(r,dtype=np.int64)
        for i,(c,l) in enumerate(arrays):
            x=np.sort((h*l+k*c)%r,axis=1)
            maxgap=np.maximum(np.diff(x,axis=1).max(axis=1),r+x[:,0]-x[:,-1])
            w=r-maxgap;scores+=w
            perblock[i]=min(perblock[i],int(w[1 if h==0 else 0:].min()))
        if h==0:scores[0]=r*len(blocks)
        j=int(scores.argmin());v=int(scores[j])
        if v<minR:minR=v;arg=(h,j)
    return minR,arg,perblock


def count_pins(r,B):
    F=[[a for a,c,t in block] for block in B]
    def conv(slots):
        x=np.zeros(r,dtype=np.int64);x[0]=1
        for S in slots:
            y=np.zeros_like(x)
            for a in S:y+=np.roll(x,a)
            x=y
        return x
    d=conv(F);assert int(d.sum())==prod(map(len,F)) and d.min()>0
    pins=[conv([[-a%r for a in S] for S in F]+[S for j,S in enumerate(F) if j!=i]) for i in range(len(F))]
    return (int(d.min()),int(d.max())),(min(int(x.min()) for x in pins),max(int(x.max()) for x in pins))


def run():
    from collections import Counter
    from fractions import Fraction
    for r,m,extra,seed in [(101,2,8,17),(257,3,12,22),(1009,2,10,17)]:
        D,zeta,tau,B=dictionary(r,m,extra,seed)
        A,u,lam,make,vals=selected_rows(r,m,D,zeta,tau,B)
        # Check the whole joint two-copy coverage by an independent direct sumset DP
        # at the small prime; at all primes test encoded examples.
        for c,t in product(range(r),repeat=2) if r==101 else [(x,(137*x+13)%r) for x in range(r)]:
            X=make(c,t);g,cc,tt=vals(X);assert (cc,tt)==(c,t)
        if r==101:
            state={(0,0)}
            for F in B:
                for _ in range(2):
                    state={((a+c)%r,(b+t)%r) for a,b in state for g,c,t in F}
            assert len(state)==r*r
        R,arg,perblock=min_range(r,B,tau)
        assert 2*R>=r
        fractions=[Fraction(max(Counter(c for a,c,t in F).values()),len(F)) for F in B]
        fibers,pins=count_pins(r,B)
        assert R>=1
        print(f'r={r}, m={m}, block sizes={[len(F) for F in B]}, D={D}, zeta={zeta}, tau={tau}: {len(A)} genuine rows, rational rank {len(u)}, modular rank {len(u)-1}, selected index {r}.')
        print(f'  Min summed circular width over {r*r-1} nonzero (h,k): {R} at {arg}; individual-block minima {perblock}; 2R={2*R} >= r. Largest quotient class fractions {[str(f) for f in fractions]}.')
        print(f'  Short tuple counts {fibers}; one-pin masses {pins}; no-wrap {0<pins[0]<=pins[1]<r}. Joint doubled-tuple coverage checked.')
        if r==1009:assert 0<pins[0]<=pins[1]<r

if __name__=='__main__':run()
