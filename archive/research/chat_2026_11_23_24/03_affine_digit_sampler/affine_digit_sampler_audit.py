"""Finite tests for the affine digit-block sampler and algebraic readout obstruction.
Uses small finite fields only; not a discrete-log solver or a contraction benchmark.
"""
from collections import Counter
from itertools import product, combinations
from math import comb, ceil, isqrt
from random import Random
import numpy as np


def points(p,a,c):
    return [None]+[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-a*x-c)%p==0]


def add(P,Q,p,a):
    if P is None: return Q
    if Q is None: return P
    x,y=P; u,v=Q
    if x==u and (y+v)%p==0: return None
    l=((v-y)*pow(u-x,-1,p) if x!=u else (3*x*x+a)*pow(2*y,-1,p))%p
    w=(l*l-x-u)%p
    return (w,(l*(x-w)-y)%p)


def neg(P,p): return None if P is None else (P[0],-P[1]%p)


def target_tuples(Ps,T,m,p,a):
    for pp in product(Ps[1:],repeat=m-1):
        s=None
        for P in pp: s=add(s,P,p,a)
        last=add(T,neg(s,p),p,a)
        if last is not None: yield pp+(last,)


def affine_pair_audit(p,S):
    masks=[{(a+t*s)%p for s in S} for a in range(p) for t in range(1,p)]
    singles=Counter(x for mask in masks for x in mask)
    pairs=Counter((x,y) for mask in masks for x in mask for y in mask if x!=y)
    assert all(singles[x]==len(S)*(p-1) for x in range(p))
    assert all(pairs[x,y]==len(S)*(len(S)-1) for x in range(p) for y in range(p) if x!=y)
    return len(masks),p*(p-1)


def digit_block_audit():
    # Toy surjective coefficient box: ordinary binary evaluation modulo p.
    # Split owner rows and within-word bit slices; the high-bit offset is absorbed in t.
    cases=[(13,4,2,2),(101,6,2,2),(257,8,2,3),(1031,8,2,3)]
    total=0
    for p,n,k,r in cases:
        V=4; A=n*k; radix=1<<k
        owner_groups=[list(range(g,V,r)) for g in range(r)]
        cuts=[(j*k//r,(j+1)*k//r) for j in range(r)]
        # Some slices may be empty at very small k; this is harmless.
        raw_blocks=[]; normalized=[]
        for og in owner_groups:
            inds=[i for i in range(n) if i%V in og]
            for lo,hi in cuts:
                s={sum(v*pow(radix,i,p) for v,i in zip(ds,inds))%p
                   for ds in product(range(1<<(hi-lo)),repeat=len(inds))}
                normalized.append(s)
                raw_blocks.append({x*(1<<lo)%p for x in s})
        ss={0}
        for s in raw_blocks: ss={(u+v)%p for u in ss for v in s}
        assert len(ss)==p
        prodsize=1
        for s in normalized: prodsize*=len(s)
        assert prodsize>=p and max(map(len,normalized))**(r*r)>=p
        total+=1
    return total


def second_moment_audit(p,a,c,S):
    Ps=points(p,a,c); N=len(Ps); w=len(S); seeds=p*(p-1)
    masks=[{(u+t*s)%p for s in S} for u in range(p) for t in range(1,p)]
    Fs=[[P for P in Ps[1:] if P[0] in mask] for mask in masks]
    eligible=[P for P in Ps[1:] if add(P,P,p,a) is not None]
    # Sum over all seed pairs, all targets, exact integer numerators.
    totals=Counter(); seconds=Counter(); positives=Counter()
    for A,B in product(Fs,repeat=2):
        h=Counter(add(P,Q,p,a) for P in A for Q in B)
        for T,z in h.items():
            totals[T]+=z; seconds[T]+=z*(z-1); positives[T]+=1
    # Reconstruct moment numerators independently from actual rational tuples
    # and the exact affine pair-retention numerators.
    for T in eligible:
        rows=list(target_tuples(Ps,T,2,p,a))
        assert totals[T]==len(rows)*(w*(p-1))**2
        expect=0
        for U,V in product(rows,repeat=2):
            if U==V: continue
            q=1
            for X,Y in zip(U,V):
                q*=w*(p-1) if X[0]==Y[0] else w*(w-1)
            expect+=q
        assert seconds[T]==expect
        # Cauchy-Schwarz exact rational lower bound.
        assert positives[T]*(seconds[T]+totals[T])>=totals[T]**2
    return len(eligible),seeds**2


def compositions_leq(d,k):
    if k==0:
        yield (); return
    for j in range(d+1):
        for tail in compositions_leq(d-j,k-1): yield (j,)+tail


def planted_seed_rank(p,a,c,m,S,D,seed):
    Ps=points(p,a,c)
    T=next(P for P in Ps[1:] if add(P,P,p,a) is not None)
    rows=list(target_tuples(Ps,T,m,p,a))
    exps=np.array(list(compositions_leq(D,2*m)),dtype=np.int64)
    N=len(exps); assert N==comb(D+2*m,2*m)
    rng=Random(seed); basis=[]; pivots=[]; used=set(); attempts=0
    while len(basis)<N and attempts<20*N:
        attempts+=1
        pp=rng.choice(rows); ss=[rng.choice(S) for _ in range(m)]; tt=[rng.randrange(1,p) for _ in range(m)]
        aa=[(P[0]-s*t)%p for P,s,t in zip(pp,ss,tt)]
        key=tuple(aa+tt)
        if key in used: continue
        used.add(key)
        powers=np.array([[pow(x,j,p) for j in range(D+1)] for x in key],dtype=np.int64)
        v=np.ones(N,dtype=np.int64)
        for j in range(2*m): v=v*powers[j,exps[:,j]]%p
        for u,piv in zip(basis,pivots):
            if v[piv]: v=(v-v[piv]*u)%p
        nz=np.flatnonzero(v)
        if len(nz):
            piv=int(nz[0]); v=v*pow(int(v[piv]),-1,p)%p
            basis.append(v); pivots.append(piv)
    assert len(basis)==N,(p,m,D,len(basis),N)
    return N,len(used),len(rows)


if __name__=='__main__':
    for p,S in [(13,{0,1,3,7}),(17,{0,1,4,8,9}),(37,{0,2,3,10,12,15})]:
        n,j=affine_pair_audit(p,S)
        print(f'p={p}: {n} affine masks; all {j} ordered distinct-abscissa inclusion counts exact')
    print(f'{digit_block_audit()} digit-block image/entropy checks passed')
    for p,a,c,S in [(13,4,1,{0,1,3,7}),(17,2,2,{0,1,4,8,9})]:
        t,n=second_moment_audit(p,a,c,S)
        print(f'p={p}: {n} ordered mask pairs; exact first and second moments at all {t} eligible targets')
    for p,a,c,m,S,D in [(13,4,1,2,[0,1,3,7],8),(17,2,2,2,[0,1,4,8,9],8),(13,4,1,3,[0,1,3,7],4)]:
        rank,seeds,rows=planted_seed_rank(p,a,c,m,S,D,20261124+m+p)
        print(f'p={p}, m={m}, degree={D}: planted evaluation rank {rank}/{rank}; {seeds} distinct sampled positive seeds; {rows} rational target tuples')
    print('All finite checks passed. No end-to-end contraction benchmark.')
