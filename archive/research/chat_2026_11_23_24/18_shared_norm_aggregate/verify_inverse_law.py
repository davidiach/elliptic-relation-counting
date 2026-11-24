#!/usr/bin/env python3
"""Exact finite checks of inverse planting and random-target reuse accounting."""
from itertools import product, combinations
from collections import Counter
from fractions import Fraction
from verify import make_curve, good_offsets
from verify_norm_aggregate import relation_poly, eval2

def check(p,a,c,m,w):
    r,h,lifts,logs=make_curve(p,a,c)
    pool=next(ss for ss in combinations(range(p),m*w) if good_offsets(p,a,c,ss))
    S=[pool[i*w:(i+1)*w] for i in range(m)]
    addr=list(product(*S)); W=len(addr)
    D=h*h*m*2**(m-1)
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        s=((v-y)*pow(u-x,-1,p) if x!=u else (3*x*x+a)*pow(2*y,-1,p))%p
        z=(s*s-x-u)%p;return z,(s*(x-z)-y)%p
    def mul(k,P):
        Q=None
        while k:
            if k&1:Q=add(Q,P)
            P=add(P,P);k>>=1
        return Q
    # logs are coordinates of [h]P. Pick P with log=1, then [h]P
    # is subgroup generator in this coordinate.
    G=mul(h,next(P for P,z in logs.items() if z==1))
    targets={j:mul(j,G) for j in range(1,r)}
    assert all(T is not None for T in targets.values())
    xT={j:T[0] for j,T in targets.items()}
    # Literal successful quadruples (address,translation,target,signed points).
    literal=Counter(); bycall={}; polys={s:relation_poly(p,a,c,s,h)[0] for s in addr}
    def tuples_at(s,A):return product(*(lifts[(A+x)%p] for x in s))
    for s in addr:
        for A in range(p):
            for pts in tuples_at(s,A):
                T=sum(logs[P] for P in pts)%r
                if T:literal[(s,A,T,pts)]+=1
    inverse=Counter();false_inputs=0;tested=0
    for s in addr:
        F=polys[s]
        for j in range(1,r):
            rec=[]
            for A in range(p):
                cand=eval2(F,A,xT[j],p)==0
                correct=[pts for pts in tuples_at(s,A) if sum(logs[P] for P in pts)%r==j]
                assert not correct or cand
                if cand and not correct:false_inputs+=1
                if cand:
                    for pts in correct:rec.append((s,A,j,pts))
                tested+=1
            assert len(rec)<=D
            bycall[(s,j)]=rec
            # each record receives exactly one acceptance ticket out of D
            inverse.update(rec)
    assert literal==inverse and all(v==1 for v in inverse.values())
    accepted=sum(inverse.values()); acceptance=Fraction(accepted,W*(r-1)*D)
    assert acceptance<=1
    # The distribution under the exact inverse rejection law equals H-weighting.
    H=Counter((A,j) for (s,A,j,pts) in inverse)
    H_direct=Counter()
    for A in range(p):
        slots=[[P for x in Si for P in lifts[(A+x)%p]] for Si in S]
        for pts in product(*slots):
            j=sum(logs[P] for P in pts)%r
            if j:H_direct[A,j]+=1
    assert H==H_direct
    # Factor base projected sign classes. Only nonidentity variables counted.
    base=[]
    for A in range(p):
        base.append({min(logs[P],r-logs[P]) for x in pool for P in lifts[(A+x)%p] if logs[P]})
    shifts={(x-y)%p for x in pool for y in pool}
    B=len(shifts) if h==1 else 2*h*len(pool)**2
    overlap_checks=0
    for A0 in range(p):
        U={A for A in range(p) if base[A]&base[A0]}
        assert len(U)<=B
        if h==1:assert all((A-A0)%p in shifts for A in U)
        for s in addr:
            # Exactly the uniform targets whose inverse lists overlap the old base.
            hit=sum(any(row[1] in U for row in bycall[(s,j)]) for j in range(1,r))
            assert hit<=len(U)*2**m
            overlap_checks+=1
    print(f'p={p}, r={r}, h={h}, m={m}, w={w}, offsets={pool}: {tested} root/target/translation checks, {accepted} accepted inverse tickets, exact acceptance={acceptance}, {false_inputs} algebraic extra candidate incidences, {overlap_checks} overlap checks',flush=True)
    return tested,accepted,false_inputs,overlap_checks

if __name__=='__main__':
    totals=[0]*4
    for args in [(13,4,1,2,2),(17,2,2,3,2),(29,1,0,2,2),(43,2,1,2,2)]:
        z=check(*args);totals=[a+b for a,b in zip(totals,z)]
    print('totals (incidences, accepted tickets, algebraic extras, overlap checks):',totals)
    print('PASS: exact planted inverse distribution and reuse bound checks.')
