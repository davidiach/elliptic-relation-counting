#!/usr/bin/env python3
"""Verify one-call recovery from the full truncated norm, not just its formula."""
from itertools import product
from random import Random
from math import prod
from verify import semaev
rng=Random(77261124)

def ringmul(A,B,p,K):
    out={}
    for (z,y,e),a in A.items():
        for (u,v,f),b in B.items():
            if z+u>K or y+v>K or (e and f):continue
            k=(z+u,y+v,e+f);out[k]=(out.get(k,0)+a*b)%p
    return {k:v for k,v in out.items() if v}

def ev(co,x,p):
    r=0
    for c in co[::-1]:r=(r*x+c)%p
    return r

def check(p,grid,F,K,lin):
    m=len(grid[0]);S=[s for s in grid if F[s]%p==0];H=len(S)
    P={(0,0,0):1}
    for s in grid:
        B={(0,0,0):F[s]%p}
        if K:
            B.update({(1,1,0):1,(1,0,0):-sum(a*x for a,x in zip(lin,s))%p})
            B.update({(1,0,j+1):-x%p for j,x in enumerate(s)})
        P=ringmul(P,B,p,K)
    if H>K:
        assert not P
        return 'overflow',H
    z=min(i for (i,j,k),c in P.items() if c and k==0)
    assert z==H
    C=prod(v for v in F.values() if v%p)%p
    assert P[H,H,0]==C
    if H==0:return 'empty',H
    poly=[P.get((H,i,0),0) for i in range(H+1)]
    expected=[C]
    for s in S:
        y=sum(a*x for a,x in zip(lin,s))%p;co=[0]*(len(expected)+1)
        for i,c in enumerate(expected):co[i]=(co[i]-y*c)%p;co[i+1]=(co[i+1]+c)%p
        expected=co
    assert poly==expected
    roots=[y for y in range(p) if ev(poly,y,p)==0]
    if len(roots)<H:
        assert len({sum(a*x for a,x in zip(lin,s))%p for s in S})<H
        return 'collision',H
    der=[i*poly[i]%p for i in range(1,H+1)];found=[]
    for y in roots:
        slope=ev(der,y,p);assert slope
        s=tuple(-ev([P.get((H,j,i+1),0) for j in range(H+1)],y,p)*pow(slope,-1,p)%p for i in range(m))
        assert s in F and F[s]%p==0
        found.append(s)
    assert set(found)==set(S)
    return 'recovered',H

if __name__=='__main__':
    assert __debug__
    counts={};roots=actual=0
    for trial in range(240):
        p=rng.choice([31,101,257]);m=rng.choice([2,3]);S=[rng.sample(range(p),rng.choice([2,3])) for j in range(m)];grid=list(product(*S));K=rng.choice([0,1,2,3,5]);H=rng.randrange(min(7,len(grid))+1);zero=set(rng.sample(grid,H));F={s:0 if s in zero else rng.randrange(1,p) for s in grid}
        lin=[0]*m if trial%9==0 else [rng.randrange(p) for j in range(m)]
        state,h=check(p,grid,F,K,lin);counts[state]=counts.get(state,0)+1;roots+=h if state=='recovered' else 0
    for p,a,c in [(19,2,1),(101,2,3)]:
        F=semaev(p,a,c);targets=[x for x in range(p) if (x**3+a*x+c)%p and pow((x**3+a*x+c)%p,(p-1)//2,p)==1]
        for trial in range(50):
            m=rng.choice([2,3]); pool=rng.sample(range(p),m*2);S=[pool[2*j:2*j+2] for j in range(m)];grid=list(product(*S));T=rng.choice(targets);vals={s:F(s,T) for s in grid};K=3;lin=[rng.randrange(p) for j in range(m)]
            state,h=check(p,grid,vals,K,lin);counts[state]=counts.get(state,0)+1;roots+=h if state=='recovered' else 0;actual+=1
    print(f'{sum(counts.values())} full truncated norms, including {actual} actual summation-polynomial cases: {counts}; {roots} recovered zero addresses.')
    print('PASS: valuation, projection polynomial, collision handling and all reconstructed coordinates.')
