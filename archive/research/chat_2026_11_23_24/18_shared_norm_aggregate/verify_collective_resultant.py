#!/usr/bin/env python3
"""Exact fixed-degree resultant checks and aggregate spectrum regression."""
from itertools import product
from collections import Counter
from random import Random
from math import prod
from verify import good_offsets, make_curve
from verify_norm_aggregate import relation_poly, ev, eval2


def det(M,p):
    M=[row[:] for row in M]; n=len(M); d=1
    for i in range(n):
        j=next((j for j in range(i,n) if M[j][i]%p),None)
        if j is None:return 0
        if j!=i:M[i],M[j]=M[j],M[i];d=-d
        x=M[i][i]%p;d=d*x%p;inv=pow(x,-1,p)
        for j in range(i+1,n):
            v=M[j][i]*inv%p
            if v:
                for k in range(i+1,n):M[j][k]=(M[j][k]-v*M[i][k])%p
                M[j][i]=0
    return d%p

def res(f,g,p):
    # f,g are ascending coefficient lists with declared degree (may have
    # leading zero). This computes the homogeneous resultant at that degree.
    d,e=len(f)-1,len(g)-1;n=d+e
    M=[[0]*n for _ in range(n)]
    for i in range(e):M[i][i:i+d+1]=f[::-1]
    for i in range(d):M[e+i][i:i+e+1]=g[::-1]
    return det(M,p)

def mul(f,g,p):
    a=[0]*(len(f)+len(g)-1)
    for i,x in enumerate(f):
        for j,y in enumerate(g):a[i+j]=(a[i+j]+x*y)%p
    return a

def s3(v,z,a,c,p):
    return [( (v*z-a)**2-4*c*(v+z))%p,(-2*((v+z)*(v*z+a)+2*c))%p,(v-z)**2%p]

def f_at(F,A,p):return [ev(c,A,p) for c in F]

def compare(A,B,p,ratios):
    if A==0:assert B==0;return
    assert B!=0;ratios.add(B*pow(A,-1,p)%p)

def interp(vals,p):
    # interpolate at 0,...,len(vals)-1 using polynomial Newton differences
    dif=vals[:];poly=[0]*len(vals);basis=[1];fac=1
    for j in range(len(vals)):
        x=dif[0]*pow(fac,-1,p)%p
        for k,v in enumerate(basis):poly[k]=(poly[k]+v*x)%p
        dif=[(b-a)%p for a,b in zip(dif,dif[1:])]
        basis=mul(basis,[-j,1],p);fac=fac*(j+1)%p
    return poly

def recurrence(p,a,c,ss):
    F=relation_poly(p,a,c,ss)[0];ratios=set();balanced=set();checks=bc=0
    if len(ss)==2: G=((-ss[0]%p,p-1),(1,)) # y-(A+s0)
    else:G=relation_poly(p,a,c,ss[:-1])[0]
    mid=len(ss)//2
    def subpoly(s):return ((-s[0]%p,p-1),(1,)) if len(s)==1 else relation_poly(p,a,c,s)[0]
    L,R=subpoly(ss[:mid]),subpoly(ss[mid:]); dR=len(R)-1;assert 2*dR<p
    for A in range(p):
        g=f_at(G,A,p);l=f_at(L,A,p);r=f_at(R,A,p)
        for x in range(p):
            val=eval2(F,A,x,p)
            other=res(g,s3((A+ss[-1])%p,x,a,c,p),p)
            compare(val,other,p,ratios);checks+=1
            # inner resultant in Z interpolated as polynomial in Y.
            inner=interp([res(r,s3(y,x,a,c,p),p) for y in range(2*dR+1)],p)
            for y in [0,p-1]:assert ev(inner,y,p)==res(r,s3(y,x,a,c,p),p)
            compare(val,res(l,inner,p),p,balanced);bc+=1
    assert len(ratios)==len(balanced)==1
    print(f'p={p}, m={len(ss)}: {checks} sequential and {bc} balanced homogeneous resultant comparisons; ratios {ratios}, {balanced}',flush=True)
    return checks,bc

def aggregate():
    p,a,c=19,2,1; S=[(0,1),(2,3),(4,5)];assert good_offsets(p,a,c,sum(S,()))
    polys={s:relation_poly(p,a,c,s)[0] for s in product(*S)}
    partial={s:relation_poly(p,a,c,s)[0] for s in product(*S[:2])}
    ratios=set();checks=0
    for A in range(p):
        f=[1]
        for F in partial.values():f=mul(f,f_at(F,A,p),p)
        for x in range(p):
            g=[1]
            for v in S[-1]:g=mul(g,s3((A+v)%p,x,a,c,p),p)
            truth=prod(eval2(F,A,x,p) for F in polys.values())%p
            compare(truth,res(f,g,p),p,ratios);checks+=1
    assert len(ratios)==1
    print(f'{checks} exact collective resultant/product comparisons over the full grid',flush=True)
    return checks

def spectrum():
    # Distinct factor polynomials and eigenvalue collisions at random pairs.
    p=257; a,c=next((a,c) for a in range(1,12) for c in range(1,15) if (4*a**3+27*c*c)%p and good_offsets(p,a,c,range(12)) and all((sum(1+(0 if (x**3+a*x+c)%p==0 else (1 if pow((x**3+a*x+c)%p,(p-1)//2,p)==1 else -1)) for x in range(p))+1)%d for d in range(2,18)))
    r,h,lifts,logs=make_curve(p,a,c); assert h==1
    S=[tuple(range(i*4,(i+1)*4)) for i in range(3)];polys=[relation_poly(p,a,c,s)[0] for s in product(*S)]
    W=len(polys);assert len(set(polys))==W
    rng=Random(932); pts=[P for P in logs if P is not None];K=512
    distinct=[];energies=[];pos=[]
    for _ in range(K):
        A=rng.randrange(p);x=rng.choice(pts)[0]
        counts=Counter(eval2(F,A,x,p) for F in polys)
        d=len(counts);e=sum(v*v for v in counts.values())
        assert d*e>=W*W
        distinct.append(d);energies.append(e)
        if counts[0]:pos.append(d)
    print(f'spectrum p={p}, curve a={a}, c={c}, r={r}, W={W}: {K} pair samples, {W} distinct symbolic factors; numeric distinct values min/mean/max={min(distinct)}/{sum(distinct)/K:.3f}/{max(distinct)}; mean energy={sum(energies)/K:.3f}; {len(pos)} algebraically positive samples, mean distinct on positives={sum(pos)/len(pos):.3f}',flush=True)

if __name__=='__main__':
    s=b=0
    for args in [(13,4,1,[0,1]),(17,2,2,[0,1,2]),(37,0,5,[0,1,3,4]),(43,2,1,[0,1,2,3])]:
        x,y=recurrence(*args);s+=x;b+=y
    print(f'totals: {s} sequential, {b} balanced comparisons',flush=True)
    aggregate();spectrum()
    print('PASS: homogeneous resultant recursion, collective evaluation and spectrum.',flush=True)
