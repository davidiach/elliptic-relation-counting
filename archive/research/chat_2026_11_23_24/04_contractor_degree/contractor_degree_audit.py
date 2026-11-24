#!/usr/bin/env python3
"""Finite regressions for the degree-sensitive contractor note; stdlib only."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
from random import Random


def v2(x):
    assert x != 0
    x = abs(x)
    return (x & -x).bit_length() - 1


def mul(a,b,mod=None):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    if mod: c=[x%mod for x in c]
    return c


def poly_pow(a,n,mod):
    b=[1]
    while n:
        if n&1: b=mul(b,a,mod)
        n//=2
        if n: a=mul(a,a,mod)
    return b


def collision():
    checked=expanded=coeff_count=0
    # Direct factor valuations for all centers. Explicitly expand a smaller
    # subset of the normalized polynomials as a separate modular check.
    for a in range(1,8):
        L=2**a; grid=range(-L,L); alpha=2*L
        for u in range(1,a+1):
            for M in (8*u,16*u,32*u):
                A=M//4; h=M//(4*u)
                orders={s:A+h*(u-v2(s)) for s in grid if s and v2(s)<u}
                deg=sum(orders.values()); assert deg<=M*(2*L-1)
                for t in grid:
                    order=orders.get(t,0)
                    val=0
                    for s,j in orders.items():
                        val+=j*((a+1 if s==t else v2(t-s))-v2(s))
                    assert val>=0
                    if t in orders:
                        i=v2(t); D=L//(2**i)
                        assert order>=A and 2*val>=h*D
                    else: assert val==0
                    checked+=1
                    if a<=3 and M==8*u:
                        # Expand the unit-normalized local polynomial directly
                        # in Z/2^P. Scaling out each factor's minimum valuation
                        # makes its two coefficients integral with odd denominator.
                        P=max(10,min(80,val+8)); mod=1<<P
                        poly=[1]
                        for s,j in orders.items():
                            v = a+1 if s==t else v2(t-s)
                            den=(-s)//(2**v2(s))
                            inv=pow(den,-1,mod)
                            f=[((t-s)//(2**v))*inv%mod,(alpha//(2**v))*inv%mod]
                            poly=mul(poly,poly_pow(f,j,mod),mod)
                        poly=[(x<<val)%mod for x in poly]
                        assert all(x==0 for x in poly[:order])
                        if val<P:
                            assert v2(poly[order])==val
                        else: assert not any(poly)
                        if t in orders:
                            assert all(x%(1<<min(P,(h*(L//(2**v2(t))))//2))==0 for x in poly)
                        if t==0: assert poly[0]==1
                        expanded+=1; coeff_count+=len(poly)
    print(f'Joint obstruction: {checked:,} center checks; {expanded:,} explicit local polynomials; {coeff_count:,} coefficients.')


def enum_bound():
    cases=0
    # Binary subsets of the chosen low-weight raw modes meet both budgets.
    for a in range(8,25):
        L=2**a; M=4*L
        for u in range(1,a//2+1):
            h=M//(4*u); A=M//4
            for i in range(u):
                D=L//(2**i)
                k=int(D/(32*u**0.5))
                if not k: continue
                tiles=M//L*(2**(i+1))
                K=tiles*k; W=tiles*k*(k+1)//2
                assert K<=A and 2*W<=h*D
                cases+=1
    print(f'Joint mode budgets: {cases:,} exact degree/weight inequalities.')


def low_degree():
    rng=Random(20261124)
    checks=0
    for L in (2,4,8):
        nodes=list(range(-L,L))
        c=1
        for t in nodes:
            if t: c*= -t
        W={}
        for d in range(1-L,L):
            poly=[1]
            for t in nodes:
                if t!=d: poly=mul(poly,[-t,1])
            W[d]=[Fraction(x,c) for x in poly]
        for _ in range(8):
            # Row-side parity functions on a partition, whose mod-2 lifts
            # have degree one each, even when the literal 0/1 lift has more.
            groups=[list(range(0,L//2)),list(range(L//2,L))]
            tags=[[rng.randrange(2) for _ in g] for g in groups]
            offsets=[rng.randrange(2) for _ in groups]
            F=[Fraction(0)]*(2*(2*L-1)+1)
            truth=0
            for X in product(range(2),repeat=L):
                b=rng.randrange(2)
                if not b: continue
                truth+=int(all((c0+sum(t*X[j] for t,j in zip(ts,g)))%2 for c0,ts,g in zip(offsets,tags,groups)))
                s=[]
                for x in range(L):
                    s.append([sum(X[y]*W[y-x][z] for y in range(L)) for z in range(2*L)])
                f=[Fraction(1)]
                for c0,ts,g in zip(offsets,tags,groups):
                    a=[sum(t*s[j][z] for t,j in zip(ts,g)) for z in range(2*L)];a[0]+=c0
                    f=mul(f,a)
                assert len(f)<=len(F)
                for j,x in enumerate(f):F[j]+=x
            assert F[0].denominator==1 and F[0].numerator%2==truth%2
            assert len(F)-1<=2*(2*L-1)
            checks+=1
    print(f'Low-degree lifts: {checks} exact rational deformations and parity comparisons.')


def support_degree():
    n=4; count=0
    # Exhaustive nonzero Boolean truth tables; Möbius transform in F_2.
    for f in range(1,1<<(1<<n)):
        a=[(f>>j)&1 for j in range(1<<n)]
        for i in range(n):
            for j in range(1<<n):
                if j&(1<<i): a[j]^=a[j^(1<<i)]
        d=max(j.bit_count() for j,c in enumerate(a) if c)
        assert f.bit_count()>=(1<<(n-d))
        count+=1
    print(f'Support-degree inequality: all {count:,} nonzero 4-bit truth tables.')


def ring_mul(x,y,Q):return ((x[0]*y[0]-x[1]*y[1])%Q,(x[0]*y[1]+x[1]*y[0])%Q)
def ring_pow(x,n,Q):
    r=(1,0)
    while n:
        if n&1:r=ring_mul(r,x,Q)
        x=ring_mul(x,x,Q);n//=2
    return r


def fourier_rank():
    total=0
    for l,q,B,h in ((5,19,4,5),(7,83,3,3),(11,43,2,4)):
        Q=q**h
        # q=3 mod 4; the unramified ring is Z/q^h[I]/(I^2+1).
        root=None
        for x,y in product(range(q),repeat=2):
            if not (x or y):continue
            z=ring_pow((x,y),(q*q-1)//l,q)
            if z!=(1,0):root=z;break
        assert root and ring_pow(root,l,q)==(1,0)
        # Teichmuller-lift primitive root by q-powering and removing its
        # residue exponent; q^(h-1) kills the principal-unit component.
        root=ring_pow(root,q**(h-1)*pow(q**(h-1),-1,l),Q)
        assert ring_pow(root,l,Q)==(1,0) and root!=(1,0)
        assert (root[0],-root[1]%Q)==ring_pow(root,l-1,Q)
        mid=(l-1)//4
        A=[0]+list(range(1,mid+1)); C=list(range(mid+1,(l+1)//2))
        smin=2*min(len(A)-1,len(C))
        assert Q**smin>(2*l*B)**l
        for lam in (1,2):
            theta=pow(lam,pow(l,-1,q**(h-1)*(q-1)),Q)
            assert pow(theta,l,Q)==lam
            points=[ring_mul((theta,0),ring_pow(root,j,Q),Q) for j in A+C]
            tables={}; cols=Counter()
            for a in product(range(B),repeat=l):
                vals=[]
                for z in points:
                    v=(0,0)
                    for c in reversed(a):
                        v=ring_mul(v,z,Q);v=((v[0]+c)%Q,v[1])
                    vals.append(v)
                x,y=tuple(vals[:len(A)]),tuple(vals[len(A):])
                assert x not in tables;tables[x]=y;cols[y]+=1
            assert max(cols.values())<=(B if lam==1 else 1)
            # Each occupied A-row has one entry, so rank = # nonzero columns
            # over every field, without constructing the enormous padded matrix.
            rank=len(cols);assert rank>=B**(l-1)
            print(f'Fourier row l={l}, B={B}, lambda={lam}, Q={Q}: {len(tables):,} legal rows; rank {rank:,}; max fibre {max(cols.values())}.')
            total+=len(tables)
    print(f'Fourier projection: {total:,} legal assignments checked.')


def middle_image():
    rng=Random(20261124);checks=assignments=0
    for q,h,f in ((3,2,1),(5,2,1),(7,2,1),(3,2,2),(7,2,2),(11,1,2)):
        Q=q**h
        if f==1:
            vals=[(x,0) for x in range(Q)]
        else:
            assert q%4==3
            vals=list(product(range(Q),repeat=2))
        add=lambda x,y:((x[0]+y[0])%Q,(x[1]+y[1])%Q)
        neg=lambda x:((-x[0])%Q,(-x[1])%Q)
        sub=lambda x,y:add(x,neg(y))
        mul=lambda x,y:ring_mul(x,y,Q)
        sq=lambda x:mul(x,x)
        sc=lambda k,x:((k*x[0])%Q,(k*x[1])%Q)
        def s3(u,v,w,a,c):
            return add(sub(sub(mul(sq(sub(u,v)),sq(w)),
              sc(2,mul(add(mul(add(u,v),add(mul(u,v),a)),sc(2,c)),w))),
              sc(4,mul(c,add(u,v)))),sq(sub(mul(u,v),a)))
        for _ in range(6):
            v,a,c=(rng.choice(vals) for j in range(3));w=add(v,(1,0))
            # Derive coefficient of U via evaluations at 0 and 1; the
            # leading coefficient is (v-w)^2=1, also checked pointwise.
            C=s3((0,0),v,w,a,c)
            B=sub(sub(s3((1,0),v,w,a,c),C),(1,0))
            outputs=set();num=0
            for u in vals:
                actual=s3(u,v,w,a,c)
                assert actual==add(add(sq(u),mul(B,u)),C)
                derivative=add(sc(2,u),B)
                if derivative[0]%q or derivative[1]%q:
                    outputs.add(actual);num+=1
            assert num==len(vals)*(q**f-1)//(q**f)
            assert len(outputs)*2>=num
            checks+=1;assignments+=len(vals)
    print(f'Middle quadratic images: {checks} fixtures; {assignments:,} assignments; every simple-root image bound holds.')


def xor_spectrum():
    print('XOR diagnostic (width, median/max Walsh eigenvalue valuation at t=L+1):')
    for a in (4,6,8,10,12):
        L=2**a;t=L+1
        v=[(-1)**(L-1-d)*comb(t,d)*comb(t-d-1,L-1-d) for d in range(0,L,2)]
        for j in range(a-1):
            s=2**j
            for i in range(0,len(v),2*s):
                for k in range(i,i+s): v[k],v[k+s]=v[k]+v[k+s],v[k]-v[k+s]
        assert all(v)
        z=sorted(map(v2,v)); print(f'  {L}: {z[len(z)//2]}/{z[-1]}')

if __name__=='__main__':
    assert __debug__
    collision();enum_bound();low_degree();support_degree();fourier_rank();middle_image();xor_spectrum()
