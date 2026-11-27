#!/usr/bin/env python3
"""Exact finite checks for abstract additive range-avoidance boundaries."""
from itertools import combinations, product
from math import prod
from random import Random

if not __debug__:
    raise RuntimeError('Run with assertions enabled.')


def add(a,b,p):
    return ((a[0]+b[0])%p,(a[1]+b[1])%p)


def sums(A,p):
    S={(0,0)}
    for B in A:
        S={add(a,b,p) for a in S for b in B}
    return S


def mask(S,p):
    return sum(1<<(a*p+b) for a,b in S)


def lower_cover(H,sizes,p):
    assert 0<len(H)<=sum(n-1 for n in sizes)+1
    H=list(H); origin=H[0]
    others=[((x-origin[0])%p,(y-origin[1])%p) for x,y in H[1:]]
    A=[]; i=0
    for n in sizes:
        B={(0,0)} | set(others[i:i+n-1]);i+=n-1
        A.append(B)
    assert i>=len(others)
    A[0]={add(a,origin,p) for a in A[0]}
    allpts=list(product(range(p),repeat=2))
    for B,n in zip(A,sizes):
        for x in allpts:
            if len(B)==n:break
            B.add(x)
        assert len(B)==n
    assert set(H)<=sums(A,p)
    return A


def check_universal():
    print('UNIVERSAL CANDIDATE SETS (exhaustive small cases)')
    print('prime | sizes | distinct sumsets | covered subsets checked | minimum | lower | upper (capped)')
    for p,sizes in [(3,(2,2)),(3,(2,3)),(3,(2,2,2)),(5,(2,2))]:
        Q=p*p; pts=list(product(range(p),repeat=2)); N=sum(sizes); M=prod(sizes)
        assert M<Q
        choices=[[{(0,0)}|set(t) for t in combinations(pts[1:],n-1)] for n in sizes]
        supports=set()
        for A in product(*choices):
            S=sums(A,p)
            for t in pts:
                supports.add(mask({add(x,t,p) for x in S},p))
        low=N-len(sizes)+2
        high=1
        while high<Q and Q**(N-len(sizes)+1)*M**high>=Q**high:
            high+=1
        checked=0
        for k in range(1,low):
            for H in combinations(pts,k):
                lower_cover(H,sizes,p)
                checked+=1
        optimal=None
        for k in range(low,high+1):
            for H in combinations(pts,k):
                h=mask(H,p)
                if not any(h&s==h for s in supports):
                    optimal=k;break
            if optimal is not None:break
        assert optimal is not None and low<=optimal<=high
        print(f'{p} | {sizes} | {len(supports)} | {checked} | {optimal} | {low} | {high}')


def fval(t,coefs,p):
    ans=[]
    for c in coefs:
        v=0
        for a in reversed(c):v=(v*t+a)%p
        ans.append(v)
    return tuple(ans)


def component(I,z,base,coefs,p):
    ans=(0,0)
    for k in range(len(I)+1):
        sign=(-1)**(len(I)-k)
        for J in combinations(range(len(I)),k):
            t=sum(z[j]*base**I[j] for j in J)
            x,y=fval(t,coefs,p)
            ans=((ans[0]+sign*x)%p,(ans[1]+sign*y)%p)
    return ans


def check_polynomial():
    print('\nLOW-DEGREE PARAMETER COVER')
    print('prime | degree | radix | digits | blocks | max block size | parameter checks')
    rng=Random(20261127); total=0; mixed=0
    cases=[(101,1,5,3),(65537,2,3,4),(65537,2,5,3),
           (1000003,3,3,3),(2147483647,3,2,4),(2147483647,4,2,4)]
    for p,d,q,s in cases:
        I=[I for k in range(1,min(d,s)+1) for I in combinations(range(s),k)]
        n=q**d; assert n**len(I)<p*p
        # Polynomial map with first coordinate t, and a second independent map.
        maps=[((0,1),tuple(rng.randrange(p) for _ in range(d))+(1,)),
              tuple(tuple(rng.randrange(p) for _ in range(d+1)) for _ in range(2))]
        for coefs in maps:
            C=[{component(i,z,q,coefs,p) for z in product(range(q),repeat=len(i))} for i in I]
            assert all(1<=len(B)<=n for B in C)
            constant=fval(0,coefs,p)
            C[0]={add(x,constant,p) for x in C[0]}
            for z in product(range(q),repeat=s):
                a=[]
                for j,i in enumerate(I):
                    x=component(i,tuple(z[k] for k in i),q,coefs,p)
                    if j==0:x=add(x,constant,p)
                    assert x in C[j];a.append(x)
                result=(0,0)
                for x in a:result=add(result,x,p)
                assert result==fval(sum(x*q**i for i,x in enumerate(z)),coefs,p)
                total+=1
            for k in range(d+1,s+1):
                for i in combinations(range(s),k):
                    for z in product(range(q),repeat=k):
                        assert component(i,z,q,coefs,p)==(0,0)
                        mixed+=1
            print(f'{p} | {d} | {q} | {s} | {len(I)} | {max(map(len,C))} | {q**s}')
    print(f'{total} exact parameter representations; {mixed} vanishing higher mixed differences.')


def arc_length(values,p):
    v=sorted(set(values)); assert v
    gaps=[b-a for a,b in zip(v,v[1:])]+[p+v[0]-v[-1]]
    return p-max(gaps)+1


def min_rect(sample,p):
    # Distinct projective directions, and all nonzero scalings within each.
    minima=[]; tested=0
    for a,b in [(1,k) for k in range(p)]+[(0,1)]:
        vals={ (a*x+b*y)%p for x,y in sample }
        best=p
        for c in range(1,p):
            best=min(best,arc_length((c*x%p for x in vals),p));tested+=1
        minima.append(best)
    minima.sort()
    return minima[0]*minima[1],tested


def axis_digits(p,A,i):
    return {(j*pow(A,i,p)%p,0) for j in range(A)} | {(0,j*pow(A,i,p)%p) for j in range(A)}


def check_boxes():
    print('\nFULL-PROJECTION SINGLE-BOX OBSTRUCTIONS')
    print('prime | sizes | minimal box areas | product bound | actual sumset size | affine linear forms checked')
    rng=Random(20261127);total=0
    for p,sample_size in [(101,60),(127,70),(257,90)]:
        m=2;A=1
        while A**m<p:A+=1
        pts=list(product(range(p),repeat=2)); L=[];areas=[];forms=0
        for i in range(m):
            S=set(rng.sample(pts,sample_size))|axis_digits(p,A,i)
            area,n=min_rect(S,p)
            assert 2*area>p*p,(p,area)
            areas.append(area);forms+=n;L.append(S)
        M=prod(map(len,L));assert M<p*p
        Z=sums(L,p)
        assert len(Z)<=M
        assert all((x,0) in Z and (0,x) in Z for x in range(p))
        for a,b in [(1,k) for k in range(p)]+[(0,1)]:
            assert {(a*x+b*y)%p for x,y in Z}==set(range(p))
        total+=forms
        print(f'{p} | {tuple(map(len,L))} | {tuple(areas)} | {M} | {len(Z)} | {forms}')
    print(f'{total} exact linear-form arc checks; all projective scalar projections are full.')


def check_combined():
    # NumPy accelerates the exact arc scan;
    # all operations are bounded integers and the final support is enumerated.
    import numpy as np
    p=997;A=10;q=3;m=3;I=[(0,),(1,),(0,1)]
    coefs=((0,1),(5,3,1));rng=Random(20261127)
    L=[];areas=[];scales=np.arange(1,p,dtype=np.int64)[:,None]
    for i in range(m):
        S={divmod(j,p) for j in rng.sample(range(p*p),65)}|axis_digits(p,A,i)
        C={component(I[i],z,q,coefs,p) for z in product(range(q),repeat=len(I[i]))}
        if i==0:C={add(x,fval(0,coefs,p),p) for x in C}
        S|=C;L.append(S)
        P=np.array(sorted(S),dtype=np.int64);minima=[]
        for a,b in [(1,k) for k in range(p)]+[(0,1)]:
            vals=(P[:,0]*a+P[:,1]*b)%p
            rows=np.sort((scales*vals[None,:])%p,axis=1)
            gaps=np.maximum(np.max(np.diff(rows,axis=1),axis=1),p+rows[:,0]-rows[:,-1])
            minima.append(int(np.min(p-gaps+1)))
        minima.sort();area=minima[0]*minima[1];assert 2*area>p*p;areas.append(area)
    M=prod(map(len,L));assert M<p*p
    Z=sums(L,p)
    assert all((x,0) in Z and (0,x) in Z for x in range(p))
    assert all(fval(t,coefs,p) in Z for t in range(q*q))
    print('\nCOMBINED THREE-BLOCK CONTROL')
    print(f'r={p}, sizes={tuple(map(len,L))}, product={M}, support={len(Z)}, minimum box areas={tuple(areas)}.')
    print(f'{m*(p*p-1)} exact linear-form arc checks; both axes and all {q*q} quadratic prefix points present.')


def cyclic_intervals(p):
    return {tuple(sorted((a+j)%p for j in range(k))) for a in range(p) for k in range(1,p+1)}


def check_multi_box():
    # Standalone union-of-boxes example at p=7 (not a sparse-product fixture).
    p=7;pts=list(product(range(p),repeat=2));half=p*p//2
    first={(0,0),(0,2),(0,5),(1,1),(1,4),(1,6),(2,0),(2,3),(2,5),(3,1),(3,4),(3,6),
           (4,0),(4,2),(4,5),(5,1),(5,3),(5,6),(6,0),(6,2),(6,4),(6,6)}
    S={(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(2,1),(2,2),(2,5),(3,2),(3,4),
       (3,5),(3,6),(4,2),(4,3),(4,4),(4,6),(5,0),(5,6),(6,1),(6,4),(6,5)}
    assert len(S)==len(first)==22
    # Enumerate all affine boxes of area <= half, with each box's own frame.
    intervals=cyclic_intervals(p);forms=[(a,b) for a,b in pts if (a,b)!=(0,0)]
    projected={u:{I:mask({x for x in pts if (u[0]*x[0]+u[1]*x[1])%p in I},p) for I in intervals} for u in forms}
    allboxes={0}; orientations=0
    for u,v in combinations(forms,2):
        if (u[0]*v[1]-u[1]*v[0])%p==0:continue
        orientations+=1
        for I,a in projected[u].items():
            for J,b in projected[v].items():
                if len(I)*len(J)<=half:allboxes.add(a&b)

    def check(sample):
        s=mask(sample,p); n=len(sample)
        possible=sorted(((z&s).bit_count(),z) for z in allboxes if (z&~s).bit_count()<=half-n)
        # A covering pair must have combined counts on S >= |S|. Skip every
        # smaller pair, and include identical boxes to cover the one-box case.
        from bisect import bisect_left
        counts=[k for k,z in possible]; pairs=0
        for i,(k,a) in enumerate(possible):
            for j in range(max(i,bisect_left(counts,n-k)),len(possible)):
                pairs+=1;b=possible[j][1];z=a|b
                if z&s==s and z.bit_count()<=half:
                    return (len(possible),pairs,z.bit_count())
        return len(possible),pairs,None

    prior=check(first); final=check(S)
    assert prior[2] is not None and final[2] is None
    print('\nTWO-BOX UNION CONTROL (standalone)')
    print(f'r={p}, sample={len(S)}, half-area={half}, independent form pairs={orientations}, distinct boxes={len(allboxes)}.')
    print(f'Initial patterned sample: found a two-box cover of area {prior[2]}.')
    print(f'Replacement sample: {final[0]} candidate boxes, {final[1]} eligible pairs examined; no cover of area <= half.')


if __name__=='__main__':
    check_universal();check_polynomial();check_boxes();check_combined();check_multi_box();print('\nAll checks passed.')
