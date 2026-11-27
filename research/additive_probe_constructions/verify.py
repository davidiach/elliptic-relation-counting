#!/usr/bin/env python3
"""Exact finite checks for additive probe constructions (abstract finite planes)."""
from itertools import product, permutations, combinations
from math import comb, isqrt, prod
from random import Random

if not __debug__:
    raise RuntimeError('Assertions must be enabled')


def add(a, b, p):
    return ((a[0]+b[0])%p, (a[1]+b[1])%p)


def sub(a, b, p):
    return ((a[0]-b[0])%p, (a[1]-b[1])%p)


def sumsets(sets, p):
    S = {(0,0)}
    for B in sets:
        S = {add(a,b,p) for a in S for b in B}
    return S


def rank(rows, p):
    piv = {}
    for row in rows:
        a = [x%p for x in row]
        for j,b in sorted(piv.items()):
            if a[j]:
                c = a[j]
                a = [(x-c*y)%p for x,y in zip(a,b)]
        j = next((i for i,x in enumerate(a) if x),None)
        if j is not None:
            c = pow(a[j],-1,p)
            piv[j] = [x*c%p for x in a]
    return len(piv)


def incidence(sizes):
    return [(1,)+tuple(int(x==j) for x,n in zip(a,sizes) for j in range(1,n))
            for a in product(*(range(n) for n in sizes))]


def radix_probes(sizes,p):
    D=1+sum(n-1 for n in sizes)
    C=1 if len(sizes)==2 else isqrt((len(sizes)+1)**D)
    B=C+1; k=(D+1)//2
    assert B**k<=p
    H=[(0,0)]+[(B**i,0) for i in range(k)]+[(0,B**i) for i in range(D-k)]
    assert len(set(H))==D+1
    return H


def check_radix():
    print('EXACT-MINIMUM RADIX PROBES')
    rng=Random(20261127); total=0
    for p,sizes,exhaust in [(5,(2,2),True),(11,(2,3),True),(11,(3,3),True),
                             (293,(2,2,2),True),(175633,(2,2,2,2),False)]:
        assert all(p%i for i in range(2,isqrt(p)+1))
        H=radix_probes(sizes,p); R=incidence(sizes); K=len(H)
        seqs=permutations(R,K) if exhaust else (rng.sample(R,K) for _ in range(1000))
        count=0
        for A in seqs:
            assert rank([a+h for a,h in zip(A,H)],p)>rank(A,p)
            count+=1
        total+=count
        print(f'r={p} sizes={sizes} K={K} witness arrays={count} exhaustive={exhaust}')
    print(f'{total} witness arrays ruled out exactly.')


def sidon_bound(P,Q):
    return min(P*Q,(Q+isqrt(Q*Q+4*Q*P*(P-1)))//2,
               (P+isqrt(P*P+4*P*Q*(Q-1)))//2)


def check_sidon():
    print('\nSIDON GRAPH PROBES')
    rng=Random(918); total=0
    for p,P,Q in [(5,2,2),(7,2,3),(11,3,4),(31,6,6),(101,12,12)]:
        K=sidon_bound(P,Q)+1
        assert K<=p and P*Q<p*p
        H={(t,t*t%p) for t in range(K)}
        pts=list(product(range(p),repeat=2)); count=0; maxhit=0
        # All translations can be tested by moving a normalized sumset.
        if P==Q==2:
            choices=( ([(0,0),a],[(0,0),b]) for a in pts[1:] for b in pts[1:] )
        else:
            choices=((rng.sample(pts,P),rng.sample(pts,Q)) for _ in range(250))
        for A,B in choices:
            S=sumsets([A,B],p)
            shifts=pts if P==Q==2 else [rng.choice(pts)]
            for z in shifts:
                n=sum(sub(x,z,p) in S for x in H)
                assert n<K
                maxhit=max(maxhit,n); count+=1
        total+=count
        print(f'r={p} sizes={P,Q} K={K} max hit={maxhit} checks={count}')
    print(f'{total} exact cover exclusions.')


def subset_collision(T,k,p):
    seen={}
    for I in combinations(range(len(T)),k):
        z=(sum(T[i][0] for i in I)%p,sum(T[i][1] for i in I)%p)
        if z in seen:
            A=set(seen[z]); B=set(I)
            return sorted(A-B),sorted(B-A)
        seen[z]=I
    return None


def cycle_lists(U,V,p):
    assert len(U)==len(V)>=2
    a=(0,0); A=[a]; B=[]
    for i,(u,v) in enumerate(zip(U,V)):
        b=sub(u,a,p); B.append(b); a=sub(v,b,p)
        if i+1<len(U): A.append(a)
    assert a==(0,0) and len(A)==len(B)==len(U)
    return set(A),set(B)|{(0,0)}


def threshold(p,k):
    s=k
    while comb(s,k)<=p*p: s+=1
    return s


def cover_pair(T,n,k,p):
    T=list(T); s=threshold(p,k); A={(0,0)};B={(0,0)}; cost=[0,0]; q=0
    while len(T)>=s:
        I,J=subset_collision(T,k,p)
        assert len(I)==len(J)>=2
        U=[T[i] for i in I]; V=[T[i] for i in J]
        C,D=cycle_lists(U,V,p); j=len(I)
        if q%2: A|=D;B|=C;cost[0]+=j;cost[1]+=j-1
        else: A|=C;B|=D;cost[0]+=j-1;cost[1]+=j
        q+=1; gone=set(I+J); T=[v for i,v in enumerate(T) if i not in gone]
    assert sum(cost)+len(T)<=2*(n-1) and max(cost)<=n-1
    for v in T:
        if cost[0]<n-1: A.add(v);cost[0]+=1
        else: B.add(v);cost[1]+=1
    assert len(A)<=n and len(B)<=n
    return A,B,q


def check_cycles():
    print('\nCYCLE-PACKING COVER LOWER BOUND')
    rng=Random(606); total=0; cycles=0
    for p,m,n,k in [(31,2,20,3),(101,2,50,4),(101,3,20,5),(1009,4,30,7)]:
        s=threshold(p,k); g=max(0,(2*(n-1)-s+1)//(2*k)); cap=2*(n-1)+g
        assert n**m<p*p and g>0,(p,m,n,k,s)
        K=1+m*(n-1)+(m//2)*g
        for _ in range(20):
            H=[(0,0)]+[divmod(x,p) for x in rng.sample(range(1,p*p),K-1)]
            shift=divmod(rng.randrange(p*p),p)
            sets=[]; index=1
            for j in range(m//2):
                chunk=H[index:index+cap]; index+=cap
                A,B,q=cover_pair(chunk,n,k,p)
                assert set(chunk)<=sumsets([A,B],p)
                assert q>=g;cycles+=q; sets.extend([A,B])
            if m%2: sets.append(set(H[index:])|{(0,0)})
            assert len(sets)==m and all(len(A)<=n and (0,0) in A for A in sets)
            # Each chunk is covered with every other block at zero.
            # Translate the first list; the same witnesses cover the translated H.
            sets[0]={add(a,shift,p) for a in sets[0]}
            assert all(len(A)<=n for A in sets) and shift in sets[0]
            if m==2: assert {add(a,shift,p) for a in H}<=sumsets(sets,p)
            total+=1
        print(f'r={p} m={m} n={n} k={k} threshold={s} extra/pair={g} cover size={K}')
    print(f'{total} complete constructed covers; {cycles} disjoint cycles.')


def ceil_root(x,k):
    lo,hi=0,1
    while hi**k<x: hi*=2
    while lo+1<hi:
        t=(lo+hi)//2
        if t**k>=x: hi=t
        else: lo=t
    return hi


def basis(p,h):
    q=ceil_root(p*p,h);k=h//2
    A=[{(a*q**j%p,0) for a in range(q)} for j in range(k)]
    A+=[{(0,a*q**j%p) for a in range(q)} for j in range(k)]
    if h%2:
        z=q**k; s=(p+z-1)//z
        A.append({(a*z%p,b*z%p) for a,b in product(range(s),repeat=2)})
    return A,set.union(*A)


def capped_growth(A,p,h,cap):
    S={(0,0)}
    for _ in range(h):
        T=set()
        for a in S:
            for b in A:
                T.add(add(a,b,p))
                if len(T)>cap:return None
        S=T
    return len(S)


def check_growth():
    print('\nGROWTH PROBES AND GROUPING')
    bases=0;checks=0
    for p in (7,11,23,101):
        for h in range(2,9):
            B,H=basis(p,h)
            assert len(B)==h and len(sumsets(B,p))==p*p
            assert len(H)<=h*ceil_root(p*p,h)+2*isqrt(ceil_root(p*p,h))+3
            assert comb(len(H)+h-1,h)>=p*p
            bases+=1
    print(f'{bases} exact plane-basis constructions.')
    for p,m,n,h in [(101,3,4,6),(1009,4,6,4),(2003,3,4,6)]:
        A=[{(j,j*(i+1)%p) for j in range(n)} for i in range(m)]
        exact=[capped_growth(B,p,h,100000) for B in A]
        assert all(x is not None for x in exact) and prod(exact)<p*p
        S=sumsets(A,p);_,H=basis(p,h)
        assert H-S
        print(f'r={p} sizes={(n,)*m} h={h} growth product={prod(exact)} probes={len(H)} missing={len(H-S)}')
        checks+=1
    p=101;h=4
    A=[{(j,0) for j in range(5)}]*2+[{(0,j) for j in range(5)}]*2
    ungrouped=prod(capped_growth(B,p,h,100000) for B in A)
    grouped=prod(capped_growth(sumsets(A[i:i+2],p),p,h,100000) for i in (0,2))
    assert ungrouped>=p*p>grouped and basis(p,h)[1]-sumsets(A,p)
    assert capped_growth(A[0],p,h,16) is None
    print(f'grouping rescue: r={p} ungrouped={ungrouped} grouped={grouped}')
    p=101;A0=ceil_root(p,3)
    A=[{(j*A0**i%p,0) for j in range(A0)}|{(0,j*A0**i%p) for j in range(A0)} for i in range(3)]
    S=sumsets(A,p)
    assert prod(map(len,A))<p*p and all((x,0) in S and (0,x) in S for x in range(p))
    assert len(sumsets([S,S],p))==p*p
    for a,b in [(1,2),(3,p-1),(p-1,p-1)]:
        assert len(sumsets([{(a*x%p,a*y%p) for x,y in S},{(b*x%p,b*y%p) for x,y in S}],p))==p*p
    print(f'honest fast-growth obstruction: r={p}, product={prod(map(len,A))}, |S|={len(S)}, |2S|={p*p}')
    print(f'{checks} growth certificates, one grouping rescue and one obstruction.')


if __name__=='__main__':
    check_radix();check_sidon();check_cycles();check_growth()
    print('\nAll checks passed.')
