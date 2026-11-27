"""Exact abstract finite-group checks. Standard library + SymPy for independent HNF.
Optional: supply PR36's check_joint_shadows.py to test the previous selector.
"""
from itertools import product, combinations
from math import gcd, prod, comb
from random import Random
from collections import Counter
from pathlib import Path
import sys, importlib.util
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form


def add_pivot(row, piv, q):
    row = [x % q for x in row]
    for j in sorted(piv):
        if row[j]:
            c = row[j]
            row = [(x-c*y)%q for x,y in zip(row,piv[j])]
    j=next((j for j,x in enumerate(row) if x), None)
    if j is None:return False
    c=pow(row[j],-1,q);piv[j]=[x*c%q for x in row]
    return True


def basis_indices(A,q):
    piv={};B=[]
    for j,a in enumerate(A):
        if add_pivot((1,)+a,piv,q):B.append(j)
    return B


def contains(A,x,q):
    piv={}
    for a in A:add_pivot((1,)+a,piv,q)
    return not add_pivot((1,)+x,piv,q)


def peel(A,qs,e):
    """No core labels are used. After peeling, restore inside every remaining hull."""
    R=list(A)
    for q in qs:
        for _ in range(e):
            B=set(basis_indices(R,q))
            R=[a for i,a in enumerate(R) if i not in B]
    if not R:return [[a] for a in A]
    pivots=[]
    for q in qs:
        piv={}
        for a in R:add_pivot((1,)+a,piv,q)
        pivots.append(piv)
    core=[];out=[]
    for a in A:
        inside=True
        for q,piv in zip(qs,pivots):
            if add_pivot((1,)+a,piv.copy(),q):inside=False;break
        (core if inside else out).append(a)
    return [core]+[[a] for a in out]


def crt_pair(x,y):return (10*x+6*y)%15


def planted_tests():
    rng=Random(648921);qs=(3,5);total=restored=removed=0
    for trial in range(80):
        d=6+(trial%3);e=1+trial%3
        # Independent rank-two affine spaces, lifted by CRT into [0,15)^d.
        V=[]
        for q in qs:
            while True:
                u=tuple(rng.randrange(q) for _ in range(d))
                v=tuple(rng.randrange(q) for _ in range(d))
                piv={}
                if add_pivot(u,piv,q) and add_pivot(v,piv,q):break
            a=tuple(rng.randrange(q) for _ in range(d))
            V.append([tuple((c+x*b+y*f)%q for c,b,f in zip(a,u,v))
                      for x,y in product(range(q),repeat=2)])
        core=[tuple(crt_pair(x,y) for x,y in zip(u,v)) for u,v in product(*V)]
        assert len(set(core))==225
        core=rng.sample(core,140+trial%50)
        out=[]
        while len(out)<e:
            a=tuple(rng.randrange(17) for _ in range(d))
            if a not in out and a not in core and any(not contains(core,a,q) for q in qs):out.append(a)
        A=core+out;rng.shuffle(A)
        partition=peel(A,qs,e)
        assert sorted(sum(partition,[]))==sorted(A)
        assert len(partition)-1<=len(qs)*e*(d+1)
        assert all(contains(core,a,q) for a in partition[0] for q in qs)
        restored+=set(partition[0])==set(core)
        removed+=len(partition)-1;total+=len(A)
    print(f'Planted affine cores: 80 exact partitions, {total} input words; {restored} complete core recoveries; {removed} singleton exceptions after restoration.')


def closure2(C,p,q):
    a,b=C[0];D=0
    dif=[(x-a,y-b) for x,y in C]
    for u,v in dif:
        for x,y in dif:D=gcd(D,u*y-v*x)
    lo,hi=min(x for x,y in C),max(x for x,y in C)
    lj,hj=min(y for x,y in C),max(y for x,y in C)
    out=set()
    if D:
        for x in range(lo,hi+1):
            for y in range(lj,hj+1):
                if all(((x-a)*v-(y-b)*u)%D==0 for u,v in dif):out.add((x,y))
    else:
        v=next(((x,y) for x,y in dif if x or y),None)
        if v is None:return {(a,b)}
        u,v=v;g=gcd(u,v);u//=g;v//=g;G=0
        for x,y in dif:G=gcd(G,x//u if u else y//v)
        for x in range(lo,hi+1):
            for y in range(lj,hj+1):
                X,Y=x-a,y-b
                if X*v==Y*u and (X//u if u else Y//v)%G==0:out.add((x,y))
    return out


def hnf(C):
    a=C[0];d=len(a)
    return hermite_normal_form(Matrix(d,len(C),lambda i,j:C[j][i]-a[i]))


def hnf_member(H,v):
    # HNF columns have their last nonzero entries in increasing rows.
    x=list(v)
    for j in reversed(range(H.cols)):
        col=[int(H[i,j]) for i in range(H.rows)]
        i=max(i for i,z in enumerate(col) if z)
        q,r=divmod(x[i],col[i])
        if r:return False
        x=[a-q*b for a,b in zip(x,col)]
    return not any(x)


def closure_hnf(C):
    H=hnf(C);a=C[0];d=len(a)
    ranges=[range(min(x[j] for x in C),max(x[j] for x in C)+1) for j in range(d)]
    return {x for x in product(*ranges) if hnf_member(H,[u-v for u,v in zip(x,a)])}


def closure_tests():
    # A is an exact small counterexample to two-class partition avoidance.
    A=[(10,10),(5,3),(5,6),(2,5),(7,11),(0,3),(9,0),(4,2),(5,1),(6,8),(8,11)]
    p,q=11,13;n=len(A);full=(1<<n)-1;cache={};members=0
    for mask in range(1,full+1):
        C=[a for i,a in enumerate(A) if mask>>i&1]
        K=closure2(C,p,q);assert K==closure_hnf(C);assert set(C)<=K
        for a in K:
            for Q in (2,3,5,7,11,13,17):assert contains(C,a,Q)
        members+=len(K)
        cache[mask]=K
    sizes=[len(cache[full])]
    for mask in range(1,full,2):sizes.append(len(cache[mask]|cache[full^mask]))
    assert min(sizes)==73 and 2*min(sizes)>p*q and n*n<p*q
    S={( (x+u)%p,(y+v)%q ) for x,y in A for u,v in A}
    # Both closures for any choices of partitions have size>half, so all targets
    # are present in some branch even with every prime enforced. Check representative pairs.
    rng=Random(4106);count=0
    for _ in range(100):
        a,b=(rng.randrange(1,full) for _ in range(2))
        U=cache[a]|cache[full^a];V=cache[b]|cache[full^b]
        assert len({((x+u)%p,(y+v)%q) for x,y in U for u,v in V})==p*q
        count+=1
    print(f'Two-class obstruction: {full} lattice closures checked against HNF, {7*members} modular membership checks; {len(sizes)-1} bipartitions plus unpartitioned list; minimum union {min(sizes)}/{p*q}. True sumset {len(S)}/{p*q} with {n*n} ordered choices; {count} partition-pair sumsets full.')


def generator_tests():
    rng=Random(984211);total=0;peak=0
    for trial in range(120):
        d=1+trial%5;P=3+trial%7
        A=[tuple(rng.randrange(P) for _ in range(d)) for _ in range(2+trial%20)]
        A=list(dict.fromkeys(A));B=[A[0]];H=hnf(B)
        for a in A[1:]:
            if not hnf_member(H,[x-y for x,y in zip(a,A[0])]):
                B.append(a);H=hnf(B)
        assert H==hnf(A)
        assert len(B)-1<=d+(d*P).bit_length()*d
        peak=max(peak,len(B)-1);total+=len(A)
    print(f'Lattice generators: 120 exact HNF agreements, {total} input words; largest selected difference count {peak}.')


def ranks_and_moments():
    rng=Random(7581);tests=0;comp=0;maxfiber=0
    # Exhaustive equality-pattern moment bound in an independent tiny array model.
    for p in (3,5,7):
        for m,n,h in ((2,2,2),(2,2,3),(2,3,2)):
            mass=0;allarrays=0
            for arr in product(range(p),repeat=m*n):
                Z=sum(sum(arr[i*n+j] for i,j in enumerate(t))%p==0 for t in product(range(n),repeat=m))
                mass+=Z**h;allarrays+=1
            # Bell(2)=2, Bell(3)=5, use integer-cleared bound.
            B={2:2,3:5}[h];M=n**m
            assert mass*p**h<=allarrays*B**m*max(p,M)**h
            comp+=allarrays;tests+=1
    # Balanced-cut support/rank certificate; moduli product exceeds count mass.
    ps=(5,7,11,13);G=prod(ps);d=4
    for trial in range(24):
        m=3;n=3+trial%2
        A=[[tuple(rng.randrange(p) for p in ps) for _ in range(n)] for _ in range(m)]
        S=Counter({tuple(sum(a[j] for a in X)%ps[j] for j in range(d)):0 for X in product(*A)})
        for X in product(*A):S[tuple(sum(a[j] for a in X)%ps[j] for j in range(d))]+=1
        M=n**m;assert sum(S.values())==M<G
        C=[(0,)*d,tuple(rng.randrange(p) for p in ps)];Z=Counter()
        for x,c in S.items():
            for u in C:Z[tuple((v+w)%p for v,w,p in zip(x,u,ps))]+=c
        primes=(5,7,11);assert prod(primes)>sum(Z.values())
        for J in combinations(range(d),2):
            K=tuple(j for j in range(d) if j not in J)
            proj=lambda x,I:tuple(x[j] for j in I)
            F=Counter();E=Counter()
            for x,c in S.items():F[proj(x,J)]+=c;E[proj(x,K)]+=c
            B=max(max(F.values()),max(E.values()));maxfiber=max(maxfiber,B)
            rows={}
            for x,c in Z.items():rows.setdefault(proj(x,J),{})[proj(x,K)]=c
            cols=sorted({k for v in rows.values() for k in v})
            rr=[]
            for Q in primes:
                piv={}
                for v in rows.values():add_pivot([v.get(k,0) for k in cols],piv,Q)
                rr.append(len(piv))
            assert max(rr)*len(C)*B*B*len(primes)>=M
            tests+=1
    print(f'Moment/rank bounds: {comp} exhaustive tiny input arrays; {tests-144} moment inequalities; 144 exact balanced-cut rank checks; largest raw fiber {maxfiber}.')


def selector_test(path):
    spec=importlib.util.spec_from_file_location('joint',path);j=importlib.util.module_from_spec(spec);spec.loader.exec_module(j)
    d=10;ps=(17,)*d;qs=(3,5);m=3
    core=[(a,)*5+(b,)*5 for a,b in product(range(15),repeat=2)]
    rng=Random(67214);lists=[]
    for i in range(m):
        A=core+[tuple(int(k==j) for k in range(d)) for j in (i,i+5)]
        rng.shuffle(A);lists.append(A)
    C=[(0,)*d]
    assert all(len(j.shadow(lists,q,d)[2])>2 for q in qs)
    old=j.branch(lists,ps,qs,C)
    partitions=[peel(A,qs,2) for A in lists]
    assert all(set(P[0])==set(core) for P in partitions)
    branches=[j.branch(B,ps,qs,C) for B in product(*partitions)]
    z,phi,cap,peak,checks=j.choose(branches,ps)
    assert z is not None
    # Independent exact pair-sum membership check without enumerating 227^3 tuples.
    P={tuple((x+y)%p for x,y,p in zip(a,b,ps)) for a in lists[0] for b in lists[1]}
    assert all(tuple((x-y)%p for x,y,p in zip(z,c,ps)) not in P for c in lists[2])
    print(f'Automatic PR36-style rescue: list sizes {list(map(len,lists))}, product {prod(map(len,lists))}; whole-list ranks {[len(j.shadow(lists,q,d)[2]) for q in qs]}, recovered core ranks {[len(j.shadow([P[0] for P in partitions],q,d)[2]) for q in qs]}; {len(branches)} branches; potential {phi}; cap {cap}, peak {peak}, {checks} averaging checks; absent {z}.')


if __name__=='__main__':
    planted_tests();closure_tests();generator_tests();ranks_and_moments()
    if len(sys.argv)>1:selector_test(Path(sys.argv[1]))
