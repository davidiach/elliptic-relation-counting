#!/usr/bin/env python3
"""Finite verification for replica_mixing_note.md. Exact integer/NTT checks;
complex FFT character regressions are explicitly numerical, not proofs."""
from itertools import product, combinations, permutations
from collections import Counter
from math import comb, isqrt, sqrt
from random import Random
import numpy as np


def curve(p,a,c):
    assert (4*a**3+27*c*c)%p
    pts=[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-a*x-c)%p==0]
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P; u,v=Q
        if x==u and (y+v)%p==0:return None
        s=((v-y)*pow(u-x,-1,p) if x!=u else (3*x*x+a)*pow(2*y,-1,p))%p
        z=(s*s-x-u)%p
        return z,(s*(x-z)-y)%p
    seq=[None]; P=pts[0]
    while P is not None:
        seq.append(P);P=add(P,pts[0])
        assert len(seq)<=len(pts)+1
    r=len(seq);assert r==len(pts)+1 and prime(r)
    lifts=[[] for _ in range(p)]
    for j,P in enumerate(seq[1:],1):lifts[P[0]].append(j)
    return r,lifts

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def masks(p,lifts,S):
    return Counter(tuple(sorted(P for s in S for P in lifts[(a+t*s)%p]))
                   for a in range(p) for t in range(1,p))

def signed_types(k):
    """Rows encoded (block index, sign), first occurrence has sign +1."""
    def rec(rows,d):
        if len(rows)==k:
            yield tuple(rows);return
        yield from rec(rows+[(d,1)],d+1)
        for j in range(d):
            for e in [-1,1]:yield from rec(rows+[(j,e)],d)
    return list(rec([],0))

def rep_tuple(typ, vals, r):
    return tuple(e*vals[j]%r for j,e in typ)

def d_of(typ):return 1+max(j for j,e in typ)

def classify(inds,r,x_of):
    blocks={}; reps=[]; rows=[]
    for v in inds:
        x=x_of[v]
        if x not in blocks:
            blocks[x]=len(reps);reps.append(v)
        j=blocks[x]; s=1 if reps[j]==v else -1
        assert reps[j]*s%r==v
        rows.append((j,s))
    return tuple(rows),tuple(reps)

def type_tables(p,r,lifts,S,k):
    ts=signed_types(k); dims=(r,)*k
    by_type={t:np.zeros(dims,dtype=np.int64) for t in ts}
    independent={t:np.zeros(dims,dtype=np.int64) for t in ts}
    x_of={P:x for x,L in enumerate(lifts) for P in L}
    masks0=masks(p,lifts,S)
    for F,mul in masks0.items():
        for row in product(F,repeat=k):
            t,_=classify(row,r,x_of);by_type[t][row]+=mul
        # Independent construction via retained representatives of each signed pattern.
        for t in ts:
            d=d_of(t)
            for vals in product(F,repeat=d):
                if len({x_of[v] for v in vals})==d:
                    independent[t][rep_tuple(t,vals,r)]+=mul
    for t in ts: assert np.array_equal(by_type[t],independent[t])
    return {t:v for t,v in by_type.items() if v.any()},masks0,len(ts)

def root_prime(r,start):
    q=start+((1-start)%r)
    while not prime(q):q+=r
    z=next(pow(g,(q-1)//r,q) for g in range(2,q) if pow(g,(q-1)//r,q)!=1)
    assert pow(z,r,q)==1
    return q,z

def ntt(a,r,q,z,inverse=False):
    if inverse:z=pow(z,-1,q)
    mat=np.array([[pow(z,i*j,q) for j in range(r)] for i in range(r)],dtype=np.int64)
    for ax in range(a.ndim):
        a=np.moveaxis(a,ax,0)
        a=(mat@a.reshape(r,-1)%q).reshape(a.shape)
        a=np.moveaxis(a,0,ax)
    if inverse:a=a*pow(pow(r,a.ndim,q),-1,q)%q
    return a

def admissible(types,k):
    return all(any(t[i]!=t[j] for t in types) for i in range(k) for j in range(i))

def exact_moment_case(p,a,c,supports,k):
    r,lifts=curve(p,a,c); m=len(supports); dims=(r,)*k
    tabs=[]; mm=[]; local=0
    for S in supports:
        tab,mask,n=type_tables(p,r,lifts,S,k);tabs.append(tab);mm.append(mask);local+=n
    # Ground truth from seed combinations, ordinary integer group convolution,
    # and falling-factorial products for repeated targets.
    ix=np.indices(dims); prior=[sum((ix[j]==ix[i]) for j in range(i)) if i else 0 for i in range(k)]
    ground=np.zeros(dims,dtype=np.int64); seed_classes=0
    for chosen in product(*[list(z.items()) for z in mm]):
        h=np.zeros(r,dtype=np.int64);h[0]=1;mul=1
        for F,n in chosen:
            h=sum((np.roll(h,j) for j in F),np.zeros(r,dtype=np.int64));mul*=n
        term=np.ones(dims,dtype=np.int64)
        for i in range(k):term*=h[ix[i]]-prior[i]
        ground+=mul*term;seed_classes+=1
    # Independent signed-pattern convolution, exactly in two prime fields.
    qs=[]; ns=[]; pattern_count=0
    for base in (1_000_000,2_000_000,3_000_000):
        q,z=root_prime(r,base);qs.append(q)
        nt=[{t:ntt(v,r,q,z) for t,v in tab.items()} for tab in tabs]
        total=np.zeros(dims,dtype=np.int64);cnt=0
        for ts in product(*[list(t) for t in tabs]):
            if not admissible(ts,k):continue
            v=np.ones(dims,dtype=np.int64)
            for i,t in enumerate(ts):v=v*nt[i][t]%q
            total=(total+v)%q;cnt+=1
        result=ntt(total,r,q,z,inverse=True)
        assert np.array_equal(result,ground%q)
        pattern_count=cnt
        if len(qs)>=2 and np.prod(qs,dtype=object)>int(ground.max()):break
    # A priori count bound ensures modular equality certifies integer equality.
    bound=1
    for tab in tabs:bound*=sum(int(v.sum()) for v in tab.values())
    # Here ground.max() and each pattern-sum cell are both <= bound. If needed,
    # independently verify more primes; all selected fixtures fit the product.
    if np.prod(qs,dtype=object)<=bound:
        for base in (4_000_000,5_000_000,6_000_000):
            q,z=root_prime(r,base);qs.append(q)
            nt=[{t:ntt(v,r,q,z) for t,v in tab.items()} for tab in tabs]
            total=np.zeros(dims,dtype=np.int64)
            for ts in product(*[list(t) for t in tabs]):
                if not admissible(ts,k):continue
                v=np.ones(dims,dtype=np.int64)
                for i,t in enumerate(ts):v=v*nt[i][t]%q
                total=(total+v)%q
            assert np.array_equal(ntt(total,r,q,z,True),ground%q)
            if np.prod(qs,dtype=object)>bound:break
    assert np.prod(qs,dtype=object)>bound
    print(f'exact moment p={p} r={r} m={m} k={k}: {seed_classes:,} seed classes, {local} local signed types, {pattern_count:,} admissible patterns, {r**k:,} target vectors, {len(qs)} primes')
    return local,pattern_count,r**k,seed_classes

def trim(a,p):
    a=[x%p for x in a]
    while a and a[-1]==0:a.pop()
    return a

def rem(a,b,p):
    a=trim(a,p);b=trim(b,p);iv=pow(b[-1],-1,p)
    while len(a)>=len(b):
        k=len(a)-len(b);v=a[-1]*iv%p
        for j in range(len(b)):a[j+k]=(a[j+k]-v*b[j])%p
        a=trim(a,p)
    return a

def gcddeg(a,b,p):
    while b:a,b=b,rem(a,b,p)
    return len(a)-1

def cubic(a,c,h,p):return [(h**3+a*h+c)%p,(3*h*h+a)%p,3*h%p,1]

def character_case(p,a,c,S):
    r,lifts=curve(p,a,c);d=len(S);bad=good=0;maxrat=0.;num=0;mass=0
    g2=(3*d-3)*2**(d-1) # 2g-2
    for t in range(1,p):
        fs=[cubic(a,c,t*s%p,p) for s in S]
        ok=all(gcddeg(f,h,p)==0 for f,h in combinations(fs,2))
        arr=np.zeros((r,)*d,dtype=np.int64)
        for v in range(p):
            for row in product(*[lifts[(v+t*s)%p] for s in S]):arr[row]+=1
        mass+=int(arr.sum())
        if not ok:bad+=1;continue
        good+=1
        assert abs(int(arr.sum())-p)<= (g2+2)*sqrt(p)+2**d+1
        f=np.fft.fftn(arr);f.flat[0]=0
        upper=g2*sqrt(p)+2**d
        rat=float(np.max(np.abs(f)))/upper
        assert rat<=1+1e-8
        maxrat=max(maxrat,rat);num+=r**d-1
    assert bad<=9*comb(d,2)
    print(f'character p={p} r={r} d={d}: {good} good slopes, {bad} bad slopes, {num:,} nontrivial character comparisons (numerical), max ratio to bound {maxrat:.6f}, total mass {mass}')
    return good,bad,num

def rref_key(rows,p):
    rows=[[x%p for x in row] for row in rows];j=0;n=len(rows[0])
    for i in range(n):
        k=next((k for k in range(j,len(rows)) if rows[k][i]),None)
        if k is None:continue
        rows[j],rows[k]=rows[k],rows[j];v=pow(rows[j][i],-1,p)
        rows[j]=[x*v%p for x in rows[j]]
        for k in range(len(rows)):
            if k!=j:
                v=rows[k][i];rows[k]=[(x-v*y)%p for x,y in zip(rows[k],rows[j])]
        j+=1
        if j==len(rows):break
    return tuple(tuple(v) for v in rows[:j])

def pattern_check():
    total=0;compat=0
    rng=Random(17)
    for k in range(2,7):
        keys=set();ts=signed_types(k)
        for t in ts:
            d=d_of(t); rows=[[s if j==h else 0 for j,s in t] for h in range(d)]
            key=rref_key(rows,7);assert len(key)==d and key not in keys;keys.add(key)
            total+=1
            if d<k:
                for _ in range(10):
                    vs=[rng.randrange(1,7) for h in range(d)];T=rep_tuple(t,vs,7)
                    if any(T[i]==-T[j]%7 for i in range(k) for j in range(i)):continue
                    assert any(t[i]==t[j] for i in range(k) for j in range(i))
                    compat+=1
        print(f'signed partitions k={k}: {len(ts):,} distinct subspaces')
    print(f'pattern checks: {total:,} exact rank/uniqueness; {compat:,} compatible proper-type duplicate checks')

if __name__=='__main__':
    assert __debug__
    pattern_check()
    local=patterns=targets=classes=0
    for args in [(5,2,1,[(0,1),(0,2),(0,1)],2),
                 (5,2,1,[(0,1),(0,2),(0,1)],3),
                 (5,2,1,[(0,1,3),(0,2),(0,1),(0,2,4)],3),
                 (7,0,3,[(0,1),(0,2,4),(1,3)],3)]:
        t=exact_moment_case(*args);local+=t[0];patterns+=t[1];targets+=t[2];classes+=t[3]
    print(f'exact moment totals: {local} local signed-type arrays, {patterns:,} admissible pattern tuples, {targets:,} target-vector coefficients, {classes:,} seed classes')
    c=[character_case(*x) for x in [(5,2,1,[0]),(13,0,2,[0,1]),(13,0,2,[0,1,4]),(13,0,2,[0,1,3,8]),(37,0,5,[0,1,7]),(101,1,32,[0,9])]]
    print('character totals:',tuple(sum(t[i] for t in c) for i in range(3)))
    print('PASS (finite component checks only; no exponent improvement claimed)')
