#!/usr/bin/env python3
"""Exact finite checks for shared-translation sampler. FFT checks numerical only."""
from itertools import product, combinations, permutations
from collections import Counter
from math import gcd, isqrt, sqrt
from random import Random
import numpy as np


def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def make_curve(p,a,c,r=None):
    assert (4*a**3+27*c*c)%p
    pts=[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-a*x-c)%p==0]
    N=len(pts)+1
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P; u,v=Q
        if x==u and (y+v)%p==0:return None
        s=((v-y)*pow(u-x,-1,p) if x!=u else (3*x*x+a)*pow(2*y,-1,p))%p
        z=(s*s-x-u)%p
        return z,(s*(x-z)-y)%p
    def mul(k,P):
        out=None
        while k:
            if k&1:out=add(out,P)
            P=add(P,P);k>>=1
        return out
    if r is None:r=max(d for d in range(2,N+1) if N%d==0 and prime(d))
    h=N//r
    assert r>2 and N==r*h and gcd(r,h)==1
    G=next(mul(h,P) for P in pts if mul(h,P) is not None)
    lookup={mul(j,G):j for j in range(r)}
    assert len(lookup)==r and mul(r,G) is None
    lifts=[[] for _ in range(p)]
    for P in pts:
        if P[1]!=0:lifts[P[0]].append(P)
    log={P:lookup[mul(h,P)] for P in pts};log[None]=0
    return r,h,lifts,log

def trim(a,p):
    a=[x%p for x in a]
    while a and not a[-1]:a.pop()
    return a

def rem(a,b,p):
    a=trim(a,p); b=trim(b,p); iv=pow(b[-1],-1,p)
    while len(a)>=len(b):
        j=len(a)-len(b);v=a[-1]*iv%p
        for i,x in enumerate(b):a[i+j]=(a[i+j]-v*x)%p
        a=trim(a,p)
    return a

def gcddeg(a,b,p):
    while b:a,b=b,rem(a,b,p)
    return len(a)-1

def cubic(a,c,s,p):return [(s**3+a*s+c)%p,(3*s*s+a)%p,3*s%p,1]

def good_offsets(p,a,c,S):
    return len(set(S))==len(S) and all(gcddeg(cubic(a,c,s,p),cubic(a,c,t,p),p)==0 for s,t in combinations(S,2))

def signed_types(k):
    def rec(rows,d):
        if len(rows)==k:yield tuple(rows);return
        yield from rec(rows+[(d,1)],d+1)
        for j in range(d):
            for e in [-1,1]:yield from rec(rows+[(j,e)],d)
    return list(rec([],0))

def classify(pts):
    blocks={};reps=[];rows=[]
    for P in pts:
        x,y=P
        if x not in blocks:blocks[x]=len(reps);reps.append(P)
        j=blocks[x]; rows.append((j,1 if P==reps[j] else -1))
    return tuple(rows)

def dim(t):return 1+max(j for j,e in t)

def sum_rows(vals,r,k):return tuple(sum(v[i] for v in vals)%r for i in range(k))

def block_table(p,a,c,r,lifts,log,block,k):
    assert good_offsets(p,a,c,[s for S in block for s in S])
    dims=(r,)*k; direct={}; independent={}
    for shift in range(p):
        Fs=[[P for s in S for P in lifts[(shift+s)%p]] for S in block]
        for rows in product(*(list(product(F,repeat=k)) for F in Fs)):
            ts=tuple(map(classify,rows));out=sum_rows([[log[P] for P in row] for row in rows],r,k)
            if ts not in direct:direct[ts]=np.zeros(dims,dtype=np.int64)
            direct[ts][out]+=1
    for ts in product(signed_types(k),repeat=len(block)):
        v=np.zeros(dims,dtype=np.int64)
        for offsets in product(*(list(permutations(S,dim(t))) for S,t in zip(block,ts))):
            for shift in range(p):
                Ls=[[lifts[(shift+s)%p] for s in seq] for seq in offsets]
                for reps in product(*(list(product(*L)) for L in Ls)):
                    out=sum_rows([[e*log[ps[j]]%r for j,e in t] for ps,t in zip(reps,ts)],r,k)
                    v[out]+=1
        if v.any():independent[ts]=v
    assert direct.keys()==independent.keys()
    for ts in direct:assert np.array_equal(direct[ts],independent[ts])
    return direct

def root_prime(r,start):
    q=start+((1-start)%r)
    while not prime(q):q+=r
    z=next(pow(g,(q-1)//r,q) for g in range(2,q) if pow(g,(q-1)//r,q)!=1)
    assert pow(z,r,q)==1
    return q,z

def ntt(A,r,q,z,inverse=False):
    if inverse:z=pow(z,-1,q)
    mat=np.array([[pow(z,i*j,q) for j in range(r)] for i in range(r)],dtype=np.int64)
    for ax in range(A.ndim):
        A=np.moveaxis(A,ax,0);A=(mat@A.reshape(r,-1)%q).reshape(A.shape);A=np.moveaxis(A,0,ax)
    if inverse:A=A*pow(pow(r,A.ndim,q),-1,q)%q
    return A

def admissible(types,k):return all(any(t[i]!=t[j] for t in types) for i in range(k) for j in range(i))

def test_moments(p,a,c,blocks,k):
    r,h,lifts,log=make_curve(p,a,c);dims=(r,)*k; qn=len(blocks)
    tabs=[block_table(p,a,c,r,lifts,log,B,k) for B in blocks]
    # Ground truth: for every actual translation vector, count raw ordered
    # tuples by subgroup target, including multiple different raw points with
    # the same projected logarithm. Then form falling factorials.
    masks=[]
    for block in blocks:
        hist=Counter()
        for shift in range(p):
            H=np.zeros(r,dtype=np.int64);H[0]=1
            for S in block:
                L=[log[P] for s in S for P in lifts[(shift+s)%p]]
                H=sum((np.roll(H,j) for j in L),np.zeros(r,dtype=np.int64))
            hist[tuple(int(v) for v in H)]+=1
        masks.append(hist)
    ix=np.indices(dims);prior=[sum((ix[j]==ix[i]) for j in range(i)) if i else 0 for i in range(k)]
    ground=np.zeros(dims,dtype=np.int64); cnt=0
    for choices in product(*(list(z.items()) for z in masks)):
        H=np.zeros(r,dtype=np.int64);H[0]=1;weight=1
        for H1,n in choices:
            H=sum((np.roll(H,j)*v for j,v in enumerate(H1) if v),np.zeros(r,dtype=np.int64));weight*=n
        term=np.ones(dims,dtype=np.int64)
        for i in range(k):term*=H[ix[i]]-prior[i]
        ground+=weight*term;cnt+=1
    # Independently convolve signed pattern tables; all comparisons exact.
    bound=1
    for tab in tabs:bound*=sum(int(v.sum()) for v in tab.values())
    prodq=1;npats=0
    for base in (100_000,200_000,300_000,400_000,500_000,600_000):
        q,z=root_prime(r,base);prodq*=q
        nt=[{t:ntt(v,r,q,z) for t,v in tab.items()} for tab in tabs]
        total=np.zeros(dims,dtype=np.int64);npats=0
        for ts in product(*(list(tab) for tab in tabs)):
            if not admissible([t for block in ts for t in block],k):continue
            v=np.ones(dims,dtype=np.int64)
            for i,t in enumerate(ts):v=v*nt[i][t]%q
            total=(total+v)%q;npats+=1
        assert np.array_equal(ntt(total,r,q,z,True),ground%q)
        if prodq>bound:break
    assert prodq>bound
    print(f'moments p={p}, group={r}*{h}, blocks={qn}, m={sum(map(len,blocks))}, k={k}: {cnt} seed classes, {sum(map(len,tabs))} block-pattern arrays, {npats} admissible combinations, {r**k} coefficients')
    return cnt,sum(map(len,tabs)),npats,r**k

def rref(rows,r,k):
    rows=[[int(x)%r for x in row] for row in rows];j=0
    for c in range(k):
        i=next((i for i in range(j,len(rows)) if rows[i][c]),None)
        if i is None:continue
        rows[i],rows[j]=rows[j],rows[i];v=pow(rows[j][c],-1,r);rows[j]=[x*v%r for x in rows[j]]
        for i in range(len(rows)):
            if i!=j:
                v=rows[i][c];rows[i]=[(x-v*y)%r for x,y in zip(rows[i],rows[j])]
        j+=1
        if j==len(rows):break
    return tuple(tuple(x) for x in rows[:j])

def tcols(t):return [tuple(e if j==a else 0 for j,e in t) for a in range(dim(t))]

def test_pruning():
    rng=Random(812);fixtures=0;checks=0
    for r,k in [(3,2),(5,3),(3,4),(7,3)]:
        ts=signed_types(k);fs=list(product(range(r),repeat=k))[1:];s=2*k+4
        for n in (8,17,30):
            for _ in range(20):
                B=[rref([v for i in range(rng.randint(1,3)) for v in tcols(rng.choice(ts))],r,k) for j in range(n)]
                alive=list(range(n));deleted=[];steps=0
                def active(f):return [j for j in alive if any(sum(x*y for x,y in zip(f,v))%r for v in B[j])]
                while True:
                    V=rref([v for j in alive for v in B[j]],r,k)
                    bad=None
                    for f in fs:
                        a=active(f)
                        if 0<len(a)<s:bad=a;break
                    if bad is None:break
                    deleted+=bad;alive=[j for j in alive if j not in bad];steps+=1
                    V2=rref([v for j in alive for v in B[j]],r,k)
                    assert len(V2)<len(V)
                assert steps<=k and len(deleted)<=k*(s-1)
                assert set(alive).isdisjoint(deleted) and sorted(alive+deleted)==list(range(n))
                for f in fs:
                    a=active(f);assert len(a)==0 or len(a)>=s;checks+=1
                fixtures+=1
    print(f'robust-span pruning: {fixtures} fixtures, {checks} functional checks')

def test_slope_avoidance():
    # D is the set of all nonzero base-field root differences, found here
    # independently from cubic gcds. Check t selection over the entire field.
    checked=0;bad=0;good=0;cases=0
    for p,a,c in [(13,4,1),(29,1,0),(43,0,2),(47,5,16)]:
        D={u for u in range(1,p) if gcddeg(cubic(a,c,0,p),cubic(a,c,u,p),p)>0}
        assert len(D)<=6
        for S in list(combinations(range(8),3)):
            diff={(x-y)%p for x in S for y in S if x!=y}
            excluded={v*pow(d,-1,p)%p for v in D for d in diff}
            assert len(excluded)<=6*len(diff)
            for t in range(1,p):
                actual=good_offsets(p,a,c,[t*x%p for x in S])
                assert actual==(t not in excluded)
                checked+=1;bad+=not actual;good+=actual
            cases+=1
    assert bad and good
    print(f'slope-avoidance identities: {cases} supports, {checked} slope checks, {good} allowed, {bad} excluded')

def test_geometry():
    combos=0;fibers=0;charchecks=0;avoid=0;bad=0;ma=0
    for p,a,c in [(13,4,1),(17,2,2),(19,1,0),(29,1,1),(37,3,1),(43,2,1),(43,0,2)]:
        try:r,h,lifts,log=make_curve(p,a,c)
        except (AssertionError,StopIteration):continue
        for d in (1,2,3):
            for offsets in list(combinations(range(min(p,7)),d))[:20]:
                if not good_offsets(p,a,c,offsets):bad+=1;continue
                avoid+=1
                arr=np.zeros((r,)*d,dtype=np.int64)
                for shift in range(p):
                    for ps in product(*(lifts[(shift+s)%p] for s in offsets)):
                        arr[tuple(log[P] for P in ps)]+=1
                # Geometric char-sum bound after deleting infinity and finite
                # branch fibers. Deliberately loose on these toy curves.
                omitted=(3*d+1)*2**d;g2=2**(d-1)*(3*d-3)
                assert abs(int(arr.sum())-p)<=(g2+2)*sqrt(p)+omitted+1
                f=np.fft.fftn(arr);f.flat[0]=0
                assert float(np.max(np.abs(f)))<=g2*sqrt(p)+omitted+1e-8
                charchecks+=r**d-1
                # Degree formula is a bound on rational fibers of [h]sum c_j P_j.
                for coeffs in product((-2,-1,0,1,2),repeat=d):
                    if not any(coeffs):continue
                    deg=h*h*2**(d-1)*sum(c*c for c in coeffs)
                    counts=[0]*r
                    for ind in np.argwhere(arr):
                        counts[sum(int(i)*c for i,c in zip(ind,coeffs))%r]+=int(arr[tuple(ind)])
                    # For a nonconstant map to E, every geometric fiber has at most degree many points.
                    assert max(counts)<=deg
                    combos+=1;fibers+=r;ma=max(ma,max(counts))
        print(f'geometry p={p}, r={r}, h={h}: cumulative good offset sets {avoid}, bad {bad}')
    print(f'geometry totals: {avoid} admissible offset sets, {bad} rejected sets, {combos} coefficient choices, {fibers} exact projected-fiber bounds, {charchecks} numerical character comparisons; max observed fiber {ma}')

if __name__=='__main__':
    assert __debug__
    # Choose small fixtures with disjoint branch sets within each shared block.
    tests=[(5,2,1,[[(0,1),(2,3)]],3),
           (13,4,1,[[(0,1),(2,3)],[(0,1)]],2),
           (13,4,1,[[(0,),(1,)],[(0,)],[(2,)]],3),
           (19,1,0,[[(0,),(1,)],[(0,1)]],2),
           (5,2,1,[[(0,),(1,)],[(2,)],[(1,)],[(0,)],[(1,)]],2),
           (43,0,2,[[(0,),(1,2)],[(0,1)]],2)]
    sums=[0]*4
    for p,a,c,bs,k in tests:
        # Find a common slope that gives distinct branch loci in each block.
        for t in range(1,p):
            B=[[tuple(t*s%p for s in S) for S in block] for block in bs]
            if all(good_offsets(p,a,c,[s for S in block for s in S]) for block in B):break
        else:raise AssertionError('no good slope')
        z=test_moments(p,a,c,B,k);sums=[a+b for a,b in zip(sums,z)]
    print('exact moment totals (seed classes, block-pattern arrays, admissible combinations, coefficients):',sums)
    test_pruning();test_slope_avoidance();test_geometry()
    print('PASS: exact component identities and numerical character checks; not an ECDLP runtime test.')
