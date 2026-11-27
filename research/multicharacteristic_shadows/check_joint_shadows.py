"""Exact synthetic tests for simultaneous affine shadows. Python 3.10+, stdlib;
NumPy is used only for the bounded small unfolding-rank fixture.
"""
from itertools import product, combinations
from collections import Counter
from fractions import Fraction
from math import prod, lcm
from random import Random


def ceildiv(a, b):
    return -(-a // b)


def rref(rows, q, d):
    pivots = {}
    for row in rows:
        row = [v % q for v in row]
        for j in sorted(pivots):
            if row[j]:
                c = row[j]
                row = [(v-c*w) % q for v,w in zip(row,pivots[j])]
        j = next((j for j,v in enumerate(row) if v), None)
        if j is None:
            continue
        t = pow(row[j], -1, q)
        row = [v*t % q for v in row]
        for h in list(pivots):
            if pivots[h][j]:
                c = pivots[h][j]
                pivots[h] = [(v-c*w) % q for v,w in zip(pivots[h],row)]
        pivots[j] = row
    I = sorted(pivots)
    return tuple(I), tuple(tuple(pivots[j]) for j in I)


def shadow(lists, q, d):
    ref = [min(A) for A in lists]
    off = tuple(sum(x[j] for x in ref)%q for j in range(d))
    I, R = rref((tuple(a[j]-v[j] for j in range(d))
                 for A,v in zip(lists,ref) for a in A),q,d)
    return q, off, I, R


def branch(lists, ps, qs, corrections):
    d=len(ps)
    lows = tuple(sum(min(a[j] for a in A) for A in lists) for j in range(d))
    highs= tuple(sum(max(a[j] for a in A) for A in lists) for j in range(d))
    shadows=tuple(shadow(lists,q,d) for q in qs)
    weights=[]; maxima=[]
    for j,p in enumerate(ps):
        w=highs[j]-lows[j]+1
        P=prod(q for q,off,I,R in shadows if j in I)
        D=prod(q for q,off,I,R in shadows if j not in I)
        weights.append(min(Fraction(P), Fraction(ceildiv(w,D),p)))
        maxima.append(min(P,ceildiv(w,lcm(p,D))))
    suffix=[Fraction(1)]*(d+1)
    for j in reversed(range(d)):
        suffix[j]=suffix[j+1]*weights[j]
    states=[(c, tuple(() for _ in qs)) for c in corrections]
    return {'lo':lows,'hi':highs,'sh':shadows,'a':weights,'B':maxima,
            'suf':suffix,'states':states,'L':len(corrections)}


def children(br, j, p, z, states):
    result=[]
    lo,hi=br['lo'][j],br['hi'][j]
    for c,coeffs in states:
        first=lo+((z-c[j]-lo)%p)
        options=set()
        for w in range(first,hi+1,p):
            nxt=[]
            for (q,off,I,R),coeff in zip(br['sh'],coeffs):
                if j in I:
                    assert len(coeff)==I.index(j)
                    nxt.append(coeff+((w-off[j])%q,))
                else:
                    assert len(coeff)==sum(i<j for i in I)
                    pred=(off[j]+sum(a*row[j] for a,row in zip(coeff,R)))%q
                    if pred!=w%q:
                        break
                    nxt.append(coeff)
            else:
                options.add(tuple(nxt))
        assert len(options)<=br['B'][j]
        result.extend((c,s) for s in sorted(options))
    return result


def choose(branches, ps, cap=None):
    live = [br['states'] for br in branches]
    start=sum(len(states)*br['suf'][0] for br,states in zip(branches,live))
    bound=sum(br['L']*prod(br['B']) for br in branches)
    if start>=1 or (cap is not None and bound>cap):
        return None,start,bound,0,0
    prefix=[];peak=sum(map(len,live)); checks=0
    potential=start
    for j,p in enumerate(ps):
        best=None
        allcounts=[0]*len(branches)
        mass=0
        for z in range(p):
            ns=[children(br,j,p,z,states) for br,states in zip(branches,live)]
            val=sum(len(st)*br['suf'][j+1] for st,br in zip(ns,branches))
            for i,st in enumerate(ns):allcounts[i]+=len(st)
            mass+=val
            if best is None or val<best[0]:best=(val,z,ns)
        for br,count,states in zip(branches,allcounts,live):
            assert Fraction(count,p)<=len(states)*br['a'][j]
            checks+=1
        assert mass/p<=potential
        potential,z,live=best
        assert potential<1
        prefix.append(z)
        peak=max(peak,sum(map(len,live)))
        assert peak<=bound
    assert not any(live)
    return tuple(prefix),start,bound,peak,checks


def check_selector_reuse_and_cap():
    ps = (3,)
    lists = [[(0,)]]
    corrections = [(0,), (1,)]
    br = branch(lists, ps, (2,), corrections)
    actual = set(convolve(lists, ps, corrections))
    assert choose([br], ps, cap=1) == (None, Fraction(2, 3), 2, 0, 0), 'cap 1 must reject before searching'
    first = choose([br], ps, cap=2)
    assert first[0] == (2,) and first[0] not in actual
    second = choose([br], ps, cap=2)
    assert second == first, 'reusing a branch changed the selector result'
    assert {z for z in product(range(3)) if feasible(br, ps, z)} == actual
    print(f'Bound {first[2]} at potential {first[1]} rejected cap 1. Two searches at cap 2 returned {first[0]} and preserved branch membership.')


def brute_states(br, ps, z):
    """Independent direct Cartesian lift test, only for tiny inputs."""
    cnt=0
    for c,_ in br['states']:
        ranges=[range(lo+((x-y-lo)%p),hi+1,p)
                for x,y,lo,hi,p in zip(z,c,br['lo'],br['hi'],ps)]
        for w in product(*ranges):
            if all(all((w[j]-off[j]-sum((w[i]-off[i])*row[j]
                         for i,row in zip(I,R)))%q==0 for j in range(len(ps)))
                   for q,off,I,R in br['sh']):
                cnt+=1
    return cnt


def feasible(br, ps, z):
    states=br['states']
    for j,(p,x) in enumerate(zip(ps,z)):
        states=children(br,j,p,x,states)
    return bool(states)


def convolve(lists, ps, corrections=None):
    S=Counter({(0,)*len(ps):1})
    if corrections is not None:lists=list(lists)+[corrections]
    for A in lists:
        T=Counter()
        for x,c in S.items():
            for y in A:
                T[tuple((a+b)%p for a,b,p in zip(x,y,ps))]+=c
        S=T
    assert sum(S.values())==prod(map(len,lists))
    return S


def run_small():
    rng=Random(140571)
    comparisons=0;inclusions=0;success=0;failed=0;avg=0;different=0
    for t in range(100):
        ps=(3,5,3) if t%2 else (5,3,4)
        qs=(2,3) if t%3 else (2,5)
        d=len(ps)
        points=list(product(*map(range,ps)))
        lists=[sorted(set(rng.sample(points,rng.randrange(1,4)))) for _ in range(2)]
        C=rng.sample(points,rng.randrange(1,3))
        br=branch(lists,ps,qs,C)
        different+=br['sh'][0][2]!=br['sh'][1][2]
        S=convolve(lists,ps,C)
        for z in points:
            assert feasible(br,ps,z)==bool(brute_states(br,ps,z))
            comparisons+=1
            if z in S:
                assert feasible(br,ps,z)
                inclusions+=1
        z,g,B,peak,a=choose([branch(lists,ps,qs,C)],ps)
        avg+=a
        if z is None:failed+=1
        else:
            success+=1;assert z not in S
    print(f'Random tiny mixed groups: {comparisons} exact direct-lift comparisons, {inclusions} represented-point inclusions; {different} different pivot-set cases; {success} certified absences, {failed} honest guard failures; {avg} exact averaging checks.')


def run_union():
    # Each class consists of a single point. All nine branches can be
    # avoided simultaneously. Different intervals/pivot sets tested below.
    ps=(7,11,13);qs=(2,3,5); d=len(ps)
    A=[[(0,0,0)],[(2,4,6)],[(6,10,12)]]
    B=[[(0,0,0)],[(1,2,3)],[(4,3,2)]]
    C=[(0,0,0),(2,1,1)]
    br=[branch([a,b],ps,qs,C) for a,b in product(A,B)]
    z,g,cap,peak,checks=choose(br,ps)
    S=convolve([sum(A,[]),sum(B,[])],ps,C)
    assert z is not None and z not in S
    print(f'Affine-cover union: {len(br)} branches, {len(C)} corrections, initial potential {g}, cap {cap}, peak {peak}, absent {z}; {checks} averaging checks.')
    # Different pivot sets and nonzero intervals, verified against direct sums.
    rng=Random(111)
    totals=[0,0,0]
    ps=(7,7,7,7,7);qs=(2,3)
    for t in range(30):
        e=lambda j,v: tuple(v if k==j else 0 for k in range(5))
        cover=[ [[(0,)*5,e(t%5,1)],[(2,)*5,tuple(2+v for v in e((t+1)%5,1))]], [[(1,)*5]] ]
        lists=[sorted(set(sum(classes,[]))) for classes in cover]
        C=[(0,)*5, tuple(rng.randrange(7) for _ in range(5))]
        br=[branch([a,b],ps,qs,C) for a,b in product(*cover)]
        z,g,cap,peak,checks=choose(br,ps)
        assert z is not None and z not in convolve(lists,ps,C)
        totals[0]+=1; totals[1]+=checks; totals[2]+=len({tuple(s[2] for s in x['sh']) for x in br})>1
    print(f'Varying affine-cover fixtures: {totals[0]} successes, {totals[1]} exact averaging checks, {totals[2]} fixtures with different branch pivot patterns.')


def run_repetition():
    ps=(7,)*5; qs=(2,3); m=4;d=len(ps)
    A=[(0,)*d,(5,)*d]+[tuple(6 if k==j else 0 for k in range(d)) for j in range(d)]
    lists=[A]*m;C=[(0,)*d]
    S=convolve(lists,ps)
    assert prod(map(len,lists))<prod(ps)
    br=branch(lists,ps,qs,C)
    # Each prime individually is compatible with every target: U>=q*p-1,
    # gcd(q,p)=1. Every other auxiliary prime has full affine rank.
    for q,off,I,R in br['sh']:
        assert len(I)==1 and all(hi>=q*p-1 for hi,p in zip(br['hi'],ps))
        single=branch(lists,ps,(q,),C)
        assert all(feasible(single,ps,z) for z in product(*map(range,ps)))
    for q in [5,7,11,13,17,19,23]:
        assert len(shadow(lists,q,d)[2])==d
    assert all(len({x[j] for x in S})==7 for j in range(d))
    z,g,cap,peak,c=choose([br],ps)
    assert z is not None and z not in S
    print(f'Joint-only repetition example: ambient {prod(ps)}, product {7**4}, support {len(S)}, single-shadow support {prod(ps)} for both primes; all singleton projections full. Joint initial potential {g}, candidate cap {cap}, actual peak {peak}, absent {z}.')
    # Adding standard units destroys every global affine-shadow obstruction.
    A2=sorted(set(A+[tuple(1 if k==j else 0 for k in range(d)) for j in range(d)]))
    lists2=[A2]+[A]*(m-1)
    assert prod(map(len,lists2))<prod(ps)
    for q in [2,3,5,7,11]:assert len(shadow(lists2,q,d)[2])==d
    br2=branch(lists2,ps,(2,3,5,7,11),C)
    assert choose([br2],ps)[0] is None
    # Every target has common lifts w=z; this gives joint compatibility.
    for z in product(*map(range,ps)):
        assert all(lo<=a<=hi for a,lo,hi in zip(z,br2['lo'],br2['hi']))
    print(f'All-field full-rank negative: product {prod(map(len,lists2))} < {prod(ps)}, joint affine relaxation still covers the ambient space.')


def crt(residues, primes):
    P=prod(primes)
    return sum(x*(P//p)*pow(P//p,-1,p) for x,p in zip(residues,primes))%P


def unfolding_rank(S,ps,J,q):
    import numpy as np
    J=tuple(J);K=tuple(i for i in range(len(ps)) if i not in J)
    def idx(x,I):
        v=0
        for j in I:v=v*ps[j]+x[j]
        return v
    M=np.zeros((prod(ps[j] for j in J),prod(ps[j] for j in K)),dtype=np.int64)
    for x,c in S.items():M[idx(x,J),idx(x,K)]=c%q
    # Only rank, eliminate rows below pivot.
    r=0
    for j in range(M.shape[1]):
        ix=np.flatnonzero(M[r:,j])
        if not len(ix):continue
        k=r+int(ix[0]);M[[r,k],:]=M[[k,r],:]
        M[r,j:]=M[r,j:]*pow(int(M[r,j]),-1,q)%q
        c=M[r+1:,j].copy()
        nz=np.flatnonzero(c)+r+1
        M[nz,j:]=(M[nz,j:]-c[nz-r-1,None]*M[r,j:])%q
        r+=1
        if r==M.shape[0]:break
    return r


def run_code_fixture():
    # Small analogue of the product-code construction. The Q0 coefficient
    # groups can be empty in the toy, unlike the asymptotic proof.
    ps=(31,)*4;d=len(ps);qs=(2,3);q0=5;L=6
    lists=[]
    for i in range(4):
        core=[]
        for a,b,c in product(range(q0) if i<2 else range(1),range(2),range(3)):
            core.append(tuple(crt((b,c,(a*pow(j+1,i,q0))%q0),(2,3,q0)) for j in range(d)))
        assert len(set(core))==len(core)
        A=sorted(set(core+[tuple(L if k==j else 0 for k in range(d)) for j in range(d)]))
        lists.append(A)
    assert prod(map(len,lists))<prod(ps)
    C=[(0,)*d]
    br=branch(lists,ps,qs,C)
    for q,off,I,R in br['sh']:
        assert len(I)==1 and all(hi>=q*p-1 for hi,p in zip(br['hi'],ps))
    for q in [5,7,11,31]:assert len(shadow(lists,q,d)[2])==d
    z,g,cap,peak,c=choose([br],ps)
    S=convolve(lists,ps)
    assert z is not None and z not in S
    assert all(len({x[j] for x in S})==31 for j in range(d))
    ranks=[unfolding_rank(S,ps,J,101) for J in [(0,1),(0,2),(0,3)]]
    print(f'Joint-only CRT-code fixture: list sizes {list(map(len,lists))}, product {prod(map(len,lists))}, support {len(S)}, ambient {prod(ps)}; all singleton projections full. All three balanced unfolding ranks over F101: {ranks}. Joint potential {g}, cap {cap}, peak {peak}, absent {z}.')


if __name__=='__main__':
    run_small();run_union();run_repetition();run_code_fixture()
    check_selector_reuse_and_cap()
