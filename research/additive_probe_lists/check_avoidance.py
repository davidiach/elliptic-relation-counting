"""Exact finite checks of additive range avoidance and probe constructions.
All inputs are explicit synthetic sets in finite planes; no curve or key inputs.
"""
from bisect import bisect_left
from collections import Counter
from itertools import combinations, product
from math import prod
from random import Random
from fractions import Fraction


def add(a,b,p):
    return ((a[0]+b[0])%p,(a[1]+b[1])%p)


def is_plane_point(point, p):
    return len(point) == 2 and all(type(c) is int and 0 <= c < p for c in point)


def direct_tuple_sums(lists, p):
    for row in product(*lists):
        yield (sum(x for x, y in row) % p, sum(y for x, y in row) % p)


def literal_sumset(lists,p):
    s={(0,0)}
    for L in lists:
        s={add(a,b,p) for a in s for b in L}
    return s


def cyclic_interval(a,b,shift,p):
    """Half-open [a,b) shifted by -shift modulo p, disjoint intervals."""
    w=b-a
    if not w:return []
    if w==p:return [(0,p)]
    s=(a-shift)%p
    return [(s,s+w)] if s+w<=p else [(s,p),(0,s+w-p)]


def rectangle_count(A,B,R,p):
    """Offline range counts, preserving multiplicity. Word-RAM finite verifier."""
    lo,hi,bot,top=R
    assert 0<=lo<=hi<=p and 0<=bot<=top<=p
    events=[]
    for x,y in A:
        for l,h in cyclic_interval(lo,hi,x,p):
            for s,t in cyclic_interval(bot,top,y,p):
                events.extend(((h,t,1),(l,t,-1),(h,s,-1),(l,s,1)))
    B=sorted(B);ys=sorted({y for _,y in B});F=[0]*(len(ys)+1)
    def put(i):
        i+=1
        while i<len(F):F[i]+=1;i+=i&-i
    def pref(i):
        total=0
        while i:total+=F[i];i-=i&-i
        return total
    j=0;ans=0
    for x,y,s in sorted(events):
        while j<len(B) and B[j][0]<x:
            put(bisect_left(ys,B[j][1]));j+=1
        ans+=s*pref(bisect_left(ys,y))
    return ans


def partial_sums(lists,p):
    A=[(0,0)]
    for L in lists:A=[add(a,b,p) for a in A for b in L]
    return A


def avoid(lists,p):
    M=prod(map(len,lists));assert 0<M<p*p
    L=[[],[]];P=[1,1]
    for s in sorted(lists,key=len,reverse=True):
        j=0 if P[0]<=P[1] else 1
        L[j].append(s);P[j]*=len(s)
    assert max(P)<=min(P)*max(map(len,lists))
    A,B=(partial_sums(a,p) for a in L)
    R=(0,p,0,p);c=M;queries=0
    def area(R):return (R[1]-R[0])*(R[3]-R[2])
    while area(R)>1:
        x,X,y,Y=R
        if X-x>1:
            mid=(x+X)//2;R0=(x,mid,y,Y);R1=(mid,X,y,Y)
        else:
            mid=(y+Y)//2;R0=(x,X,y,mid);R1=(x,X,mid,Y)
        c0=rectangle_count(A,B,R0,p);queries+=1
        assert 0<=c0<=c
        if c0<area(R0):R,c=R0,c0
        else:R,c=R1,c-c0
        assert c<area(R)
    assert c==0
    return (R[0],R[2]),queries,max(P)


def range_tests():
    rng=Random(87342);tests=0;cases=0;maxhalf=0;sumqueries=0
    for p in [3,5,7,11,17,31,101]:
        G=list(product(range(p),repeat=2))
        for _ in range(24):
            m=rng.randrange(1,8);sizes=[];M=1
            for i in range(m):
                s=rng.randint(1,min(6,(p*p-1)//M));sizes.append(s);M*=s
            # Repetitions are intentional: tuple multiplicities matter.
            L=[[rng.choice(G) for j in range(s)] for s in sizes]
            j=rng.randrange(m+1);A=partial_sums(L[:j],p);B=partial_sums(L[j:],p)
            C = Counter(direct_tuple_sums(L, p))
            for k in range(7):
                x,X=sorted((rng.randrange(p+1),rng.randrange(p+1)))
                y,Y=sorted((rng.randrange(p+1),rng.randrange(p+1)))
                c=sum(v for (a,b),v in C.items() if x<=a<X and y<=b<Y)
                assert rectangle_count(A,B,(x,X,y,Y),p)==c
                tests+=1
            z,q,h=avoid(L,p)
            assert is_plane_point(z, p), f'Noncanonical avoidance point for p={p}: {z}'
            assert z not in C, f'Returned point occurs in the sumset: {z}'
            cases+=1;sumqueries+=q;maxhalf=max(maxhalf,h)
    print(f'Exact modular rectangle counts: {tests}; complete avoidance instances: {cases}; bisection queries: {sumqueries}; largest enumerated half: {maxhalf}.')


def normalized_cover(E,sizes,p):
    """Every E of size <=1+sum(n_i-1) is coverable by legal lists."""
    G=list(product(range(p),repeat=2));E=list(dict.fromkeys(E))
    assert 0<len(E)<=1+sum(s-1 for s in sizes)
    z=E[0];D=[((x-z[0])%p,(y-z[1])%p) for x,y in E[1:]]
    A=[];j=0
    for s in sizes:
        L={(0,0),*D[j:j+s-1]};j+=s-1
        for v in G:
            if len(L)==s:break
            L.add(v)
        assert len(L)==s
        A.append(L)
    A[0]={add(z,a,p) for a in A[0]}
    assert set(E)<=literal_sumset(A,p)
    return A


def universe_minimum(p,sizes):
    """Exhaust exact universal probe minimum for tiny finite planes."""
    G=list(product(range(p),repeat=2));idx={a:i for i,a in enumerate(G)}
    sets=[[((0,0),*s) for s in combinations(G[1:],n-1)] for n in sizes]
    support=set()
    for L in product(*sets):
        S=literal_sumset(L,p)
        for t in G:support.add(sum(1<<idx[add(a,t,p)] for a in S))
    k0=1+sum(n-1 for n in sizes)
    assert prod(sizes)<p*p
    # Lower bound is also checked by directly verifying containment for every size-k0 probe set.
    for I in combinations(range(p*p),k0):
        mask=sum(1<<j for j in I)
        assert any(mask&s==mask for s in support)
    for k in range(k0+1,p*p+1):
        for I in combinations(range(p*p),k):
            mask=sum(1<<j for j in I)
            if not any(mask&s==mask for s in support):
                return k,len(support),[G[j] for j in I]
    raise AssertionError('whole plane must avoid')


def universal_tests():
    rng=Random(1729);nchecks=0
    for p,ss in [(3,[2,2]),(3,[2,3]),(3,[2,2,2]),(5,[2,2]),(7,[2,3,4])]:
        G=list(product(range(p),repeat=2));D=1+sum(n-1 for n in ss)
        for _ in range(25):
            probes = rng.sample(G, rng.randrange(1,D+1))
            blocks = normalized_cover(probes, ss, p)
            assert [len(block) for block in blocks] == ss, f'Incorrect returned block sizes for caps {ss}'
            for block in blocks:
                assert all(is_plane_point(point, p) for point in block)
                assert len(set(block)) == len(block)
            assert set(probes) <= set(direct_tuple_sums(blocks, p)), 'Returned blocks do not cover all probes'
            nchecks += 1
    print(f'Constructive lower-bound covers: {nchecks}.')
    for p,ss in [(3,[2,2]),(3,[2,3]),(3,[2,2,2]),(5,[2,2])]:
        k,S,E=universe_minimum(p,ss)
        D=1+sum(n-1 for n in ss)
        M=prod(ss);upper=1
        while p**(2*D)*M**upper>=p**(2*upper):upper+=1
        assert D<k<=upper
        print(f'Universal probes p={p}, sizes={ss}: {S} distinct sumsets, minimum={k}, counting bounds [{D+1},{min(upper,p*p)}], witness={E}.')



def probe_reuse_tests():
    # Tiny abstract prime-order group map. Quotient coordinates are computed
    # independently by exact enumeration, never used as a membership oracle.
    p=3;trials=0;increments=0;probes=0;feasible_progress=0
    for m in (2,3):
        _,_,E=universe_minimum(p,[2]*m)
        V=list(product(range(p),repeat=m));T=list(product(range(2),repeat=m));zero=(0,)*m
        unit=[tuple(int(i==j) for i in range(m)) for j in range(m)]
        def lin(cs,bs):return tuple(sum(c*a[j] for c,a in zip(cs,bs))%p for j in range(m))
        def span(bs):return {lin(cs,bs) for cs in product(range(p),repeat=len(bs))}
        for lam in product((1,2),repeat=m):
            val=lambda x:sum(a*b for a,b in zip(x,lam))%p
            fs={g:[x for x in T if val(x)==g] for g in range(p)};assert all(fs.values())
            for reverse in (False,True):
                U=[]
                while len(U)<m-1:
                    basis=U.copy()
                    for e in unit:
                        if e not in span(basis):basis.append(e)
                    coords={lin(cs,basis):cs[len(U):] for cs in product(range(p),repeat=m)}
                    assert len(coords)==p**m
                    S={coords[x][:2] for x in T} # diagnostic only
                    assert any(z not in S for z in E)
                    order=list(reversed(E)) if reverse else E
                    for z in order:
                        x=lin(z,basis[len(U):len(U)+2]);assert coords[x][:2]==z
                        X=fs[val(x)][-1 if reverse else 0]
                        v=tuple((a-b)%p for a,b in zip(x,X));assert val(v)==0
                        probes+=1
                        if coords[v][:2]!=(0,0):
                            assert v not in span(U);U.append(v);increments+=1
                            feasible_progress+=int(z in S)
                            break
                    else:raise AssertionError('universal list failed')
                assert span(U)=={x for x in V if val(x)==0};trials+=1
    print(f'Universal-list reuse: {trials} synthetic calibrated completions, {increments} verified rank increments, {probes} decomposition probes; {feasible_progress} increments at probes actually in the projected support.')

def peval(f,x,p):
    z=0
    for a in reversed(f):z=(z*x+a)%p
    return z


def lagrange_weights(nodes,q,p):
    assert len(set(x%p for x in nodes))==len(nodes)
    ans=[]
    for x in nodes:
        a=1;b=1
        for y in nodes:
            if x!=y:a=a*(q-y)%p;b=b*(x-y)%p
        ans.append(a*pow(b,-1,p)%p)
    return ans


def digits(x,q,k):
    D=[]
    for j in range(k):x,t=divmod(x,q);D.append(t)
    assert not x
    return D


def digit_eval(a,z):
    x=0
    for v in reversed(a):x=z*x+v
    return x


def polynomial_cover(p,q,k,m,curves,explicit=True):
    d=max(len(f)-1 for F in curves for f in F)
    L=d*(k-1)+1;assert L<=m and L<=p
    nodes=list(range(L));w=lagrange_weights(nodes,q,p)
    ranges=[(q-1)*sum(z**j for j in range(k)) for z in nodes]
    bound=len(curves)*(max(ranges)+1)
    # This is an upper bound on cardinality, not a claim all points are distinct.
    assert bound**m<p*p
    if explicit:
        A=[{tuple(c*peval(f,x,p)%p for f in F) for x in range(R+1) for F in curves} for c,R in zip(w,ranges)]
        assert max(map(len,A))<=bound
        A += [{(0,0)} for _ in range(m-L)]
    else:A=None
    return nodes,w,ranges,bound,A


def curve_tests():
    rng=Random(1919);exhaustive=0;large=0
    # Prime and size choices are small, with all probe parameters checked.
    fixtures=[(10000019,16,3,5,[([0,1],[0,0,1])]),
              (2**61-1,10,3,7,[([1,2,3,4],[8,7,1,3]),([2,1],[3,0,1])]),
              (1009,4,2,2,[([1,3],[4,2])]),
              (1000000007,7,3,5,[([0,1,1],[3,2,4]),([3,7,2],[1,3,2])])]
    for p,q,k,m,curves in fixtures:
        nodes,w,ranges,bound,A=polynomial_cover(p,q,k,m,curves)
        assert prod(map(len,A))<p*p
        for t in range(q**k):
            ds=digits(t,q,k)
            for F in curves:
                parts=[]
                for i,(z,c,R) in enumerate(zip(nodes,w,ranges)):
                    v=digit_eval(ds,z);assert 0<=v<=R
                    a=tuple(c*peval(f,v,p)%p for f in F);assert a in A[i];parts.append(a)
                y=tuple(sum(a[j] for a in parts)%p for j in (0,1))
                assert y==tuple(peval(f,t,p) for f in F);exhaustive+=1
        print(f'Polynomial probe cover: p={p}, m={m}, degree={max(len(f)-1 for F in curves for f in F)}, {len(curves)} curves, {q**k} parameters each, max list size {max(map(len,A))}, product {prod(map(len,A))}.')
    # Large arithmetic-only check: no gigantic list or Cartesian product enumerated.
    p=2**521-1;q=2**12;k=9;m=17;curves=[([5,3],[1,7,2])]
    nodes,w,ranges,bound,_=polynomial_cover(p,q,k,m,curves,False)
    for _ in range(400):
        t=rng.randrange(q**k);ds=digits(t,q,k);F=curves[0]
        vv=[digit_eval(ds,z) for z in nodes]
        assert all(0<=v<=R for v,R in zip(vv,ranges))
        assert tuple(sum(c*peval(f,v,p) for c,v in zip(w,vv))%p for f in F)==tuple(peval(f,t,p) for f in F)
        large+=1
    print(f'Polynomial probe interpolation: {exhaustive} exhaustive small probe checks, {large} large arithmetic checks; large segment length 2^{(q**k).bit_length()-1}, upper list-size bits {bound.bit_length()}, product bound below field-plane size.')



def multivariate_and_line_tests():
    # Total-degree-3 polynomial pair in two parameters. Verify full grid.
    p=2**61-1;q=6;k=2;d=3;m=d*(k-1)+1
    nodes=list(range(m));weights=lagrange_weights(nodes,q,p)
    def F(x,y):return ((x*y+x*x*y+2*y*y+1)%p,(x*x*x+3*x*y*y+2*x+7)%p)
    R=[(q-1)*sum(z**j for j in range(k)) for z in nodes]
    A=[{tuple(c*a%p for a in F(x,y)) for x in range(t+1) for y in range(t+1)} for c,t in zip(weights,R)]
    cap=(max(R)+1)**2;assert max(map(len,A))<=cap and cap**m<p*p
    n=0
    for x,y in product(range(q**k),repeat=2):
        dx,dy=digits(x,q,k),digits(y,q,k);parts=[]
        for j,(z,c) in enumerate(zip(nodes,weights)):
            a=tuple(c*t%p for t in F(digit_eval(dx,z),digit_eval(dy,z)))
            assert a in A[j];parts.append(a)
        assert tuple(sum(a[j] for a in parts)%p for j in (0,1))==F(x,y);n+=1
    print(f'Multivariate polynomial cover: degree {d}, {m} lists, {n} parameter pairs; max list size {max(map(len,A))}; product below plane size.')
    rng=Random(982);tests=0
    for p,m,q,ell in [(257,2,17,3),(1009,3,11,9),(4093,3,16,8)]:
        assert q**m>=p and (ell*q)**m<p*p
        lines=[]
        for _ in range(ell):
            a=(rng.randrange(p),rng.randrange(p));v=(1,rng.randrange(p));lines.append((a,v))
        A=[]
        for i in range(m):
            A.append({((int(i==0)*a[0]+j*q**i*v[0])%p,(int(i==0)*a[1]+j*q**i*v[1])%p) for a,v in lines for j in range(q)})
        assert max(map(len,A))<=ell*q
        for a,v in lines:
            for t in range(p):
                ds=digits(t,q,m);parts=[((int(i==0)*a[0]+j*q**i*v[0])%p,(int(i==0)*a[1]+j*q**i*v[1])%p) for i,j in enumerate(ds)]
                assert all(part in A[i] for i,part in enumerate(parts))
                assert tuple(sum(x[j] for x in parts)%p for j in (0,1))==((a[0]+t*v[0])%p,(a[1]+t*v[1])%p)
                tests+=1
        print(f'Whole-line cover: p={p}, {ell} lines, {m} lists; max size {max(map(len,A))}; product {prod(map(len,A))} < {p*p}.')
    print(f'Whole-line identities: {tests}.')

def incidence_rows(q,m,t):
    return [tuple(peval(f,j,q) for j in range(m)) for f in product(range(q),repeat=t+1)]


def has_singleton(T):
    return any(1 in Counter(v).values() for v in zip(*T))


def girth_tests():
    rng=Random(521);pairs=0;subs=0;analytic=0
    for q,m,d in [(5,5,1),(7,7,2),(11,9,2)]:
        T=incidence_rows(q,m,d);assert len(T)==q**(d+1) and len(set(T))==len(T)
        for row in T:
            assert len(row) == m, f'Incorrect incidence row width for m={m}'
            assert all(type(x) is int and 0 <= x < q for x in row), f'Invalid incidence label for q={q}'
        # Evaluate the coefficient formula directly, independently of peval's Horner loop.
        expected_rows = {
            tuple(sum(c * j**i for i, c in enumerate(coeffs)) % q for j in range(m))
            for coeffs in product(range(q), repeat=d+1)
        }
        assert set(T) == expected_rows, f'Incorrect polynomial evaluation rows for q={q}, m={m}, d={d}'
        # Pair agreements <= degree, exactly the hypothesis of the singleton proof.
        if len(T)<=350:ps=combinations(T,2)
        else:ps=(rng.sample(T,2) for _ in range(2500))
        for a,b in ps:assert sum(x==y for x,y in zip(a,b))<=d;pairs+=1
        smax=(m-1)//d+1
        for s in range(2,smax+1):
            samples=combinations(T,s) if len(T)<=25 and s<=4 else (rng.sample(T,s) for _ in range(120))
            for E in samples:assert has_singleton(E);subs+=1
    # Exact integer verification of the elementary probability bound for sample parameters.
    for n in (16,64,256,1024):
        for m in (8,16,24):
            K=n**(m//8);s0=int(n**0.25)
            v=sum(Fraction(K**s * s**(m*s),n**(m*s//2)) for s in range(2,s0+1))
            assert v<1;analytic+=1
    print(f'Incidence-code checks: {pairs} pair-agreement comparisons, {subs} singleton-subset comparisons, {analytic} exact finite union-bound estimates.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('assertions required')
    range_tests();universal_tests();probe_reuse_tests();curve_tests();multivariate_and_line_tests();girth_tests()
