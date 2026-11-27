"""Exact synthetic tests for structured remote syndromes; Python standard library."""
from collections import Counter
from itertools import product, combinations
from fractions import Fraction
from math import comb, prod, isqrt
from random import Random


def add(a,b,q): return tuple((x+y)%q for x,y in zip(a,b))
def sub(a,b,q): return tuple((x-y)%q for x,y in zip(a,b))
def scale(a,c,q): return tuple(c*x%q for x in a)
def dot(a,b,q): return sum(x*y for x,y in zip(a,b))%q

def rref(rows,q,d):
    a=[list(map(lambda x:x%q,row)) for row in rows if any(row)]; piv=[]; i=0
    for j in range(d):
        z=next((z for z in range(i,len(a)) if a[z][j]),None)
        if z is None:continue
        a[i],a[z]=a[z],a[i]; c=pow(a[i][j],-1,q); a[i]=[x*c%q for x in a[i]]
        for z in range(len(a)):
            if z==i or not a[z][j]:continue
            c=a[z][j]; a[z]=[(x-c*y)%q for x,y in zip(a[z],a[i])]
        piv.append(j); i+=1
        if i==len(a):break
    return tuple(tuple(row) for row in a[:i]),tuple(piv)

def basis(rows,q,d): return rref(rows,q,d)[0]
def rank(rows,q,d): return len(basis(rows,q,d))
def matrank(v,q,s): return rank([v[i*s:(i+1)*s] for i in range(s)],q,s)
def normal(v,B,piv,q):
    a=list(v)
    for b,j in zip(B,piv):
        c=a[j]
        if c:a=[(x-c*y)%q for x,y in zip(a,b)]
    return tuple(a)

def nullspace(rows,q,d):
    B,piv=rref(rows,q,d); out=[]
    for j in range(d):
        if j in piv:continue
        v=[0]*d;v[j]=1
        for row,i in zip(B,piv):v[i]=-row[j]%q
        out.append(tuple(v))
    return tuple(out)

def intersection(U,V,q,d):return nullspace((*nullspace(U,q,d),*nullspace(V,q,d)),q,d)
def lincomb(B,coeff,q,d):
    out=[0]*d
    for b,c in zip(B,coeff):
        if c:out=[(x+c*y)%q for x,y in zip(out,b)]
    return tuple(out)

def span(B,q,d):return tuple(lincomb(B,a,q,d) for a in product(range(q),repeat=len(B)))
def units(d):return [tuple(int(i==j) for j in range(d)) for i in range(d)]
def extend(first,rest,q,d):
    out=list(first)
    for v in rest:
        if rank([*out,v],q,d)>len(out):out.append(v)
    return out

def coordinate_map(B,q,d):
    """Rows B are a square basis; solve x = coeff*B using exact inverse."""
    assert len(B)==d
    rows=[list(v)+list(e) for v,e in zip(B,units(d))]
    R,piv=rref(rows,q,d);assert piv==tuple(range(d))
    inverse=[x[d:] for x in R]
    return lambda x:lincomb(inverse,x,q,d)

def random_basis(d,q,rng):
    E=units(d)
    for _ in range(8*d):
        i,j=rng.sample(range(d),2); E[i]=add(E[i],scale(E[j],rng.randrange(q),q),q)
    assert rank(E,q,d)==d
    return E

def hidden_ball(q,s,k):
    d=s*s;zero=(0,)*d
    if k==1:
        # Canonical first factor: first nonzero coordinate one.
        U=[u for u in product(range(q),repeat=s) if any(u) and next(x for x in u if x)==1]
        V=[v for v in product(range(q),repeat=s) if any(v)]
        return {zero,*[tuple(x*y%q for x in u for y in v) for u in U for v in V]}
    assert q**d<=100000
    return {v for v in product(range(q),repeat=d) if matrank(v,q,s)<=k}


def overlaps(B,q,d):
    """Count ordered differences, deduplicating the input set."""
    B = {tuple(x) for x in B}
    if q==2:
        def pack(v):return sum(x<<i for i,x in enumerate(v))
        D=[pack(v) for v in B];C=Counter(x^y for x in D for y in D if x!=y)
        return {tuple(v>>i&1 for i in range(d)):n for v,n in C.items()}
    return Counter(sub(x,y,q) for x in B for y in B if x!=y)


def recover_grid(R,q,s):
    d=s*s; z=(0,)*d; A=set(R)|{z};h=min(R)
    def near(a):return {x for x in A if add(x,a,q) in A}
    def line(a):return {scale(a,c,q) for c in range(q)}
    N=near(h);a=min(N-line(h));L={x for x in N if add(x,a,q) in A};U=(N-L)|line(h)
    assert len(L)==len(U)==q**s and rank(L,q,d)==rank(U,q,d)==s
    ar=extend([h],sorted(L),q,d);br=extend([h],sorted(U),q,d)
    assert len(ar)==len(br)==s
    G=[[z]*s for _ in range(s)]
    for j in range(s):G[0][j]=ar[j]
    for i in range(s):G[i][0]=br[i]
    # The two factor spaces at a grid point are identified by adjacency.
    Lj=[None]+[(near(br[i])-U)|line(br[i]) for i in range(1,s)]
    Uj=[None]+[(near(ar[j])-L)|line(ar[j]) for j in range(1,s)]
    for i in range(1,s):
        for j in range(1,s):
            inter=Lj[i]&Uj[j];assert len(inter)==q and z in inter
            base=add(add(h,ar[j],q),br[i],q)
            candidates=[x for x in inter if any(x) and add(base,x,q) in A]
            assert len(candidates)==1;G[i][j]=candidates[0]
    G=[x for row in G for x in row];assert rank(G,q,d)==d
    return G


def rank1_counts_formula(q,s,k):
    def n(h):return prod((q**s-q**i)**2 for i in range(h))//prod(q**h-q**i for i in range(h))
    Rk=n(k);tail=sum(n(h) for h in range(k));
    return Rk,tail


def analytic_checks():
    tests=0
    for q in (7,11,17,31):
        for k in range(1,8):
            for s in (8*k,8*k+1,12*k+1):
                R,tail=rank1_counts_formula(q,s,k)
                assert Fraction(tail,R)<=Fraction(8,q**(2*(s-k)+1))
                low=Fraction(q**k-1,q**s-1)
                upper=Fraction(12,q**(2*(s-3*k)))+Fraction(tail,R)
                assert upper<low/16
                # Rank-k conditional factor distribution and all lower ranks.
                pfull=1-Fraction(2,q**(s-k)*(q-1))
                hyper=(Fraction(1,q)+Fraction(q-1,q**(k+1)))/pfull+Fraction(tail,R)
                assert hyper<Fraction(1,3)
                tests+=1
    print(f'Analytic integer/rational guards: {tests} parameter triples; autocorrelation separation <1/16 and hyperplane bound <1/3.')


def rank_recovery_checks():
    rng=Random(2807);cases=0;preserves=0;noisy=0
    B0=hidden_ball(2,2,1); C0=overlaps(B0,2,4)
    assert C0 == {v: 6 for v in product(range(2), repeat=4) if any(v)}, "incorrect binary overlap counts"
    print("  Retained negative control F_2 2x2: all 15 nonzero shifts have overlap 6; rank-one recovery fails as expected outside the guard.")
    for q,s,k,w in [(2,3,1,2),(3,2,1,1),(7,2,1,1),(2,3,2,1),(2,5,1,2)]:
        d=s*s;B=hidden_ball(q,s,k);hidden=random_basis(d,q,rng);T=lambda v:lincomb(hidden,v,q,d);inv=coordinate_map(hidden,q,d)
        D={T(v) for v in B};R={T(v) for v in B if matrank(v,q,s)==1}
        C=overlaps(D,q,d);mx=max(C.values());assert {v for v,n in C.items() if n==mx}==R
        if (q,s,k) in [(2,3,1),(3,2,1)]:
            counts = {}
            for shift in product(range(q), repeat=d):
                if any(shift):
                    count = sum(add(x, shift, q) in D for x in D)
                    if count:
                        counts[shift] = count
            repeated = [list(x) for x in [*D, min(R), min(R), (0,)*d]]
            assert overlaps(repeated, q, d) == C == counts, "duplicate atoms changed overlap counts"
            print(f'  Duplicate atoms: exact overlap counts unchanged at all {q**d-1} nonzero shifts over F_{q}.')
        G=recover_grid(R,q,s);co=coordinate_map(G,q,d)
        for x in D:
            assert matrank(co(x),q,s)==matrank(inv(x),q,s);preserves+=1
        I=tuple(int(i//s==i%s) for i in range(d));y=lincomb(G,I,q,d)
        assert matrank(inv(y),q,s)==s and w*k<s
        # Exact capped sparse spans, independent of the rank reasoning.
        if q==2:
            pack=lambda v:sum(x<<i for i,x in enumerate(v))
            X={pack(x) for x in D};Y={0}
            for _ in range(w):Y={x^z for x in X for z in Y}
            assert pack(y) not in Y
        else:
            X={scale(v,c,q) for v in D for c in range(q)};Y={(0,)*d}
            for _ in range(w):Y={add(x,z,q) for x in X for z in Y}
            assert y not in Y
        cases+=1
        print(f'  F_{q}, {s}x{s}, k={k}, w={w}: {len(D)} distinct atoms, {len(R)} rank-one vectors recovered, full-rank remote witness.')
        # Noise controls: exact rank-one recovery from a thinned set and outliers.
        if q==2 and s==5:
            C1=mx; drop=set(rng.sample(sorted(D),4));
            pool=[T(tuple(rng.randrange(q) for _ in range(d))) for _ in range(20)]
            outside={v for v in pool if matrank(inv(v),q,s)>k};outs=set(sorted(outside)[:2]);assert len(outs)==2
            E=(D-drop)|outs;C2=overlaps(E,q,d);M=max(C2.values());R2={v for v,n in C2.items() if 2*n>M};assert R2==R
            G2=recover_grid(R2,q,s);co2=coordinate_map(G2,q,d)
            assert {x for x in E if matrank(co2(x),q,s)>k}==outs
            codeword,cc,nn=outlier_remote(q,s,k,w,[co2(x) for x in outs]);assert cc<nn
            y2=lincomb(G2,codeword,q,d);X={pack(x) for x in E};Y={0}
            for _ in range(w):Y={x^z for x in X for z in Y}
            assert pack(y2) not in Y
            print(f'  Noisy F_2 5x5: {len(drop)} deletions, {len(outs)} arbitrary outliers; all {len(R2)} rank-one differences recovered, {cc} distinct outlier centers, {nn} code probes, remote witness verified.')
            noisy+=1
    print(f'Unknown-coordinate checks: {cases} exact ball recoveries; {preserves} individual rank equalities; {noisy} noisy outlier recoveries.')


def polynomial_remainder(a,b,q):
    a=list(a)
    for i in range(len(a)-1,len(b)-2,-1):
        c=a[i]%q
        for j in range(len(b)):a[i-len(b)+1+j]=(a[i-len(b)+1+j]-c*b[j])%q
    return tuple(a[:len(b)-1])

def irreducible(q,s):
    if s==1:return (0,1)
    ds=[(*z,1) for j in range(1,s//2+1) for z in product(range(q),repeat=j) if z[0]]
    for z in product(range(q),repeat=s):
        if not z[0]:continue
        f=(*z,1)
        if all(any(polynomial_remainder(f,g,q)) for g in ds):return f
    raise AssertionError('no irreducible')

def field_ops(q,s):
    f=irreducible(q,s);z=(0,)*s;one=(1,)+(0,)*(s-1)
    def mul(a,b):
        out=[0]*(2*s-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%q
        return polynomial_remainder(out,f,q)
    def power(a,n):
        ans=one
        while n:
            if n&1:ans=mul(ans,a)
            a=mul(a,a);n//=2
        return ans
    return mul,power

def mr_code(q,s,ell,count):
    d=s*s;mul,power=field_ops(q,s);U=units(s);powers=[[power(x,q**i) for x in U] for i in range(ell)]
    assert count<=q**(s*ell)
    for n in range(count):
        a=[]
        for i in range(s*ell):n,x=divmod(n,q);a.append(x)
        coeff=[tuple(a[i*s:(i+1)*s]) for i in range(ell)];cols=[]
        for j in range(s):
            v=(0,)*s
            for i in range(ell):v=add(v,mul(coeff[i],powers[i][j]),q)
            cols.append(v)
        yield tuple(x for row in zip(*cols) for x in row)

def sparse_centers(O,q,w,d):
    X={(0,)*d}
    for i in range(1,min(w,len(O))+1):
        for I in combinations(O,i):
            for a in product(range(1,q),repeat=i):X.add(lincomb(I,a,q,d))
    return X

def outlier_centers(outliers, q, w, d):
    """Cover the sparse centers using the smaller enumeration bound."""
    generators = basis(outliers, q, d)
    span_size = q**len(generators)
    count = len(outliers)
    sparse_bound = sum(comb(count, j) * (q-1)**j for j in range(min(w, count)+1))
    if span_size < sparse_bound:
        return set(span(generators, q, d))
    return sparse_centers(outliers, q, w, d)

def center_checks():
    # These dependent atoms give 1,941 coefficient choices but only 16 centers.
    dense = list(product(range(2), repeat=4))
    assert outlier_centers(dense[1:], 2, 4, 4) == set(dense)

    # Here the full span strictly contains the sparse centers. The two routes
    # give observably different covers, even after deduplication.
    atoms = [v+(0,0) for v in product(range(2), repeat=3) if any(v)] + units(5)[3:]
    choices = {(0,)*5, *atoms}
    reachable = {add(a, b, 2) for a in choices for b in choices}
    full = outlier_centers(atoms, 2, 2, 5)
    assert len(reachable) == 25 and reachable < full
    assert full == set(product(range(2), repeat=5))
    sparse = outlier_centers(units(5), 2, 1, 5)
    assert sparse == {(0,)*5, *units(5)}

    # Embed the same dependent atoms as nonzero even-weight diagonals (rank >= 2).
    def diagonal(v):
        entries = (*v, sum(v)%2)
        return tuple(entry if i == j else 0 for i, entry in enumerate(entries) for j in range(6))
    witness, count, probes = outlier_remote(2, 6, 1, 2, [diagonal(v) for v in atoms])
    assert count == len(full) and probes == count+1, "constructor did not use the cheaper span"
    assert all(matrank(sub(witness, diagonal(v), 2), 2, 6) > 2 for v in reachable)
    print(f'Center enumeration: 16 span elements for 1,941 coefficient choices; '
          f'{len(full)} span elements cover {len(reachable)} sparse centers; '
          f'{probes} probes certify the matrix witness; sparse route visits {len(sparse)} choices for a 32-element span.')

def outlier_remote(q,s,k,w,O):
    d=s*s;R=w*k;ell=s-2*R;assert ell>=1
    centers=outlier_centers(O,q,w,d);n=len(centers)+1;assert n<=q**(s*ell)
    probes=list(mr_code(q,s,ell,n));assert len(set(probes))==n
    assert all(matrank(sub(a,b,q),q,s)>=2*R+1 for a,b in combinations(probes,2))
    missing=[x for x in probes if all(matrank(sub(x,c,q),q,s)>R for c in centers)]
    assert missing;return missing[0],len(centers),n

def outlier_checks():
    rng=Random(118);checks=0
    fixtures = [(2,3,1,1,3),(3,3,1,1,4),(2,5,1,1,9),(2,5,1,2,3),(3,5,1,2,2),
                (2,4,1,1,16)]
    for q,s,k,w,L in fixtures:
        d=s*s;O=[]
        while len(O)<L:
            x=tuple(rng.randrange(q) for _ in range(d))
            if matrank(x,q,s)>k and x not in O:O.append(x)
        x,c,n=outlier_remote(q,s,k,w,O);checks+=n*(n-1)//2
        # Build the true sparse centers one slot at a time, independently of
        # the constructor's choice between a full span and sparse combinations.
        zero = (0,)*d
        choices = {scale(atom, a, q) for atom in O for a in range(q)} | {zero}
        actual_centers = {zero}
        for _ in range(w):
            actual_centers = {add(v, a, q) for v in actual_centers for a in choices}
        assert len(actual_centers) <= c and n == c+1
        assert len(x) == d and all(0 <= a < q for a in x)
        assert all(matrank(sub(x, v, q), q, s) > w*k for v in actual_centers), "invalid returned syndrome"
        if (q,s,k,w,L) == (2,4,1,1,16):
            assert n == q**s + 2
            probes = list(mr_code(q, s, s-2*w*k, n))
            square, mixed = probes[q**s:]
            mul, _ = field_ops(q, s)
            for v in product(range(q), repeat=s):
                squared = mul(v, v)
                for matrix, value in [(square, squared), (mixed, add(v, squared, q))]:
                    image = tuple(dot(matrix[j:j+s], v, q) for j in range(0, d, s))
                    assert image == value, "incorrect higher coefficient"
            print(f'  Second coefficient: probes {q**s} and {q**s+1} evaluate to x^2 and x+x^2 on all {q**s} extension elements.')
        print(f'  Rank-metric outlier test F_{q}, s={s}, w={w}: {L} outliers, {c} centers, {n} probes, certified missing syndrome.')
    print(f'Rank-metric pair-distance checks: {checks} across {len(fixtures)} outlier fixtures; all returned witnesses independently checked.')




def affine_noise_check():
    # Translation leaves ordered-difference counts unchanged.  Recenter at a
    # point with a maximal radius-two neighborhood; all other differences
    # are charged as outliers.  The common offset is one extra center direction.
    q,s,k,w=2,5,1,1;d=s*s;rng=Random(177);B=hidden_ball(q,s,k)
    hidden=random_basis(d,q,rng);T=lambda v:lincomb(hidden,v,q,d);inv=coordinate_map(hidden,q,d)
    shift=tuple(rng.randrange(q) for _ in range(d));D={add(T(v),shift,q) for v in B}
    pool=[tuple(rng.randrange(q) for _ in range(d)) for _ in range(40)]
    outside=[v for v in pool if v not in D];E=(D-set(rng.sample(sorted(D),3)))|set(outside[:2])
    C=overlaps(E,q,d);M=max(C.values());R={v for v,n in C.items() if 2*n>M}
    assert R=={T(v) for v in B if matrank(v,q,s)==1}
    G=recover_grid(R,q,s);co=coordinate_map(G,q,d)
    Q={co(x):x for x in E};Qkeys=sorted(Q)
    # There is an inlier, so stopping when a sufficiently large neighborhood
    # is found suffices; the theorem allows a full maximization instead.
    h0=None;O=None
    for a in Qkeys:
        far=[sub(x,a,q) for x in Qkeys if matrank(sub(x,a,q),q,s)>2*k]
        if len(far)<=2:h0=a;O=far;break
    assert h0 is not None and O is not None
    # Conservative center set: all field multiples of the shared offset and
    # all combinations of at most w high-rank difference atoms.
    centers = outlier_centers(O, q, w, d)
    Z={add(scale(h0,a,q),c,q) for a in range(q) for c in centers}
    ell=s-4*w*k;assert ell>=1 and len(Z)<q**(s*ell)
    probes=list(mr_code(q,s,ell,len(Z)+1))
    ys=[x for x in probes if all(matrank(sub(x,c,q),q,s)>2*w*k for c in Z)]
    assert ys; y=lincomb(G,ys[0],q,d)
    assert y not in E and any(y) # w=1 in F_2
    print(f'Hidden affine/noisy ball: {len(E)} supplied atoms, 3 deletions, 2 additions; tensor layout recovered; {len(O)} recentered outliers, {len(Z)} centers, {len(probes)} probes, exact one-atom exclusion.')

def extra_rank_boundaries():
    # All hyperplanes, up to a scalar, for a finite q>=7 rank-one example.
    q,s=7,2;d=s*s;B=hidden_ball(q,s,1)
    normals=[v for v in product(range(q),repeat=d) if any(v) and next(x for x in v if x)==1]
    mx=max(sum(dot(v,x,q)==0 for x in B) for v in normals)
    assert 3*mx<len(B)
    print(f'Exact hyperplane control F_7 2x2: {len(normals)} projective hyperplanes, maximum {mx}/{len(B)} < 1/3.')
    # Every invertible linear reshape of the binary 4-dimensional space.
    # The weight<=2 dictionary defeats all rank<=1 certificates for w=1.
    q,d=2,4;H=[v for v in product(range(q),repeat=d) if sum(v)<=2];n=0
    vs=[v for v in product(range(q),repeat=d) if any(v)]
    for a in vs:
        for b in vs:
            if b==a:continue
            V={a,b,add(a,b,q)}
            for c in vs:
                if c in V:continue
                W=V|{add(x,c,q) for x in V}|{c}
                for e in vs:
                    if e in W:continue
                    G=[a,b,c,e];assert any(matrank(lincomb(G,x,q,d),q,2)==2 for x in H);n+=1
    assert n==20160
    print(f'Hamming-ball rank-certificate obstruction: all {n} invertible F_2 reshapes of dimension four checked.')
    rng=Random(142);q=2;s=5;d=s*s;B=hidden_ball(q,s,1);R=B-{(0,)*d}; successes=0
    # These toy samples lie outside the stated asymptotic q,s guard;
    # failures are counted rather than suppressed.
    for _ in range(3):
        S={x for x in sorted(B) if rng.random()<.95};C=overlaps(S,q,d);mx=max(C.values());R0={v for v,n in C.items() if 2*n>mx}
        if R0==R:
            G=recover_grid(R0,q,s);co=coordinate_map(G,q,d)
            assert all(matrank(co(x),q,s)<=1 for x in S);successes+=1
    print(f'Random-thinning toy controls: {successes}/3 exact rank-one recoveries at retention 0.95 (binary controls outside asymptotic guard).')


# Binary-tree cost-profile dynamic programming.
class Node:
    def __init__(self,i=None,left=None,right=None):self.i=i;self.left=left;self.right=right

def tree(lo,hi):
    if hi-lo==1:return Node(i=lo)
    mid=(lo+hi)//2;return Node(left=tree(lo,mid),right=tree(mid,hi))

def exact_cost_profiles(tables,q,d,w,verify=False,*,root=None):
    """Use the supplied leaf-labeled tree, or a balanced tree in table order."""
    if not tables:
        raise ValueError('expected at least one leaf table')
    if root is None:
        root = tree(0, len(tables))
    cap=w+1;stats={'width':0,'max_profiles':0,'pair_tests':0,'leaf_candidates':0,'node_profiles':0}
    nodes = [root]
    seen = set()
    remaining = set(range(len(tables)))
    for n in nodes:
        if not isinstance(n, Node) or n in seen:
            raise ValueError('expected a binary tree of distinct nodes')
        seen.add(n)
        if n.i is None:
            nodes.extend((n.left, n.right))
        elif n.left is not None or n.right is not None or type(n.i) is not int or n.i not in remaining:
            raise ValueError('expected one leaf per table index')
        else:
            remaining.remove(n.i)
    if remaining:
        raise ValueError('missing table indices in the supplied tree')
    for n in reversed(nodes):
        if n.i is not None:n.V=basis(tables[n.i],q,d)
        else:n.V=basis((*n.left.V,*n.right.V),q,d)
    outside = {root: ()}
    for n in nodes:
        out = outside.pop(n)
        n.S=intersection(n.V,out,q,d);n.si=span(n.S,q,d);stats['width']=max(stats['width'],len(n.S))
        if n.i is None:
            outside[n.left] = basis((*out, *n.right.V), q, d)
            outside[n.right] = basis((*out, *n.left.V), q, d)
    def missing_coset(V,S,support):
        B,piv=rref(S,q,d);Q,qp=rref([normal(v,B,piv,q) for v in V],q,d)
        def idx(v):
            nv=normal(v,B,piv,q);return sum(nv[j]*q**i for i,j in enumerate(qp))
        seen=sorted({idx(x) for x in support});a=0
        for x in seen:
            if x!=a:break
            a+=1
        if a>=q**len(Q):return None
        co=[]
        for _ in Q:a,c=divmod(a,q);co.append(c)
        return lincomb(Q,co,q,d)
    def combine(n):
        if n.i is not None:
            T={x:c for x,c in tables[n.i].items() if c<=w};ys={sub(x,s,q) for x in T for s in n.si};stats['leaf_candidates']+=len(ys)
            P={}
            for y in sorted(ys):P.setdefault(tuple(T.get(add(y,s,q),cap) for s in n.si),y)
            miss=missing_coset(n.V,n.S,T)
            if miss is not None:P.setdefault((cap,)*len(n.si),miss)
            n.P=P
            if verify:n.actual=T
        else:
            pindex={x:i for i,x in enumerate(n.si)};pairs=[]
            for i,x in enumerate(n.left.si):
                for j,y in enumerate(n.right.si):
                    s=add(x,y,q)
                    if s in pindex:pairs.append((i,j,pindex[s]))
            P={}
            for A,a in n.left.P.items():
                for B,b in n.right.P.items():
                    C=[cap]*len(n.si)
                    for i,j,k in pairs:C[k]=min(C[k],A[i]+B[j])
                    P.setdefault(tuple(C),add(a,b,q));stats['pair_tests']+=1
            n.P=P
            if verify:
                T={}
                for a,i in n.left.actual.items():
                    for b,j in n.right.actual.items():
                        if i+j<=w:
                            x=add(a,b,q);T[x]=min(i+j,T.get(x,cap))
                n.actual=T
        stats['max_profiles']=max(stats['max_profiles'],len(n.P));stats['node_profiles']+=len(n.P)
        if verify:
            actual={tuple(n.actual.get(add(y,s,q),cap) for s in n.si) for y in span(n.V,q,d)}
            assert set(n.P)==actual
            for P,y in n.P.items():assert P==tuple(n.actual.get(add(y,s,q),cap) for s in n.si)
    for n in reversed(nodes):
        combine(n)
    return root.P,stats

def profile_checks():
    rng=Random(9302);checks=0;profiles=0
    for q,d,K in [(2,4,5),(3,3,4),(2,5,7),(3,4,5)]:
        for trial in range(12):
            tables=[]
            for i in range(K):
                if trial%2:
                    a=tuple(rng.randrange(q) for _ in range(d));T={scale(a,c,q):int(c!=0) for c in range(q)};T[(0,)*d]=0
                else:
                    T={tuple(rng.randrange(q) for _ in range(d)):rng.randrange(3) for j in range(4)}
                tables.append(T)
            P,S=exact_cost_profiles(tables,q,d,3,verify=True);checks+=1;profiles+=S['node_profiles']
    print(f'Exact tree profile comparisons: {checks} random weighted/sparse inputs, {profiles} distinct stored node profiles independently checked.')
    # Keep table indices fixed. The two trees put equal columns on opposite
    # sides of the root, or together in sibling leaves.
    d = 4
    tables = [{(0,)*d: 0, v: 1} for v in units(d)*2]
    pairs = [Node(left=Node(i=i), right=Node(i=i+d)) for i in range(d)]
    paired = Node(left=Node(left=pairs[0], right=pairs[1]),
                  right=Node(left=pairs[2], right=pairs[3]))
    widths = []
    for supplied, width in [(tree(0, 8), 4), (paired, 1)]:
        root_profiles, stats = exact_cost_profiles(tables, 2, d, 2, verify=True, root=supplied)
        assert stats['width'] == width, "incorrect supplied tree width"
        assert set(root_profiles) == {(0,), (1,), (2,), (3,)}
        for (cost,), witness in root_profiles.items():
            assert cost == min(3, sum(witness))
        widths.append(stats['width'])
    print(f'Same eight labeled leaves: separated copies have width {widths[0]}, '
          f'sibling copies have width {widths[1]}; both realize capped costs 0 through 3.')
    # 64 one-choice edge lists, each atom covers at most one coordinate.
    q=2;d=65;E=units(d);z=(0,)*d;tables=[{z:0,E[i]:0,E[i+1]:0} for i in range(d-1)]
    P,S=exact_cost_profiles(tables,q,d,0);assert (1,) in P and sum(P[(1,)])==65
    print(f'Large colored chain: {d} dimensions, {len(tables)} lists, 3^{len(tables)} tuple choices, width {S["width"]}, max profiles {S["max_profiles"]}, {S["pair_tests"]} profile pairs; missing all-ones vector found.')
    # Overlapping groups with four private coordinates and two boundary coordinates.
    q=3;g=12;h=4;d=g*h+g-1;E=units(d);z=(0,)*d;tables=[]
    for i in range(g):
        ix=list(range(i*h,(i+1)*h))
        if i>0:ix.append(g*h+i-1)
        if i+1<g:ix.append(g*h+i)
        local=span([E[j] for j in ix],q,d);tables.append({x:int(any(x)) for x in local})
    P,S=exact_cost_profiles(tables,q,d,g);assert max(x[0] for x in P)==g
    y=P[(g,)];assert all(any(y[i*h:(i+1)*h]) for i in range(g))
    K=sum(len(T)-1 for T in tables)
    print(f'Overlapping full local dictionaries: q={q}, d={d}, K={K} nonzero atoms, exact covering radius {g}, width {S["width"]}, max profiles {S["max_profiles"]}, {S["pair_tests"]} profile pairs; no enumeration of ambient 3^{d}.')

if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    analytic_checks();rank_recovery_checks();center_checks();outlier_checks();affine_noise_check();extra_rank_boundaries();profile_checks()
