#!/usr/bin/env python3
"""Exact finite-geometry checks. No cryptographic instances or computations."""
from collections import Counter
from itertools import combinations, product
from math import comb
from operator import index
from unittest import TestCase
from random import Random
import numpy as np

if not __debug__:
    raise RuntimeError('Run with assertions enabled.')


# Modular products can overflow int64 before reduction, even in a 2x2 frame.
# operator.index also converts NumPy integer scalars to Python integers.
as_integers = np.frompyfunc(index, 1, 1)


def V(a,p): return as_integers(a) % p

def mul(a,b,p): return (V(a,p) @ V(b,p)) % p

def rref(a,p):
    a=V(a,p).copy(); k=0; piv=[]
    for j in range(a.shape[1]):
        cand=np.flatnonzero(a[k:,j])
        if not len(cand):continue
        i=k+int(cand[0]);a[[k,i]]=a[[i,k]]
        a[k]=a[k]*pow(int(a[k,j]),-1,p)%p
        for i in range(len(a)):
            if i!=k and a[i,j]:a[i]=(a[i]-a[i,j]*a[k])%p
        piv.append(j);k+=1
        if k==len(a):break
    return a,piv

def rank(a,p):return len(rref(a,p)[1])

def inv(a,p):
    n=len(a);assert a.shape==(n,n)
    b,piv=rref(np.hstack([V(a,p),identity(n)]),p)
    if piv[:n]!=list(range(n)):raise ValueError('singular')
    return b[:,n:]

def invertible(a,p):
    return rank(a,p)==len(a)

def tup(a):return tuple(map(int,a))

def plus(a,b,p):return tuple((x+y)%p for x,y in zip(a,b))

def sub(a,b,p):return tuple((x-y)%p for x,y in zip(a,b))

def identity(n):return np.eye(n,dtype=np.int64)

def mpow(a,k,p):
    b=identity(len(a))
    while k:
        if k&1:b=mul(b,a,p)
        a=mul(a,a,p);k//=2
    return b


def ball_volume(L,Q,t):
    return sum(comb(L,i)*(Q-1)**i for i in range(max(0,min(L,t)+1)))


def overlap_formula(L,Q,t,j):
    if j>L:return 0
    z=0
    for a in range(j+1):
        for b in range(j-a+1):
            c=j-a-b; h=min(t-j+a,t-j+b)
            if h>=0:
                z+=comb(j,a)*comb(j-a,b)*(Q-2)**c*ball_volume(L-j,Q,h)
    return z


def shell_gap(L,Q,t,j):
    # Difference between shifts of weights j and j+1, by equal-weight shells
    z=0
    for a in range(j//2+1):
        c=j-2*a; h=t-j+a
        if 0<=h<=L-1-j:
            z+=comb(j,a)*comb(j-a,a)*(Q-2)**c*comb(L-1-j,h)*(Q-1)**h
    return z


def hamming_checks():
    print('HAMMING BALL INTERSECTION IDENTITIES')
    tests=0;strict=0;binary_ties=0
    for Q,L,t in [(2,5,2),(3,5,2),(3,7,3),(4,5,2),(5,5,2),(3,3,1),(7,4,1)]:
        B=[x for x in product(range(Q),repeat=L) if sum(v!=0 for v in x)<=t]; S=set(B)
        assert len(B)==ball_volume(L,Q,t)
        seq=[]
        for j in range(L+1):
            h=(1,)*j+(0,)*(L-j)
            n=sum(plus(x,h,Q) in S for x in B)
            assert n==overlap_formula(L,Q,t,j);seq.append(n);tests+=1
        for j in range(L):
            assert seq[j]-seq[j+1]==shell_gap(L,Q,t,j)>=0
            if Q>2 and L>=2*t+1 and j<=2*t:
                assert seq[j]>seq[j+1];strict+=1
            if Q==2 and seq[j]==seq[j+1] and seq[j]>0:binary_ties+=1
        print(f'Q={Q}, L={L}, t={t}: overlaps={seq}')
    print(f'{tests} exact overlap counts; {strict} strict gaps; {binary_ties} positive binary ties.')


class Extension:
    def __init__(self,p,poly):
        self.p=p;self.a=len(poly);self.Q=p**self.a;self.poly=poly
        self.words=[np.array([(x//p**i)%p for i in range(self.a)],dtype=np.int64) for x in range(self.Q)]
    def mat(self,x):
        a=self.a;p=self.p;M=np.zeros((a,a),dtype=np.int64)
        for j in range(a):
            w=[0]*(2*a-1)
            for i,v in enumerate(self.words[x]):w[i+j]=int(v)
            for i in range(2*a-2,a-1,-1):
                for k,c in enumerate(self.poly):w[i-a+k]=(w[i-a+k]-w[i]*c)%p
            M[:,j]=w[:a]
        return M%p


def random_inv(n,p,rng):
    # Product of exact elementary row operations; no random search required.
    M=identity(n)
    for _ in range(5*n):
        i,j=rng.sample(range(n),2);M[i]=(M[i]+rng.randrange(1,p)*M[j])%p
    assert invertible(M,p);return M


def eval_groups(nodes,s,a,p,rng):
    """Nodes are matrices in a commutative algebra, None denotes infinity."""
    groups=[];d=s*a
    T=random_inv(d,p,rng)
    for t in nodes:
        B=(np.vstack([np.zeros((a,a),dtype=np.int64)]*(s-1)+[identity(a)]) if t is None
           else np.vstack([mpow(t,j,p) for j in range(s)]))
        groups.append(mul(mul(T,B,p),random_inv(a,p,rng) if a>1 else identity(a),p))
    rng.shuffle(groups)
    return groups


def matrix_algebra(gens,p,a):
    """Small incremental linear-span closure; only used for field diagnostics."""
    B=[]
    def push(x):
        if rank(np.stack([m.ravel() for m in B+[x]]),p)>len(B):B.append(x);return True
        return False
    push(identity(a))
    for x in gens:push(x)
    i=0
    while i<len(B):
        for x in gens:push(mul(B[i],x,p))
        i+=1
    F=np.stack([mpow(x,p,p).ravel() for x in B]);D=np.stack([x.ravel() for x in B]);
    r=rank(F,p);fix=len(B)-rank(F-D,p)
    return len(B),r,fix



def span_mats(rows,p,a):
    basis={}
    for X in rows:
        x=V(X,p).reshape(-1)
        for j,b in sorted(basis.items()):
            if x[j]:x=(x-x[j]*b)%p
        nz=np.flatnonzero(x)
        if len(nz):
            j=int(nz[0]);basis[j]=x*pow(int(x[j]),-1,p)%p
    return [b.reshape(a,a) for b in basis.values()]


def ideal_test(gens,p,a):
    # Quotient by the two-sided ideal of generator commutators.
    # It is commutative by construction. Nilpotence is an exact certificate
    # of simultaneous triangularizability over the algebraic closure.
    gens=span_mats([identity(a)]+gens,p,a)
    A=gens.copy()
    while True:
        B=span_mats(A+[mul(x,y,p) for x in A for y in gens],p,a)
        if len(B)==len(A):break
        A=B
    J=span_mats([(mul(x,y,p)-mul(y,x,p))%p for x,y in combinations(gens,2)],p,a)
    while J:
        B=span_mats(J+[mul(x,y,p) for x in J for y in gens]+[mul(y,x,p) for x in J for y in gens],p,a)
        if len(B)==len(J):break
        J=B
    P=J
    for _ in range(1,a):
        if not P:break
        P=span_mats([mul(x,y,p) for x in P for y in J],p,a)
    if P:raise ValueError('commutator ideal is not nilpotent')
    return A,J

def recognize_arc(groups,p):
    """Exact certificate, independent of how the candidate groups were obtained."""
    p=index(p)
    groups=[V(B,p) for B in groups]
    g=len(groups);d,a=groups[0].shape
    if d%a:raise ValueError('ambient dimension')
    s=d//a
    if s<2 or g<s:raise ValueError('not enough reference groups')
    assert all(B.shape==(d,a) for B in groups)
    if any(rank(B,p)!=a for B in groups):raise ValueError('group rank')
    M=np.hstack(groups[:s]);Mi=inv(M,p)
    if g==s:return s,a,None
    anchor=mul(Mi,groups[s],p)
    align=[inv(anchor[i*a:(i+1)*a],p) for i in range(s)]
    if g==s+1:return s,a,None
    G=[];I=V(identity(a),p);O=np.zeros_like(I)
    for B in groups[s:]:
        C=mul(Mi,B,p);C=[mul(align[i],C[i*a:(i+1)*a],p) for i in range(s)]
        C0=inv(C[0],p);G.append([inv(mul(x,C0,p),p) for x in C])
    if any(np.any(x!=I) for x in G[0]):raise ValueError('anchor')
    # First nonbasis group after the anchor defines the other reference nodes.
    t=G[1][1];nodes=[O]+[mul(inv((I-z)%p,p),(t-z)%p,p) for z in G[1][2:]]
    pars=[Z[1] for Z in G];finite=nodes+pars
    A,J=ideal_test([z for row in G for z in row],p,a)
    for x,y in combinations(finite,2):
        if not invertible((x-y)%p,p):
            raise ValueError('parameters are not unit separated')
    nz=0
    for Z,t in zip(G,pars):
        if np.any(Z[0]!=I):raise ValueError('first coordinate')
        for z,a0 in zip(Z[1:],nodes):
            R=(mul(z,(I-a0)%p,p)-(t-a0))%p
            if np.any(mpow(R,a,p)):raise ValueError('ratio identity fails on a diagonal')
            nz+=bool(np.any(R))
    # Field diagnostics optional; recognition itself does not require a field.
    diag=matrix_algebra(finite,p,a) if not J else ('triangular',len(A),len(J),nz)
    return s,a,diag


def certificate_regressions():
    for p, dtype in [(2**61 - 1, np.int64), (2**89 - 1, object)]:
        v = np.array([[2], [2]], dtype=dtype)
        opposite = np.array([[p-2], [p-2]], dtype=dtype)
        assert rank(np.hstack([v, opposite]), p) == 1, 'overflow in modular elimination'
        with TestCase().assertRaisesRegex(ValueError, '^singular$'):
            recognize_arc([v, opposite], p)
        # A nearby independent frame and nonreference groups must still certify.
        transform = np.array([[p-1, p-2], [p-3, p-5]], dtype=dtype)
        assert mul(transform, inv(transform, p), p).tolist() == [[1, 0], [0, 1]]
        columns = [(1, 0), (0, 1), (1, 1), (1, 2), (1, 3)]
        arc = [mul(transform, np.array([col]).T, p) for col in columns]
        assert recognize_arc(arc, p)[:2] == (2, 1)
    print('Over the primes 2^61-1 and 2^89-1, dependent frames were rejected and transformed five-group arcs certified.')

    # The g=s+1 shortcut must check every projection of the anchor.
    frame = identity(2)
    axes = [frame[:, i:i+1] for i in range(2)]
    with TestCase().assertRaisesRegex(ValueError, '^singular$'):
        recognize_arc(axes + [axes[0]], 5)

    # All other guards pass on these two fixtures. The indicated subset
    # ranks witness dependence independently of the certificate's decision.
    p = 7
    I = identity(2)
    O = np.zeros_like(I)
    axes = [np.vstack([I if i == j else O for i in range(3)]) for j in range(3)]
    beta = np.diag([3, 5])
    z = np.array([[0, 1], [6, 3]])
    groups = axes + [np.vstack([I, I, I]), np.vstack([I, 2*I, 3*I]),
                     np.vstack([I, inv(beta, p), inv(z, p)])]
    assert rank(np.hstack([groups[i] for i in (0, 3, 5)]), p) == 5
    with TestCase().assertRaisesRegex(ValueError, '^commutator ideal is not nilpotent$'):
        recognize_arc(groups, p)

    frame = identity(3)
    axes = [frame[:, i:i+1] for i in range(3)]
    triple = [np.array([v]).T for v in [(1, 1, 1), (1, 2, 3), (1, 3, 5)]]
    assert rank(np.hstack(triple), p) == 2
    with TestCase().assertRaisesRegex(ValueError, '^ratio identity fails on a diagonal$'):
        recognize_arc(axes + triple, p)
    print('Rejected singular anchor, nonnilpotent commutator ideal (rank 5/6), and failed ratio residual (rank 2/3).')


def arc_checks():
    print('\nHIDDEN COMMUTING-OPERATOR ARCS')
    rng=Random(271126);tests=0;subsets=0
    cases=[]
    for p,poly,s,g in [(3,(1,0),4,10),(2,(1,1,0),4,9),(5,(2,0),3,10),(7,(1,),3,8)]:
        F=Extension(p,poly)
        assert all(invertible(F.mat(x),p) for x in range(1,F.Q))
        nodes=[F.mat(x) for x in range(g-1)]+[None]
        cases.append(('field',p,F.a,s,nodes))
    p=7;I=identity(2);nil=np.array([[0,1],[0,0]],dtype=np.int64)
    cases.append(('dual numbers',p,2,3,[(i*I+(i*i%p)*nil)%p for i in range(p)]+[None]))
    p=5;cases.append(('product ring',p,2,3,[np.diag([i,(i*i*i)%p]) for i in range(p)]+[None]))
    # The cube map is a permutation of F_5, so each difference is a unit.
    for label,p,a,s,nodes in cases:
        for _ in range(3):
            G=eval_groups(nodes,s,a,p,rng);assert recognize_arc(G,p)[:2]==(s,a)
            diagnostic=recognize_arc(G,p)[2];n=0
            for I in combinations(range(len(G)),s):
                assert rank(np.hstack([G[i] for i in I]),p)==s*a;n+=1
            subsets+=n;tests+=1
        print(f'{label}: p={p}, block dimension={a}, s={s}, g={len(nodes)}, algebra(dim,Frobenius rank,fixed)={diagnostic}, {3*n} exhaustive maximal-subset checks')
    # Independent rational-normal curves on diagonal layers, arbitrary
    # upper-triangular off-diagonal lifts, and completely hidden coordinates.
    p=7;a=3;s=3;g=8;flag=0;flag_subsets=0;nonzero=0
    for _ in range(8):
        perms=[rng.sample(range(p),p) for j in range(a)]
        T=random_inv(s*a,p,rng);G=[]
        for i in range(g):
            blocks=[]
            for k in range(s):
                diag=([pow(perm[i],k,p) for perm in perms] if i<p else [int(k==s-1)]*a)
                B=np.diag(diag).astype(np.int64)
                for j in range(a):
                    for h in range(j+1,a):B[j,h]=rng.randrange(p)
                blocks.append(B)
            G.append(mul(mul(T,np.vstack(blocks),p),random_inv(a,p,rng),p))
        rng.shuffle(G)
        ans=recognize_arc(G,p);assert ans[:2]==(s,a) and ans[2][0]=='triangular'
        nonzero+=ans[2][-1]
        for I in combinations(range(g),s):
            assert rank(np.hstack([G[i] for i in I]),p)==s*a;flag_subsets+=1
        flag+=1
    assert nonzero>0
    print(f'{flag} noncommutative flag-lifted arcs; {flag_subsets} maximal-subset checks; {nonzero} nonzero nilpotent ratio residuals')
    negative=0
    for _ in range(20):
        F=Extension(3,(1,0));G=eval_groups([F.mat(x) for x in range(9)]+[None],4,2,3,rng)
        G[-1]=G[-2].copy()
        try:recognize_arc(G,3)
        except (ValueError,AssertionError):negative+=1
        else:raise AssertionError('duplicate group accepted')
    # A noncommutative graph perturbation in normalized coordinates.
    p=5;a=2;s=3;I=identity(a);Z=np.zeros((a,a),dtype=np.int64)
    axes=[np.vstack([I if i==j else Z for i in range(s)]) for j in range(s)]
    G=axes+[np.vstack([I,I,I]),np.vstack([I,2*I,3*I]),np.vstack([I,np.array([[1,1],[0,1]]),np.array([[1,0],[1,1]])])]
    try:recognize_arc(G,p)
    except (ValueError,AssertionError):negative+=1
    else:raise AssertionError('noncommuting perturbation accepted')
    print(f'{tests+flag} recognized layouts; {subsets+flag_subsets} maximal-subset checks; {negative} negative controls rejected.')


def make_ball(groups,p,t):
    a=groups[0].shape[1];d=groups[0].shape[0];zero=(0,)*d
    cores=[[tup(mul(B,np.array(x,dtype=np.int64),p)) for x in product(range(p),repeat=a) if any(x)] for B in groups]
    out={zero:()}
    for k in range(1,t+1):
        for I in combinations(range(len(groups)),k):
            for vs in product(*(cores[i] for i in I)):
                x=zero
                for v in vs:x=plus(x,v,p)
                assert x not in out
                out[x]=tuple(sorted(vs))
    return out,cores


def counts(S,p):
    # The binary case packs vectors into a machine integer, XOR is exact.
    d=len(S[0]);K=len(S)
    if p==2 and d<=32:
        def pack(x):return sum(v<<j for j,v in enumerate(x))
        A=np.array([pack(x) for x in S],dtype=np.uint32)
        keys,vals=np.unique(np.bitwise_xor(A[:,None],A[None,:]),return_counts=True)
        return {tuple((int(k)>>j)&1 for j in range(d)):int(c) for k,c in zip(keys,vals)}
    return Counter(sub(x,y,p) for x in S for y in S)


def split_core(R,p):
    z=(0,)*len(next(iter(R)));left=set(R);groups=[]
    while left:
        v=min(left);U={x for x in left if plus(x,v,p) in R or plus(x,v,p)==z}
        rows,piv=rref(np.array(sorted(U)),p);a=len(piv);B=rows[:a].T
        if len(U)!=p**a-1:raise ValueError('not a complete nonzero subspace')
        groups.append(B);left-=U
    return groups


def recover_blocks(S,p,L=None,Q=None,t=None,e=0,retention=(1,1)):
    C=counts(S,p);d=len(S[0]);z=(0,)*d
    if L is None:
        peak=max(v for x,v in C.items() if x!=z)
        R={x for x,v in C.items() if x!=z and v==peak}
        groups=split_core(R,p)
        L=len(groups);a=groups[0].shape[1];Q=p**a
        assert Q>=3 and all(B.shape[1]==a for B in groups)
        t=next((i for i in range(1,(L-1)//2+1) if ball_volume(L,Q,i)==len(S)),None)
        if t is None:raise ValueError('unrecognized ball size')
    else:
        c1=overlap_formula(L,Q,t,1);c2=overlap_formula(L,Q,t,2)
        a,b=retention
        R={x for x,v in C.items() if x!=z and 2*v*b*b>a*a*(c1+c2)}
        if len(R)!=L*(Q-1):raise ValueError('wrong core size')
        groups=split_core(R,p)
        if len(groups)!=L or any(p**B.shape[1]!=Q for B in groups):raise ValueError('wrong core component size')
    return C,R,groups,t


def factor(x,C,R,p,t):
    cur=x;z=(0,)*len(x);pieces=[]
    for _ in range(t):
        if cur==z:break
        v=cur if cur in R else max(sorted(R),key=lambda v:C.get(sub(cur,v,p),0))
        cur=sub(cur,v,p);pieces.append(v)
    if cur!=z:return None
    return tuple(sorted(pieces))


def check_shift_counts(observed, cores, p, t):
    """Enumerate formal supports in the planted blocks, independently of pair counts."""
    L = len(cores)
    Q = len(cores[0]) + 1
    zero = (0,) * len(cores[0][0])
    expected = [overlap_formula(L, Q, t, j) for j in range(2*t + 1)]
    remaining = dict(observed)

    def visit(shift, first, weight):
        assert remaining.pop(shift, None) == expected[weight], 'incorrect overlap for formal shift'
        if weight < 2*t:
            for i in range(first, L):
                for v in cores[i]:
                    visit(plus(shift, v, p), i + 1, weight + 1)

    visit(zero, 0, 0)
    assert not remaining, 'unexpected shifts in difference table'


def recovery_regressions():
    L, p = 7, 3
    axes = [identity(L)[:, i:i+1] for i in range(L)]
    ball, _ = make_ball(axes, p, 3)
    counts3, core, _, radius = recover_blocks(list(ball), p)
    x = (1, 1, 1, 0, 0, 0, 0)
    assert radius == 3, 'incorrect inferred radius'
    assert factor(x, counts3, core, p, radius) == ball[x]
    assert factor(x, counts3, core, p, radius-1) is None
    print('Inferred radius 3 from the complete dictionary and recovered a point requiring three blocks.')

    ball, cores = make_ball([identity(5)[:, i:i+1] for i in range(5)], p, 2)
    observed = counts(list(ball), p)
    swapped = observed.copy()
    for sign in (1, 2):
        u = (sign, sign, sign, 0, 0)
        v = (sign, sign, sign, sign, 0)
        assert (observed[u], observed[v]) == (12, 6)
        swapped[u], swapped[v] = observed[v], observed[u]
    assert Counter(swapped.values()) == Counter(observed.values())
    with TestCase().assertRaisesRegex(AssertionError, '^incorrect overlap for formal shift$'):
        check_shift_counts(swapped, cores, p, 2)
    print('Rejected counts 12 and 6 exchanged on opposite shift pairs with an unchanged histogram.')


def recovery_checks():
    print('\nOVERLAPPING HIDDEN BLOCK-BALL RECOVERY')
    rng=Random(106);total=0;counts_total=0
    fixtures=[]
    for p,L,t in [(3,5,2),(5,5,2),(3,9,2)]:
        s=L if L<9 else L-1
        A=random_inv(s,p,rng);groups=[A[:,j:j+1] for j in range(s)]
        if L>s:groups.append(mul(A,np.ones((s,1),dtype=np.int64),p))
        fixtures.append((f'ordinary p={p}',p,L,t,groups))
    F=Extension(3,(1,0));fixtures.append(('projected GF(9)',3,10,1,eval_groups([F.mat(x) for x in range(9)]+[None],4,2,3,rng)))
    F=Extension(2,(1,1,0));fixtures.append(('projected GF(8)',2,9,2,eval_groups([F.mat(x) for x in range(8)]+[None],8,3,2,rng)))
    for name,p,L,t,groups in fixtures:
        out,cores=make_ball(groups,p,t);S=list(out);Q=len(cores[0])+1
        C,R,G,radius=recover_blocks(S,p)
        assert radius==t, 'incorrect inferred radius'
        assert R==set(x for core in cores for x in core)
        assert {frozenset(core) for core in cores}=={frozenset(tup(mul(B,np.array(x),p)) for x in product(range(p),repeat=B.shape[1]) if any(x)) for B in G}
        check_shift_counts(C,cores,p,t)
        for x,v in out.items():assert factor(x,C,R,p,radius)==v;total+=1
        if name.startswith('projected GF'):
            s,a,diag=recognize_arc(G,p);assert s>=4*t and (diag is None or diag[0]==a)
        # Earlier rich-line criterion cannot separate the blocks when t>=2 and p>=3:
        # all directions on a line joining two distinct core blocks are present.
        if t>=2 and p>=3:
            x,y=cores[0][0],cores[1][0]
            assert all(plus(x,tuple(v*j%p for v in y),p) in out for j in range(p))
        print(f'{name}: L={L}, Q={Q}, t={t}, rank={rank(np.hstack(groups),p)}, K={len(S)}, difference shifts={len(C)}, exact decompositions={len(S)}')
        counts_total+=len(C)
    print(f'{total} exact block decompositions used inferred radii. All {counts_total} overlap counts matched their shifts.')


def perturbed_checks():
    print('\nDETERMINISTIC PERTURBATION GUARDS')
    rng=Random(201);n=0
    for p,L,t in [(3,7,2),(5,7,2),(3,8,3),(7,5,1)]:
        A=random_inv(L,p,rng);groups=[A[:,j:j+1] for j in range(L)];out,cores=make_ball(groups,p,t)
        last_gap = max(1, t-1)
        gap = min(shell_gap(L, p, t, j) for j in range(1, last_gap+1))
        e = min(3, (gap-1)//4)
        assert 1 <= e and 4*e < gap
        if t == 1:
            assert (gap, e) == (5, 1), 'incorrect radius-one edit budget'
        for _ in range(3):
            S=set(out);deleted=rng.sample(sorted(S),e-1);S.difference_update(deleted)
            outside=next(x for x in (tup([rng.randrange(p) for _ in range(L)]) for _ in range(200)) if x not in out)
            S.add(outside);C,R,G,radius=recover_blocks(list(S),p,L,p,t,e=e)
            assert R==set(x for core in cores for x in core)
            for x in S:
                v=factor(x,C,R,p,radius)
                if x in out:assert v==out[x];n+=1
                else:assert v is None
        print(f'p={p}, L={L}, t={t}, gap={gap}, edits={e}: 3 recovery trials, outliers rejected')
    print(f'{n} exact inlier decompositions after deterministic perturbations.')


def recover_center(S,B,p,e):
    # B is the exact recovered unshifted ball in dictionary coordinates.
    S=set(S);B=set(B);z=next(iter(B));assert e<len(S)
    candidates={sub(x,y,p) for x in sorted(S)[:e+1] for y in B}
    answers=[]
    for c in sorted(candidates):
        # Count overlap by hash lookup, rather than allocate a translated set.
        overlap=sum(plus(c,x,p) in S for x in B)
        if len(S)+len(B)-2*overlap<=e:answers.append(c)
    if len(answers)!=1:raise ValueError('no unique center')
    return answers[0]


def binary_layout(S,L,t):
    p=2;C=counts(S,p);z=(0,)*len(S[0])
    threshold_twice=overlap_formula(L,2,t,2)+overlap_formula(L,2,t,3)
    R={x for x,v in C.items() if x!=z and 2*v>threshold_twice}
    assert len(R)==comb(L+1,2), 'incorrect binary edge set size'
    e=min(R);f=next(x for x in sorted(R) if plus(x,e,p) in R)
    U={e,f}|{x for x in R if x!=e and x!=f and plus(x,e,p) in R and plus(x,f,p) in R}
    assert len(U)==L+1
    star=[x for x in U if sum(y!=x and plus(x,y,p) in R for y in U)>2]
    assert len(star)==L
    return [np.array(x,dtype=np.int64).reshape(-1,1) for x in sorted(star)]


def binary_edit_regression():
    L, t = 9, 2
    original, _ = make_ball([identity(L)[:, i:i+1] for i in range(L)], 2, t)
    missing = (1,) + (0,) * (L-1)
    sample = set(original) - {missing}
    assert 4 < shell_gap(L, 2, t, 2) == 12
    observed = counts(list(sample), 2)
    del observed[(0,) * L]
    peak = max(observed.values())
    assert sum(v == peak for v in observed.values()) == 28
    midpoint = (overlap_formula(L, 2, t, 2) + overlap_formula(L, 2, t, 3)) // 2
    edges = {x for x, v in observed.items() if v > midpoint}
    assert edges == {x for x in product(range(2), repeat=L) if 1 <= sum(x) <= 2}
    assert len(edges) == 45
    star = binary_layout(list(sample), L, t)
    ball, _ = make_ball(star, 2, t)
    center = recover_center(sample, ball, 2, 1)
    assert center == (0,) * L and set(ball) == set(original)
    print('After one binary radius-two deletion, 28 shifts attained the peak. The midpoint recovered all 45 edge shifts and the full ball.')


def affine_binary_checks():
    print('\nAFFINE CENTER AND BINARY STAR RECOVERY')
    rng=Random(1717);centers=0;stars=0;points=0
    for p,L,t,e in [(3,5,2,1),(5,5,2,2),(2,7,2,0),(2,13,3,0),(2,9,2,0)]:
        d=L if L<9 else L-1
        A=random_inv(d,p,rng);G=[A[:,i:i+1] for i in range(d)]
        if L>d:G.append(mul(A,np.ones((d,1),dtype=np.int64),p))
        B,_=make_ball(G,p,t);zero=(0,)*d
        for _ in range(3):
            c=tuple(rng.randrange(p) for i in range(d));S={plus(x,c,p) for x in B}
            if e:
                # Remove the actual center, then add arbitrary points if allowed.
                S.remove(c)
                for i in range(e-1):
                    v=next(x for x in (tuple(rng.randrange(p) for j in range(d)) for _ in range(1000)) if x not in S and sub(x,c,p) not in B)
                    S.add(v)
            if p==2:
                G2=binary_layout(list(S),L,t);stars+=1
            else:
                C,R,G2,_=recover_blocks(list(S),p,L,p,t,e=e)
            B2,_=make_ball(G2,p,t)
            center=recover_center(S,B2,p,e)
            assert len(S.symmetric_difference({plus(x,center,p) for x in B2}))<=e
            if p!=2:assert center==c
            else:
                # With an arbitrary star, even radii preserve center; odd radii
                # permit a recentered coordinate description.
                if t%2==0:assert center==c
            centers+=1;points+=len(S)
        print(f'p={p}, L={L}, t={t}, rank={d}, edits={e}: 3 exact affine recoveries')
    print(f'{centers} center recoveries, {stars} recovered binary star layouts, {points} input points checked.')



def recover_thinned_center(sample, ball, p):
    """Maximize the number of retained points explained by a translate."""
    ball = set(ball)
    scores = Counter(sub(x, b, p) for x in set(sample) for b in ball)
    leaders = scores.most_common(2)
    if not leaders or (len(leaders) == 2 and leaders[0][1] == leaders[1][1]):
        raise ValueError('no unique center')
    return leaders[0][0]


def translated_thinning_check():
    p, L, t = 5, 5, 2
    center = (1, 4, 2, 0, 3)
    original, _ = make_ball([identity(L)[:, i:i+1] for i in range(L)], p, t)
    rng = Random(0)
    sample = {plus(x, center, p) for x in sorted(original) if rng.randrange(4) < 3}
    assert len(sample) == 143 and center not in sample
    C, R, groups, radius = recover_blocks(list(sample), p, L, p, t, retention=(3, 4))
    ball, _ = make_ball(groups, p, radius)
    recovered = recover_thinned_center(sample, ball, p)
    assert recovered == center, 'incorrect thinned center'
    # Direct intersections over this small ambient space check both leading scores.
    rivals = (c for c in product(range(p), repeat=L) if c != recovered)
    best = sum(sub(x, recovered, p) in ball for x in sample)
    runner_up = max(sum(sub(x, c, p) in ball for x in sample) for c in rivals)
    assert (best, runner_up) == (143, 67)
    for x in sample:
        residual = sub(x, recovered, p)
        assert factor(residual, C, R, p, radius) == original[residual]
    with TestCase().assertRaisesRegex(ValueError, '^no unique center$'):
        recover_thinned_center([min(sample)], ball, p)
    print(f'Translated 3/4 thinning recovered the absent center {recovered} by {best} to {runner_up}. All {len(sample)} retained points factored. Tied centers were rejected.')


def thinning_examples():
    print('\nSEEDED THINNING CONTROLS (NO CLAIM THE ANALYTICAL GUARD HOLDS)')
    rng=Random(622);p=2;F=Extension(2,(1,1,0));t=2;L=9;Q=8
    G=eval_groups([F.mat(x) for x in range(8)]+[None],8,3,p,rng);B,cores=make_ball(G,p,t)
    R0=set(x for core in cores for x in core);n=0;fail=0
    for a,b in [(3,4),(1,2),(1,4),(1,8)]:
        S=[x for x in sorted(B) if rng.randrange(b)<a]
        try:
            C,R,G2,radius=recover_blocks(S,p,L,Q,t,retention=(a,b))
            assert R==R0
            assert all(factor(x,C,R,p,radius)==B[x] for x in S)
        except (ValueError,AssertionError):
            print(f'retention={a}/{b}, sample size={len(S)}: guard not met; recognition failed');fail+=1
        else:
            print(f'retention={a}/{b}, sample size={len(S)}: every retained vector recovered');n+=1
    print(f'{n} successful and {fail} failed thinning controls.')

def thinning_guards():
    from math import log, exp
    print('\nANALYTICAL THINNING GUARDS (NO LARGE BALL ENUMERATION)')
    for L,Q,t,d,q,theta in [(100,3,2,100,3,0.5),(10000,3,4,10000,3,0.5),(2000,7,5,2000,7,0.5),(2000,4,5,4000,2,0.5)]:
        J=max(1,t-1);Cs=[overlap_formula(L,Q,t,j) for j in range(J+2)]
        m=min((Cs[j]-Cs[j+1])**2/Cs[j] for j in range(1,J+1))
        logfail=log(6*J)+d*log(q)-theta**2*m/640
        print(f'L={L}, Q={Q}, t={t}, retention={theta}: log failure bound={logfail:.3f}, guard={logfail<0}')


if __name__=='__main__':
    hamming_checks()
    arc_checks()
    certificate_regressions()
    recovery_checks()
    recovery_regressions()
    perturbed_checks()
    affine_binary_checks()
    binary_edit_regression()
    thinning_examples()
    translated_thinning_check()
    thinning_guards()
    print('\nAll checks passed.')
