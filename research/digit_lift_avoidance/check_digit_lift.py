"""Exact synthetic verification of carry lifts and polynomial separation.
No elliptic-curve computation, asymptotic contractor, or remote-point solver.
Requires Python 3.10+ standard library, assertions enabled.
"""
from itertools import product, combinations
from math import prod, isqrt
from random import Random
from fractions import Fraction


def digits(x,q,k):
    out=[]
    for _ in range(k):
        x,a=divmod(x,q);out.append(a)
    return tuple(out)


def decode(v,q):
    return sum(a*q**j for j,a in enumerate(v))


def k_floor(r,q):
    k=0
    while q**(k+1)<=r:k+=1
    return k


def lift(p,q,k):
    return digits(p[0],q,k)+digits(p[1],q,k)


def lower_point(v,q,k):
    return decode(v[:k],q),decode(v[k:],q)


def addv(a,b,q):
    return tuple((x+y)%q for x,y in zip(a,b))


def sums(lists,q,d):
    S={(0,)*d}
    for L in lists:S={addv(s,x,q) for s in S for x in L}
    return S


def carry_word(xs,r,q,k):
    m=len(xs);z=sum(xs)%r;t=sum(xs)//r
    assert 0<=t<m
    xds=[digits(x,q,k) for x in xs];zd=digits(z,q,k);rd=digits(r,q,k)
    c=0;cs=[]
    for j in range(k):
        cs.append(c)
        s=sum(x[j] for x in xds)-t*rd[j]-zd[j]
        assert (s+c)%q==0
        assert (sum(x[j] for x in xds)-t*rd[j]+c)%q==zd[j]
        c=(s+c)//q
        assert -m<=c<=m
    return z,t,tuple(cs)


def carry_lists(r,q,k,m,ell):
    rd=digits(r,q,k);R={tuple((-t*x)%q for x in rd)+tuple((-u*x)%q for x in rd) for t in range(m) for u in range(m)}
    coords=[j for j in range(2*k) if j not in (0,k)]
    lists=[R]
    for a in range(0,len(coords),ell):
        block=coords[a:a+ell];L=set()
        for v in product(range(-m,m+1),repeat=len(block)):
            x=[0]*(2*k)
            for i,c in zip(block,v):x[i]=c%q
            L.add(tuple(x))
        assert len(L)<= (2*m+1)**len(block)
        lists.append(L)
    return lists


def check_carries():
    total=0
    for r,q,m in [(7,2,2),(11,3,3),(13,5,4),(17,3,3),(31,5,3)]:
        k=k_floor(r,q)
        for xs in product(range(r),repeat=m):
            carry_word(xs,r,q,k);total+=1
    rng=Random(8207);random_tests=0
    for r,q,m in [(127,11,2),(733,7,5),(1009,17,8),(4093,13,7),(65537,31,11)]:
        k=k_floor(r,q)
        for _ in range(1000):
            carry_word([rng.randrange(r) for _ in range(m)],r,q,k);random_tests+=1
    print(f'Scalar carry identities: {total} exhaustive and {random_tests} seeded random tuples passed.')


def check_plane():
    rng=Random(492);fixtures=0;sum_checks=0;outside_checks=0;factor_checks=0
    for r,q,m,ell in [(19,7,2,1),(127,11,2,1),(733,7,2,2),(1009,17,3,2)]:
        k=k_floor(r,q);d=2*k
        carries=carry_lists(r,q,k,m,ell)
        R=carries[0];C_fact=sums(carries,q,d)
        direct=set()
        for v in product(range(-m,m+1),repeat=2*(k-1)):
            e=(0,)+v[:k-1]+(0,)+v[k-1:]
            for w in R:direct.add(addv(e,w,q))
        assert C_fact==direct;factor_checks+=1
        for _ in range(3):
            blocks=[list({(rng.randrange(r),rng.randrange(r)) for j in range(3 if r!=733 else 2)}) for _ in range(m)]
            raw=sums(blocks,r,2);lifted=[{lift(a,q,k) for a in B} for B in blocks]
            sup=sums([*lifted,*carries],q,d)
            for z in raw:
                assert lift(z,q,k) in sup;sum_checks+=1
            assert prod(map(len,[*lifted,*carries]))<q**d
            # Check every absent lower digit vector, not just a sampled output.
            for v in product(range(q),repeat=d):
                if v not in sup:
                    x=lower_point(v,q,k)
                    assert x[0]<r and x[1]<r and x not in raw
                    outside_checks+=1
            fixtures+=1
    print(f'Plane lifts: {fixtures} complete instances; {factor_checks} exact carry factorization identities; {sum_checks} original-support inclusions; {outside_checks} safe absent-vector lifts.')


def rank(A,q,n):
    piv={}
    for row in A:
        v=[x%q for x in row]
        for i in sorted(piv):
            if v[i]:
                c=v[i];v=[(x-c*y)%q for x,y in zip(v,piv[i])]
        i=next((i for i,x in enumerate(v) if x),None)
        if i is not None:
            t=pow(v[i],-1,q);piv[i]=[x*t%q for x in v]
            if len(piv)==n:return n
    return len(piv)


def exponents(q,d,k):
    return [e for e in product(range(q),repeat=d) if sum(e)<=k]


def value(p,e,q):
    z=1
    for x,a in zip(p,e):z=z*pow(x,a,q)%q
    return z


def evaluations(points,exps,q):
    return [[value(p,e,q) for e in exps] for p in sorted(points)]


def anchored(blocks,q,d,k):
    base=tuple(sum(B[0][j] for B in blocks)%q for j in range(d))
    D=[{tuple((x-y)%q for x,y in zip(v,B[0])) for v in B} for B in blocks]
    out={base}
    for s in range(1,min(k,len(blocks))+1):
        for J in combinations(range(len(blocks)),s):
            for xs in product(*(D[j] for j in J)):
                out.add(tuple((base[j]+sum(x[j] for x in xs))%q for j in range(d)))
    return out


def check_truncation():
    rng=Random(778);tests=0;strict=0
    for q,d,m in [(3,3,4),(3,4,4),(5,3,4),(7,3,5)]:
        for _ in range(12):
            blocks=[[tuple(rng.randrange(q) for _ in range(d)) for _ in range(3)] for _ in range(m)]
            S=sums(blocks,q,d)
            for k in range(4):
                sk=anchored(blocks,q,d,k);exps=exponents(q,d,k)
                assert sk<=S
                assert rank(evaluations(sk,exps,q),q,len(exps))==rank(evaluations(S,exps,q),q,len(exps))
                strict+=len(sk)<len(S);tests+=1
    print(f'Anchored truncation: {tests} exact evaluation-span identities; {strict} strictly smaller restricted supports.')


def ball(q,d,t):
    S=set()
    for j in range(t+1):
        for I in combinations(range(d),j):
            for v in product(range(1,q),repeat=j):
                x=[0]*d
                for i,a in zip(I,v):x[i]=a
                S.add(tuple(x))
    return S


def check_balls():
    fixtures=0
    for q,d,m,t in [(3,4,2,1),(3,5,3,1),(5,4,2,1),(3,5,2,2)]:
        h=m*t;assert h<d
        B=ball(q,d,t);S=sums([B]*m,q,d)
        assert S==ball(q,d,h)
        exps=[e for e in product(range(q),repeat=d) if sum(a!=0 for a in e)<=h]
        assert len(exps)==len(S)
        assert rank(evaluations(S,exps,q),q,len(exps))==len(exps)
        z=(1,)*(h+1)+(0,)*(d-h-1)
        assert z not in S and all(prod(x[:h+1])%q==0 for x in S)
        fixtures+=1
        print(f'  F_{q}^{d}, m={m}, t={t}: |B_t|={len(B)}, |mB_t|={len(S)}, reduced support-degree <= {h} evaluation rank {len(exps)}, first annihilator degree {h+1}.')
    print(f'Hamming ball obstructions: {fixtures} exact support and full interpolation-rank checks.')


def check_sparse_syndrome():
    rng=Random(45);tests=0;nonfull=0
    for q,d,N,w in [(3,5,6,1),(3,6,7,2),(5,5,6,1)]:
        for _ in range(10):
            cols=[tuple(rng.randrange(q) for _ in range(d)) for i in range(N)]
            # One sparse coefficient vector for each reachable syndrome.
            hit={}
            for s in range(w+1):
                for I in combinations(range(N),s):
                    for coeff in product(range(1,q),repeat=s):
                        v=[0]*N
                        for i,c in zip(I,coeff):v[i]=c
                        z=tuple(sum(v[j]*cols[j][i] for j in range(N))%q for i in range(d))
                        hit[z]=v
            assert len(hit)<q**d
            z=next(z for z in product(range(q),repeat=d) if z not in hit)
            # Direct independent enumeration of all vectors in the promised radius.
            for v in product(range(q),repeat=N):
                if sum(a!=0 for a in v)<=w:
                    assert tuple(sum(v[j]*cols[j][i] for j in range(N))%q for i in range(d))!=z
            nonfull+=rank(cols,q,d)<d;tests+=1
    print(f'Sparse syndrome fixtures: {tests} exhaustive missing-syndrome checks ({nonfull} non-full-row-rank dictionaries).')



def kernel_vector(A,q,n):
    piv={}
    for row in A:
        v=[x%q for x in row]
        for i in sorted(piv):
            if v[i]:
                c=v[i];v=[(x-c*y)%q for x,y in zip(v,piv[i])]
        i=next((j for j,x in enumerate(v) if x),None)
        if i is not None:
            t=pow(v[i],-1,q);piv[i]=[x*t%q for x in v]
            if len(piv)==n:return None
    f=next(i for i in range(n) if i not in piv);u=[0]*n;u[f]=1
    for i in sorted(piv,reverse=True):u[i]=-sum(x*y for x,y in zip(piv[i][i+1:],u[i+1:]))%q
    assert any(u) and all(sum(x*y for x,y in zip(row,u))%q==0 for row in A)
    return u


def greedy_nonzero(coeffs,exps,q):
    terms={e:c%q for e,c in zip(exps,coeffs) if c%q};assert terms
    out=[]
    for j in range(len(exps[0])):
        for a in range(q):
            D={}
            for e,c in terms.items():D[e[1:]]=(D.get(e[1:],0)+c*pow(a,e[0],q))%q
            D={e:c for e,c in D.items() if c}
            if D:terms=D;out.append(a);break
        else:raise AssertionError('reduced polynomial vanished at every residue')
    return tuple(out)


def check_algebraic_selector():
    # A product of two nonlinear graph blocks has every nonzero scalar
    # projection surjective, yet a quadratic polynomial separates it.
    for q in (3,5,7):
        d=6;B=[(x,y,x*y%q) for x in range(q) for y in range(q)]
        A=[v+(0,0,0) for v in B];C=[(0,0,0)+v for v in B]
        S=sums([A,C],q,d)
        L=exponents(q,d,1);assert rank(evaluations(S,L,q),q,len(L))==len(L)
        L=exponents(q,d,2);sk=anchored([A,C],q,d,2)
        u=kernel_vector(evaluations(sk,L,q),q,len(L));assert u is not None
        z=greedy_nonzero(u,L,q);assert z not in S
        assert all(sum(c*value(x,e,q) for c,e in zip(u,L))%q==0 for x in S)
        assert sum(c*value(z,e,q) for c,e in zip(u,L))%q!=0
        # Verify every scalar projection of one graph; a nonzero functional
        # on the product has a nonzero restriction to at least one graph.
        for a in product(range(q),repeat=3):
            if any(a):assert len({sum(x*y for x,y in zip(a,v))%q for v in B})==q
        print(f'  F_{q} graph product: {len(S)} support points, full affine span, all scalar projections surjective, quadratic selector succeeded.')
    # Carry thickening itself makes low-degree annihilators impossible on
    # a valid lifted input with full low-coordinate projection.
    for r,q in [(127,11),(173,13)]:
        k=2;m=2;A=[(a,b) for a in range(q) for b in range(q)];C=[(0,0)]
        lists=[[lift(x,q,k) for x in A],[lift(x,q,k) for x in C],*carry_lists(r,q,k,m,1)]
        assert prod(map(len,lists))<q**(2*k)
        S=sums(lists,q,2*k);assert len(S)<q**(2*k)
        L=exponents(q,2*k,2*m);assert rank(evaluations(S,L,q),q,len(L))==len(L)
        print(f'  Carry-box obstruction: r={r}, q={q}, {len(S)} lifted points, full rank on all {len(L)} monomials of total degree <= {2*m}.')
    print('Algebraic selector: 3 quadratic successes with surjective scalar projections, 2 exact carry-box negative controls.')


def check_parameters():
    # Conservative integer bit bounds, not floating point asymptotics.
    # q may be any prime in [B,2B], r any b-bit prime; Bertrand supplies q.
    for b,m,sigma,c in [(2**20,4,Fraction(1,3),16),(2**20,32,Fraction(1,3),16),(2**22,64,Fraction(1,2),24),(2**24,128,Fraction(3,4),48),(2**25,256,Fraction(1,2),24)]:
        B=(2*m+1)**c;ql=B.bit_length()-1;qu=(2*B).bit_length();a=(2*m+1).bit_length()
        kmax=b//ql;ell=(b//m)//a
        assert ell>0
        chunks=(2*max(0,kmax-1)+ell-1)//ell
        mp=m+1+chunks
        logM=(1+sigma)*b+2*m.bit_length()+2*max(0,kmax-1)*a
        # If each main block has <= ceil(2^((1+sigma)b/m)) entries,
        # N' <= (m+chunks+1)*2^(ceil((1+sigma)b/m)+2) for these bounds.
        x=(1+sigma)*b/m
        logN=(x.numerator+x.denominator-1)//x.denominator+2+mp.bit_length()
        logspan=mp*(1+qu+logN)
        lower=2*(b-1-qu)
        assert logM<lower
        if m==4:
            assert logspan>=lower
            print(f'  b={b}, m={m}, c={c}: volume guard passes; conservative uncolored sparse-span guard fails (finite rounding).')
        else:
            assert logspan<lower
            print(f'  b={b}, m={m}, c={c}: at most {mp} lifted lists, carry chunks <= {chunks}; volume and uncolored sparse-span guards hold by integer bit bounds.')
    print('Parameter guards: 4 passing conservative large-bitlength cases and 1 finite failure retained; no large prime or factor-base enumeration.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    check_carries();check_plane();check_truncation();check_balls();check_sparse_syndrome();check_algebraic_selector();check_parameters()
