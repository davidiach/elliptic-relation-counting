"""Exact synthetic group and integer matrix checks; no cryptographic-system tests."""
from collections import defaultdict
from itertools import product
from random import Random
from math import prod
from sympy import Matrix, ZZ, zeros
from sympy.matrices.normalforms import smith_normal_form


def rank_mod(rows,r,n):
    basis={}
    for row in rows:
        a=[int(x)%r for x in row]
        for j in sorted(basis):
            if a[j]:
                c=a[j]; a=[(x-c*y)%r for x,y in zip(a,basis[j])]
        j=next((j for j,x in enumerate(a) if x),None)
        if j is not None:
            c=pow(a[j],-1,r);basis[j]=[x*c%r for x in a]
    return len(basis)


def tuple_fibers(r,blocks,cpb):
    slots=[i for i in range(len(blocks)) for _ in range(cpb)]
    fibers=defaultdict(set)
    for X in product(*(range(len(blocks[i])) for i in slots)):
        g=sum(blocks[i][x] for i,x in zip(slots,X))%r
        c=tuple(sum(i==j and x==k for i,x in zip(slots,X)) for j,F in enumerate(blocks) for k in range(1,len(F)))
        fibers[g].add(c)
    return fibers


def rows_at(fibers,r,t):
    return sorted({tuple(y-x for x,y in zip(X,Y)) for g,U in fibers.items() for X in U for Y in fibers.get((g+t)%r,())})


def snf_diagonal(A):
    M=smith_normal_form(Matrix(A),domain=ZZ)
    return [abs(int(M[j,j])) for j in range(min(M.rows,M.cols))]


def lattice_checks():
    rng=Random(8371); ntests=0; target_tests=0; maxrows=0
    for r,m,k,cpb in [(5,2,2,2),(7,2,3,1),(11,3,3,1),(13,3,3,1),(13,2,3,2),(17,3,2,2),(19,3,3,1),(23,3,3,2)]:
        done=0
        for _ in range(14):
            blocks=[rng.sample(range(1,r),k) for _ in range(m)]
            fibers=tuple_fibers(r,blocks,cpb)
            if len(fibers)<r:continue
            A=rows_at(fibers,r,1);D=sum(len(F)-1 for F in blocks);maxrows=max(maxrows,len(A))
            assert snf_diagonal(A)==[1]*D
            aug=[list(a)+[-1] for a in A]
            assert snf_diagonal(aug)==[1]*D+[r]
            C=[tuple(x-y for x,y in zip(a,A[0])) for a in A]
            assert snf_diagonal(C)==[1]*(D-1)+[r]
            # Target family alone has same raw rank; its differences span
            # the zero-label subspace. Add one genuine unit-shift row to force normalization.
            for t in [1,r-1,rng.randrange(1,r)]:
                T=rows_at(fibers,r,t)
                assert rank_mod(T,r,D)==D
                M=[list(a)+[-1] for a in T]+[list(A[0])+[0]]
                assert rank_mod(M,r,D+1)==D+1
                target_tests+=1
            ntests+=1;done+=1
            if done==8:break
    print(f'Integral-lattice theorem: {ntests} full-support fixtures; all three Smith forms exact; {target_tests} target-augmentation rank checks; maximum distinct shift rows {maxrows}.')


def target_obstruction():
    r=13;blocks=[[8,2,3],[6,12,3],[8,5,7]];D=6
    B=[[-1,1,0,-1,1,0],[0,-1,1,0,-1,0],[1,-1,0,-1,0,1],[0,1,1,0,0,0],[0,1,0,-1,-1,0],[0,0,0,-1,1,-1]]
    u=(1,3,10,9,7,11); lam=tuple((x-F[0])%r for F in blocks for x in F[1:])
    assert abs(int(Matrix(B).det()))==r and rank_mod(B,r,D)==D-1
    assert all(sum(x*y for x,y in zip(a,u))%r==0 for a in B)
    assert all(sum(x*y for x,y in zip(a,lam))%r==1 for a in B)
    # Actual E/F7, y^2=x^3+3, with group order 13.
    p=7
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;z,w=Q
        if x==z and (y+w)%p==0:return None
        t=((3*x*x)*pow(2*y,-1,p) if P==Q else (w-y)*pow(z-x,-1,p))%p
        X=(t*t-x-z)%p;return X,(t*(x-X)-y)%p
    pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-3)%p==0]
    P=pts[0];G=[None]
    for _ in range(r-1):G.append(add(G[-1],P))
    assert len(pts)==r-1 and set(G)=={None,*pts} and add(G[-1],P) is None
    for cpb in (1,2):
        slots=[i for i in range(3) for _ in range(cpb)]
        f=defaultdict(set)
        for X in product(range(3),repeat=len(slots)):
            Q=None
            for i,x in zip(slots,X):Q=add(Q,G[blocks[i][x]])
            g=G.index(Q); assert g==sum(blocks[i][x] for i,x in zip(slots,X))%r
            c=tuple(sum(i==j and x==k for i,x in zip(slots,X)) for j in range(3) for k in (1,2))
            f[g].add(c)
        assert len(f)==r and set(map(tuple,B))<=set(rows_at(f,r,1))
        assert snf_diagonal(rows_at(f,r,1))==[1]*D
        missing=[];counts=[]
        for q in range(1,r):
            A=rows_at(f,r,q); scores={sum(x*y for x,y in zip(a,u))%r for a in A}
            assert all(sum(x*y for x,y in zip(a,lam))%r==q for a in A)
            if 0 not in scores:missing.append(q);counts.append(len(A))
        assert 5 in missing
        print(f'Order-13 curve, copies per block={cpb}: all individual target rows ambiguous for shifts {missing}; distinct target-row counts {counts}; rational basis determinant ±13, modular rank 5/6.')


def gram_weighted(rows,n):
    G=zeros(n,n); s=zeros(n,1); total=0
    for a,k in rows:
        if not k:continue
        v=Matrix(a);G+=k*v*v.T;s+=k*v;total+=k
    return G,s,total


def gram_trade(r):
    # b+3 coordinates: b prefix digits, one auxiliary, one balancing
    # coordinate z, one constant column c. Actual log labels all equal one.
    b=r.bit_length();n=b+3;aux=b;z=b+1;c=b+2
    def E(j):return Matrix([int(i==j) for i in range(n)])
    e0=E(0);v=e0-E(z); ec=E(c)
    bits=[int(x) for x in bin(r)[3:]];assert len(bits)==b-1
    u=[1]
    for x in bits:u.append((2*u[-1]+x)%r)
    u=Matrix(u+[1,0,0]); assert u[b-1]==0
    raw=[]
    for j,bit in enumerate(bits,1):
        a=2*E(j-1)+bit*E(aux)-E(j)-(1+bit)*E(z)
        raw.append((a,2**(b-1-j)))
    raw.extend([(E(0)-E(aux),r-2**(b-1)),(E(b-1)-E(z),1),(zeros(n,1),r)])
    assert sum(w for _,w in raw)==2*r
    G,S,M=gram_weighted(raw,n);assert S==r*v
    A0=[(a-v,w) for a,w in raw];G1,S1,M1=gram_weighted(A0,n)
    assert G1==G and S1==-r*v
    B=[];A=[]
    for family,out in [(raw,B),(A0,A)]:
        for a,w in family:
            for sgn in [1,-1]:out.append((tuple(ec+sgn*a),w))
    G0,S0,M0=gram_weighted(B,n);G1,S1,M1=gram_weighted(A,n)
    assert G0==G1 and S0==S1 and M0==M1==4*r
    assert all(sum(row)==1 for family in [B,A] for row,w in family if w)
    assert all(sum(x*x for x in row)<=11 and max(map(abs,row))<=2 for family in [B,A] for row,w in family if w)
    assert all(sum(x*y for x,y in zip(row,u))%r==0 for row,w in B if w)
    assert all(sum(x*y for x,y in zip(row,u))%r in [1,r-1] for row,w in A if w)
    anchor=tuple(E(z)); rb=[a for a,w in B if w]+[anchor];ra=[a for a,w in A if w]+[anchor]
    assert Matrix(rb).rank()==Matrix(ra).rank()==n
    assert rank_mod(rb,r,n)==n-1 and rank_mod(ra,r,n)==n
    return n,B,A,anchor


def gram_tests():
    from sympy import nextprime
    rs=[3,5,7,13,31,101,1009,65537]+[int(nextprime(2**k)) for k in (20,32,64)]
    for r in rs:
        n,B,A,z=gram_trade(r)
    print(f'Exact Gram trade: {len(rs)} primes through {rs[-1].bit_length()} bits (multiplicities compressed); identical Gram and RHS, full rational rank, modular ranks n−1 versus n, squared row norm <=11.')
    # Finite oracle simulation with dyadic padding and a known common anchor.
    r=13;n,B,A,z=gram_trade(r)
    def expand(rows):return [tuple(int(x) for x in a) for a,w in rows for _ in range(w)]
    B=expand(B);A=expand(A);L=1<<(len(B)-1).bit_length();c=tuple(int(i==n-1) for i in range(n))
    B+= [c]*(L-len(B));A+=[c]*(L-len(A))
    K=8
    from collections import Counter
    def oracle(hidden,lo,hi):
        rows=Counter(a for j in range(K) for a in (A if j==hidden else B)[max(0,lo-j*L):max(0,min(L,hi-j*L))])
        return gram_weighted(list(rows.items()),n)[0]
    # Whole aligned blocks cannot identify the hidden block; a strict prefix
    # can cut at most one block, and any interval at most two.
    Gbase=oracle(-1,0,K*L)
    Gblock=oracle(-1,0,L)
    for h in range(K):
        assert oracle(h,0,K*L)==Gbase
        assert oracle(h,h*L,(h+1)*L)==Gblock
        for j in range(K):
            if j!=h: assert oracle(h,j*L,j*L+1)==oracle(-1,j*L,j*L+1)
    print(f'Oracle fixture: {K} hidden-block possibilities of {L} rows; exact whole-block and untouched-prefix replies agree.')


if __name__=='__main__':
    lattice_checks();target_obstruction();gram_tests()
