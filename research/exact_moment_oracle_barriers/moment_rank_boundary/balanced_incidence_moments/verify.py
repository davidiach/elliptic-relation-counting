#!/usr/bin/env python3
"""Exact checks for balanced integer/ternary interpolation and cycle trades."""
from collections import Counter
from itertools import product, permutations, combinations, combinations_with_replacement
from math import comb, factorial, prod

if not __debug__:
    raise RuntimeError('Run with assertions enabled')


def simplex(d, w):
    if d == 0:
        yield ()
        return
    for x in range(w+1):
        for a in simplex(d-1,w-x):
            yield (x,)+a


def mass(a):
    return sum(max(x,0) for x in a)


def support(a):
    return tuple(i for i,x in enumerate(a) if x)


def balanced(n,w,unit=False):
    # Enumerate signed disjoint supports from weak compositions of the mass.
    out={(0,)*n}
    for k in range(1,w+1):
        comps=[x for x in simplex(n,k) if sum(x)==k and (not unit or max(x)<=1)]
        for a in comps:
            for b in comps:
                if not any(x*y for x,y in zip(a,b)):
                    out.add(tuple(x-y for x,y in zip(a,b)))
    return sorted(out,key=lambda a:(-len(support(a)),-mass(a),a))


def separator(a,w,unit=False):
    """Raw triangular separator, integer numerator; normalization is diagonal."""
    S=support(a); s=len(S)
    if s==0:
        return 0,lambda x:1
    if s==2 and not unit:
        i,j=S; v=a[i]
        roots=tuple(t for t in range(-w,w+1) if t not in (0,v))
        return 2+len(roots),lambda x:x[i]*x[j]*prod(x[i]-t for t in roots)
    assert s>=2
    P=tuple(i for i in S if a[i]>0)
    N=tuple(i for i in S if a[i]<0)
    R=P if len(P)<=len(N) else N
    alpha=tuple(abs(v) for v in a)
    def f(x):
        return (prod(x[i]*prod(x[i]*x[i]-j*j for j in range(1,alpha[i])) for i in S)
                *prod(alpha[i]+(1 if a[i]>0 else -1)*x[i] for i in R))
    degree=2*sum(alpha)-s+len(R)
    return degree,f


def check_triangular():
    print('BALANCED BLOCK TRIANGULAR SEPARATORS')
    print('n | w | unit | points | max degree | exact triangular zeros')
    total=0
    cases=[(2,1,False),(3,1,False),(5,1,False),(4,2,False),
           (6,2,False),(5,3,False),(4,4,False),(5,2,True),
           (7,3,True),(8,4,True)]
    for n,w,unit in cases:
        X=balanced(n,w,unit)
        degree=0; count=0
        for j,a in enumerate(X):
            k,f=separator(a,w,unit)
            assert f(a)!=0
            degree=max(degree,k)
            for x in X[j+1:]:
                assert f(x)==0,(a,x)
                count+=1
        bound=3*w if unit else (3 if w==1 else 4*w-2)
        assert degree<=bound
        total+=count
        print(f'{n} | {w} | {int(unit)} | {len(X)} | {degree} | {count}')
    print(f'{total} exact triangular zero checks.')


def det(A):
    A=[list(x) for x in A]; n=len(A); sign=1;last=1
    if n==0:return 1
    for k in range(n-1):
        i=next((i for i in range(k,n) if A[i][k]),None)
        if i is None:return 0
        if i!=k:A[k],A[i]=A[i],A[k];sign=-sign
        z=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                a,b=divmod(A[i][j]*z-A[i][k]*A[k][j],last)
                assert b==0;A[i][j]=a
        for i in range(k+1,n):A[i][k]=0
        last=z
    return sign*A[-1][-1]


def inv_parity(a):
    return (-1)**sum(x>y for i,x in enumerate(a) for y in a[i+1:])


def moment(C,e):
    return sum(v*prod(a**b for a,b in zip(x,e)) for x,v in C.items())


def alternant(k):
    vals=tuple(range(-k,k+1))
    return {a:inv_parity(a) for a in permutations(vals)}


def alt_moment(k,e):
    vals=tuple(range(-k,k+1))
    return det([[v**j for v in vals] for j in e])


def check_alternants():
    print('\nNORMALIZED VANDERMONDE TRADES')
    print('k | coordinates | mass | first moment degree | rows | exhaustive low-degree moments')
    count=0
    for k in range(1,4):
        n=2*k+1; w=k*(k+1)//2; t=n*(n-1)//2
        C=alternant(k)
        assert len(C)==factorial(n) and sum(C.values())==0
        assert all(sum(a)==0 and mass(a)==w and abs(v)==1 for a,v in C.items())
        assert sum(v==1 for v in C.values())==factorial(n)//2
        es=tuple(range(n))
        expected=prod(j-i for i in range(n) for j in range(i+1,n))
        assert alt_moment(k,es)==moment(C,es)==expected!=0
        checked=0
        if k<=2:
            for e in simplex(n,t-1):
                assert alt_moment(k,e)==moment(C,e)==0
                checked+=1
        else:
            # An arbitrary low-degree exponent vector must repeat an exponent.
            # Check one exhaustive smaller simplex and all repeated-pair variants
            # of the witness. The general argument is the determinant identity.
            for e in simplex(n,5):
                assert alt_moment(k,e)==0
                checked+=1
            for i in range(n):
                for j in range(i):
                    e=list(es);e[i]=e[j]
                    assert alt_moment(k,e)==0
                    checked+=1
        count+=checked
        print(f'{k} | {n} | {w} | {t} | {len(C)} | {checked}')
    print(f'{count} zero alternant moments checked.')
    print('w | general lower degree | general upper degree | signed-unit degree')
    for w in range(1,21):
        k=1
        while (k+1)*(k+2)//2<=w:k+=1
        s=w-k*(k+1)//2
        lower=k*(2*k+1)+3*s
        upper=3 if w==1 else 4*w-2
        assert 0<=s<=k and lower==4*w-k-s<=upper
        assert lower>=4*w-2*k
        print(f'{w} | {lower} | {upper} | {3*w}')
    # Small exact Cartesian products, with all monomials below the sum tested
    # via factorization. No need to list the exponentially larger full table.
    A=alternant(1); B=alternant(1)
    C={a+b:x*y for a,x in A.items() for b,y in B.items()}
    assert len(C)==36 and all(mass(a)==2 for a in C)
    for e in simplex(6,5):
        assert moment(C,e)==moment(A,e[:3])*moment(B,e[3:])==0
    witness=(0,1,2)*2
    assert moment(C,witness)==4
    print('Two triangle factors: all degree <=5 moments zero, degree-6 witness = 4.')


def rank(rows,p):
    piv={}
    for row in rows:
        a=[x%p for x in row]
        for j,b in sorted(piv.items()):
            if a[j]:
                z=a[j];a=[(x-z*y)%p for x,y in zip(a,b)]
        j=next((j for j,x in enumerate(a) if x),None)
        if j is not None:
            z=pow(a[j],-1,p);piv[j]=[x*z%p for x in a]
            if len(piv)==len(a):return len(piv)
    return len(piv)


def cycle(p):
    E=[]
    for i in range(p):
        v=[0]*p;v[i]=-1;v[(i+1)%p]=1
        E.append(tuple(v[1:])) # remove label-zero reference
    return E


def neg(a):return tuple(-x for x in a)
def dot(a,b):return sum(x*y for x,y in zip(a,b))


def sparse_moments(rows,t):
    C=Counter()
    for x in rows:
        S=support(x)
        for d in range(t+1):
            for I in combinations_with_replacement(S,d):
                C[I]+=prod(x[i] for i in I)
    return {I:v for I,v in C.items() if v}


def check_cycles():
    print('\nSIGNED-UNIT WEIGHT-FOUR MODULAR RANK TRADES')
    print('r | rows per tagged family | reduced dimension | deficient/full ranks | rational certificates')
    for p in (5,7,11,13,17,31,101):
        E=cycle(p);F=list(map(neg,E));d=p-1;u=tuple(range(1,p))
        assert len(set(E+F))==2*p
        assert sparse_moments(E,2)==sparse_moments(F,2)
        # Both orientations of each edge are distinct, cubic orientation witness.
        assert sparse_moments(E,3)!=sparse_moments(F,3)
        assert all(dot(x,u)%p==1 for x in E)
        assert all(dot(x,u)%p==p-1 for x in F)
        # A small explicit spanning subset of each family: fixed second or first
        # edge; all pairwise changes and the quotient direction are present.
        H=[a+F[0] for a in E]+[E[0]+b for b in F]+[F[0]+E[0]]
        J=[a+E[0] for a in E]+[E[0]+b for b in E]+[F[0]+F[0]]
        assert rank(H,p)==rank(J,p)==2*d-1
        assert rank(H+J,p)==2*d
        assert rank(H,1009)==rank(J,1009)==2*d
        z=u+u
        assert all(dot(a,z)%p==0 for a in H)
        assert all(dot(a,z)%p!=0 for a in J)
        # Distinct tagged rank-deficient vs full-rank systems. H,J are centrally
        # symmetric, so tag and constant directions can also be recovered.
        X=[a+(0,1) for a in H]+[a+(1,1) for a in H]
        Y=[a+(0,1) for a in H]+[a+(1,1) for a in J]
        assert rank(X,p)==2*d+1 and rank(Y,p)==2*d+2
        assert rank(X,1009)==rank(Y,1009)==2*d+2
        if p<=11:
            h=[a+b for a in E for b in F]+[a+b for a in F for b in E]
            j=[a+b for a in E for b in E]+[a+b for a in F for b in F]
            assert len(set(h+j))==4*p*p
            assert all(dot(x,z)%p==0 for x in h)
            assert all(dot(x,z)%p!=0 for x in j)
            x=[a+(0,1) for a in h]+[a+(1,1) for a in h]
            y=[a+(0,1) for a in h]+[a+(1,1) for a in j]
            assert len(set(x))==len(set(y))==4*p*p
            assert sparse_moments(x,5)==sparse_moments(y,5)
            assert sparse_moments(x,6)!=sparse_moments(y,6)
            assert rank(x,p)==2*d+1 and rank(y,p)==2*d+2
        print(f'{p} | {4*p*p} | {2*d+2} | {2*d+1}/{2*d+2} | mod 1009')
    print('7 rank fixtures; complete tagged moment tables checked for r=5,7,11.')


if __name__=='__main__':
    check_triangular();check_alternants();check_cycles()
    print('\nAll checks passed.')
