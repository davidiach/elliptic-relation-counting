#!/usr/bin/env python3
"""Exact checks for abstract moment/rank boundaries. Standard library only."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, prod
import random

if not __debug__:
    raise RuntimeError("Run the verifier with assertions enabled.")

def rank(rows, p):
    piv = {}
    for row in rows:
        a = [x % p for x in row]
        for j, b in sorted(piv.items()):
            if a[j]:
                z = a[j]
                a = [(x-z*y) % p for x,y in zip(a,b)]
        j = next((j for j,x in enumerate(a) if x), None)
        if j is not None:
            z = pow(a[j],-1,p)
            piv[j] = [(x*z)%p for x in a]
            if len(piv)==len(a): return len(piv)
    return len(piv)


def bound(q,t):
    if t<=0:return 1
    a,b=divmod(t,q-1)
    return q**a*(b+1)


def monomials(d,t):
    return [a for a in product(range(t+1),repeat=d) if sum(a)<=t]


def mon(x,e):return prod(v**j for v,j in zip(x,e))


def check_distance():
    print('CARTESIAN DUAL SUPPORT BOUNDS')
    print('alphabet | dimension | degree+1 | bound | sub-bound supports | modulus')
    cases=[(2,3,t,5,True) for t in range(1,4)]
    cases += [(3,2,t,7,True) for t in range(1,5)]
    cases += [(5,2,t,11,True) for t in range(1,4)]
    cases += [(3,3,t,7,False) for t in range(3,7)]
    cases += [(5,2,t,11,False) for t in range(4,9)]
    rng=random.Random(20261127); count=0
    for q,d,t,p,exhaust in cases:
        X=list(product(range(q),repeat=d)); es=monomials(d,t-1)
        E=[tuple(mon(x,e)%p for e in es) for x in X]
        B=bound(q,t); n=0
        if exhaust:
            tests=(I for k in range(1,B) for I in combinations(range(len(X)),k))
        else:
            tests=(rng.sample(range(len(X)),rng.randrange(1,B)) for _ in range(100))
        for I in tests:
            assert rank([E[i] for i in I],p)==len(I)
            n+=1
        # Product finite-difference vector at the exact support bound.
        a,b=divmod(t,q-1); sizes=[q]*a+([b+1] if b else [])
        S=[(tuple(x)+(0,)*(d-len(sizes)),prod((-1)**j*comb(k-1,j) for k,j in zip(sizes,x)))
           for x in product(*(range(k) for k in sizes))]
        assert len(sizes)<=d and len(S)==B and all(w%p for _,w in S)
        assert all(sum(w*mon(x,e) for x,w in S)%p==0 for e in es)
        print(f'{q} | {d} | {t} | {B} | {n} | {p}')
        count+=n
    print(f'{count} sub-bound support rank checks; {len(cases)} sharp unrestricted-weight examples.')


def signs(alpha):
    S=[i for i,v in enumerate(alpha) if v]
    for e in product((-1,1),repeat=len(S)):
        x=list(alpha)
        for i,s in zip(S,e):x[i]*=s
        yield tuple(x)


def add(P,Q,c=F(1)):
    R=dict(P)
    for e,a in Q.items():
        R[e]=R.get(e,F(0))+c*a
        if not R[e]:del R[e]
    return R


def mul(P,Q):
    R={}
    for e,a in P.items():
        for f,b in Q.items():
            g=tuple(x+y for x,y in zip(e,f))
            R[g]=R.get(g,F(0))+a*b
    return {e:a for e,a in R.items() if a}


def unit(d,i,k):return tuple(k if j==i else 0 for j in range(d))
def evaluate(P,x):return sum(a*mon(x,e) for e,a in P.items())


@lru_cache(None)
def interpolation(d,W):
    B=[x for x in product(range(W+1),repeat=d) if sum(x)<=W];zero=(0,)*d
    @lru_cache(None)
    def indicator(alpha,R):
        S=[i for i,v in enumerate(alpha) if v]
        assert R and len(R)%2==0 and set(R)<=set(S)
        P={zero:F(1)}
        for i in S:
            P=mul(P,{unit(d,i,1 if i in R else 2):F(1)})
            for j in range(1,alpha[i]):
                P=mul(P,{unit(d,i,2):F(1),zero:F(-j*j)})
        c=evaluate(P,alpha);assert c
        Q=dict(P)
        for beta in B:
            if beta!=alpha and all(y>=x for x,y in zip(alpha,beta)):
                Q=add(Q,indicator(beta,R),-evaluate(P,beta))
        Q={e:a/c for e,a in Q.items()}
        assert all(sum(e)<=2*W-2 for e in Q)
        return Q
    return B,indicator


def check_certificates():
    print('\nNONNEGATIVE SPARSE RANK CERTIFICATES')
    print('d | W | prime | directions | grid size | max degree | max terms')
    cases=[(2,2,5),(2,3,5),(2,3,7),(2,4,7),(3,3,7),(3,4,11),(3,5,11),(4,3,11),(4,4,13)]
    rng=random.Random(727); count=0; points=0; indicators=0
    for d,W,p in cases:
        B,I=interpolation(d,W)
        X=[x for a in B for x in signs(a)]
        for a in B:
            S=[i for i,v in enumerate(a) if v]
            for k in range(2,len(S)+1,2):
                for R in combinations(S,k):
                    P=I(a,R)
                    for x in X:
                        expected=prod(1 if x[i]>0 else -1 for i in R) if tuple(map(abs,x))==a else 0
                        assert evaluate(P,x)==expected
                    indicators+=1
        if d==2: dirs=list(product(range(p),repeat=d))[1:]
        else:
            dirs=[tuple(rng.randrange(p) for _ in range(d)) for _ in range(8)]
            dirs += [tuple(int(j==i) for j in range(d)) for i in range(d)]
        max_degree=0;max_terms=0
        for z in dirs:
            if not any(z):continue
            P={unit(d,i,2):F(1) for i in range(d) if z[i]};targets={}
            for a in B:
                base=sum(v*v for v,w in zip(a,z) if w)
                S=[i for i,v in enumerate(a) if v]
                A=list(signs(a));bad=[x for x in A if sum(v*w for v,w in zip(x,z))%p]
                assert bool(bad)==bool(base)
                if not bad:
                    targets.update((x,F(0)) for x in A);continue
                assert 2*len(bad)>=len(A)
                badset=set(bad);g=F(base*len(A),len(bad))
                targets.update((x,g if x in badset else F(0)) for x in A)
                for k in range(2,len(S)+1,2):
                    for R in combinations(S,k):
                        c=F(base,len(bad))*sum(prod(1 if x[i]>0 else -1 for i in R) for x in bad)
                        if c:P=add(P,I(a,R),c)
            assert all(sum(e)<=2*W-2 for e in P)
            for x in X:
                f=evaluate(P,x)
                assert f==targets[x] and f>=0
                assert (f==0)==(sum(v*w for v,w in zip(x,z))%p==0)
            max_degree=max(max_degree,max(map(sum,P),default=0));max_terms=max(max_terms,len(P))
            count+=1;points+=len(X)
        print(f'{d} | {W} | {p} | {len(dirs)} | {len(X)} | {max_degree} | {max_terms}')
    print(f'{indicators} exact sign-character interpolation checks on complete balls.')
    print(f'{count} rational polynomial certificates; {points} exact point evaluations.')


def gram(A):
    T=list(zip(*A));return tuple(tuple(sum(x*y for x,y in zip(a,b)) for b in T) for a in T)


def moments(A,t):
    out=Counter()
    for x in A:
        supp=[i for i,v in enumerate(x) if v]
        for k in range(t+1):
            for I in combinations_with_replacement(supp,k):
                out[I]+=prod(x[i] for i in I)
    return {I:v for I,v in out.items() if v}


def det(A):
    A=[list(x) for x in A]; n=len(A); sign=1;last=1
    for k in range(n-1):
        i=next((i for i in range(k,n) if A[i][k]),None)
        if i is None:return 0
        if i!=k:A[i],A[k]=A[k],A[i];sign=-sign
        z=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                a,b=divmod(A[i][j]*z-A[i][k]*A[k][j],last)
                assert b==0;A[i][j]=a
        for i in range(k+1,n):A[i][k]=0
        last=z
    return sign*A[-1][-1]


def null_one(A,p):
    A=[[v%p for v in a] for a in A];n=len(A[0]);k=0;piv=[]
    for j in range(n):
        i=next((i for i in range(k,len(A)) if A[i][j]),None)
        if i is None:continue
        A[k],A[i]=A[i],A[k];z=pow(A[k][j],-1,p);A[k]=[(v*z)%p for v in A[k]]
        for i in range(len(A)):
            if i!=k and A[i][j]:
                z=A[i][j];A[i]=[(v-z*w)%p for v,w in zip(A[i],A[k])]
        k+=1;piv.append(j)
        if k==len(A):break
    free=[j for j in range(n) if j not in piv];assert len(free)==1
    x=[0]*n;x[free[0]]=1
    for i,j in enumerate(piv):x[j]=-A[i][free[0]]%p
    return tuple(x)


def fan(k):
    U=[];V=[]
    for i in range(k):
        a=[0]*k;a[i]=1; b=[0]*(k-1)
        if i<k-1:b[i]=1
        if i:b[i-1]=-1
        U.append(a);V.append(b)
    for i in range(k-1):
        a=[0]*k;a[i]=-1;a[i+1]=1;b=[0]*(k-1);b[i]=1
        U.append(a);V.append(b)
    return U,V


def check_fans():
    print('\nTHREE-ACTIVE-COORDINATE CUBIC COUNTEREXAMPLES')
    print('fan k | prime | spanning trees | rows | dimension | ranks | max active norm_1 | distinct violations')
    cases=[(11,89),(11,199),(13,233),(13,521),(17,1597),(17,3571)]
    for k,p in cases:
        assert all(p%i for i in range(2,int(p**.5)+1))
        U,V=fan(k);n=2*k-1
        tau=det(gram(U));assert tau%p==0
        f,g=0,1
        for _ in range(2*k):f,g=g,f+g
        assert tau==f
        assert all(sum(v*w for v,w in zip(a,b))==0 for a in zip(*U) for b in zip(*V))
        A=[tuple(a+b) for a,b in zip(U,V)]; B=[tuple(a+[-x for x in b]) for a,b in zip(U,V)]
        assert abs(det(A))==abs(det(B))==tau
        assert rank(A,p)==rank(B,p)==n-1 and rank(A+B,p)==n
        z=null_one(A,p)+(0,0)
        H=[row+(1,) for a in A for row in (a,tuple(-v for v in a))]
        J=[row+(1,) for a in B for row in (a,tuple(-v for v in a))]
        X=sorted(a+(0,) for a in H)+sorted(a+(1,) for a in H)
        Y=sorted(a+(0,) for a in H)+sorted(a+(1,) for a in J)
        assert len(X)==len(Y)==4*n<p and len(set(X))==len(X) and len(set(Y))==len(Y)
        assert moments(X,3)==moments(Y,3) and moments(X,4)!=moments(Y,4)
        assert rank(X,p)==n+1 and rank(Y,p)==n+2
        assert rank(X,65537)==rank(Y,65537)==n+2
        active=max(sum(abs(v) for v,w in zip(a,z) if w) for a in X+Y)
        assert active<=3
        assert all(sum(v*w for v,w in zip(a,z))%p==0 for a in X)
        bad=sum(sum(v*w for v,w in zip(a,z))%p!=0 for a in Y)
        assert bad>0
        print(f'{k} | {p} | {tau} | {len(X)} | {n+2} | {n+1}/{n+2} | {active} | {bad}')
    print(f'{len(cases)} exact full-rational-rank, no-wrap, distinct-row fan examples.')

if __name__=='__main__':
    check_distance();check_certificates();check_fans();print('\nAll checks passed.')
