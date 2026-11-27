#!/usr/bin/env python3
"""Exact finite-field tests for abstract hidden determinantal geometry.
NumPy int64 is used only for the bounded small-field fixtures in this file.
"""
from itertools import combinations_with_replacement, product
from collections import Counter
from math import comb,isqrt
import numpy as np

if not __debug__:
    raise RuntimeError('Assertions required')


def rref(A,p):
    A=np.array(A,dtype=np.int64,copy=True)%p
    n,d=A.shape;k=0;P=[]
    for j in range(d):
        I=np.flatnonzero(A[k:,j])
        if not len(I):continue
        i=k+int(I[0]);A[[k,i]]=A[[i,k]]
        A[k]=A[k]*pow(int(A[k,j]),-1,p)%p
        c=A[:,j].copy();c[k]=0;I=np.flatnonzero(c)
        A[I]=(A[I]-c[I,None]*A[k])%p
        P.append(j);k+=1
        if k==n:break
    return A[:k],P


def basis(A,p):return rref(A,p)[0]
def rank(A,p):return len(basis(A,p))


def null(A,p):
    A=np.asarray(A,dtype=np.int64);d=A.shape[1]
    R,P=rref(A,p);free=[j for j in range(d) if j not in P]
    X=np.zeros((len(free),d),dtype=np.int64)
    for i,j in enumerate(free):
        X[i,j]=1;X[i,P]=-R[:,j]%p
    return X


def solve(A,b,p):
    A=np.asarray(A,dtype=np.int64);b=np.asarray(b,dtype=np.int64).reshape(-1,1)
    R,P=rref(np.hstack([A,b]),p);d=A.shape[1]
    if d in P:return None
    x=np.zeros(d,dtype=np.int64);x[P]=R[:,-1]
    assert np.all((A@x-b[:,0])%p==0)
    return x


def inv(A,p):
    d=len(A);R,P=rref(np.hstack([A,np.eye(d,dtype=np.int64)]),p)
    assert P==list(range(d))
    return R[:,d:]


def contains(U,x,p):return not np.any(null(U,p)@x%p)
def intersect(U,V,p):return null(np.vstack([null(U,p),null(V,p)]),p)


def extend(x,U,p):
    S=[np.asarray(x,dtype=np.int64)]
    for y in U:
        if rank(np.array(S+[y]),p)>len(S):S.append(y)
    return np.array(S)


def mons(d,t):return list(combinations_with_replacement(range(d),t))


def evaluate(X,E,p):
    X=np.asarray(X,dtype=np.int64);out=np.ones((len(X),len(E)),dtype=np.int64)
    for j,e in enumerate(E):
        for i in e:out[:,j]=out[:,j]*X[:,i]%p
    return out


def partial_tensor(F,d,t,p):
    E=mons(d,t);B=mons(d,t-1);index={e:i for i,e in enumerate(B)}
    R=np.zeros((len(F),len(B),d),dtype=np.int64)
    for j,e in enumerate(E):
        for i,c in Counter(e).items():
            z=list(e);z.remove(i);R[:,index[tuple(z)],i]+=c*F[:,j]
    return R%p


def learn(X,p,a,b,k):
    d=a*b;t=k+1;E=sum((mons(d,j) for j in range(t+1)),[])
    B=comb(a,t)*comb(b,t);I=null(evaluate(X,E,p),p)
    if len(I)!=B:return None
    h=len(mons(d,t));g=len(mons(d,t-1));F=I[:,-h:];G=I[:,-h-g:-h]
    if rank(F,p)!=B:return None
    T=partial_tensor(F,d,t,p).reshape(-1,d)
    c=solve(T,-G.ravel(),p)
    if c is None or rank(T,p)!=d:return None
    for j in range(t,2,-1):
        D=partial_tensor(F,d,j,p)
        F=basis(D.transpose(0,2,1).reshape(-1,D.shape[1]),p)
    if len(F)!=comb(a,2)*comb(b,2):return None
    return c,F,len(E)-B


def sym_forms(Q,d,p,dual=False):
    M=np.zeros((len(Q),d,d),dtype=np.int64);half=pow(2,-1,p)
    for j,(i,k) in enumerate(mons(d,2)):
        c=Q[:,j] if i==k or dual else half*Q[:,j]%p
        M[:,i,k]=c;M[:,k,i]=c
    return M


def zero(Q,x,p):return not np.any(((Q@x)@x)%p)
def tangent(Q,x,p):return null((Q@x)%p,p)


def stabilizer(Q,p,a,b):
    d=a*b;M=sym_forms(Q,d,p);H=sym_forms(null(Q,p),d,p,dual=True)
    # Tr(H(X^T M+M X))=2 <M H,X> for symmetric M,H.
    C=np.vstack([(M@h).reshape(-1,d*d)%p for h in H])
    g=null(C,p).reshape(-1,d,d)
    assert len(g)==a*a+b*b-1
    return g


def split_lie(g,p,a,b):
    d=a*b
    G,P=rref(np.array([(x@y-y@x).ravel()%p for i,x in enumerate(g) for y in g[i+1:]]),p)
    s=len(G);assert s==a*a+b*b-2;A=G.reshape(s,d,d)
    adj=[]
    for x in A:
        B=(x@A-A@x).reshape(s,d*d)%p
        coef=B[:,P];assert np.array_equal(coef@G%p,B);adj.append(coef.T)
    C=np.eye(s*s,dtype=np.int64);I=np.eye(s,dtype=np.int64)
    for ad in adj:
        eq=(np.kron(I,ad.T)-np.kron(ad,I))%p
        C=null(eq@C.T%p,p)@C%p
        if len(C)==2:break
    assert len(C)==2
    for ad in adj:
        for M in C.reshape(-1,s,s):assert not np.any((M@ad-ad@M)%p)
    L=next(M for M in C.reshape(-1,s,s) if np.any((M-M[0,0]*I)%p))
    parts=[]
    for z in range(p):
        K=null((L-z*I)%p,p)
        if len(K):parts.append((K@G%p).reshape(-1,d,d))
    assert sorted(map(len,parts))==sorted([a*a-1,b*b-1])
    assert all(not np.any((x@y-y@x)%p) for x in parts[0] for y in parts[1])
    return [np.concatenate([T,np.eye(d,dtype=np.int64)[None]]) for T in parts]


def orbit(A,x,p):return basis(A@x%p,p)


def projective_coeffs(q,d):
    for j in range(d):
        for tail in product(range(q),repeat=d-j-1):yield (0,)*j+(1,)+tail


def rank_one(Q,x,p,algs):
    A,B=algs;a=isqrt(len(A));b=isqrt(len(B));assert a*a==len(A) and b*b==len(B)
    checks=0;path=[]
    while True:
        U=orbit(A,x,p);V=orbit(B,x,p);h=len(U)//a
        assert h>=1 and len(U)==a*h and len(V)==b*h
        path.append(h)
        if h==1:
            assert zero(Q,x,p);return x,path,checks
        J=intersect(U,V,p);assert len(J)==h*h
        for coeff in projective_coeffs(p,h+1):
            y=np.array(coeff)@J[:h+1]%p;checks+=1
            if len(orbit(A,y,p))<len(U):x=y;break
        else:raise AssertionError('rank reduction failed')


def pure_parts(Q,y,W,p):
    if zero(Q,y,p):return [y]
    K=null(((Q@y)@W.T)%p,p);assert len(K)==1
    z=K[0]@W%p
    out=[(y+c*z)%p for c in range(1,p) if zero(Q,(y+c*z)%p,p)]
    assert len(out)==2 and all(np.any(x) for x in out)
    return out


def rulings(Q,x,p):
    assert np.any(x) and zero(Q,x,p)
    T=tangent(Q,x,p);W=extend(x,T,p)[1:]
    y=pure_parts(Q,W[0],W,p)[0];L=intersect(T,tangent(Q,y,p),p)
    v=next(v for v in W if not contains(L,v,p))
    z=next(z for z in pure_parts(Q,v,W,p) if not contains(L,z,p))
    R=intersect(T,tangent(Q,z,p),p)
    assert len(intersect(L,R,p))==1
    return L,R


def layout(Q,x,p,a,b):
    L,R=rulings(Q,x,p)
    if len(L)!=b:L,R=R,L
    assert (len(L),len(R))==(b,a)
    A=extend(x,L,p);B=extend(x,R,p)
    T=np.zeros((a,b,a*b),dtype=np.int64);T[0]=A;T[:,0]=B
    for i in range(1,a):
        V=next(U for U in rulings(Q,B[i],p) if len(intersect(U,R,p))==1)
        for j in range(1,b):
            U=next(U for U in rulings(Q,A[j],p) if len(intersect(U,L,p))==1)
            C=intersect(U,V,p);assert len(C)==1
            c=C[0];y=(x+A[j]+B[i])%p
            z=[s for s in range(1,p) if zero(Q,(y+s*c)%p,p)]
            assert len(z)==1;T[i,j]=z[0]*c%p
    T=T.reshape(a*b,a*b).T;assert rank(T,p)==a*b
    return T


def random_invertible(rng,d,p):
    while True:
        U=rng.integers(0,p,(d,d),dtype=np.int64)
        if rank(U,p)==d:return U


def random_rank(rng,a,b,k,p):
    while True:
        X=rng.integers(0,p,(a,k),dtype=np.int64)@rng.integers(0,p,(k,b),dtype=np.int64)%p
        if rank(X,p)==k:return X.ravel()


def hidden_fixture(rng,a,b,k,p,shifted):
    d=a*b;U=random_invertible(rng,d,p);Ui=inv(U,p)
    c=rng.integers(0,p,d,dtype=np.int64) if shifted else np.zeros(d,dtype=np.int64)
    t=k+1;E=sum((mons(d,j) for j in range(t+1)),[]);N=2*len(E)
    # Rank-exact samples test that no rank-one anchor is supplied for k>1.
    raw=np.array([random_rank(rng,a,b,k,p) for _ in range(N)]);X=(raw@U.T+c)%p
    out=learn(X,p,a,b,k);assert out is not None
    center,F,e=out;assert np.array_equal(center,c);Q=sym_forms(F,d,p)
    x=(X[0]-center)%p
    if k>1:
        assert not zero(Q,x,p)
        g=stabilizer(F,p,a,b);algs=split_lie(g,p,a,b)
        x,path,checks=rank_one(Q,x,p,algs)
        assert rank((Ui@x%p).reshape(a,b),p)==1 and path[0]==k and path[-1]==1
    else:path=[1];checks=0
    T=layout(Q,x,p,a,b);Ti=inv(T,p);R=(X-center)@Ti.T%p
    assert all(rank(y.reshape(a,b),p)==k for y in R)
    nchecks=0
    for j in range(min(a,b)+1):
        for _ in range(25):
            y=random_rank(rng,a,b,j,p) if j else np.zeros(d,dtype=np.int64)
            h=Ti@U@y%p;assert rank(h.reshape(a,b),p)==j;nchecks+=1
    z=T@np.eye(a,b,dtype=np.int64).ravel()%p
    assert rank((Ui@z%p).reshape(a,b),p)==min(a,b)>k
    B=np.array([(U@random_rank(rng,a,b,1,p)+c)%p for _ in range(12)])
    bmax=max(rank((Ti@(y-c)%p).reshape(a,b),p) for y in B);budget=k+bmax
    if budget<min(a,b):
        target=(2*c+z)%p
        assert rank((Ui@(target-2*c)%p).reshape(a,b),p)>budget
        for v in X[:12]:
            for w in B:assert np.any((v+w-target)%p)
    bad=(c+z)%p;assert rank((Ti@(bad-c)%p).reshape(a,b),p)>k
    assert learn(X[:1],p,a,b,k) is None
    print(f'{a}x{b} k={k} p={p} affine={shifted} samples={N} eval-rank={e} rank-path={path} search={checks} rank-tests={nchecks} two-list={budget<min(a,b)}',flush=True)
    return N,nchecks,checks


def main():
    rng=np.random.default_rng(20261127);count=tests=search=0
    print('HIDDEN DETERMINANTAL GEOMETRY FROM POLYNOMIALLY COMPLETE SAMPLES',flush=True)
    cases=[(2,3,1,7,False),(3,3,1,7,True),(3,4,1,11,True),(3,3,2,13,False),(3,3,2,13,True),(3,4,2,13,True),(4,4,2,13,True)]
    for c in cases:
        x,y,z=hidden_fixture(rng,*c);count+=x;tests+=y;search+=z
    print(f'{len(cases)} recovered hidden layouts; {count} training vectors; {tests} independent exact rank checks; {search} projective rank-reduction candidates.')
    print('All affine centers exact; insufficient-sample and outlier-cap controls checked.')
    print('All checks passed.')

if __name__=='__main__':main()
