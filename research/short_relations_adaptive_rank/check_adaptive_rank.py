"""Finite exact tests of adaptive r-adic rank extraction. No tensor contractor."""
from itertools import product
from random import Random
from math import prod


def rref(rows,r,n):
    a=[list(map(lambda v:v%r,row)) for row in rows]; j=0; piv=[]
    for i in range(n):
        k=next((k for k in range(j,len(a)) if a[k][i]),None)
        if k is None:continue
        a[j],a[k]=a[k],a[j]; c=pow(a[j][i],-1,r);a[j]=[x*c%r for x in a[j]]
        for k in range(len(a)):
            if k!=j and a[k][i]:
                c=a[k][i];a[k]=[(x-c*y)%r for x,y in zip(a[k],a[j])]
        piv.append(i);j+=1
    return a[:j],piv


def kernel(rows,r,n):
    a,piv=rref(rows,r,n);U=[]
    for i in range(n):
        if i in piv:continue
        u=[0]*n;u[i]=1
        for row,j in zip(a,piv):u[j]=(-row[i])%r
        U.append(u)
    return U


def dot(a,b):return sum(x*y for x,y in zip(a,b))

def gram(rows,n):return [[sum(a[i]*a[j] for a in rows) for j in range(n)] for i in range(n)]

def extract(rows,weight,mod):
    assert sum(map(weight,rows))%mod
    lo,hi=0,len(rows)
    while hi-lo>1:
        mid=(lo+hi)//2
        if sum(map(weight,rows[lo:mid]))%mod:hi=mid
        else:lo=mid
    assert weight(rows[lo])%mod
    return rows[lo]


def adaptive(families,r,n,initial=()):
    W=list(initial);count=0;Gs=[gram(F,n) for F in families]
    while True:
        U=kernel(W,r,n)
        if not U:return W,count
        selected=None
        for F,G in zip(families,Gs):
            for u in U:
                a=[dot(row,u)%r for row in G]
                j=next((j for j,v in enumerate(a) if v),None)
                if j is not None:
                    selected=extract(F,lambda row:row[j]*dot(row,u),r);break
            if selected is not None:break
            for i,u in enumerate(U):
                v=next((v for v in [u]+[[x+y for x,y in zip(u,z)] for z in U[i+1:]] if dot(v,[dot(row,v) for row in G])%(r*r)),None)
                if v is not None:
                    selected=extract(F,lambda row:dot(row,v)**2,r*r);break
            if selected is not None:break
            assert all(dot(u,[dot(row,v) for row in G])%(r*r)==0 for u in U for v in U)
        if selected is None:return W,count
        assert any(dot(selected,u)%r for u in U)
        W.append(selected);count+=1


def fixture(r,blocks,target):
    vals=[v for B in blocks for v in B]+[1,target];n=len(vals);ix=[];j=0
    for B in blocks:
        ix.extend([range(j,j+len(B))]*2);j+=len(B)
    T=[];by={}
    for x in product(*ix):
        t=sum(vals[i] for i in x)%r;a=[x.count(i) for i in range(n)];T.append((t,a));by.setdefault(t,a)
    if len(by)<r:return None
    E=[[int(i==n-2) for i in range(n)]]+[[int(i==j)-a for i,a in enumerate(by[v])] for j,v in enumerate(vals)]
    S=[[y-x for x,y in zip(a,b)] for t,a in T for s,b in T if s==(t+1)%r]
    return E,S,n


def main():
    import sympy as sp
    r=7;E,S,n=fixture(r,[(1,4,5),(1,4,5)],6)
    d=int(sp.Matrix(gram(E+S,n)).det());assert d%(r*r)==0
    W,k=adaptive([E+S],r,n)
    V,l=adaptive([E,S],r,n,E)
    print(f'C7 determinant mass fails modulo 49; adaptive combined-family rank {n-len(kernel(W,r,n))}/{n}; partitioned rank {n-len(kernel(V,r,n))}/{n}.')
    rng=Random(847);cases=0;rescued=0;stuck=0
    for r in [3,5,7,11]:
        for _ in range(100):
            n=3;A=[[rng.randrange(-3,4) for _ in range(n)] for i in range(rng.randrange(3,18))]
            if kernel(A,r,n):continue
            G=sp.Matrix(gram(A,n));d=int(G.det());v=0
            while d%r==0:d//=r;v+=1
            W,k=adaptive([A],r,n);D=len(kernel(W,r,n));assert 2*D<=v
            if v>=2 and D==0:rescued+=1
            if D:stuck+=1
            cases+=1
    # All rows identical repeated r^2 times: every quadratic test fails despite rank one.
    W,k=adaptive([[[1]]*25],5,1);assert not W
    print(f'Adaptive exact-rank tests: {cases} full-rank arrays, {rescued} valuation-at-least-two rescues, {stuck} stalled; rank-deficiency bound and repetition obstruction passed.')

if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    main()
