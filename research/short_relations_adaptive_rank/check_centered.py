"""Exact small-field tests of the centered collision bound and pin-rank estimates."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product, combinations
from math import prod
from random import Random
import numpy as np
import sympy as sp


def rank(rows,r,n):
    B={}
    for row in rows:
        u=[int(x)%r for x in row]
        for j in sorted(B):
            if u[j]:u=[(a-u[j]*b)%r for a,b in zip(u,B[j])]
        j=next((i for i,a in enumerate(u) if a),None)
        if j is not None:B[j]=[a*pow(u[j],-1,r)%r for a in u]
    return len(B)


def partitions(h):
    def go(a):
        if len(a)==h:yield tuple(a);return
        for j in range(max(a)+2):yield from go(a+[j])
    return list(go([0]))


def pat(t):
    d={};return tuple(d.setdefault(x,len(d)) for x in t)


def isolated(ps,i):
    return all(p.count(p[i])==1 for p in ps)


def incidence(ps,h):
    return [[int(p[i]==j) for p in ps for j in range(max(p)+1)] for i in range(h)]


def partition_checks():
    tested=nonisolated=0
    for h,m in [(2,2),(3,3),(4,2),(4,3),(5,2)]:
        for ps in product(partitions(h),repeat=m):
            if any(isolated(ps,i) for i in range(h)):continue
            A=incidence(ps,h);ell=len(A[0]);C=m*h-ell
            for r in (3,5,19):
                R=rank(A,r,ell);D=h-R
                s=sum(rank(A[:i]+A[i+1:],r,ell)==R for i in range(h))
                assert 2*C>=m*max(2*D,s)+(h-s)
                assert C>=m*D
                if 4*D<h:assert 4*(C-m*D)>=h
                tested+=1
            nonisolated+=1
    print(f'Partition bound: {nonisolated} nonisolated patterns, {tested} modular rank/support checks passed.')


def expectation_checks():
    # Independent Bernoulli half-retention, all seeds enumerated (stronger than h-wise independence).
    cases=0;patterns=0;seeds=0
    for r,m,h in [(5,2,2),(5,2,4),(5,3,4)]:
        T=[[] for _ in range(r)]
        for x in product(range(1,r),repeat=m):T[sum(x)%r].append(x)
        moments=[Fraction(0) for _ in range(r)];theta=Fraction(1,2)
        mu=[len(t)*theta**m for t in T]
        total=1<<(m*(r-1));seeds+=total
        for mask in range(total):
            d=[sum(all((mask>>(i*(r-1)+a-1))&1 for i,a in enumerate(x)) for x in tt) for tt in T]
            for g in range(r):moments[g]+=(d[g]-mu[g])**h/total
        qp=r*theta;mup=qp**m/r
        for g in range(r):
            pc=Counter(tuple(pat(col) for col in zip(*rows)) for rows in product(T[g],repeat=h))
            answer=Fraction(0);upper=Fraction(0)
            for ps,count in pc.items():
                e=0
                for bits in product((0,1),repeat=h):
                    keep=[i for i,j in enumerate(bits) if j];q=sum(len({p[i] for i in keep}) for p in ps)
                    e+=(-1)**(h-len(keep))*theta**(m*(h-len(keep))+q)
                A=incidence(ps,h);ell=len(A[0]);R=rank(A,r,ell)
                assert count<=r**(ell-R)
                assert abs(e)<=2**h*theta**ell
                if any(isolated(ps,i) for i in range(h)):assert not e
                else:
                    upper+=2**h*qp**ell/r**R
                    assert 4*(h-R)>=h or 4*(m*h-ell-m*(h-R))>=h
                answer+=count*e;patterns+=1
            assert answer==moments[g] and 0<=answer<=upper
            cases+=1
    print(f'Exact centered moment expansion: {seeds} independent-mask seeds, {cases} target checks, {patterns} realized partition patterns.')


def fixture_tables(blocks,r):
    vals=[a for B in blocks for a in B];n=len(vals);ids=[];off=0
    for B in blocks:ids.append(range(off,off+len(B)));off+=len(B)
    T=[[] for _ in range(r)]
    for X in product(*ids):T[sum(vals[i] for i in X)%r].append(X)
    return vals,ids,T


def rows_at(T,t,r,n):
    rows=[];H=[]
    for g,XX in enumerate(T):
        for X in XX:
            for Y in T[(g+t)%r]:
                row=[0]*n;pin=[0]*(2*n)
                for j in X:row[j]-=1;pin[j]=1
                for j in Y:row[j]+=1;pin[n+j]=1
                rows.append(row);H.append(pin)
    return np.array(rows,dtype=np.int64),np.array(H,dtype=np.int64)


def valuation(n,r):
    assert n>0;v=0
    while n%r==0:n//=r;v+=1
    return v


def kernel(rows,r,n):
    A=np.array(rows,dtype=np.int64).reshape(-1,n)%r;k=0;piv=[]
    for j in range(n):
        ix=np.flatnonzero(A[k:,j]);
        if not len(ix):continue
        i=k+int(ix[0]);A[[k,i]]=A[[i,k]];A[k]=A[k]*pow(int(A[k,j]),-1,r)%r
        for i in range(len(A)):
            if i!=k:A[i]=(A[i]-int(A[i,j])*A[k])%r
        piv.append(j);k+=1
    U=[]
    for j in range(n):
        if j in piv:continue
        u=np.zeros(n,dtype=np.int64);u[j]=1
        for i,a in enumerate(piv):u[a]=-A[i,j]%r
        U.append(u)
    return U


def adaptive(A,r):
    G=A.T@A;W=[];n=A.shape[1]
    while True:
        U=kernel(W,r,n)
        if not U:return len(W)
        chosen=None
        for u in U:
            g=G@u;j=next((j for j,v in enumerate(g) if v%r),None)
            if j is not None:chosen=(A[:,j]*(A@u),r,u);break
        if chosen is None:
            for i,u in enumerate(U):
                for v in [u]+[u+w for w in U[i+1:]]:
                    if int(v@G@v)%(r*r):chosen=((A@v)**2,r*r,v);break
                if chosen is not None:break
        if chosen is None:return len(W)
        weights,mod,u=chosen
        lo=0;hi=len(A);assert int(sum(weights))%mod
        while hi-lo>1:
            mid=(lo+hi)//2
            if int(sum(weights[lo:mid]))%mod:hi=mid
            else:lo=mid
        assert int(A[lo]@u)%r
        W.append(A[lo])


def rank_checks():
    rng=Random(2026);passed=singular=rescues=0;maxdef=0;example=None
    for r,m,k in [(7,2,3),(19,3,3),(31,3,4),(43,3,4),(61,3,5)]:
        for _ in range(500):
            B=[sorted(rng.sample(range(1,r),k)) for i in range(m)]
            vals,ids,T=fixture_tables(B,r);n=len(vals)
            if not all(T):continue
            A,H=rows_at(T,1,r,n);mass=H.sum(axis=0)
            if min(mass)<=0 or max(mass)>=r:continue
            keep=[j for I in ids for j in list(I)[1:]];A=A[:,keep];D=len(keep)
            assert rank(A,r,D)==D
            V=H.T@A;R=mass[0]*V-mass[:,None]*V[0];rr=rank(R,r,D)
            G=A.T@A;det=int(sp.Matrix(G.tolist()).det());vv=valuation(det,r)
            assert det<=prod(int(x) for x in G.diagonal())
            assert rr>=rank(G,r,D)-1
            got=adaptive(A,r);assert 2*(D-got)<=vv
            singular+=vv>0;rescues+=(vv>0 and got==D);maxdef=max(maxdef,D-1-rr)
            if example is None and vv>0 and got==D:example=(r,B,D,vv,got,int(max(mass)))
            passed+=1
            if passed>=150:break
        if passed>=150:break
    assert passed>0
    print(f'Full-coverage/no-wrap rank: {passed} fixed-shift fixtures, {singular} singular Gram matrices, {rescues} rescued, maximum contrast pin-average rank loss {maxdef}.')
    print('First singular-Gram rescue:',example)


def curve_fixture():
    p,a,b=5,2,1;r=7
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        z=((3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
        X=(z*z-x-u)%p;return X,(z*(x-X)-y)%p
    pts=[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-b)%p==0]
    assert len(pts)+1==r and (4*a**3+27*b*b)%p
    P=(0,1);G=[None]
    for j in range(1,r):G.append(add(G[-1],P))
    assert set(G)=={None,*pts} and add(G[-1],P) is None
    blocks=[[4,5,6],[2,4,6]];vals,ids,T=fixture_tables(blocks,r)
    for g,tt in enumerate(T):
        for X in tt:
            z=None
            for j in X:z=add(z,G[vals[j]])
            assert z==G[g]
    A,H=rows_at(T,1,r,len(vals));mass=H.sum(axis=0)
    for row in A:
        z=None
        for j,c in enumerate(row):
            if c:z=add(z,G[(int(c)*vals[j])%r])
        assert z==P
    keep=[j for I in ids for j in list(I)[1:]];A=A[:,keep];D=len(keep)
    Gram=A.T@A;vv=valuation(int(sp.Matrix(Gram.tolist()).det()),r)
    V=H.T@A;R=mass[0]*V-mass[:,None]*V[0]
    assert all(T) and min(mass)>0 and max(mass)<r
    assert rank(A,r,D)==D and rank(R,r,D)==D-1 and vv==1 and adaptive(A,r)==D
    print(f'Curve E/F5, y^2=x^3+2x+1, order 7: {len(A)} actual group-law relations, fiber counts {min(map(len,T))}..{max(map(len(T)))}, pin masses {min(mass)}..{max(mass)}, Gram rank {rank(Gram,r,D)}/{D}, adaptive rank {D}/{D}.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    partition_checks();expectation_checks();rank_checks();curve_fixture()
