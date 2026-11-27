#!/usr/bin/env python3
"""Finite exact checks for support-complex recognition and CRT jet certificates."""
from itertools import product, combinations
from math import comb
from random import Random
import numpy as np

if not __debug__:
    raise RuntimeError('Assertions must be enabled')

rng=Random(271126)

def rref(A,p,d=None):
    A=np.asarray(A,dtype=np.int64)
    if A.size==0:
        return np.zeros((0,d if d is not None else A.shape[-1]),dtype=np.int64),[]
    A=A.copy()%p;k=0;P=[]
    for j in range(A.shape[1]):
        ids=np.flatnonzero(A[k:,j])
        if not len(ids):continue
        i=k+int(ids[0]);A[[k,i]]=A[[i,k]]
        A[k]=A[k]*pow(int(A[k,j]),-1,p)%p
        c=A[:,j].copy();c[k]=0;A=(A-c[:,None]*A[k])%p
        P.append(j);k+=1
        if k==len(A):break
    return A[:k],P

def basis(A,p,d):return rref(A,p,d)[0]
def key(A):return tuple(map(tuple,A.tolist()))
def null(A,p,d):
    R,P=rref(A,p,d);F=[i for i in range(d) if i not in P];N=np.zeros((len(F),d),dtype=np.int64)
    for k,j in enumerate(F):
        N[k,j]=1
        if P:N[k,P]=-R[:,j]%p
    return N

def inside(x,B,p):
    y=np.asarray(x,dtype=np.int64).copy()%p
    for b in B:
        j=int(np.flatnonzero(b)[0]);y=(y-y[j]*b)%p
    return not np.any(y)

def contains(A,B,p):return all(inside(x,B,p) for x in A)
def intersect(A,B,p,d):return null(np.vstack([null(A,p,d),null(B,p,d)]),p,d)

def rndmat(n,p):
    while True:
        A=np.array([[rng.randrange(p) for _ in range(n)] for i in range(n)],dtype=np.int64)
        if len(basis(A,p,n))==n:return A

def inv(A,p):
    n=len(A);R,P=rref(np.hstack([A,np.eye(n,dtype=np.int64)]),p)
    assert P==list(range(n));return R[:,n:]

def vs(B,p):
    d=B.shape[1];return {tuple(np.asarray(c,dtype=np.int64)@B%p) for c in product(range(p),repeat=len(B))} if len(B) else {(0,)*d}

def add(x,y,p):return tuple((a+b)%p for a,b in zip(x,y))
def sub(x,y,p):return tuple((a-b)%p for a,b in zip(x,y))

def flat_key(x,B,p):
    x=np.array(x,dtype=np.int64).copy();B=basis(B,p,len(x))
    for b in B:
        j=int(np.flatnonzero(b)[0]);x=(x-x[j]*b)%p
    return tuple(x),key(B)

def recover_complete(S,p):
    # Every candidate is explicitly checked for full affine-flat containment.
    S=set(S);ordered=sorted(S);covered=set();facets={};attempts=0
    for x in ordered:
        if x in covered:continue
        P={x};directions=[]
        for y in ordered:
            if y in P:continue
            z=sub(y,x,p);attempts+=1;new=set(P);ok=True
            for c in range(1,p):
                cz=tuple(c*a%p for a in z)
                for v in P:
                    h=add(v,cz,p)
                    if h not in S:ok=False;break
                    new.add(h)
                if not ok:break
            if ok:P=new;directions.append(z)
        B=basis(directions,p,len(x));assert len(P)==p**len(B)
        facets[flat_key(x,B,p)]=P;covered|=P
    assert covered==S
    return facets,attempts

def layers(facets,p,d):
    # Affine center may not be unique. Any common point works.
    C=set.intersection(*map(set,facets.values()));assert C
    center=min(C);F=[np.array(k[1],dtype=np.int64).reshape(-1,d) for k in facets]
    I={key(B):B for B in F};todo=list(I.values())
    for A in todo:
        for B in F:
            Z=basis(intersect(A,B,p,d),p,d)
            if len(Z) and key(Z) not in I:I[key(Z)]=Z;todo.append(Z)
    records=[]
    for A in sorted(I.values(),key=lambda B:(len(B),key(B))):
        smaller=[B for B in I.values() if len(B)<len(A) and contains(B,A,p)]
        Z=basis(np.vstack(smaller) if smaller else [],p,d);start=len(Z);chosen=[]
        for v in A:
            if not inside(v,Z,p):chosen.append(v);Z=basis(np.vstack([Z,v]),p,d)
        assert len(Z)==len(A) and len(chosen)==len(A)-start
        if chosen:records.append((A,basis(chosen,p,d)))
    # Verify direct sums and actual decompositions in every facet.
    local=[]
    for A in F:
        ids=[i for i,(I,B) in enumerate(records) if contains(I,A,p)]
        Z=np.vstack([records[i][1] for i in ids]);assert len(Z)==len(A)==len(basis(Z,p,d))
        assert contains(Z,A,p);local.append((ids,Z))
    n=0
    for x in set.union(*facets.values()):
        x0=sub(x,center,p);answers=[]
        for ids,Z in local:
            if not inside(x0,basis(Z,p,d),p):continue
            # Full-row-rank Z: use a pivot column information set.
            J=rref(Z,p)[1];coeff=np.asarray(x0)[J]@inv(Z[:,J],p)%p
            assert tuple(coeff@Z%p)==x0
            out={};k=0
            for i in ids:
                B=records[i][1];v=tuple(coeff[k:k+len(B)]@B%p);k+=len(B)
                if any(v):out[i]=v
            answers.append(out)
        assert answers and all(a==answers[0] for a in answers);n+=1
    return records,n,len(I)

def fixture(dims,F,p,projection='none'):
    D=sum(dims);off=np.cumsum([0]+dims);I=np.eye(D,dtype=np.int64)
    if projection=='one':
        A=np.vstack([np.eye(D-1,dtype=np.int64),-np.ones(D-1,dtype=np.int64)])%p
    elif projection=='vander':
        assert all(a==1 for a in dims);A=np.array([[pow(i,j,p) for j in range(8)] for i in range(D)],dtype=np.int64)
    else:A=I
    d=A.shape[1];A=A@rndmat(d,p)%p;c=tuple(rng.randrange(p) for _ in range(d))
    groups=[A[off[i]:off[i+1]] for i in range(len(dims))]
    spaces=[basis(np.vstack([groups[i] for i in J]),p,d) for J in F]
    sets=[{add(c,x,p) for x in vs(B,p)} for B in spaces]
    return c,A,groups,spaces,sets

def check_complete():
    print('COMPLETE ARBITRARY SUPPORT COMPLEXES')
    cases=[(5,[1,1,1,1],[(0,1),(1,2),(3,)],'none'),
           (5,[1,2,1,1],[(0,1),(1,2),(0,3)],'none'),
           (5,[1]*9,[(i,(i+1)%9) for i in range(9)],'one'),
           (13,[1]*13,[(i,(i+1)%13) for i in range(13)],'vander'),
           (7,[1]*5,[(0,1,2),(0,3),(4,)],'none'),
           (5,[1]*5,[(0,1),(2,3,4)],'none')]
    total=pts=0
    for p,dims,F,proj in cases:
        c,A,G,B,sets=fixture(dims,F,p,proj);S=set.union(*sets);d=A.shape[1];t=max(map(len,F));assert p>t+1
        recovered,attempts=recover_complete(S,p)
        assert set(recovered)=={flat_key(c,V,p) for V in B}
        R,n,I=layers(recovered,p,d);sig={}
        for i,a in enumerate(dims):
            J=tuple(j for j,f in enumerate(F) if i in f);sig[J]=sig.get(J,0)+a
        assert sorted(len(B) for _,B in R)==sorted(sig.values())
        total+=len(recovered);pts+=n
        print(f'q={p} blocks={len(dims)} rank={d} points={len(S)} facets={len(F)} layers={len(R)} intersections={I} extension attempts={attempts}')
    print(f'{total} exact facet recoveries; {pts} point decompositions consistent across every incident facet.')
    # Sharp field-size failures: affine q=t+1; linear q=t.
    for p in (2,3,5,7):
        affine=[tuple((a-j)%p for j in range(p)) for a in range(p)]
        assert all(sum(x!=0 for x in v)==p-1 for v in affine)
        assert all(any(v[i]!=0 for v in affine) for i in range(p))
        linear=[tuple((x+j*y)%p for j in range(p))+(y,) for x,y in product(range(p),repeat=2)]
        assert all(sum(a!=0 for a in v)<=p for v in linear)
        assert all(any(v[i]!=0 for v in linear) for i in range(p+1))
    print('8 sharp small-field obstruction families checked.')

def simplex(d,D):
    if not d:yield ();return
    for x in range(D+1):
        for v in simplex(d-1,D-x):yield (x,)+v

def evalmat(X,E,p):
    X=np.asarray(X,dtype=np.int64);R=np.ones((len(X),len(E)),dtype=np.int64)
    for j in range(X.shape[1]):
        ex=np.array([e[j] for e in E]);T=np.ones((len(X),int(max(ex))+1),dtype=np.int64)
        for k in range(1,T.shape[1]):T[:,k]=T[:,k-1]*X[:,j]%p
        R=R*T[:,ex]%p
    return R

def ideal(X,p,D):
    d=len(X[0]);E=list(simplex(d,D));I=null(evalmat(X,E,p),p,len(E))
    return E,I

def tangent(x,E,I,p):
    d=len(x);Z=np.zeros((len(E),d),dtype=np.int64)
    for j in range(d):
        e2=[tuple(v-int(i==j and v>0) for i,v in enumerate(e)) for e in E]
        Z[:,j]=evalmat([x],e2,p)[0]*np.array([e[j] for e in E])%p
    return basis(null(I@Z%p,p,d),p,d)

def recover_sample(X,p,D):
    X=sorted(set(X));E,I=ideal(X,p,D);T={}
    for x in X:
        B=tangent(x,E,I,p);T[flat_key(x,B,p)]=(x,B)
    candidates=[]
    for x,B in sorted(T.values(),key=lambda z:len(z[1])):
        if any(len(C)<len(B) and inside(sub(y,x,p),B,p) and contains(C,B,p) for y,C in candidates):continue
        candidates.append((x,B))
    return candidates,E,I

def sample_grid(c,A,p,D):
    return {add(c,tuple(np.array(z,dtype=np.int64)@A%p),p) for z in simplex(len(A),D)}

def sample_complete(X,x,B,p,D):
    B=basis(B,p,len(x));J=rref(B,p)[1];C=inv(B[:,J],p);L=[]
    for y in X:
        v=sub(y,x,p)
        if inside(v,B,p):L.append(tuple(np.asarray(v)[J]@C%p))
    E=list(simplex(len(B),D));return bool(L) and len(basis(evalmat(L,E,p),p,len(E)))==len(E)

def check_samples():
    print('\nTHIN SAMPLE TANGENT RECOVERY')
    total=points=thresholds=0
    cases=[(17,[1]*6,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,5)],'none',3),
           (101,[1,2,1,1,1],[(0,1,2),(1,3),(4,)],'none',4),
           (101,[1]*8,[(i,j) for i in range(4) for j in range(4,8)],'none',3),
           (13,[1]*13,[(i,) for i in range(13)],'vander',3)]
    for p,dims,F,proj,D in cases:
        # Avoid enumerating complete facets: the simplex samples alone suffice.
        n=sum(dims);off=np.cumsum([0]+dims)
        A=np.eye(n,dtype=np.int64) if proj=='none' else np.array([[pow(i,j,p) for j in range(8)] for i in range(n)],dtype=np.int64)
        d=A.shape[1];A=A@rndmat(d,p)%p;c=tuple(rng.randrange(p) for _ in range(d));G=[A[off[i]:off[i+1]] for i in range(len(dims))]
        B=[np.vstack([G[i] for i in f]) for f in F];X=set.union(*(sample_grid(c,V,p,D) for V in B));R,E,I=recover_sample(X,p,D)
        truth={flat_key(c,V,p) for V in B};assert {flat_key(x,V,p) for x,V in R}==truth
        assert all(sample_complete(X,x,V,p,D) for x,V in R)
        if proj=='none':
            for delta in range(1,D+1):
                Ex,Ix=ideal(sorted(X),p,delta)
                for f in F:
                    x=add(c,tuple(np.vstack([G[i] for i in f]).sum(axis=0)%p),p)
                    # Full-support anchor: every coordinate in each active block is one.
                    visible=set(f)
                    for j in range(len(dims)):
                        if j in f:continue
                        sizes=[len(J)+1 for k in range(len(f)+1) for J in combinations(f,k) if not any(set(J)|{j}<=set(g) for g in F)]
                        if min(sizes)>delta:visible.add(j)
                    V=basis(np.vstack([G[i] for i in sorted(visible)]),p,d)
                    assert key(tangent(x,Ex,Ix,p))==key(V);thresholds+=1
        total+=len(R);points+=len(X)
        print(f'q={p} ambient={d} facets={len(F)} D={D} monomials={len(E)} samples={len(X)} sum facet sizes={sum(p**len(V) for V in B)}')
    # A missing component cannot be inferred from samples from other components.
    X=[(x,0) for x in range(5)];R,_,_=recover_sample(X,5,2)
    assert len(R)==1 and len(R[0][1])==1
    print(f'{total} exact facet recoveries from {points} samples; {thresholds} exact tangent-threshold comparisons; one absent-component negative control.')

def check_general_arrangements():
    print('\nBOUNDED-REDUNDANCY MAXIMAL-RANK ARRANGEMENTS')
    total=0;facets=0;points=0
    for p,d,dims,random_sample in [(101,5,[2,2,1,2,1,2,1],False),
                                  (101,5,[2]*7,True),
                                  (101,6,[1,2,3,2,2,1],False)]:
        # No common coordinate blocks, curve, jet or flag promise is imposed.
        for _ in range(100):
            G=[basis([[rng.randrange(p) for _ in range(d)] for j in range(r)],p,d) for r in dims]
            if [len(B) for B in G]!=dims:continue
            if all(len(basis(np.vstack(J),p,d))==min(d,sum(map(len,J)))
                   for k in range(1,len(G)+1) for J in combinations(G,k)):break
        else:raise AssertionError('No maximal-rank fixture found')
        total+=2**len(G)-1
        # Next-fit batches of other groups of total dimension < d.
        counts=[]
        for i in range(len(G)):
            batches=[];current=[];load=0
            for j,r in enumerate(dims):
                if j==i:continue
                if load+r>=d:batches.append(current);current=[];load=0
                current.append(j);load+=r
            if current:batches.append(current)
            for B in batches:
                Z=basis(np.vstack([G[j] for j in B]),p,d)
                assert len(Z)<d and not contains(G[i],Z,p)
            counts.append(len(batches))
        D=1+max(counts);assert p>D
        c=tuple(rng.randrange(p) for _ in range(d));X=set()
        for B in G:
            if random_sample:
                M=comb(len(B)+D,D)
                for _ in range(4*M):X.add(add(c,tuple(np.array([rng.randrange(p) for _ in B])@B%p),p))
            else:X|=sample_grid(c,B,p,D)
        R,E,I=recover_sample(X,p,D)
        assert {flat_key(x,B,p) for x,B in R}=={flat_key(c,B,p) for B in G}
        assert all(sample_complete(X,x,B,p,D) for x,B in R)
        facets+=len(G);points+=len(X)
        print(f'q={p} ambient={d} total rank={sum(dims)} facets={len(G)} D={D} samples={len(X)} random sample={random_sample}')
    print(f'{total} exact subset rank checks; {facets} generic maximal-rank facets recovered from {points} samples.')


def trim(a):
    a=list(map(int,a))
    while a and not a[-1]:a.pop()
    return a

def gcd(a,b,p):
    a,b=trim(a),trim(b)
    while b:
        while a and len(a)>=len(b):
            c=a[-1]*pow(b[-1],-1,p)%p;k=len(a)-len(b)
            for j,v in enumerate(b):a[j+k]=(a[j+k]-c*v)%p
            a=trim(a)
        a,b=b,a
    return [v*pow(a[-1],-1,p)%p for v in a] if a else []

def pmul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return trim(c)

def jet(f,r,d,p):
    rows=[]
    for j in range(d-r):
        a=[0]*d;a[j:j+len(f)]=f;rows.append(a)
    return basis(null(rows,p,d),p,d)

def certify_jets(V,p):
    d=V[0].shape[1];D=d-1;forms=[]
    for B in V:
        r=len(basis(B,p,d));N=null(B,p,d)
        if not (0<r<d):return None
        g=[];mx=0
        for a in N:
            b=trim(a);mx=max(mx,len(b)-1);g=gcd(g,b,p)
        infinity=D-mx
        if len(g)-1+infinity!=r:return None
        forms.append((g,infinity))
    for (a,i),(b,j) in combinations(forms,2):
        if (i and j) or len(gcd(a,b,p))>1:return None
    return forms

def frame_from_lines(V,p):
    d=V[0].shape[1];lines=[B[0] for B in V if len(B)==1]
    if len(lines)<d+2:return None
    B=np.array(lines[:d]).T
    if len(basis(B,p,d))!=d:return None
    BI=inv(B,p);a=BI@lines[d]%p;b=BI@lines[d+1]%p
    if not all(a) or not all(b):return None
    t=[int(y)*pow(int(x),-1,p)%p for x,y in zip(a,b)]
    if len(set(t))!=d:return None
    C=[]
    for i,u in enumerate(t):
        lam=1
        for j,v in enumerate(t):
            if i!=j:lam=lam*(-v)*pow((u-v)%p,-1,p)%p
        C.append([lam*pow(int(a[i]),-1,p)*pow(u,k,p)%p for k in range(d)])
    T=np.array(C,dtype=np.int64).T@BI%p
    return T if len(basis(T,p,d))==d else None

def check_jets():
    print('\nHOMOGENEOUS CRT JET CERTIFICATES')
    total=0;frames=0
    # Mixed group dimensions; repeated roots and infinity are intentional.
    cases=[(2,6,[([0,1],1),([1,1,1],0),([1,1,0,1],0)]),
           (7,6,[([1],2),([1,5,1],0),([1,0,1],0),([0,1],0)]),
           (23,6,[([(-i)%23,1],0) for i in range(8)]+[(pmul([(-9)%23,1],[(-9)%23,1],23),0),([1],2)]),
           (31,5,[([(-i)%31,1],0) for i in range(7)]+[(pmul([(-8)%31,1],[(-9)%31,1],31),0)])]
    for p,d,F in cases:
        V=[jet(f,len(f)-1+i,d,p) for f,i in F];assert certify_jets(V,p)
        if sum(len(B)==1 for B in V)>=d+2:
            M=rndmat(d,p);W=[basis(B@M%p,p,d) for B in V];T=frame_from_lines(W,p);assert T is not None
            V=[basis(B@T.T%p,p,d) for B in W];assert certify_jets(V,p);frames+=1
        n=0
        for k in range(1,len(V)+1):
            for J in combinations(V,k):
                assert len(basis(np.vstack(J),p,d))==min(d,sum(map(len,J)));n+=1
        total+=n;print(f'q={p} ambient={d} ranks={tuple(map(len,V))} verified subsets={n}')
    p=7;d=6
    assert certify_jets([jet([1,1],1,d,p),jet(pmul([1,1],[2,1],p),2,d,p)],p) is None
    V=jet([1,1,1],2,d,p);V[0,3]=(V[0,3]+1)%p;assert certify_jets([V],p) is None
    print(f'{total} exact subset-rank comparisons; {frames} hidden-frame recoveries; 2 negative certificates rejected.')

if __name__=='__main__':
    check_complete();check_samples();check_general_arrangements();check_jets();print('\nAll checks passed.')
