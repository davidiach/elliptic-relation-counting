#!/usr/bin/env python3
"""Exact tests for affine rank-ball center avoidance in abstract matrix spaces."""
from itertools import product
import numpy as np
from verify import rank,basis,null,inv,random_invertible,random_rank

if not __debug__:
    raise RuntimeError('Assertions required')


def canonical(C,p):
    D=C.copy()%p;a,b=D.shape;P=np.eye(a,dtype=np.int64);Q=np.eye(b,dtype=np.int64)
    t=0
    while t<min(a,b):
        ij=np.argwhere(D[t:,t:])
        if not len(ij):break
        i,j=map(int,ij[0]+t)
        D[[t,i]]=D[[i,t]];P[[t,i]]=P[[i,t]]
        D[:,[t,j]]=D[:,[j,t]];Q[:,[t,j]]=Q[:,[j,t]]
        s=pow(int(D[t,t]),-1,p);D[t]=s*D[t]%p;P[t]=s*P[t]%p
        for i in range(a):
            if i!=t:
                s=D[i,t];D[i]=(D[i]-s*D[t])%p;P[i]=(P[i]-s*P[t])%p
        for j in range(b):
            if j!=t:
                s=D[t,j];D[:,j]=(D[:,j]-s*D[:,t])%p;Q[:,j]=(Q[:,j]-s*Q[:,t])%p
        t+=1
    J=np.zeros_like(D);J[:t,:t]=np.eye(t,dtype=np.int64)
    assert np.array_equal(D,J) and np.array_equal(P@C@Q%p,D)
    return t,P,Q


def rootless_block(k,p):
    assert k in (2,3)
    for a in product(range(p),repeat=k):
        if all((pow(x,k,p)+sum(v*pow(x,i,p) for i,v in enumerate(a)))%p for x in range(p)):
            M=np.zeros((k,k),dtype=np.int64)
            for i in range(k-1):M[i+1,i]=1
            M[:,-1]=-np.array(a)%p
            assert all(rank((M-x*np.eye(k,dtype=np.int64))%p,p)==k for x in range(p))
            return M
    raise AssertionError('irreducible block missing')


def line_avoid(C,p):
    a,b=C.shape
    if a>b:return line_avoid(C.T,p).T
    assert a>=2
    t,P,Q=canonical(C,p)
    if t<a:
        B=np.roll(np.eye(a,dtype=np.int64),1,axis=0)
    else:
        parts=[2]*(a//2) if a%2==0 else [3]+[2]*((a-3)//2)
        B=np.zeros((a,a),dtype=np.int64);i=0
        for k in parts:B[i:i+k,i:i+k]=rootless_block(k,p);i+=k
    Z=np.zeros_like(C);Z[:,:a]=B
    X=inv(P,p)@Z@inv(Q,p)%p
    return X


def extend(U,n,p):
    S=[x.copy() for x in U]
    for v in np.eye(n,dtype=np.int64):
        if rank(np.array(S+[v]),p)>len(S):S.append(v)
    return np.array(S,dtype=np.int64)


def common_space_avoid(W,p):
    # All W entries have the same shape, and W is nonempty.
    a,b=W[0].shape
    if a>b:
        Z=common_space_avoid([w.T for w in W],p)
        return None if Z is None else Z.T
    U=basis(np.hstack(W).T,p);K=null(np.vstack(W),p);t=len(U)
    if len(K)<t:return None
    D=extend(K[:t],b,p).T;R=extend(U,a,p).T
    F=np.zeros((a,b),dtype=np.int64);F[:,:a]=R
    return F@inv(D,p)%p


def main():
    print('AFFINE CENTER AVOIDANCE: FULL-RANK MATRIX PENCILS')
    rng=np.random.default_rng(27112026);fixtures=checks=0
    for p in (2,3,5):
        for c in product(range(p),repeat=4):
            C=np.array(c,dtype=np.int64).reshape(2,2);Z=line_avoid(C,p)
            for lam in range(p):assert rank((Z-lam*C)%p,p)==2;checks+=1
            fixtures+=1
    for p in (2,3,7,13):
        for a,b in ((2,3),(3,2),(3,3),(4,5),(5,4),(7,7)):
            for t in range(min(a,b)+1):
                for _ in range(4):
                    C=random_rank(rng,a,b,t,p).reshape(a,b);Z=line_avoid(C,p)
                    for lam in range(p):assert rank((Z-lam*C)%p,p)==min(a,b);checks+=1
                    fixtures+=1
    print(f'{fixtures} center matrices, {checks} full-rank scalar-shift checks.')
    print('COMMON KERNEL/IMAGE CENTER SPACES')
    cases=checks=0
    for p,a,b,t,c in ((3,4,4,2,2),(5,4,5,3,2),(3,5,5,2,3),(7,3,5,2,3),(3,4,3,1,2),(5,4,4,0,2)):
        for _ in range(8):
            # Work with a<=b, then transpose the rectangular reverse fixture.
            reverse=a>b;u,v=(b,a) if reverse else (a,b)
            assert v-c>=t
            P=random_invertible(rng,u,p);Q=random_invertible(rng,v,p);W=[]
            for j in range(3):
                A=np.zeros((u,v),dtype=np.int64);A[:t,:c]=rng.integers(0,p,(t,c))
                B=P@A@Q%p;W.append(B.T if reverse else B)
            Z=common_space_avoid(W,p);assert Z is not None
            for alpha in product(range(p),repeat=len(W)):
                C=sum(int(x)*w for x,w in zip(alpha,W))%p
                assert rank((Z-C)%p,p)==min(a,b);checks+=1
            cases+=1
    assert common_space_avoid([np.eye(3,dtype=np.int64)],5) is None
    print(f'{cases} center spaces, {checks} full-rank checks; insufficient common kernel control rejected.')
    print('All checks passed.')

if __name__=='__main__':main()
