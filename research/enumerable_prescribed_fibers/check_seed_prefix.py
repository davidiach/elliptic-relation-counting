"""Exact finite conditional expectation on affine masks; small curve only."""
from itertools import product, combinations
from fractions import Fraction as F
from math import prod
from functools import lru_cache
import numpy as np
from check_fibers import conv, fiber_rank, exact_pin_data, covariance_check

# E/F5, y^2=x^3+2x+1, of prime order 7.
p=5;r=7;m=4;T=4;anchors=(0,1,3)

def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;z,w=Q
    if x==z and (y+w)%p==0:return None
    s=((3*x*x+2)*pow(2*y,-1,p) if P==Q else (w-y)*pow(z-x,-1,p))%p
    X=(s*s-x-z)%p;return X,(s*(x-X)-y)%p
pts=[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-2*x-1)%p==0]
assert len(pts)==r-1
G=[None];P=pts[0]
for _ in range(r-1):G.append(add(G[-1],P))
assert set(G)=={None,*pts} and add(G[-1],P) is None
q=F((r-1)*T,p);mu=q**m/r
# Each parity's completion list. Condition independently on one affine word.
allseeds=tuple(product(range(p),repeat=2));prefix=[[allseeds,allseeds] for _ in range(m)]
Isets=[I for s in (0,1,2) for I in combinations(range(m),s)]

@lru_cache(None)
def one_mask(par,S):
    selected=[]
    for a,c in S:selected.append(tuple(g for g in range(1,r) if G[g][1]%2==par and (a*G[g][0]+c)%p<T))
    F1=[sum(x in V for V in selected) for x in range(r)]
    F2=[[sum(x in V and y in V for V in selected) for y in range(r)] for x in range(r)]
    return F1,F2

def data(a,b):
    A,C=one_mask(0,a);B,D=one_mask(1,b);den=len(a)*len(b)
    H=np.zeros((r,r),dtype=object);one=[0]*r
    for x in range(1,r):
        px=G[x][1]%2;one[x]=A[x]*len(b) if px==0 else B[x]*len(a)
        for y in range(1,r):
            py=G[y][1]%2
            H[x,y]=C[x][y]*len(b) if px==py==0 else D[x][y]*len(a) if px==py==1 else A[x]*B[y] if px==0 else B[x]*A[y]
    return H,den,sum(one),int(H.sum())

def conv2(A,B):
    Z=np.zeros((r,r),dtype=object)
    for i in range(r):
        for j in range(r):
            if B[i,j]:Z+=int(B[i,j])*np.roll(A,(i,j),axis=(0,1))
    return Z

def expectation(pref):
    Ds=[data(a,b) for a,b in pref];den=prod(a[1] for a in Ds)
    cost=m**4*sum((F(n2,d)-2*q*F(n1,d)+q*q)/(q*q) for H,d,n1,n2 in Ds)
    for I in Isets:
        Z=np.zeros((r,r),dtype=object);Z[0,0]=1
        for j,(H,d,n1,n2) in enumerate(Ds):
            A=np.diag(np.diag(H)) if j in I else H
            Z=conv2(Z,A)
        for t in anchors:
            e2=F(r*r*int(Z[t,t])-r*(int(Z[t,:].sum())+int(Z[:,t].sum()))+int(Z.sum()),r*r*den)
            assert e2>=0
            cost+=e2/(mu*mu)*(q*q if I else 1)
    return cost

def direct(Bs):
    cost=m**4*sum((F(len(B))-q)**2/(q*q) for B in Bs)
    for t in anchors:
        for I in Isets:
            J=[j for j in range(m) if j not in I];W=prod(len(Bs[j]) for j in J);D=conv(r,[Bs[j] for j in J]);e2=F(0)
            for pp in product(*(Bs[i] for i in I)):
                e2+=(F(int(D[(t-sum(pp))%r]))-F(W,r))**2
            cost+=e2/(mu*mu)*(q*q if I else 1)
    return cost

# Split each seed domain into two nonempty halves; conditional expectation
# weighted-average identity must hold exactly at every prefix.
old=expectation(prefix);start=old;steps=0
for i in range(m):
    for par in (0,1):
        while len(prefix[i][par])>1:
            S=prefix[i][par];k=len(S)//2;left,right=S[:k],S[k:]
            prefix[i][par]=left;v=expectation(prefix)
            prefix[i][par]=right;w=expectation(prefix)
            assert old*len(S)==v*len(left)+w*len(right)
            prefix[i][par]=left if v<=w else right;new=min(v,w);assert new<=old;old=new;steps+=1
Bs=[]
for a,b in prefix:
    A,C=a[0];B,D=b[0]
    Bs.append([g for g in range(1,r) if ((A*G[g][0]+C) if G[g][1]%2==0 else (B*G[g][0]+D))%p<T])
assert direct(Bs)==old
print(f'Actual E/F5, order 7: {steps} exact conditional-expectation prefix identities; initial potential {float(start):.6f}; final {float(old):.6f}; direct equality exact.')
print(f'  Selected block sizes {[len(B) for B in Bs]}; anchors {anchors}; finite guard {old<F(1,64)**2} (asymptotic K not available here).')
for t in anchors:
    D=sum(len(B)-1 for B in Bs);d=int(conv(r,Bs)[t]);rank,_=fiber_rank(r,Bs,t)
    # Verify the actual group law for every tuple in each tested fiber.
    ct=0
    for X in product(*Bs):
        if sum(X)%r!=t:continue
        A=None
        for g in X:A=add(A,G[g])
        assert A==G[t];ct+=1
    assert ct==d
    print(f'  t={t}: {d} genuine short tuples; modular contrast rank {rank}/{D-1}.')
