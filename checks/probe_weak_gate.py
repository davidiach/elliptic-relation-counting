"""Holdout for one fixed field-valued Gram readout on planted and matched seeds."""
from itertools import combinations, product
from math import log, sqrt

import numpy as np

from fixtures import curve_cycle, determinant, evaluation


def hankel_det(s,indices,p):
    if not indices:
        return 1
    def minor(rows,cols):
        if len(rows)==1:
            return s[:,rows[0]+cols[0]]
        return sum((-1)**j*s[:,rows[0]+v]*minor(rows[1:],cols[:j]+cols[j+1:])
                   for j,v in enumerate(cols)) % p
    return minor(indices,indices)


def gram_score(sums,p,a):
    count,m,_=sums.shape
    even=list(range((m+1)//2+1))
    odd=list(range(m//2))
    total=np.zeros(count,dtype=np.int64)
    for n in range(1,m+1):
        for subset in combinations(range(m),n):
            s=sums[:,subset,:].sum(axis=1)%p
            o=(s[:,3:]+a*s[:,1:-2]+s[:,:-3])%p
            value=(hankel_det(s,even[1:],p)*hankel_det(o,odd,p)
                   +hankel_det(s,even,p)*hankel_det(o,odd[1:],p))
            total+=(-1)**(m-n)*value
    return total%p


def check_score():
    p,a=13,4
    curve=curve_cycle(p,a)
    for m in (4,6):
        for shift in (0,2):
            slots=[]
            for j in range(m):
                ids=[1+(j+shift)%8]
                if j==1:
                    ids.append(1+(j+shift+3)%8)
                slots.append([curve[i] for u in ids for i in (u,len(curve)-u)])
            sums=np.array([[sum(pow(x,d,p) for x,y in slot)%p for d in range(m+2)]
                           for slot in slots])
            z=int(gram_score(sums[None,:,:],p,a)[0])
            vectors=[[evaluation(P,m+1,p) for P in slot] for slot in slots]
            t=evaluation(curve[1],m+1,p)
            assert z==sum(determinant([t,*rows],p)**2 for rows in product(*vectors))%p
    print("Gram readout agrees with direct signed determinant sums for m=4 and m=6")


def run(p,a,widths,draws,seed):
    curve=curve_cycle(p,a)
    N=len(curve)
    m=len(widths)
    M=1<<p.bit_length()
    xs=np.array([-1]+[pt[0] for pt in curve[1:]])
    assert curve[1]==(0,1)
    powers=np.zeros((M,m+2),dtype=np.int64)
    for x,y in curve[1:]:
        powers[x]+=np.array([pow(x,j,p) for j in range(m+2)])
    powers%=p
    inv=np.array([pow(u,-1,M) for u in range(1,M,2)])
    squares=np.zeros(p,dtype=bool)
    squares[(np.arange(1,p)**2)%p]=True
    rates=[]
    for planted in (False,True):
        rng=np.random.default_rng(seed+int(planted))
        yes=0
        for start in range(0,draws,512):
            count=min(draws-start,512)
            labels=rng.integers(1,N,size=(count,m))
            if planted:
                labels[:,-1]=(1-labels[:,:-1].sum(axis=1))%N
                bad=labels[:,-1]==0
                while bad.any():
                    labels[bad,:-1]=rng.integers(1,N,size=(int(bad.sum()),m-1))
                    labels[bad,-1]=(1-labels[bad,:-1].sum(axis=1))%N
                    bad=labels[:,-1]==0
            alphas=rng.integers(0,M//2,size=(count,m))
            slots=[]
            for i,w in enumerate(widths):
                addr=rng.integers(0,w,size=count)
                beta=(addr-(2*alphas[:,i]+1)*xs[labels[:,i]])%M
                ids=((np.arange(w)-beta[:,None])*inv[alphas[:,i],None])%M
                slots.append(powers[ids].sum(axis=1)%p)
            scores=gram_score(np.stack(slots,axis=1),p,a)
            yes+=int(squares[scores].sum())
        rates.append(yes/draws)
    band=sqrt(log(2/.05)/draws)
    print(f"p={p}, m={m}, n={draws} per law, matched {rates[0]:.5f}, planted {rates[1]:.5f}, "
          f"gap {rates[1]-rates[0]:+.5f}, 95% Hoeffding half-width {band:.5f}")


if __name__=="__main__":
    check_score()
    for mwidths,seed in (((8,8,16,16),720641),((4,4,4,4,8,8),720661)):
        run(65521,44,mwidths,1000000,seed)
