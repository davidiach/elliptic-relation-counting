#!/usr/bin/env python3
"""Exact symbolic mixed-addition norm checks for inverse shared translations.
Root search in the finite regression uses exhaustive field evaluation; the
asymptotic algorithm uses standard univariate finite-field root factoring.
"""
from itertools import product
from random import Random
from functools import lru_cache
import numpy as np
from verify import good_offsets

class Ring:
    def __init__(self,p,a,c,offsets):
        self.p=p;self.d=len(offsets)
        self.f=[tuple(v%p for v in (s**3+a*s+c,3*s*s+a,3*s,1)) for s in offsets]
    def trim(self,v):
        v=[int(x)%self.p for x in v]
        while v and not v[-1]:v.pop()
        return tuple(v)
    def pa(self,a,b):
        return self.trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
    def pm(self,a,b):
        if not a or not b:return ()
        return self.trim(np.convolve(np.array(a,dtype=np.int64),np.array(b,dtype=np.int64)))
    @lru_cache(None)
    def fac(self,mask):
        if not mask:return (1,)
        q=mask&-mask
        return self.pm(self.fac(mask^q),self.f[q.bit_length()-1])
    def C(self,x):return {0:(x%self.p,)} if x%self.p else {}
    def A(self,s):return {0:(s%self.p,1)}
    def Y(self,j):return {1<<j:(1,)}
    def add(self,a,b):
        out=dict(a)
        for j,v in b.items():
            w=self.pa(out.get(j,()),v)
            if w:out[j]=w
            elif j in out:del out[j]
        return out
    def scale(self,a,s):return {j:self.trim([x*s for x in v]) for j,v in a.items() if self.trim([x*s for x in v])}
    def sub(self,a,b):return self.add(a,self.scale(b,-1))
    def mul(self,a,b):
        out={}
        for j,v in a.items():
            for k,w in b.items():
                z=self.pm(self.pm(v,w),self.fac(j&k));idx=j^k
                zz=self.pa(out.get(idx,()),z)
                if zz:out[idx]=zz
                elif idx in out:del out[idx]
        return out
    def sq(self,a):return self.mul(a,a)
    def norm(self,a):
        for j in range(self.d):
            cj={k:(self.trim([-x for x in v]) if k>>j&1 else v) for k,v in a.items()}
            a=self.mul(a,cj)
            assert all(not (k>>j&1) for k in a)
        assert set(a)=={0},'generic rational numerator must be nonzero'
        return a[0]
    def evaluate(self,a,x,ys):
        v=0
        for mask,poly in a.items():
            t=0
            for c in reversed(poly):t=(t*x+c)%self.p
            for j,y in enumerate(ys):
                if mask>>j&1:t=t*y%self.p
            v+=t
        return v%self.p
    def evalpoly(self,a,x):
        v=0
        for c in reversed(a):v=(v*x+c)%self.p
        return v

def symbolic_sum(R,offsets):
    X,Y,Z=R.A(offsets[0]),R.Y(0),R.C(1)
    for j,s in enumerate(offsets[1:],1):
        z2=R.sq(Z);z3=R.mul(Z,z2)
        H=R.sub(R.mul(R.A(s),z2),X)
        D=R.sub(R.mul(R.Y(j),z3),Y)
        H2=R.sq(H);H3=R.mul(H,H2);V=R.mul(X,H2)
        X1=R.sub(R.sub(R.sq(D),H3),R.scale(V,2))
        Y1=R.sub(R.mul(D,R.sub(V,X1)),R.mul(Y,H3))
        Z1=R.mul(Z,H)
        X,Y,Z=X1,Y1,Z1
    return X,Y,Z

def run(p,a,c,offsets):
    assert good_offsets(p,a,c,offsets)
    R=Ring(p,a,c,offsets);X,Y,Z=symbolic_sum(R,offsets)
    Z2=R.sq(Z)
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        s=((v-y)*pow(u-x,-1,p) if x!=u else (3*x*x+a)*pow(2*y,-1,p))%p
        z=(s*s-x-u)%p
        return z,(s*(x-z)-y)%p
    lifts=[[y for y in range(1,p) if (y*y-x*x*x-a*x-c)%p==0] for x in range(p)]
    pts=[(x,y) for x in range(p) for y in lifts[x]]
    targets=sorted(set(pts)); rng=Random(781)
    targets=[None]+rng.sample(targets,min(5,len(targets)))
    # direct group sums and independent symbolic projective evaluation
    truth={T:set() for T in targets}; tuples=0;exceptional=0
    allrows=[]
    for shift in range(p):
        xs=[(shift+s)%p for s in offsets]
        for ys in product(*(lifts[x] for x in xs)):
            Q=None
            for x,y in zip(xs,ys):Q=add(Q,(x,y))
            tuples+=1
            zv=R.evaluate(Z,shift,ys)
            if zv:
                xv=R.evaluate(X,shift,ys)*pow(zv,-2,p)%p
                yv=R.evaluate(Y,shift,ys)*pow(zv,-3,p)%p
                assert Q==(xv,yv)
            else:exceptional+=1
            if Q in truth:truth[Q].add(shift)
            allrows.append((shift,ys,Q))
    checks=0;maxdeg=0;roots=0;hits=0
    for T in targets:
        # Include every bad addition denominator by retaining the final Z.
        N=Z if T is None else R.mul(Z,R.sub(X,R.scale(Z2,T[0])))
        P=R.norm(N);maxdeg=max(maxdeg,len(P)-1)
        cand={s for s in range(p) if R.evalpoly(P,s)==0}
        # Verify the norm identity independently at every rational lift.
        for shift,ys,Q in allrows:
            # Its value is the product over the independent sign conjugates.
            v=1
            for sg in product((-1,1),repeat=len(offsets)):
                v=v*R.evaluate(N,shift,[a*b%p for a,b in zip(ys,sg)])%p
            assert v==R.evalpoly(P,shift)
            checks+=1
        assert truth[T]<=cand
        got={s for s,ys,Q in allrows if s in cand and Q==T}
        assert got==truth[T]
        roots+=len(cand);hits+=len(got)
    print(f'inverse p={p}, m={len(offsets)}: {tuples} lift tuples, {exceptional} exceptional paths, {len(targets)} targets, {checks} norm evaluations, {roots} candidate translations, {hits} true translations, max polynomial degree {maxdeg}')
    return tuples,exceptional,len(targets),checks,roots,hits

if __name__=='__main__':
    totals=[0]*6
    for args in [(13,4,1,[0,1]),(17,2,2,[0,1,2]),(37,0,5,[0,1,3,4]),(43,2,1,[0,1,2,3])]:
        z=run(*args);totals=[a+b for a,b in zip(totals,z)]
    print('inverse totals (tuples, exceptional paths, targets, norm checks, candidate translations, true translations):',totals)
    print('PASS: symbolic inverse translations, exact toy root checks.')
