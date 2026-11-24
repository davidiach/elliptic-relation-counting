#!/usr/bin/env python3
"""Finite tests of affine-sampler Fourier identities. Integer histogram checks,
exact rational L2 checks, floating-point SVD checks clearly labeled. No large
ECDLP runs, no claim of end-to-end speedup.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import prod, isqrt
import random
import numpy as np

rng = random.Random(20261124)

def prime(n): return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def curve(p):
    for a,c in product(range(p), repeat=2):
        if (4*a**3+27*c*c)%p==0: continue
        pts=[None]+[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-c)%p==0]
        r=len(pts)
        if r<=3 or not prime(r) or r==2: continue
        def add(u,v):
            if u is None:return v
            if v is None:return u
            x,y=u; z,t=v
            if x==z:
                if (y+t)%p==0:return None
                s=(3*x*x+a)*pow(2*y,-1,p)%p
            else:s=(t-y)*pow(z-x,-1,p)%p
            xx=(s*s-x-z)%p
            return xx,(s*(x-xx)-y)%p
        g=pts[1]; group=[None]
        for j in range(1,r): group.append(add(group[-1],g))
        assert len(set(group))==r and add(group[-1],g) is None
        unused=next((x for x in range(p) if all(P is None or P[0]!=x for P in group)),None)
        if unused is None:continue
        xs=np.array([unused if P is None else P[0] for P in group],dtype=np.int64)
        assert all(xs[j]==xs[(-j)%r] for j in range(r))
        return a,c,r,xs
    raise ValueError('no fixture')

def seeds(p, zero=False):return list(product(range(p),range(0 if zero else 1,p)))

def incidence(p,S,zero=False):
    # Multiplicity at t=0 is essential.
    A=np.zeros((len(seeds(p,zero)),p),dtype=np.int64)
    for i,(a,t) in enumerate(seeds(p,zero)):
        for z in S:A[i,(a+t*z)%p]+=1
    return A

def support(p,w):return tuple(sorted(rng.sample(range(p),w)))

def tensor_num(p,xs,ss,T,zero=False,aff=False):
    r=len(xs);m=len(ss)
    mats=[incidence(p,S,zero)[:,xs] for S in ss]
    N=mats[0].shape[0]
    points=range(1 if aff else 0,r)
    H=np.zeros((N,)*m,dtype=np.int64)
    R=np.ones_like(H)
    def outer(vecs):
        ans=vecs[0]
        for v in vecs[1:]:ans=np.multiply.outer(ans,v)
        return ans
    for inds in product(points,repeat=m-1):
        last=(T-sum(inds))%r
        if aff and last==0:continue
        H+=outer([M[:,j] for M,j in zip(mats,inds+(last,))])
    R=outer([M[:,1 if aff else 0:].sum(axis=1) for M in mats])
    if not aff: return r*H-R,H,R
    nsum=sum(1 for inds in product(points,repeat=m-1) if (T-sum(inds))%r!=0)
    return (r-1)**m*H-nsum*R,H,R

def rational_norm(Nmat,p,r,ss):
    m=len(ss); W=prod(map(len,ss)); scale=F(p**m,r**m*W)
    exact=F(int(np.sum(Nmat*Nmat)),Nmat.size)*scale*scale
    d=prod(F(p*(p-len(S)),r*len(S)*(p-1)) for S in ss)
    expected=d*((r-2)*F(r-1,r)**m+(3-r)*F(-1,r)**m)
    assert exact==expected, (exact,expected)
    return exact,d

def hist(values,weights,p):
    out=np.zeros(p,dtype=np.int64)
    np.add.at(out,values.ravel()%p,weights.ravel())
    return out

def convolve(a,b,p):
    z=np.zeros(p,dtype=np.int64)
    for i in range(p):z+=a[i]*np.roll(b,i)
    return z

def cycl_equal(a,b):
    # Since p is prime, equality at a primitive p-th root iff the integer
    # coefficient difference is constant, i.e. a multiple of Phi_p.
    c=a-b;return bool(np.all(c==c[0]))

def group_hist(xs,xi,T,p):
    r=len(xs); m=len(xi)
    a=np.indices((r,)*m,dtype=np.int64).reshape(m,-1)
    phase=np.sum(np.array(xi)[:,None]*xs[a],axis=0)%p
    weights=r*(np.sum(a,axis=0)%r==T).astype(np.int64)-1
    return hist(phase,weights,p)

# I. Exact incidence Gram matrices for all widths and many shapes.
num_gram=0;entries=0;curves=[]
for p in (5,7,11,13,17,19):
    a,c,r,xs=curve(p);curves.append((p,a,c,r,xs))
    for w in range(1,p):
        choices=list(combinations(range(p),w)) if p<=7 else [support(p,w) for _ in range(8)]
        for S in choices:
            A=incidence(p,S); B=A.T@A
            want=np.full((p,p),w*(w-1),dtype=np.int64)
            np.fill_diagonal(want,(p-1)*w)
            assert np.array_equal(B,want)
            assert np.all(A.sum(axis=1)==w)
            num_gram+=1;entries+=p*p
print('Curve fixtures (p, a, c, r):',[(p,a,c,r) for p,a,c,r,_ in curves],flush=True)
print('Exact incidence Gram matrices:',num_gram,'matrix entries:',entries,flush=True)

# II. Full extended additive Fourier spectrum.
full_checks=0; active=0; killed=0
for p,m in ((5,2),(7,2),(5,3),(7,3)):
    a,c,r,xs=curve(p); ss=[support(p,2+(i%2)) for i in range(m)]
    T=1; E,_,_=tensor_num(p,xs,ss,T,zero=True)
    aa,tt=np.array(seeds(p,True),dtype=np.int64).T
    idx=np.indices(E.shape,dtype=np.int64).reshape(m,-1)
    Eflat=E.ravel();A=aa[idx];Tt=tt[idx]
    if m==2:
        freqs=list(product(range(p),repeat=2*m))
    else:
        freqs=[]
        for _ in range(150):
            xi=[rng.randrange(1,p) for _ in range(m)]
            eta=[x*rng.choice(S)%p for x,S in zip(xi,ss)]
            freqs.append(tuple(xi+eta))
            bad=eta.copy();j=rng.randrange(m)
            bad[j]=next(y for y in range(p) if y*pow(xi[j],-1,p)%p not in ss[j])
            freqs.append(tuple(xi+bad))
            xi[j]=0;eta[j]=0;freqs.append(tuple(xi+eta))
    cache={}
    for v in freqs:
        xi,eta=v[:m],v[m:]
        phase=(np.array(xi)@A+np.array(eta)@Tt)%p
        lhs=hist(phase,Eflat,p)
        if all(x and y*pow(x,-1,p)%p in S for x,y,S in zip(xi,eta,ss)):
            if xi not in cache:cache[xi]=group_hist(xs,xi,T,p)
            rhs=p**m*cache[xi];active+=1
        else:rhs=np.zeros(p,dtype=np.int64);killed+=1
        assert cycl_equal(lhs,rhs), (p,m,v,lhs,rhs)
        full_checks+=1
    print('Extended spectrum fixture:',p,m,'frequencies:',len(freqs),flush=True)
print('Exact cyclotomic coefficient checks:',full_checks,'on lines:',active,'off lines:',killed,flush=True)

# III. Translation Fourier for arbitrary revealed nonzero multipliers.
cond=0
for p,m in ((5,2),(7,2),(7,3),(11,3),(5,4)):
    a,c,r,xs=curve(p);ss=[support(p,2+(i%2)) for i in range(m)]
    E,_,_=tensor_num(p,xs,ss,1)
    for _ in range(100):
        ts=[rng.randrange(1,p) for _ in range(m)]
        xi=[rng.randrange(p) for _ in range(m)]
        ids=[slice(t-1,None,p-1) for t in ts]
        v=E[tuple(ids)]
        assert v.shape==(p,)*m
        idx=np.indices(v.shape,dtype=np.int64).reshape(m,-1)
        left=hist(np.array(xi)@idx,v.ravel(),p)
        right=group_hist(xs,xi,1,p)
        for x,t,S in zip(xi,ts,ss):
            u=np.zeros(p,dtype=np.int64)
            for z in S:u[-x*t*z%p]+=1
            right=convolve(right,u,p)
        assert np.array_equal(left,right), (p,m,ts,xi,left,right)
        cond+=1
print('Exact conditional-translation coefficient histograms:',cond,flush=True)

# IV. Exact likelihood L2, and numerical singular spectra under support changes.
# Shapes small enough to compute full integer tensors, SVD only on selected cuts.
normchecks=0; svdchecks=0; singular_pairs=0
for p,m in ((5,2),(7,2),(11,2),(5,3),(7,3),(5,4)):
    a,c,r,xs=curve(p);widths=[2+i%2 for i in range(m)]
    for T in range(1,min(r,4)):
        svs=[]
        for rep in range(2):
            ss=[support(p,w) for w in widths]
            E,_,_=tensor_num(p,xs,ss,T)
            val,d=rational_norm(E,p,r,ss);normchecks+=1
            # SVD uses a one-slot vs rest split; exact identities above do
            # not depend on any numerical tolerance.
            N=p*(p-1);scale=float(F(p**m,r**m*prod(widths)))
            mat=E.reshape(N,-1).astype(float)*(scale/(N**(m/2)))
            sv=np.linalg.svd(mat,compute_uv=False)
            n=(r-1)//2
            assert sv[0]<=2*np.sqrt(float(d))+1e-10
            assert np.all(sv[:n]>1e-10) and np.max(sv[n:],initial=0)<1e-10
            assert abs(np.sum(sv**2)-float(val))<1e-9
            svs.append(sv[:n]);svdchecks+=1
        assert np.max(np.abs(svs[0]-svs[1]))<1e-10
        singular_pairs+=1
print('Exact rational L2 formulas:',normchecks,flush=True)
print('Numerical cut SVD checks:',svdchecks,'support-shape spectrum pairs:',singular_pairs,flush=True)

# V. Full-group sentinel extension versus actual affine-only marking.
conditioning=0
for p,m in ((5,2),(7,2),(5,3),(7,3),(5,4)):
    a,c,r,xs=curve(p);ss=[support(p,2) for _ in range(m)];W=prod(map(len,ss));L=(p*(p-1))**m
    for T in (1,2):
        _,H,R=tensor_num(p,xs,ss,T)
        _,Ha,Ra=tensor_num(p,xs,ss,T,aff=True)
        nsum=sum(1 for inds in product(range(1,r),repeat=m-1) if (T-sum(inds))%r!=0)
        if not nsum:continue
        # Compare integer numerators of seed probability measures, not densities.
        nq=r**(m-1); nr=r**m; na=(r-1)**m
        tvq=F(int(np.abs(H*nsum-Ha*nq).sum()),2*nq*nsum*((p-1)**m)*W)
        tvr=F(int(np.abs(R*na-Ra*nr).sum()),2*nr*na*((p-1)**m)*W)
        assert tvq<=F(m,r) and tvr<=F(m,r)
        conditioning+=2
print('Exact conditioning-TV checks:',conditioning,flush=True)
print('All assertions passed.',flush=True)

# VI. Exact average retained energy for random rank-(m-1) projections.
# For each line in F_p^m, take the coset sums along that line (the kernel
# of the projection). Compare the total retained L2 energy over all lines.
projection_checks=0;coset_partitions=0
for p,m in ((5,2),(7,2),(5,3),(7,3)):
    a,c,r,xs=curve(p);ss=[support(p,2) for _ in range(m)]
    E,_,_=tensor_num(p,xs,ss,1)
    lines=[]
    for dvec in product(range(p),repeat=m):
        j=next((i for i,x in enumerate(dvec) if x),None)
        if j is not None and dvec[j]==1:lines.append(dvec)
    assert len(lines)==(p**m-1)//(p-1)
    for _ in range(6):
        ts=[rng.randrange(1,p) for _ in range(m)]
        v=E[tuple(slice(t-1,None,p-1) for t in ts)]
        assert int(v.sum())==0
        full=int(np.sum(v*v))
        vals=0
        for dvec in lines:
            j=next(i for i,x in enumerate(dvec) if x)
            total=0
            for rest in product(range(p),repeat=m-1):
                base=list(rest[:j]) + [0] + list(rest[j:])
                u=sum(int(v[tuple((x+k*d)%p for x,d in zip(base,dvec))]) for k in range(p))
                total+=u*u
            vals+=total
            coset_partitions+=1
        # Retained energy = sum of squared coset sums / p; full energy
        # = sum of squared entries. Their ratio averages to (p^(m-1)-1)/(p^m-1).
        assert F(vals,p*len(lines))==F(p**(m-1)-1,p**m-1)*full
        projection_checks+=1
print('Exact random-projection energy identities:',projection_checks,'coset partitions:',coset_partitions,flush=True)
print('Final verifier passed.',flush=True)

# VII. Exact high-order retention for intervals (continued-fraction
# rational reconstruction; complexity independent of interval length).
from math import gcd, lcm

def bounded_ratio(t,p,w):
    """Unique reduced A/B with |A|, B <= w-1, B>0, if any."""
    if t==0:return (0,1)
    n,d=t,p;h0,h1=0,1;k0,k1=1,0
    while d:
        q,n,d=n//d,d,n%d
        h0,h1=h1,q*h1+h0;k0,k1=k1,q*k1+k0
        if 1<=k1<=w-1:
            A=t*k1-h1*p
            if abs(A)<=w-1 and gcd(A,k1)==1:return A,k1
    return None

def interval_count(xs,p,w):
    xx=list(dict.fromkeys(xs));assert len(xx)>=2
    inv=pow((xx[1]-xx[0])%p,-1,p)
    ratios=[(0,1),(1,1)];B=1
    for x in xx[2:]:
        v=bounded_ratio((x-xx[0])*inv%p,p,w)
        if v is None:return 0
        ratios.append(v);B=lcm(B,v[1])
        if B>=w:return 0
    a=[a*(B//b) for a,b in ratios]
    span=max(a)-min(a)
    if span>=w:return 0
    k=(w-1)//span
    return 2*k*w-span*k*(k+1)

def literal_count(xs,p,w):
    xx=list(dict.fromkeys(xs));assert len(xx)>=2
    inv=pow((xx[1]-xx[0])%p,-1,p)
    ts=[(x-xx[0])*inv%p for x in xx[2:]]
    return sum(all((s+z*(t-s))%p<w for z in ts)
               for s,t in product(range(w),repeat=2) if s!=t)

ratio_checks=0;interval_checks=0;nonzero=0
for p,w in ((13,3),(29,4),(101,7),(257,11),(1031,23),(4093,45),(65521,181)):
    assert prime(p) and 2*(w-1)**2<p
    tab={}
    for a in range(-(w-1),w):
        for b in range(1,w):
            if gcd(a,b)==1:
                z=a*pow(b,-1,p)%p
                assert z not in tab
                tab[z]=(a,b)
    for z in range(p):
        assert bounded_ratio(z,p,w)==tab.get(z),(p,w,z,bounded_ratio(z,p,w),tab.get(z))
        ratio_checks+=1
    trials=50 if w<=45 else 20
    for trial in range(trials):
        k=rng.randrange(2,11)
        if trial%2==0:
            a=rng.randrange(p);t=rng.randrange(1,p)
            ss=[rng.randrange(w) for _ in range(k)]
            if ss[1]==ss[0]:ss[1]=(ss[0]+1)%w
            xs=[(a+t*z)%p for z in ss]
        else:
            xs=rng.sample(range(p),k)
        expected=literal_count(xs,p,w)
        assert interval_count(xs,p,w)==expected,(p,w,xs,interval_count(xs,p,w),expected)
        nonzero+=bool(expected);interval_checks+=1
print('Exact rational-reconstruction comparisons:',ratio_checks,flush=True)
print('High-order interval incidence checks:',interval_checks,'nonzero:',nonzero,flush=True)
print('All final affine-spectrum and interval assertions passed.',flush=True)
