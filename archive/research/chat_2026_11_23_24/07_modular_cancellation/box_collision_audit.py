"""Exact finite box-code and collision diagnostics for toy normalized residuals."""
from itertools import product
import numpy as np
from modular_cancellation_audit import Residual, curve, chain_hist, rank_binary, parity, isprime


def boxes_audit():
    rng=np.random.default_rng(18149); decompositions=0; projected=0; fibers=0; recon=0
    for p in (13,17,37):
        a,c,pts,_=curve(p); R=Residual(p,a,c); b=R.b; mask=(1<<b)-1
        for w in range(p):
            mat=R.get(np.arange(p)[:,None],np.arange(p)[None,:],w)
            for u in range(p):
                h=np.bincount(mat[u],minlength=p)
                assert h[1:].max()<=4; fibers+=1
                for v in range(u):
                    assert np.count_nonzero((mat[u]==mat[v])&(mat[u]!=0))<=4; recon+=1
        for m in (2,3,4) if p<=17 else (2,3):
          r=m-1; width=2
          for _ in range(12):
            slots=[(int(rng.integers(p))+int(rng.integers(1,p))*np.arange(width))%p for j in range(m)]
            T=pts[int(rng.integers(1,len(pts)))][0]; hist=set(chain_hist(R,slots,T)); rebuilt=set(); boxes=[]
            for q in hist:
                cols=[q& (mask<<(b*j)) for j in range(r)]
                if all(cols):
                    boxes.append(cols)
                    for bits in product((0,1),repeat=r):
                        v=sum(x for x,t in zip(cols,bits) if t)
                        rebuilt.symmetric_difference_update((v,))
            assert rebuilt==hist
            if hist: assert len(hist)>=1<<r
            decompositions+=1
            for k in range(r+3):
              for _ in range(4):
                A=[int(rng.integers(1<<(b*r))) for j in range(k)]
                def project(v): return sum(parity(v&row)<<j for j,row in enumerate(A))
                direct=np.zeros(1<<k,np.uint8); via=np.zeros_like(direct)
                for v in hist: direct[project(v)]^=1
                for col in boxes:
                    cc=[project(v) for v in col]
                    if rank_binary(cc)<r: continue
                    for bits in product((0,1),repeat=r):
                        h=0
                        for v,t in zip(cc,bits):
                            if t: h^=v
                        via[h]^=1
                assert np.array_equal(direct,via)
                if direct.any():
                    assert k>=r and direct.sum()>=1<<r
                if k<r: assert not direct.any()
                else:
                    mob=direct.copy()
                    for i in range(k):
                        for h in range(1<<k):
                            if h>>i&1: mob[h]^=mob[h^(1<<i)]
                    assert all(not mob[h] or h.bit_count()<=k-r for h in range(1<<k))
                projected+=1
    print(f'Box code: {decompositions} exact decompositions, {projected:,} linear projections; {fibers:,} bridge fiber checks, {recon:,} fixed-output reconvergence checks',flush=True)


def collisions(p,m,w,n):
    assert isprime(p) and (4*4**3+27)%p; a=4;c=1; T=0
    rng=np.random.default_rng(34891+p+m)
    inv=np.array([0]+[pow(j,-1,p) for j in range(1,p)],np.int64)
    field=np.arange(p,dtype=np.int64); f=(field**3+a*field+c)%p
    def g(u,v,z):
        return (((u-v)**2%p*z-(u+v)*(u*v+a)-2*c)%p)**2*inv[4*f[u]*f[v]%p]*(u!=v)%p
    mids=np.indices([p]*(m-2),dtype=np.int64).reshape(m-2,-1); B0=w**m*p**(m-2)
    invalid=0;same=0;other=0;active=0;single=0; mins=1.
    for _ in range(n):
        slots=[(int(rng.integers(p))+int(rng.integers(1,p))*np.arange(w))%p for j in range(m)]
        vals=[]; same_one=0
        for xs in product(*slots):
            path=(xs[0],*mids,T); key=np.zeros(mids.shape[1],np.int64); good=np.ones(len(key),dtype=bool)
            for j in range(1,m):
                q=g(path[j-1],xs[j],path[j]); good&=q!=0; key+=q*p**(j-1)
            key=key[good]; vals.append(key)
            _,cnt=np.unique(key,return_counts=True); same_one+=int((cnt*(cnt-1)).sum())
        allq=np.concatenate(vals); _,cnt=np.unique(allq,return_counts=True); total=int((cnt*(cnt-1)).sum()); B=int((cnt&1).sum()); U=int((cnt==1).sum())
        I=B0-len(allq); assert B0-I-total<=U<=B<=B0 and total>=same_one
        invalid+=I;same+=same_one;other+=total-same_one;active+=B;single+=U;mins=min(mins,B/B0)
    denom=n*B0; lam=w**m/p; bound=4**(m-1)*lam*(1+1/w)**m
    print(f'Collision audit p={p}, m={m}, w={w}, n={n}, nominal lambda={lam:.8g}, candidate paths={denom:,}',flush=True)
    print(f'  invalid fraction={invalid/denom:.8g}; same-leaf ordered collision pairs/path={same/denom:.8g}; other-leaf ordered collision pairs/path={other/denom:.8g}',flush=True)
    print(f'  singleton fiber fraction={single/denom:.8g}; active full-support boxes/path={active/denom:.8g}; min root ratio={mins:.8g}; analytic other-leaf expectation bound={bound:.8g}',flush=True)

if __name__=='__main__':
    boxes_audit()
    for args in ((257,3,2,64),(4093,3,2,32),(65521,3,2,12),(101,4,2,20)):
        collisions(*args)
