"""Finite audits for normalized chord residual parity; toy fields only."""
from itertools import product
from math import isqrt, log, sqrt
import numpy as np


def isprime(n): return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def add(P,Q,p,a):
    if P is None: return Q
    if Q is None: return P
    x,y=P; u,v=Q
    if x==u and (y+v)%p==0: return None
    slope=((v-y)*pow((u-x)%p,-1,p) if x!=u else (3*x*x+a)*pow(2*y%p,-1,p))%p
    w=(slope*slope-x-u)%p
    return w,(slope*(x-w)-y)%p

def curve(p):
    roots=[[] for _ in range(p)]
    for y in range(p): roots[y*y%p].append(y)
    for a in range(30):
      for c in range(1,30):
        if (4*a**3+27*c*c)%p==0: continue
        N=1+sum(len(roots[(x**3+a*x+c)%p]) for x in range(p))
        if isprime(N) and N>3:
            G=next((x,roots[(x**3+a*x+c)%p][0]) for x in range(p) if roots[(x**3+a*x+c)%p])
            pts=[None]
            for _ in range(N-1): pts.append(add(pts[-1],G,p,a))
            assert len(set(pts))==N and add(pts[-1],G,p,a) is None
            return a,c,pts,roots
    raise ValueError(p)

class Residual:
    def __init__(self,p,a,c):
        self.p=p; self.a=a; self.c=c; self.b=p.bit_length()
        x=np.arange(p,dtype=np.int64); u=x[:,None]; v=x[None,:]
        f=(x*x*x+a*x+c)%p
        inv=np.array([0]+[pow(j,-1,p) for j in range(1,p)],dtype=np.int64)
        self.d=(u-v)**2%p
        self.h=((u+v)*(u*v+a)+2*c)%p
        self.fac=inv[4*f[:,None]*f[None,:]%p]*(u!=v)%p
    def get(self,u,v,w):
        p=self.p
        return ((self.d[u,v]*w-self.h[u,v])%p)**2*self.fac[u,v]%p
    def matrix(self,leaves):
        p=self.p; u=np.arange(p)[:,None,None]; w=np.arange(p)[None,None,:]
        return self.get(u,np.asarray(leaves)[None,:,None],w)

def rank_binary(rows):
    basis={}
    for x in rows:
        x=int(x)
        while x:
            i=x.bit_length()-1
            if i in basis: x^=basis[i]
            else: basis[i]=x; break
    return len(basis)

def parity(x): return int(x).bit_count()&1

def two_check_formula_audit():
    rng=np.random.default_rng(7281); cases=0; values=0
    for n in range(1,7):
      for _ in range(30):
        supp=np.flatnonzero(rng.integers(0,2,size=1<<n)).tolist()
        rows=[0]*n
        for v in supp:
            for i in range(n):
                if v>>i&1: rows[i]^=v
        r=rank_binary(rows); expect=(-1)**((len(supp)+r)&1)*(1<<(2*n-r))
        signed=0
        for a in range(1<<n):
          for b in range(1<<n):
            z=sum(1 for v in supp if not parity(a&v) and not parity(b&v))&1
            signed+=1-2*z; values+=1
        assert signed==expect; cases+=1
    print(f'Two-check exact bias identity: {cases} histograms, {values:,} ordered check pairs',flush=True)


def chain_hist(R,slots,target):
    p=R.p;m=len(slots); hist=0
    for xs in product(*slots):
      for mids in product(range(p),repeat=m-2):
        path=(xs[0],*mids,target); mask=0
        for j in range(1,m):
            mask |= int(R.get(path[j-1],xs[j],path[j])) << (R.b*(j-1))
        hist^=1<<mask
    return [i for i in range(hist.bit_length()) if hist>>i&1]

def direct_guarded(p,a,c,roots,slots,T):
    signed=[[(int(x),y) for x in S for y in roots[(int(x)**3+a*int(x)+c)%p]] for S in slots]
    count=0
    for tup in product(*signed):
        P=tup[0]
        if not P[1]: continue
        for Q in tup[1:]:
            if P is None or P[0]==Q[0] or not P[1] or not Q[1]: break
            P=add(P,Q,p,a)
            if P is None or not P[1]: break
        else: count+=P==T
    return count

def toy_audit():
    rng=np.random.default_rng(9941); rows=0; chains=0; hist_checks=0
    for p in (13,17,37,101):
        a,c,pts,roots=curve(p); R=Residual(p,a,c)
        for u,v in product(range(p),repeat=2):
            vals=R.get(u,v,np.arange(p)); h=np.bincount(vals,minlength=p)%2
            assert h[0]==1 and not h[1:].any(); rows+=1
        for m in (2,3,4):
          for w in (2,3,4):
            for _ in range(12):
                slots=[(int(rng.integers(p))+int(rng.integers(1,p))*np.arange(w))%p for j in range(m)]
                T=pts[int(rng.integers(1,len(pts)))]; vec=np.zeros(p,dtype=np.int64); vec[slots[0]]=1
                for S in slots[1:]: vec=vec@np.sum(R.matrix(S)==1,axis=1)
                assert int(vec[T[0]])==direct_guarded(p,a,c,roots,slots,T); chains+=1
                if p<=17 and m<=3:
                    hist=chain_hist(R,slots,T[0]); assert (int(vec[T[0]])&1)==int(sum(1<<(R.b*j) for j in range(m-1)) in hist)
                    if w%2==0:
                      for k in range(m):
                        for _ in range(5):
                            A=[int(rng.integers(1<<(R.b*(m-1)))) for j in range(k)]
                            out=[0]*(1<<k)
                            for v in hist: out[sum(parity(row&v)<<j for j,row in enumerate(A))]^=1
                            assert len(set(out))==1
                            if k<m-1: assert not out[0]
                            hist_checks+=len(out)
                    for _ in range(30):
                        mask=int(rng.integers(1<<(R.b*(m-1))))
                        support=[j for j in range(m-1) if (mask>>(R.b*j))&((1<<R.b)-1)]
                        if support and support!=list(range(support[0],m-1)):
                            assert sum((v&mask)==mask for v in hist)%2==0; hist_checks+=1
    print(f'Normalized chord identities: {rows:,} complete field rows; {chains} guarded chain counts; {hist_checks:,} histogram/syndrome/moment comparisons',flush=True)


def sample_masks(rng,p,pts,m,w,n,planted):
    N=len(pts); xs=np.array([-1]+[P[0] for P in pts[1:]],dtype=np.int64)
    if planted:
        idx=rng.integers(1,N,size=(n,m-1)); last=(1-idx.sum(axis=1))%N; bad=last==0
        while bad.any():
            idx[bad]=rng.integers(1,N,size=(int(bad.sum()),m-1)); last=(1-idx.sum(axis=1))%N; bad=last==0
        idx=np.column_stack((idx,last))
    else: idx=rng.integers(1,N,size=(n,m))
    t=rng.integers(1,p,size=(n,m)); s=rng.integers(0,w,size=(n,m)); a=(xs[idx]-t*s)%p
    masks=(a[:,:,None]+t[:,:,None]*np.arange(w))%p
    assert np.all(masks[np.arange(n)[:,None],np.arange(m),s]==xs[idx])
    return masks


def cross_rank3(R,masks,target):
    p=R.p; n=len(masks); b=R.b; u=np.arange(p)[None,None,None,:]
    out=np.empty(n,dtype=np.int16)
    for start in range(0,n,128):
        x=masks[start:start+128]; N=len(x)
        A=np.bitwise_xor.reduce(R.get(x[:,0,:,None,None],x[:,1,None,:,None],u).reshape(N,-1,p),axis=1)
        B=np.bitwise_xor.reduce(R.get(np.arange(p)[None,:,None],x[:,2,None,:],target),axis=2)
        for l in range(N):
            ar=A[l]; br=B[l]; rows=[int(np.bitwise_xor.reduce(br[((ar>>j)&1)!=0],initial=0)) for j in range(b)]
            out[start+l]=rank_binary(rows)
    return out


def vertex_scores(R,masks,target):
    p=R.p; n,m,w=masks.shape; out=np.empty(n,np.uint8)
    inv=np.array([0]+[pow(j,-1,p) for j in range(1,p)])
    vm=R.h*inv[R.d]%p
    for start in range(0,n,256):
        xx=masks[start:start+256]; nn=len(xx); vec=np.zeros((nn,p),dtype=np.uint8)
        vec[np.arange(nn)[:,None],xx[:,0]]=1
        for j in range(1,m):
            res=np.zeros_like(vec)
            for t in range(w):
                leaf=xx[:,j,t]; uu=np.arange(p)[None,:]; vv=leaf[:,None]
                inds=vm[uu,vv]; vals=vec.copy(); deg=uu==vv; vals[deg]=0
                np.bitwise_xor.at(res,(np.arange(nn)[:,None],inds),vals)
                res^=vec[np.arange(nn),leaf][:,None]
            vec=res
        out[start:start+nn]=vec[:,target]
    return out


def probes(p,w,n):
    a,c,pts,_=curve(p); R=Residual(p,a,c); target=pts[1][0]; rng=np.random.default_rng(90177+p)
    assert all((x*x*x+a*x+c)%p for x in range(p))
    data=[]
    for plant in (True,False):
        masks=sample_masks(rng,p,pts,3,w,n,plant)
        ranks=cross_rank3(R,masks,target)
        vertex=vertex_scores(R,masks,target)
        for x,v,r in zip(masks[:6],vertex[:6],ranks[:6]):
            vec=np.zeros(p,dtype=np.int64); vec[x[0]]=1
            for S in x[1:]: vec=vec@np.sum(R.matrix(S)==0,axis=1)%2
            assert int(vec[target])==int(v)
            A=np.bitwise_xor.reduce(R.get(x[0,:,None,None],x[1,None,:,None],np.arange(p)[None,None,:]).reshape(-1,p),axis=0)
            B=np.bitwise_xor.reduce(R.get(np.arange(p)[:,None],x[2,None,:],target),axis=1)
            rows=[0]*R.b
            for u in range(p):
                for i in range(R.b):
                    if int(A[u])>>i&1: rows[i]^=int(B[u])
            assert rank_binary(rows)==r
        means=(1-np.exp2(-2.0*ranks))/2
        data.append((ranks,means,vertex))
    gap=float(data[0][1].mean()-data[1][1].mean()); se=float(np.sqrt(data[0][1].var(ddof=1)/n+data[1][1].var(ddof=1)/n))
    vg=float(data[0][2].mean()-data[1][2].mean()); vr=sqrt(log(40)/n)
    print(f'Probe p={p}, N={len(pts)}, m=3, w={w}, n/law={n:,}, lambda_nom={w**3/p:.6f}',flush=True)
    for label,(r,means,v) in zip(('Q','R'),data):
        print(f'  {label}: rank histogram {np.bincount(r,minlength=R.b+1).tolist()}; two-check rate={means.mean():.10f}; vertex rate={v.mean():.6f}',flush=True)
    print(f'  Q-R two-check gap={gap:+.10f}, sample SE={se:.10f}, rigorous pointwise 95% half-width={vr/2:.8f}; vertex gap={vg:+.6f}, 95% half-width={vr:.8f}',flush=True)

if __name__=='__main__':
    two_check_formula_audit(); toy_audit()
    for p,w,n in ((101,4,16000),(257,4,16000),(509,6,16000),(1021,8,12000)):
        probes(p,w,n)
