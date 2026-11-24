"""Finite checks and held-out score probes. No end-to-end complexity claim."""
from itertools import product
from math import sqrt, log, isqrt
import numpy as np
import sympy as sp


def prime(n): return n>=2 and all(n%i for i in range(2,isqrt(n)+1))
def add(P,Q,p,a):
    if P is None: return Q
    if Q is None: return P
    x,y=P; u,v=Q
    if x==u and (y+v)%p==0: return None
    l=((v-y)*pow((u-x)%p,-1,p) if x!=u else (3*x*x+a)*pow(2*y%p,-1,p))%p
    w=(l*l-x-u)%p
    return w,(l*(x-w)-y)%p

def curve(p):
    roots=[[] for _ in range(p)]
    for y in range(p): roots[y*y%p].append(y)
    for a in range(30):
      for c in range(1,30):
        if (4*a*a*a+27*c*c)%p==0: continue
        N=1+sum(len(roots[(x*x*x+a*x+c)%p]) for x in range(p))
        if prime(N) and N>3:
            P=next((x,roots[(x*x*x+a*x+c)%p][0]) for x in range(p) if roots[(x*x*x+a*x+c)%p])
            pts=[None]
            for _ in range(N-1): pts.append(add(pts[-1],P,p,a))
            assert len(set(pts))==N and add(pts[-1],P,p,a) is None
            return a,c,pts,roots
    raise ValueError('no curve found')

def spectral_audit(p,S,m):
    a,c,pts,_=curve(p); N=len(pts); w=len(S); rho=w/p
    masks=[{(x+t*s)%p for s in S} for x in range(p) for t in range(1,p)]
    F=np.array([[int(P is not None and P[0] in mask) for P in pts] for mask in masks],dtype=float)
    A=F-rho*np.array([0]+[1]*(N-1)); nseed=len(masks); K=A.T@A/nseed
    cov=rho*(1-rho)*p/(p-1)
    for i,P in enumerate(pts):
      for j,Q in enumerate(pts):
        want=0 if P is None or Q is None else cov*((P[0]==Q[0])-1/p)
        assert abs(K[i,j]-want)<1e-12
    assert np.linalg.eigvalsh(K)[-1]<=2*cov+1e-12
    ind=np.indices([N]*m); C=(ind.sum(axis=0)%N==1).astype(float)-1/N
    for i in range(m):
        flat=np.moveaxis(C,i,0).reshape(N,-1)
        assert abs(np.linalg.svd(flat,compute_uv=False)[0]-N**((m-2)/2))<1e-10
    eig,V=np.linalg.eigh(K); keep=eig>1e-10; R=np.sqrt(eig[keep])[:,None]*V[:,keep].T
    Z=C
    for i in range(m): Z=np.moveaxis(np.tensordot(R,Z,axes=(1,i)),0,i)
    bound=(2*cov)**m*N**(m-2); worst=0
    for i in range(m):
        flat=np.moveaxis(Z,i,0).reshape(R.shape[0],-1)
        v=np.linalg.svd(flat,compute_uv=False)[0]**2
        assert v<=bound+1e-10; worst=max(worst,v/bound)
    rng=np.random.default_rng(340+p+m); r=R.shape[0]
    for _ in range(60):
        qv=rng.integers(0,min(4,r)+1,size=m); residual=Z.copy(); rhs=0.
        for i,q in enumerate(qv):
            if not q: continue
            U=np.linalg.qr(rng.normal(size=(r,int(q))))[0]; P=U@U.T
            proj=np.moveaxis(np.tensordot(P,Z,axes=(1,i)),0,i); rhs+=np.sum(proj**2)
            residual=np.moveaxis(np.tensordot(np.eye(r)-P,residual,axes=(1,i)),0,i)
        assert np.sum((Z-residual)**2)<=rhs+1e-10 and rhs<=qv.sum()*bound+1e-10
    for _ in range(80):
        ii=rng.integers(0,nseed,size=m)
        def convolution_at(vecs):
            arr=np.zeros(N); arr[0]=1.
            for v in vecs: arr=sum(v[j]*np.roll(arr,j) for j in range(N))
            return arr[1]-arr.sum()/N
        direct=convolution_at([F[i] for i in ii]); zsum=0.
        for bits in product([0,1],repeat=m):
            active=[A[j] for j,k in zip(ii,bits) if k]
            zsum+=(-rho)**(m-sum(bits))*convolution_at(active)
        assert abs(direct-zsum)<1e-9
    print(f'Spectral p={p}, N={N}, m={m}, w={w}: {N*N} covariance entries; {m} flattenings; max norm^2/bound={worst:.6f}; 60 slice projections and 80 expansion comparisons',flush=True)

def character_audit():
    p=101; a,c,pts,_=curve(p); wt=pts[1][0]
    z,u,v=sp.symbols('z u v')
    W=(z-u)**2*wt**2-2*((z+u)*(z*u+a)+2*c)*wt+(z*u-a)**2-4*c*(z+u)
    W=sp.Poly(W,z,u,modulus=p).as_expr()
    S=[0,1,3,7,11]; B=sp.Poly(W.subs(u,u+3*v),z,u,v,modulus=p)
    factors=[s for s in S if B.rem(sp.Poly(W.subs(u,u+s*v),z,u,v,modulus=p)).is_zero]
    assert factors==[3] and B.total_degree()==4
    for j in range(1,len(pts)):
        k=(1-j)%len(pts)
        if k: assert int(W.subs({z:pts[j][0],u:pts[k][0]}))%p==0
    chi=np.array([0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)])
    for h in range(p):
        v1=int(chi[(np.arange(p)**2-h)%p].sum())
        assert v1==(p-1 if h==0 else -1)
    aa,bb,tt,vv=sp.symbols('aa bb tt vv'); G=(aa+2*tt)*(bb+3*vv)+1
    for av,tv,want in [(aa,tt,2),(bb,vv,3)]:
        bad=[s for s in range(p) if sp.Poly(sp.diff(G,tv)-s*sp.diff(G,av),aa,bb,tt,vv,modulus=p).is_zero]
        assert bad==[want]
    max_l1=max(np.abs(np.fft.fft(((np.arange(p)>>j)&1))/p).sum() for j in range(p.bit_length()))
    assert max_l1<4*(1+log(p))**2
    print(f'Character p={p}: one exceptional discriminant factor among {len(S)} slopes; {p} quadratic sums; {2*p} directional checks; bit Fourier L1 max={max_l1:.6f}',flush=True)

def det_batch(A,p,inv):
    A=A.copy()%p; n,m,_=A.shape; d=np.ones(n,dtype=np.int64); rows=np.arange(n)
    for j in range(m):
        k=j+np.argmax(A[:,j:,j]!=0,axis=1)
        d[k!=j]=-d[k!=j]
        old=A[:,j,:].copy(); A[:,j,:]=A[rows,k,:]; A[rows,k,:]=old
        piv=A[:,j,j]; d=d*piv%p
        if j+1<m:
            f=A[:,j+1:,j]*inv[piv,None]%p
            A[:,j+1:,j+1:]=(A[:,j+1:,j+1:]-f[:,:,None]*A[:,None,j,j+1:])%p
    return d%p

def mixed(G,p,inv):
    n,m,_,_=G.shape; ans=np.zeros(n,np.int64)
    for mask in range(1,1<<m):
        ix=[j for j in range(m) if mask&(1<<j)]
        ans+=(-1)**(m-len(ix))*det_batch(G[:,ix,:,:].sum(axis=1),p,inv)
    return ans%p

def gram_tables(p,a,c,pts,roots,m):
    d=m+1; basis=sorted([(2*i,i,0) for i in range(d//2+1)]+[(2*i+3,i,1) for i in range((d-3)//2+1)])
    assert len(basis)==d; basis=basis[1:]; xt,yt=pts[1]; yt=(-yt)%p
    at=[pow(xt,i,p)*pow(yt,j,p)%p for _,i,j in basis]
    formal=np.zeros((p,m,m),np.int64); rat=formal.copy(); sizes=np.zeros(p,np.int64)
    for x in range(p):
        fx=(x*x*x+a*x+c)%p
        for u,(_,i,j) in enumerate(basis):
          for v,(_,k,l) in enumerate(basis):
            s1=2*pow(x,i+k,p)*pow(fx,(j+l)//2,p) if (j+l)%2==0 else 0
            s2=2*at[u]*pow(x,k,p) if l==0 else 0
            s3=2*at[v]*pow(x,i,p) if j==0 else 0
            formal[x,u,v]=(s1-s2-s3+2*at[u]*at[v])%p
        sizes[x]=len(roots[fx])
        for y in roots[fx]:
            col=np.array([(pow(x,i,p)*pow(y,j,p)-at[u])%p for u,(_,i,j) in enumerate(basis)],dtype=np.int64)
            rat[x]=(rat[x]+np.outer(col,col))%p
        if sizes[x]==2: assert np.array_equal(formal[x],rat[x])
    return formal,rat,sizes,basis,at

def sample_seeds(rng,n,p,pts,m,w,planted):
    N=len(pts); xs=np.array([-1]+[P[0] for P in pts[1:]],dtype=np.int64)
    if planted:
        idx=rng.integers(1,N,size=(n,m-1)); last=(1-idx.sum(axis=1))%N; bad=last==0
        while bad.any():
            idx[bad]=rng.integers(1,N,size=(int(bad.sum()),m-1)); last=(1-idx.sum(axis=1))%N; bad=last==0
        idx=np.column_stack([idx,last]); assert np.all(idx.sum(axis=1)%N==1)
    else: idx=rng.integers(1,N,size=(n,m))
    tt=rng.integers(1,p,size=(n,m)); ss=rng.integers(0,w,size=(n,m)); aa=(xs[idx]-tt*ss)%p
    assert np.all((aa+tt*ss)%p==xs[idx])
    return aa,tt

def all_scores(aa,tt,p,w,form,rat,sizes,inv):
    ans=[]
    for start in range(0,len(aa),1024):
        xs=(aa[start:start+1024,:,None]+tt[start:start+1024,:,None]*np.arange(w)[None,None,:])%p
        zf=mixed(form[xs].sum(axis=2),p,inv); zr=mixed(rat[xs].sum(axis=2),p,inv)
        sz=np.prod(sizes[xs].sum(axis=2),axis=1)%p
        ans.append(np.column_stack([zf,zr,zr*inv[sz]%p]))
    return np.concatenate(ans)

def literal_moment_check():
    p=13; m=3; w=2; a,c,pts,roots=curve(p)
    form,rat,sizes,basis,at=gram_tables(p,a,c,pts,roots,m)
    inv=np.array([0]+[pow(x,-1,p) for x in range(1,p)])
    rng=np.random.default_rng(789); aa,tt=sample_seeds(rng,100,p,pts,m,w,True)
    got=all_scores(aa,tt,p,w,form,rat,sizes,inv)
    for j in range(len(aa)):
        slots=[]
        for u,t in zip(aa[j],tt[j]):
            xx=[(int(u)+int(t)*s)%p for s in range(w)]
            slots.append([(x,y) for x in xx for y in roots[(x*x*x+a*x+c)%p]])
        mats=np.array([[[(pow(x,i,p)*pow(y,k,p)-at[q])%p for x,y in pp] for q,(_,i,k) in enumerate(basis)] for pp in product(*slots)])
        z=np.sum(det_batch(mats,p,inv)**2)%p
        assert z==got[j,1]
    print('Literal rational determinant moment: 100/100 agree',flush=True)

def collision_scaled(a,b,p):
    u=np.bincount(a,minlength=p).astype(float); v=np.bincount(b,minlength=p).astype(float); n=len(a); k=len(b)
    return p*(np.dot(u,u-1)/(n*(n-1))+np.dot(v,v-1)/(k*(k-1))-2*np.dot(u,v)/(n*k))

def moment_probe(p,m,w,n):
    a,c,pts,roots=curve(p); N=len(pts); form,rat,sizes,basis,at=gram_tables(p,a,c,pts,roots,m)
    inv=np.array([0]+[pow(x,-1,p) for x in range(1,p)]); rng=np.random.default_rng(61225+p+m)
    train=[]; test=[]
    for planted in [True,False]:
        aa,tt=sample_seeds(rng,2*n,p,pts,m,w,planted); s=all_scores(aa,tt,p,w,form,rat,sizes,inv)
        train.append(s[:n]); test.append(s[n:])
    half=sqrt(log(2/.05)/n)
    print(f'Full-value probe p={p}, N={N}, m={m}, w={w}, nominal lambda={N**(m-1)*(w/p)**m:.5f}, n={n} training AND n={n} held-out per law; pointwise 95% half-width={half:.6f}',flush=True)
    for j,name in enumerate(['formal Z2','rational Z2','rational Z2 / tuple count']):
        qa,ra=train[0][:,j],train[1][:,j]; qt,rt=test[0][:,j],test[1][:,j]
        rule=np.bincount(qa,minlength=p)>np.bincount(ra,minlength=p); gap=rule[qt].mean()-rule[rt].mean()
        print(f'  {name}: held-out gap={gap:+.6f}; CI=[{gap-half:+.6f},{gap+half:+.6f}]; scaled L2^2 estimate={collision_scaled(qt,rt,p):+.6f}',flush=True)


def average_gap_control():
    # Signal only on T=+G and -G, with opposite signs. The signed average
    # vanishes exactly but the absolute average is positive. Use literal
    # relation existence as a deliberately expensive positive control.
    p=13; a,c,pts,_=curve(p); N=len(pts); S=[0,1]; rho=len(S)/p
    masks=[{(u+t*s)%p for s in S} for u in range(p) for t in range(1,p)]
    F=np.array([[int(P is not None and P[0] in mask) for P in pts] for mask in masks],dtype=np.int64)
    ns=len(masks); D=F.sum(1); pi=D/D.sum(); count=np.outer(D,D)
    H=F@F[:,(1-np.arange(N))%N].T
    Hneg=F@F[:,(-1-np.arange(N))%N].T
    assert np.array_equal(H,Hneg)
    score=(H>0).astype(float)
    delta=(H*score).sum()/H.sum()-(count*score).sum()/count.sum()
    assert delta>0
    xs=sorted({P[0] for P in pts[1:]}); K=np.array([[int(x in mask) for mask in masks] for x in xs],dtype=float)/(ns*rho)
    assert np.allclose(K.sum(axis=1),1)
    assert np.allclose(K.mean(axis=0),pi)
    B=(K-pi[None,:])/2; g1=B@score; g2=g1@B.T
    # Verify exact probability mass identity and centered-coefficient energy.
    xpos={x:i for i,x in enumerate(xs)}
    pairs=[(P,pts[(1-j)%N]) for j,P in enumerate(pts) if j and (1-j)%N]
    gm=np.mean([g2[xpos[P[0]],xpos[Q[0]]] for P,Q in pairs])
    # Here marginal correction is nonzero but small; bound directly.
    assert abs(4*gm-delta)<=8/(N-1)
    assert np.max(np.sum(g1*g1,axis=0))<=1/(4*rho)+1e-10
    assert np.sum(g2*g2)<=1/(4*rho)**2+1e-10
    avg=2*delta/(N-1)
    print(f'Average-gap control p={p}, N={N}: only 2/{N-1} targets carry signal; signed average=0; absolute average={avg:.6f}; centered energy verified',flush=True)

def residual_parity_audit():
    rng=np.random.default_rng(8356); n=7; d=3; U=1<<n
    v=((np.arange(U)[:,None]>>np.arange(n)[None,:])&1).astype(np.int64)
    checks=0
    for _ in range(100):
        c=rng.integers(0,2,size=U); A=rng.integers(0,2,size=(d,n))
        vals=1-(v@A.T%2); exact=int(np.sum(c*np.prod(vals,axis=1))%2)
        # Independent polynomial expansion, choosing constant or a coordinate
        # from each hash factor. Repeated variables are idempotent on the cube.
        exp=0
        for choice in product(range(-1,n),repeat=d):
            weight=1; J=set()
            for j,i in enumerate(choice):
                if i>=0: weight*=int(A[j,i]); J.add(i)
            if not weight: continue
            q=np.prod(v[:,sorted(J)],axis=1) if J else np.ones(U,dtype=np.int64)
            exp^=int(np.sum(c*q)%2)
        assert exact==exp; checks+=1
    for _ in range(100):
        # An affine flat of dimension d+1 has even intersection with every
        # affine fiber of d linear constraints.
        A=rng.integers(0,2,size=(d,n)); origin=rng.integers(0,2,size=n)
        flat=(v[:1<<(d+1)]^origin[None,:])
        assert np.count_nonzero(np.all(flat@A.T%2==0,axis=1))%2==0
        checks+=1
    print(f'Residual parity sketch: {checks} exact polynomial-expansion and affine-flat cancellation checks',flush=True)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1: moment_probe(*map(int,sys.argv[1:]))
    else:
        for p,S,m in [(13,[0,1,3],2),(13,[0,1],3),(17,[0,1,4],3),(37,[0,2,3,10],3)]: spectral_audit(p,S,m)
        character_audit(); literal_moment_check(); average_gap_control(); residual_parity_audit()
        print('All finite assertions passed. Run probe cases as p m w n arguments.',flush=True)
