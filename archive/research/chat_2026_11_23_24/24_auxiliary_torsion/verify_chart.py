#!/usr/bin/env python3
"""Exact homogeneous-chart and compact circuit checks for auxiliary-torsion fibers.
Determinants here use toy Gaussian elimination / subset DP. Production division-free
polynomial-size determinants can use Berkowitz; no collective trace is evaluated.
"""
from collections import Counter
from itertools import product
from random import Random
from previous_fiber_verifier import Ring,roots_poly,pe,interp,res,s3
from verify import Quotient,division,aux,FIX,Scalar,sum_poly,projection_poly
rng=Random(2026112422)
class K2:
    def __init__(self,p):
        self.p=p;self.nr=next(x for x in range(2,p) if pow(x,(p-1)//2,p)==p-1);self.zero=(0,0);self.one=(1,0)
    def c(self,n):return (n%self.p,0)
    def add(self,a,b):return ((a[0]+b[0])%self.p,(a[1]+b[1])%self.p)
    def neg(self,a):return (-a[0]%self.p,-a[1]%self.p)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):return ((a[0]*b[0]+self.nr*a[1]*b[1])%self.p,(a[0]*b[1]+a[1]*b[0])%self.p)
    def pw(self,a,n):
        x=self.one
        while n:
            if n&1:x=self.mul(x,a)
            n//=2
            if n:a=self.mul(a,a)
        return x
    def inv(self,a):
        assert a!=self.zero
        t=pow((a[0]*a[0]-self.nr*a[1]*a[1])%self.p,-1,self.p);return (a[0]*t%self.p,-a[1]*t%self.p)
class ProductK(Ring):
    def __init__(self,p,n):self.K=K2(p);self.p=p;self.n=n;self.z=(self.K.zero,)*n;self.o=(self.K.one,)*n
    def c(self,a):return (self.K.c(a),)*self.n
    def const(self,a):return (a,)*self.n
    def add(self,a,b):return tuple(self.K.add(x,y) for x,y in zip(a,b))
    def neg(self,a):return tuple(self.K.neg(x) for x in a)
    def mul(self,a,b):return tuple(self.K.mul(x,y) for x,y in zip(a,b))
    def pw(self,a,k):return tuple(self.K.pw(x,k) for x in a)

def homog(n,X,Z,a,c,R):
    from functools import lru_cache
    def mul(*args):
        out=R.one
        for v in args:out=R.mul(out,v)
        return out
    def add(*args):
        out=R.zero
        for v in args:out=R.add(out,v)
        return out
    def sq(x):return R.mul(x,x)
    def cube(x):return R.mul(sq(x),x)
    xs=[R.one];zs=[R.one]
    for _ in range(6):xs.append(R.mul(xs[-1],X));zs.append(R.mul(zs[-1],Z))
    D=R.mul(R.c(4),add(xs[3],mul(R.c(a),X,zs[2]),mul(R.c(c),zs[3])));D2=sq(D)
    base={0:R.zero,1:R.one,2:R.one,
          3:add(mul(R.c(3),xs[4]),mul(R.c(6*a),xs[2],zs[2]),mul(R.c(12*c),X,zs[3]),mul(R.c(-a*a),zs[4])),
          4:mul(R.c(2),add(xs[6],mul(R.c(5*a),xs[4],zs[2]),mul(R.c(20*c),xs[3],zs[3]),mul(R.c(-5*a*a),xs[2],zs[4]),mul(R.c(-4*a*c),X,zs[5]),mul(R.c(-8*c*c-a**3),zs[6])))}
    @lru_cache(None)
    def phi(k):
        if k<=4:return base[k]
        j=k//2
        if not k%2:return mul(phi(j),R.sub(mul(phi(j+2),sq(phi(j-1))),mul(phi(j-2),sq(phi(j+1)))))
        L=mul(phi(j+2),cube(phi(j)));Q=mul(phi(j-1),cube(phi(j+1)))
        return R.sub(mul(D2,L) if not j%2 else L,mul(D2,Q) if j%2 else Q)
    return phi(n)

def det(K,M):
    M=[list(row) for row in M];n=len(M);v=K.one
    for j in range(n):
        i=next((i for i in range(j,n) if M[i][j]!=K.zero),None)
        if i is None:return K.zero
        if i!=j:M[i],M[j]=M[j],M[i];v=K.neg(v)
        x=M[j][j];v=K.mul(v,x);iv=K.inv(x)
        for k in range(j+1,n):
            t=K.mul(M[k][j],iv)
            for l in range(j+1,n):M[k][l]=K.sub(M[k][l],K.mul(t,M[j][l]))
    return v

def chart_matrices(p,ell,a,c,fs,d,alpha,flt):
    R=ProductK(p,len(fs));K=R.K
    # F'(V)=V^d F(alpha+1/V).  Preserve declared degree, including roots at infinity.
    from math import comb
    o=[R.z]*(d+1);F=list(zip(*[f+[0]*(d+1-len(f)) for f in fs]));V=[K.one]
    for _ in range(d):V.append(K.mul(V[-1],alpha))
    for j,f in enumerate(F):
        for k in range(j+1):o[d-j+k]=R.add(o[d-j+k],R.mul(tuple(K.c(x*comb(j,k)) for x in f),R.const(V[k])))
    assert all(x!=K.zero for x in o[-1]);iv=R.pw(o[-1],p*p-2);o=[R.mul(x,iv) for x in o];assert o[-1]==R.o
    Q=Quotient(R,o);Z=Q.x;A0=(R.const(alpha),)+(R.z,)*(d-1);X=Q.add(Q.one,Q.mul(A0,Z))
    psi=homog(ell,X,Z,a,c,Q)
    rational=Q.sub(Q.pw(X,p),Q.mul(X,Q.pw(Z,p-1)))
    A,t,aa,cc=flt;L=Q.add(Q.mul(Q.c(A),Z),Q.mul(Q.c(t),X));Z2=Q.mul(Z,Z)
    cubic=Q.add(Q.add(Q.mul(Q.mul(L,L),L),Q.mul(Q.mul(Q.c(aa),L),Z2)),Q.mul(Q.c(cc),Q.mul(Z,Z2)))
    lift=Q.sub(Q.pw(cubic,(p-1)//2),Q.pw(Z,3*((p-1)//2)))
    matrices=[]
    for f in (psi,rational,lift):
        cols=[];q=Q.one
        for j in range(d):cols.append(Q.mul(f,q));q=Q.mul(q,Z)
        matrices.append([[[cols[j][i][z] for j in range(d)] for i in range(d)] for z in range(len(fs))])
    return R.K,matrices

def homogeneous_tests():
    C=Counter()
    for p,ell in FIX:
        a,c,N,H,S,tr=aux(p,ell);K=K2(p)
        for _ in range(60):
            n=rng.choice((ell,3,4,7,32,127,129));X=(rng.randrange(p),rng.randrange(p));Z=(rng.randrange(p),rng.randrange(p));dd=(n*n-(1 if n%2 else 4))//2
            got=homog(n,X,Z,a,c,K)
            if Z!=K.zero:
                ref=division(n,K.mul(X,K.inv(Z)),a,c,K)[0];assert got==K.mul(K.pw(Z,dd),ref);C['homogeneous']+=1
            if n%2:
                assert homog(n,X,K.zero,a,c,K)==K.mul(K.c(n),K.pw(X,dd));C['infinity']+=1
    print('Homogeneous recurrence:',dict(C),flush=True)

def chart_tests():
    C=Counter()
    for _ in range(80):
        p,ell=rng.choice(FIX[:5]);a,c,N,H,S,tr=aux(p,ell);d=rng.randrange(1,8);n=rng.randrange(2,7)
        K=K2(p);fs=[]
        for j in range(n):
            deg=rng.randrange(d+1);roots=[rng.choice(S) if rng.randrange(2) else rng.randrange(p) for _ in range(deg)];f=roots_poly(p,roots)
            fs.append(f);C['degree_drop']+=deg<d;C['repeated']+=len(set(roots))<len(roots)
        # Choose a public alpha; the theorem charges the small bad-chart probability.
        while True:
            alpha=(rng.randrange(p),rng.randrange(1,p));
            def evalk(f,x):
                y=K.zero
                for v in reversed(f):y=K.add(K.mul(y,x),K.c(v))
                return y
            if all(evalk(f,alpha)!=K.zero for f in fs):break
            C['rejected_chart']+=1
        A,t=rng.randrange(p),rng.randrange(1,p);aa,cc=a,c;flt=(A,t,aa,cc)
        truth=[any(pe(f,s,p)==0 and pow(((A+t*s)**3+aa*(A+t*s)+cc)%p,(p-1)//2,p)==1 for s in S) for f in fs]
        K,Ms=chart_matrices(p,ell,a,c,fs,d,alpha,flt)
        for _ in range(12):
            eta=(rng.randrange(p),rng.randrange(p));eta2=K.mul(eta,eta)
            for z in range(n):
                M=[[K.add(Ms[0][z][i][j],K.add(K.mul(eta,Ms[1][z][i][j]),K.mul(eta2,Ms[2][z][i][j]))) for j in range(d)] for i in range(d)]
                hit=det(K,M)==K.zero;assert not truth[z] or hit;C['determinants']+=1;C['true_positives']+=truth[z];C['false_positives']+=hit and not truth[z]
        C['fibers']+=n
    # Exhaust all eta in a small extension, including roots at infinity and degree drops.
    p,ell=13,5;a,c,N,H,S,tr=aux(p,ell);K=K2(p)
    for _ in range(25):
        d=rng.randrange(1,5);f=roots_poly(p,[rng.randrange(p) for _ in range(rng.randrange(d+1))]);alpha=(rng.randrange(p),1)
        flt=(rng.randrange(p),rng.randrange(1,p),a,c);A,t,aa,cc=flt
        truth=any(pe(f,s,p)==0 and pow(((A+t*s)**3+aa*(A+t*s)+cc)%p,(p-1)//2,p)==1 for s in S)
        K,Ms=chart_matrices(p,ell,a,c,[f],d,alpha,flt);bad=0
        for eta in product(range(p),repeat=2):
            eta2=K.mul(eta,eta);M=[[K.add(Ms[0][0][i][j],K.add(K.mul(eta,Ms[1][0][i][j]),K.mul(eta2,Ms[2][0][i][j]))) for j in range(d)] for i in range(d)]
            hit=det(K,M)==K.zero;assert not truth or hit;bad+=hit;C['exhaustive_eta']+=1
        if not truth:assert bad<=2*d;C['root_bound_checks']+=1
    assert C['true_positives'] and C['root_bound_checks']
    print('One-determinant chart checks:',dict(C),flush=True)

def det_ring(R,M):
    # Division-free subset determinant for tiny test matrices.
    D={0:R.o};n=len(M)
    for row in M:
        E={}
        for mask,a in D.items():
            for j,b in enumerate(row):
                if mask>>j&1:continue
                x=R.mul(a,b)
                if (mask>>(j+1)).bit_count()&1:x=R.neg(x)
                k=mask|1<<j;E[k]=R.add(E.get(k,R.z),x)
        D=E
    return D[(1<<n)-1]
def res_ring(R,f,g):
    d,e=len(f)-1,len(g)-1;M=[[R.z]*(d+e) for _ in range(d+e)]
    for j in range(e):M[j][j:j+d+1]=f[::-1]
    for j in range(d):M[e+j][j:j+e+1]=g[::-1]
    return det_ring(R,M)
def s3_ring(R,u,v,a,c):
    z=R.mul(u,v);s=R.add(u,v);q=R.sub(z,R.c(a));e=R.sub(u,v)
    return [R.sub(R.mul(q,q),R.mul(R.c(4*c),s)),R.mul(R.c(-2),R.add(R.mul(s,R.add(z,R.c(a))),R.c(2*c))),R.mul(e,e)]
def interp_ring(R,vals):
    n=len(vals);D=vals[:];co=[R.z]*n;basis=[1];fac=1;p=R.p
    for j in range(n):
        k=R.mul(D[0],R.c(pow(fac,-1,p)))
        for i,v in enumerate(basis):co[i]=R.add(co[i],R.mul(k,R.c(v)))
        D=[R.sub(v,u) for u,v in zip(D,D[1:])]
        if j+1<n:
            b=[0]*(len(basis)+1)
            for i,v in enumerate(basis):b[i]=(b[i]-j*v)%p;b[i+1]=(b[i+1]+v)%p
            basis=b;fac=fac*(j+1)%p
    return co
def partial(R,xs,a,c):
    F=[R.neg(xs[0]),R.o]
    for x in xs[1:]:
        d=2*(len(F)-1);F=interp_ring(R,[res_ring(R,F,s3_ring(R,x,R.c(j),a,c)) for j in range(d+1)])
    return F

def circuit_tests():
    C=Counter()
    for _ in range(40):
        p=rng.choice((101,127,257));a,c=rng.randrange(p),rng.randrange(p)
        if (4*a**3+27*c*c)%p==0:continue
        m=rng.choice((2,3));h=rng.choice((1,2));d=h*h*2**(m-1);T=rng.randrange(p);A,t=rng.randrange(p),rng.randrange(1,p);z=[tuple(rng.randrange(p) for _ in range(m-1)) for j in range(3)]
        R=Ring(p,len(z));xs=list(zip(*z));Frest=partial(R,xs,a,c)
        P=[R.c(x) for x in projection_poly(p,a,c,h,T)];vals=[]
        for u in range(d+1):
            x=R.c(A+t*u)
            if h==1:vals.append(res_ring(R,Frest,s3_ring(R,x,R.c(T),a,c)))
            else:
                dv=2*(len(Frest)-1);F=interp_ring(R,[res_ring(R,Frest,s3_ring(R,x,R.c(j),a,c)) for j in range(dv+1)])
                vals.append(res_ring(R,F,P))
        F=interp_ring(R,vals)
        for j,z0 in enumerate(z):
            for u in (0,1,2,3,17,29):
                x=(A+t*u)%p
                ref=res(sum_poly(z0+(x,),p,a,c),projection_poly(p,a,c,h,T),p)
                assert pe([v[j] for v in F],u,p)==ref;C['values']+=1
        C['coefficient_lanes']+=len(z)*(d+1);C['fixtures']+=1
    print('Ring-valued summation polynomial construction:',dict(C),flush=True)
if __name__=='__main__':
    assert __debug__
    homogeneous_tests();chart_tests();circuit_tests()
    print('All assertions passed. No collective trace was computed.')
