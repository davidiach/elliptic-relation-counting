#!/usr/bin/env python3
"""Exact tests of smooth isogeny kernel circuits, not a collective evaluator."""
from random import Random
from collections import Counter
from math import prod,comb
from functools import lru_cache
rng=Random(202611242043); C=Counter()

def points(p,a,c):
    ys={}
    for y in range(p):ys.setdefault(y*y%p,[]).append(y)
    return [(x,y) for x in range(p) for y in ys.get((x*x*x+a*x+c)%p,())]
def add(P,Q,M,a):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if (y+v)%M==0:return None
        assert y==v
        t=(3*x*x+a)*pow(2*y,-1,M)%M
    else:t=(v-y)*pow(u-x,-1,M)%M
    z=(t*t-x-u)%M;return z,(t*(x-z)-y)%M
def mul(P,n,M,a):
    Q=None
    while n:
        if n&1:Q=add(Q,P,M,a)
        n//=2
        if n:P=add(P,P,M,a)
    return Q

def pe(f,x,M):
    out=0
    for v in reversed(f):out=(out*x+v)%M
    return out
def plus(f,g,M):
    h=[0]*max(len(f),len(g))
    for i,v in enumerate(f):h[i]+=v
    for i,v in enumerate(g):h[i]+=v
    return [v%M for v in h]
def scale(f,c,M):return [v*c%M for v in f]
def pm(f,g,M):
    h=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        for j,b in enumerate(g):h[i+j]=(h[i+j]+a*b)%M
    return h
def roots_poly(roots,M):
    h=[1]
    for x in roots:h=pm(h,[-x%M,1],M)
    return h
def rem(f,g,M):
    d=len(g)-1;assert g[-1]==1
    h=list(f)+[0]*max(0,d-len(f))
    for j in range(len(h)-1,d-1,-1):
        v=h[j]
        for k in range(d):h[j-d+k]=(h[j-d+k]-v*g[k])%M
    return h[:d]
def deriv(f,M):return [(i*v)%M for i,v in enumerate(f) if i]

def velu(a,c,H,M,ell):
    us=[mul(H,j,M,a)[0] for j in range(1,(ell+1)//2)]
    assert len(set(us))==len(us) and mul(H,ell,M,a) is None
    B=roots_poly(us,M);A=[0]+pm(B,B,M)
    for i,u in enumerate(us):
        b=roots_poly(us[:i]+us[i+1:],M);b2=pm(b,b,M)
        f=(u**3+a*u+c)%M;fp=(3*u*u+a)%M
        A=plus(A,scale(pm([-u%M,1],b2,M),2*fp,M),M)
        A=plus(A,scale(b2,4*f,M),M)
    assert len(A)==ell+1 and A[-1]==B[-1]==1
    aa=(a-10*sum(3*u*u+a for u in us))%M
    cc=(c-14*sum(5*u**3+3*a*u+2*c for u in us))%M
    return A,B,aa,cc

def image(P,A,B,M):
    if P is None:return None
    x,y=P;u=pe(A,x,M);v=pe(B,x,M)
    if v==0:return None
    iv=pow(v,-1,M);xx=u*iv*iv%M
    yy=y*(pe(deriv(A,M),x,M)*v-2*u*pe(deriv(B,M),x,M))*iv**3%M
    return xx,yy

def tower(a,c,H,M,ell,k):
    out=[]
    for j in range(k):
        K=mul(H,ell**(k-j-1),M,a)
        A,B,aa,cc=velu(a,c,K,M,ell)
        out.append((A,B));H=image(H,A,B,M)
        if H is not None:assert (H[1]**2-H[0]**3-aa*H[0]-cc)%M==0
        a,c=aa,cc
    assert H is None
    return out,a,c

class Ring:
    def __init__(self,M,g=None):self.M=M;self.g=g;self.n=1 if g is None else len(g)-1
    def c(self,x):return (x%self.M,)+(0,)*(self.n-1)
    def add(self,a,b):return tuple((x+y)%self.M for x,y in zip(a,b))
    def mul(self,a,b):
        if self.g is None:return (a[0]*b[0]%self.M,)
        return tuple(rem(pm(a,b,self.M),self.g,self.M))
    def hom(self,f,X,Z):
        v=self.c(f[-1]);power=self.c(1)
        for a in f[-2::-1]:
            power=self.mul(power,Z);v=self.add(self.mul(v,X),self.mul(self.c(a),power))
        return v

def kernel_eval(stages,R,U,V):
    X,Z,K=U,V,R.c(1)
    for A,B in stages:
        Bv=R.hom(B,X,Z);X=R.hom(A,X,Z)
        Z=R.mul(Z,R.mul(Bv,Bv));K=R.mul(K,Bv)
        assert Z==R.mul(V,R.mul(K,K))
    return K,X,Z

class Dual:
    """First derivative propagated through the homogeneous kernel circuit."""
    def __init__(self,R):self.R=R
    def c(self,x):return self.R.c(x),self.R.c(0)
    def add(self,a,b):return self.R.add(a[0],b[0]),self.R.add(a[1],b[1])
    def mul(self,a,b):
        R=self.R
        return R.mul(a[0],b[0]),R.add(R.mul(a[0],b[1]),R.mul(a[1],b[0]))
    def hom(self,f,X,Z):
        v=self.c(f[-1]);power=self.c(1)
        for a in f[-2::-1]:
            power=self.mul(power,Z);v=self.add(self.mul(v,X),self.mul(self.c(a),power))
        return v

def unit_inverse(q,g,M,p):
    d=len(g)-1
    A=[list(v)+[int(i==0)] for i,v in enumerate(zip(*(rem([0]*j+list(q),g,M) for j in range(d))))]
    for i in range(d):
        j=next(j for j in range(i,d) if A[j][i]%p)
        A[i],A[j]=A[j],A[i];iv=pow(A[i][i],-1,M);A[i]=[v*iv%M for v in A[i]]
        for j in range(d):
            if j==i:continue
            c=A[j][i];A[j]=[(x-c*y)%M for x,y in zip(A[j],A[i])]
    return tuple(row[-1] for row in A)

def rational_trace(stages,g,P,M,p):
    R=Ring(M,g);x=tuple(rem([0,1],g,M));D=Dual(R)
    (Q,Qp),_,_=kernel_eval(stages,D,(x,R.c(1)),D.c(1))
    v=R.mul(tuple(P),R.mul(Qp,unit_inverse(Q,g,M,p)))
    return -v[-1]%M

def power_sums(stages,jmax,M,w):
    R=Ring(M,[0]*(jmax+1)+[1]);t=tuple(rem([0,1],R.g,M))
    Q,_,_=kernel_eval(stages,R,R.c(1),t)
    assert Q[0]==1
    P=[w%M]
    for j in range(1,jmax+1):P.append((-j*Q[j]-sum(Q[k]*P[j-k] for k in range(1,j)))%M)
    return P

def det(A,M):
    n=len(A);dp={0:1}
    for i in range(n):
        nxt={}
        for S,v in dp.items():
            for j in range(n):
                if S>>j&1:continue
                s=S|1<<j;sign=(-1)**((S>>(j+1)).bit_count())
                nxt[s]=(nxt.get(s,0)+sign*v*A[i][j])%M
        dp=nxt
    return dp[(1<<n)-1]%M

def muldet(q,g,M):
    d=len(g)-1;cols=[rem([0]*j+list(q),g,M) for j in range(d)]
    return det(list(zip(*cols)),M)

def norm(stages,L,f,M,alpha):
    d=len(f)-1;w=(L-1)//2
    if d==0:return pow(f[0],w,M)
    g=[0]*(d+1)
    for j,v in enumerate(f):
        for i in range(j+1):g[d-j+i]=(g[d-j+i]+v*comb(j,i)*alpha**i)%M
    l=g[-1];g=[v*pow(l,-1,M)%M for v in g]
    R=Ring(M,g);X=tuple(rem([0,1],g,M));U=R.add(R.c(1),R.mul(R.c(alpha),X))
    q,_,_=kernel_eval(stages,R,U,X)
    return pow(l,w,M)*muldet(q,g,M)%M

@lru_cache(None)
def fixture(p,ell,k):
    L=ell**k
    for tries in range(1,4000):
        a,c=rng.randrange(p),rng.randrange(p)
        if (4*a**3+27*c*c)%p==0:continue
        ps=points(p,a,c);N=len(ps)+1
        if N%L:continue
        v=N
        while v%ell==0:v//=ell
        for P in ps[:min(100,len(ps))]:
            H=mul(P,v,p,a)
            if H is None:continue
            while mul(H,L,p,a) is not None:H=mul(H,ell,p,a)
            if H is not None and mul(H,L//ell,p,a) is not None:return a,c,H,N,tries
    raise AssertionError((p,ell,k))

def psi(n,x,a,c,M):
    D=4*(x**3+a*x+c)%M;D2=D*D%M
    base=[0,1,1,(3*x**4+6*a*x*x+12*c*x-a*a)%M,2*(x**6+5*a*x**4+20*c*x**3-5*a*a*x*x-4*a*c*x-8*c*c-a**3)%M]
    @lru_cache(None)
    def t(k):
        if k<5:return base[k]
        j=k//2
        if k%2==0:return t(j)*(t(j+2)*t(j-1)**2-t(j-2)*t(j+1)**2)%M
        x0=t(j+2)*t(j)**3;x1=t(j-1)*t(j+1)**3
        return ((x0*D2 if j%2==0 else x0)-(x1*D2 if j%2 else x1))%M
    return t(n)
def lift(H,L,p,e,a,c):
    x,y=H;step=p
    for j in range(1,e):
        M=step*p;q=psi(L,x,a,c,M);assert q%step==0
        dq=((psi(L,x+step,a,c,M)-q)//step)%p;assert dq
        x+=step*((-q//step)*pow(dq,-1,p)%p)
        u=(x**3+a*x+c-y*y)%M;assert u%step==0
        y+=step*(u//step*pow(2*y,-1,p)%p);step=M
    assert psi(L,x,a,c,step)==0 and (y*y-x**3-a*x-c)%step==0
    assert mul((x,y),L,step,a) is None
    return x,y

def run():
    configs=[(101,7,1),(101,7,2),(257,7,2),(509,7,2),(1009,7,2),(347,7,3),(127,3,3),(257,5,2)]
    for p,ell,k in configs:
        a,c,H,N,tr=fixture(p,ell,k);L=ell**k
        print('fixture',p,ell,k,'#E',N,'curve',a,c,'H',H,'proposals',tr,flush=True)
        for e in (1,2,3):
            M=p**e;J=H if e==1 else lift(H,L,p,e,a,c)
            stages,aa,cc=tower(a,c,J,M,ell,k)
            roots=[mul(J,j,M,a)[0] for j in range(1,(L+1)//2)]
            Q=roots_poly(roots,M);assert len(set(x%p for x in roots))==len(roots);R=Ring(M)
            for x in list(roots)+[rng.randrange(M) for _ in range(30)]:
                q,X,Z=kernel_eval(stages,R,(x,),(1,));assert q[0]==pe(Q,x,M);C['kernel evaluations']+=1
            for j in range(40):
                d=rng.randrange(1,7);f=[rng.randrange(M) for _ in range(d+1)]
                if j%4==0:f[-1]=0;C['declared degree drops']+=1
                if j%3==0:f=pm(f[:-1],[-rng.choice(roots)%M,1],M);C['planted zeros']+=1
                if not any(v%p for v in f):f[0]=1
                alpha=next(x for x in range(p) if pe(f,x,p));true=prod(pe(f,x,M) for x in roots)%M
                got=norm(stages,L,f,M,alpha);assert got==true,(p,e,f,alpha,got,true)
                C['scalar norm checks']+=1
                alt=next((x for x in range(alpha+1,p) if pe(f,x,p)),None)
                if alt is not None:assert norm(stages,L,f,M,alt)==true;C['chart comparisons']+=1
                if j%4==0:
                    g=[rng.randrange(M) for _ in range(d)]+[1];R2=Ring(M,g);U=tuple(rem([0,1],g,M))
                    qr,_,_=kernel_eval(stages,R2,U,R2.c(1));assert list(qr)==rem(Q,g,M);C['modular kernel checks']+=1
            assert power_sums(stages,12,M,len(roots))==[sum(pow(x,j,M) for x in roots)%M for j in range(13)]
            C['power sum checks']+=13
            outside=[x for x in range(p) if x not in {u%p for u in roots}]
            for j in range(20):
                d=rng.randrange(1,7)
                poles=[rng.choice(outside)+p*rng.randrange(M//p) for _ in range(d)]
                if j%3==0 and d>1:poles[-1]=poles[0];C['repeated-pole controls']+=1
                g=roots_poly(poles,M);P=[rng.randrange(M) for _ in range(d)]
                actual=sum(pe(P,x,M)*pow(pe(g,x,M),-1,M) for x in roots)%M
                assert rational_trace(stages,g,P,M,p)==actual
                C['rational trace checks']+=1
            if e==1:
                ps=points(p,a,c)
                def push(P):
                    for A,B in stages:P=image(P,A,B,p)
                    return P
                for _ in range(60):
                    P,Qp=rng.choice(ps),rng.choice(ps);U,V=push(P),push(Qp)
                    for z in (U,V):
                        if z:assert (z[1]**2-z[0]**3-aa*z[0]-cc)%p==0
                    assert push(add(P,Qp,p,a))==add(U,V,p,aa);C['isogeny homomorphism checks']+=1
            C['towers']+=1
    print(dict(C),flush=True)
    print('All assertions passed. Exact finite component tests; no collective trace benchmark.')
if __name__=='__main__':
    assert __debug__;run()
