#!/usr/bin/env python3
"""Exact toy-grid checks. Uses explicit products and schoolbook ring arithmetic."""
from itertools import product
from collections import Counter,defaultdict
from functools import lru_cache
from math import comb, prod,gcd
from random import Random
from previous_mixed_verifier import Zmod,Ext,irr,split,vp
rng=Random(2026112405)
def sub(R,a,b):return R.add(a,R.neg(b))
def muls(R,*xs):
    v=R.one
    for x in xs:v=R.mul(v,x)
    return v
def pm(R,f,g):
    o=[R.zero]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        if a!=R.zero:
            for j,b in enumerate(g):
                if b!=R.zero:o[i+j]=R.add(o[i+j],R.mul(a,b))
    return o
def pe(R,f,x):
    a=R.zero
    for v in reversed(f):a=R.add(R.mul(a,x),v)
    return a
def det(R,M):
    # Small division-free determinant; exponential DP is ONLY for tests.
    D={0:R.one};n=len(M)
    for row in M:
        E={}
        for mask,a in D.items():
            if a==R.zero:continue
            for j,b in enumerate(row):
                if mask>>j&1 or b==R.zero:continue
                v=R.mul(a,b)
                if (mask>>(j+1)).bit_count()&1:v=R.neg(v)
                k=mask|1<<j;E[k]=R.add(E.get(k,R.zero),v)
        D=E
    return D.get((1<<n)-1,R.zero)
def res(R,f,g):
    d,e=len(f)-1,len(g)-1;n=d+e;M=[[R.zero]*n for _ in range(n)]
    for i in range(e):M[i][i:i+d+1]=f[::-1]
    for i in range(d):M[e+i][i:i+e+1]=g[::-1]
    return det(R,M)
def interp(R,vals,p,e):
    assert len(vals)<=p;out=[R.zero]*len(vals);D=vals[:];B=[R.one];fac=1
    for j in range(len(vals)):
        v=R.mul(D[0],R.c(pow(fac,-1,p**e)))
        for i,b in enumerate(B):out[i]=R.add(out[i],R.mul(v,b))
        D=[sub(R,b,a) for a,b in zip(D,D[1:])];B=pm(R,B,[R.c(-j),R.one]);fac=fac*(j+1)%p**e
    return out
def s3(R,u,v,a,c):
    uv=R.mul(u,v);s=R.add(u,v);d=sub(R,uv,R.c(a));e=sub(R,u,v)
    return [sub(R,R.mul(d,d),R.mul(R.c(4*c),s)),R.mul(R.c(-2),R.add(R.mul(s,R.add(uv,R.c(a))),R.c(2*c))),R.mul(e,e)]
def poly(R,xs,p,e,a,c):
    F=[R.neg(xs[0]),R.one]
    for x in xs[1:]:
        d=2*(len(F)-1);F=interp(R,[res(R,F,s3(R,x,R.c(j),a,c)) for j in range(d+1)],p,e)
    return F
def proj(R,X,a,c,h):
    if h==1:return [R.neg(X),R.one]
    assert h==2
    return [sub(R,R.c(a*a),R.mul(R.c(4*c),X)),sub(R,R.c(-8*c),R.mul(R.c(4*a),X)),R.c(-2*a),R.mul(R.c(-4),X),R.one]
def value(R,xs,X,p,e,a,c,h):
    if h==1:
        if len(xs)==1:return sub(R,X,xs[0])
        if len(xs)==2:return pe(R,s3(R,*xs,a,c),X)
        return res(R,poly(R,xs[:-1],p,e,a,c),s3(R,xs[-1],X,a,c))
    return res(R,poly(R,xs,p,e,a,c),proj(R,X,a,c,h))
def grid(R,xs,X,p,e,a,c,h):
    v=R.one
    for z in product(*xs):v=R.mul(v,value(R,z,X,p,e,a,c,h))
    return v
@lru_cache(None)
def nr(p):return next(i for i in range(2,p) if pow(i,(p-1)//2,p)==p-1)
class Curve:
    def __init__(self,p,a,c):
        assert (4*a**3+27*c*c)%p
        self.p=p;self.a=a;self.c=c;self.F=Ext(Zmod(p),[(-nr(p))%p,0]);R=self.F
        D=defaultdict(list)
        for y in product(range(p),repeat=2):D[R.mul(y,y)].append(y)
        self.roots={x:tuple(y for y in D[((x**3+a*x+c)%p,0)] if y!=R.zero) for x in range(p)}
        self.pts=[((x,0),y) for x,L in self.roots.items() for y in L if y[1]==0]
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        R=self.F;x,y=P;u,v=Q
        if x==u:
            if R.add(y,v)==R.zero:return None
            assert y==v
            z=R.mul(R.add(R.mul(R.c(3),R.mul(x,x)),R.c(self.a)),R.inv(R.mul(R.c(2),y)))
        else:z=R.mul(sub(R,v,y),R.inv(sub(R,u,x)))
        w=sub(R,sub(R,R.mul(z,z),x),u);t=sub(R,R.mul(z,sub(R,x,w)),y)
        assert R.mul(t,t)==R.add(R.add(R.mul(w,R.mul(w,w)),R.mul(R.c(self.a),w)),R.c(self.c))
        return w,t
    def truth(self,S,A,t,h):
        D=defaultdict(set);C=Counter();p=self.p
        for s in product(*S):
            xs=[(A[i]+t[i]*v)%p for i,v in enumerate(s)]
            assert all(self.roots[x] for x in xs)
            for ys in product(*(self.roots[x] for x in xs)):
                Q=None
                for x,y in zip(xs,ys):Q=self.add(Q,((x,0),y))
                C['sum_at_infinity']+=Q is None
                if h==2:Q=self.add(Q,Q)
                C['projected_infinity']+=Q is None;C['sign_choices']+=1
                if Q is not None and Q[0][1]==Q[1][1]==0 and Q[1]!=(0,0):D[Q].add((s,ys))
        return D,C
class Jet:
    def __init__(self,p,n):self.p=p;self.n=n;self.zero=(0,)*(n+1);self.one=self.c(1)
    def c(self,x):return (x%self.p,)+(0,)*self.n
    def var(self,x,j):return (x%self.p,)+tuple(int(i==j) for i in range(self.n))
    def add(self,x,y):return tuple((a+b)%self.p for a,b in zip(x,y))
    def neg(self,x):return tuple(-a%self.p for a in x)
    def mul(self,x,y):return (x[0]*y[0]%self.p,)+tuple((x[0]*b+y[0]*a)%self.p for a,b in zip(x[1:],y[1:]))
def derivative(C,S,A,t,h,T,truth):
    p=C.p;m=len(S);R=Jet(p,2*m+1);X=R.var(T[0][0],2*m)
    xs=[[R.add(R.var(A[i],2*i),R.mul(R.var(t[i],2*i+1),R.c(s))) for s in L] for i,L in enumerate(S)]
    g=grid(R,xs,X,p,1,C.a,C.c,h);H=len(truth)
    assert (g[0]==0)==bool(H)
    if H>1:assert g==R.zero
    if H==1:
        s,ys=next(iter(truth));assert g[-1]!=0
        for i in range(m):
            da,dt=g[1+2*i:3+2*i];assert da and ys[i][1]==0
            assert dt*pow(da,-1,p)%p==s[i]
            assert (-h*T[1][0]*g[-1]*pow(da,-1,p))%p==ys[i][0]
    return H
@lru_cache(None)
def irr0(p,d):return tuple(irr(Zmod(p),d))
@lru_cache(None)
def tower(p,m):
    # cap=2, m=2 or 3: coprime extension degrees.
    assert gcd(2,2*m+1)==gcd(3,2*(2*m+1))==1
    k0=Ext(Zmod(p),[(-nr(p))%p,0]);R0=Ext(Zmod(p,3),k0.f)
    f=[k0.c(c) for c in irr0(p,2*m+1)];k=Ext(k0,f);R=Ext(R0,f)
    L=Ext(R,[k.c(c) for c in irr0(p,3)])
    return k0,k,R,L
def packed(C,S,A,t,h,T,truth):
    p=C.p;m=len(S);k0,k,R,L=tower(p,m);v=L.emb(R.gen());theta=L.gen();V=[L.one]
    for i in range(2*m):V.append(L.mul(V[-1],v))
    xs=[[L.add(L.c(A[i]+t[i]*s),L.mul(L.c(p),L.add(V[2*i+1],L.mul(L.c(s),V[2*i+2])))) for s in L0] for i,L0 in enumerate(S)]
    X=L.add(L.c(T[0][0]),L.mul(L.c(p),theta));g=grid(L,xs,X,p,3,C.a,C.c,h)
    H=min(vp(a,p,3) for a in L.flat(g));N=len(truth);assert H==min(N,3),(p,m,h,H,N)
    if H>2:return 'overflow',N
    P=[tuple(tuple(a//p**H%p for a in v) for v in u) for u in g]
    assert all(a==k.zero for a in P[H+1:]) and P[H]!=k.zero
    assert all(a==k0.zero for a in P[H][1:]);iv=k.inv(P[H]);P=[k.mul(a,iv) for a in P[:H+1]]
    got=set();hy=k0.mul(k0.c(h),T[1])
    for r in split(P,k):
        ss=[];ys=[]
        for i in range(m):
            u,v=r[2*i+1:2*i+3];assert u!=k0.zero
            j=k0.mul(v,k0.inv(u));assert j[1]==0
            ss.append(j[0]);ys.append(k0.mul(hy,k0.inv(u)))
        got.add((tuple(ss),tuple(ys)))
    assert got==truth,(got,truth)
    return 'ok',N

def packed_single(C,S,A,t,h,T,truth):
    p=C.p;m=len(S);k=Ext(Zmod(p),irr0(p,2*m+2));R=Ext(Zmod(p,2),k.f)
    v=R.gen();V=[R.one]
    for _ in range(2*m+1):V.append(R.mul(V[-1],v))
    xs=[[R.add(R.c(A[i]+t[i]*s),R.mul(R.c(p),R.add(V[2*i+1],R.mul(R.c(s),V[2*i+2])))) for s in L] for i,L in enumerate(S)]
    X=R.add(R.c(T[0][0]),R.mul(R.c(p),V[-1]))
    z=grid(R,xs,X,p,2,C.a,C.c,h);H=min(vp(a,p,2) for a in z);assert H==min(len(truth),2)
    if H==1:
        P=[a//p%p for a in z];iv=pow(P[-1],-1,p);P=[a*iv%p for a in P]
        s=tuple(P[2*i+2]*pow(P[2*i+1],-1,p)%p for i in range(m))
        ys=tuple((-h*T[1][0]*pow(P[2*i+1],-1,p)%p,0) for i in range(m))
        assert {(s,ys)}==truth
    return H

def chart(R,f,a):
    # Homogeneous change (U,V)=(alpha*z+1,z), det=-1.
    d=len(f)-1;o=[R.zero]*(d+1);V=[R.one]
    for i in range(d):V.append(R.mul(V[-1],a))
    for i,c in enumerate(f):
        for j in range(i+1):o[d-i+j]=R.add(o[d-i+j],muls(R,c,R.c(comb(i,j)),V[j]))
    return o
def rem(R,f,g):
    d=len(g)-1;assert d and g[-1]==R.one;f=f[:]
    for i in range(len(f)-1,d-1,-1):
        a=f[i]
        for j in range(d):f[i-d+j]=sub(R,f[i-d+j],R.mul(a,g[j]))
    return (f[:d]+[R.zero]*d)[:d]
def tree(R,gs):
    if len(gs)==1:return gs[0],None,None
    j=len(gs)//2;L=tree(R,gs[:j]);Q=tree(R,gs[j:]);return pm(R,L[0],Q[0]),L,Q
def walk(R,f,T):
    f=rem(R,f,T[0]);L,Q=T[1:]
    return [f] if L is None else walk(R,f,L)+walk(R,f,Q)
def unit(R,a,p):return any(v%p for v in R.flat(a))
def half(C,S,A,t,h,T):
    p=C.p;e=2;R=Ext(Zmod(p,e),[(-nr(p))%p,0]);m=len(S)
    xs=[[R.add(R.c(A[i]+t[i]*s),(0,p*rng.randrange(p))) for s in L] for i,L in enumerate(S)]
    X=R.add(R.c(T[0][0]),(0,p*rng.randrange(p)))
    # Balanced for h=1; one-slot left half for the small h=2 tests.
    mid=1 if h==2 else m//2
    Qs=[poly(R,z,p,e,C.a,C.c) for z in product(*xs[mid:])]
    if h==1:
        Bs=[]
        for z in product(*xs[:mid]):
            F=poly(R,z,p,e,C.a,C.c);d0=2*(len(F)-1)
            Bs.append(interp(R,[res(R,F,s3(R,R.c(j),X,C.a,C.c)) for j in range(d0+1)],p,e))
    else:Bs=[interp(R,[res(R,s3(R,z,R.c(j),C.a,C.c),proj(R,X,C.a,C.c,h)) for j in range(9)],p,e) for z in xs[0]]
    d=sum(len(f)-1 for f in Qs);assert sum(len(f)-1 for f in Bs)<p*p
    alpha=next(z for z in product(range(p),repeat=2) if all(unit(R,pe(R,B,z),p) for B in Bs))
    Q=[R.one]
    for f in Qs:Q=pm(R,Q,chart(R,f,alpha))
    gs=[];lc=[]
    for B in Bs:
        g=chart(R,B,alpha);l=g[-1];iv=R.pw(l,(R.q-1)*p**(e-1)-1);assert R.mul(l,iv)==R.one
        gs.append([R.mul(x,iv) for x in g]);lc.append(l)
    fast=R.one;literal=R.one;rs=walk(R,Q,tree(R,gs))
    for g,r,l in zip(gs,rs,lc):
        assert r==rem(R,Q,g);fast=muls(R,fast,R.pw(l,d),res(R,g,r))
    for f in Qs:
        for B in Bs:literal=R.mul(literal,res(R,f,B))
    assert fast==literal,(fast,literal)
    return len(Bs),len(Qs)*len(Bs),any(not unit(R,B[-1],p) for B in Bs)

def fixtures():
    curves=[Curve(13,4,1),Curve(17,2,2),Curve(19,1,0),Curve(23,1,1),Curve(29,1,0),Curve(37,1,1)]
    quota=Counter();want={0:8,1:12,2:8,3:5};items=[];mixed=0
    for _ in range(1000):
        C=rng.choice(curves);p=C.p;m=rng.choice((2,3));h=rng.choice((1,1,2));w=rng.choice((1,2,2,3))
        av=[x for x,L in C.roots.items() if L]
        if len(av)<m*w:continue
        pool=rng.sample(av,m*w);A=[rng.randrange(p) for i in range(m)];t=[rng.randrange(1,p) for i in range(m)]
        S=[tuple((x-A[i])*pow(t[i],-1,p)%p for x in pool[i*w:(i+1)*w]) for i in range(m)]
        D,stat=C.truth(S,A,t,h);T=rng.choice(C.pts);H=len(D.get(T,set()));key=min(H,3)
        if quota[key]<want[key]:items.append((C,S,A,t,h,T,D.get(T,set()),stat));quota[key]+=1
        if mixed<2:
            opts=[(T,Z) for T,Z in D.items() if 0<len(Z)<=2 and any(any(y[1] for y in ys) for _,ys in Z)]
            if opts:
                T,Z=rng.choice(opts);items.append((C,S,A,t,h,T,Z,stat));mixed+=1
        if all(quota[k]>=v for k,v in want.items()) and mixed>=2:break
    assert all(quota[k]>=v for k,v in want.items()),quota
    print('Selected fixture quotas:',dict(quota),'extra nonrational controls:',mixed,flush=True)
    return items
if __name__=='__main__':
    items=fixtures();C=Counter()
    for i,(E,S,A,t,h,T,truth,stat) in enumerate(items):
        N=derivative(E,S,A,t,h,T,truth);C.update(stat);C['jet_grids']+=1;C['singleton_jet_coordinates']+=2*len(S)*(N==1);C['projected_grids']+=h==2
        C['grid_addresses']+=prod(map(len,S))
        H1=packed_single(E,S,A,t,h,T,truth);C['single_extension_grids']+=1;C['single_extension_recovered']+=H1==1
        # Packed tests on all m=2 and selected m=3, to limit slow schoolbook arithmetic.
        if len(S)==2 or (max(map(len,S))<=2 and h==1):
            kind,N=packed(E,S,A,t,h,T,truth);C['packed_grids']+=1
            if kind=='overflow':C['overflow']+=1
            else:
                C['empty']+=N==0;C['singleton']+=N==1;C['multiple']+=N>1;C['recovered_signed_addresses']+=N
                C['recovered_nonrational']+=sum(any(y[1] for y in ys) for s,ys in truth)
        if max(map(len,S))<=2 and (h==1 or len(S)==2):
            a,b,d=half(E,S,A,t,h,T);C['half_grids']+=1;C['remainder_comparisons']+=a;C['pairwise_resultants']+=b;C['degree_drop_grids']+=d
        print('Completed',i+1,'/',len(items),flush=True)
    for _ in range(8):
        E=Curve(17,2,2);p=E.p;m=4
        av=[x for x,L in E.roots.items() if L];pool=rng.sample(av,8)
        S=[tuple(pool[2*i:2*i+2]) for i in range(m)];A=[0]*m;t=[1]*m;T=rng.choice(E.pts)
        a,b,d=half(E,S,A,t,1,T);C['balanced_four_slot_grids']+=1;C['remainder_comparisons']+=a;C['pairwise_resultants']+=b;C['degree_drop_grids']+=d
    print('Final counts:',dict(C),flush=True)
    print('PASS. Explicit finite grid products only; no faster ECDLP exponent demonstrated.',flush=True)
