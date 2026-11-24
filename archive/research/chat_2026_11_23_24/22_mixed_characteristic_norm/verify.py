#!/usr/bin/env python3
"""Exact finite checks for mixed-characteristic norm readouts.
No fast collective norm evaluator is assumed or implemented; grid products explicit.
"""
from itertools import product
from random import Random
from collections import Counter
from math import prod
from functools import lru_cache
from previous_fiber_verifier import semaev
rng=Random(2026112402)

class Zmod:
    def __init__(self,p,e=1):self.p=p;self.mod=p**e;self.q=p;self.zero=0;self.one=1
    def c(self,x):return x%self.mod
    def add(self,x,y):return (x+y)%self.mod
    def neg(self,x):return -x%self.mod
    def mul(self,x,y):return x*y%self.mod
    def pw(self,x,n):return pow(x,n,self.mod)
    def inv(self,x):return pow(x,-1,self.mod)
    def rand(self):return rng.randrange(self.q)
    def flat(self,x):return (x,)

class Ext:
    def __init__(self,R,coeff):
        self.R=R;self.f=tuple(coeff);self.d=len(coeff);self.q=R.q**self.d
        self.zero=(R.zero,)*self.d;self.one=(R.one,)+(R.zero,)*(self.d-1)
    def c(self,x):return (self.R.c(x),)+(self.R.zero,)*(self.d-1)
    def emb(self,x):return (x,)+(self.R.zero,)*(self.d-1)
    def add(self,x,y):return tuple(self.R.add(a,b) for a,b in zip(x,y))
    def neg(self,x):return tuple(self.R.neg(a) for a in x)
    def mul(self,x,y):
        R=self.R;d=self.d;z=[R.zero]*(2*d-1)
        for i,a in enumerate(x):
            if a==R.zero:continue
            for j,b in enumerate(y):
                if b!=R.zero:z[i+j]=R.add(z[i+j],R.mul(a,b))
        for j in range(2*d-2,d-1,-1):
            a=z[j]
            if a==R.zero:continue
            for i,b in enumerate(self.f):z[j-d+i]=R.add(z[j-d+i],R.neg(R.mul(a,b)))
        return tuple(z[:d])
    def pw(self,x,n):
        r=self.one
        while n:
            if n&1:r=self.mul(r,x)
            n//=2
            if n:x=self.mul(x,x)
        return r
    def inv(self,x):assert x!=self.zero;return self.pw(x,self.q-2) # only used in residue fields
    def rand(self):return tuple(self.R.rand() for _ in range(self.d))
    def flat(self,x):return tuple(a for v in x for a in self.R.flat(v))
    def gen(self):assert self.d>=2;return (self.R.zero,self.R.one)+(self.R.zero,)*(self.d-2)

def trim(a,R):
    a=list(a)
    while a and a[-1]==R.zero:a.pop()
    return a

def sub(a,b,R):
    return trim([R.add(a[i] if i<len(a) else R.zero,R.neg(b[i] if i<len(b) else R.zero)) for i in range(max(len(a),len(b)))],R)

def mul(a,b,R):
    if not a or not b:return []
    c=[R.zero]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=R.add(c[i+j],R.mul(x,y))
    return trim(c,R)

def divrem(a,b,R):
    a=trim(a,R);b=trim(b,R);assert b
    q=[R.zero]*max(0,len(a)-len(b)+1);iv=R.inv(b[-1])
    while len(a)>=len(b):
        k=len(a)-len(b);v=R.mul(a[-1],iv);q[k]=v
        for i,c in enumerate(b):a[i+k]=R.add(a[i+k],R.neg(R.mul(v,c)))
        a=trim(a,R)
    return trim(q,R),a

def pgcd(a,b,R):
    while b:a,b=b,divrem(a,b,R)[1]
    if not a:return []
    iv=R.inv(a[-1]);return [R.mul(x,iv) for x in a]

def ppow(a,k,f,R):
    v=[R.one]
    while k:
        if k&1:v=divrem(mul(v,a,R),f,R)[1]
        k//=2
        if k:a=divrem(mul(a,a,R),f,R)[1]
    return v

def irr(R,d):
    if d==1:return [R.zero]
    for _ in range(1000):
        f=[R.rand() for i in range(d)]+[R.one]
        if f[0]==R.zero:continue
        x=[R.zero,R.one];v=x;good=True
        for i in range(1,d+1):
            v=ppow(v,R.q,f,R)
            if i<=d//2 and len(pgcd(f,sub(v,x,R),R))!=1:good=False;break
        if good and v==x:return f[:-1]
    raise AssertionError('irreducible search cap')

def split(f,R):
    """Cantor-Zassenhaus splitting on promised split squarefree polynomials."""
    n=len(f)-1
    if n==0:return []
    if n==1:return [R.neg(R.mul(f[0],R.inv(f[1])))]
    for _ in range(200):
        a=[R.rand() for i in range(n)]
        g=pgcd(f,sub(ppow(a,(R.q-1)//2,f,R),[R.one],R),R)
        if 1<len(g)<len(f):
            q,r=divrem(f,g,R);assert not r
            return split(g,R)+split(q,R)
    raise AssertionError('split search cap')

@lru_cache(None)
def spec(p,m,K):
    k=Ext(Zmod(p),irr(Zmod(p),m));h=irr(k,K+1)
    # Verify degree K+1 monic polynomial has no root in k via Frobenius/gcd.
    f=h+[k.one];x=[k.zero,k.one]
    assert pgcd(f,sub(ppow(x,k.q,f,k),x,k),k)==[k.one]
    return k,h

def vp(n,p,cap):
    k=0
    while k<cap and n%p==0:n//=p;k+=1
    return k

def mixed_norm(p,m,K,grid,values):
    """One scalar norm over an unramified tower, modulo p^(K+1)."""
    k,h=spec(p,m,K);R=Ext(Zmod(p,K+1),k.f);T=Ext(R,h);th=T.gen();out=T.one
    for s,F in zip(grid,values):
        a=T.add(T.c(F*F),T.mul(T.c(p),T.add(th,T.neg(T.emb(tuple(s))))))
        out=T.mul(out,a)
    H=min(vp(n,p,K+1) for n in T.flat(out))
    if H>K:return 'overflow',H,None
    P=[tuple((a//p**H)%p for a in v) for v in out]
    assert all(v==k.zero for v in P[H+1:]);assert P[H]!=k.zero
    C=P[H];iv=k.inv(C);f=[k.mul(c,iv) for c in P[:H+1]]
    assert f[-1]==k.one
    return set(split(f,k)),H,C

def singleton_norm(p,m,grid,values):
    """One extension, degree m+1, precision p^2."""
    k,_=spec(p,m+1,1);R=Ext(Zmod(p,2),k.f);N=R.one
    for s,F in zip(grid,values):N=R.mul(N,R.add(R.c(F*F),R.mul(R.c(p),tuple(s)+(1,))))
    H=min(vp(n,p,2) for n in N)
    if H==0:return set()
    if H==2:return 'overflow'
    v=[(a//p)%p for a in N];assert v[-1]
    return {tuple(a*pow(v[-1],-1,p)%p for a in v[:-1])}

def tag_norms(p,K,grid,values,tags):
    M=p**(K+2);N0=prod((F*F+p)%M for F in values)%M
    H=vp(N0,p,K+2)
    if H>K:return 'overflow',H
    base=N0//p**H;iv=pow(base,-1,p*p);ans=[]
    for tag in tags:
        N=prod((F*F+p+p*p*tag(s))%M for s,F in zip(grid,values))%M
        assert vp(N,p,K+2)==H
        v=(N//p**H)*iv%(p*p);assert v%p==1
        ans.append((v-1)//p)
    return ans,H

class Tape:
    """Arithmetic-circuit evaluation with reverse-mode differentiation over Z/p^e."""
    def __init__(self,M):self.M=M;self.val=[];self.ops=[]
    def node(self,v,op):self.val.append(v%self.M);self.ops.append(op);return len(self.val)-1
    def const(self,v):return self.node(v,('c',))
    def inp(self,v):return self.node(v,('i',))
    def add(self,i,j):return self.node(self.val[i]+self.val[j],('+',i,j))
    def mul(self,i,j):return self.node(self.val[i]*self.val[j],('*',i,j))
    def reverse(self,o):
        g=[0]*len(self.val);g[o]=1;M=self.M
        for k in range(len(g)-1,-1,-1):
            op=self.ops[k]
            if op[0] in ('c','i'):continue
            _,i,j=op
            if op[0]=='+':g[i]=(g[i]+g[k])%M;g[j]=(g[j]+g[k])%M
            else:g[i]=(g[i]+g[k]*self.val[j])%M;g[j]=(g[j]+g[k]*self.val[i])%M
        return g

def reverse_tags(p,K,grid,values,tags):
    T=Tape(p**(K+1));xs=[T.inp(0) for _ in tags];factors=[]
    for s,F in zip(grid,values):
        f=T.const(F*F+p)
        for x,tag in zip(xs,tags):f=T.add(f,T.mul(x,T.const(p*tag(s))))
        factors.append(f)
    # Balanced scalar norm circuit; the same reverse rule applies to a fast norm circuit.
    while len(factors)>1:
        factors=[T.mul(factors[i],factors[i+1]) if i+1<len(factors) else factors[i] for i in range(0,len(factors),2)]
    o=factors[0] if factors else T.const(1);H=vp(T.val[o],p,K+1)
    if H>K:return 'overflow',H,len(T.val)
    g=T.reverse(o);C=T.val[o]//p**H%p;iv=pow(C,-1,p)
    assert all(g[x]%p**H==0 for x in xs)
    return [g[x]//p**H*iv%p for x in xs],H,len(T.val)

def trial(p,m,K,grid,vals,cnt):
    truth={s for s,v in zip(grid,vals) if v%p==0};H=len(truth)
    a,h,c=mixed_norm(p,m,K,grid,vals);cnt['tower_norms']+=1
    if H>K:assert a=='overflow' and h==K+1;cnt['overflow']+=1
    else:
        assert a==truth and h==H
        C=prod((v*v)%p for v in vals if v%p)%p
        assert c==(C,)+(0,)*(m-1)
        cnt['recovered_addresses']+=H;cnt['empty']+=H==0;cnt['multiple']+=H>1
        vals2=[v+p*rng.randint(-9,9) for v in vals]
        assert mixed_norm(p,m,K,grid,vals2)==(a,h,c);cnt['lift_invariance']+=1
    z=singleton_norm(p,m,grid,vals);assert z==(truth if H<=1 else 'overflow');cnt['singleton_norms']+=1
    tags=[lambda s:1]+[lambda s,j=j:s[j] for j in range(m)]+[lambda s:sum(s)**2,lambda s:prod(s)]
    tags.extend([lambda s,j=j:(sum((i+1)*v for i,v in enumerate(s)))**j for j in (3,4)])
    r,h=tag_norms(p,K,grid,vals,tags);cnt['tag_batches']+=1
    if H>K:assert r=='overflow' and h>=K+1
    else:
        assert h==H and r==[sum(t(s) for s in truth)%p for t in tags];cnt['tag_comparisons']+=len(tags)
    r,h,nodes=reverse_tags(p,K,grid,vals,tags);cnt['reverse_batches']+=1;cnt['tape_nodes']+=nodes
    if H>K:assert r=='overflow' and h==K+1
    else:assert h==H and r==[sum(t(s) for s in truth)%p for t in tags];cnt['reverse_comparisons']+=len(tags)
    cnt['grid_addresses']+=len(grid)

if __name__=='__main__':
    C=Counter()
    for q in range(120):
        p=rng.choice((5,7,11,13));m=rng.choice((2,3));K=rng.randrange(1,4)
        grid=list(product(*[rng.sample(range(p),rng.randrange(1,4)) for i in range(m)]))
        H=rng.randrange(min(len(grid),K+2)+1);roots=set(rng.sample(grid,H));vals=[(0 if s in roots else rng.randrange(1,p))+p*rng.randint(-7,7) for s in grid]
        trial(p,m,K,grid,vals,C)
    print('Arbitrary explicit grids:',dict(C),flush=True)
    C=Counter()
    for p,a,c,T in ((13,4,1,0),(17,2,2,0),(19,1,0,0),(37,1,1,0)):
        assert (4*a**3+27*c*c)%p
        F=semaev(p,a,c)
        for q in range(18):
            m=rng.choice((2,3));K=rng.randrange(1,4)
            pool=rng.sample(range(p),2*m);grid=list(product(*[pool[2*i:2*i+2] for i in range(m)]))
            target=rng.randrange(p);vals=[F(s,target)+p*rng.randint(-4,4) for s in grid]
            trial(p,m,K,grid,vals,C)
    print('Actual Semaev grids (unprojected):',dict(C),flush=True)
    C=Counter()
    for p,m in ((5,2),(7,3)):
        grid=list(product(range(2),repeat=m))
        for vals in ([0]*len(grid),[p]*len(grid),[1]*len(grid),[0]+[1]*(len(grid)-1)):
            trial(p,m,2,grid,vals,C)
    print('Whole-grid zero and exact overflow controls:',dict(C),flush=True)
    print('PASS. Explicit finite norm products only; no fast collective evaluator implemented.')
