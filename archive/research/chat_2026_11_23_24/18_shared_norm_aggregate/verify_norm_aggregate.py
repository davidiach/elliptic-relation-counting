#!/usr/bin/env python3
"""Exact toy checks of primitive shared-translation relation polynomials.
Uses previous checkpoint's squarefree multiquadratic arithmetic on small p.
Checks include nonrational lifts via F_{p^2}, torsion and exceptional sums.
No asymptotic runtime claim is measured by these fixtures.
"""
from itertools import product, combinations
from collections import Counter
from random import Random
from inverse_translation_verify import Ring, symbolic_sum
from verify import good_offsets


def divrem(a,b,p):
    a=list(a); b=list(b); out=[0]*max(0,len(a)-len(b)+1)
    while a and len(a)>=len(b):
        k=len(a)-len(b); c=a[-1]*pow(b[-1],-1,p)%p; out[k]=c
        for j,v in enumerate(b): a[j+k]=(a[j+k]-c*v)%p
        while a and not a[-1]: a.pop()
    while out and not out[-1]:out.pop()
    return tuple(out),tuple(a)

def gcd(a,b,p):
    while b:a,b=b,divrem(a,b,p)[1]
    return tuple(v*pow(a[-1],-1,p)%p for v in a) if a else ()

def conj(R,a,j):
    return {mask:R.trim([-v for v in c]) if (mask>>j)&1 else c for mask,c in a.items()}

def addJ(R,P,Q):
    X,Y,Z=P; U,V,W=Q
    z2=R.sq(Z); w2=R.sq(W)
    x1=R.mul(X,w2);x2=R.mul(U,z2)
    y1=R.mul(Y,R.mul(W,w2));y2=R.mul(V,R.mul(Z,z2))
    H=R.sub(x2,x1); D=R.sub(y2,y1);H2=R.sq(H); H3=R.mul(H,H2)
    B=R.mul(x1,H2); X3=R.sub(R.sub(R.sq(D),H3),R.scale(B,2))
    Y3=R.sub(R.mul(D,R.sub(B,X3)),R.mul(y1,H3)); Z3=R.mul(R.mul(Z,W),H)
    assert Z3
    return X3,Y3,Z3

def doubleJ(R,P,a):
    X,Y,Z=P; y2=R.sq(Y);S=R.scale(R.mul(X,y2),4)
    M=R.add(R.scale(R.sq(X),3),R.scale(R.sq(R.sq(Z)),a))
    XX=R.sub(R.sq(M),R.scale(S,2)); YY=R.sub(R.mul(M,R.sub(S,XX)),R.scale(R.sq(y2),8));ZZ=R.scale(R.mul(Y,Z),2)
    assert ZZ
    return XX,YY,ZZ

def mulJ(R,P,a,h):
    if h==1:return P
    Q=mulJ(R,P,a,h//2);Q=doubleJ(R,Q,a)
    return addJ(R,Q,P) if h%2 else Q

def relation_poly(p,a,c,ss,h=1):
    assert good_offsets(p,a,c,ss)
    R=Ring(p,a,c,ss); P=mulJ(R,symbolic_sum(R,ss),a,h)
    X,Y,Z=P
    coeff=[R.scale(X,-1),R.sq(Z)] # polynomial in target coordinate z
    # These coefficients are invariant under simultaneous negation.
    assert all(all(mask.bit_count()%2==0 for mask in x) for x in coeff)
    for j in range(1,len(ss)):
        cc=[conj(R,x,j) for x in coeff]; out=[{} for _ in range(2*len(coeff)-1)]
        for i,x in enumerate(coeff):
            for k,y in enumerate(cc):out[i+k]=R.add(out[i+k],R.mul(x,y))
        coeff=out
        assert all(all(mask>>j&1==0 for mask in x) for x in coeff)
    assert all(set(x)<= {0} for x in coeff)
    coeff=[x.get(0,()) for x in coeff];content=()
    for x in coeff:content=gcd(content,x,p)
    primitive=[]
    for x in coeff:
        q,rem=divrem(x,content,p);assert not rem;primitive.append(q)
    # Normalize by first nonzero coefficient of highest target power.
    u=pow(primitive[-1][-1],-1,p)
    out=tuple(tuple(v*u%p for v in x) for x in primitive)
    assert len(out)-1==2**(len(ss)-1)
    assert max(len(x)-1 for x in out)==h*h*len(ss)*2**(len(ss)-1)
    return out,max(len(x)-1 for x in coeff),len(content)-1

def ev(poly,x,p):
    s=0
    for a in reversed(poly):s=(s*x+a)%p
    return s

def eval2(F,A,z,p):return ev([ev(x,A,p) for x in F],z,p)

class Curve2:
    def __init__(self,p,a,c):
        self.p=p;self.a=a;self.c=c;self.nr=next(i for i in range(2,p) if pow(i,(p-1)//2,p)==p-1)
        self.zero=(0,0);self.one=(1,0); self.sqroots={}
        for i in range(p):
            self.sqroots.setdefault(self.mul((i,0),(i,0))[0],set()).add((i,0))
            self.sqroots.setdefault(self.mul((0,i),(0,i))[0],set()).add((0,i))
    def add(self,x,y):return ((x[0]+y[0])%self.p,(x[1]+y[1])%self.p)
    def neg(self,x):return (-x[0]%self.p,-x[1]%self.p)
    def sub(self,x,y):return self.add(x,self.neg(y))
    def mul(self,x,y):return ((x[0]*y[0]+self.nr*x[1]*y[1])%self.p,(x[0]*y[1]+x[1]*y[0])%self.p)
    def inv(self,x):
        d=pow((x[0]**2-self.nr*x[1]**2)%self.p,-1,self.p);return (x[0]*d%self.p,-x[1]*d%self.p)
    def scale(self,x,n):return (x[0]*n%self.p,x[1]*n%self.p)
    def padd(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and self.add(y,v)==self.zero:return None
        if x==u:k=self.mul(self.add(self.scale(self.mul(x,x),3),(self.a,0)),self.inv(self.scale(y,2)))
        else:k=self.mul(self.sub(v,y),self.inv(self.sub(u,x)))
        z=self.sub(self.sub(self.mul(k,k),x),u)
        return z,self.sub(self.mul(k,self.sub(x,z)),y)
    def pmul(self,h,P):
        Q=None
        while h:
            if h&1:Q=self.padd(Q,P)
            P=self.padd(P,P);h>>=1
        return Q
    def points_at(self,x):return [((x,0),y) for y in self.sqroots[(x**3+self.a*x+self.c)%self.p]]
    def sums(self,A,ss,h=1):
        for pts in product(*(self.points_at((A+s)%self.p) for s in ss)):
            S=None
            for P in pts:S=self.padd(S,P)
            yield pts,self.pmul(h,S)

def fixture(p,a,c,ss,h):
    F,olddeg,content=relation_poly(p,a,c,ss,h)
    E=Curve2(p,a,c); checks=alg=rat=bad=branch=inf=0
    for A in range(p):
        truths=set();real=set()
        for pts,T in E.sums(A,ss,h):
            if T is None:inf+=1;continue
            if T[0][1]:continue
            truths.add(T[0][0]);alg+=1
            if all(y[1]==0 and y[0]!=0 for x,y in pts) and T[1][1]==0:
                real.add(T[0][0]);rat+=1
            if any(y==(0,0) for x,y in pts):branch+=1
        for z in range(p):
            assert (eval2(F,A,z,p)==0)==(z in truths),(p,ss,h,A,z,F,truths)
            checks+=1
        bad+=len(truths-real)
    print(f'p={p}, m={len(ss)}, h={h}: degrees (A,X)=({max(len(x)-1 for x in F)},{len(F)-1}); raw degree={olddeg}, removed content={content}; {checks} exact incidence checks; {alg} algebraic signed finite hits, {rat} rational nonbranch hits, {bad} extra coordinate hits, {branch} branch-path hits, {inf} infinity paths',flush=True)
    return checks,alg,rat,bad,branch,inf


def orbits_check():
    # For disjoint interval slots, canonical representative has min digit zero.
    tests=0;covered=0
    for m in range(2,6):
        for w in range(1,8):
            reps=[d for d in product(range(w),repeat=m) if min(d)==0]
            assert len(reps)==w**m-(w-1)**m
            grid=Counter(tuple(x+j for x in d) for d in reps for j in range(w-max(d)))
            assert len(grid)==w**m and all(v==1 for v in grid.values())
            covered+=len(grid);tests+=1
    print(f'{tests} exact interval orbit partitions, {covered} covered address tuples',flush=True)
    # Verify bivariate shift equivalence, distinct normalized polynomials and
    # equality of batched root lookup vs literal product on a complete grid.
    p,a,c,m,w=19,2,1,3,2
    ss=list(range(m*w))
    assert good_offsets(p,a,c,ss)
    E=Curve2(p,a,c)
    def offsets(d):return [i*w+v for i,v in enumerate(d)]
    reps=[d for d in product(range(w),repeat=m) if min(d)==0]
    polys={d:relation_poly(p,a,c,offsets(d))[0] for d in reps}
    allpoly={d:relation_poly(p,a,c,offsets(d))[0] for d in product(range(w),repeat=m)}
    assert len(set(allpoly.values()))==w**m
    roots={(d,z):{A for A in range(p) if eval2(F,A,z,p)==0} for d,F in polys.items() for z in range(p)}
    batched=0;shiftchecks=0
    for d,F in allpoly.items():
        j=min(d);rep=tuple(v-j for v in d)
        # Normalization may change a scalar, but the zero sets must coincide.
        for A,z in product(range(p),repeat=2):
            assert (eval2(F,A,z,p)==0)==(eval2(polys[rep],(A+j)%p,z,p)==0)
            shiftchecks+=1
    for A,z in product(range(p),repeat=2):
        expected={d for d,F in allpoly.items() if eval2(F,A,z,p)==0}
        got={tuple(v+j for v in d) for d in reps for j in range(w-max(d)) if (A+j)%p in roots[d,z]}
        assert expected==got
        batched+=1
    print(f'{len(allpoly)} distinct normalized polynomials, {shiftchecks} shift-incidence comparisons, {batched} complete batched forward queries',flush=True)

if __name__=='__main__':
    totals=[0]*6
    for args in [(13,4,1,[0,1],1),(17,2,2,[0,1,2],1),(37,0,5,[0,1,3,4],1),(43,2,1,[0,1,2,3],1),(19,2,1,[0,1],2),(19,2,1,[0,1,2],2),(19,2,1,[0,1],3),(29,1,0,[0,1],4),(43,2,1,[0,1],8)]:
        q=fixture(*args);totals=[x+y for x,y in zip(totals,q)]
    print('totals (incidence, algebraic signed hits, rational nonbranch signed hits, extra coordinates, branch paths, infinity paths):',totals,flush=True)
    orbits_check()
    print('PASS: exact primitive image polynomials and diagonal batching.',flush=True)
