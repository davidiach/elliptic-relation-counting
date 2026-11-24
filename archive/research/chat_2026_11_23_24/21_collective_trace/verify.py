#!/usr/bin/env python3
"""Exact component regressions for projected norms and rational traces.
Explicit finite product-ring lanes; NOT an implementation of a fast collective trace.
The previous checkpoint's Semaev/interpolation helpers are reused with attribution.
"""
from collections import Counter
from itertools import product
from math import prod, isqrt
from random import Random
from previous_fiber_verifier import det, interp, pe, roots_poly, semaev, solve

rng = Random(241126)

class Fp2:
    """a+b*i encoded as a+p*b, i^2=nonsquare."""
    def __init__(self,p):
        self.p=p
        self.ns=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
    def add(self,x,y):
        p=self.p;return (x%p+y%p)%p+p*((x//p+y//p)%p)
    def neg(self,x):
        p=self.p;return (-x%p)+p*(-(x//p)%p)
    def sub(self,x,y):return self.add(x,self.neg(y))
    def mul(self,x,y):
        p=self.p;a,b=x%p,x//p;c,d=y%p,y//p
        return (a*c+self.ns*b*d)%p+p*((a*d+b*c)%p)
    def pw(self,x,k):
        r=1
        while k:
            if k&1:r=self.mul(r,x)
            x=self.mul(x,x);k//=2
        return r
    def inv(self,x):
        assert x
        return self.pw(x,self.p*self.p-2)
    def sm(self,xs):
        r=0
        for x in xs:r=self.add(r,x)
        return r
    def pe(self,f,x):
        r=0
        for a in reversed(f):r=self.add(self.mul(r,x),a)
        return r

class ProductDual:
    """Product of N square-zero extensions of Fp2, K formal tags."""
    def __init__(self,p,n,k):
        self.F=Fp2(p);self.p=p;self.n=n;self.k=k
        self.z=self.c(0);self.o=self.c(1)
    def c(self,a):return ((a%self.p,)+(0,)*self.k,)*self.n
    def add(self,a,b):
        return tuple(tuple(self.F.add(x,y) for x,y in zip(u,v)) for u,v in zip(a,b))
    def neg(self,a):return tuple(tuple(self.F.neg(x) for x in u) for u in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):
        F=self.F
        return tuple((F.mul(u[0],v[0]),)+tuple(F.add(F.mul(u[0],y),F.mul(v[0],x)) for x,y in zip(u[1:],v[1:])) for u,v in zip(a,b))
    def pm(self,a,b):
        r=[self.z]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):r[i+j]=self.add(r[i+j],self.mul(x,y))
        return r
    def shift(self,f,j):
        r=[f[-1]]
        for a in f[-2::-1]:
            r=self.pm(r,[self.c(j),self.o]);r[0]=self.add(r[0],a)
        return r
    def rem(self,f,g):
        assert g[-1]==self.o
        d=len(g)-1
        if not d:return []
        r=list(f)+[self.z]*max(0,d-len(f))
        for j in range(len(r)-1,d-1,-1):
            for i in range(d):r[j-d+i]=self.sub(r[j-d+i],self.mul(r[j],g[i]))
        return r[:d]
    def ev(self,f,x):
        r=self.z
        for a in reversed(f):r=self.add(self.mul(r,self.c(x)),a)
        return r
    def balanced_product(self,fs):
        if not fs:return [self.o]
        while len(fs)>1:fs=[self.pm(fs[j],fs[j+1]) if j+1<len(fs) else fs[j] for j in range(0,len(fs),2)]
        return fs[0]
    def tree(self,xs):
        if len(xs)==1:return ([self.c(-xs[0]),self.o],None,None)
        mid=len(xs)//2;l=self.tree(xs[:mid]);r=self.tree(xs[mid:]);return (self.pm(l[0],r[0]),l,r)
    def multipoint(self,f,xs):
        def walk(g,t):
            r=self.rem(g,t[0]);a,b=t[1:]
            return [r[0]] if a is None else walk(r,a)+walk(r,b)
        return walk(f,self.tree(xs))
    def block_norm(self,f,w,H=None):
        assert w>=1
        d=max(1,len(f)-1)
        if H is None:H=max(1,min(w,isqrt(max(1,w//d))))
        J=w//H
        B=self.balanced_product([self.shift(f,j) for j in range(H)])
        acc=self.o
        for v in self.multipoint(B,[q*H for q in range(J)]):acc=self.mul(acc,v)
        for j in range(J*H,w):acc=self.mul(acc,self.ev(f,j))
        return acc


def residue_from_distinct_poles(F,vals,weights):
    """Construct common-denominator numerator at t=0 independently.
    Residue N(0)/Q'(0), if 0 is a pole of Q; cancellation allowed."""
    S=set(vals)
    if 0 not in S:return 0
    p=F.p;prod_nonzero=prod(v for v in S if v)%p
    assert prod_nonzero
    N=0
    for v,a in zip(vals,weights):
        q_over_factor=prod(u for u in S if u!=v)%p
        N=(N+a*q_over_factor)%p
    return N*pow(prod_nonzero,-1,p)%p


def block_and_residue_tests():
    counts=Counter()
    for trial in range(220):
        p=rng.choice([11,17,29,43]);n=rng.randrange(2,6);d=rng.randrange(1,6);w=rng.randrange(1,min(p,24)+1);k=3
        f=[]
        for j in range(n):
            deg=rng.randrange(d+1)
            if rng.randrange(6)==0:g=[0];counts['zero_fibers']+=1
            elif rng.randrange(2)==0:g=[rng.randrange(p) for _ in range(deg+1)]
            else:
                roots=[rng.randrange(p) for _ in range(deg)];g=roots_poly(p,roots)
                counts['repeated_input_roots']+=len(roots)!=len(set(roots))
                # Keep the polynomial monic in this generated-root branch.
            f.append(g+[0]*(d+1-len(g)))
        tags=[[[rng.randrange(p) for _ in range(d+1)] for j in range(n)] for _ in range(k)]
        H=rng.randrange(1,w+1)
        R=ProductDual(p,n,0);v=R.block_norm([tuple((row[i],) for row in f) for i in range(d+1)],w,H)
        counts['blocks']+=1;counts['tail']+=bool(w%H)
        for j in range(n):
            v0=prod(pe(f[j],u,p) for u in range(w))%p
            assert v[j][0]==v0
            assert (1-pow(v0,p-1,p))%p==int(any(pe(f[j],u,p)==0 for u in range(w)))
            counts['projected_masks']+=1
        # Work at an actual extension point outside the base field.
        R=ProductDual(p,n,k);F=R.F;zeta=rng.randrange(p)+p*rng.randrange(1,p)
        coeff=[]
        for i in range(d+1):
            coeff.append(tuple((F.add(row[i],zeta) if i==0 else row[i],)+tuple(tags[t][j][i] for t in range(k)) for j,row in enumerate(f)))
        D=R.block_norm(coeff,w,H)
        counts['dual_norms']+=1
        for j,row in enumerate(f):
            vals=[F.add(pe(row,u,p),zeta) for u in range(w)]
            direct=1
            for a in vals:direct=F.mul(direct,a)
            assert direct==D[j][0] and direct
            for t in range(k):
                s=F.sm(F.mul(pe(tags[t][j],u,p),F.inv(v)) for u,v in enumerate(vals))
                assert F.mul(D[j][t+1],F.inv(direct))==s
                counts['reciprocal_tags']+=1
        vals=[pe(row,u,p) for row in f for u in range(w)]
        for t in range(k):
            wt=[pe(tags[t][j],u,p) for j in range(n) for u in range(w)]
            assert residue_from_distinct_poles(F,vals,wt)==sum(a for a,vv in zip(wt,vals) if vv==0)%p
            counts['global_residues']+=1
    print('Block and rational trace tests:',dict(counts))
    return counts


def prem(f,g,p):
    d=len(g)-1;assert g[-1]==1
    r=list(f)+[0]*max(0,d-len(f))
    for j in range(len(r)-1,d-1,-1):
        for i in range(d):r[j-d+i]=(r[j-d+i]-r[j]*g[i])%p
    return r[:d]

def pm(f,g,p):
    r=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        for j,b in enumerate(g):r[i+j]=(r[i+j]+a*b)%p
    return r

def projective_tests():
    C=Counter()
    for _ in range(300):
        p=rng.choice([11,17,29,43]);d=rng.randrange(1,7);w=rng.randrange(1,p+1);alpha=rng.randrange(p)
        f=[rng.randrange(p) for _ in range(rng.randrange(d+1)+1)]
        if rng.randrange(10)==0:f=[0]
        f=f+[0]*(d+1-len(f));a=pe(f,alpha,p)
        from math import comb
        # V^d f(alpha+1/V)
        t=[0]*(d+1)
        for j in range(d+1):
            for k in range(j+1):t[d-j+k]=(t[d-j+k]+f[j]*comb(j,k)*pow(alpha,k,p))%p
        assert t[d]==a
        g=[t[j]*pow(a,d-1-j,p)%p for j in range(d)]+[1]
        R=[1]
        for u in range(w):R=prem(pm(R,[a,(alpha-u)%p],p),g,p)
        M=[]
        for j in range(d):M.append(prem([0]*j+R,g,p))
        got=det(list(map(list,zip(*M))),p)
        expected=pow(a,(d-1)*w,p)*prod(pe(f,u,p) for u in range(w))%p
        assert got==expected
        C['projective_determinants']+=1;C['degree_drops']+=f[-1]==0;C['a_zero']+=a==0
        if a and alpha>=w:
            assert (got==0)==any(pe(f,u,p)==0 for u in range(w));C['valid_chart_zero_tests']+=1
    print('Projective determinant identity:',dict(C));return C


def recover_projected(p,rest,J,cs):
    """Newton/Vandermonde recovery using only projected-mask traces."""
    n=len(rest[0]);H=sum(J);assert H<p
    if not H:return set()
    yy=[sum(a*b for a,b in zip(cs,z))%p for z in rest]
    M=[sum(j*pow(y,k,p) for y,j in zip(yy,J))%p for k in range(H+1)]
    tag=[[sum(j*z[i]*pow(y,k,p) for z,y,j in zip(rest,yy,J))%p for k in range(H)] for i in range(n)]
    coeff=[1]
    for k in range(1,H+1):coeff.append(-sum(coeff[k-j]*M[j] for j in range(1,k+1))*pow(k,-1,p)%p)
    ys=[v for v in range(p) if pe(coeff[::-1],v,p)==0]
    if len(ys)!=H:return None
    V=[[pow(y,k,p) for y in ys] for k in range(H)]
    return set(zip(*(solve(V,t,p) for t in tag)))


def ec_add(P,Q,p,a):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and (y+v)%p==0:return None
    t=((3*x*x+a)*pow(2*y%p,-1,p) if x==u else (v-y)*pow((u-x)%p,-1,p))%p
    z=(t*t-x-u)%p
    return (z,(t*(x-z)-y)%p)

def actual_multi_hit_controls():
    C=Counter()
    p,a,c=13,4,1
    points=[(x,y) for x in range(p) for y in range(p) if y*y%p==(x*x*x+a*x+c)%p]
    assert len(points)+1==19
    F=semaev(p,a,c)
    for T in points:
        for Q in points:
            if T[0]==Q[0]:continue
            x=ec_add(T,Q,p,a);y=ec_add(T,(Q[0],-Q[1]%p),p,a)
            if x is None or y is None or x[0]==y[0]:continue
            base=x[0];slope=(y[0]-base)%p
            row=interp([F(((base+slope*u)%p,Q[0]),T[0]) for u in range(3)],p)
            R=ProductDual(p,1,0);D=R.block_norm([((v,),) for v in row],2)
            assert D[0][0]==0
            assert all(pe(row,u,p)==0 for u in [0,1])
            # Independently verify each selected abscissa has a signed point completing Q or -Q to T.
            for u in [0,1]:
                lifts=[P for P in points if P[0]==(base+slope*u)%p]
                assert any(ec_add(P,V,p,a)==T for P in lifts for V in (Q,(Q[0],-Q[1]%p)))
                C['actual_two_root_addresses']+=1
            C['actual_two_root_fibers']+=1
    print('Controlled genuine elliptic multiple-hit fibers:',dict(C));return C

def elliptic_projection_tests():
    C=Counter()
    for p,a,c,T in [(101,2,3,3),(257,1,1,0),(263,4,1,0)]:
        F=semaev(p,a,c)
        for m in [2,3,4]:
            d=2**(m-1)
            for it in range(12):
                w=rng.randrange(2,6);rest=list(product(*(rng.sample(range(p),rng.randrange(2,4)) for _ in range(m-1))))
                b=rng.randrange(p);t=rng.randrange(1,p)
                rows=[interp([F(((b+t*u)%p,)+z,T) for u in range(d+1)],p) for z in rest]
                root_sets=[set(u for u in range(w) if F(((b+t*u)%p,)+z,T)==0) for z in rest]
                R=ProductDual(p,len(rest),0);D=R.block_norm([tuple((f[i],) for f in rows) for i in range(d+1)],w)
                J=[(1-pow(v[0],p-1,p))%p for v in D]
                assert J==[int(bool(x)) for x in root_sets]
                C['elliptic_grids']+=1;C['fibers']+=len(rest);C['addresses']+=w*len(rest)
                C['multi_hit_fibers']+=sum(len(x)>1 for x in root_sets)
                C['positive_fibers']+=sum(J);C['whole_slot_fibers']+=sum(len(x)==w for x in root_sets)
                true={(u,)+z for z,rs in zip(rest,root_sets) for u in rs}
                # Two readouts: a random projection and deliberate all-zero challenge.
                for cs in [[rng.randrange(p) for _ in range(m-1)],[0]*(m-1)]:
                    rec=recover_projected(p,rest,J,cs)
                    if rec is None:
                        assert len({sum(a*b for a,b in zip(cs,z))%p for z,j in zip(rest,J) if j})<sum(J)
                        C['projection_collisions']+=1;continue
                    assert rec=={z for z,j in zip(rest,J) if j}
                    # Only the recovered fibers are inspected here; tiny fixture factors by scanning.
                    got={(u,)+z for z in rec for u in range(w) if pe(rows[rest.index(z)],u,p)==0}
                    assert got==true
                    C['successful_readouts']+=1;C['recovered_addresses']+=len(got)
    # Explicit edge cases including an all-zero polynomial fiber and multiple roots in one fiber.
    p=17;rest=[(0,),(1,),(2,)];w=5;rows=[[0,0,0],roots_poly(p,[1,3]),[1,0,0]]
    R=ProductDual(p,3,0);D=R.block_norm([tuple((f[i],) for f in rows) for i in range(3)],w)
    J=[(1-pow(v[0],p-1,p))%p for v in D];assert J==[1,1,0]
    rec=recover_projected(p,rest,J,[1]);assert rec=={(0,),(1,)}
    C['explicit_edge_addresses']=sum(1 for z in rec for u in range(w) if pe(rows[z[0]],u,p)==0);assert C['explicit_edge_addresses']==7
    print('Actual elliptic projected traces:',dict(C));return C


def rank(A,p):
    if not A:return 0
    A=[list(a) for a in A];r=0
    for j in range(len(A[0])):
        k=next((k for k in range(r,len(A)) if A[k][j]%p),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];iv=pow(A[r][j]%p,-1,p);A[r]=[x*iv%p for x in A[r]]
        for k in range(len(A)):
            if k!=r:
                q=A[k][j];A[k]=[(a-q*b)%p for a,b in zip(A[k],A[r])]
        r+=1
        if r==len(A):break
    return r


def oracle_collision(A,p):
    """Solve A n=0, n0=1 by exact row reduction, free entries set zero."""
    n=len(A[0]);a=[list(row)+[0] for row in A]+[[1]+[0]*(n-1)+[1]];r=0;pivs=[]
    for j in range(n):
        k=next((k for k in range(r,len(a)) if a[k][j]%p),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r];iv=pow(a[r][j]%p,-1,p);a[r]=[x*iv%p for x in a[r]]
        for k in range(len(a)):
            if k!=r:
                q=a[k][j];a[k]=[(v-q*w)%p for v,w in zip(a[k],a[r])]
        pivs.append(j);r+=1
        if r==len(a):break
    assert not any(not any(row[:-1]) and row[-1] for row in a)
    x=[0]*n
    for k,j in enumerate(pivs):x[j]=a[k][-1]
    assert x[0]==1 and all(sum(u*v for u,v in zip(row,x))%p==0 for row in A)
    return x


def rational_oracle_tests():
    C=Counter()
    for _ in range(500):
        p=rng.choice([17,29,43,101]);F=Fp2(p);s=rng.randrange(3,min(p,13));S=[0]+rng.sample(range(1,p),s-1);A=[[1]*s]
        budget=0;queries=[]
        # A mixture of base/extension values and derivative (Hasse-jet) queries.
        while budget<s-2 and rng.randrange(5)!=0:
            ext=rng.randrange(2)
            e=2 if ext else 1
            jmax=min(rng.randrange(3),(s-2-budget)//e-1)
            if jmax<0:break
            if ext:t=rng.randrange(p)+p*rng.randrange(1,p)
            else:t=rng.choice([a for a in range(p) if (-a)%p not in S])
            for j in range(jmax+1):
                row=[F.mul((-1)**j%p,F.pw(F.inv(F.add(v,t)),j+1)) for v in S]
                A.append([x%p for x in row]);
                if ext:A.append([x//p for x in row])
                C['jet_rows']+=e
            queries.append((t,jmax));budget+=e*(jmax+1)
        assert budget<s-1
        r=rank(A,p);assert rank(A+[[1]+[0]*(s-1)],p)==r+1
        x=oracle_collision(A,p)
        for t,jmax in queries:
            for j in range(jmax+1):
                assert F.sm(F.mul(vv,F.mul((-1)**j%p,F.pw(F.inv(F.add(v,t)),j+1))) for vv,v in zip(x,S))==0
                C['independent_zero_transcripts']+=1
        C['dual_collisions']+=1
        # Matching upper control: constant plus s-1 distinct admissible base-field resolvents.
        if p-s>=s-1:
            ts=[t for t in range(p) if (-t)%p not in S][:s-1]
            assert rank([[1]*s]+[[pow((t+v)%p,-1,p) for v in S] for t in ts],p)==s
            C['full_rank_controls']+=1
    print('Rational readout obstruction:',dict(C));return C


def isprime(n):return n>1 and all(n%k for k in range(2,isqrt(n)+1))
def prime_curve(p):
    for a in range(1,p):
        for c in range(1,min(p,9)):
            if (4*a*a*a+27*c*c)%p==0:continue
            N=1;xs=[]
            for x in range(p):
                y=(x*x*x+a*x+c)%p;leg=pow(y,(p-1)//2,p)
                if not y:N+=1
                elif leg==1:N+=2;xs.append(x)
            if isprime(N) and N%2:return a,c,N,xs
    raise AssertionError('curve not found')


def spectrum_tests():
    C=Counter();samples=[]
    # Actual specialized Semaev polynomials in A,X, unprojected prime-order curves.
    # We test nonassociation by distinct polynomial signatures at a sufficient grid.
    for p in [101,137]:
        a,c,r,xs=prime_curve(p);F=semaev(p,a,c);m=3;S=[list(range(2*j,2*j+2)) for j in range(m)]
        for elim in [(0,),(0,1)]:
            other=[i for i in range(m) if i not in elim];I=list(product(*(S[i] for i in elim)));Z=list(product(*(S[i] for i in other)))
            DA=len(I)*m*2**(m-1);DX=len(I)*2**(m-1);assert max(DA,DX)<p
            def addr(u,z,A):
                s=[0]*m
                for i,v in zip(elim,u):s[i]=(v+A)%p
                for i,v in zip(other,z):s[i]=(v+A)%p
                return tuple(s)
            def D(A,X,z):return prod(F(addr(u,z,A),X) for u in I)%p
            signatures=[tuple(D(A,X,z) for A in range(DA+1) for X in range(DX+1)) for z in Z]
            assert len(set(signatures))==len(Z)
            C['distinct_polynomial_pairs']+=len(Z)*(len(Z)-1)//2
            C['signature_evaluations']+=len(Z)*(DA+1)*(DX+1)
            # Compare ordered collision expectation with the exact bivariate degree bound.
            E=0;nobs=0;success=0
            for A in range(p):
                for X in xs:
                    vs=[D(A,X,z) for z in Z];energy=sum(k*k for k in Counter(vs).values())
                    assert len(Z)**2<=energy*len(set(vs));E+=energy;nobs+=1;success+=0 in vs
            N=len(Z)
            # p*len(xs) sample pairs, epsilon=DA/p+DX/len(xs).
            assert E<=N*nobs+N*(N-1)*(DA*len(xs)+DX*p)
            C['spectrum_observations']+=nobs;C['positive_algebraic_observations']+=success
            samples.append((p,len(I),len(Z),round(E/nobs,5)))
    print('Projected spectrum checks:',dict(C));print('Mean ordered collision energies (p, eliminated width, lanes, energy):',samples)
    return C

if __name__=='__main__':
    assert __debug__
    block_and_residue_tests();projective_tests();elliptic_projection_tests();actual_multi_hit_controls();rational_oracle_tests();spectrum_tests()
    print('All exact assertions passed. No fast collective evaluator or end-to-end timing is implemented.')
