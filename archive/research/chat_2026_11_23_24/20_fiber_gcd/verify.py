#!/usr/bin/env python3
"""Exact component checks; product-ring lanes are explicit, not a fast trace oracle."""
from itertools import product
from math import prod, isqrt, comb
from random import Random
from functools import lru_cache
from collections import Counter
rng=Random(20261124)

class Ring:
    def __init__(self,p,n):
        self.p=p;self.n=n;self.z=(0,)*n;self.o=(1,)*n
    def c(self,a):return (a%self.p,)*self.n
    def add(self,a,b):return tuple((x+y)%self.p for x,y in zip(a,b))
    def neg(self,a):return tuple(-x%self.p for x in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def mul(self,a,b):return tuple(x*y%self.p for x,y in zip(a,b))
    def pw(self,a,k):return tuple(pow(x,k,self.p) for x in a)
    def sm(self,vals):
        out=self.z
        for v in vals:out=self.add(out,v)
        return out
    def pm(self,f,g):
        out=[self.z]*(len(f)+len(g)-1)
        for i,a in enumerate(f):
            for j,b in enumerate(g):out[i+j]=self.add(out[i+j],self.mul(a,b))
        return out
    def shift(self,f,c): # f(X+c); schoolbook Horner, fast version uses Taylor convolution
        out=[f[-1]]
        for a in f[-2::-1]:
            out=self.pm(out,[self.c(c),self.o]);out[0]=self.add(out[0],a)
        return out
    def rem(self,f,g):
        d=len(g)-1;assert g[-1]==self.o
        if d==0:return []
        a=list(f)+[self.z]*max(0,d-len(f))
        for j in range(len(a)-1,d-1,-1):
            c=a[j]
            for i in range(d):a[j-d+i]=self.sub(a[j-d+i],self.mul(c,g[i]))
        return a[:d]
    def mask(self,f):
        """Return degree masks in ascending order and all-zero mask."""
        out=[self.z]*len(f);left=self.o
        for i in reversed(range(len(f))):
            nz=self.pw(f[i],self.p-1)
            out[i]=self.mul(left,nz);left=self.mul(left,self.sub(self.o,nz))
        return out,left
    def monic(self,f,k):
        iv=self.pw(f[k],self.p-2)
        return [self.mul(a,iv) for a in f[:k]]+[self.o]
    def modproduct(self,roots,f):
        g=[self.o]
        for x in roots:g=self.rem(self.pm(g,[self.c(-x),self.o]),f)
        return g
    def block(self,f,w,H=None):
        """Reduce falling factorial Q(X) modulo monic f using shifted moduli.
        A true subproduct/remainder tree; polynomial products are schoolbook here."""
        d=len(f)-1;assert d>=1
        if H is None:H=min(w,max(1,isqrt(d*w)))
        def poly_roots(L):
            out=[self.o]
            for j in range(L):out=self.pm(out,[self.c(-j),self.o])
            return out
        P=poly_roots(H);J=w//H
        mods=[self.shift(f,q*H) for q in range(J)]
        def tree(L):
            if len(L)==1:return (L[0],None,None)
            h=len(L)//2;a=tree(L[:h]);b=tree(L[h:]);return (self.pm(a[0],b[0]),a,b)
        def walk(g,t):
            r=self.rem(g,t[0]);a,b=t[1:]
            if a is None:return [r]
            return walk(r,a)+walk(r,b)
        rs=walk(P,tree(mods));out=[self.o]
        for q,r in enumerate(rs):out=self.rem(self.pm(out,self.shift(r,-q*H)),f)
        if w%H:out=self.rem(self.pm(out,self.shift(poly_roots(w%H),-J*H)),f)
        return out+[self.z]*(d-len(out))
    def gcd(self,f,r):
        """Branch-free Euclidean algorithm on explicit product ring, fixed degree bound."""
        d=len(f)-1;A=list(f);B=list(r)+[self.z]*(d+1-len(r))
        for _ in range(d+1):
            masks,zero=self.mask(B)
            # In each lane either freeze A or select exactly one monic divisor.
            nextA=[self.mul(zero,a) for a in A];nextB=[self.z]*(d+1)
            for k,e in enumerate(masks):
                g=self.monic(B,k);rem=self.rem(A,g)
                for j,x in enumerate(g):nextA[j]=self.add(nextA[j],self.mul(e,x))
                for j,x in enumerate(rem):nextB[j]=self.add(nextB[j],self.mul(e,x))
            A,B=nextA,nextB
        assert all(x==self.z for x in B)
        return A
    def all_gcd(self,F,w):
        """Arbitrary degree including zero lanes. Return low-degree gcds, plus
        mask for whole-slot fibers (their moments are handled separately)."""
        d=len(F)-1;mask,zero=self.mask(F)
        G=[mask[0]]+[self.z]*d
        for k in range(1,d+1):
            f=self.monic(F,k);q=self.block(f,w);g=self.gcd(f,q)
            for j,x in enumerate(g):G[j]=self.add(G[j],self.mul(mask[k],x))
        return G,zero
    def moments(self,G,zero,w,K):
        """Power sums in index coordinate, including zero polynomial fibers."""
        d=len(G)-1;masks,_=self.mask(G)
        # Reverse a variable-degree monic polynomial into 1+O(t).
        H=self.sm([self.mul(self.c(j),e) for j,e in enumerate(masks)])
        C=[self.o]+[self.sm(self.mul(masks[j],G[j-i]) for j in range(i,d+1)) for i in range(1,d+1)]
        P=[H]
        for j in range(1,K+1):
            # -t C'(t)/C(t) with constant coefficient 1.
            c=self.mul(self.c(-j),C[j]) if j<=d else self.z
            for i in range(1,min(d,j-1)+1):c=self.sub(c,self.mul(C[i],P[j-i]))
            P.append(c)
        assert K+1<self.p
        W=[]
        for j in range(K+1):
            W.append((pow(w,j+1,self.p)-sum(comb(j+1,k)*W[k] for k in range(j)))*pow(j+1,-1,self.p)%self.p)
        return [self.add(v,self.mul(zero,self.c(W[j]))) for j,v in enumerate(P)]

def pe(f,x,p):
    s=0
    for a in reversed(f):s=(s*x+a)%p
    return s

def roots_poly(p,L):
    f=[1]
    for x in L:
        g=[0]*(len(f)+1)
        for i,c in enumerate(f):g[i]=(g[i]-x*c)%p;g[i+1]=(g[i+1]+c)%p
        f=g
    return f

def ring_tests():
    checks=lanes=tail=var=whole=repeated=0
    for _ in range(200):
        p=rng.choice([11,17,31,101]);d=rng.randrange(1,min(p,7));n=rng.randrange(2,7);w=rng.randrange(1,p+1)
        R=Ring(p,n);H=rng.randrange(1,w+1)
        rows=[[rng.randrange(p) for _ in range(d)]+[1] for _ in range(n)]
        F=list(zip(*rows));q=R.block(F,w,H);direct=R.modproduct(range(w),F)
        assert q==direct;checks+=1;tail+=w%H!=0
        # Independent mixed degree inputs, with duplicate roots and zero polynomial lanes.
        fs=[]
        for j in range(n):
            deg=rng.randrange(d+1)
            if rng.randrange(5)==0:ff=[0];whole+=1
            else:
                L=[rng.randrange(p) for k in range(deg)];repeated+=len(set(L))<len(L)
                ff=roots_poly(p,L);c=rng.randrange(1,p);ff=[v*c%p for v in ff]
            fs.append(ff+[0]*(d+1-len(ff)))
        G,z=R.all_gcd(list(zip(*fs)),w);P=R.moments(G,z,w,9)
        var+=len({max((i for i,v in enumerate(f) if v),default=-1) for f in fs})>1
        for j,f in enumerate(fs):
            true=[x for x in range(w) if pe(f,x,p)==0];lanes+=1
            if not z[j]:assert [x[j] for x in G]==roots_poly(p,true)+[0]*(d-len(true))
            else:assert not any(f) and len(true)==w
            assert [x[j] for x in P]==[sum(pow(x,k,p) for x in true)%p for k in range(10)]
    print(f'Block remainders: {checks} exact product-ring comparisons ({tail} nonempty tails).')
    print(f'Masked gcd/moments: {lanes} scalar lanes, {var} mixed-degree fixtures, {whole} whole-slot lanes, {repeated} repeated-root input lanes.')

# independent declared-degree resultant and finite-field interpolation, for actual Semaev tests.
def det(M,p):
    M=[list(r) for r in M];n=len(M);out=1
    for j in range(n):
        k=next((k for k in range(j,n) if M[k][j]%p),None)
        if k is None:return 0
        if k!=j:M[k],M[j]=M[j],M[k];out=-out
        v=M[j][j]%p;out=out*v%p;iv=pow(v,-1,p)
        for k in range(j+1,n):
            q=M[k][j]*iv%p
            for l in range(j+1,n):M[k][l]=(M[k][l]-q*M[j][l])%p
            M[k][j]=0
    return out%p

def res(f,g,p):
    d,e=len(f)-1,len(g)-1;M=[[0]*(d+e) for _ in range(d+e)]
    for i in range(e):M[i][i:i+d+1]=f[::-1]
    for i in range(d):M[e+i][i:i+e+1]=g[::-1]
    return det(M,p)

def s3(v,z,a,c,p):return [( (v*z-a)**2-4*c*(v+z))%p,(-2*((v+z)*(v*z+a)+2*c))%p,(v-z)**2%p]

def interp(vals,p):
    n=len(vals);diff=list(vals);co=[0]*n;basis=[1];fac=1
    for j in range(n):
        k=diff[0]*pow(fac,-1,p)%p
        for i,v in enumerate(basis):co[i]=(co[i]+k*v)%p
        diff=[(v-u)%p for u,v in zip(diff,diff[1:])]
        if j+1<n:
            b=[0]*(len(basis)+1)
            for i,v in enumerate(basis):b[i]=(b[i]-j*v)%p;b[i+1]=(b[i+1]+v)%p
            basis=b;fac=fac*(j+1)%p
    return co

def semaev(p,a,c):
    @lru_cache(None)
    def partial(xs):
        if len(xs)==1:return ((-xs[0])%p,1)
        if len(xs)==2:return tuple(s3(*xs,a,c,p))
        prev=partial(xs[:-1]);d=2**(len(xs)-1)
        return tuple(interp([res(prev,s3(xs[-1],y,a,c,p),p) for y in range(d+1)],p))
    return lambda xs,T:pe(partial(tuple(xs)),T,p)

def solve(M,y,p):
    n=len(y);a=[list(row)+[v] for row,v in zip(M,y)]
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j]%p);a[k],a[j]=a[j],a[k]
        u=pow(a[j][j],-1,p);a[j]=[x*u%p for x in a[j]]
        for k in range(n):
            if k==j:continue
            c=a[k][j];a[k]=[(x-c*z)%p for x,z in zip(a[k],a[j])]
    return [r[-1] for r in a]

def projected_recover(roots,p,cs,P_by_lane,rest,w,K):
    m=len(cs);H=sum(P_by_lane[0])%p
    assert H==len(roots)%p
    if H>K:return 'overflow'
    if not H:return 'empty'
    # powers of projected address c_0 u + sum c_j z_j
    def local(j,zidx,extra):
        z=rest[zidx];c=sum(a*b for a,b in zip(cs[1:],z))%p
        return sum(comb(j,k)*pow(cs[0],k,p)*pow(c,j-k,p)*P_by_lane[k+extra][zidx] for k in range(j+1))%p
    M=[sum(local(j,k,0) for k in range(len(rest)))%p for j in range(H+1)]
    tag=[[sum(local(j,k,1) for k in range(len(rest)))%p for j in range(H)]]
    for i in range(m-1):tag.append([sum(rest[k][i]*local(j,k,0) for k in range(len(rest)))%p for j in range(H)])
    assert M==[sum(pow(sum(a*b for a,b in zip(cs,s))%p,j,p) for s in roots)%p for j in range(H+1)]
    for i in range(m):assert tag[i]==[sum(s[i]*pow(sum(a*b for a,b in zip(cs,s))%p,j,p) for s in roots)%p for j in range(H)]
    # Newton identities, coefficient order from leading coefficient down.
    coef=[1]
    for j in range(1,H+1):coef.append(-sum(coef[j-k]*M[k] for k in range(1,j+1))*pow(j,-1,p)%p)
    pr=[x for x in range(p) if pe(coef[::-1],x,p)==0]
    if len(pr)!=H:return 'collision'
    V=[[pow(x,j,p) for x in pr] for j in range(H)]
    recovered=set(zip(*(solve(V,t,p) for t in tag)))
    assert recovered==set(roots)
    return 'ok'

def elliptic_tests():
    outcomes=Counter();slices=degree_drops=addresses=checks=0;hits=0
    for p,a,c,T in [(101,2,3,3),(257,1,1,0),(263,4,1,0)]:
        assert (4*a**3+27*c*c)%p and pow((T**3+a*T+c)%p,(p-1)//2,p)==1
        F=semaev(p,a,c)
        for m in (2,3,4):
            d=2**(m-1)
            for q in range(15):
                w=rng.choice([2,3,4,5]);rest=list(product(*[rng.sample(range(p),rng.choice([2,3])) for i in range(m-1)]));base=rng.randrange(p);slope=rng.randrange(1,p)
                # F is viewed in interval index u, with arbitrary actual affine leaf abscissa base+slope*u.
                fs=[interp([F(((base+slope*u)%p,)+z,T) for u in range(d+1)],p) for z in rest]
                R=Ring(p,len(rest));G,zero=R.all_gcd(list(zip(*fs)),w);assert zero==R.z
                P=R.moments(G,zero,w,9);roots=[]
                for k,z in enumerate(rest):
                    slices+=1;degree_drops+=fs[k][-1]==0
                    true=[u for u in range(w) if F(((base+slope*u)%p,)+z,T)==0]
                    for u in range(w):assert pe(fs[k],u,p)==F(((base+slope*u)%p,)+z,T)
                    roots.extend((u,)+z for u in true);addresses+=w
                    assert [v[k] for v in P]==[sum(pow(u,j,p) for u in true)%p for j in range(10)]
                assert len(roots)<p
                # Also force some collisions in the projection to test failure detection.
                for cs in ([rng.randrange(p) for i in range(m)],[0]*m):
                    outcome=projected_recover(roots,p,cs,P,rest,w,8);outcomes[outcome]+=1
                    if outcome=='ok':hits+=len(roots)
                checks+=1
    print(f'Elliptic grid projections: {checks} actual Semaev grids, {slices} slices ({degree_drops} degree drops), {addresses} literal address comparisons.')
    print(f'Trace recovery: {dict(outcomes)}; {hits} recovered addresses (with multiplicity across tests).')

# Certificate-degree test: exact Lagrange reconstruction from any grid values.
def lagrange(p,S,x):
    g=roots_poly(p,[v for v in S if v!=x]);iv=pow(prod(x-v for v in S if v!=x)%p,-1,p)
    return [v*iv%p for v in g]

def weights(S,grid,p):return [pow(prod(prod(x-v for v in Si if v!=x) for x,Si in zip(s,S))%p,-1,p) for s in grid]

def interp_grid(p,S,grid,vals):
    E=list(product(*[range(len(Si)) for Si in S]));ls=[[lagrange(p,Si,x) for x in Si] for Si in S]
    index=[{x:j for j,x in enumerate(Si)} for Si in S]
    return {e:sum(v*prod(ls[i][index[i][s[i]]][e[i]] for i in range(len(S))) for s,v in zip(grid,vals))%p for e in E}

def cert_tests():
    checks=dual=singular=generic=nonunit=0
    for _ in range(100):
        p=101;m=rng.choice([2,3]);S=[rng.sample(range(p),rng.choice([2,3,4])) for i in range(m)];grid=list(product(*S));W=len(grid);s=rng.choice(grid);k0=grid.index(s);D=sum(len(x)-1 for x in S);wt=weights(S,grid,p)
        F=[0 if x==s else rng.randrange(1,p) for x in grid]
        def R(f):return [[sum(wt[k]*(x[i]-s[i])*(x[j]-s[j])*pow(f[k],-1,p) for k,x in enumerate(grid) if k!=k0)%p for j in range(m)] for i in range(m)]
        r=R(F);singular+=not det(r,p)
        for c in [[rng.randrange(p) for j in range(m)] for _ in range(5)]+[[int(i==j) for j in range(m)] for i in range(m)]:
            v=[0 if k==k0 else sum(a*(xj-yj) for a,xj,yj in zip(c,x,s))*pow(F[k],-1,p)%p for k,x in enumerate(grid)]
            v[k0]=-sum(w*u for w,u in zip(wt,v))*pow(wt[k0],-1,p)%p
            C=interp_grid(p,S,grid,v);de=max((sum(e) for e,a in C.items() if a),default=-1)
            assert de<=D-1
            assert (de<=D-2)==all(sum(a*b for a,b in zip(row,c))%p==0 for row in r)
            checks+=1
        # Construct the preconditioning nonzero determinant witness by a formal
        # pole limit: isolate the m single-coordinate neighbors plus s.
        t=[next(v for v in Si if v!=si) for si,Si in zip(s,S)]
        neighbors={tuple(t[i] if i==j else s[i] for i in range(m)) for j in range(m)}
        forms=lambda x:[(x[i]-s[i])*(x[i]-t[i])%p for i in range(m)]+[(x[i]-s[i])*(x[j]-s[j])%p for i in range(m) for j in range(i+1,m)]
        assert {x for x in grid if not any(forms(x))}==neighbors|{s}
        pole=[[sum(wt[k]*(x[i]-s[i])*(x[j]-s[j])*pow(F[k],-1,p) for k,x in enumerate(grid) if x in neighbors)%p for j in range(m)] for i in range(m)]
        assert det(pole,p);dual+=1
        # Actual random quadratic unit (one trial; failures not suppressed).
        mons=[e for e in product(range(3),repeat=m) if sum(e)<=2];coef=[rng.randrange(p) for e in mons]
        G=[sum(c*prod(pow(xi,ei,p) for xi,ei in zip(x,e)) for c,e in zip(coef,mons))%p for x in grid]
        if not all(G):nonunit+=1;continue
        generic+=bool(det(R([f*g%p for f,g in zip(F,G)]),p))
    print(f'Image certificates: {checks} exact degree/reciprocal-matrix comparisons; {singular} singular unpreconditioned matrices.')
    print(f'Unit preconditioning: {dual} exact diagonal pole certificates; {generic} nonsingular random unit trials, {nonunit} nonunit trials (not retried).')

def actual_certificate_tests():
    tests=attempts=0
    for p,m in [(257,2),(1009,3)]:
        F0=semaev(p,2,1);d=2**(m-1);w=d+2;done=0
        while done<10 and attempts<2000:
            attempts+=1;S=[rng.sample(range(p),w) for i in range(m)];grid=list(product(*S));F=[F0(s,0) for s in grid]
            if F.count(0)!=1:continue
            z=F.index(0);s=grid[z];wt=weights(S,grid,p)
            for i in range(m):
                v=[0 if k==z else (x[i]-s[i])*pow(F[k],-1,p)%p for k,x in enumerate(grid)]
                v[z]=-sum(w*u for w,u in zip(wt,v))*pow(wt[z],-1,p)%p
                C=interp_grid(p,S,grid,v)
                for j in range(m):assert max(e[j] for e,a in C.items() if a)>=w-d
                tests+=1
            done+=1
        assert done==10
    print(f'Actual Semaev image certificates: {tests} affine generators on 20 singleton grids, all separate-degree bounds verified ({attempts} attempted grids).')

if __name__=='__main__':
    assert __debug__
    ring_tests();elliptic_tests();cert_tests();actual_certificate_tests()
