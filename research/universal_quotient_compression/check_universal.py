"""Exact abstract finite-group checks of quotient geometry and integral transport.

No elliptic-curve solver or tensor-contraction implementation. The finite
range-avoidance oracle in the 2D test exhausts the toy support; it does not
implement the unresolved efficient avoidance algorithm.
"""
from itertools import product
from math import prod, ceil
from random import Random
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form


def add(a,b,p): return tuple((x+y)%p for x,y in zip(a,b))
def sub(a,b,p): return tuple((x-y)%p for x,y in zip(a,b))
def scale(k,a,p): return tuple(k*x%p for x in a)
def dot(a,b,p): return sum(x*y for x,y in zip(a,b))%p

def rref(rows,p,n):
    B={}
    for rr in rows:
        a=[x%p for x in rr]
        for j in sorted(B):
            if a[j]:
                c=a[j];a=[(x-c*y)%p for x,y in zip(a,B[j])]
        j=next((i for i,x in enumerate(a) if x),None)
        if j is None: continue
        c=pow(a[j],-1,p);a=[x*c%p for x in a]
        for k,b in list(B.items()):
            if b[j]:
                c=b[j];B[k]=[(x-c*y)%p for x,y in zip(b,a)]
        B[j]=a
    return B

def normal(rows,p,n):
    B=rref(rows,p,n)
    j=next((i for i in range(n) if i not in B),None)
    assert j is not None
    a=[0]*n;a[j]=1
    for k,b in B.items():a[k]=-b[j]%p
    assert all(dot(a,b,p)==0 for b in rows)
    return tuple(a)

def quotient(U,p,n):
    B=rref(U,p,n);free=[j for j in range(n) if j not in B]
    def project(v):
        return tuple((v[j]-sum(v[i]*a[j] for i,a in B.items()))%p for j in free)
    def lift(v):
        x=[0]*n
        for i,j in enumerate(free):x[j]=v[i]%p
        return tuple(x)
    return len(free),project,lift

class Toy:
    def __init__(self,p,blocks,w):
        self.p=p;self.B=blocks;self.m=len(blocks);self.w=w
        self.n=sum(len(F)-1 for F in blocks)
        self.lam=tuple((a-F[0])%p for F in blocks for a in F[1:]);self.g0=sum(w*F[0] for F in blocks)%p
        j=0;self.cols=[]
        for F in blocks:
            A=[(0,)*self.n]
            for _ in F[1:]:
                a=[0]*self.n;a[j]=1;j+=1;A.append(tuple(a))
            self.cols.append(A)
        self.slots=[i for i in range(self.m) for _ in range(w)]
        self.rep=self.represent(self.slots)
        assert len(self.rep)==p
    def represent(self,slots):
        rep={0:()}
        for i in slots:
            nxt={}
            for a,T in rep.items():
                for j,v in enumerate(self.B[i]):nxt.setdefault((a+v)%self.p,T+(j,))
            rep=nxt
        return rep
    def inc(self,X,slots=None):
        if slots is None:slots=self.slots
        a=[0]*self.n
        for i,j in zip(slots,X):
            if j:a[self.cols[i][j].index(1)]+=1
        return tuple(a)
    def val(self,X,slots=None):
        if slots is None:slots=self.slots
        return sum(self.B[i][j] for i,j in zip(slots,X))%self.p
    def calibrated(self):
        a=sub(self.inc(self.rep[1]),self.inc(self.rep[0]),self.p)
        assert dot(a,self.lam,self.p)==1
        return a
    def row(self,z,lift):
        x=lift(z);g=(self.g0+dot(x,self.lam,self.p))%self.p
        c=self.inc(self.rep[g]);row=sub(x,c,self.p)
        assert dot(row,self.lam,self.p)==0
        return row

def partition_cover(toy,t):
    groups=[];k=0
    for C in toy.cols:
        size=ceil(len(C)/t)
        G=[C[i:i+size] for i in range(0,len(C),size)]
        groups.append(G);k+=toy.w*max(len(g)-1 for g in G)
    return groups,k,prod(len(G)**toy.w for G in groups)

def find_outside(groups,toy,h,project):
    p=toy.p;eq=[]
    for choice in product(*(groups[i] for i in toy.slots)):
        origin=(0,)*h;dirs=[]
        for G in choice:
            base=project(G[0]);origin=add(origin,base,p)
            dirs.extend(sub(project(v),base,p) for v in G[1:])
        a=normal(dirs,p,h);eq.append((a,dot(a,origin,p)))
    assert len(eq)<p
    ending=[[] for _ in range(h)]
    for a,b in eq: ending[max(i for i,x in enumerate(a) if x)].append((a,b))
    z=[]
    for i,E in enumerate(ending):
        bad={(b-sum(a[j]*z[j] for j in range(i)))*pow(a[i],-1,p)%p for a,b in E}
        z.append(next(j for j in range(len(bad)+1) if j not in bad))
    assert all(dot(z,a,p)!=b for a,b in eq)
    return tuple(z)

def compression():
    rng=Random(822)
    for p,m,n,w,t in [(101,2,11,2,3),(127,2,12,2,3),(257,2,16,2,4),(1009,3,9,2,2)]:
        toy=Toy(p,[rng.sample(range(1,p),n) for _ in range(m)],w)
        groups,k,C=partition_cover(toy,t); assert C<p
        a=toy.calibrated();U=[];initial=toy.n-1;steps=0
        while True:
            h,proj,lift=quotient(U,p,toy.n)
            if h<=max(1,k):break
            z=find_outside(groups,toy,h,proj);row=toy.row(z,lift)
            assert proj(row)!=(0,)*h
            U.append(row);assert len(rref(U,p,toy.n))==len(U)
            steps+=1
        h,_,_=quotient(U,p,toy.n)
        assert h-1<=max(0,k-1) and steps==initial-(h-1)
        assert len(rref([a,*U],p,toy.n))==1+len(U)
        print(f'Partition cover C{p}: D={toy.n}, {C} flats, dimension <= {k}; corank {initial}->{h-1}, {steps} verified increases.')

def projected_sum(toy,proj,dim):
    S={(0,)*dim}
    for i in toy.slots:
        C=[proj(a)[:dim] for a in toy.cols[i]]
        S={add(a,b,toy.p) for a in S for b in C}
    return S

def plane_reduction():
    rng=Random(73)
    for p,m,n in [(17,3,4),(31,3,5),(101,3,8)]:
        # Draw until the finite toy has full coverage. This is NOT the seed algorithm.
        while True:
            try:toy=Toy(p,[rng.sample(range(1,p),n) for _ in range(m)],1);break
            except AssertionError:pass
        assert n**m<p*p
        a=toy.calibrated();U=[];Mmax=0;steps=0
        while True:
            h,proj,lift=quotient(U,p,toy.n)
            if h==1:break
            S=projected_sum(toy,proj,2);Mmax=max(Mmax,len(S))
            z=next(z for z in product(range(p),repeat=2) if z not in S)
            z=z+(0,)*(h-2);row=toy.row(z,lift)
            assert proj(row)[:2]!=(0,0)
            U.append(row);steps+=1
        assert len(rref([a,*U],p,toy.n))==toy.n and steps==toy.n-1
        print(f'2D reduction C{p}: {steps} independent increases; maximum toy projected support {Mmax} < {p*p}. (Exhaustive toy avoidance.)')

def transport():
    rng=Random(41)
    p=17;m=3;n=5
    while True:
        try:
            toy=Toy(p,[rng.sample(range(1,p),n) for _ in range(m)],1)
            omitted=[toy.represent([j for j in range(m) if j!=i]) for i in range(m)]
            if all(len(R)==p for R in omitted):break
        except AssertionError:pass
    refs=(0,)*m;rows=[]
    def Z(g):return toy.inc(toy.rep[g%p])
    def v(g):return tuple(a-b for a,b in zip(Z(g),Z(0)))
    def contrast(X,Y):
        assert toy.val(X)==toy.val(Y)
        a=tuple(x-y for x,y in zip(toy.inc(X),toy.inc(Y)))
        assert dot(a,toy.lam,p)==0
        rows.append(a);return a
    def pinned(g,i,j):
        Y=list(omitted[i][(g-toy.B[i][j])%p]);Y.insert(i,j);return tuple(Y)
    def plus(g,h):
        """Construct actual short contrasts whose sum is v(g)+v(h)-v(g+h)."""
        g%=p;h%=p;gstart,hstart=g,h;X0=toy.rep[0];Xh=toy.rep[h]
        witness=[0]*toy.n
        for i,(a,b) in enumerate(zip(X0,Xh)):
            X=pinned(g,i,a);Y=pinned(h,i,b)
            Xp=list(X);Yp=list(Y);Xp[i]=b;Yp[i]=a;Xp=tuple(Xp);Yp=tuple(Yp)
            gn=toy.val(Xp);hn=toy.val(Yp)
            assert (gn+hn)%p==(g+h)%p
            R=[contrast(toy.rep[g],X),contrast(toy.rep[h],Y),contrast(Xp,toy.rep[gn]),contrast(Yp,toy.rep[hn])]
            for r in R:witness=[x+y for x,y in zip(witness,r)]
            g,h=gn,hn
        assert g==(gstart+hstart)%p and h==0
        assert tuple(witness)==tuple(x+y-z for x,y,z in zip(v(gstart),v(hstart),v(gstart+hstart)))
    # Compare actual singleton column change to v(group difference).
    for i,F in enumerate(toy.B):
        for j,A in enumerate(F[1:],1):
            X=list(refs);X[i]=j;X=tuple(X);delta=(A-F[0])%p
            contrast(X,toy.rep[(toy.g0+delta)%p]);contrast(refs,toy.rep[toy.g0]);plus(toy.g0,delta)
            # Binary addition chain for multiplication by p (known integer).
            y=delta
            for bit in bin(p)[3:]:
                plus(y,y);y=2*y%p
                if bit=='1':plus(y,delta);y=(y+delta)%p
            assert y==0
    H=hermite_normal_form(Matrix(rows).T)
    assert H.shape==(toy.n,toy.n)
    inv=H.inv()
    assert all(x.q==1 for x in inv*p)
    SNF=smith_normal_form(H,domain=ZZ)
    diag=[abs(int(SNF[j,j])) for j in range(toy.n)]
    assert all(d in (1,p) for d in diag) and p in diag
    print(f'Short-fiber transport C{p}: {len(rows)} genuine contrasts, {toy.n} columns; Smith invariants {diag}; exact inclusion {p}Z^D verified by HNF.')

if __name__=='__main__':
    compression();plane_reduction();transport()
