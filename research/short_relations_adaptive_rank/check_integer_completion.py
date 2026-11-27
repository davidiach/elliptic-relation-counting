"""Exact finite checks of integer completion; not a tensor-contractor implementation."""
from itertools import combinations, product
from math import gcd, lcm, prod
from random import Random
import sympy as sp


def rank(rows, p, n):
    B={}
    for row in rows:
        u=[int(x)%p for x in row]
        for j in sorted(B):
            if u[j]:u=[(x-u[j]*y)%p for x,y in zip(u,B[j])]
        j=next((j for j,x in enumerate(u) if x),None)
        if j is not None:B[j]=[x*pow(u[j],-1,p)%p for x in u]
    return len(B)


def gram(A,n):
    return sp.Matrix(n,n,lambda i,j:sum(a[i]*a[j] for a in A))


def mass(G,W):
    n=G.rows;k=len(W)
    if not k:return int(G.det())
    W=sp.Matrix(W)
    Z=G.row_join(W.T).col_join(W.row_join(sp.zeros(k)))
    return (-1)**k*int(Z.det())


def valuation(x,p):
    assert x!=0;v=0
    while x%p==0:v+=1;x//=p
    return v


def choose_real(A,W,n):
    """Use only full/prefix Grams to select a real-independent completion."""
    B=list(W);G=gram(A,n)
    while len(B)<n:
        ns=sp.Matrix(B).nullspace() if B else [sp.eye(n)[:,0]]
        v=ns[0];u=[int(x*lcm(*(int(y.q) for y in v))) for x in v]
        u=sp.Matrix(u);assert (u.T*G*u)[0]>0
        lo=0;hi=len(A)
        while hi-lo>1:
            mid=(lo+hi)//2
            if (u.T*gram(A[lo:mid],n)*u)[0]>0:hi=mid
            else:lo=mid
        assert sum(x*y for x,y in zip(A[lo],u))!=0
        B.append(A[lo])
    assert sp.Matrix(B).det()!=0
    return B


def choose_mass(A,W,n,p):
    """Preserve squared-minor mass modulo p^(v+1), via suffix queries."""
    B=list(W);lo=0;C=mass(gram(A,n),W);assert C>0
    v=valuation(C,p);mod=p**(v+1)
    while len(B)<n:
        L=lo;R=len(A)
        assert mass(gram(A[L:],n),B)%mod
        while R-L>1:
            mid=(L+R)//2
            q=(mass(gram(A[L:],n),B)-mass(gram(A[mid:],n),B))%mod
            if q:R=mid
            else:L=mid
        B.append(A[L]);lo=L+1
        assert mass(gram(A[lo:],n),B)%mod
    d=int(sp.Matrix(B).det());assert d and 2*valuation(d,p)<=v
    return B,v


def bounds(A,W,p,n):
    G=gram(A,n);s=max(sum(x!=0 for x in a) for a in A);L=s*max(int(G[i,i]) for i in range(n))
    k=len(W);d=n-k;delta=int((sp.Matrix(W)*sp.Matrix(W).T).det()) if k else 1
    C=mass(G,W);assert 0<C<=delta*L**d
    B=choose_real(A,W,n);det=int(sp.Matrix(B).det())
    assert det*det<=delta*s**d
    rnk=rank(B,p,n);assert n-rnk<=valuation(det,p)
    return C,delta,s,L,rnk


def fixture():
    p=13;lam=[3,4,5,8,6,7];n=len(lam)
    A=[[0,-1,1,0,1,1],[-1,0,1,0,1,-1],[1,0,1,0,1,0],[1,-1,0,1,-1,0],[0,0,-1,0,1,0],[0,1,-1,1,0,1],[0,0,0,0,-1,1],[-1,0,0,-1,1,-1],[-1,-1,0,1,0,0],[-1,1,0,0,0,0],[0,-1,1,0,0,0],[0,0,-1,-1,-1,1]]
    W=[A[i] for i in (3,8,9,10,11)];u=sp.Matrix([8,8,8,3,3,1]);G=gram(A,n)
    assert len({tuple(a) for a in A})==len(A)<p and max(G.diagonal())<p
    assert all(sum(x*y for x,y in zip(a,lam))%p==1 for a in A)
    assert rank(A,p,n)==n and rank(W,p,n)==n-1
    assert all(x%p==0 for x in sp.Matrix(W)*u) and all(x%p==0 for x in G*u)
    assert int((u.T*G*u)[0])%(p*p)==0
    C,delta,s,L,rnk=bounds(A,W,p,n)
    assert (C,delta,s,rnk)==(169,36,4,6) and delta*s<p*p
    # Every row is a valid signed tuple difference on E/F7 of order 13.
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;v,w=Q
        if x==v and (y+w)%7==0:return None
        t=((3*x*x)*pow(2*y,-1,7) if P==Q else (w-y)*pow(v-x,-1,7))%7
        X=(t*t-x-v)%7;return X,(t*(x-X)-y)%7
    pts=[(x,y) for x,y in product(range(7),repeat=2) if (y*y-x*x*x-3)%7==0]
    P=pts[0];cycle=[None]
    for _ in range(12):cycle.append(add(cycle[-1],P))
    assert len(pts)==12 and set(cycle)=={None,*pts} and add(cycle[-1],P) is None
    for a in A:
        X=Y=None
        for c,l in zip(a,lam):
            X=add(X,cycle[1+l if c==-1 else 1]);Y=add(Y,cycle[1+l if c==1 else 1])
        assert Y==add(X,P)
    print(f'Order-13 curve subset: {len(A)} genuine fixed-shift rows, rank {n}/{n}, initial rank {len(W)}. Diagonals {list(G.diagonal())}; quadratic mass {int((u.T*G*u)[0])}; completion mass {C}; squared W volume {delta}. Integer completion rank {rnk}/{n}.')


def random_checks():
    rng=Random(812);checks=0;mass_checks=0;direct=0;rescued=0
    for p,n in [(5,3),(7,4),(13,5),(17,6)]:
        for _ in range(12):
            A=[tuple(rng.choice([-1,0,0,1]) for _ in range(n)) for j in range(n+3)]
            if rank(A,p,n)!=n:continue
            W=[]
            for a in A:
                if len(W)>=n//2:break
                if rank([*W,a],p,n)>len(W):W.append(a)
            C,delta,s,L,rnk=bounds(A,W,p,n);checks+=1
            if n<=4:
                q=sum(int(sp.Matrix([*W,*[A[j] for j in ids]]).det())**2 for ids in combinations(range(len(A)),n-len(W)))
                assert q==C;direct+=1
            B,v=choose_mass(A,W,n,p);mass_checks+=1
            assert n-rank(B,p,n)<=v//2
            if C%(p*p)==0 and rnk==n:rescued+=1
    print(f'Exact sparse arrays: {checks} volume bounds and real completions; {direct} direct Cauchy-Binet checks; {mass_checks} modular suffix completions. {rescued} real completions full rank despite mass divisible by p^2.')


def group_checks():
    rng=Random(61);checks=0
    for p,m,k in [(7,2,3),(11,2,4),(13,3,3),(19,3,3)]:
        for _ in range(12):
            Bs=[rng.sample(range(1,p),k) for i in range(m)];n=m*(k-1)
            Ts=[[] for _ in range(p)]
            for X in product(range(k),repeat=m):Ts[sum(Bs[i][a] for i,a in enumerate(X))%p].append(X)
            if not all(Ts):continue
            A=[]
            for g in range(p):
                for X in Ts[g]:
                    for Y in Ts[(g+1)%p]:
                        A.append([int(Y[i]==j)-int(X[i]==j) for i in range(m) for j in range(1,k)])
            assert rank(A,p,n)==n
            B=choose_real(A,[],n);det=int(sp.Matrix(B).det());s=max(sum(x*x for x in a) for a in A)
            assert det*det<=s**n and n-rank(B,p,n)<=valuation(det,p)
            checks+=1
    print(f'Complete fixed-shift cyclic fixtures: {checks} real-basis and modular rank bounds checked.')


def bad_real_basis():
    p=13;Bs=[[8,2,3],[6,12,3],[8,5,7]];n=6
    Ts=[[] for _ in range(p)]
    for X in product(range(3),repeat=3):Ts[sum(Bs[i][a] for i,a in enumerate(X))%p].append(X)
    assert all(Ts)
    A=[]
    for g in range(p):
        for X in Ts[g]:
            for Y in Ts[(g+1)%p]:
                A.append([int(Y[i]==j)-int(X[i]==j) for i in range(3) for j in (1,2)])
    B=[[-1,1,0,-1,1,0],[0,-1,1,0,-1,0],[1,-1,0,-1,0,1],[0,1,1,0,0,0],[0,1,0,-1,-1,0],[0,0,0,-1,1,-1]]
    assert all(a in A for a in B) and rank(A,p,n)==n
    assert abs(int(sp.Matrix(B).det()))==p and rank(B,p,n)==n-1
    C=choose_real(B+[a for a in A if a not in B],[],n)
    assert abs(int(sp.Matrix(C).det()))==p
    print(f'Complete C13 obstruction: {len(A)} fixed-shift rows, full rank 6; a possible rational-basis ordering returns determinant +/-13 and modular rank 5. This fixture has wrap and is not a canonical-order counterexample.')


def doubled_checks():
    # Complete doubled-tuple relations, same ordinary reduced incidence matrix.
    rng=Random(152);cases=0
    for p,m,k in [(7,2,2),(11,2,3),(13,3,2),(19,3,2)]:
        for _ in range(5):
            Bs=[rng.sample(range(1,p),k) for i in range(m)];n=m*(k-1)
            Ts=[[] for _ in range(p)];loc=[i for i in range(m) for _ in (0,1)]
            for X in product(range(k),repeat=2*m):
                g=sum(Bs[i][j] for i,j in zip(loc,X))%p
                v=[sum(int(j==a) for h,j in zip(loc,X) if h==i) for i in range(m) for a in range(1,k)]
                Ts[g].append(v)
            if not all(Ts):continue
            A=[[y-x for x,y in zip(X,Y)] for g in range(p) for X in Ts[g] for Y in Ts[(g+1)%p]]
            assert rank(A,p,n)==n and max(sum(x*x for x in a) for a in A)<=8*m
            B=choose_real(A,[],n);det=int(sp.Matrix(B).det());assert det*det<=(8*m)**n
            cases+=1
    print(f'Doubled-tuple fixtures: {cases} complete-group relation matrices verified, squared row norm <=8m and integer determinant bound passed.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    fixture();random_checks();group_checks();bad_real_basis();doubled_checks()
