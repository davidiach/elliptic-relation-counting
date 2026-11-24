"""Finite checks for the short multiplier, shifted masks and mixed-twist bound.
Not an end-to-end contraction or cryptographic-scale timing experiment.
Requires sympy and the previous polynomial_ring_verifier.py.
"""
import itertools, math, random
from collections import Counter
from fractions import Fraction
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from polynomial_ring_verifier import (construct, fixture_prime, fold, ev, add, scale, lll_exact, add_ec)


def toy_construct(p,n=5):
    both=(p-1)%n==0;lam=1 if both else 2;d=n-1 if both else n
    if both:
        gamma=next(g for a in range(2,100) if (g:=pow(a,(p-1)//n,p))!=1)
    else:gamma=pow(2,pow(n,-1,p-1),p)
    B=sp.zeros(d,d);B[0,0]=p
    for i in range(1,d):B[i,0]=-pow(gamma,i,p);B[i,i]=1
    f=lll_exact(B.tolist())[0]+[0]*(n-d)
    def shift(a,j):
        v=[0]*n;v[j]=1;return fold(a,v,lam)
    cols=[shift(f,i) for i in range(d)]
    if both:cols.append([1]*n)
    H=sp.Matrix.hstack(*map(sp.Matrix,cols));S=max(abs(v) for col in cols for v in col)
    C=8*n*n*S
    v=sp.Matrix([C]+[0]*(n-1))
    M=list(map(int,H*sp.Matrix([sp.floor(t+sp.Rational(1,2)) for t in H.inv()*v])))
    A=sp.Matrix.hstack(*[sp.Matrix(shift(M,i)) for i in range(n)])
    D=abs(int(A.det()));assert D%p==0 and D<p*p
    e=(D//p)*pow(D//p,-1,p);Ai=A.inv()
    def lift(x):
        v=sp.Matrix([e*(x%p)]+[0]*(n-1))
        return list(map(int,v-A*sp.Matrix([sp.floor(z+sp.Rational(1,2)) for z in Ai*v])))
    def quotient(v):
        t=Ai*sp.Matrix(v)
        return list(map(int,t)) if all(z.q==1 for z in t) else None
    return gamma,lam,M,D,lift,quotient


def projector_checks():
    rng=random.Random(20261123)
    for bits,both in ((96,False),(128,True)):
        p=int(sp.nextprime(2**(bits-1)+314159))
        while ((p-1)%5==0)!=both:p=int(sp.nextprime(p))
        gamma,lam,M,D,lift,quotient=toy_construct(p)
        n=len(M); c=D//p; e=c*pow(c,-1,p)
        def shift(a,j):
            v=[0]*n;v[j]=1
            return fold(a,v,lam)
        A=sp.Matrix.hstack(*[sp.Matrix(shift(M,j)) for j in range(n)])
        Ai=A.inv()
        J=hermite_normal_form(A.row_join(e*sp.eye(n)))
        assert abs(int(J.det()))==c
        E=lll_exact(J.T.tolist())[0]
        assert quotient(scale(p,E)) is not None
        assert ev(E,gamma,p)!=0
        assert max(map(abs,E)) < M[0]-2*sum(map(abs,M[1:]))
        X=[-gamma]+[1]+[0]*(n-2)
        def arbitrary(x):
            v=lift(x)
            v=add(v,scale(rng.randrange(-4,5),X))
            a=sp.Matrix(v); t=Ai*a
            a=a-A*sp.Matrix([sp.floor(z+sp.Rational(1,2)) for z in t])
            v=list(map(int,a)); assert ev(v,gamma,p)==x%p
            return v
        good=bad=unprojected_failure=0
        for _ in range(60):
            x,y,z=[rng.randrange(p) for i in range(3)]
            for delta in (0,1):
                t=(x*x*y*z+delta)%p
                U,V,W,T=map(arbitrary,(x,y,z,t))
                F=add(fold(fold(fold(U,U,lam),V,lam),W,lam),scale(-1,T))
                raw=quotient(F)
                Q=quotient(fold(E,F,lam))
                assert (Q is not None)==(delta==0)
                if delta: bad+=1
                else:
                    good+=1;unprojected_failure+=raw is None
        print(f'{bits}-bit p, n={n}, lambda={lam}: bit length c={c.bit_length()}, '
              f'bit length max |E_i|={max(map(abs,E)).bit_length()}, '
              f'bit length M0={M[0].bit_length()}, valid {good}, rejected {bad}, '
              f'valid but failed without E {unprojected_failure}')
        assert unprojected_failure>0


def pointsets(p,a,c):
    # Work in F_p on the curve and on a quadratic twist d*y^2=f(x).
    # Use y'=d*y and x'=d*x: the twist is y'^2=x'^3+d^2*a*x'+d^3*c.
    d=next(d for d in range(2,p) if pow(d,(p-1)//2,p)==p-1)
    G=[None]+[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-a*x-c)%p==0]
    Gm=[None]+[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-d*d*a*x-d**3*c)%p==0]
    return d,G,Gm


def energy_checks():
    # Exact distribution over all independent translations of a small support.
    # The energy identity is checked by direct full target counts.
    for p,a,c,S in [(13,4,1,{0,1,4}),(17,2,2,{0,2,5}),(19,1,0,{0,1,3})]:
        d,G,Gm=pointsets(p,a,c); N=len(G);m=3; rho=Fraction(len(S),p)
        bases=[[P for P in G[1:] if (P[0]-shift)%p in S] for shift in range(p)]
        A={D:Fraction(sum(sum(P in base and add_ec(P,D,p,a) in base for P in G) for base in bases),p) for D in G}
        diff={(x-y)%p for x in S for y in S}
        assert all(A[D]<=4*rho*len(diff) for D in G if D is not None)
        assert sum(A.values())==Fraction(sum(len(F)**2 for F in bases),p)
        def neg(P): return None if P is None else (P[0],-P[1]%p)
        energy=sum(A[D]*A[E]*A[neg(add_ec(D,E,p,a))] for D,E in itertools.product(G,repeat=2))/N
        squares=0; positive=0
        for slots in itertools.product(bases,repeat=m):
            counts=Counter(add_ec(add_ec(P,Q,p,a),R,p,a) for P,Q,R in itertools.product(*slots))
            squares+=sum(v*v for v in counts.values());positive+=len(counts)
        direct=Fraction(squares,p**m*N); assert direct==energy
        bound=(Fraction(N-1)*rho)**m/N+Fraction(m,N)*4*rho*len(diff)*(4*len(S)**2)**(m-1)
        assert energy<=bound
        # Verify the exact group-theoretic mixed bound via two group sum lists.
        # For each proper assignment of twist slots and each rational target,
        # sum of the anti-rational part must be one of the shared 2-torsion points.
        shared_x={x for x,y in G[1:] if y==0}
        shared=[None]+[(x,0) for x in shared_x]
        anti=[P for P in Gm[1:] if P[1]!=0]
        def sumpts(xs,aa):
            r=None
            for x in xs:r=add_ec(r,x,p,aa)
            return r
        Nm=len(Gm); checks=0
        for T in G[1:]:
            if T[1]==0:continue
            total=0
            for k in (1,2):
                num=0
                for U in shared:
                    Um=None if U is None else (d*U[0]%p,0)
                    R=add_ec(T,U,p,a) # U=-U
                    rc=sum(sumpts(seq,a)==R for seq in itertools.product(G[1:],repeat=k))
                    ac=sum(sumpts(seq,d*d*a%p)==Um for seq in itertools.product(anti,repeat=m-k))
                    num+=rc*ac
                assert num<=len(shared)*N**(k-1)*Nm**(m-k-1)
                total+=math.comb(m,k)*num
            assert total<=4*(2**m-2)*max(N,Nm)**(m-2);checks+=1
        print(f'p={p}: energy {float(energy):.6f}, exact average support {positive/(p**m*N):.6f}, mixed-count target checks {checks}')


def partition_checks():
    # Toy surjective digit map, sorted in owner order; all partitions checked.
    p=31; gamma=3; n=4; width=2
    weights=[(pow(gamma,i,p)*2**j)%p for i in range(n) for j in range(width)]
    def supports(weights):
        s={0}
        for w in weights:s|={(x+w)%p for x in s}
        return s
    assert len(supports(weights))==p
    for m in range(2,len(weights)+1):
        chunks=[weights[j*len(weights)//m:(j+1)*len(weights)//m] for j in range(m)]
        counts=[len(supports(x)) for x in chunks]
        assert math.prod(counts)>=p and max(counts)**m>=p
    print('Entropy partition: all 7 partitions of a surjective 8-bit digit map pass.')

if __name__=='__main__':
    projector_checks(); energy_checks(); partition_checks()
