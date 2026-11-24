"""Finite regression for the proposed auxiliary polynomial-ring arithmetic.

This is NOT an ECDLP implementation or an end-to-end contraction benchmark.
The asymptotic proof uses exact SVP; these finite fixtures use LLL and check
its result. Requires Python 3.10+ and sympy.
"""
from math import gcd
from fractions import Fraction
from random import Random
import sympy as sp


def fold(a, b, lam):
    n = len(a)
    c = [0] * n
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            k = i+j
            c[k % n] += x*y*(lam if k >= n else 1)
    return c


def ev(a, gamma, p):
    v=0
    for x in reversed(a): v=(v*gamma+x)%p
    return v


def add(*args): return [sum(t) for t in zip(*args)]
def scale(k,a): return [k*x for x in a]
def round_rat(x): return int(sp.floor(x+sp.Rational(1,2)))


def lll_exact(rows):
    # Small exact rational LLL for fixtures; not used for asymptotic SVP.
    B=[list(map(int,row)) for row in rows];n=len(B)
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    def gs():
        ortho=[];norm=[];mu=[[Fraction(0)]*n for _ in B]
        for i,b in enumerate(B):
            v=list(map(Fraction,b))
            for j,w in enumerate(ortho):
                t=dot(b,w)/norm[j];mu[i][j]=t
                v=[x-t*y for x,y in zip(v,w)]
            ortho.append(v);norm.append(dot(v,v))
        return mu,norm
    mu,norm=gs();k=1
    while k<n:
        for j in range(k-1,-1,-1):
            q=(mu[k][j]+Fraction(1,2)).__floor__()
            if q:
                B[k]=[a-q*b for a,b in zip(B[k],B[j])]
                for l in range(j):mu[k][l]-=q*mu[j][l]
                mu[k][j]-=q
        if norm[k]>=(Fraction(3,4)-mu[k][k-1]**2)*norm[k-1]:
            k+=1
        else:
            B[k],B[k-1]=B[k-1],B[k]
            mu,norm=gs();k=max(1,k-1)
    return B


def fixture_prime(bits, both):
    p=int(sp.nextprime(2**(bits-1)+314159265358979323846))
    while ((p-1)%15==0 if both else gcd(p-1,15)==1) is False:
        p=int(sp.nextprime(p))
    return p


def construct(p, n=15):
    X=sp.Symbol('X'); qs=(3,5); both=(p-1)%n==0
    if both:
        lam=1; d=(qs[0]-1)*(qs[1]-1)
        for a in range(2,10000):
            gamma=pow(a,(p-1)//n,p)
            if all(pow(gamma,n//q,p)!=1 for q in qs): break
        else: raise AssertionError('root not found')
        F=sp.Poly(sp.cyclotomic_poly(n,X),X)
    else:
        assert gcd(p-1,n)==1
        lam=2; d=n; gamma=pow(lam,pow(n,-1,p-1),p)
        F=None
    assert pow(gamma,n,p)==lam
    basis=sp.zeros(d,d);basis[0,0]=p
    for i in range(1,d): basis[i,0]=-pow(gamma,i,p);basis[i,i]=1
    # This is only a finite fixture: verify the LLL vector is small enough.
    f=lll_exact(basis.tolist())[0]+[0]*(n-d)
    S=d*2**((p.bit_length()+d-1)//d)
    assert max(map(abs,f))<=S and ev(f,gamma,p)==0
    def shift(a,j):
        e=[0]*n;e[j]=1
        return fold(a,e,lam)
    cols=[shift(f,j) for j in range(d)]
    if both:
        fv=[int(F.nth(i)) for i in range(n)]
        cols += [shift(fv,j) for j in range(n-d)]
    H=sp.Matrix.hstack(*map(sp.Matrix,cols))
    assert H.det()!=0 and all(ev(c,gamma,p)==0 for c in cols)
    # Approximate a large scalar by a lattice vector.
    actual=max(abs(x) for c in cols for x in c)
    for multiplier in range(1,16):
        C=8*n*n*actual*multiplier
        target=sp.Matrix([C]+[0]*(n-1))
        M=list(map(int,H*sp.Matrix([round_rat(t) for t in H.inv()*target])))
        assert ev(M,gamma,p)==0
        assert M[0]>2*sum(map(abs,M[1:]))
        A=sp.Matrix.hstack(*[sp.Matrix(shift(M,j)) for j in range(n)])
        N=abs(int(A.det()));assert N%p==0
        if gcd(N//p,p)==1:break
    else:raise AssertionError('could not find simple p factor')
    c=N//p; e=c*pow(c,-1,p)
    assert (e*e-e)%N==0 and e%p==1
    Ai=A.inv()
    def lift(x):
        v=sp.Matrix([e*(x%p)]+[0]*(n-1))
        t=Ai*v
        r=list(map(int,v-A*sp.Matrix([round_rat(z) for z in t])))
        assert ev(r,gamma,p)==x%p
        assert 2*max(map(abs,r))<=M[0]+2*sum(map(abs,M[1:]))
        return r
    def quotient(v):
        t=Ai*sp.Matrix(v)
        if not all(z.q==1 for z in t):return None
        q=list(map(int,t));assert fold(M,q,lam)==v
        return q
    return gamma,lam,M,N,lift,quotient


def add_ec(P,Q,p,a):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if (y+v)%p==0:return None
        s=(3*x*x+a)*pow(2*y,-1,p)%p
    else:s=(v-y)*pow(u-x,-1,p)%p
    w=(s*s-x-u)%p
    return w,(s*(x-w)-y)%p


def chord_residuals(u,v,w,s,y,a,even):
    z=w if even else u
    return [s*s-u-v-w,
            a+(2 if even else -2)*s*y+2*z*(u+v+w)-(u*v+u*w+v*w)]


def check_chords():
    total=double=0
    for p,a,c in [(13,4,1),(17,2,2),(19,1,0),(37,3,1)]:
        points=[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-a*x-c)%p==0]
        # Exhaustively check all line intersections, including repeated roots.
        for s in range(p):
            for d in range(p):
                xs=[x for x in range(p) if ((s*x+d)**2-x**3-a*x-c)%p==0]
                for u in xs:
                    for v in xs:
                        w=(s*s-u-v)%p
                        if w not in xs:continue
                        P=(u,(s*u+d)%p);Q=(v,(s*v+d)%p);R=(w,(-s*w-d)%p)
                        # Repeated roots must actually have multiplicity: check X coefficient.
                        if (a-2*s*d-u*v-u*w-v*w)%p:continue
                        for even in (False,True):
                            assert not any(t%p for t in chord_residuals(u,v,w,s,R[1] if even else P[1],a,even))
                        assert add_ec(P,Q,p,a)==R
                        total+=1;double+=u==v
    print(f'Chord identities: {total} verified line triples, including {double} doubling triples.')


def transform_checks(p,gamma,lam,M,lift,quotient):
    n=len(M);rng=Random(p)
    # Random quadratic equations with known zero/nonzero scalar residuals.
    divisible=nonzero=0
    for _ in range(100):
        x,y,z=[rng.randrange(p) for _ in range(3)]
        t=(x*y+z)%p
        for delta in (0,1):
            A,B,C,T=map(lift,(x,y,z,(t+delta)%p))
            residual=add(fold(A,B,lam),C,scale(-1,T))
            q=quotient(residual)
            assert (q is not None)==(delta==0)
            assert ev(residual,gamma,p)==(-delta)%p
            if q is not None: divisible+=1
            else:nonzero+=1
    print(f'  Polynomial quotient checks: {divisible} exact, {nonzero} rejected.')
    # Exact twisted DFT over an unramified quadratic extension, with conjugacy.
    # Find q0 == -1 mod 8n; lift primitive n-th root by exponentiation.
    q0=8*n-1
    while not sp.isprime(q0):q0+=8*n
    # A quadratic, not cubic, coefficient modulus suffices.
    capM=max(map(abs,M))
    capA=(M[0]+2*sum(map(abs,M[1:]))+1)//2
    gap=M[0]-2*sum(map(abs,M[1:]))
    capD=2*n*capA**2+2*capA
    capQ=(capD+gap-1)//gap
    bound=capD+2*n*capM*capQ
    mod=q0
    while mod<=2*bound:mod*=q0
    def E(a):return (a%mod,0)
    def ea(a,b):return ((a[0]+b[0])%mod,(a[1]+b[1])%mod)
    def em(a,b):return ((a[0]*b[0]-a[1]*b[1])%mod,(a[0]*b[1]+a[1]*b[0])%mod)
    def ep(a,k):
        r=E(1)
        while k:
            if k&1:r=em(r,a)
            a=em(a,a);k//=2
        return r
    def conj(a):return (a[0],-a[1]%mod)
    # Unit group order in the extension: (q0^2-1)*q0^(2(h-1));
    # take Teichmuller n-torsion via raising to (q0^2-1)/n * q0^(h-1).
    h=0;v=mod
    while v>1:v//=q0;h+=1
    exponent=((q0*q0-1)//n)*(q0**(h-1))
    for _ in range(1000):
        w=ep((rng.randrange(mod),rng.randrange(mod)),exponent)
        if ep(w,n)==E(1) and ep(w,n//3)!=E(1) and ep(w,n//5)!=E(1):break
    else:raise AssertionError('root search failed')
    assert conj(w)==ep(w,n-1)
    beta=pow(lam,pow(n,-1,(mod//q0)*(q0-1)),mod)
    assert pow(beta,n,mod)==lam
    def dft(a):
        out=[]
        for j in range(n):
            t=E(0)
            for i,x in enumerate(a):t=ea(t,em(E(x*pow(beta,i,mod)),ep(w,(i*j)%n)))
            out.append(t)
        assert all(out[-j%n]==conj(out[j]) for j in range(n))
        return out
    def row_dft(a):
        L,V=3,5
        U=[]
        for j in range(V):
            row=[]
            for s in range(L):
                v=E(0)
                for h in range(L):
                    i=j+V*h
                    v=ea(v,em(E(a[i]*pow(beta,i,mod)),ep(w,V*h*s)))
                row.append(v)
            assert row[2]==conj(row[1]) and row[0][1]==0
            U.append(row)
        return U
    def full_from_rows(U):
        out=[]
        for f in range(n):
            s=f%3
            v=E(0)
            for j in range(5):v=ea(v,em(U[j][s],ep(w,j*f)))
            out.append(v)
        return out
    for _ in range(8):
        a,b=lift(rng.randrange(p)),lift(rng.randrange(p))
        assert full_from_rows(row_dft(a))==dft(a)
        assert full_from_rows(row_dft(b))==dft(b)
        assert dft(fold(a,b,lam))==[em(x,y) for x,y in zip(dft(a),dft(b))]
        t=lift(ev(a,gamma,p)*ev(b,gamma,p)%p)
        residual=add(fold(a,b,lam),scale(-1,t));Q=quotient(residual);assert Q is not None
        assert max(map(abs,Q))<=capQ
        # Middle processors need only the input transforms: no output array.
        F,G,T,H,P=[full_from_rows(row_dft(v)) for v in (a,b,t,M,Q)]
        assert all(ea(em(f,g),E(0))==ea(tt,em(h,pp))
                   for f,g,tt,h,pp in zip(F,G,T,H,P))
        assert [em(x,y) for x,y in zip(dft(M),dft(Q))]==dft(residual)
    print(f'  Twisted Fourier: 8 product and exact quotient identities; two-stage row factorization and conjugate compression verified.')


def main():
    for bits,both in [(256,False),(384,True)]:
        p=fixture_prime(bits,both);g,lam,M,N,lift,quo=construct(p)
        print(f'{p.bit_length()}-bit prime, n=15, lambda={lam}, det bit length={N.bit_length()}, v_p(det)=1.')
        transform_checks(p,g,lam,M,lift,quo)
    check_chords()

if __name__=='__main__': main()
