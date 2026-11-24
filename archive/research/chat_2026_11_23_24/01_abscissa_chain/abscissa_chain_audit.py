"""Toy regression of quartic chain soundness and torsion-pinning enumeration.

Requires Python 3 and sympy. Imports the previous polynomial_ring_verifier.py
for arithmetic component fixtures. Not an end-to-end contraction benchmark.
"""
from itertools import product
from random import Random
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from polynomial_ring_verifier import construct, fixture_prime, fold, ev, add, scale, add_ec, round_rat


def S(u,v,w,a,c):
    return (u-v)**2*w*w - 2*((u+v)*(u*v+a)+2*c)*w + (u*v-a)**2 - 4*c*(u+v)


def check_chain(p,a,c,m):
    points = [(x,y) for x,y in product(range(p),repeat=2) if (y*y-x*x*x-a*x-c)%p==0]
    torsion = {x for x,y in points if y==0}
    targets = {x for x,y in points if y!=0}
    # Independent signed group addition establishes the allowed abscissa tuples.
    truth={t:set() for t in targets}
    def actual(prefix, state):
        if len(prefix)==m:
            if state is not None and state[1] and state[0] in targets:
                truth[state[0]].add(tuple(prefix))
            return
        for Q in points:
            R=add_ec(state,Q,p,a)
            if R is None: continue
            if 1<len(prefix)+1<m and R[1]==0: continue
            actual(prefix+[Q[0]],R)
    actual([],None)
    pred={t:set() for t in targets}
    unrestricted={t:set() for t in targets}
    # All field abscissas, including non-liftable values, are allowed in the candidate chain.
    prev={w:[(u,v) for u,v in product(range(p),repeat=2) if S(u,v,w,a,c)%p==0] for w in range(p)}
    def reverse(w,k,guard):
        if k==1: return {(w,)}
        out=set()
        for u,v in prev[w]:
            if guard and k>2 and u in torsion: continue
            out.update(t+(v,) for t in reverse(u,k-1,guard))
        return out
    for t in targets:
        pred[t]=reverse(t,m,True)
        unrestricted[t]=reverse(t,m,False)
        assert pred[t]==truth[t], (p,a,c,m,t,pred[t]^truth[t])
    hits=sum(map(len,truth.values()))
    bad=sum(len(unrestricted[t]-truth[t]) for t in targets)
    nonlift={x for x in range(p) if not any(u==x for u,y in points)}
    false_nonlift=sum(any(x in nonlift for x in tup) for t in targets for tup in unrestricted[t]-truth[t])
    print(f'p={p:2}, curve a={a:2}, c={c:2}, m={m}, targets={len(targets):2}, guarded supports={hits:6}, unguarded extra={bad:6}, containing non-liftable leaf={false_nonlift:6}',flush=True)
    return hits,bad,false_nonlift


def check_small_cosets():
    # n=3, X^3=2, gamma=7 mod 11, M=4+X, determinant 66.
    p=11; gamma=7; lam=2; M=[4,1,0]; n=3; cap=3
    A=sp.Matrix.hstack(*[sp.Matrix(fold(M,[int(i==j) for i in range(n)],lam)) for j in range(n)])
    Ainv=A.inv()
    B=sp.Matrix([[p,-gamma,-gamma**2%p],[0,1,0],[0,0,1]])
    # fix sign of third basis entry: should represent -gamma^2 modulo p
    assert all(ev(list(B[:,j]),gamma,p)==0 for j in range(n))
    C=B.inv()*A
    assert all(x.q==1 for x in C)
    H=hermite_normal_form(C.applyfunc(int))
    diagonal=[int(H[i,i]) for i in range(n)]
    assert abs(int(A.det()))//p==sp.prod(diagonal)==6
    checked=0
    for x in range(p):
        got=set()
        for z in product(*[range(h) for h in diagonal]):
            v=sp.Matrix([x,0,0])+B*sp.Matrix(z)
            r=v-A*sp.Matrix([round_rat(t) for t in Ainv*v])
            # |q_i| <= ||A^-1||_infinity (cap + ||r||_infinity).
            Rmax=max(abs(t) for t in r)
            norm=max(sum(abs(Ainv[i,j]) for j in range(n)) for i in range(n))
            qmax=int(sp.floor(norm*(cap+Rmax)))
            for q in product(range(-qmax,qmax+1),repeat=n):
                v2=tuple(map(int,r+A*sp.Matrix(q)))
                if max(map(abs,v2))<=cap:
                    assert v2 not in got and ev(v2,gamma,p)==x
                    got.add(v2)
        want={v for v in product(range(-cap,cap+1),repeat=n) if ev(v,gamma,p)==x}
        assert got==want,(x,got^want)
        checked+=len(got)
    print(f'Ideal-coset enumeration: {checked} bounded representatives, all {p} field values, {sp.prod(diagonal)} cosets per value; exact, disjoint, exhaustive.')


def quartic_poly(U,V,W,A,C,lam):
    def mul(*xs):
        x=xs[0]
        for y in xs[1:]: x=fold(x,y,lam)
        return x
    d=add(U,scale(-1,V));z=add(mul(U,V),scale(-1,A))
    return add(mul(d,d,W,W),scale(-2,mul(add(mul(add(U,V),add(mul(U,V),A)),scale(2,C)),W)),mul(z,z),scale(-4,mul(C,add(U,V))))


def check_quartic_lifts():
    for bits,both in [(256,False),(384,True)]:
        p=fixture_prime(bits,both); g,lam,M,N,lift,quo=construct(p)
        rng=Random(p+11); pos=neg=0
        # Choose u,v,w,a and solve S3=0 for c (linear unless u+v+w=0).
        for _ in range(40):
            u,v,w,a=[rng.randrange(p) for _ in range(4)]
            denom=4*(u+v+w)%p
            if denom==0: continue
            c=S(u,v,w,a,0)*pow(denom,-1,p)%p
            for c0 in (c,(c+1)%p):
                residual=quartic_poly(*map(lift,(u,v,w,a,c0)),lam)
                assert ev(residual,g,p)==S(u,v,w,a,c0)%p
                q=quo(residual)
                assert (q is not None)==(c==c0)
                if q is None: neg+=1
                else: pos+=1; assert fold(M,q,lam)==residual
        print(f'{bits}-bit polynomial ring lambda={lam}: {pos} divisible quartic lifts and {neg} nonzero residuals rejected.',flush=True)

if __name__=='__main__':
    totals=[0,0,0]
    for p,a,c in [(7,0,3),(11,1,0),(13,4,1),(17,2,0),(19,1,0)]:
        assert (4*a**3+27*c*c)%p
        for m in (2,3,4):
            totals=[x+y for x,y in zip(totals,check_chain(p,a,c,m))]
    print(f'Aggregate guarded support comparisons: {totals[0]}; unguarded extra supports: {totals[1]}; non-liftable extras: {totals[2]}.')
    check_small_cosets()
    check_quartic_lifts()
