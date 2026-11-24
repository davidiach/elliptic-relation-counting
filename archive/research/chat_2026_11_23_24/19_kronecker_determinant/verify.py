#!/usr/bin/env python3
"""Exact finite regressions; not a benchmark of fast polynomial multiplication."""
from itertools import product
from math import prod, isqrt
from random import Random
from functools import lru_cache
from sympy import Matrix, ZZ, nextprime
from sympy.matrices.normalforms import smith_normal_form
rng=Random(20261124)

# 1. Block norm, using componentwise arithmetic in a split finite F_p algebra.
def block_tests():
    checks=zero_div=tail=0
    for _ in range(600):
        p=rng.choice([7,11,19,31,101]); r=rng.randrange(1,5)
        d=rng.randrange(0,min(8,p)); w=rng.randrange(1,p+1)
        a,t=rng.randrange(p),rng.randrange(1,p)
        zero=(0,)*r; one=(1,)*r
        add=lambda x,y:tuple((u+v)%p for u,v in zip(x,y))
        mul=lambda x,y:tuple(u*v%p for u,v in zip(x,y))
        sc=lambda x:(x%p,)*r
        def pm(f,g):
            out=[zero]*(len(f)+len(g)-1)
            for i,x in enumerate(f):
                for j,y in enumerate(g):out[i+j]=add(out[i+j],mul(x,y))
            return out
        def pe(f,x):
            z=zero
            for c in f[::-1]:z=add(mul(z,sc(x)),c)
            return z
        F=[tuple(rng.randrange(p) for _ in range(r)) for j in range(d+1)]
        if rng.randrange(4)==0: # force some zero divisors at a slot
            at=pe(F,a);F[0]=((F[0][0]-at[0])%p,)+F[0][1:]
        H=max(1,isqrt(max(1,w//max(1,d))))
        if rng.randrange(3)==0:H=rng.randrange(1,w+1)
        B=[one]
        for j in range(H):
            # Compose F with the affine linear expression a+t(X+j), by Horner.
            g=[zero]
            for c in F[::-1]:g=pm(g,[sc(a+t*j),sc(t)]);g[0]=add(g[0],c)
            B=pm(B,g)
        actual=one; expected=one
        for q in range(w//H):actual=mul(actual,pe(B,q*H))
        for j in range((w//H)*H,w):actual=mul(actual,pe(F,a+t*j))
        for j in range(w):expected=mul(expected,pe(F,a+t*j))
        assert actual==expected
        checks+=1;zero_div+=0 in actual and any(actual);tail+=w%H!=0
    print(f'Block norms: {checks} exact product-ring comparisons, {zero_div} nonzero zero-divisor outputs, {tail} tail cases.',flush=True)

# 2. Direct sparse integer polynomial tests of the (p,F)^t leading-term theorem.
def add(f,g):
    h=f.copy()
    for e,c in g.items():
        h[e]=h.get(e,0)+c
        if not h[e]:del h[e]
    return h

def mul(f,g):
    h={}
    for a,c in f.items():
        for b,d in g.items():
            e=tuple(x+y for x,y in zip(a,b));h[e]=h.get(e,0)+c*d
    return {e:c for e,c in h.items() if c}

def scale_shift(f,c,a):
    return {tuple(x+y for x,y in zip(e,a)):v*c for e,v in f.items() if v*c}

def lifting_tests():
    checks=decomps=cancels=0
    for _ in range(400):
        m=rng.choice([2,3]); d=rng.choice([1,2]); p=rng.choice([5,7,11,13]);t=rng.randrange(1,5)
        corner=(d,)*m;zero=(0,)*m
        F={e:rng.randrange(-p,p+1) for e in product(range(d+1),repeat=m)}
        F={e:c for e,c in F.items() if c};F[corner]=1
        powers=[{zero:1}]
        for j in range(t+3):powers.append(mul(powers[-1],F))
        g={}
        for j in range(rng.randrange(2,9)):
            k=rng.randrange(t+2);alpha=tuple(rng.randrange(3) for i in range(m));c=rng.choice([-3,-2,-1,1,2,3])
            g=add(g,scale_shift(powers[k],c*p**max(t-k,0),alpha))
        # Insert two cancelling shifts with arbitrarily large pivot coefficients.
        x=scale_shift(powers[t],13,(3,)*m);g=add(add(g,x),scale_shift(x,-1,zero));cancels+=1
        if not g:continue
        e=max(g);j=min(ei//d for ei in e)
        assert g[e] % p**max(t-j,0)==0
        X=tuple(rng.randrange(2,5) for i in range(m));X=(max(X[0],p),)+X[1:]
        assert prod(x**d for x in X)>=p
        assert abs(g[e])*prod(x**a for x,a in zip(X,e))>=p**t
        checks+=1
        # Exact decomposition in F-adic monomial basis. Lex leading term strictly decreases.
        rem=g.copy();rebuild={}
        while rem:
            e=max(rem);c=rem[e];j=min(ei//d for ei in e)
            beta=tuple(ei-j*d for ei in e)
            while len(powers)<=j:powers.append(mul(powers[-1],F))
            h=scale_shift(powers[j],c,beta)
            assert max(h)==e and h[e]==c
            assert c % p**max(t-j,0)==0
            rem=add(rem,scale_shift(h,-1,zero));rebuild=add(rebuild,h)
            decomps+=1
        assert rebuild==g
    print(f'Integer lifting: {checks} nonzero ideal elements, {decomps} exact F-adic coefficient checks, {cancels} cancellation controls.',flush=True)

# 3. Exact ideal index in the integer grid quotient, via an independent Smith form.
def grid_tests():
    checks=singleton=neg=multi=0
    for _ in range(60):
        p=rng.choice([7,11,13]);widths=rng.choice([(2,2),(2,3),(3,3),(2,2,2)])
        S=[rng.sample(range(p),w) for w in widths];grid=list(product(*S)); exps=list(product(*[range(w) for w in widths]));W=len(grid)
        E=Matrix([[prod(x**e for x,e in zip(s,a)) for a in exps] for s in grid])
        assert int(E.det())%p
        # Arbitrary polynomial with chosen zero patterns, by integer CRT modulo p.
        H0=rng.choice([0,1,1,2,W]); zeros=set(rng.sample(range(W),H0))
        target=Matrix([0 if i in zeros else rng.randrange(1,p) for i in range(W)])
        # Interpolate coefficient vector modulo p (rational inverse denominator is a p-unit).
        c=E.inv()*target
        coeff=Matrix([int(v.p)*pow(int(v.q),-1,p)%p for v in c])
        vals=E*coeff;assert {i for i,v in enumerate(vals) if int(v)%p==0}==zeros
        A=E.inv()*Matrix.diag(*vals)*E
        assert all(x.q==1 for x in A)
        for t in [1,2,3]:
            gens=Matrix.hstack(*[p**(t-j)*A**j for j in range(t+1)])
            D=smith_normal_form(gens,domain=ZZ)
            assert abs(prod(int(D[i,i]) for i in range(W)))==p**(t*H0)
            checks+=1
        if H0==1:
            s=grid[next(iter(zeros))]
            # The coefficients of the tensor-product Lagrange polynomial.
            col=E.inv()[:,next(iter(zeros))]
            scaled=Matrix([int(v.p)*pow(int(v.q),-1,p)%p for v in col])
            assert scaled[-1]%p!=0
            assert all(int(v)%p==0 for v in A*scaled)
            singleton+=1
        elif H0==0:neg+=1
        else:multi+=1
    print(f'Grid quotient: {checks} Smith-normal-form index comparisons; {singleton} singleton annihilators, {neg} empty and {multi} multiple-zero fixtures.',flush=True)

# 4. Actual Semaev polynomials, checked numerically with declared-degree resultants
# and exact finite-field interpolation. Helpers are independent of the block norm.
def det(M,p):
    M=[list(row) for row in M];n=len(M);out=1
    for j in range(n):
        k=next((k for k in range(j,n) if M[k][j]%p),None)
        if k is None:return 0
        if k!=j:M[j],M[k]=M[k],M[j];out=-out
        v=M[j][j]%p;out=out*v%p;iv=pow(v,-1,p)
        for k in range(j+1,n):
            q=M[k][j]*iv%p
            if q:
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
    # Newton forward differences at consecutive nodes, retaining all declared degrees.
    n=len(vals);diff=vals[:];co=[0]*n;basis=[1];fac=1
    for j in range(n):
        k=diff[0]*pow(fac,-1,p)%p
        for i,v in enumerate(basis):co[i]=(co[i]+k*v)%p
        diff=[(v-u)%p for u,v in zip(diff,diff[1:])]
        if j+1<n:
            b=[0]*(len(basis)+1)
            for i,v in enumerate(basis):b[i]=(b[i]-j*v)%p;b[i+1]=(b[i+1]+v)%p
            basis=b;fac=fac*(j+1)%p
    return co

def ev(f,x,p):
    z=0
    for v in f[::-1]:z=(z*x+v)%p
    return z

def semaev(p,a,c):
    @lru_cache(None)
    def partial(xs):
        if len(xs)==1:return ((-xs[0])%p,1)
        if len(xs)==2:return tuple(s3(*xs,a,c,p))
        prev=partial(xs[:-1]);d=2**(len(xs)-1)
        return tuple(interp([res(prev,s3(xs[-1],y,a,c,p),p) for y in range(d+1)],p))
    def F(xs,T):return ev(partial(tuple(xs)),T,p)
    return F

def degree_tests():
    checks=agg=corners=0
    for p,a,c,T in [(101,2,3,3),(257,1,1,0),(263,4,1,0)]:
        assert (4*a**3+27*c*c)%p and (T**3+a*T+c)%p and pow((T**3+a*T+c)%p,(p-1)//2,p)==1
        F=semaev(p,a,c)
        for m in [2,3,4]:
            d=2**(m-1);ratios=set();cor=1
            for j in range(m):
                for _ in range(8):
                    other=tuple(rng.randrange(p) for _ in range(m-1))
                    f=interp([F(other[:j]+(z,)+other[j:],T) for z in range(d+1)],p)
                    want=F(other,T)**2%p
                    if want:ratios.add(f[-1]*pow(want,-1,p)%p)
                    else:assert f[-1]==0
                    assert ev(f,p-1,p)==F(other[:j]+(p-1,)+other[j:],T)
                    checks+=1
            assert len(ratios)==1;cor=next(iter(ratios))
            # The full corner coefficient is the product of the constants in
            # the recursive leading-term formula; m=1 corner is 1.
            assert cor!=0;corners+=1
            for k in range(1,m-1):
                S=[rng.sample(range(p),2) for i in range(k)]
                grid=list(product(*S));degree=d*len(grid)
                for _ in range(3):
                    other=tuple(rng.randrange(p) for _ in range(m-k-1))
                    def norm(z):return prod(F(s+other+(z,),T) for s in grid)%p
                    lead=prod(cor*F(s+other,T)**2 for s in grid)%p
                    f=interp([norm(z) for z in range(degree+1)],p)
                    assert f[-1]==lead and ev(f,p-1,p)==norm(p-1)
                    agg+=1
    print(f'Semaev degrees: {checks} exact leading-coefficient comparisons, {agg} collective degree/leading comparisons, {corners} arity/curve checks.',flush=True)

# 5. Balanced generic Kronecker reduction from multicolored clique.
def clique_tests():
    checks=assignments=isolated=empty=0
    for trial in range(22):
        k=3;n=2;w=n*n;m=k*(k+1)//2
        pairs=[(i,j) for i in range(k) for j in range(i+1,k)]
        codes=rng.sample(range(1,20),k*n);V=[codes[i*n:(i+1)*n] for i in range(k)];B=max(codes)+w+1
        E={pair:{(u,v) for u in V[pair[0]] for v in V[pair[1]] if trial==0 or (trial!=1 and rng.randrange(3)!=0)} for pair in pairs}
        cl=[a for a in product(*V) if all((a[i],a[j]) in E[i,j] for i,j in pairs)]
        slots=[v+list(range(B-w,B-w+(w-n))) for v in V]
        assert all(len(set(s))==w for s in slots)
        for i,j in pairs:
            es=[u+B*v for u,v in E[i,j]];slots.append(es+list(range(B*B,B*B+w-len(es))))
        p=int(nextprime(max(w**m,8*m*B**4)))
        for tau in [None, min(map(sum,cl)) if cl else sum(v[0] for v in V)]:
            def f(xs):
                out=sum((xs[k+l]-xs[i]-B*xs[j])**2 for l,(i,j) in enumerate(pairs))
                if tau is not None:out+=(sum(xs[:k])-tau)**2
                assert out<p
                return out%p
            roots=[xs for xs in product(*slots) if f(xs)==0];assignments+=w**m
            expected=[a for a in cl if tau is None or sum(a)==tau]
            assert {x[:k] for x in roots}==set(expected) and len(roots)==len(expected)
            checks+=1;isolated+=len(roots)==1;empty+=not roots
    print(f'Clique reduction: {checks} complete grid comparisons, {assignments} grid assignments, {isolated} singleton and {empty} empty instances.',flush=True)

if __name__=='__main__':
    assert __debug__
    block_tests();lifting_tests();grid_tests();degree_tests();clique_tests()
    print('PASS: all exact component assertions.')
