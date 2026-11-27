#!/usr/bin/env python3
"""Exact arithmetic checks of abstract Newton certificates and moment trades."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, prod

if not __debug__:
    raise RuntimeError('Run with assertions enabled')


def simplex(d, W):
    if d == 0:
        yield ()
        return
    for a in range(W+1):
        for rest in simplex(d-1, W-a):
            yield (a,) + rest


def signs(alpha):
    S = [i for i,v in enumerate(alpha) if v]
    for eps in product((-1,1), repeat=len(S)):
        x = list(alpha)
        for i,s in zip(S,eps):
            x[i] *= s
        yield tuple(x)


def dot(x,u):
    return sum(a*b for a,b in zip(x,u))


def desired(alpha,u,r,c):
    B = sum(a*a for a in alpha) if c == 0 else 1
    A = tuple(signs(alpha))
    bad = tuple(x for x in A if (dot(x,u)-c)%r)
    if not bad:
        assert B == 0
        return F(0), ()
    assert 2*len(bad) >= len(A) or not any(alpha)
    return F(B*len(A),len(bad)), bad


def newton(u,r,W,c):
    d = len(u)
    assert all(v%r for v in u) and 0 <= c < r and r > W
    terms = {}
    @lru_cache(None)
    def fourier(alpha,R):
        height,bad = desired(alpha,u,r,c)
        if not bad:
            return F(0)
        # Divide the Boolean Fourier coefficient by product(alpha_i for i in R).
        return height * F(sum(prod(1 if x[i]>0 else -1 for i in R) for x in bad),
                          2**sum(a!=0 for a in alpha)*prod(alpha[i] for i in R))
    for n in range(1,min(d,W)+1):
        if c == 0 and n%2:
            continue
        for R in combinations(range(d),n):
            eta = tuple(int(i in R) for i in range(d))
            for beta in simplex(d,W-n):
                coeff = F(0)
                for gamma in product(*(range(k+1) for k in beta)):
                    alpha = tuple(a+b for a,b in zip(gamma,eta))
                    denom = prod((a+s)**2-(j+s)**2 for a,s,k in zip(gamma,eta,beta)
                                 for j in range(k+1) if j != a)
                    coeff += fourier(alpha,R)/denom
                if coeff:
                    terms[R,beta] = coeff
    assert len(terms) <= 2**W*comb(d+W,W)
    return terms


def evaluate(x,terms,c):
    result = F(sum(a*a for a in x) if c == 0 else 1)
    for (R,beta),coeff in terms.items():
        factor = prod(x[i] for i in R)
        if not factor: continue
        factor *= prod(v*v-(j+int(i in R))**2 for i,(v,k) in enumerate(zip(x,beta)) for j in range(k))
        if factor: result += coeff*factor
    return result


def check_newton():
    print('DIRECT NEWTON CERTIFICATES')
    print('d | W | r | directions | targets | point checks | max terms | max degree')
    cases=[(1,1,3,[(1,)]),(1,2,3,[(1,),(2,)]),
           (2,2,5,[(1,1),(1,2),(2,3)]),(2,3,5,[(1,2),(2,3)]),
           (3,3,7,[(1,2,3),(1,1,2),(3,4,5)]),
           (3,4,7,[(1,2,3),(1,1,1),(2,3,4)]),
           (3,5,7,[(1,2,3)]),(4,3,11,[(1,2,3,4),(1,1,2,3)]),
           (4,4,11,[(1,2,3,4)])]
    total=0; points=0
    for d,W,r,dirs in cases:
        X=[x for a in simplex(d,W) for x in signs(a)]
        targets=(0,1,r-1)
        mx=degree=0
        for u in dirs:
            for c in targets:
                T=newton(u,r,W,c)
                deg=max([2 if c==0 else 0]+[len(R)+2*sum(beta) for R,beta in T])
                assert deg <= (max(2,2*W-2) if c==0 else 2*W-1)
                degree=max(degree,deg); mx=max(mx,len(T))
                # Check the exact prescribed function, not just its zero set.
                values={}
                for a in simplex(d,W):
                    height,bad=desired(a,u,r,c); bad=set(bad)
                    for x in signs(a): values[x]=height if x in bad else F(0)
                for x in X:
                    y=evaluate(x,T,c)
                    assert y==values[x] and y>=0
                    assert (y==0)==((dot(x,u)-c)%r==0)
                total+=1; points+=len(X)
        print(f'{d} | {W} | {r} | {len(dirs)} | {len(targets)} | {len(X)*len(dirs)*len(targets)} | {mx} | {degree}')
    print(f'{total} certificates; {points} exact rational point checks.')


def orbit(a,p=7):
    # Multiplication of the residue labels 1,2,3 by k acts as a signed permutation.
    ret=set()
    for k in range(1,p):
        x=[0]*len(a)
        for i,v in enumerate(a,1):
            s=k*i%p;j=min(s,p-s)
            x[j-1]=v if s==j else -v
        ret.add(tuple(x))
    return ret


def expand(spec):
    C=Counter()
    for a,n in spec:
        assert type(n) is int and n > 0, f'Orbit {a} requires a positive integer multiplicity, got {n!r}'
        for x in orbit(a): C[x]+=n
    return C


def moments(C,d,t):
    return tuple(sum(n*prod(v**e for v,e in zip(x,a)) for x,n in C.items())
                 for a in simplex(d,t))


def rank(A,p):
    piv={}
    for row in A:
        x=[v%p for v in row]
        for j,b in sorted(piv.items()):
            if x[j]:
                c=x[j];x=[a-c*z for a,z in zip(x,b)];x=[a%p for a in x]
        j=next((j for j,v in enumerate(x) if v),None)
        if j is not None:
            c=pow(x[j],-1,p);piv[j]=[v*c%p for v in x]
    return len(piv)


H4=[((-2,-1,-1),1),((-2,1,0),2)]
J4=[((-2,0,-1),2),((-2,1,1),1)]
H5=[((-4,0,-1),7718),((-3,-2,0),8568),((-3,0,1),26400),
    ((-2,-1,-1),55440),((-2,1,0),4620),((-1,-1,1),166320)]
J5=[((-5,0,0),1379),((-4,-1,0),945),((-4,0,0),24),((-4,0,1),288),
    ((-4,1,0),67),((-3,-1,-1),11130),((-3,-1,0),864),((-3,-1,1),1260),
    ((-3,0,-2),56),((-3,0,-1),17744),((-3,0,0),1719),((-3,0,2),3696),
    ((-3,1,-1),1344),((-3,1,1),42),((-2,-2,-1),3906),((-2,-2,0),17808),
    ((-2,-2,1),1050),((-2,-1,-2),2772),((-2,-1,0),5460),((-2,-1,1),7392),
    ((-2,0,-2),280),((-2,0,-1),420),((-2,0,1),41580),((-1,-1,0),36960)]


def check_trades():
    print('\nSHARP HOMOGENEOUS TRADES (compressed multiplicities, r=7)')
    print('W | N per population | distinct H/J | equal moment degree | ranks H/J/union mod 7')
    for W,H0,J0 in [(4,H4,J4),(5,H5,J5)]:
        H,J=expand(H0),expand(J0)
        assert max(sum(map(abs,x)) for x in H|J)==W
        assert all(dot(x,(1,2,3))%7==0 for x in H)
        assert all(dot(x,(1,2,3))%7!=0 for x in J)
        assert sum(H.values())==sum(J.values())
        assert moments(H,3,2*W-3)==moments(J,3,2*W-3)
        assert moments(H,3,2*W-2)!=moments(J,3,2*W-2)
        assert rank(H,1009)==rank(J,1009)==3
        assert rank(H,7)==2 and rank(H|J,7)==3
        print(f'{W} | {sum(H.values())} | {len(H)}/{len(J)} | {2*W-3} | 2/{rank(J,7)}/3')
    H=Counter({(-2,):1,(1,):2});J=Counter({(-1,):2,(2,):1})
    assert all(x[0]%3==1 for x in H) and all(x[0]%3!=1 for x in J)
    assert moments(H,1,2)==moments(J,1,2) and moments(H,1,3)!=moments(J,1,3)
    print('Affine W=2, r=3: equal moments through 2, different at 3.')


def parity_tags(t,h,code):
    k=t+1
    groups=[]
    for i in range(h):
        groups.append([x for x in product((0,1),repeat=k) if sum(x)%2==(code>>i)&1])
    return [sum(v,()) for v in product(*groups)]


def check_tags():
    # Test finite tag-marginal equality without expanding the massive sharpness examples.
    print('\nINACTIVE TAG DISTINCTIFICATION')
    n=0
    for t,h in [(2,2),(3,2),(5,1),(7,1)]:
        tags=[parity_tags(t,h,c) for c in range(2**h)];d=h*(t+1)
        assert len(set(x for S in tags for x in S))==sum(map(len,tags))
        ref=tags[0]
        for S in tags:
            for tag in S:
                assert len(tag) == d, f'Expected a tag of length {d}, got {tag!r}'
                assert all(type(bit) is int and bit in (0, 1) for bit in tag), f'Expected a binary tag, got {tag!r}'
            assert rank([x+(1,) for x in S],7)==d+1
            assert rank([x+(1,) for x in S],1009)==d+1
            for j in range(t+1):
                for I in combinations(range(d),j):
                    assert sum(prod(x[i] for i in I) for x in S)==sum(prod(x[i] for i in I) for x in ref)
                    n+=1
        print(f't={t}, h={h}: {len(tags)} disjoint classes, {len(ref)} tags/class, affine rank {d+1}')
    print(f'{n} exact tag marginal comparisons.')


if __name__=='__main__':
    check_newton();check_trades();check_tags();print('\nAll checks passed.')
