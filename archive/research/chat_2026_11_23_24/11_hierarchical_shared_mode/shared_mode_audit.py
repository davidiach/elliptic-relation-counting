#!/usr/bin/env python3
"""Exact shared-mode Newton-interpolation checks. Standard library only."""
from itertools import product
from random import Random
from fractions import Fraction
from functools import lru_cache
from math import prod


def rem2(a,b):
    while a and a.bit_length()>=b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a

def irreducible(h):
    for f in range((1<<h)|1,1<<(h+1),2):
        if all(rem2(f,g) for d in range(1,h//2+1) for g in range(1<<d,1<<(d+1))):return f
    raise AssertionError(h)

class Ring:
    def __init__(self,h,P):
        self.h=h;self.P=P;self.mod=1<<P
        f=irreducible(h);self.f=tuple((f>>i)&1 for i in range(h))
        self.zero=(0,)*h;self.one=(1,)+(0,)*(h-1)
    def scalar(self,n):return (n%self.mod,)+(0,)*(self.h-1)
    def rational(self,x):
        assert x.denominator%2
        return self.scalar(x.numerator*pow(x.denominator,-1,self.mod))
    def bits(self,x):return tuple((x>>i)&1 for i in range(self.h))
    def add(self,a,b):return tuple((x+y)%self.mod for x,y in zip(a,b))
    def sub(self,a,b):return tuple((x-y)%self.mod for x,y in zip(a,b))
    @lru_cache(maxsize=100000)
    def mul(self,a,b):
        h=self.h;c=[0]*(2*h-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):c[i+j]+=x*y
        for j in range(2*h-2,h-1,-1):
            if c[j]:
                for k,f in enumerate(self.f):
                    if f:c[j-h+k]-=c[j]
        return tuple(x%self.mod for x in c[:h])
    @lru_cache(None)
    def inv(self,a):
        aa=tuple(x&1 for x in a);assert any(aa)
        for j in range(1,1<<self.h):
            z=self.bits(j)
            if tuple(x&1 for x in self.mul(aa,z))==self.one:break
        else:raise AssertionError('no inverse')
        for _ in range((self.P-1).bit_length()):z=self.mul(z,self.sub(self.scalar(2),self.mul(a,z)))
        assert self.mul(a,z)==self.one;return z
    def product(self,values):
        out=self.one
        for x in values:out=self.mul(out,x)
        return out
    def sum(self,values):
        out=self.zero
        for x in values:out=self.add(out,x)
        return out


def lower_set(weights,P):
    def rec(j,budget,prefix):
        if j==len(weights):yield prefix;return
        for e in range(budget//weights[j]+1):yield from rec(j+1,budget-e*weights[j],prefix+(e,))
    return list(rec(0,P-1,()))


def compare(R,weights,source,forms,local,cross,table):
    """One X and one Y port per group, arbitrary local weight. Shared modes can have many targets."""
    g=len(local);n=len(weights);P=R.P
    assert all(cross[i][j]==R.zero or all(v%(1<<weights[j])==0 for v in cross[i][j])
               for i in range(g) for j in range(n))
    assert all(cross[source[j]][j]==R.zero for j in range(n))
    def row(i,x,y,mode_vals):
        s=R.mul(local[i],R.scalar(x))
        s=R.add(s,R.sum(R.mul(c,z) for c,z in zip(cross[i],mode_vals)))
        return R.mul(table[i][2*x+y],s if y else R.sub(R.one,s))
    expected=R.zero
    for X in product((0,1),repeat=g):
        Z=[forms[j][X[source[j]]] for j in range(n)]
        for Y in product((0,1),repeat=g):expected=R.add(expected,R.product(row(i,X[i],Y[i],Z) for i in range(g)))
    # Unit-separated grid and tensor-product divided differences on a lower set.
    K=lower_set(weights,P);knots=[R.bits(i) for i in range(P)]
    den={};newton={}
    for b in range(P):
        for a in range(b+1):
            den[b,a]=R.inv(R.product(R.sub(knots[a],knots[h]) for h in range(b+1) if h!=a))
        for j in range(n):
            for x in (0,1):newton[j,b,x]=R.product(R.sub(forms[j][x],knots[h]) for h in range(b))
    target={};src={}
    for a in K:
        z=[knots[e] for e in a]
        target[a]=[[row(i,x,y,z) for x,y in product((0,1),repeat=2)] for i in range(g)]
        src[a]=[[R.product(newton[j,a[j],x] for j in range(n) if source[j]==i) for x in (0,1)] for i in range(g)]
    got=R.zero;terms=0
    for beta in K:
        for alpha in product(*(range(b+1) for b in beta)):
            w=R.product(den[b,a] for b,a in zip(beta,alpha))
            Z=R.product(R.sum(R.mul(target[alpha][i][2*x+y],src[beta][i][x])
                              for x,y in product((0,1),repeat=2)) for i in range(g))
            got=R.add(got,R.mul(w,Z));terms+=1
    assert got==expected
    assert terms==sum(prod(b+1 for b in beta) for beta in K)
    fanout=max(sum(cross[i][j]!=R.zero for i in range(g)) for j in range(n))
    return len(K),terms,fanout,any(expected)


def random_tests():
    rng=Random(2026112403);nodes=terms=fanout=nonzero=0
    for trial in range(40):
        P=rng.randrange(3,6);R=Ring(3,P);g=rng.randrange(3,5);n=rng.randrange(2,5)
        weights=[rng.randrange(1,3) for _ in range(n)];source=[rng.randrange(g) for _ in range(n)]
        forms=[(R.zero,R.scalar(rng.randrange(1,8))) for _ in range(n)]
        local=[R.scalar(rng.randrange(16)) for _ in range(g)]
        cross=[[R.zero if source[j]==i else R.scalar((1<<weights[j])*rng.randrange(1,1<<(P-weights[j]),2))
                for j in range(n)] for i in range(g)]
        table=[[R.scalar(rng.randrange(1<<P)) for _ in range(4)] for i in range(g)]
        k,t,f,nz=compare(R,weights,source,forms,local,cross,table)
        nodes+=k;terms+=t;fanout=max(fanout,f);nonzero+=nz
    return nodes,terms,fanout,nonzero


def W_value(L,T,delta):
    S=range(-L,L);C=prod(T-s for s in S);den=prod(-s for s in S if s)
    return Fraction(C,den*(T-delta))


def cauchy_tests():
    # Full hierarchical packet of 4 positions. At level 1 the two low-bit
    # prefixes use disjoint group pairs. Each source tile is a singleton.
    # At T=2, local edges stay inside each group. Nearby they cross to 3 groups.
    rng=Random(2026112404);a=(0,2,1,3);b=(1,3,0,2);L=4;t=2
    assert all(a[x]==b[y] for x in range(L) for y in range(L) if abs(y-x)==2)
    nodes=terms=fanout=nonzero=0
    for trial in range(12):
        P=rng.randrange(3,6);R=Ring(3,P);T=t+2*L*(1+trial)
        rows={a[x]:x for x in range(L)};cols={b[y]:y for y in range(L)}
        source=list(range(L));forms=[(R.zero,R.one)]*L;weights=[1]*L
        local=[R.rational(W_value(L,T,cols[i]-rows[i])) for i in range(L)]
        cross=[[R.zero if i==j else R.rational(W_value(L,T,cols[j]-rows[i])) for j in range(L)] for i in range(L)]
        table=[[R.scalar(rng.randrange(1<<P)) for _ in range(4)] for i in range(L)]
        # Direct row values above are actual Cauchy entries, not synthetic modes.
        k,s,f,nz=compare(R,weights,source,forms,local,cross,table)
        nodes+=k;terms+=s;fanout=max(fanout,f);nonzero+=nz
    return nodes,terms,fanout,nonzero


if __name__=='__main__':
    if not __debug__:raise RuntimeError('assertions must be enabled')
    print('40 random shared-mode joins: grid nodes, Newton pairs, max fanout, nonzero outputs:',*random_tests(),flush=True)
    print('12 actual hierarchical Cauchy joins: grid nodes, Newton pairs, max fanout, nonzero outputs:',*cauchy_tests(),flush=True)
    print('All assertions passed.',flush=True)
