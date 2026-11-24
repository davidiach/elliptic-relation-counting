#!/usr/bin/env python3
"""Exact finite checks for the numerical-sketch checkpoint.
No cryptographic-scale benchmark; all probability comparisons use integers.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from functools import lru_cache
from math import prod, isqrt
import random


def rank(cols, p, d):
    a = [list(x) for x in cols]
    t = 0
    for j in range(d):
        pivot = next((i for i in range(t, len(a)) if a[i][j] % p), None)
        if pivot is None:
            continue
        a[t], a[pivot] = a[pivot], a[t]
        v = pow(a[t][j] % p, -1, p)
        a[t] = [(x*v)%p for x in a[t]]
        for i in range(t+1,len(a)):
            v=a[i][j]%p
            if v: a[i]=[(x-v*y)%p for x,y in zip(a[i],a[t])]
        t+=1
    return t


def packed_bases(columns, p, d, t, k):
    if t==0: return [()] * k
    bases=[s for s in combinations(range(len(columns)),t)
           if rank([columns[i][1] for i in s],p,d)==t]
    masks=[sum(1<<i for i in s) for s in bases]
    @lru_cache(None)
    def rec(start, used, left):
        if not left: return ()
        if len(columns)-used.bit_count()<left*t: return None
        for j in range(start,len(bases)):
            if masks[j]&used: continue
            tail=rec(j+1,used|masks[j],left-1)
            if tail is not None: return (bases[j],)+tail
        return None
    return rec(0,0,k)


def peel(blocks,p,d,k):
    cols=[(i,c) for i,block in enumerate(blocks) for c in block]
    B=set(); original=rank([c for _,c in cols],p,d)
    while True:
        t=rank([c for _,c in cols],p,d)
        if not t: return B,cols,[()]*k,original
        bad=None
        for h in range(t):
            for ids in combinations(range(len(cols)),h):
                basis=[cols[i][1] for i in ids]
                if rank(basis,p,d)!=h: continue
                outside=[i for i,(_,c) in enumerate(cols) if rank(basis+[c],p,d)>h]
                if len(outside)<k*(t-h):
                    bad=(outside,h); break
            if bad is not None: break
        if bad is None:
            bases=packed_bases(cols,p,d,t,k)
            assert bases is not None
            return B,cols,bases,original
        outside,h=bad
        bnew={cols[i][0] for i in outside}
        assert bnew and len(bnew)<k*(t-h)
        B.update(bnew)
        cols=[x for x in cols if x[0] not in bnew]
        assert rank([c for _,c in cols],p,d)<=h


def check_peeling():
    rng=random.Random(614201)
    checks=0; dense=0; removed=0
    for _ in range(75):
        p=rng.choice([3,5,7]); d=rng.choice([1,2,3]); k=rng.choice([2,3])
        blocks=[[tuple(rng.randrange(p) for _ in range(d)) for _ in range(rng.choice([1,2]))]
                for j in range(rng.randrange(2,8))]
        B,C,bases,t0=peel(blocks,p,d,k)
        assert len(B)<=k*t0
        t=rank([c for _,c in C],p,d)
        assert len(bases)==k
        ids=[i for b in bases for i in b]
        assert len(ids)==len(set(ids))
        assert all(len(b)==t and rank([C[i][1] for i in b],p,d)==t for b in bases)
        # The sum of active scalar columns of the packed bases is at most
        # block width times the number of active blocks, for every dual vector.
        for v in product(range(p),repeat=d):
            h=sum(any(sum(x*y for x,y in zip(v,c))%p for c in block)
                  for j,block in enumerate(blocks) if j not in B)
            used=sum(bool(sum(x*y for x,y in zip(v,C[i][1]))%p) for i in ids)
            assert used<=max(map(len,blocks))*h
            checks+=1
        dense+=bool(t); removed+=len(B)
    print('Packing/peeling fixtures:',75,'nonzero surviving spaces:',dense,'peeled blocks:',removed)
    print('Exact dual-support inequalities:',checks)


P=13; A=4; C=1
O=None

def add(u,v):
    if u is None:return v
    if v is None:return u
    x,y=u;z,t=v
    if x==z:
        if (y+t)%P==0:return None
        s=(3*x*x+A)*pow(2*y,-1,P)%P
    else:s=(t-y)*pow(z-x,-1,P)%P
    q=(s*s-x-z)%P
    return q,(s*(x-q)-y)%P

pts=[O]+[(x,y) for x in range(P) for y in range(P) if (y*y-x*x*x-A*x-C)%P==0]
R=len(pts); assert R==19 and all(R%j for j in range(2,isqrt(R)+1))
g=pts[1]; group=[None]
for _ in range(1,R): group.append(add(group[-1],g))
assert len(set(group))==R and add(group[-1],g) is None
coords=[(0,0) if u is None else u for u in group]

# Work with distributions on F_p^2 indexed by x*p+y, exact integer masses.
S=P*P
shift=[[((i//P+j//P)%P)*P+(i%P+j%P)%P for j in range(S)] for i in range(S)]
def convolution(a,b):
    ans=[0]*S
    aa=[(i,v) for i,v in enumerate(a) if v]
    bb=[(i,v) for i,v in enumerate(b) if v]
    for i,v in aa:
        row=shift[i]
        for j,w in bb: ans[row[j]]+=v*w
    return tuple(ans)
unit=(1,)+(0,)*(S-1)
# A type is the choice of the second scalar coordinate of a measurement.
# The first coordinate is always x.
features={
    'Y':tuple(x*P+y for x,y in coords),
    'X':tuple(x*P+x for x,y in coords),
    '0':tuple(x*P for x,y in coords),
    'Q':tuple(x*P+(x*x)%P for x,y in coords),
}
hists={k:tuple(Counter(v).get(j,0) for j in range(S)) for k,v in features.items()}
@lru_cache(None)
def power(k,e):
    if e==0:return unit
    if e==1:return hists[k]
    a=power(k,e//2)
    a=convolution(a,a)
    return convolution(a,hists[k]) if e%2 else a

@lru_cache(None)
def independent(counts):
    out=unit
    for k,e in counts: out=convolution(out,power(k,e))
    return out

def counted(cnt):return tuple(sorted((k,v) for k,v in cnt.items() if v))

@lru_cache(None)
def selected_distribution(types,signs,target):
    result=[0]*S
    for i,j in product(range(R),repeat=2):
        z=(signs[2]*(target-signs[0]*i-signs[1]*j))%R
        a,b,c=(features[k][x] for k,x in zip(types,(i,j,z)))
        result[shift[shift[a][b]][c]]+=1
    assert sum(result)==R**2
    return tuple(result)


def check_adaptive():
    w=32; n=3*w; target=5; W=(2*w)**3
    # Second query determined by first answer modulo three.
    branches=[([Counter({'Y':w}) for i in range(3)],set()),
              ([Counter({'X':2,'0':w-2}) for i in range(3)],{'X'}),
              ([Counter({'Y':w}),Counter({'0':w}),Counter({'0':w})],{'Y'})]
    # Each branch has a known peeled set and enough untouched x or (x,y)
    # inputs to be nearly uniform. Use an explicit tolerance; no floating
    # approximation or assumed character bound in these finite comparisons.
    eps=Fraction(1,10**10)
    ptrans=[0]*S;qtrans=[0]*S
    good_checks=0; relative_checks=0; type_cases=0; good_trans=[0]*S
    pden=R**n;qden=W*R**(n-1)
    for branch,(slots,Btypes) in enumerate(branches):
        full=Counter()
        for slot in slots: full.update(slot)
        fresh=independent(counted(full)); assert sum(fresh)==pden
        if branch==0:
            ref=[Fraction(1,S)]*S
        else:
            bdist=independent(counted({k:v for k,v in full.items() if k in Btypes}))
            bden=R**sum(v for k,v in full.items() if k in Btypes)
            ycount=[sum(bdist[x*P+y] for x in range(P)) for y in range(P)]
            ref=[Fraction(ycount[j%P],bden*P) for j in range(S)]
        for j in range(S):
            v=Fraction(fresh[j],pden)
            assert (1-eps)*ref[j]<=v<=(1+eps)*ref[j]
            relative_checks+=1
        mixed=[0]*S; good=[0]*S; good_mult=0
        for types in product(*(tuple(s) for s in slots)):
            mult=prod(slot[k] for slot,k in zip(slots,types))
            cnt=full.copy();cnt.subtract(types)
            rest=independent(counted(cnt))
            covered=all(k in Btypes for k in types)
            for signs in product((-1,1),repeat=3):
                cur=convolution(rest,selected_distribution(types,signs,target))
                assert sum(cur)==R**(n-1)
                type_cases+=1
                for j,v in enumerate(cur): mixed[j]+=mult*v
                if not covered:
                    good_mult+=mult
                    for j,v in enumerate(cur):
                        val=Fraction(v,R**(n-1))
                        assert (1-eps)*ref[j]<=val<=(1+eps)*ref[j]
                        good[j]+=mult*v; good_checks+=1
        assert sum(mixed)==qden
        assert good_mult==W-prod(2*sum(v for k,v in slot.items() if k in Btypes) for slot in slots)
        for j in range(S):
            if (j//P)%3==branch:
                ptrans[j]=fresh[j];qtrans[j]=mixed[j];good_trans[j]=good[j]
                assert Fraction(good[j],qden)>=(Fraction(good_mult,W)*(1-eps)/(1+eps))*Fraction(fresh[j],pden)
    assert sum(ptrans)==pden and sum(qtrans)==qden
    tv=sum(abs(Fraction(a,pden)-Fraction(b,qden)) for a,b in zip(ptrans,qtrans))/2
    bad=1-Fraction(sum(good_trans),qden)
    theta=Fraction(2,w)**3
    bound=theta+2*eps/(1+eps)
    assert tv<=bound and bad<=bound
    print('Adaptive two-query curve: p=',P,'r=',R,'m=3, w=',w,sep='')
    print('Fresh relative comparisons:',relative_checks,'good-plant comparisons:',good_checks,'signed type cases:',type_cases)
    print('Adaptive transcript probabilities:',S,'all normalized exactly')
    print('Observed TV (rounded for display only):',float(tv))
    print('Covered-address probability (rounded):',float(bad))
    print('Certified TV and covered-address upper bound:',float(bound))


def check_moment_smoothing():
    print('Full two-power-sum vector (x, x^2), maximum relative deviation from uniform:')
    for w in [4,8,16,32,64,128]:
        law=power('Q',w-1);den=R**(w-1)
        err=max(abs(Fraction(S*v,den)-1) for v in law)
        # Conditioning on any marked row only translates the distribution.
        for v in set(features['Q']):
            shifted=[0]*S
            for j,c in enumerate(law):shifted[shift[v][j]]=c
            assert max(abs(Fraction(S*c,den)-1) for c in shifted)==err
        print(' w=',w,'relative deviation=',float(err),sep='')
    assert err<Fraction(1,10**9)
    print('All point-conditioned moment translations verified exactly.')

def check_summary_bound():
    w=3
    # Each of three arbitrary target-independent summaries emits one bit.
    # The summaries intentionally include nonsymmetric, nonlinear functions.
    rules=[lambda v: coords[v[0]][0]%2,
           lambda v: sum(coords[j][0] for j in v)%2,
           lambda v: int(v[0]==v[1] or coords[v[2]][1]%2==1)]
    ps=[];mark=[];deviations=[]
    den=2*w*R**w
    for rule in rules:
        c=[[0]*R for _ in range(2)]
        for v in product(range(R),repeat=w):
            z=rule(v)
            for a in v:
                c[z][a]+=1;c[z][(-a)%R]+=1
        assert sum(map(sum,c))==den
        p=[Fraction(sum(row),den) for row in c]
        joint=[[Fraction(v,den) for v in row] for row in c]
        delta=sum(abs(joint[z][a]-p[z]/R) for z in range(2) for a in range(R))
        assert delta**2<=Fraction(2,w) # use ln(2)<1 for a rational bound
        ps.append(p);mark.append(joint);deviations.append(delta)
    def group_conv(a,b):
        return [sum(a[i]*b[(j-i)%R] for i in range(R)) for j in range(R)]
    tv=Fraction(0);total=Fraction(0)
    for z in product(range(2),repeat=3):
        law=group_conv(group_conv(mark[0][z[0]],mark[1][z[1]]),mark[2][z[2]])
        pz=prod(ps[i][z[i]] for i in range(3))
        total+=sum(law)
        tv+=sum(abs(x-pz/R) for x in law)/2
    assert total==1 and tv<=prod(deviations)/2
    assert (2*tv)**2<=Fraction(2,w)**3
    print('Nonlinear local-summary fixtures: 3 one-bit summaries; joint target/sketch cells:',R*8)
    print('Exact joint-law TV (rounded):',float(tv),'product bound:',float(prod(deviations)/2))

if __name__=='__main__':
    check_peeling();check_adaptive();check_moment_smoothing();check_summary_bound();print('ALL CHECKS PASSED')
