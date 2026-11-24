#!/usr/bin/env python3
"""Exact finite checks for the seed-local degree/edge tradeoff.
No external packages. This does not benchmark a contractor or prove an ECDLP bound.
"""
from itertools import product, combinations
from math import comb, isqrt, log2
from random import Random
from collections import Counter
from functools import lru_cache


def ec_add(P,Q,p,a):
    if P is None: return Q
    if Q is None: return P
    x,y=P; u,v=Q
    if x==u and (y+v)%p==0: return None
    s=((v-y)*pow(u-x,-1,p) if x!=u else (3*x*x+a)*pow(2*y,-1,p))%p
    w=(s*s-x-u)%p
    return (w,(s*(x-w)-y)%p)

def points(p,a,c):
    return [None]+[(x,y) for x in range(p) for y in range(p) if (y*y-x*x*x-a*x-c)%p==0]

def is_prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def gf2_rank(rows):
    basis={}
    for x in rows:
        while x:
            i=x.bit_length()-1
            if i in basis: x^=basis[i]
            else:
                basis[i]=x; break
    return len(basis)

# Polynomial arithmetic over GF(2), for explicit irreducible extension construction.
def pmul(a,b):
    z=0
    while b:
        if b&1:z^=a
        a<<=1;b>>=1
    return z

def prem(a,b):
    while a and a.bit_length()>=b.bit_length():
        a^=b<<(a.bit_length()-b.bit_length())
    return a

def pgcd(a,b):
    while b:a,b=b,prem(a,b)
    return a

def irreducible(f,d):
    x=2; y=x
    for j in range(1,d+1):
        y=prem(pmul(y,y),f)
        if j<=d//2 and pgcd(y^x,f)!=1: return False
    return y==x

class GF:
    def __init__(self,d):
        self.d=d; self.size=1<<d
        self.f=next(f for f in range(self.size+1,2*self.size,2) if irreducible(f,d))
    @lru_cache(maxsize=200000)
    def mul(self,a,b):return prem(pmul(a,b),self.f)
    def pow(self,a,k):
        z=1
        while k:
            if k&1:z=self.mul(z,a)
            a=self.mul(a,a);k>>=1
        return z
    def root(self,r):
        assert (self.size-1)%r==0 and is_prime(r)
        for a in range(2,self.size):
            z=self.pow(a,(self.size-1)//r)
            if z!=1:
                assert self.pow(z,r)==1
                return z
        raise AssertionError('no root')
    def rank(self,rows):
        rows=[list(row) for row in rows]; rank=0
        for j in range(len(rows[0])):
            k=next((k for k in range(rank,len(rows)) if rows[k][j]),None)
            if k is None:continue
            rows[rank],rows[k]=rows[k],rows[rank]
            inv=self.pow(rows[rank][j],self.size-2)
            rows[rank]=[self.mul(x,inv) for x in rows[rank]]
            for k in range(rank+1,len(rows)):
                a=rows[k][j]
                if a:rows[k]=[x^self.mul(a,y) for x,y in zip(rows[k],rows[rank])]
            rank+=1
            if rank==len(rows):break
        return rank

def check_curves():
    mask_cases=mask_images=diag_checks=fixed_checks=0
    fixtures=[]
    for p,a,c in [(5,2,1),(7,0,3),(13,4,1),(17,2,2)]:
        if (4*a*a*a+27*c*c)%p==0:continue
        G=points(p,a,c);r=len(G)
        if not is_prime(r) or r==2:continue
        fixtures.append((p,a,c,r))
    assert len(fixtures)>=3, fixtures
    for p,a,c,r in fixtures:
        G=points(p,a,c); gen=G[1]; P=None; seq=[]
        for j in range(r):seq.append(P);P=ec_add(P,gen,p,a)
        assert P is None and set(seq)==set(G)
        d=(r-1)//2; xs=[seq[j][0] for j in range(1,d+1)]
        assert len(set(xs))==d and d<p
        # Every nonempty proper support S (exhaustive for p<=7, otherwise fixed random fixtures).
        rng=Random(p)
        if p<=7: supports=[set(s) for k in range(1,p) for s in combinations(range(p),k)]
        else: supports=[set(rng.sample(range(p),k)) for k in range(1,p) for _ in range(2)]
        for S in supports:
            full=[];restricted=[]
            for aa in range(p):
                for t in range(1,p):
                    A={(aa+t*s)%p for s in S}
                    full.append(sum(1<<x for x in A))
                    restricted.append(sum(1<<i for i,x in enumerate(xs) if x in A))
                    mask_images+=1
            assert gf2_rank(full)==p-(len(S)%2==0)
            assert gf2_rank(restricted)==d
            mask_cases+=1
        f=1
        while pow(2,f,r)!=1:f+=1
        F=GF(f);z=F.root(r);powers=[F.pow(z,j) for j in range(r)]
        H=[[powers[(i*j)%r]^powers[(-i*j)%r] for j in range(1,d+1)] for i in range(1,d+1)]
        assert F.rank(H)==d
        # Direct signed point counts modulo 2 versus the paired Fourier diagonal.
        for m in (2,3,4):
            for ids in (product(range(d),repeat=m) if d**m<=1000 else [tuple(rng.randrange(d) for _ in range(m)) for _ in range(1000)]):
                count=Counter(sum((i+1)*s for i,s in zip(ids,signs))%r for signs in product((-1,1),repeat=m))
                coeff=[]
                for j in range(d):
                    v=1
                    for i in ids:v=F.mul(v,H[j][i])
                    coeff.append(v)
                for t in range(1,d+1):
                    got=0
                    for j,v in enumerate(coeff,1):got^=F.mul(v,H[j-1][t-1])
                    assert got==(count[t]%2), (p,m,ids,t,got,count[t])
                    diag_checks+=1
        # Nonempty mask Fourier support bound, and restriction by fixing inactive masks.
        for S in supports[:min(len(supports),25)]:
            A=next({(aa+t*s)%p for s in S} for aa in range(p) for t in range(1,p) if any(x in {(aa+t*s)%p for s in S} for x in xs))
            mask=[int(x in A) for x in xs];w=sum(mask)*2
            fh=[0]*d
            for j in range(d):
                for i,bit in enumerate(mask):
                    if bit:fh[j]^=H[j][i]
            live=sum(x!=0 for x in fh)
            assert 2*live*w>=r
            for i,k in product(range(d),repeat=2):
                signs=[(s1*(i+1)+s2*(k+1)+s3*(h+1))%r for h,bit in enumerate(mask) if bit for s1,s2,s3 in product((-1,1),repeat=3)]
                got=0
                for j in range(d):got^=F.mul(F.mul(F.mul(H[j][i],H[j][k]),fh[j]),H[j][0])
                assert got==(signs.count(1)%2)
                fixed_checks+=1
        print(f'Curve F_{p}, group order {r}: mask span {d}, Fourier extension degree {f}, paired transform invertible')
    print('Affine support families:',mask_cases,'; affine images:',mask_images)
    print('Exact paired Fourier parity comparisons:',diag_checks)
    print('Exact fixed-mask restrictions:',fixed_checks)

# Test the monomial-to-slice identity on random small seed-local Boolean networks.
# Row/column tensors are functions of (local port bits, own seed), with some
# edges incident to two distinct seed-dependent vertices. All arithmetic mod 2.
def parity_moment(table, mask, deg, seed):
    return sum(table[seed][x] for x in range(1<<deg) if x&mask==mask)&1

def check_slice_decomposition():
    rng=Random(20261124); ncases=assignments=terms=0
    for _ in range(100):
        m=3; edges=[(i,j) for i in range(3) for j in range(3)];M=len(edges)
        rports=[[k for k,(a,b) in enumerate(edges) if a==i] for i in range(3)]
        cports=[[k for k,(a,b) in enumerate(edges) if b==i] for i in range(3)]
        # Both sides may be seed dependent, with different slot labels.
        rseed=[rng.randrange(m+1) for _ in range(3)];cseed=[rng.randrange(m+1) for _ in range(3)]
        # m means seed-free, other labels correspond to actual slot seeds.
        rpoly=[]; ctables=[]; cap=[]
        for i in range(3):
            d=rng.randrange(3);cap.append(d)
            rpoly.append([[rng.randrange(2) if j.bit_count()<=d else 0 for j in range(8)] for _ in range(2)])
            ctables.append([[rng.randrange(2) for _ in range(8)] for _ in range(2)])
        D=sum(cap)
        E=[{k for k,(a,b) in enumerate(edges) if rseed[a]==i or cseed[b]==i} for i in range(m)]
        # Any set J of degree <=D has a slot with <=floor(2D/m) used ports.
        bound=2*D//m
        def localmask(J,ports):return sum(1<<j for j,e in enumerate(ports) if J>>e&1)
        def f_row(j,loc,seed):return rpoly[j][seed][loc]
        def f_col(j,loc,seed):return parity_moment(ctables[j],loc,3,seed)
        def features(J,seed):
            vals=[1]*(m+1)
            for j in range(3):
                i=rseed[j];vals[i]&=f_row(j,localmask(J,rports[j]),seed[i] if i<m else 0)
                i=cseed[j];vals[i]&=f_col(j,localmask(J,cports[j]),seed[i] if i<m else 0)
            return vals
        def rowval(j,x,seed):
            return sum(rpoly[j][seed][S] for S in range(8) if S&x==S)&1
        buckets={}
        for J in range(1<<M):
            if J.bit_count()>D or any(localmask(J,rports[j]).bit_count()>cap[j] for j in range(3)):continue
            choice=next(i for i in range(m) if sum(J>>e&1 for e in E[i])<=bound)
            local=tuple(e for e in sorted(E[choice]) if J>>e&1)
            buckets.setdefault((choice,local),[]).append(J);terms+=1
        assert len(buckets)<=sum(sum(comb(len(e),j) for j in range(min(bound,len(e))+1)) for e in E)
        for seeds in product(range(2),repeat=m):
            def get(i):return seeds[i] if i<m else 0
            direct=0
            for X in range(1<<M):
                v=1
                for j in range(3):
                    v&=rowval(j,localmask(X,rports[j]),get(rseed[j]))
                    v&=ctables[j][get(cseed[j])][localmask(X,cports[j])]
                direct^=v
            got=0
            for (i,local),Js in buckets.items():
                representative=features(Js[0],seeds)[i]
                rest=0
                for J in Js:
                    fs=features(J,seeds)
                    assert fs[i]==representative
                    v=1
                    for k,z in enumerate(fs):
                        if k!=i:v&=z
                    rest^=v
                got^=representative&rest
            assert got==direct
            assignments+=1
        ncases+=1
    print('Seed-local networks:',ncases,'; exact seed-tuple contractions:',assignments,'; grouped monomials:',terms)

# Elementary counting inequality for all n,k in a finite range.
def check_entropy():
    checks=0
    for n in range(1,257):
        total=0
        for k in range(n+1):
            total+=comb(n,k)
            assert log2(total)<=2*(n*k)**.5+1e-12
            checks+=1
    print('Finite binomial/degree inequalities:',checks)

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    check_curves();check_slice_decomposition();check_entropy()
    print('All checks passed. Component tests only; no improved ECDLP exponent claimed.')
