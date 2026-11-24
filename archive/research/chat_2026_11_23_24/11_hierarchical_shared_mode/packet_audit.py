#!/usr/bin/env python3
"""Exact regressions for hierarchical packing, entropy symmetry, and variable precision."""
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import product, combinations_with_replacement
from math import ceil, log2
from random import Random


def v2(x):
    assert x
    x=abs(x)
    return (x & -x).bit_length()-1


def direct_local(a,b,i):
    return all(a[x]==b[y] for x in range(len(a)) for y in range(len(a))
               if x!=y and v2(y-x)==i)


def prefix_swap(a,b,i):
    step=1<<(i+1)
    for x in range(len(a)):
        low=x%step
        if a[x]!=a[low] or b[x]!=a[low^(1<<i)]:
            return False
    return True


def locality_tests():
    cnt=0
    for a in (1,2,3):
        L=1<<a
        # Precompute all locality pairs for efficient exhaustive checks.
        pairs=[[ (x,y) for x in range(L) for y in range(L)
                if x!=y and v2(y-x)==i] for i in range(a)]
        for bits in product((0,1), repeat=2*L):
            r,c=bits[:L],bits[L:]
            for i in range(a):
                assert all(r[x]==c[y] for x,y in pairs[i])==prefix_swap(r,c,i)
                cnt+=1
    return cnt


def seq_counts(tuples):
    return Counter(tuples)


def transitions(counts,i):
    out=defaultdict(Counter)
    for s,n in counts.items():
        out[s[:i]][s[i]]+=n
    return out


def product_symmetry(A,B,u):
    for i in range(u):
        for ac in transitions(A,i).values():
            for bc in transitions(B,i).values():
                if ac[0]*bc[1]!=ac[1]*bc[0]:
                    return False
    return True


def independent_same(A,B,u):
    na,nb=sum(A.values()),sum(B.values())
    for i in range(u):
        qa=sum(n for s,n in A.items() if s[i])
        qb=sum(n for s,n in B.items() if s[i])
        if qa*nb!=qb*na:
            return False
        for C,q,N in ((A,qa,na),(B,qb,nb)):
            if any(c[1]*N != q*sum(c.values()) for c in transitions(C,i).values()):
                return False
    return True


def entropy_tests():
    tested=passed=0
    for u,n in ((2,4),(3,3)):
        labels=list(product((0,1),repeat=u))
        sets=[seq_counts(x) for x in combinations_with_replacement(labels,n)]
        for A in sets:
            for B in sets:
                good=product_symmetry(A,B,u)
                assert good==independent_same(A,B,u)
                tested+=1
                passed+=good
    return tested,passed


def min_packets(counts,a,u,depth=0):
    """Minimum subpacket count, fixed order of u label levels, each of 2^a slots."""
    if not counts or not any(counts.values()): return 0
    if depth==u:
        return (sum(counts.values())+(1<<(a-u))-1)//(1<<(a-u))
    groups=defaultdict(Counter)
    for s,n in counts.items():
        if n:
            groups[s[0]][s[1:]]+=n
    needs={k:min_packets(v,a,u,depth+1) for k,v in groups.items()}
    out=0
    for x,y in {tuple(sorted(k)) for k in needs}:
        if x==y: out+=(needs.get((x,x),0)+1)//2
        else: out+=max(needs.get((x,y),0),needs.get((y,x),0))
    return out


def all_patterns(u,r=2):
    if u==0: return [((),)]
    prev=all_patterns(u-1,r)
    result=set()
    for x,y in product(range(r),repeat=2):
        for p,q in product(prev,repeat=2):
            result.add(tuple(sorted(tuple(((x,y),)+z for z in p)+tuple(((y,x),)+z for z in q))))
    return sorted(result)


def packing_tests():
    rng=Random(20261124)
    u=a=2
    universe=list(product(list(product(range(2),repeat=2)),repeat=u))
    patterns=[]
    for p in all_patterns(u):
        c=Counter(p); patterns.append(tuple(c[x] for x in universe))
    @lru_cache(None)
    def exact(d):
        if not any(d):return 0
        i=next(j for j,x in enumerate(d) if x)
        return 1+min(exact(tuple(max(0,x-y) for x,y in zip(d,p))) for p in patterns if p[i])
    tests=0
    for _ in range(250):
        c=Counter(rng.choices(universe,k=rng.randrange(1,10)))
        d=tuple(c[x] for x in universe)
        assert min_packets(c,a,u)==exact(d)
        tests+=1
    records=[]
    for u in range(1,7):
        labels=list(product((0,1),repeat=u))
        c=Counter({tuple(zip(x,y)):1<<u for x,y in product(labels,repeat=2)})
        real=sum(c.values()); slots=(1<<u)*min_packets(c,u,u)
        assert real==slots
        for r in (2,3,5):
            c=Counter({((x,y),)*u:1<<u for x,y in product(range(r),repeat=2)})
            real=sum(c.values());slots=(1<<u)*min_packets(c,u,u)
            ratio=Fraction(slots,real)
            assert ratio==Fraction(1,r)+Fraction(r-1,r)*(1<<(u-1))
            records.append((u,r,ratio))
    return tests,len(patterns),records


def flat_min_packets(counts,a):
    groups=defaultdict(Counter)
    for seq,n in counts.items():
        key=tuple(tuple(sorted(pair)) for pair in seq)
        groups[key][seq]+=n
    total=0
    for key, C in groups.items():
        multiplicity=1<<(a-sum(x!=y for x,y in key))
        total+=max((n+multiplicity-1)//multiplicity for n in C.values())
    return total


def hierarchy_examples():
    records=[]
    for u in range(1,9):
        c=Counter()
        for bits in product((0,1),repeat=u):
            prefix=0;seq=[]
            for bit in bits:
                aa,bb=2*prefix,2*prefix+1
                seq.append((bb,aa) if bit else (aa,bb));prefix=2*prefix+bit
            c[tuple(seq)]+=1
        h=min_packets(c,u,u);f=flat_min_packets(c,u)
        assert h==1 and f==1<<(u-1)
        records.append((u,h,f))
    return records


def materialize(counts,a,u,depth=0):
    if depth==u:
        return [tuple([()]*(1<<(a-u)))]*min_packets(counts,a,u,depth)
    groups=defaultdict(Counter)
    for s,n in counts.items():
        if n:groups[s[0]][s[1:]]+=n
    sub={k:materialize(v,a,u,depth+1) for k,v in groups.items()}
    dummy=tuple([((0,0),)*(u-depth-1)]*(1<<(a-depth-1)))
    out=[]
    for x,y in sorted({tuple(sorted(k)) for k in sub}):
        l=sub.get((x,y),[]); r=sub.get((y,x),[])
        if x==y:
            for i in range(0,len(l),2):
                z=l[i+1] if i+1<len(l) else dummy
                out.append(tuple(((x,x),)+t for t in l[i]+z))
        else:
            for i in range(max(len(l),len(r))):
                L=l[i] if i<len(l) else dummy;R=r[i] if i<len(r) else dummy
                out.append(tuple(((x,y),)+t for t in L)+tuple(((y,x),)+t for t in R))
    return out


def H(counts):
    N=sum(counts.values())
    return -sum(n/N*log2(n/N) for n in counts.values() if n)


def padding_entropy_tests():
    rng=Random(2026112402);tests=0;largest_gap=0.
    from math import log
    for _ in range(150):
        u=rng.randrange(1,5);r=rng.randrange(2,5)
        na=rng.randrange(2,11);nb=rng.randrange(2,11)
        A=[tuple(rng.randrange(r) for _ in range(u)) for i in range(na)]
        B=[tuple(rng.randrange(r) for _ in range(u)) for i in range(nb)]
        C=Counter(tuple(zip(a,b)) for a,b in product(A,B))
        packets=materialize(C,u,u);Q=Counter(s for packet in packets for s in packet)
        assert len(packets)==min_packets(C,u,u)
        assert all(Q[s]>=n for s,n in C.items())
        for i in range(u):
            child=defaultdict(Counter)
            for s,n in Q.items():child[s[:i]][s[i]]+=n
            assert all(cc[x,y]==cc[y,x] for cc in child.values() for x,y in cc)
        N=na*nb;NQ=sum(Q.values());pad=NQ/N
        KL=sum(n/N*log2((n/N)/(Q[s]/NQ)) for s,n in C.items())
        assert -1e-10<=KL<=log2(pad)+1e-10
        ha=[H(Counter(s[i] for s in A)) for i in range(u)]
        # Test the direct entropy-loss bound before imposing balanced marginals.
        eps=[]
        for i in range(u):
            P=Counter(s[i] for s in A);part=defaultdict(Counter)
            for s in A:part[s[:i]][s[i]]+=1
            ep=sum(sum(abs(c[j]*na-P[j]*sum(c.values())) for j in range(r))
                   /(2*na*na) for c in part.values())
            eps.append(ep)
        assert sum(x*x for x in eps)<=8*log(2)*KL+1e-9
        info=sum(ha)-H(Counter(A));upper=sum(x*log2(r)+1 for x in eps)
        assert info<=upper+1e-10
        largest_gap=max(largest_gap,info)
        tests+=1
    return tests,largest_gap


def q_from_precision(M,L,p):
    h=[(2*P+(L>>i)-1)//(L>>i) for i,P in enumerate(p)]
    assert sum(h)<=M
    return h,[sum(h[i:]) for i in range(len(h))]


def exact_val(t,a,q):
    L=1<<a
    if t==0:return 0
    i=v2(t); qi=q[i] if i<len(q) else 0
    return (a+1-i)*qi+sum((q[v2(s)] if v2(s)<len(q) else 0)*(v2(t-s)-v2(s))
                           for s in range(-L,L) if s!=0 and s!=t)


def formula_val(t,a,q):
    if t==0:return 0
    L=1<<a; i=v2(t); qi=q[i] if i<len(q) else 0; D=L>>i
    return (a+1-i)*qi+sum((D>>j)*j*(qi-(q[i+j] if i+j<len(q) else 0))
                          for j in range(1,a-i+1))


def frac_v2(x):
    return v2(x.numerator)-v2(x.denominator) if x else 10**9


def local_poly(t,a,q):
    L=1<<a; coeff=[Fraction(1)]
    for s in range(-L,L):
        if s==0 or v2(s)>=len(q):continue
        aa=Fraction(s-t,s);bb=Fraction(-2*L,s)
        for _ in range(q[v2(s)]):
            out=[Fraction(0)]*(len(coeff)+1)
            for j,c in enumerate(coeff):
                out[j]+=aa*c;out[j+1]+=bb*c
            coeff=out
    return coeff


def precision_tests():
    rng=Random(2026112401); tested=expanded=coef_count=0
    for a in range(1,9):
        L=1<<a
        for _ in range(40):
            M=rng.randrange(max(8,a),max(20,4*a)+40)
            u=rng.randrange(1,a+1)
            h=[0]*u
            for _ in range(M):
                if rng.random()<0.75:h[rng.randrange(u)]+=1
            p=[(L>>i)*x//2 for i,x in enumerate(h)]
            hh,q=q_from_precision(M,L,p)
            deg=sum((L>>i)*x for i,x in enumerate(q))
            assert deg<=M*(2*L-1)
            for t in range(-L,L):
                V=exact_val(t,a,q)
                assert V==formula_val(t,a,q) and V>=0
                if t and v2(t)<u:assert V>=p[v2(t)]
                tested+=1
    for a in (1,2,3):
        L=1<<a
        for _ in range(12):
            M=4;u=rng.randrange(1,a+1);h=[0]*u
            for _ in range(M):h[rng.randrange(u)]+=1
            p=[(L>>i)*x//2 for i,x in enumerate(h)]
            hh,q=q_from_precision(M,L,p)
            for t in range(-L,L):
                coeff=local_poly(t,a,q)
                v=min(frac_v2(c) for c in coeff)
                assert v==formula_val(t,a,q)
                assert v>=0
                if t==0:assert coeff[0]==1
                elif v2(t)<u:assert v>=p[v2(t)]
                expanded+=1;coef_count+=len(coeff)
    # Positive taper: check loss of Hermite prefactor directly (no interpolation oracle).
    taper=0
    for a in range(2,9):
        L=1<<a
        for _ in range(20):
            u=rng.randrange(1,a+1)
            h=sorted([rng.randrange(1,20) for _ in range(u)],reverse=True)
            q=[sum(h[i:]) for i in range(u)]
            for t in range(-L,L):
                if not t or v2(t)>=u:continue
                i=v2(t);D=L>>i
                loss=sum(q[v2(s)]*(v2(t-s)-v2(s)) for s in range(-L,L)
                         if s and s!=t and v2(s)<u)
                formula=sum((D>>j)*j*(q[i]-(q[i+j] if i+j<u else 0)) for j in range(1,a-i+1))
                assert loss==formula and 0<=loss<=6*D*h[i]
                taper+=1
    return tested,expanded,coef_count,taper


# Matrix-packet diagnostic. Binary matrices represented by tuples of row words.
def matmul(X,Z,k):
    return tuple(xor_rows(Z,x) for x in X)

def xor_rows(Z,x):
    ans=0
    for i,row in enumerate(Z):
        if (x>>i)&1:ans^=row
    return ans

def madd(X,Y):return tuple(x^y for x,y in zip(X,Y))

def matrix_tests():
    zcount=packets=checks=ranktests=0
    for k,ell in ((1,2),(1,4),(2,2),(2,3),(2,4)):
        xs=list(product(range(1<<k),repeat=k));zs=list(product(range(1<<ell),repeat=k))
        I=tuple(1<<i for i in range(k)); goodz=[]
        Ts=[]
        for T in xs:
            v=next((v for v in range(1,1<<k) if xor_rows(madd(T,I),v)==0),0)
            if v:Ts.append((T,v))
        for Z in zs:
            if any(xor_rows(Z,v)==0 for v in range(1,1<<k)):continue
            zcount+=1;goodz.append(Z)
            H={matmul(x,Z,k) for x in xs};assert len(H)==len(xs)
            unseen=set(zs)
            while unseen:
                a0=next(iter(unseen));coset={madd(a0,h) for h in H}
                assert coset<=unseen
                unseen-=coset;packets+=1
            for T,v in Ts:
                for X in xs:
                    A=matmul(X,Z,k);B=madd(matmul(madd(X,T),Z,k),Z)
                    assert xor_rows(A,v)==xor_rows(B,v)
                    checks+=1
        # Union bound for non-full-row-rank differences, exactly checked.
        assert (len(zs)-len(goodz))*(1<<ell) <= ((1<<k)-1)*len(zs)
    # Verify the Cauchy residue crossing rank for k=1,2 using GF(2^(k*k)).
    for k,poly in ((1,0b11),(2,0b10011)):
        d=k*k;q=1<<d
        def mul(a,b):
            c=0
            while b:
                if b&1:c^=a
                b>>=1;a<<=1
                if a&q:a^=poly
            return c
        def pw(a,n):
            c=1
            while n:
                if n&1:c=mul(c,a)
                a=mul(a,a);n>>=1
            return c
        def inv(a):
            assert a
            z=pw(a,q-2);assert mul(a,z)==1
            return z
        def enc(X):return sum(x<<(i*k) for i,x in enumerate(X))
        xs=list(product(range(1<<k),repeat=k));I=tuple(1<<i for i in range(k))
        for T in xs:
            v=next((v for v in range(1,1<<k) if xor_rows(madd(T,I),v)==0),0)
            if not v:continue
            # A row group parity uses one column of v^T X; column label adds v^T I.
            R=[enc(X) for X in xs if xor_rows(X,v)&1==0]
            C=[enc(Y) for Y in xs if (xor_rows(Y,v)^v)&1==1]
            W=[[inv(enc(T)^x^y) for y in C] for x in R]
            n=len(W);assert n==q//2 and len(C)==n
            for j in range(n):
                pivot=next(i for i in range(j,n) if W[i][j]);W[j],W[pivot]=W[pivot],W[j]
                z=inv(W[j][j]);W[j]=[mul(a,z) for a in W[j]]
                for i in range(j+1,n):
                    z=W[i][j]
                    W[i]=[a^mul(z,b) for a,b in zip(W[i],W[j])]
            ranktests+=1
    return zcount,packets,checks,ranktests


if __name__=='__main__':
    if not __debug__:raise RuntimeError('assertions must be enabled')
    print('Locality/prefix-swap comparisons:',locality_tests(),flush=True)
    print('Independent label-law comparisons, symmetric cases:',*entropy_tests(),flush=True)
    n,pat,rows=packing_tests();print('Exact packing optimum comparisons:',n,'packet patterns:',pat,flush=True)
    print('Repeated-label padding ratios, r=3:',[(u,str(ratio)) for u,r,ratio in rows if r==3],flush=True)
    print('Hierarchical vs flat synthetic packet counts:',hierarchy_examples(),flush=True)
    print('Padded entropy regressions, maximum observed entropy defect:',*padding_entropy_tests(),flush=True)
    print('Precision valuations, explicit polynomial expansions, coefficients, upper taper checks:',*precision_tests(),flush=True)
    print('Matrix full-rank differences, coset packets, localization checks, crossing rank tests:',*matrix_tests(),flush=True)
    print('All assertions passed.',flush=True)
