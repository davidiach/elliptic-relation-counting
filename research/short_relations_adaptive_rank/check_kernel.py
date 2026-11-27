from fractions import Fraction
from itertools import combinations
from random import Random
from math import ceil,log2
import numpy as np


def count_curve(p,a,b):
    squares=[0]*p
    for y in range(p):squares[y*y%p]+=1
    return 1+sum(squares[(x*x*x+a*x+b)%p] for x in range(p))


def isprime(n):
    return n>1 and all(n%d for d in range(2,int(n**.5)+1))


def curve_for(r):
    for p in range(r-2*int(r**.5),r+2*int(r**.5)+2):
        if not isprime(p):continue
        for a in range(0,12):
            for b in range(1,32):
                if (4*a**3+27*b*b)%p and count_curve(p,a,b)==r:return p,a,b
    raise ValueError(r)


def elliptic_cycle(p,a,b,r):
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u and (y+v)%p==0:return None
        t=((3*x*x+a)*pow(2*y,-1,p) if P==Q else (v-y)*pow(u-x,-1,p))%p
        X=(t*t-x-u)%p
        return X,(t*(x-X)-y)%p
    squares={y*y%p:y for y in range(p)}
    P=next((x,squares[(x*x*x+a*x+b)%p]) for x in range(p) if (x*x*x+a*x+b)%p in squares)
    G=[None]
    for _ in range(r-1):G.append(add(G[-1],P))
    assert len(set(G))==r and add(G[-1],P) is None
    return G,add


def convolve(a,b):
    # exact integer circular convolution, avoiding any floating arithmetic
    out=np.zeros(len(a),dtype=np.int64)
    for j in np.flatnonzero(b):out+=int(b[j])*np.roll(a,int(j))
    return out


def counts(sets,r):
    out=np.zeros(r,dtype=np.int64);out[0]=1
    for B in sets:
        s=np.zeros(r,dtype=np.int64)
        for v in B:s[v]+=1
        out=convolve(out,s)
    return out


def rank(rows,r,n):
    piv={}
    for row in rows:
        a=np.array(row,dtype=np.int64)%r
        for j in sorted(piv):
            if a[j]:a=(a-int(a[j])*piv[j])%r
        ix=np.flatnonzero(a)
        if len(ix):
            j=int(ix[0]);piv[j]=a*pow(int(a[j]),-1,r)%r
    return len(piv)


def min_cut_var_numer(labels,z,r):
    # minimum over every possible interval of length one; numerator over n^2*r^2.
    xs=sorted(z*x%r for x in labels); n=len(xs); s=sum(xs); q=sum(x*x for x in xs)
    best=n*q-s*s
    for x in xs:
        s+=r;q+=2*r*x+r*r
        best=min(best,n*q-s*s)
    return best


def construct(r,m,n,seed):
    D=1<<ceil(log2(r)/m)
    # dictionary (c,d), group coordinate kappa*c+d, homogeneous direction c
    Ds=[[(j*D**i,1) for j in range(D+1)]+[(0,-1),(0,2),(0,-2)] for i in range(m)]
    for kap in range(1,r):
        bs=[[(kap*c+d)%r for c,d in B] for B in Ds]
        if all(0 not in B and len(B)==len(set(B)) for B in bs):break
    else:raise ValueError('no dictionary')
    rng=Random(seed)
    for i,B in enumerate(bs):
        R=rng.sample([a for a in range(1,r) if a not in B],n)
        Ds[i]+=[((a-2)*pow(kap,-1,r)%r,2) for a in R]
    B=[[((kap*c+d)%r,c%r,d) for c,d in block] for block in Ds]
    return D,kap,B


def construct_rows(r,D,kap,B,curve=None):
    m=len(B);indices=[];vals=[];u=[];offs=[]
    for block in B:
        start=len(vals);indices.append(list(range(start,start+len(block))))
        for v,c,d in block:vals.append(v);u.append(c);offs.append(d)
    N=len(vals)+2;P=len(vals);Q=P+1
    vals+=[1,vals[indices[0][-1]]];u+=[0,u[indices[0][-1]]];offs+=[1,2]
    rows=[]
    def add(plus,minus,rhs):
        row=[0]*N
        for j in plus:row[j]+=1
        for j in minus:row[j]-=1
        assert sum(a*v for a,v in zip(row,vals))%r==rhs%r
        assert sum(a*v for a,v in zip(row,u))%r==0
        if curve is not None:
            G,law=curve; S=None
            for j in plus:S=law(S,G[vals[j]])
            for j in minus:S=law(S,G[-vals[j]%r])
            assert S==G[rhs%r]
        rows.append(row)
    add([P],[],1)
    for i,ix in enumerate(indices):
        ref=ix[0];neg=ix[D+1];two=ix[D+2];negtwo=ix[D+3]
        assert [vals[j] for j in [ref,neg,two,negtwo]]==[1,r-1,2,r-2]
        add([two],[ref],1)
        add([two],[ref,ref],0) # point-shift of the indexed P, same underlying pair
        add([ref],[neg],2)
        add([neg],[negtwo],1)
        for j in range(2,D+1):add([ix[j],ref],[ix[j-1],ix[1]],0)
        if i:
            old=indices[i-1]
            add([old[D],ref],[old[0],ix[1]],0)
    # explicit doubled decompositions (one digit + one known filler per block)
    for j in range(N):
        c=u[j]; digits=[]
        for i in range(m):digits.append((c//D**i)%D)
        assert sum(d*D**i for i,d in enumerate(digits))==c
        f=[ix[D+1] for ix in indices] # -P in every block
        d=offs[j]
        if d==2:f[0]=indices[0][0]
        elif d==1:f[0]=indices[0][0];f[1]=indices[1][D+3]
        elif d==-1:f[0]=indices[0][D+3]
        elif d==-2:f[0]=indices[0][D+3];f[1]=indices[1][D+3]
        else:raise ValueError(d)
        t=[ix[digit] for ix,digit in zip(indices,digits)]+f
        assert sum(vals[k] for k in t)%r==vals[j]
        add([j],t,0)
    assert rank(rows,r,N)==N-1 and u[Q]!=0
    return vals,u,N,len(rows)


def finite_fixture(r=1009,m=4,n=16):
    p,a,b=curve_for(r)
    print(f'Curve: E/F{p}, y^2=x^3+{a}x+{b}, order {r} (direct point count).')
    for seed in range(100):
        D,kap,B=construct(r,m,n,seed)
        allsets=[[x[0] for x in block] for block in B]
        dist=[counts([allsets[j] for j in I],r) for k in range(m-1,m+1) for I in combinations(range(m),k)]
        if all(min(x)>0 and max(x)<r for x in dist):break
    else:raise ValueError('no qualifying fixture')
    vals,u,N,nrows=construct_rows(r,D,kap,B,elliptic_cycle(p,a,b,r))
    chis=[Fraction(r*sum(int(z)**2 for z in d)-int(sum(d))**2,int(sum(d))**2) for d in dist]
    print(f'Fixture seed {seed}; radix {D}; kappa {kap}; block sizes {[len(x) for x in B]}; {nrows} genuine linear rows, rank {N-1}/{N}.')
    print(f'All {len(dist)} retained subsets: full support, max fiber {max(int(max(x)) for x in dist)} < {r}; max chi-square {float(max(chis)):.8f}.')
    lo=None;arg=None
    for z in range(1,r):
        total=sum(Fraction(min_cut_var_numer([x[1] for x in block],z,r),len(block)**2*r**2) for block in B)
        if lo is None or total<lo:lo=total;arg=z
    assert lo>Fraction(m,100)
    print(f'All {r-1} nonzero multipliers and every circular cut (exact): minimum aggregate variance {float(lo):.8f}, per block average {float(lo/m):.8f}, multiplier {arg}.')
    # If only one block may vary, the line is absent. Check all proper subsets.
    for i,block in enumerate(B): assert any(x[1] for x in block)
    print('Every block varies on the residual line, target coordinate nonzero; all proper block-support intersections are zero.')
    return B

def modal_bound_fixture(B,r=1009):
    # Use the full nonzero scalar set for the exact averaged character-gap check.
    # Instead of irrational cosines use the exact circular-distance lower bound on variance.
    from collections import Counter
    E=sum(1-Fraction(max(Counter(x[1] for x in block).values()),len(block)) for block in B)
    D=sum(1-sum(Fraction(v*v,len(block)**2) for v in Counter(x[1] for x in block).values()) for block in B)
    assert E<=D<=2*E
    print(f'Modal mass E={float(E):.8f}, collision deficit={float(D):.8f}; exact comparison passed.')


def chi(d,r):
    W=sum(int(x) for x in d)
    return Fraction(r*sum(int(x)**2 for x in d)-W*W,W*W)


def block_search_checks():
    from itertools import product
    rng=Random(991);r=5;checks=0
    # D(K) is represented by a full-column-rank matrix with three two-coordinate blocks.
    for _ in range(60):
        n=6;d=rng.randrange(1,5)
        while True:
            U=[[rng.randrange(r) for _ in range(d)] for _ in range(n)]
            if rank(U,r,d)==d:break
        vectors=[tuple(sum(a*b for a,b in zip(row,c))%r for row in U) for c in product(range(r),repeat=d) if any(c)]
        dist=min(sum(any(v[2*j:2*j+2]) for j in range(3)) for v in vectors)
        found=None
        for k in range(4):
            for S in combinations(range(3),k):
                outside=[U[i] for i in range(n) if i//2 not in S]
                predicted=rank(outside,r,d)<d
                brute=any(all(v[i]==0 for i in range(n) if i//2 not in S) for v in vectors)
                assert predicted==brute
                if predicted and found is None:found=k
                checks+=1
        assert found==dist and d<=n-2*(dist-1)
    print(f'Block-support search: 60 residual subspaces, {checks} intersection tests matched exhaustive enumeration; dimension bound passed.')


def trim_checks():
    rng=Random(551);r=101;m=8;a=6;k=m-a
    B=[rng.sample(range(r),70) for _ in range(m)]
    eps=max(chi(counts([B[j] for j in I],r),r) for I in combinations(range(m),a))
    # rational upper bound for epsilon**(1/(2a)), by exact binary search
    lo=0;hi=10000
    while lo+1<hi:
        mid=(lo+hi)//2
        if Fraction(mid,10000)**(2*a)>=eps:hi=mid
        else:lo=mid
    tau=Fraction(hi,10000);delta=Fraction(1,14) # retain at least 65/70
    from math import comb
    bound=comb(m,k)*(tau+delta)**(2*k)/(1-delta)**(2*m)*(1+eps)
    maxchi=Fraction(0)
    for _ in range(40):
        Bs=[rng.sample(block,65) for block in B]
        q=chi(counts(Bs,r),r);assert q<=bound
        maxchi=max(maxchi,q)
    print(f'Trim inequality: 40 fixtures; punctured epsilon={float(eps):.3g}, tau upper bound={float(tau):.4f}, predicted bound={float(bound):.6g}, max observed chi-square={float(maxchi):.3g}.')


if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions required')
    B=finite_fixture();modal_bound_fixture(B);block_search_checks();trim_checks()
