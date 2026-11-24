"""Finite exact regressions for high-precision lifting and mode graph accounting.
Does not implement the asymptotic packet contractor.
Python standard library only.
"""
from collections import defaultdict
from itertools import product, combinations
from random import Random
from math import ceil, log2, pi, sqrt, log, prod

RNG = Random(20261124)

# F_2[x] polynomial helpers, then an unramified Z/(2^B) extension.
def rem2(a,b):
    while a and a.bit_length()>=b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a

def irreducible(h):
    for f in range((1<<h)|1, 1<<(h+1), 2):
        if all(rem2(f,g) for d in range(1,h//2+1)
               for g in range(1<<d,1<<(d+1))):
            return f
    raise AssertionError(h)

class Ring:
    def __init__(self,h,B):
        self.h=h; self.B=B; self.mod=1<<B
        f=irreducible(h); self.f=tuple((f>>i)&1 for i in range(h))
        self.zero=(0,)*h; self.one=(1,)+(0,)*(h-1)
    def scalar(self,n):return (n%self.mod,)+(0,)*(self.h-1)
    def bits(self,x):return tuple((x>>i)&1 for i in range(self.h))
    def add(self,a,b):return tuple((x+y)%self.mod for x,y in zip(a,b))
    def sub(self,a,b):return tuple((x-y)%self.mod for x,y in zip(a,b))
    def mul(self,a,b):
        h=self.h; mod=self.mod; c=[0]*(2*h-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):
                    if y:c[i+j]+=x*y
        for j in range(2*h-2,h-1,-1):
            if c[j]:
                for k,f in enumerate(self.f):
                    if f:c[j-h+k]-=c[j]
        return tuple(x%mod for x in c[:h])
    def inv(self,a):
        # Compute inverse mod 2 by exhaustive field search (small fixtures), then Newton.
        aa=tuple(x&1 for x in a)
        assert any(aa)
        for j in range(1,1<<self.h):
            z=self.bits(j)
            if tuple(x&1 for x in self.mul(aa,z))==self.one:break
        else:raise AssertionError('no inverse')
        for _ in range((self.B-1).bit_length()):
            z=self.mul(z,self.sub(self.scalar(2),self.mul(a,z)))
        assert self.mul(a,z)==self.one
        return z
    def eval(self,coeffs,x):
        z=self.zero
        for c in reversed(coeffs):z=self.add(self.mul(z,x),c)
        return z
    def interpolate(self,xs,ys):
        # Incremental Newton interpolation using only unit inverses.
        c=[self.zero]; basis=[self.one]
        for x,y in zip(xs,ys):
            b=self.eval(basis,x); assert any(v&1 for v in b)
            a=self.mul(self.sub(y,self.eval(c,x)),self.inv(b))
            if len(c)<len(basis):c += [self.zero]*(len(basis)-len(c))
            for i,d in enumerate(basis):c[i]=self.add(c[i],self.mul(a,d))
            nb=[self.zero]*(len(basis)+1)
            for i,d in enumerate(basis):
                nb[i]=self.sub(nb[i],self.mul(x,d)); nb[i+1]=self.add(nb[i+1],d)
            basis=nb
        return c

# Ring of original outputs: Z/(2^P)[U,V]/(U^e+... , V^ell).
# For fixture e=2 use U^2+U+1; for e=1 U is absent.
def old_mul(a,b,P,e,ell):
    c=[[0]*(2*e-1) for _ in range(ell)]
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            if j+l<ell:c[j+l][i+k]+=x*y
    if e==2:
        for j in range(ell):
            c[j][0]-=c[j][2]; c[j][1]-=c[j][2]
    return {(i,j):c[j][i]%(1<<P) for i in range(e) for j in range(ell)}

def old_add(a,b,P):return {k:(a.get(k,0)+b.get(k,0))%(1<<P) for k in a.keys()|b.keys()}

def assignments(edges):return product((0,1),repeat=len(edges))

def key(x,edges,v):return tuple(x[i] for i,uv in enumerate(edges) if v in uv)

def poly_for(a,P,e,ell,baseT,baseU):
    # T is the binary digit variable; U,V are raw (not yet reduced).
    coeff={}
    for (i,j),x in a.items():
        for k in range(P):
            if (x>>k)&1:coeff[k+baseT*(i+baseU*j)]=1
    return coeff

def precision_case(edges,P,e,ell):
    V=1+max(max(uv) for uv in edges); S=len(edges)
    tables=[]
    for v in range(V):
        deg=sum(v in uv for uv in edges)
        tables.append({bits:{(i,j):RNG.randrange(1<<P) for i in range(e) for j in range(ell)}
                       for bits in product((0,1),repeat=deg)})
    direct={(i,j):0 for i in range(e) for j in range(ell)}
    for x in assignments(edges):
        z={(0,0):1}
        for v in range(V):z=old_mul(z,tables[v][key(x,edges,v)],P,e,ell)
        direct=old_add(direct,z,P)
    baseT=V*(P-1)+1; baseU=V*(e-1)+1; baseV=V*(ell-1)+1
    D=baseT*baseU*baseV-1
    # Strict height bound: each coefficient <= 2^S*(P*e*ell)^V.
    B=( (1<<S)*(P*e*ell)**V ).bit_length()
    R=Ring(max(1,D.bit_length()),B)
    enc=[{bits:poly_for(a,P,e,ell,baseT,baseU) for bits,a in table.items()} for table in tables]
    # Build an independent exact integer coefficient vector by sparse convolution.
    want=defaultdict(int)
    for x in assignments(edges):
        z={0:1}
        for v in range(V):
            nz=defaultdict(int)
            for i,a in z.items():
                for j,b in enc[v][key(x,edges,v)].items():nz[i+j]+=a*b
            z=nz
        for i,a in z.items():want[i]+=a
    assert max(want.values(),default=0)<(1<<B)
    xs=[R.bits(i) for i in range(D+1)]
    ys=[]
    for q in xs:
        powers=[R.one]
        for i in range(D):powers.append(R.mul(powers[-1],q))
        vals=[]
        for table in enc:
            tv={}
            for bits,p in table.items():
                z=R.zero
                for i in p:z=R.add(z,powers[i])
                tv[bits]=z
            vals.append(tv)
        ans=R.zero
        for x in assignments(edges):
            z=R.one
            for v in range(V):z=R.mul(z,vals[v][key(x,edges,v)])
            ans=R.add(ans,z)
        # This comparison also catches any bad extension-ring implementation.
        literal=R.zero
        for i,c in want.items():literal=R.add(literal,R.mul(R.scalar(c),powers[i]))
        assert ans==literal
        ys.append(ans)
    coeff=R.interpolate(xs,ys)
    assert all(c==R.scalar(want[i]) for i,c in enumerate(coeff))
    lifted={(i,j):0 for i in range(e) for j in range(ell)}
    for idx,c in enumerate(coeff):
        k=idx%baseT; rem=idx//baseT; i=rem%baseU; j=rem//baseU
        if j>=ell:continue
        # Independently reduce U^i by repeated multiplication.
        u={(0,0):1}
        for _ in range(i):u=old_mul(u,{(1,0):1},P,e,ell)
        for (ii,jj),a in u.items():
            if j+jj<ell:lifted[ii,j+jj]=(lifted[ii,j+jj]+a*c[0]*(1<<k))%(1<<P)
    assert lifted==direct
    print(f'precision: edges={S}, tensors={V}, P={P}, e={e}, ell={ell}, digit precision={B}, degree={D}, exact coefficients={D+1}',flush=True)
    return D+1

def mode_counts(c,J,P):
    out=[1]+[0]*P
    for j in range(1,J+1):
        for _ in range(c):
            for h in range(j,P+1):out[h]+=out[h-j]
    return out

def entropy_tests():
    tests=0; K=pi*sqrt(2/3)/log(2)
    for c in range(1,8):
        for J in (3,5,9):
            N=mode_counts(c,J,120)
            for h,n in enumerate(N):
                if n:assert log2(n)<=K*sqrt(c*h)+1e-9
                tests+=1
    # Pair budgets: bound the total width and every degree by Cauchy--Schwarz.
    cases=0
    for _ in range(300):
        r=RNG.randrange(3,10); es=list(combinations(range(r),2))
        caps=[RNG.randrange(1,8) for e in es]; h=[RNG.randrange(0,35) for e in es]
        bits=[(mode_counts(c,9,p)[p]-1).bit_length() for c,p in zip(caps,h)]
        assert sum(bits)<=K*sqrt(sum(caps)*sum(h))+len(es)+1e-8
        for v in range(r):
            I=[i for i,e in enumerate(es) if v in e]
            assert sum(bits[i] for i in I)<=K*sqrt(sum(caps[i] for i in I)*sum(h[i] for i in I))+len(I)+1e-8
            cases+=1
    # Dense-budget witness: all pair alphabets simultaneously have large size.
    witnesses=[]
    for r,c,J in [(4,2,10),(5,3,12),(6,4,16),(8,5,18)]:
        # Subset (0/1 exponent) distribution by exact weight.
        maxh=c*J*(J+1)//2; z=[1]+[0]*maxh
        for j in range(1,J+1):
            for _ in range(c):
                for h in range(maxh,j-1,-1):z[h]+=z[h-j]
        assert sum(z)==1<<(c*J)
        h=max(range(len(z)),key=z.__getitem__)
        assert z[h]*(maxh+1)>=1<<(c*J)
        n=mode_counts(c,J,h)[h]; assert n>=z[h]
        l=(n-1).bit_length(); S=r*(r-1)//2*l
        P=1+r*(r-1)//2*maxh
        assert r*(r-1)//2*h<P
        # Every balanced cut in a uniform bundle-complete graph is large.
        for a in range(ceil(r/3),2*r//3+1):assert a*(r-a)*l>=S*4/9
        witnesses.append((r,c,J,l,S,h,P))
    print(f'entropy: {tests:,} partition-count inequalities; {cases:,} incident-degree inequalities; {len(witnesses)} dense-budget witnesses')
    for t in witnesses:print('  r,c,J,bundle_bits,total_bits,budget_per_pair,P =',t)

def v2(n):return (n&-n).bit_length()-1

def antipodal_tests():
    cyclic=0; linear=0; cuts=0
    for L in (2,4,8,16):
        mask=(1<<L)-1
        for bits in range(1<<L):
            good=[t for t in range(1,L) if (((bits<<t)|(bits>>(L-t)))&mask)==(bits^mask)]
            assert len({v2(t) for t in good})<=1
            cyclic+=1
    from math import gcd
    for L in range(2,13):
        for bits in range(1<<L):
            good=[t for t in range(1,L) if ((bits^(bits>>t))&((1<<(L-t))-1))==((1<<(L-t))-1)]
            for t,s in combinations(good,2):
                if v2(t)!=v2(s):assert L<t+s-gcd(t,s)+1
            linear+=1
    for r in range(2,9):
        full=(1<<r)-1
        # Complementary cuts are identical; fix the first vertex on one side.
        masks=[a for a in range(1,full) if a&1]
        for a,b in combinations(masks,2):
            assert any(((a>>i)^(a>>j))&1 and ((b>>i)^(b>>j))&1 for i,j in combinations(range(r),2))
            cuts+=1
    print(f'pair geometry: {cyclic:,} cyclic colorings; {linear:,} linear colorings; {cuts:,} intersecting complete-graph cut pairs')

if __name__=='__main__':
    if not __debug__:raise RuntimeError('assertions must be enabled')
    total=0
    for args in [([(0,1)],64,1,1), ([(0,1),(0,1)],48,1,1),
                 ([(0,1),(1,2),(2,0)],18,1,1), ([(0,1),(1,2)],12,1,1),
                 ([(0,1)],4,2,2), ([(0,1),(1,2)],3,2,2),
                 ([(0,1)],8,2,1), ([(0,1)],6,1,2)]:
        total+=precision_case(*args)
    print(f'precision total: 8 networks; {total:,} exact evaluations and interpolated coefficients')
    entropy_tests()
    antipodal_tests()
