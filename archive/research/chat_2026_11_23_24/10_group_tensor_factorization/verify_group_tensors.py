"""Exact component tests for group tensor splitting, weighted packet moments,
and unit-separated XOR Cauchy blocks. Standard library only."""
from fractions import Fraction
from itertools import product
from math import factorial
from random import Random


def v2(x):
    if not x:
        return 10**9
    x = abs(x)
    return (x & -x).bit_length()-1


def poly_mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c


def sparse_mul(a,b):
    out={}
    for u,x in a.items():
        for v,y in b.items():
            w=tuple(i+j for i,j in zip(u,v)); out[w]=out.get(w,0)+x*y
    return {w:x for w,x in out.items() if x}


def dot(a,b): return sum(x*y for x,y in zip(a,b))


def moment(j,bits): return sum((i**j)*x for i,x in enumerate(bits))


def eval_mon(beta,bits):
    out=1
    for j,e in enumerate(beta): out*=moment(j,bits)**e
    return out


def indicator_expansion(x):
    # Lagrange polynomials at the integer locations 0,...,k-1.
    k=len(x); d=factorial(k-1); zero=(0,)*k
    out={zero:1}
    for t,xt in enumerate(x):
        p=[1]; den=1
        for j in range(k):
            if j!=t:
                p=poly_mul(p,[-j,1]); den*=t-j
        assert d%den==0
        p=[v*(d//den) for v in p]
        factor={zero:0 if xt else d}
        for j,v in enumerate(p):
            e=tuple(int(h==j) for h in range(k)); factor[e]=v if xt else -v
        out=sparse_mul(out,factor)
    return out,d**k


def test_moments():
    n_eval=n_fixed=n_gram=0
    for k in range(1,7):
        xs=list(product((0,1),repeat=k)); exps=[]; c=None
        H=k*k
        for x in xs:
            e,c=indicator_expansion(x); exps.append(e)
            assert all(sum((j+1)*v for j,v in enumerate(beta))<=H and sum(beta)<=k for beta in e)
        betas=sorted(set().union(*[set(e) for e in exps]))
        vals={beta:[eval_mon(beta,x) for x in xs] for beta in betas}
        # Include exact-weight sentinel padding: two independent fixed-one moments of weights k,1.
        for beta in betas:
            w=sum((j+1)*v for j,v in enumerate(beta)); q,r=divmod(H-w,k)
            assert w+k*q+r==H and sum(beta)+q+r<=3*k
            for x,v in zip(xs,vals[beta]):
                assert v==eval_mon(beta,x)*(1**q)*(1**r); n_fixed+=1
        reconstructed=[]
        for i,e in enumerate(exps):
            row=[sum(a*vals[beta][j] for beta,a in e.items()) for j in range(len(xs))]
            for j,val in enumerate(row):
                assert val==c*(i==j); n_eval+=1
            reconstructed.append(row)
        # Gram matrix transformed on both sides: A (V V^t) A^t = c^2 I.
        # Use direct sums over the copied source bits, without materializing a huge raw Gram matrix.
        precision=2*v2(c)+2; modulus=1<<precision
        for i,j in product(range(len(xs)),repeat=2):
            assert dot(reconstructed[i],reconstructed[j])%modulus==(c*c*(i==j))%modulus; n_gram+=1
        print(f'k={k}: {len(xs)} source states, {len(betas)} used mode monomials, fixed weight {H}, precision {precision}')
    print(f'Moment indicator identities: {n_eval}; fixed-budget comparisons: {n_fixed}; transformed Gram comparisons: {n_gram}')


# Binary extension fields F_{2^d}, with known irreducible polynomials for the test degrees.
class GF:
    def __init__(self,d,f): self.d=d; self.f=f; self.q=1<<d
    def mul(self,a,b):
        v=0
        while b:
            if b&1: v^=a
            b>>=1; a<<=1
            if a&self.q: a^=self.f
        return v
    def pow(self,a,n):
        v=1
        while n:
            if n&1: v=self.mul(v,a)
            a=self.mul(a,a); n>>=1
        return v
    def inv(self,a):
        assert a; return self.pow(a,self.q-2)
    def rank(self,mat):
        rows=[row[:] for row in mat]; k=0
        for j in range(len(rows[0])):
            i=next((i for i in range(k,len(rows)) if rows[i][j]),None)
            if i is None: continue
            rows[i],rows[k]=rows[k],rows[i]
            c=self.inv(rows[k][j]); rows[k]=[self.mul(x,c) for x in rows[k]]
            for i in range(k+1,len(rows)):
                c=rows[i][j]
                if c: rows[i]=[x^self.mul(c,y) for x,y in zip(rows[i],rows[k])]
            k+=1
            if k==len(rows): break
        return k


def test_xor():
    count=entries=0
    for d,f in [(2,0b111),(3,0b1011),(4,0b10011),(5,0b100101),(6,0b1000011)]:
        field=GF(d,f); q=field.q
        assert all(field.mul(x,field.inv(x))==1 for x in range(1,q))
        # Every coordinate subspace and every disjoint affine coset.
        for mask in range(1,q):
            H=[x for x in range(q) if x&mask==0]
            # Representatives supported on the fixed coordinates, nonzero means disjoint.
            for t in range(1,q):
                if t&~mask or not t: continue
                inv=[field.inv(t^h) for h in H]; total=0
                for x in inv: total^=x
                assert total
                # Matrix coefficient after dividing the crossing derivative by common uniformizer.
                mat=[[field.inv(t^x^y) for y in H] for x in H]
                assert field.rank(mat)==len(H); count+=1; entries+=len(H)**2
                # Derivative/product identity of the additive root polynomial.
                p=1; der=1
                for h in H:
                    p=field.mul(p,t^h)
                    if h: der=field.mul(der,h)
                assert total==field.mul(der,field.inv(p))
        print(f'F_{{2^{d}}}: all tested off-coset inverse blocks nonsingular')
    print(f'XOR crossing blocks: {count}; field matrix entries: {entries}')


def test_symbolic():
    rng=Random(20261124); checks=0; noncanonical=0
    # For each parent index, use signed reconstruction sums with denominators powers of 2.
    # Independent local residues can be lifted arbitrarily before distributing products.
    for _ in range(200):
        P=rng.randrange(2,10); r=rng.randrange(2,5); mod=1<<P
        losses=[rng.randrange(1,7) for _ in range(r)]; B=sum(losses); big=1<<(P+B)
        expected=got=0
        for beta in range(rng.randrange(1,6)):
            xs=[rng.randrange(mod) for _ in range(r)]; expected=(expected+__import__('math').prod(xs))%mod
            local=[]
            for a,x in zip(losses,xs):
                small=1<<(P+a); terms=[rng.randrange(-small,small) for _ in range(3)]
                terms.append(((1<<a)*x-sum(terms))%small)
                # Noncanonical lifts simulate child-network products at their creation precision.
                lifted=[v+small*rng.randrange(-3,4) for v in terms]
                noncanonical+=sum(v!=w for v,w in zip(terms,lifted)); local.append(lifted)
            for terms in product(*local): got=(got+__import__('math').prod(terms))%big; checks+=1
        assert got%(1<<B)==0 and (got>>B)==expected
    print(f'Symbolic scaled-join summands: {checks}; changed entry lifts: {noncanonical}; 200 exact joins')


if __name__=='__main__':
    test_moments(); test_xor(); test_symbolic(); print('All exact checks passed.')
