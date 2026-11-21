"""Exact integer seed covariance and numerical singular spectra under the full hash law."""
from fractions import Fraction
from itertools import product
from math import prod

import numpy as np

from fixtures import curve_cycle, hash_families


def slot_gram(curve,M,w):
    masks,weight = hash_families(curve,M,w)
    N = len(curve)
    indices = np.arange(1,(N+1)//2)
    xs = np.array([curve[j][0] for j in indices])
    mask = np.array([[j in body for j in indices] for body in masks],dtype=np.int64)
    masses = np.array(weight,dtype=np.int64)
    seed_count = int(masses.sum())
    # The integer identity includes diagonal entries and all nonempty residue classes.
    second = mask.T@(masses[:,None]*mask)
    blocks = (xs[:,None]%w == xs[None,:]%w).astype(np.int64)
    identity = np.eye(len(xs),dtype=np.int64)
    assert np.array_equal(second*M*M,seed_count*(M*w*identity+w*w*(1-blocks)))
    char = 2*np.cos(2*np.pi*np.outer(indices,indices)/N)
    rho = w/M
    normalized = char.T@second@char/(seed_count*N*rho)
    formula = identity-(2-rho)*np.ones_like(identity)/N-rho*(char.T@blocks@char)/N
    assert np.allclose(normalized,formula)
    deficit = np.linalg.eigvalsh(identity-normalized)
    assert deficit.min()>-1e-9 and deficit.max()<1+1e-9
    assert sum(deficit>1e-8)<=w+1
    return normalized,masks,masses,mask,char


def square_root(matrix):
    e,v = np.linalg.eigh(matrix)
    assert e.min()>-1e-8
    return (v*np.sqrt(e.clip(0)))@v.T


def study(p,a,widths,cut,direct=False):
    curve = curve_cycle(p,a)
    N = len(curve)
    M = 1<<p.bit_length()
    slots = [slot_gram(curve,M,w) for w in widths]
    grams = [s[0] for s in slots]
    left = np.multiply.reduce(grams[:cut])
    right = np.multiply.reduce(grams[cut:])
    lam = N**(len(widths)-1)*prod(Fraction(w,M) for w in widths)
    root_left,root_right = square_root(left),square_root(right)
    modes = np.arange(1,(N+1)//2)
    for t in (1,2,5):
        cos = 2*np.cos(2*np.pi*modes*t/N)
        actual = np.linalg.svd((root_left*cos)@root_right,compute_uv=False)**2
        cap = np.sort(cos**2)[::-1]
        assert np.all(actual <= cap+1e-8)
        assert -1e-8 < (cap-actual).sum() <= 4*sum(w+1 for w in widths)+1e-8
        if direct:
            # Independent witness counts by literal addition of signed labels.
            assert len(widths)==3 and cut==1
            N1,N2,N3 = [len(s[1]) for s in slots]
            counts = np.empty((N1,N2,N3))
            sizes = [np.array(list(map(len,s[1]))) for s in slots]
            thirds = [set(body) for body in slots[2][1]]
            for i,A in enumerate(slots[0][1]):
                for j,B in enumerate(slots[1][1]):
                    needed = [(t-u-v)%N for u,v in product(A,B)]
                    counts[i,j] = [sum(z in S for z in needed) for S in thirds]
            K = counts-np.einsum("i,j,k->ijk",*sizes)/N
            for axis,s in enumerate(slots):
                shape = [1,1,1]
                shape[axis] = len(s[2])
                K *= np.sqrt(s[2]/s[2].sum()).reshape(shape)
            svals = np.linalg.svd(K.reshape(N1,-1),compute_uv=False)**2
            assert np.allclose(svals[:len(actual)],actual*float(lam)/N)
            assert np.allclose(svals[len(actual):],0)
        if t==1:
            R = len(cos)//2
            print(f"N={N}, widths={widths}, cut={cut}, modes={len(cos)}, "
                  f"best residual of nonconstant count at rank {R} / lambda = {actual[R:].sum()/N:.6f}, "
                  f"full nonconstant energy / lambda = {actual.sum()/N:.6f}")


def cofactor_check(p):
    """Full rational curve y²=x³-x, with group C_n x C_2 and prime target order n/2."""
    points = [(x,y) for x,y in product(range(p),repeat=2) if (y*y-x**3+x)%p==0]
    N = len(points)+1
    n = N//2
    def add(P,Q):
        if P is None or Q is None:
            return Q if P is None else P
        x,y=P
        z,w=Q
        if x==z and (y+w)%p==0:
            return None
        num,den = (3*x*x-1,2*y) if P==Q else (w-y,z-x)
        slope=num*pow(den,-1,p)%p
        X=(slope*slope-x-z)%p
        return X,(slope*(x-X)-y)%p
    for G in points:
        orbit=[None,G]
        while (P:=add(orbit[-1],G)) is not None:
            orbit.append(P)
        if len(orbit)==n:
            break
    J=next(P for P in points if P[1]==0 and P not in orbit)
    group = [add(P,H) for H in (None,J) for P in orbit]
    assert set(group)==set(points+[None]) and len(group)==N
    modes = [(l,k) for l in range(1,n//2) for k in (0,1)]
    real = [(l,k) for l in (0,n//2) for k in (0,1)]
    xs=sorted({x for x,y in points})
    fibers=[[i for i,P in enumerate(group) if P and P[0]==x] for x in xs]
    def evaluations(characters):
        return np.array([[sum(np.cos(2*np.pi*l*(i%n)/n)*(-1)**(k*(i//n))
                              for i in fiber) for l,k in characters] for fiber in fibers])
    F=evaluations(modes)
    single=F[[len(f)==1 for f in fibers]]
    one=np.ones((len(modes),len(modes)))
    assert np.allclose(F.T@F,N*np.eye(len(modes))-2*one-single.T@single)
    M=1<<p.bit_length()
    widths=(2,2,4)
    grams=[]
    real_features=[]
    configs=[]
    for w in widths:
        masks,mass=hash_families(group,M,w)
        weights=np.array(mass)/sum(mass)
        Z=np.array([[fiber[0] in S for fiber in fibers] for S in masks],dtype=float)
        rho=w/M
        Q=(Z@F).T@(weights[:,None]*(Z@F))/(N*rho)
        I=np.eye(len(modes))
        eig=np.linalg.eigvalsh(I-Q)
        assert eig.min()>-1e-8 and eig.max()<1+1e-8
        assert eig.sum() <= (2+len(single))*len(modes)/N+rho*len(xs)+1e-8
        grams.append(Q)
        real_features.append(Z@evaluations(real))
        configs.append((masks,weights))
    t=2
    d=2*np.cos(2*np.pi*np.array([l for l,k in modes])*t/n)
    normalized=(square_root(grams[0])*d)@square_root(grams[1]*grams[2])
    cap=np.sort(d*d)[::-1]
    svals=np.linalg.svd(normalized,compute_uv=False)**2
    assert np.all(svals<=cap+1e-8)
    if p==19:
        mats=[len(m[0]) for m in configs]
        H=np.empty(mats)
        third=[set(body) for body in configs[2][0]]
        for i,A in enumerate(configs[0][0]):
            for j,B in enumerate(configs[1][0]):
                required=[(t-u%n-v%n)%n+n*((u//n+v//n)%2) for u,v in product(A,B)]
                H[i,j]=[sum(z in body for z in required) for body in third]
        factors=np.array([np.cos(2*np.pi*l*t/n) for l,k in real])
        residual=H-np.einsum("il,jl,kl,l->ijk",*real_features,factors)/N
        for axis,(_,weights) in enumerate(configs):
            shape=[1,1,1]
            shape[axis]=len(weights)
            residual*=np.sqrt(weights).reshape(shape)
        literal=np.linalg.svd(residual.reshape(mats[0],-1),compute_uv=False)**2
        lam=N**2*prod(widths)/M**3
        assert np.allclose(literal[:len(svals)],svals*lam/N)
        assert np.allclose(literal[len(svals):],0)
    print(f"Noncyclic group order {N}, target order {n//2}, "
          f"{len(real)} real modes, {len(modes)} opposite pairs verified")


def margin_check():
    xs=np.linspace(-1,1,20001)
    for gap in (0.5,0.2,0.1,0.04):
        degree=int(np.ceil(8*np.pi/gap))
        k=np.arange(1,degree)
        coeff=np.r_[0.5,2/np.pi*(1-k/degree)*np.sin(k*np.pi/2)/k]
        vals=np.polynomial.chebyshev.chebval(xs,coeff)
        assert vals.min()>=0 and vals.max()<=1
        assert vals[xs<=-gap].max()<=1/8 and vals[xs>=gap].min()>=7/8
    print("Bounded polynomial margin filter verified at four gap sizes")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    study(13,4,(2,2,2),1,direct=True)
    study(23,5,(2,2,2),1,direct=True)
    study(101,17,(4,4,4),1)
    study(103,33,(4,8,4,8),2)
    cofactor_check(19)
    cofactor_check(43)
    margin_check()
