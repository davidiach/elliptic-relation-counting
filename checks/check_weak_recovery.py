"""Finite-law check of prefix correlations and path recovery on the original hash slots."""
from itertools import product
from math import prod

import numpy as np

from probe_target_jets import curve_cycle, gram_data, hash_families, moment_polynomial


def finite(p,a,width,m):
    curve=curve_cycle(p,a)
    N,M=len(curve),1<<p.bit_length()
    masks,mass=hash_families(curve,M,width)
    P=np.array(mass,dtype=float)/sum(mass)
    rho=width/M
    xs=[curve[j][0] for j in range(1,(N+1)//2)]
    inc=np.array([[j in mask for mask in masks] for j in range(1,len(xs)+1)])
    K=inc*P/rho
    R=K.mean(axis=0)
    center=(K-R)/2
    sigma=np.linalg.svd(center/np.sqrt(P),compute_uv=False)[0]
    assert 4*rho*sigma**2<=1+1e-10
    # All ordered signed tuples at T; their canonical-abscissa IDs are retained.
    def ab(i):
        return min(i,N-i)-1
    marks=np.array([tuple(map(ab,(*t,z))) for t in product(range(1,N),repeat=m-1)
                    if (z:=(1-sum(t))%N)])
    # Direct integer oracle H in the synthetic noisy-score control only.
    by_label=np.ones((len(masks),N),dtype=int)
    for j,slot in enumerate(masks):
        by_label[j]=[int(i in slot) for i in range(N)]
    def mix(A,B,i):
        return np.moveaxis(np.tensordot(A,B,axes=(1,i)),0,i)
    grid=np.indices((len(masks),)*m)
    point_fft=np.fft.fft(by_label,axis=1)
    spectra=np.prod(np.stack([point_fft[i] for i in grid]),axis=0)
    H=np.rint(np.fft.ifft(spectra,axis=-1).real[...,1]).astype(int)
    prob=prod(P[i] for i in grid)
    matched=prod(R[i] for i in grid)
    mu=np.sum(prob*H)
    planted=prob*H/mu
    assert np.isclose(mu,((N-1)**m-(-1)**m)/N*rho**m)
    noisy=.42+.16*(H>0)+.12*((grid[0]%2)==1)
    tests={"noisy witness control":noisy}
    if m==3:
        gr=np.array([gram_data([curve[i] for i in body],p) for body in masks])
        coeff=moment_polynomial([gr[i.ravel()] for i in grid],p,a)
        # A fixed arithmetic score; no training and no count labels as input.
        z2=coeff[:,0]
        tests["square character of Z2"]=np.array(
            [pow(int(x),(p-1)//2,p)==1 for x in z2],dtype=float).reshape(H.shape)
    for name,f in tests.items():
        delta=float(np.sum((planted-matched)*f))
        arrays=[f]
        for i in range(m):
            arrays.append(mix(center,arrays[-1],i))
        top=arrays[-1][tuple(marks.T)].mean()
        assert abs(delta-(2**m)*top)<=2*((1+1/(N-1))**m-1)+1e-10
        # Exact expected number of squared prefix coefficients in the root,
        # summed over prefix abscissas, with the suffix averaged over the mixture.
        for k in range(1,m+1):
            energy=rho**k*(arrays[k]**2).sum(axis=tuple(range(k)))
            assert np.max(energy)<=4**(-k)+1e-10
        # Compute actual path-survival probability exactly, backwards over suffixes.
        # Start at the full prefix. On revealing the next mixed seed, impose its barrier.
        threshold=abs(top)/2
        surv=0.
        for mark in marks:
            mixed=(K[mark]+R)/2
            path=1.
            for k in range(m,0,-1):
                alive=np.abs(arrays[k][tuple(mark[:k])])>=threshold
                path=alive*path
                if k>1:
                    path=np.expand_dims(path,0)*mixed[k-1].reshape((-1,)+(1,)*(m-k))
            surv+=np.sum(path)/len(marks)
        assert surv+1e-10>=abs(top)-threshold
        print(f"p={p}, m={m}, {name}: Q-R={delta:.6f}, mixed gap={2**m*top:.6f}, "
              f"planted prefix survival={surv:.5f}, guaranteed >= {abs(top)-threshold:.6f}")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    finite(13,4,2,3)
    finite(13,4,1,4)
