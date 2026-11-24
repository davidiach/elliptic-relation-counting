"""Exact finite checks of pairwise mode-index tensorization.
This is a combinatorial identity check, not a benchmark of the contractor.
"""
from random import Random
from itertools import product
from functools import lru_cache
from math import ceil, log2, sqrt


def states(weights,h):
    if not weights:
        if h==0:yield ()
        return
    w=weights[0]
    for z in range(h//w+1):
        for tail in states(weights[1:],h-w*z):yield (z,)+tail


@lru_cache(None)
def count(weights,h):
    if h<0:return 0
    if not weights:return int(h==0)
    return sum(count(weights[1:],h-k*weights[0]) for k in range(h//weights[0]+1))


def unrank(weights,h,j):
    assert 0<=j<count(weights,h)
    if not weights:return ()
    for k in range(h//weights[0]+1):
        n=count(weights[1:],h-k*weights[0])
        if j<n:return (k,)+unrank(weights[1:],h-k*weights[0],j)
        j-=n
    raise AssertionError


def run():
    rng=Random(20261124);decoded=0;checks=0;networks=0;bad_codes=0
    for _ in range(80):
        r=rng.randrange(2,5);p=rng.randrange(3,8)
        possible=[(i,j) for i in range(r) for j in range(r) if i!=j]
        rng.shuffle(possible)
        edges=possible[:min(len(possible),rng.randrange(2,5))]
        W=[tuple(rng.randrange(1,4) for _ in range(rng.randrange(1,4))) for _ in edges]
        inc=[[i for i,e in enumerate(edges) if g in e] for g in range(r)]
        # A reproducible arbitrary local tensor on the incident exponent tuples.
        def f(g,keys):
            z=g+1
            for v in keys:
                for k in v:z=(z*19+k+1)%251
            return z%13-6
        def term(beta):
            z=1
            for g in range(r):z*=f(g,[beta[i] for i in inc[g]])
            return z
        all_states=[[(h,x) for h in range(p) for x in states(w,h)] for w in W]
        direct=0;all_count=0
        for terms in product(*all_states):
            hs,vs=zip(*terms)
            if sum(hs)>=p:continue
            direct+=2**sum(hs)*term(vs);all_count+=1
        reencoded=0
        for hs in product(range(p),repeat=len(W)):
            if sum(hs)>=p:continue
            N=[count(w,h) for w,h in zip(W,hs)]
            if not all(N):continue
            assert __import__('math').prod(N)<=all_count
            bits=[(n-1).bit_length() for n in N]
            assert sum(bits)<=log2(all_count)+len(W)+1e-8
            for j,(w,h,n) in enumerate(zip(W,hs,N)):
                v=list(states(w,h));assert len(v)==n
                assert [unrank(w,h,k) for k in range(n)]==v
                decoded+=n
            # Sum over all bit codes, rejecting invalid labels at their endpoints.
            for js in product(*[range(1<<b) for b in bits]):
                if any(j>=n for j,n in zip(js,N)):
                    bad_codes+=1;continue
                beta=[unrank(w,h,j) for w,h,j in zip(W,hs,js)]
                reencoded+=2**sum(hs)*term(beta);checks+=1
        assert reencoded==direct
        assert reencoded%(2**p)==direct%(2**p)
        networks+=1
    print(f'{networks} mode networks; {decoded:,} exact rank/unrank checks; {checks:,} valid code assignments; {bad_codes:,} invalid padded codes rejected.')

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    run()
