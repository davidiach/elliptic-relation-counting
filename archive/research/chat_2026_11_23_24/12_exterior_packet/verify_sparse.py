from itertools import product
from math import comb,log2,e
from random import Random
rng=Random(4921)

def checks():
    cases=evals=groups=0
    # Each row's capacity restricts the total Hamming support.  Public baselines
    # may differ per case.  Local functions depend on one ternary seed or no seed.
    for trial in range(120):
        m=3;M=rng.randint(4,9);n=4
        ends=[(rng.randrange(2),2+rng.randrange(2)) for _ in range(M)]
        ports=[[i for i,p in enumerate(ends) if v in p] for v in range(n)]
        who=[rng.randrange(-1,m) for _ in range(n)]
        cap=[rng.randrange(len(ports[v])+1) for v in range(2)];H=sum(cap)
        base=[rng.randrange(2) for _ in range(M)]
        tables=[]
        for v in range(n):
            tb={}
            for s in range(3):
                for bits in product(range(2),repeat=len(ports[v])):
                    z=sum(x^base[i] for x,i in zip(bits,ports[v]))
                    tb[s,bits]=0 if v<2 and z>cap[v] else rng.randrange(2)
            tables.append(tb)
        E=[set(i for v in range(n) if who[v]==j for i in ports[v]) for j in range(m)]
        ks=[len(q) for q in E]
        # Threshold reductions are allowed only while the pigeonhole property holds.
        # If all local edge sets are small, keep a fully assigned mode.
        if sum(k+1 for k in ks)>2*H:
            order=list(range(m));rng.shuffle(order)
            for j in order:
                while ks[j]>0 and sum(k+1 for k in ks)-1>2*H:ks[j]-=1
        else:
            # A single unrestricted seed gives a valid slice decomposition.
            ks=[len(q) for q in E]
        blocks={}
        for bits in product(range(2),repeat=M):
            dev={i for i in range(M) if bits[i]!=base[i]}
            if len(dev)>H:continue
            j=next(i for i in range(m) if len(dev&E[i])<=ks[i])
            key=(j,tuple(sorted(dev&E[j])))
            blocks.setdefault(key,[]).append(bits)
        assert len(blocks)<=sum(sum(comb(len(E[i]),j) for j in range(ks[i]+1)) for i in range(m))
        # In each slice block, the seed-i factor is identical for all assignments.
        def val(v,bits,s):return tables[v][s if who[v]>=0 else 0,tuple(bits[i] for i in ports[v])]
        for seed in product(range(3),repeat=m):
            direct=0
            for bits in product(range(2),repeat=M):
                term=1
                for v in range(n):term*=val(v,bits,seed[who[v]])
                direct^=term
            got=0
            for (j,idx),assigns in blocks.items():
                u=1
                for v in range(n):
                    if who[v]==j:u*=val(v,assigns[0],seed[j])
                z=0
                for bits in assigns:
                    uu=1;vv=1
                    for v in range(n):
                        if who[v]==j:uu*=val(v,bits,seed[j])
                        else:vv*=val(v,bits,seed[who[v]])
                    assert uu==u;z^=vv
                got^=u*z
            assert got==direct;evals+=1
        cases+=1;groups+=len(blocks)
    print('Sparse seed-local networks:',cases)
    print('Exact contraction / slice decomposition comparisons:',evals)
    print('Slice blocks checked:',groups)
    nchecks=0
    for _ in range(2000):
        m=rng.randrange(2,12);cs=[rng.randrange(1,150) for _ in range(m)];B=rng.randrange(1,min(cs)+1)
        if any(c<=B for c in cs):continue
        ks=[]
        for c in cs:
            total=1;k=0
            while k<c and total+comb(c,k+1)<=2**B:k+=1;total+=comb(c,k)
            assert total<=2**B and total+comb(c,k+1)>2**B
            ks.append(k)
        T=sum(k+1 for k in ks);M=(sum(cs)+1)//2;H=(T+1)//2
        assert H<=M
        assert m*B<2*H*log2(e*M/H)
        assert M*H>B*B*m*m/16-1e-9
        nchecks+=1
    print('Product and strengthened entropy inequality checks:',nchecks)
    print('ALL ASSERTIONS PASSED')
if __name__=='__main__':checks()
