"""Check the mixed mask law, the geometric counting bound and specialized C38 ports."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import log, log2, pi, prod, sqrt
from random import Random

from strided_compiler import arith, free_positions, sc, specialize


def curve(p,a,c):
    pts=[None]+[(x,y) for x,y in product(range(p),repeat=2) if (y*y-x**3-a*x-c)%p==0]
    ids={P:i for i,P in enumerate(pts)}
    add=[[ids[arith.point_add(P,Q,p,a)] for Q in pts] for P in pts]
    return pts,add


def geometry(p,a,c,positions):
    pts,add=curve(p,a,c)
    M=1<<p.bit_length()
    free=sum(1<<j for j in positions)
    subs=[i for i in range(M) if not (i & ~free)]
    diff={u-v for u in subs for v in subs}
    runs=sum(j-1 not in positions for j in positions)
    assert len(diff)<=len(subs)*(1<<runs)
    peak=0
    for D in range(1,len(pts)):
        fibers=Counter()
        same=0
        for i,P in enumerate(pts[1:],1):
            Q=pts[add[i][D]]
            if Q is None:
                continue
            fibers[(Q[0]-P[0])%p]+=1
            if (Q[0]&~free)==(P[0]&~free):
                assert Q[0]-P[0] in diff
                same+=1
        peak=max(peak,max(fibers.values(),default=0))
        assert same<=4*len(diff)
    assert peak<=4
    print(f"Geometry p={p}, N={len(pts)}, runs={runs}, difference-set size {len(diff)}, "
          f"maximum fiber {peak}")
    return pts,add


def exact_law(p,a,c,m,positions):
    pts,add=geometry(p,a,c,positions)
    N=len(pts)
    M=1<<p.bit_length()
    width=1<<len(positions)
    free=sum(1<<j for j in positions)
    universal=Counter(tuple(i for i,P in enumerate(pts[1:],1) if (alpha*P[0]+beta)%M<width)
                      for alpha in range(1,M,2) for beta in range(M))
    boxes=Counter(tuple(i for i,P in enumerate(pts[1:],1) if P[0]&~free==code)
                  for code in range(M) if code&free==0)
    laws=[universal]+[boxes]*(m-1)
    total=prod(sum(C.values()) for C in laws)
    first,second,guarded_first,covered=[[0]*N for i in range(4)]
    for masks in product(*laws):
        mass=prod(C[A] for C,A in zip(laws,masks))
        h,guarded=Counter(),Counter()
        for choices in product(*masks):
            S=choices[0]
            valid=True
            for v in choices[1:]:
                valid=valid and bool(S) and pts[S][0]!=pts[v][0]
                S=add[S][v]
            h[S]+=1
            if valid and S:
                guarded[S]+=1
        for t,count in h.items():
            first[t]+=mass*count
            second[t]+=mass*count*count
        for t,count in guarded.items():
            guarded_first[t]+=mass*count
            covered[t]+=mass
    rho=Fraction(width,M)
    lam=N**(m-1)*rho**m
    mu=Fraction((N-1)**m-(-1)**m,N)*rho**m
    runs=sum(j-1 not in positions for j in positions)
    limit=rho**2*(2*width)**(2*m-2)*(1+Fraction(2*(m-1)*(1<<runs),width))+2*lam
    assert all(Fraction(s,total)==mu for s in first[1:])
    assert all(Fraction(s,total)<=limit for s in second[1:])
    if N>3*m:
        expected_guarded=(1-Fraction(3*m,N))*lam
        assert all(Fraction(s,total)>=expected_guarded for s in guarded_first[1:])
        bound=expected_guarded**2/limit
        assert all(Fraction(s,total)>=bound for s in covered[1:])
    print(f"Exact m={m}, {total:,} seed combinations, all {N-1} affine targets, "
          f"guarded coverage {min(covered[1:])/total:.4f}..{max(covered[1:])/total:.4f}")


def accepts(net,labels):
    ring=sc.Ring(1)
    return set(labels)==set(net.edges) and all(
        t.entry(tuple(labels[e] for e in t.ports),ring)==ring.one for t in net.tensors)


def layout_checks():
    cases=((127,3,4,37),(1024,7,8,155),(16003,25,32,2001),(65537,30,32,5891))
    for b,k,V,ell in cases:
        bits=set(free_positions(b,k,V,ell))
        active={(i//k)%V for i in bits}
        runs=sum(i-1 not in bits for i in bits)
        minimum=k*(b//(k*V))
        assert len(bits)==ell and len(active)<=(ell+minimum-1)//minimum
        assert 1<runs<=(b+k*V-1)//(k*V)
        assert any(i not in bits and (i//k)%V in active for i in range(b))
        for start in range(0,b,k*V):
            block=sorted(i for i in range(start,min(b,start+k*V)) if i in bits)
            assert block==list(range(start,start+len(block)))
    print(f"Strided layout verified on {len(cases)} nontrivial shapes up to {cases[-1][0]:,} bits")


def graph_checks():
    rng=Random(202611201)
    for b,L,offset,m,dimension in ((4,8,3,3,1),(7,16,27,5,2)):
        data,values,_=arith.make_instance(b,L,offset,m,rng,"general")
        par=data["par"]
        positions=free_positions(b,par["k"],par["L2"],data["hash_ell"])
        mask=sum(1<<j for j in positions)
        patterns={f"x{i}":values[f"x{i}"]&~mask for i in range(2,m+1)}
        base,net,pins=specialize(data,dimension,positions,patterns)
        full=base.encode_witness(arith.witness(base.data,values))
        assert base.accepts(full) and all(full[e]==v for e,v in pins.items())
        labels={e:full[e] for e in net.edges}
        assert accepts(net,labels)
        active={(i//par["k"])%par["L2"] for i in positions}
        saved=sum(base.paths[e][0]=="hs" for e in pins)
        assert saved==(m-1)*(par["L2"]-len(active))*par["L1"]*base.codec.width
        assert sum(base.paths[e][0]=="outs" for e in pins)==(m-1)*par["L"]*base.codec.width
        # Vary a retained private word and a fixed bit in an active owner.
        e=next(e for e in net.edges if base.paths[e][:2]==("hs","x2"))
        corrupted=dict(labels)
        corrupted[e]=(labels[e]+1)%base.codec.Q
        assert not accepts(net,corrupted)
        fixed=[j for j in range(b) if not mask>>j&1]
        bad_bit=next((j for j in fixed if (j//par["k"])%par["L2"] in active),fixed[0])
        changed={**patterns,"x2":patterns["x2"]^(1<<bad_bit)}
        _,bad_net,_=specialize(data,dimension,positions,changed)
        assert not accepts(bad_net,labels)
        print(f"C38 b={b}, m={m}, alphabet={base.codec.Q}: "
              f"{len(base.net.edges):,} -> {len(net.edges):,} edges, "
              f"{saved} fixed input symbols, corrupted cases rejected")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    exact_law(13,4,1,3,(0,2))
    exact_law(13,4,1,4,(1,))
    exact_law(19,-1,0,3,(0,2))
    geometry(101,17,1,(1,2,5))
    layout_checks()
    graph_checks()
    K=2*pi*sqrt(14*log2(90))/(9*log(2))
    print(f"conditional time constant {2*sqrt(15*K):.12f}, space constant {sqrt(15*K):.12f}")
