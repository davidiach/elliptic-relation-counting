"""Independent checks and fresh draws for the weighted dense-pool decision."""
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, product
from math import prod
from random import Random

from flint import nmod_poly

from dense_weighted import complete_locator, decide, decide_slots, duplicate_probes, join, moment_channels, translate, visible_locator

from fixtures import curve_cycle


def unit_checks():
    p,a,c = 37,27,2
    X = nmod_poly([0,1],p)
    slots = [[2,32],[1,21],[4,5,6,16,17,18,22,23,27]]
    rng = Random(449)
    for size in range(1,18):
        xs = [0]+rng.sample(range(1,p),size-1)
        poly = prod((X-u for u in xs),start=nmod_poly([1],p))
        for t in (0,1,2,13):
            expected = prod(((t-u)**2*X**2-2*((t+u)*(u*t+a)+2*c)*X
                             +(u*t-a)**2-4*c*(u+t) for u in xs),start=nmod_poly([1],p))
            assert translate(poly,t,a,c,p) == expected
        for _ in range(4):
            tags = [rng.randrange(1,p) for x in xs]
            moments = [sum(h*pow(x,j,p) for h,x in zip(tags,xs)) % p for j in range(2*size)]
            assert complete_locator(visible_locator(moments,p),size,p)==poly
            moments = [sum(h*pow(x,j,p) for h,x in zip(tags[1:],xs[1:])) % p for j in range(2*size)]
            assert complete_locator(visible_locator(moments,p),size,p) is None
    for third,truth in ((slots[2],True),([3,6,8,9,14,17,22,23,36],False)):
        slots[2] = third
        tags = [[rng.randrange(1,p) for x in s] for s in slots]
        channels = moment_channels(slots,tags,p)
        for i,seq in enumerate(channels):
            brute = [0]*len(seq)
            for indices in product(*(range(len(s)) for s in slots)):
                u,v,w = xs = [s[k] for s,k in zip(slots,indices)]
                tagged = prod(tags[j][indices[j]] for j in range(3) if j!=i)
                amp = tagged*((u-v)*(v-w)*(w-u))**2
                for l in range(len(seq)):
                    brute[l] += amp*pow(xs[i],l,p)
            assert seq == [s % p for s in brute]
        for attempt in range(20):
            found = decide(list(map(len,slots)),channels,13,a,c,p)
            if found is not None:
                break
            tags = [[rng.randrange(1,p) for x in s] for s in slots]
            channels = moment_channels(slots,tags,p)
        assert found == truth


def fresh(p,a,widths,draws,seed):
    curve = curve_cycle(p,a)
    order,mask = len(curve),(1<<p.bit_length())-1
    by_x = defaultdict(list)
    for j,(x,y) in enumerate(curve[1:],1):
        by_x[x].append(j)
    seed_rng,tag_rng = Random(seed),Random(seed+1)
    counts = Counter()
    for _ in range(draws):
        slots = []
        seeds = []
        for width in widths:
            alpha=seed_rng.randrange(1,mask+1,2)
            inv = pow(alpha,-1,mask+1)
            beta = seed_rng.randrange(mask+1)
            seeds.append((alpha,beta,width))
            xs = [((u-beta)*inv)&mask for u in range(width)]
            slots.append([x for x in xs if x in by_x])
        signed = [[j for x in s for j in by_x[x]] for s in slots]
        last = set(signed[2])
        hit_count = sum((1-i-j)%order in last for i in signed[0] for j in signed[1])
        counts["all hits"] += hit_count>0
        flat = [x for s in slots for x in s]
        regular = len(set(flat)) == len(flat) and curve[1][0] not in flat
        if not regular or not all(slots):
            counts["outside hits"] += hit_count>0
        dense = min(map(len,slots))>=2 and regular
        counts["eligible"] += regular and all(slots)
        counts["dense"] += dense
        counts["dense hits"] += dense and hit_count>0
        counts["dense single hits"] += dense and hit_count==1
        counts["dense negatives"] += dense and hit_count==0
        sizes = list(map(len,slots))
        first = None
        seen=[nmod_poly([1],p) for s in slots]
        for repeat in range(8):
            tags = [[tag_rng.randrange(1,p) if repeat else 1 for x in s] for s in slots]
            ans = decide(sizes,moment_channels(slots,tags,p),curve[1][0],a,1,p,seen_locators=seen)
            if repeat==0:
                first = ans
            if ans is not None:
                break
        counts["retry inputs"] += first is None
        counts["final abstentions"] += ans is None
        if ans is not None:
            assert ans == (hit_count>0)
            counts["dense single certified"] += dense and hit_count==1
        else:
            counts["abstention positives"] += hit_count>0
            counts["regular abstentions"] += regular
            sums=[(sum(s)%p,sum(x*x for x in s)%p) for s in slots]
            repaired=decide(sizes,moment_channels(slots,tags,p),curve[1][0],a,1,p,sums,seen)
            if repaired is not None:
                assert repaired == (hit_count>0)
                counts["two-sum repairs"] += 1
                continue
            counts["after repair abstentions"] += 1
            counts["after repair positive"] += hit_count>0
            pair=next((i,j) for i,j in combinations(range(3),2)
                      if sizes[i]==sizes[j]==1 and sums[i][0]==sums[j][0])
            i,j=pair
            probes=duplicate_probes(sums[i][0],curve[1][0],a,1,p)
            alpha,beta,width=seeds[3-i-j]
            assert any(((alpha*z+beta)&mask)<width for z in probes)==(hit_count>0)
    print(p,widths,draws,dict(counts))


def check_overlaps():
    p,a=13,4
    curve=curve_cycle(p,a)
    N=len(curve)
    X=nmod_poly([0,1],p)
    labels=[s for k in (1,2) for s in combinations(range(1,5),k)]
    polys=[prod((X-curve[i][0] for i in s),start=nmod_poly([1],p)) for s in labels]
    signed=[(*s,*(N-i for i in s)) for s in labels]
    rng=Random(1007)
    checked=0
    for ids in combinations_with_replacement(range(len(labels)),3):
        truths={sum(tup)%N for tup in product(*(signed[i] for i in ids))}
        slots=[tuple(curve[q][0] for q in labels[h]) for h in ids]
        for t in range(1,(N+1)//2):
            assert join([polys[i] for i in ids],curve[t][0],a,1,p)==(t in truths)
            assert decide_slots(slots,curve[t][0],a,1,p,rng,lambda i,x:x in slots[i])==(t in truths)
            for i,j in combinations(range(3),2):
                if ids[i]==ids[j] and len(labels[ids[i]])==1:
                    u=curve[labels[ids[i]][0]][0]
                    probes=duplicate_probes(u,curve[t][0],a,1,p)
                    third={curve[q][0] for q in labels[ids[3-i-j]]}
                    assert bool(probes & third)==(t in truths)
            checked+=1
    print(f"{checked} overlapping or disjoint support/target cases checked against literal sums")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    unit_checks()
    check_overlaps()
    for params in ((103,33,(8,8,8),3000,665580),(1031,9,(8,16,16),3000,665590),
                   (65521,44,(16,32,32),3000,665600)):
        fresh(*params)
