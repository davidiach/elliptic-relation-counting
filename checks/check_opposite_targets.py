"""Compare mixed target evaluations with literal determinants at every toy target."""
from itertools import combinations, product
from math import prod
from random import Random

from flint import nmod_mat, nmod_poly

from probe_target_jets import curve_cycle
from fixtures import evaluation


def cross_polynomial(slots, weights, curve, p, a):
    m = len(slots)
    d = m+1
    pieces = []
    for slot,w in zip(slots,weights):
        rows = [evaluation(curve[i],d,p) for i in slot]
        pieces.append([[sum(weight*r[j]*r[k] for weight,r in zip(w,rows)) % p
                        for k in range(d)] for j in range(d)])
    X = nmod_poly([0,1],p)
    f = X**3+a*X+1
    even = [X**j for j in range(d//2+1)]
    odd = [X**j for j in range((d-1)//2)]
    result = nmod_poly([],p)
    for n in range(1,m+1):
        for subset in combinations(pieces,n):
            gram = [[sum(piece[i][j] for piece in subset) % p for j in range(d)] for i in range(d)]
            for offset,basis,scale in [(0,even,1),(len(even),odd,-f)]:
                for i,h in enumerate(basis,offset):
                    for j,l in enumerate(basis,offset):
                        minor = [r[:j]+r[j+1:] for k,r in enumerate(gram) if k!=i]
                        cofactor = (-1)**(i+j)*int(nmod_mat(minor,p).det())
                        result += (-1)**(m-n)*scale*cofactor*h*l
    assert result.degree() <= d
    return result


def check(shape, p=23, a=5):
    curve = curve_cycle(p,a)
    order = len(curve)
    slots = []
    start = 1
    for width in shape:
        slots.append([i for j in range(start,start+width) for i in (j,order-j)])
        start += width
    assert start-1 <= (order-1)//2
    m = len(slots)
    rng = Random(4983+m)
    channels = [[[rng.randrange(p) for i in slot] for slot in slots] for _ in range(2)]
    polys = [cross_polynomial(slots,weights,curve,p,a) for weights in channels]
    forced = {curve[slot[0]][0] for slot in slots if len(slot)==2}
    found = set()
    all_vanish = set()
    for t in range(1,(order+1)//2):
        sums = [0,0]
        individual = set()
        for indices in product(*(range(len(slot)) for slot in slots)):
            points = [slot[i] for slot,i in zip(slots,indices)]
            rows = [evaluation(curve[i],m+1,p) for i in points]
            signed = [int(nmod_mat(rows+[evaluation(curve[j],m+1,p)],p).det()) for j in (t,order-t)]
            val = prod(signed) % p
            individual.add(val)
            for j,w in enumerate(channels):
                sums[j] += val*prod(v[i] for v,i in zip(w,indices))
        x = curve[t][0]
        assert [int(poly(x)) for poly in polys] == [s % p for s in sums]
        if individual == {0}:
            all_vanish.add(x)
        if all(s % p == 0 for s in sums):
            found.add(x)
    assert all_vanish == forced and forced <= found
    print(f"m={m}, shape {shape}, forced {sorted(forced)}, extra two-channel roots {sorted(found-forced)}")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    for shape in ((2,2),(1,2,3),(2,2,2),(1,2,1,2),(2,2,2,2),
                  (1,2,1,2,1),(1,1,1,1,2,1)):
        check(shape)
