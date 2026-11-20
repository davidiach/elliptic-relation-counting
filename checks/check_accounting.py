"""Independent scalar, short-product, offset and quadrature checks."""
from itertools import product
from math import log, log2, pi, prod, sqrt
from random import Random

from flint import nmod_mat

from fixtures import curve_cycle, determinant, evaluation


def blocked_preimage(a, t, c, bits, limb, block):
    radix = 1 << limb
    n = (bits + limb-1)//limb
    coefficients = [0]*n
    for start in range(0, n, block):
        for j in range(start, min(start+block, n)):
            aj = (a >> (j*limb)) % radix
            for k in range(min(n-j, (t.bit_length()+limb-1)//limb)):
                coefficients[j+k] += aj*((t >> (k*limb)) % radix)
    carry = c
    x = 0
    for j, v in enumerate(coefficients):
        carry, digit = divmod(carry+v, radix)
        x += digit << (limb*j)
    return x % (1 << bits)


def check_hash():
    for alpha, beta, x in product(range(1, 16, 2), range(16), range(16)):
        inv = pow(alpha, -1, 16)
        residue = (alpha*x+beta) % 16
        for ell in range(4):
            accepted = residue < (1 << ell)
            preimages = [t for t in range(1 << ell)
                         if blocked_preimage(inv, t, -inv*beta % 16, 4, 1, 2) == x]
            assert preimages == ([residue] if accepted else [])
    rng = Random(6123)
    for bits in (13, 35, 61, 128):
        for _ in range(40):
            a = rng.getrandbits(bits) | 1
            c = rng.getrandbits(bits)
            t = rng.getrandbits(bits//3)
            assert blocked_preimage(a,t,c,bits,5,3) == (a*t+c) % (1 << bits)
    print("Short preimages: 8,192 exhaustive tests, 160 larger products.")


def check_offset():
    order = len(curve_cycle(101,17))
    slots = [(1,order-1), (1,order-1)]
    upper = min(map(len, slots))
    assert (2*len(slots)+7)*upper < order
    checked = 0
    for target in range(order):
        h, z = 0, 0
        for row in product(*slots):
            if sum(row) % order != target:
                continue
            h += 1
            for offset in range(1, order):
                partial = offset
                for addend in (*row, *((1,) if target == 0 else ()), -offset):
                    if (partial-addend) % order == 0 or (partial+addend) % order == 0:
                        break
                    partial = (partial+addend) % order
                else:
                    assert partial == (target or 1)
                    z += 1
        loss = 2*len(slots)+7
        assert 0 <= order*h-z <= loss*h
        assert (z+order-1)//order == h
        checked += 1
    assert checked > 0
    print(f"Offset multiplicity bounds: all {order} targets, {checked} strict rounding checks.")


def check_quadrature():
    p = 17
    # The other fixtures use constant term one. This check uses constant two.
    points = [(x,y) for x in range(p) for y in range(p)
              if (y*y-x**3-2*x-2) % p == 0]
    target = (0,6)
    dim, k = 4, 2
    n = dim*k
    nodes, rows = [], []
    for point in points:
        row = evaluation(point,n,p)
        if nmod_mat(rows+[row],p).rank() > len(rows):
            rows.append(row)
            nodes.append(point)
        if len(nodes) == n:
            break
    transform = nmod_mat(rows,p).inv()
    fixed = evaluation((target[0],-target[1] % p),dim,p)
    def powers(grid):
        out = [0]*k
        for weighted in product(*grid):
            choices, weights = zip(*weighted)
            delta = determinant([*(evaluation(P,dim,p) for P in choices),fixed],p)
            for j in range(k):
                out[j] = (out[j]+prod(weights)*pow(delta,j+1,p)) % p
        return out
    rng = Random(11239)
    for _ in range(24):
        original, compressed = [], []
        for _ in range(3):
            alpha, beta = rng.randrange(1,32,2), rng.randrange(32)
            slot = [P for P in points if (alpha*P[0]+beta) % 32 < 8]
            column_sums = [sum(evaluation(P,n,p)[i] for P in slot) % p for i in range(n)]
            weights = nmod_mat([column_sums],p)*transform
            original.append([(P,1) for P in slot])
            compressed.append([(P,int(weights[0,i])) for i,P in enumerate(nodes)])
        assert powers(original) == powers(compressed)
    print("Quadrature: both raw moments agreed on all 24 independent three-slot masks.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    K = 2*pi*sqrt(14*log2(90))/(9*log(2))
    old, new = 2*sqrt(19*K), 2*sqrt(18*K)
    assert abs(new-26.2933775155911) < 1e-12
    print(f"Conditional leading coefficients: {old:.10f} -> {new:.10f}, space {new/2:.10f}.")
    check_hash()
    check_offset()
    check_quadrature()
