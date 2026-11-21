"""Compare direct scalar resultants with independent integer/field references."""

from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import prod
from pathlib import Path
from random import Random

import numpy as np

from scalar import curve_jet, fixed_resultant, form, full_jet, polynomial, split_norm, two_slot_jet


def previous_verifier():
    path = Path(__file__).resolve().parents[1] / "joint_parity" / "verify.py"
    spec = spec_from_file_location("joint_norm_fixture", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


old = previous_verifier()


def dense_jet(sets, target, p, a, c):
    """Determinant and cofactors of the old product-algebra multiplication matrix."""
    if not all(sets):
        return 1, [0]*(len(sets)+2)
    algebra = old.SlotAlgebra(sets, p, a, c)
    coeff = algebra.coeff
    matrix = sum(pow(target, i, p)*f for i, f in enumerate(coeff)) % p
    dz = sum((i*pow(target, i-1, p) % p)*coeff[i] for i in range(1, len(coeff))) % p
    cofactors = np.zeros_like(matrix)
    for i, j in np.ndindex(matrix.shape):
        minor = np.delete(np.delete(matrix, i, axis=0), j, axis=1)
        cofactors[i, j] = (-1)**(i+j)*old.determinant(minor, p) % p
    return old.determinant(matrix, p), [int(np.sum(cofactors*f)) % p
                                        for f in (dz, algebra.identity, *algebra.coordinates)]


def raw_product(sets, target, p, a, c):
    """Literal scalar-product reference. Only the cross-check visits the full grid."""
    total, gradient = 1, [0]*(len(sets)+2)
    for row in product(*sets):
        factor = form(row, p, a, c)
        value = int(factor(target))
        directions = [int(factor.derivative()(target)), 1, *row]
        gradient = [(x*value+total*y) % p for x, y in zip(gradient, directions)]
        total = total*value % p
    return total, gradient


def small():
    p, a, c = 13, 4, 1
    curve = old.points(p, a, c)
    targets = [x for x, y in curve if 0 < y < p-y]
    counts = Counter()
    root = ((3, 4, 2), (7, 9, 2), (1, 8, 4))
    for box in product(*(old.children(s) for s in root)):
        sets = [sorted({x for x, y in slot}) for slot in old.hash_slots(curve, 16, box)]
        for target in targets:
            ref = dense_jet(sets, target, p, a, c)
            assert full_jet(sets, target, p, a, c) == ref
            assert two_slot_jet(sets, target, p, a, c) == ref
            v, (dz, c0, *_) = ref
            counts["nonzero" if v else "simple" if dz else "one-repeated" if c0 else "several"] += 1
    print(f"All 63 subwindows, {sum(counts.values())} exact scalar and gradient comparisons: "
          f"{dict(counts)}.")
    for sets in ([[0, 2]], [[0, 2]]*4, [[0, 2]]*5):
        for target in targets:
            ref = raw_product(sets, target, p, a, c)
            assert full_jet(sets, target, p, a, c) == ref == two_slot_jet(sets, target, p, a, c)
    # One abscissa tuple has two completions at the target. Its marker derivative survives.
    sets = [[0]]*3
    v, (dz, ds, *_) = full_jet(sets, 0, p, a, c)
    assert (v, dz, ds) == (0, 0, 1)
    print("Single-slot, repeated factors, one multiple zero, and degree exceeding the characteristic checked.")


def homogeneous():
    p = 13
    rng = Random(23)
    checks = 0
    for _ in range(100):
        # Known finite and infinite roots, including coincident or repeated roots.
        finite = [rng.randrange(4) for _ in range(rng.randrange(5))]
        infinity = rng.randrange(3)
        lead = rng.randrange(1, p)
        A = polynomial([lead], p)
        for x in finite:
            A *= polynomial([-x, 1], p)
        h = rng.randrange(1, 6)
        B = polynomial([rng.randrange(p) for _ in range(rng.randrange(1, h+2))], p)
        V = polynomial([rng.randrange(p) for _ in range(h+1)], p)
        if not B:
            B = polynomial([1], p)
        val, diff = fixed_resultant(A, len(finite)+infinity, B, h, [V], p)
        ref, slope = (-1)**(infinity*h)*pow(lead, h, p) % p, 0
        for f, g in [(int(B[x]), int(V[x])) for x in [h]*infinity] + [
                (int(B(x)), int(V(x))) for x in finite]:
            ref, slope = ref*f % p, (slope*f+ref*g) % p
        assert (val, diff) == (ref, [slope])
        checks += 1
    print(f"{checks} homogeneous resultant and tangent checks, including infinite roots.")


def sibling_mean():
    valid = {x for x, y in old.points(13, 4, 1) if y}
    M = 16
    for width in (2, 4, 8):
        half = width//2
        bins = Counter(x % width for x in valid)
        predicted = half*width*sum(bins[j]*bins[(j+half) % width] for j in range(width))
        observed = 0
        for A in range(1, M, 2):
            for C in range(M):
                observed += sum((A*j+C) % M in valid and (A*(j+half)+C) % M in valid
                                for j in range(half))
        assert 2*observed == predicted
    print("Dyadic sibling expectation checked at three widths over all 128 seeds.")


def torsion():
    p, a, c = 7, -1, 0
    curve = old.points(p, a, c)
    point_slots = old.hash_slots(curve, 8, ((1, 0, 4), (3, 1, 4), (5, 3, 4)))
    lookup = {x: (x, y) for x, y in curve}
    sets = [sorted({x for x, _ in slot}) for slot in point_slots]
    checked = 0
    for target in (q for q in curve if q[1]):
        value, gradient = 1, [0]*(len(sets)+2)
        for row in product(*sets):
            shifted, ordinary = target, []
            for x in row:
                if lookup[x][1]:
                    ordinary.append(x)
                else:
                    shifted = old.point_add(shifted, lookup[x], p, a)
            if not ordinary:
                continue
            f = form(ordinary, p, a, c)
            v = int(f(shifted[0]))
            changes = [int(f.derivative()(shifted[0])), 1, *row]
            gradient = [(old_v*v+value*new) % p for old_v, new in zip(gradient, changes)]
            value = value*v % p
        assert curve_jet(point_slots, target, p, a, c, old.point_add) == (value, gradient)
        checked += 1
    print(f"Order 8 with full rational two-torsion: all {checked} eligible signed targets.")


def large():
    p, a, c, M = 65521, 44, 1, 65536
    rng = Random(18112026)
    squares = {y*y % p: y for y in range(1, p)}
    rhs = lambda x: (x**3+a*x+c) % p
    assert 1 + sum(2 if rhs(x) in squares else int(rhs(x) == 0) for x in range(p)) == 65419
    widths = (8, 16, 16, 16)
    for _ in range(2):
        sets = []
        for width in widths:
            alpha, beta = rng.randrange(1, M, 2), rng.randrange(M)
            xs = [pow(alpha, -1, M)*(i-beta) % M for i in range(width)]
            sets.append([x for x in xs if x < p and rhs(x) in squares])
    dims = tuple(map(len, sets))
    assert dims == (4, 7, 10, 10)
    val, dz, geometric = split_norm(sets, 0, p, a, c)
    assert not val and dz
    xs = [scale*pow(shift, -1, p) % p for shift, scale in geometric]
    assert all(x in slot for x, slot in zip(xs, sets))
    jet = full_jet(sets, 0, p, a, c)
    assert jet == raw_product(sets, 0, p, a, c)
    assert xs == [d*pow(jet[1][1], -1, p) % p for d in jet[1][2:]]
    print(f"Order 65419: full grid {prod(dims)}, half tuple counts {prod(dims[:2])} and "
          f"{prod(dims[2:])}, resultant degrees {2*prod(dims[:2])} and {4*prod(dims[2:])}.")
    print(f"At target (0,1): norm {val}, derivative in z {dz}, additive derivatives {jet[1][1:]}.")
    # Inspect actual index siblings at a tuned two-slot width on the same curve.
    width, alpha, beta = 256, 7, 37
    preimage = [pow(alpha, -1, M)*(j-beta) % M for j in range(width)]
    steps = Counter()
    pairs = 0
    for x, y in zip(preimage[:width//2], preimage[width//2:]):
        if max(x, y) >= p or rhs(x) not in squares or rhs(y) not in squares:
            continue
        pairs += 1
        ordinates = (squares[rhs(x)], squares[rhs(y)])
        candidates = {old.point_add((x, s*ordinates[0] % p), (y, t*ordinates[1] % p), p, a)
                      for s, t in product((-1, 1), repeat=2)}
        steps.update(candidates)
    peak = max(steps.values())
    assert pairs > 8 and peak <= 8 and None not in steps
    print(f"Width 256: {pairs} valid sibling pairs, at most {peak} covered by one curve translation.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    homogeneous()
    sibling_mean()
    small()
    torsion()
    large()
