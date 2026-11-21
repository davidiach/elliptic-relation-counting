"""Compare the joint section product with literal summation forms on finite fixtures."""

from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from random import Random

from flint import nmod_poly

from joint import (CurveTable, SeedProduct, batch_remainders, chord_scale, evaluations, form,
                   polynomial, product_tree, research)

reference_path = research / "scalar_norm" / "verify.py"
spec = spec_from_file_location("scalar_check", reference_path)
reference = module_from_spec(spec)
spec.loader.exec_module(reference)


def table_for(p, a, r):
    return CurveTable(p, a, 1, r, (0, 1), reference.old.point_add)


def check_scales(table, exhaustive):
    p, r = table.p, table.r
    rng = Random(711)
    pairs = product(range(r), repeat=2) if exhaustive else [(rng.randrange(r), rng.randrange(r)) for _ in range(2000)]
    for i, j in pairs:
        assert (chord_scale(table.points[i], table.points[j], p)
                * table.gamma[(i+j) % r]*table.gamma[(i-j) % r]) % p == (table.gamma[i]*table.gamma[j])**2 % p
    assert all(pow(g, (p-1)//2, p) == 1 for g in table.gamma)


def direct_poly(table, slots):
    f = polynomial([1], table.p)
    for row in product(*slots):
        f *= form(row, table.p, table.a, table.c)
    return f


def strip(f, t, p):
    linear = polynomial([-t, 1], p)
    h = 0
    while not f(t):
        f = f // linear
        h += 1
    return h, int(f(t))


def compare(table, slots, counter):
    p = table.p
    seed = SeedProduct(table, slots)
    targets = list(table.index)
    direct = direct_poly(table, slots)
    for t, (h, unit, dz) in zip(targets, seed.target_units(targets)):
        assert (h, unit) == strip(direct, t, p)
        assert dz == int(direct.derivative()(t))
        if h:
            expected = reference.raw_product(slots, t, p, table.a, table.c)[1]
            assert seed.zero_gradient(t, h, unit, dz) == expected
            counter["zero gradients"] += 1
            counter["single tuple, multiple signs"] += h > 1 and expected[1] != 0
            counter["multiple zero tuples"] += expected[1] == 0
        counter["target jets"] += 1
    return seed


def small():
    table = table_for(13, 4, 19)
    check_scales(table, True)
    counts = Counter()
    parent = ((3, 4, 2), (7, 9, 2), (1, 8, 4))
    for box in product(*(reference.old.children(seed) for seed in parent)):
        slots = table.windows(16, box)
        ordinary = reference.old.hash_slots(table.points[1:], 16, box)
        assert list(map(set, slots)) == [{x for x, y in slot} for slot in ordinary]
        compare(table, slots, counts)
    for slots in ([[0, 2]], [[0]]*6, [[0, 2]]*4):
        compare(table, slots, counts)
    exceptional = SeedProduct(table, [[0]]*6)
    assert any(h >= table.p for h, _, _ in exceptional.target_units(table.index))
    print(f"Order 19, all 63 dyadic boxes and three additional inputs: {dict(counts)}.")


def bit_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            col = row.bit_length()
            if col not in pivots:
                pivots[col] = row
                break
            row ^= pivots[col]
    return len(pivots)


def sign_rank(table, t, width, direct):
    p, r = table.p, table.r
    bits = []
    for point in table.points:
        u = -1 if point is None else 1 if point[0] == t else t-point[0]
        bits.append(int(pow(u, (p-1)//2, p) != 1))
    binary = nmod_poly(bits, 2)
    modulus = nmod_poly([1]+[0]*(r-1)+[1], 2)
    rank = (r-binary.gcd(modulus).degree()-sum(bits) % 2)//2
    n = (r-1)//2
    if direct:
        rows = [sum((bits[(j-i) % r]^bits[(j+i) % r]) << (j-1) for j in range(1, n+1))
                for i in range(1, n+1)]
        assert bit_rank(rows) == rank
        M = 1 << p.bit_length()
        masks = table.windows(M, [(7, beta, width) for beta in range(M)])
        mask_rows = [sum(1 << (table.index[x]-1) for x in slot) for slot in masks]
        assert bit_rank(mask_rows) == n
        rng = Random(29)
        for _ in range(24):
            seeds = [(rng.randrange(1, M, 2), rng.randrange(M), width) for _ in range(2)]
            left, right = table.windows(M, seeds)
            right_mask = sum(1 << (table.index[x]-1) for x in right)
            bit = sum((rows[table.index[x]-1] & right_mask).bit_count() for x in left) % 2
            _, unit = strip(direct_poly(table, [left, right]), t, p)
            assert bit == (pow(unit, (p-1)//2, p) != 1)
    print(f"Order {r}, target x={t}: cofactor sign kernel rank {rank} out of {n}.")
    return rank


def large():
    table = table_for(65521, 44, 65419)
    check_scales(table, False)
    rng = Random(18112026)
    for _ in range(2):
        seeds = [(rng.randrange(1, 65536, 2), rng.randrange(65536), width) for width in (8, 16, 16, 16)]
    slots = table.windows(65536, seeds)
    assert tuple(map(len, slots)) == (4, 7, 10, 10)
    joint = SeedProduct(table, slots)
    targets = list(table.index)
    results = joint.target_units(targets)
    # The independent reference builds the normalized polynomial from all 2,800 tuple forms.
    one = polynomial([1], table.p)
    leaves = [(form(row, table.p, table.a, table.c), [one, *(u*one for u in row)])
              for row in product(*slots)]
    direct, markers = product_tree(leaves)
    x = polynomial([0, 1], table.p)
    moduli = [(x-t)**max(h+1, 2) for t, (h, _, _) in zip(targets, results)]
    for t, (h, u, dz), rem in zip(targets, results, batch_remainders(direct, moduli)):
        assert u != 0
        assert rem == (u*(x-t)**h if h else u+dz*(x-t))
        if h:
            assert dz == (u if h == 1 else 0)
        if t == 0:
            gradient = joint.zero_gradient(t, h, u, dz)
            assert gradient == [43446, 16953, 13292, 47523, 12119, 9163]
            print(f"2,800-tuple fixture: order {h}, target unit {u}, additive cofactor {gradient[1]}, gradient {gradient}.")
    print(f"Checked all {len(targets)} target abscissas against the literal product by local polynomial remainders.")
    zeros = [(t, jet) for t, jet in zip(targets, results) if jet[0]]
    linear_moduli = [x-t for t, _ in zeros]
    expected = zip(*(evaluations(marker, linear_moduli) for marker in markers))
    for (t, (h, u, dz)), gradient in zip(zeros, expected):
        assert joint.zero_gradient(t, h, u, dz) == [dz, *gradient]
    print(f"All {len(zeros)} full gradients at zeros agreed with the five independent marker polynomials.")
    sign_rank(table, 0, 16, False)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    small()
    medium = table_for(1031, 9, 1019)
    check_scales(medium, False)
    sign_rank(medium, 0, 32, True)
    large()
