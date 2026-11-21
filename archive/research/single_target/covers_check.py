"""Check power-class and moment count reductions by literal guarded offsets on F_211."""

from collections import Counter
from itertools import product
from math import prod

from target import chord_scale, lift_windows, local_value, point_add
from arithmetic_alternating_chord import points


p, a, c, r, modulus = 211, 1, 1, 223, 256
degrees = (2, 3, 5, 7)
affine = points(p, a, c)
anchor = (0, 1)
assert len(affine)+1 == r


def neg(P):
    return P[0], -P[1] % p


def add(P, Q):
    return point_add(P, Q, p, a)


def offsets(body, endpoint):
    good = 0
    for D in affine:
        partial = D
        for P in (*body, neg(D)):
            if partial[0] == P[0]:
                break
            partial = add(partial, P)
        else:
            assert partial == endpoint
            good += 1
    return good


def subboxes(seeds):
    slack = 2*len(seeds)+11
    if slack*prod(2*w for _, _, w in seeds) < r:
        yield seeds
        return
    i = max(range(len(seeds)), key=lambda j: seeds[j][2])
    alpha, beta, width = seeds[i]
    assert width > 1
    for shift in (0, width//2):
        yield from subboxes((*seeds[:i], (alpha, beta-shift, width//2), *seeds[i+1:]))


def feature_trace(seeds, T, prefix):
    """All counts here are independent exhaustive reference implementations of graph pins."""
    b = p.bit_length()
    power_counts = {d: Counter(pow(y, d, p) for y in range(1, p)) for d in degrees}
    class_totals, moments, log_bits = Counter(), Counter(), Counter()
    for box in subboxes(seeds):
        slots = lift_windows(p, a, c, modulus, box)
        signed = [[Q for P in slot for Q in (P, neg(P))] for slot in slots]
        guarded, ordinary = Counter(), Counter()
        for row in product(*signed):
            partial = None
            path = row[:-1] if prefix else row
            for P in path:
                partial = add(partial, P)
            if partial is None:
                weight = 1 if prefix else -1
                body = (*path, anchor)
            else:
                weight = chord_scale(partial, row[-1], p) if prefix else (T[0]-partial[0]) % p or 1
                body = (*path, neg(partial), anchor)
            count = offsets(body, anchor)
            for d in degrees:
                for j in range(1, d):
                    witnesses = power_counts[d][weight*pow(2, -j, p) % p]
                    assert witnesses in (0, d)
                    if witnesses:
                        guarded[d, j] += witnesses*count
                        ordinary[d, j] += 1
            # Simulate the unique bounded exponent witness. No log oracle feeds the count.
            for e in range(p-1):
                if pow(2, e, p) == weight % p:
                    for k in range((p-2).bit_length()):
                        if e & (1 << k):
                            key = ('power-bit', k)
                            guarded[key] += count
                            ordinary[key] += 1
            if not prefix and partial == T:
                hits = offsets(row, T)
                bits = [(i, k) for i, P in enumerate(row) for k in range(b) if P[0] & (1 << k)]
                for key in [(), *[(x,) for x in bits], *[(x, y) for j, x in enumerate(bits) for y in bits[j+1:]]]:
                    guarded[key] += hits
                    ordinary[key] += 1
        for key in ordinary:
            if key and key[0] == 'power-bit':
                recovered = (guarded[key]+r-1)//r
                log_bits[key[1]] += recovered
            elif len(key) == 2 and isinstance(key[0], int):
                d, j = key
                assert guarded[key] % d == 0
                recovered = (guarded[key]//d+r-1)//r
                class_totals[key] += recovered
            else:
                recovered = (guarded[key]+r-1)//r
                moments[key] += recovered
            assert recovered == ordinary[key]
    exponents = {}
    for d in degrees:
        exponent = 0
        for j in range(1, d):
            tally = class_totals[d, j]
            symmetry = 4 if prefix else 2
            assert tally % symmetry == 0
            exponent += j*(tally//symmetry)
        exponents[d] = exponent % d
    symmetry = 4 if prefix else 2
    assert all(tally % symmetry == 0 for tally in log_bits.values())
    exponents[p-1] = sum((1 << k)*(tally//symmetry) for k, tally in log_bits.items()) % (p-1)
    first = second = 0
    for key, count in moments.items():
        if not key:
            continue
        weights = [p**i*(1 << k) for i, k in key]
        if len(key) == 1:
            first += count*weights[0]
            second += count*weights[0]**2
        else:
            second += 2*count*prod(weights)
    return exponents, (moments[()], first, second)


def check(seeds, T):
    slots = lift_windows(p, a, c, modulus, seeds)
    expected = local_value(p, a, T, slots)
    exponents, moments = feature_trace(seeds, T, False)
    assert moments == (expected.order, expected.first_moment, expected.second_moment)
    for j in range(2, len(seeds)+1):
        prefix, _ = feature_trace(seeds[:j], T, True)
        multiplicity = prod(2*len(slot) for slot in slots[j:])
        for d in (*degrees, p-1):
            exponents[d] = (exponents[d]+multiplicity*prefix[d]) % d
    # These four orders are coprime and their product is p-1 on this fixture.
    for d in degrees:
        assert pow(expected.unit, (p-1)//d, p) == pow(2, exponents[d]*((p-1)//d), p)
    exponent = sum(exponents[d]*((p-1)//d)*pow((p-1)//d, -1, d) for d in degrees) % (p-1)
    assert pow(2, exponent, p) == expected.unit
    assert exponent == exponents[p-1]
    return expected


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    P, Q = affine[0], affine[8]
    S = add(P, Q)
    atoms = [(alpha, -alpha*point[0] % modulus, 4)
             for alpha, point in zip((5, 15, 23), (P, Q, S))]
    fixtures = [(tuple(atoms[:2]), S), (tuple(atoms), anchor), ((atoms[0],)*3, P),
                ((atoms[0],)*2, add(P, P))]
    for seeds, T in fixtures:
        answer = check(seeds, T)
        print(f"Widths {[s[2] for s in seeds]}, target {T}: order {answer.order}, "
              f"unit {answer.unit} by CRT and direct exponent-bit counts; both integer moments matched.")
    P = anchor
    T = add(P, add(P, P))
    cases = [((1, 0, 2), (1, 205, 8)), ((1, 0, 2), (11, 185, 8))]
    results = [check(seeds, T) for seeds in cases]
    assert [v.unit for v in results] == [104, 200]
    for v in results:
        assert (v.order, v.first_moment, v.second_moment) == (1, 11183, 125059489)
        assert pow(v.unit, 7, p) == 156
    print("Two seeds at the same isolated target share both moments and the order-30 character, with units 104 and 200.")
