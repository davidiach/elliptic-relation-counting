"""Exact odd-divisor counts for additive quadruples of signed moment sections."""
from collections import Counter
from itertools import product
from random import Random

from flint import nmod_mat, nmod_poly

from check_address_parity import transform
from fixtures import evaluation
from probe_signed_moment import Probe, elimination


def polynomial_key(poly):
    return tuple(map(int, poly.coeffs()))


def odd_divisor(a, b, cubic):
    """Odd geometric places of a(x)+b(x)y, encoded by Frobenius-orbit bits."""
    if not a and not b:
        return None
    pole = max(2*a.degree() if a else 0, 2*b.degree()+3 if b else 0)
    places = {"infinity"} if pole % 2 else set()
    content = a.gcd(b)
    a, b = a//content, b//content
    for factor, power in content.factor()[1]:
        if power % 2 and cubic % factor:
            places.add(("both", polynomial_key(factor)))
    norm = a*a-cubic*b*b
    for factor, power in norm.factor()[1]:
        if not power % 2:
            continue
        key = polynomial_key(factor)
        if not cubic % factor:
            places.add(("branch", key))
        else:
            # Exactly one of the two conjugate y values is a zero of the primitive part.
            unit, inverse, _ = b.xgcd(factor)
            assert unit == 1
            y = (-a*inverse) % factor
            assert (y*y-cubic) % factor == 0 and y != 0
            places.symmetric_difference_update({("low", key)})
            if polynomial_key(y) > polynomial_key((-y) % factor):
                places.symmetric_difference_update({("both", key)})
    return frozenset(places)


def check_divisors(p, cubic):
    x = nmod_poly([0, 1], p)
    zero = nmod_poly([], p)
    one = nmod_poly([1], p)
    expressions = [(x, zero), (one, one), (zero, one), (x*x+1, x), (-one, one)]
    assert odd_divisor(zero, zero, cubic) is None
    for a, b in expressions:
        for c, d in expressions:
            left = odd_divisor(a, b, cubic)
            right = odd_divisor(c, d, cubic)
            assert odd_divisor(a*c+cubic*b*d, a*d+b*c, cubic) == left ^ right


def fixture(m):
    p = 1031
    probe = Probe(p, 9, (0,))
    rng = Random(90201+m)
    # Draw actual mixed masks until the root has m disjoint full four-point slots.
    # This selects one small dense geometry for each arity, not a fresh-law statistic.
    slots, used = [], set()
    while len(slots) < m:
        mul = 1+2*rng.randrange(probe.modulus//2) if not slots else 1
        offset = rng.randrange(probe.modulus) if not slots else -2*rng.randrange(probe.modulus//2)
        mask = []
        for j in range(2):
            x = ((j-offset)*pow(mul, -1, probe.modulus)) % probe.modulus
            for index, sign in probe.by_x.get(x, ()):
                mask.append((index, 2*j+sign))
        if len(mask) == 4 and not used.intersection(i for i, _ in mask):
            slots.append(mask)
            used.update(i for i, _ in mask)
    d = m+1
    dim = 2*d
    cube = 4**m
    points = [[index for index, _ in slot] for slot in slots]
    addresses = [[a for _, a in slot] for slot in slots]
    entries = []
    for indices in product(range(4), repeat=m):
        address = 0
        rows = []
        for i, choice in enumerate(indices):
            address = 4*address+addresses[i][choice]
            rows.append(evaluation(probe.curve[points[i][choice]], d, p))
        # Determinant linear form in the target evaluation row.
        cofactors = [(-1)**(m+j)*elimination([row[:j]+row[j+1:] for row in rows], p)
                     for j in range(d)]
        entries.append((address, cofactors))

    def observations(i):
        endpoint = evaluation(probe.curve[i], d, p)
        determinants = [0]*cube
        for address, cofactor in entries:
            determinants[address] = sum(c*x for c, x in zip(cofactor, endpoint)) % p
        first = transform(determinants)
        second = transform([v*v % p for v in determinants])
        return [v % p for v in second], [(v+u*u) % p for u, v in zip(first, second)]

    sample_rows, training = [], []
    for i in range(1, probe.n):
        row = evaluation(probe.curve[i], dim, p)
        if nmod_mat(sample_rows+[row], p).rank() > len(sample_rows):
            sample_rows.append(row)
            training.append(i)
        if len(training) == dim:
            break
    assert len(training) == dim
    inverse = nmod_mat(sample_rows, p).inv()
    training_values = [observations(i) for i in training]
    cubic = nmod_poly([1, 9, 0, 1], p)
    check_divisors(p, cubic)
    for kind, name in enumerate(("M2", "M2+M1^2")):
        matrix = inverse*nmod_mat([v[kind] for v in training_values], p)
        ids, tags = {}, []
        even = d+1
        for r in range(cube):
            a = nmod_poly([int(matrix[j, r]) for j in range(even)], p)
            b = nmod_poly([int(matrix[j, r]) for j in range(even, dim)], p)
            places = odd_divisor(a, b, cubic)
            if places is None:
                tags.append(None)
                continue
            fingerprint = 0
            for place in places:
                fingerprint |= 1 << ids.setdefault(place, len(ids))
            tags.append(fingerprint)
        exceptional = nonzero = 0
        for h in range(cube):
            buckets = Counter(tags[r] ^ tags[r ^ h] for r in range(cube)
                              if tags[r] is not None and tags[r ^ h] is not None)
            exceptional += sum(v*v for v in buckets.values())
            nonzero += sum(buckets.values())**2
        # Check all interpolants at two additional public targets using fresh direct row sums.
        for i in (probe.n//3, probe.n//2):
            calculated = nmod_mat([evaluation(probe.curve[i], dim, p)], p)*matrix
            assert [int(calculated[0, r]) for r in range(cube)] == observations(i)[kind]
        active = sum(t is not None for t in tags)
        diagonals = 3*active**2-2*active
        assert 0 < exceptional == diagonals <= nonzero <= cube**3
        field_values = nmod_mat([evaluation(pt, dim, p) for pt in probe.curve[1:]], p)*matrix
        fourth = 0
        for i in range(probe.n-1):
            coefficients = transform([probe.legendre[int(field_values[i, r])] for r in range(cube)])
            fourth += sum(c**4 for c in coefficients)
        print(f"p={p}, m={m}, {name}, {cube} masks: {nonzero:,} nonzero products, "
              f"{exceptional:,} with even divisor (all diagonal), out of {cube**3:,} additive "
              f"quadruples; target-averaged fourth energy {fourth/((probe.n-1)*cube**4):.6f}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    fixture(3)
    fixture(4)
