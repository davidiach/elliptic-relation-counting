"""Exact small-field rank distribution for privately reweighted signed moments."""
from collections import Counter
from itertools import product

from check_address_parity import character
from fixtures import determinant
from probe_signed_moment import Probe


def fixture(probe, k):
    short = [[(point, sign) for point, sign in probe.by_x[x]]
             for x in range(16) if len(probe.by_x[x]) == 2][:2]
    for start in range(16, probe.p-15, 16):
        wide = [(point, 2*offset+sign) for offset in range(16)
                for point, sign in probe.by_x[start+offset]]
        if len(wide) <= 4*k:
            continue
        slots = [*short, wide]
        abscissas = {probe.curve[i][0] for slot in slots for i, _ in slot}
        for target in range(1, probe.n):
            if probe.curve[target][0] in abscissas:
                continue
            hits = [(a, b, c) for (a,_), (b,_), (c,_) in product(*slots)
                    if (a+b+c) % probe.n == target]
            if len(hits) == 1:
                a, b, c = hits[0]
                partial = (a+b) % probe.n
                if partial and probe.curve[partial][0] != probe.curve[c][0]:
                    return slots, target, hits[0]
    raise AssertionError("The fixture must contain a regular singleton.")


def row_rank(rows, p):
    matrix = [list(row) for row in rows]
    pivots = []
    for col in range(len(matrix[0])):
        j = len(pivots)
        i = next((i for i in range(j, len(matrix)) if matrix[i][col] % p), None)
        if i is None:
            continue
        matrix[j], matrix[i] = matrix[i], matrix[j]
        pivot = matrix[j]
        inverse = pow(pivot[col] % p, -1, p)
        for i in range(j+1, len(matrix)):
            scale = matrix[i][col]*inverse % p
            if scale:
                matrix[i] = [(x-scale*y) % p for x, y in zip(matrix[i], pivot)]
        pivots.append(col)
        if len(pivots) == len(matrix):
            break
    return pivots


def main():
    probe = Probe(101, 17, (0,))
    p, k = probe.p, 3
    (first, second, wide), target, witness = fixture(probe, k)
    endpoint = probe.vectors[-target]
    minors, powers = [], []
    for (a, _), (b, _) in product(first, second):
        rows = [probe.vectors[a], probe.vectors[b], endpoint]
        # The d=4 pole coordinate is x^2, at column index 2.
        leading = determinant([[v for j, v in enumerate(row) if j != 2] for row in rows], p)
        if (a, b) == witness[:2]:
            assert leading
        minors.append([pow(leading, j, p) for j in range(k+1)])
        delta = [determinant(rows[:2]+[probe.vectors[c], endpoint], p) for c, _ in wide]
        powers.append([[pow(x, j, p) for x in delta] for j in range(k+1)])
    # The uniform distribution of a pair of weights is a mass-one zero
    # vector plus p+1 projective directions of equal multiplicity p-1.
    directions = [((0,0), 1), ((0,1), p-1)] + [((1,t), p-1) for t in range(p)]
    assert sum(mass for _, mass in directions) == p*p
    ranks = Counter()
    bad_leads = 0
    test_projection = True
    for (a, mass_a), (b, mass_b) in product(directions, repeat=2):
        weights = [v*w % p for v, w in product(a, b)]
        # Arbitrary fixed Walsh signs in the two short slots are absorbed
        # by the uniform projective enumeration. Retain them in the wide slot.
        leading = [sum(w*row[j] for w, row in zip(weights, minors)) % p for j in range(k+1)]
        rows = [[character(address, 13)*sum(w*mat[j][i] for w, mat in zip(weights, powers)) % p
                 for i, (_, address) in enumerate(wide)] for j in range(k+1)]
        rank = len(row_rank(rows, p))
        ranks[rank] += mass_a*mass_b
        if all(leading):
            assert rank == k+1
        else:
            bad_leads += mass_a*mass_b
        if test_projection and rank == k+1:
            i, j = row_rank(rows[1:3], p)
            image = {tuple((v*row[i]+w*row[j]) % p for row in rows[1:3])
                     for v, w in product(range(p), repeat=2)}
            assert len(image) == p*p
            test_projection = False
    total = p**4
    assert sum(ranks.values()) == total and 0 < total-ranks[k+1] <= bad_leads
    assert p*bad_leads <= (k+1)*2*total and not test_projection
    print(f"p={p}, m=3, target=[{target}]G, regular singleton, wide slot {len(wide)} points.")
    print(f"  Of {total:,} weight choices in the first two slots, {ranks[k+1]:,} give rank four "
          f"for orders 0..3; {bad_leads:,} have a missing leading coefficient.")
    print(f"  Two varying final-slot weights produced all {p*p:,} pairs of orders one and two.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
