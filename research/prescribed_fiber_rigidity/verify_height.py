"""Check the small augmented lift and two-witness certificates in a public C101 fixture.

The dynamic program is for finite verification. The theorem uses two-pin prefix
counts for asymptotic costs. No elliptic curve arithmetic or keys are used here.
"""
import json
from collections import defaultdict
from pathlib import Path

from sympy import Matrix, ones
from verify import modular_basis


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def kernel_line(rows, prime):
    basis = modular_basis(rows, prime)
    pivots = {next(j for j, x in enumerate(row) if x): row for row in basis}
    free, = set(range(len(rows[0]))) - pivots.keys()
    u = [0]*len(rows[0])
    u[free] = 1
    for j in sorted(pivots, reverse=True):
        u[j] = -dot(pivots[j], u) % prime
    assert all(dot(row, u) % prime == 0 for row in rows)
    return u


def distribution(blocks, labels, copies, prime, modular_score=False):
    # Retain multiplicities and one witness path per attainable (value, score).
    states = {(0, 0): (1, ())}
    for i, block in enumerate(blocks*copies):
        nxt = {}
        for (g, s), (n, witness) in states.items():
            for j, a in enumerate(block):
                score = s + labels[i % len(blocks)][j]
                key = ((g+a) % prime, score % prime if modular_score else score)
                if key in nxt:
                    old_count, path = nxt[key]
                    nxt[key] = old_count+n, path
                else:
                    nxt[key] = n, witness+(j,)
        states = nxt
    return states


def run():
    data = json.loads(Path(__file__).with_name('pointwise_fixture.json').read_text())
    r, blocks = data['prime'], data['blocks']
    coords = [(i, j) for i, block in enumerate(blocks) for j in range(1, len(block))]
    d = len(coords)
    rows = []
    for description in data['rows']:
        row = [0]*d
        for j, a in description:
            row[j] = a
        # Reinsert reference changes. Pad source and successor by identical
        # references to get exactly two entries per block on both sides.
        for i, block in enumerate(blocks):
            part = [row[k] for k, (b, _) in enumerate(coords) if b == i]
            part.append(-sum(part))
            assert sum(x for x in part if x > 0) <= 2
        assert sum(x*(blocks[i][j]-blocks[i][0]) for x, (i, j) in zip(row, coords)) % r == 1
        rows.append(row)

    matrix = Matrix(rows)
    assert abs(matrix.det()) % r == 0
    lifted = matrix.inv()*ones(d, 1)
    assert all(x.q == 1 for x in lifted)
    v = list(map(int, lifted))
    u = kernel_line(rows, r)

    def indexed(x):
        labels = [[0]*len(block) for block in blocks]
        for value, (i, j) in zip(x, coords):
            labels[i][j] = value
        return labels

    vslots = indexed(v)
    width = sum(max(arr)-min(arr) for arr in vslots)
    assert 0 < 2*width < r
    signature = distribution(blocks, indexed(u), 2, r, modular_score=True)
    by_signature = defaultdict(set)
    for g, score in signature:
        by_signature[score].add(g)
    single = {(x-y) % r for group in by_signature.values() for x in group for y in group}
    assert len(single) == 49

    half = distribution(blocks, vslots, 1, r)
    double = distribution(blocks, vslots, 2, r)
    options = defaultdict(list)
    for (g, score), (count, path) in double.items():
        options[g].append((score, count, path))
    moments = {}
    for g, entries in options.items():
        n = sum(c for _, c, _ in entries)
        first = sum(s*c for s, c, _ in entries)
        second = sum(s*s*c for s, c, _ in entries)
        moments[g] = 2*(n*second-first*first)

    def incidence(path):
        return tuple(sum(x == j for x in path[i::len(blocks)]) for i, j in coords)

    target_pairs = {}
    for t in range(r):
        total = sum(c*moments[(t-z) % r] for (z, _), (c, _) in half.items())
        assert total > 0
        for (z, _), (_, common) in half.items():
            entries = options[(t-z) % r]
            if len(entries) < 2:
                continue
            lo = min(entries)[2]+common
            hi = max(entries)[2]+common
            a, b = incidence(lo), incidence(hi)
            assert all(sum(blocks[i % len(blocks)][j] for i, j in enumerate(path)) % r == t
                       for path in (lo, hi))
            assert 0 < dot(b, v)-dot(a, v) <= 2*width
            target_pairs[t] = a, b
            break
    origin = target_pairs[0][0]
    for t, (a, b) in target_pairs.items():
        x, y = ([j-i for j, i in zip(target, origin)] for target in (a, b))
        dx, dy = dot(x, u) % r, dot(y, u) % r
        alpha = -dy * pow(dx-dy, -1, r) % r
        certificate = [(alpha*j+(1-alpha)*k) % r for j, k in zip(x, y)]
        assert dot(certificate, u) % r == 0
        # Bv=1, so pairing with the rationally recovered lift transports RHS.
        assert dot(certificate, v) % r == t
    assert len(target_pairs) == r
    print(f"C{r}: genuine rational basis {d}/{d}, modular rank {d-1}/{d}.")
    print(f"Recovered integer solution range [{min(v)}, {max(v)}]; sum of block widths {width}.")
    print(f"Identifiable individual doubled-pair displacements: {len(single)}/{r}.")
    print(f"Two common-third target triples certified all {len(target_pairs)} targets against one identity triple.")


if __name__ == '__main__':
    run()
