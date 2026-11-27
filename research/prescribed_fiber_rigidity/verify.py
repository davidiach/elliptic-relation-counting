"""Exact, finite checks of prescribed-sum fibers in abstract cyclic groups.

Run with Python 3.9+ and SymPy. Labels of group elements are input residues.
"""
from collections import defaultdict
from itertools import combinations, combinations_with_replacement, product
from sympy import Matrix, ones, zeros
from sympy.matrices.normalforms import hermite_normal_form


def modular_basis(rows, prime):
    pivots = {}
    for row in rows:
        v = [int(x) % prime for x in row]
        for j in range(len(v)):
            x = v[j]
            if not x:
                continue
            if j not in pivots:
                inv = pow(x, -1, prime)
                pivots[j] = [(inv * y) % prime for y in v]
                break
            v = [(a - x * b) % prime for a, b in zip(v, pivots[j])]
    return list(pivots.values())


def minus(a, b):
    return tuple(x - y for x, y in zip(a, b))


def half_tuples(blocks, prime):
    labels = [(i, x) for i, block in enumerate(blocks) for x in block[1:]]
    return [(sum(xs) % prime, tuple(int(xs[i] == x) for i, x in labels))
            for xs in product(*blocks)]


def triple_fibers(halves, prime):
    fibers = defaultdict(set)
    # Permutations of the three copies have the same incidence vector.
    for triple in combinations_with_replacement(halves, 3):
        g = sum(x[0] for x in triple) % prime
        v = tuple(map(sum, zip(*(x[1] for x in triple))))
        fibers[g].add(v)
    return {g: sorted(vectors) for g, vectors in fibers.items()}


def check_lattice(prime, blocks):
    halves = half_tuples(blocks, prime)
    assert 3 * len(set(g for g, _ in halves)) > 2 * prime + 1
    fibers = triple_fibers(halves, prime)
    assert len(fibers) == prime
    double = defaultdict(set)
    values = {g for g, _ in halves}
    for (x, vx), (y, vy) in combinations_with_replacement(halves, 2):
        double[(x+y) % prime].add(tuple(a+b for a, b in zip(vx, vy)))
    contrasts = {}
    for h, vectors in double.items():
        first, *rest = sorted(vectors)
        contrasts[h] = [minus(v, first) for v in rest]
    for t, points in fibers.items():
        rows = [minus(v, points[0]) for v in points[1:]]
        d = len(points[0])
        H = hermite_normal_form(Matrix(rows).T)
        assert H.shape == (d, d) and abs(H.det()) == prime
        assert len(modular_basis(rows, prime)) == d - 1
        shared = [v for z in values for v in contrasts[(t-z) % prime]]
        assert hermite_normal_form(Matrix(shared).T) == H
    return len(fibers)


def check_scores():
    total = 0
    for r in (3, 5, 7):
        for n in range(1, r + 1):
            if 3 * n <= 2 * r:
                continue
            for support in combinations(range(r), n):
                # Normalize the first value, then exhaust binary-valued scores over the integers.
                for tail in product((0, 1), repeat=n-1):
                    if not any(tail):
                        continue
                    f = dict(zip(support, (0,) + tail))
                    scores = defaultdict(set)
                    for x, y in combinations_with_replacement(support, 2):
                        scores[(x+y) % r].add(f[x]+f[y])
                    mixed = sum(len(values) > 1 for values in scores.values())
                    assert mixed >= 2*n-r-1
                    total += 1
    return total


def check_intervals():
    count = 0
    for r in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        for n in range((2*r)//3+1, r+1):
            values = defaultdict(set)
            for x, y in product(range(n), repeat=2):
                values[(x+y) % r].add(x+y)
            assert sum(len(v)>1 for v in values.values()) == 2*n-r-1
            count += 1
        n = (2*r+1)//3
        totals = {x+y+z for x, y, z in product(range(n), repeat=3)
                  if (x+y+z) % r == r-1}
        assert totals == {r-1}
    return count


def check_grams(r):
    F = list(range(1, r))
    by_pin = defaultdict(list)
    for x, y in product(F, repeat=2):
        z = (-x-y) % r
        if z:
            by_pin[z].append([int(x==i)+int(y==i)+int(z==i) for i in F])
    def gram(rows):
        M = Matrix(rows)
        v = M.T * ones(M.rows, 1)
        return 2 * (M.rows * M.T * M - v*v.T)
    all_rows = sum(by_pin.values(), [])
    K = gram(all_rows)
    J = ones(r-1)
    R = Matrix([[int((i+j)%r == 0) for j in F] for i in F])
    predicted = -24 * (Matrix.eye(r-1) + R) - 48 * J
    assert (K-predicted).applyfunc(lambda x: x%r) == zeros(r-1)
    rank = len(modular_basis(K.tolist(), r))
    stack = [row for rows in by_pin.values() for row in gram(rows).tolist()]
    pinned = len(modular_basis(stack, r))
    raw = len(modular_basis([minus(x, all_rows[0]) for x in all_rows], r))
    assert all(len(rows)==r-2 for rows in by_pin.values())
    assert (rank, pinned, raw) == ((r-3)//2, r-3, r-3)
    return rank, pinned, raw


def check_projected_hulls():
    r, blocks = 7, [(1, 2, 3), (2, 4, 6)]
    halves = half_tuples(blocks, r)
    assert len(set(g for g, _ in halves)) == r
    # Explicit genuine shift row, between half tuples of successive values.
    by_value = dict(halves)
    B = [minus(by_value[1], by_value[0])]
    fibers = triple_fibers(halves, r)
    origin = fibers[0][0]
    d = len(origin)-1
    no_single = []
    for t, points in fibers.items():
        cols = [tuple((x-y)%r for x, y in zip(v, origin)) for v in points]
        directions = [minus(v, cols[0]) for v in cols[1:]]
        assert len(modular_basis(B+directions, r)) == len(origin)
        chosen = [cols[0]]
        basis = list(B)
        for v in cols[1:]:
            augmented = modular_basis(basis+[minus(v, cols[0])], r)
            if len(augmented) > len(basis):
                chosen.append(v)
                basis = augmented
        assert len(chosen) == d+1
        # Exhaust affine coefficients on these small fixtures, including signed residues.
        span = {tuple(a*z%r for z in B[0]) for a in range(r)}
        def value(weights):
            return tuple(sum(a*z for a,z in zip(weights, col))%r for col in zip(*chosen))
        assert any(value((1-sum(tail),)+tail) in span for tail in product(range(r), repeat=d))
        if all(v not in span for v in cols):
            no_single.append(t)
    return d, sorted(no_single)


def check_pin_cells_and_cuts():
    pin_cases = 0
    for r in (5, 7):
        for n in range((2*r+1)//3+1, r):
            for S in combinations(range(r), n):
                unit = {x: tuple(int(x==a) for a in S) for x in S}
                for t in range(r):
                    stack = []
                    for z in S:
                        counts = []
                        for x in S:
                            y = (t-z-x) % r
                            if y in unit:
                                counts.append([a+b for a,b in zip(unit[x], unit[y])])
                        C = Matrix(counts)
                        v = C.T * ones(C.rows,1)
                        stack.extend((C.rows*C.T*C-v*v.T).tolist())
                    assert len(modular_basis(stack,r)) == n-2
                    pin_cases += 1
    r = 7
    halves = half_tuples([(1,2,3), (2,4,6)],r)
    # This coordinate does not factor through an affine function of group value.
    sample = [(g,v[0]) for g,v in halves]
    tested = 0
    for t, a in product(range(r), repeat=2):
        pairs = [((g+h)%r, f+k) for (g,f),(h,k) in product(sample,repeat=2)]
        rows = [[(f+k-a)%r if (g+h)%r==t else 0 for h,k in pairs] for g,f in sample]
        rank = len(modular_basis(rows,r))
        block_sum = 0
        for x in range(r):
            U = {f%r for g,f in sample if g==x}
            V = {k%r for h,k in pairs if h==(t-x)%r}
            if U and V:
                if len(U)>1 and len(V)>1:
                    block_sum += 2
                else:
                    block_sum += int(any((u+v-a)%r for u in U for v in V))
        assert rank == block_sum and rank >= r-1
        tested += 1
    return pin_cases,tested


if __name__ == '__main__':
    print('Binary integer-score tables checked:', check_scores())
    print('Sharp interval counts checked:', check_intervals())
    cases = 0
    for r in (3, 5, 7):
        for n in range(1, r+1):
            if 3*n > 2*r+1:
                for S in combinations(range(r), n):
                    cases += check_lattice(r, [S])
    print('Full and common-third integral HNFs and modular rank, all dense one-block fibers on C3,C5,C7:', cases)
    cases = sum(check_lattice(r, blocks) for r, blocks in [
        (7, [(1,2,3), (2,4,6)]), (11, [(1,2), (1,3,5,7,9)]),
        (13, [(1,2,3), (1,4), (1,7)])])
    print('Full and common-third integral HNFs and modular rank, multiblock fibers:', cases)
    for r in (5, 7, 11, 13, 17):
        print('C%d unpinned, third-slot-pinned, raw ranks:' % r, check_grams(r))
    print('Projected affine hulls on all seven targets, deficiency and targets lacking a single certificate:', check_projected_hulls())
    print('Exhaustive third-slot stacks and discrepancy cut ranks:', check_pin_cells_and_cuts())
