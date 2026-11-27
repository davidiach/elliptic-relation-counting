"""Run from the repository root with python -B checks/check_fixed_shift_points.py."""
from itertools import combinations, product
from math import prod

import numpy as np

from check_batch import add
from check_fixed_shift_gram import rank_and_det
from fixed_shift_points import RowSpace, compress_stack, explicit_row_points, sparse_witness_points
from fixtures import curve_cycle


def elliptic_rows():
    cycle = curve_cycle(13, 4)
    labels = [1, 2, 3, 1, 5, 9, 1, 12]
    n = len(labels)
    tuples = []
    for digits in product(range(3), repeat=4):
        slots = tuple(d + (3 if j >= 2 else 0) for j, d in enumerate(digits))
        value = None
        for j in slots:
            value = add(value, cycle[labels[j]], 13, 4)
        address = sum(d << (2*j) for j, d in enumerate(digits))
        tuples.append((value, address, [slots.count(j) for j in range(n)]))
    by_value = {p: next(t for t in tuples if t[0] == p) for p in cycle}
    explicit = [[int(j == n-2) for j in range(n)]]
    for j, log in enumerate(labels):
        counts = by_value[cycle[log]][2]
        explicit.append([int(k == j) - v for k, v in enumerate(counts)])
    shifts = []
    bins = {p: [] for p in cycle}
    for record in tuples:
        bins[record[0]].append(record)
    for p, x, u in tuples:
        successor = add(p, cycle[1], 13, 4)
        for _, y, v in bins[successor]:
            shifts.append(((x << 8) + y, [b-a for a, b in zip(u, v)]))
    return explicit, shifts, labels


def gram_at(point, explicit, shifts, prime, radix):
    n = len(explicit[0])
    def monomial(key):
        value = 1
        for x in point:
            key, digit = divmod(key, radix)
            value = value * pow(x, digit, prime) % prime
        assert key == 0
        return value
    c = len(explicit)
    offset = monomial(c)
    rows = [(monomial(e), row) for e, row in enumerate(explicit)]
    rows.extend((offset * monomial(v) % prime, row) for v, row in shifts)
    return [[sum(w * a[j] * a[k] for w, a in rows) % prime
             for k in range(n)] for j in range(n)]


def verify_elliptic():
    explicit, shifts, labels = elliptic_rows()
    prime, radix, s = 19, 2, 16
    n, c = len(labels), len(explicit)
    points = list(explicit_row_points(prime, radix, s, c))
    def matrix(x):
        return gram_at(x, explicit, shifts, prime, radix)
    aggregate = RowSpace(prime, n)
    for z, x in enumerate(points):
        observed = [sum(col) % prime for col in zip(*matrix(x))]
        expected = [sum(pow(z, e, prime)*sum(a)*a[j] for e, a in enumerate(explicit))
                    % prime for j in range(n)]
        assert observed == expected
        aggregate.add(observed)
    required, _ = rank_and_det(explicit, prime)
    assert aggregate.rank == required
    assert not any(aggregate.add(a) for a in explicit)
    selected, space = compress_stack(points, matrix, prime, n)
    assert space.rank == n
    shift_space = RowSpace(prime, n)
    for _, a in shifts:
        shift_space.add(a)
    assert shift_space.rank == 6-2
    assert all(sum(a[j]*labels[j] for j in range(n)) % prime == 1 for _, a in shifts)
    print(f"E(F13): {len(shifts)} shifts; explicit projection rank {required}; "
          f"shift-row rank {shift_space.rank}; original stack rank {space.rank} "
          f"with {len(selected)} retained points out of {c}.")

    # Find singular matrices with complementary row spaces on this same curve.
    singular = []
    # Give six consecutive address bits the powers of one field element.
    # Report the three parameters; the evaluator uses the actual 16-digit map.
    for x in product(range(5), repeat=3):
        fine_point = tuple(pow(x[j//6], 1 << (j % 6), prime) for j in range(s))
        b = matrix(fine_point)
        rank, _ = rank_and_det(b, prime)
        if rank and rank < n:
            singular.append((x, b))
    for (x, a), (y, b) in combinations(singular, 2):
        if rank_and_det(a+b, prime)[0] == n:
            ra, rb = rank_and_det(a, prime)[0], rank_and_det(b, prime)[0]
            print(f"Two singular evaluations (six-bit groups) on the same curve: {x} rank {ra}, "
                  f"{y} rank {rb}, stacked rank {n}.")
            return
    raise AssertionError("expected a complementary pair of singular evaluations")


def verify_affine_fibers():
    cycle = curve_cycle(13, 4)
    r, threshold, slot_bits = 19, 3, 3
    blocks = []
    for alpha, beta in ((1, 0), (2, 7)):
        block = []
        for log, (x, y) in enumerate(cycle[1:], 1):
            value = (alpha*x+beta) % 13
            if value < threshold:
                block.append((log, value+threshold*(y % 2)))
        assert len({address for _, address in block}) == len(block)
        blocks.append(block)
    a, b = blocks
    sums = {(x+y) % r for x, _ in a for y, _ in b}
    assert 8*len(sums) > 7*r
    labels = [x for block in blocks for x, _ in block] + [1, 12]
    n = len(labels)
    fibers = [[] for _ in range(r)]
    for indices in product(range(len(a)), range(len(a)), range(len(b)), range(len(b))):
        columns = (indices[0], indices[1], len(a)+indices[2], len(a)+indices[3])
        slot_addresses = (a[indices[0]][1], a[indices[1]][1],
                          b[indices[2]][1], b[indices[3]][1])
        address = sum(v << (slot_bits*j) for j, v in enumerate(slot_addresses))
        group_index = sum(labels[j] for j in columns) % r
        fibers[group_index].append((address, [columns.count(j) for j in range(n)]))
    explicit = [[int(j == n-2) for j in range(n)]]
    for j, log in enumerate(labels):
        counts = fibers[log][0][1]
        explicit.append([int(k == j)-v for k, v in enumerate(counts)])
    c = len(explicit)
    shift = 1 << (4*slot_bits)
    def evaluate(point):
        z = point[0]  # Called on the Kronecker curve, giving literal key powers.
        gram = np.zeros((n, n), dtype=np.int64)
        for e, row in enumerate(explicit):
            gram += pow(z, e, r)*np.outer(row, row)
        def moments(fiber, spacing):
            addresses, rows = zip(*fiber)
            rows = np.array(rows, dtype=np.int64)
            weights = np.array([pow(z, key*spacing, r) for key in addresses])
            return int(weights.sum()) % r, weights @ rows % r, (rows.T*weights) @ rows % r
        # Local moments reduce modulo 19. The bounded arrays never approach
        # the signed 64-bit limit (each fiber contains at most 576 tuples).
        for g, fiber in enumerate(fibers):
            nx, vx, hx = moments(fiber, shift)
            ny, vy, hy = moments(fibers[(g+1) % r], 1)
            gram += pow(z, c, r)*(nx*hy+ny*hx-np.outer(vx, vy)-np.outer(vy, vx))
        return (gram % r).tolist()
    selected, stack = compress_stack(explicit_row_points(r, 2, 8*slot_bits, c),
                                     evaluate, r, n)
    assert stack.rank == n and rank_and_det(explicit, r)[0] == 10
    # Independent check of the moment expansion at a nontrivial scalar.
    flat_rows = []
    flat_weights = []
    for e, row in enumerate(explicit):
        flat_rows.append(row)
        flat_weights.append(pow(2, e, r))
    for g, fiber in enumerate(fibers):
        for x, left in fiber:
            for y, right in fibers[(g+1) % r]:
                flat_rows.append([v-u for u, v in zip(left, right)])
                flat_weights.append(pow(2, c+shift*x+y, r))
    flat = np.array(flat_rows, dtype=np.int64)
    direct = (flat.T*flat_weights) @ flat % r
    assert evaluate((2,)) == direct.tolist()
    print(f"Parity-affine curve blocks: sizes {len(a)} and {len(b)}, "
          f"{len(sums)}/{r} half sums, explicit rank {rank_and_det(explicit, r)[0]}, "
          f"fiber-sum stack rank {stack.rank}/{n} at {len(selected)} retained points.")


def verify_sparse_list():
    prime, generator = 257, 3
    assert pow(generator, 128, prime) == prime-1
    points = list(sparse_witness_points(prime, generator, 2, 3, 10))
    exponents = list(product(range(4), repeat=2))
    columns = [[prod(pow(x, a, prime) for x, a in zip(point, index)) % prime
                for index in exponents] for point in points]
    checked = 0
    for subset in combinations(range(len(exponents)), 3):
        projected = [[row[i] for i in subset] for row in columns]
        assert rank_and_det(projected, prime)[0] == 3
        checked += 1
    print(f"Sparse list: {len(points)} generated points, {checked} three-term "
          "supports tested at full column rank over F257.")
    # A degree-nine interval in each variable is dense, but two geometric
    # difference factors give a product with just four monomials.
    dense_points = list(sparse_witness_points(prime, generator, 2, 4, 10))
    for x, y in dense_points:
        value = sum(x**i for i in range(10))*sum(y**i for i in range(10))
        assert (value*(1-x)*(1-y) - (1-x**10)*(1-y**10)) % prime == 0
    assert any((1-pow(x, 10, prime))*(1-pow(y, 10, prime)) % prime
               for x, y in dense_points)
    print("Dense 100-term input: its four-term polynomial multiple was hit.")


def verify_positive_obstruction():
    prime = 5
    points = [(1, 1), (2, 3)]
    monomials = list(product(range(3), repeat=2))
    for size in range(1, len(monomials)+1):
        for subset in combinations(monomials, size):
            if all(sum(x**i*y**j for i, j in subset) % prime == 0 for x, y in points):
                print(f"Zero-sum control: {size} distinct unit-coefficient monomials "
                      f"vanish at both points over F5: {subset}.")
                return
    raise AssertionError("degree threshold guarantees a nonempty subset")


def verify_dense_fibers():
    prime = 19
    blocks = ((1, 2, 3, 4, 5), (1, 4, 7, 10, 13))
    labels = sum(blocks, ())
    half = [(a+b, (i, 5+j)) for i, a in enumerate(blocks[0])
            for j, b in enumerate(blocks[1])]
    d = len({p % prime for p, _ in half})
    q = 2*d-prime
    assert 8*d > 7*prime
    fibers = [[] for _ in range(prime)]
    for (x, slots), (y, other) in product(half, repeat=2):
        indices = slots+other
        fibers[(x+y) % prime].append([indices.count(j) for j in range(len(labels))])
    collisions = RowSpace(prime, len(labels))
    for fiber in fibers:
        assert len(fiber) >= q
        for a in fiber:
            row = [x-y for x, y in zip(a, fiber[0])]
            assert sum(x*y for x, y in zip(row, labels)) % prime == 0
            collisions.add(row)
    assert collisions.rank == len(labels)-len(blocks)-1

    # All binary label directions, and their translates by the true logarithm.
    # Histograms count all valid shifts exactly without enumerating their pairs.
    from collections import Counter
    shift_pairs = sum(len(fibers[g])*len(fibers[(g+1) % prime]) for g in range(prime))
    checked = 0
    minimum = shift_pairs
    for direction in product(range(2), repeat=len(labels)):
        if len(set(direction[:5])) == len(set(direction[5:])) == 1:
            continue
        hist = [Counter(sum(x*y for x, y in zip(a, direction)) % prime for a in f)
                for f in fibers]
        for slope in (0, 1):
            good = sum(mass*hist[(g+1) % prime].get((v-slope) % prime, 0)
                       for g, counts in enumerate(hist) for v, mass in counts.items())
            bad = shift_pairs-good
            assert 2*bad >= q*q
            minimum = min(minimum, bad)
            checked += 1
    print(f"Dense half sums: {d}/{prime} residues; collision-row rank {collisions.rank}; "
          f"{checked} directions checked, minimum {minimum}/{shift_pairs} violating shifts.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without Python optimization.")
    verify_elliptic()
    verify_affine_fibers()
    verify_sparse_list()
    verify_positive_obstruction()
    verify_dense_fibers()
