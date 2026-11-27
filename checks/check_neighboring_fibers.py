"""Exact checks of the full-grid fallback and the displacement sum."""
from itertools import product
from math import prod

import numpy as np

from check_batch import add
from check_fixed_shift_gram import rank_and_det
from fixed_shift_points import RowSpace


def curve_cycle(p, a, b):
    assert (4*a**3 + 27*b**2) % p
    affine = [(x, y) for x, y in product(range(p), repeat=2)
              if (y*y-x**3-a*x-b) % p == 0]
    generator = affine[0]
    cycle = [None]
    point = generator
    while point is not None:
        assert point not in cycle
        cycle.append(point)
        point = add(point, generator, p, a)
    assert set(cycle[1:]) == set(affine)
    return cycle


def full_tuples(blocks, cycle):
    """Each block contains (point, injective address) pairs; two slots per block."""
    points = [p for block in blocks for p, _ in block]
    n = len(points)+2
    slot_bits = max(1, max(a for block in blocks for _, a in block).bit_length())
    tuples = [[] for _ in cycle]
    r = len(cycle)
    for indices in product(*(range(len(block)) for block in blocks for _ in range(2))):
        counts = [0]*n
        address, value, offset = 0, 0, 0
        for i, selected in enumerate(indices):
            block = blocks[i//2]
            point, digit = block[selected]
            counts[offset+selected] += 1
            address = (address << slot_bits) + digit
            value += cycle.index(point)
            if i % 2:
                offset += len(block)
        tuples[value % r].append((address, counts))
    return tuples, 2*len(blocks)*slot_bits


def small_group_explicit(cycle, target):
    # LZR's decomposition of A is (A,O,O,O). Use O first, then affine
    # lexicographic order within each of the two equal full-group blocks.
    ordered = [None]+sorted(cycle[1:])
    r = len(cycle)
    labels = ordered+ordered+[cycle[1], target]
    rows = [[int(i == len(labels)-2) for i in range(len(labels))]]
    for i, point in enumerate(labels):
        row = [int(i == j) for j in range(len(labels))]
        row[ordered.index(point)] -= 1
        row[0] -= 1
        row[r] -= 2
        rows.append(row)
    return rows


def binary_monomial(point, key, r):
    return prod(point[j] for j in range(len(point)) if key & (1 << j)) % r


def moments(tuples, point, r):
    out = []
    for fiber in tuples:
        counts = np.array([row for _, row in fiber], dtype=np.int64)
        weights = np.array([binary_monomial(point, key, r) for key, _ in fiber])
        out.append((int(weights.sum()) % r, weights @ counts % r,
                    (counts.T*weights) @ counts % r))
    return out


def displacement(left, right, g, r):
    matrix = np.zeros_like(left[0][2])
    for i, (n, v, h) in enumerate(left):
        m, w, k = right[(i+g) % r]
        matrix += n*k+m*h-np.outer(v, w)-np.outer(w, v)
    return matrix % r


def matrix_at(tuples, width, explicit, point, r):
    c = len(explicit)
    left = moments(tuples, point[width:], r)
    right = moments(tuples, point[:width], r)
    matrix = binary_monomial(point, c, r)*displacement(left, right, 1, r)
    for e, row in enumerate(explicit):
        matrix += binary_monomial(point, e, r)*np.outer(row, row)
    return (matrix % r).tolist()


def raw_points(r, s, c, generator):
    # The literal displayed formulas, even when their theorem's guard fails.
    slopes = 1+(s-1)*c*(c-1)//2
    points = {tuple(pow(z, 1 << j, r) for j in range(s)) for z in range(c)}
    for t in range(1, min(slopes, r-1)+1):
        for h in range(min(c, r-1)):
            points.add(tuple(pow(generator, h*(t**j), r) for j in range(s)))
    return sorted(points)


def finite_grids():
    for a, b, r, generator in ((2, 0, 2, 1), (4, 2, 3, 2), (3, 2, 5, 2), (2, 1, 7, 3)):
        cycle = curve_cycle(5, a, b)
        assert len(cycle) == r
        assert len({pow(generator, h, r) for h in range(r-1)}) == r-1
        addresses = list(zip([None]+sorted(cycle[1:]), range(r)))
        tuples, width = full_tuples([addresses, addresses], cycle)
        ranks = []
        for target in cycle:
            explicit = small_group_explicit(cycle, target)
            n, c, s = len(explicit[0]), len(explicit), 2*width
            raw = raw_points(r, s, c, generator)
            space = RowSpace(r, n)
            for point in raw:
                for row in matrix_at(tuples, width, explicit, point, r):
                    space.add(row)
            ranks.append(space.rank)
            if r != 2:
                continue
            assert space.rank == 5
            kernel = np.array([1, 1, 0, 0, 0, 0])
            assert all(not np.any(np.array(matrix_at(tuples, width, explicit, x, r))
                                  @ kernel % r) for x in raw)
            full = RowSpace(r, n)
            for point in product(range(r), repeat=s):
                observed = np.array(matrix_at(tuples, width, explicit, point, r))
                direct = np.zeros_like(observed)
                for e, row in enumerate(explicit):
                    direct += binary_monomial(point, e, r)*np.outer(row, row)
                for g, fiber in enumerate(tuples):
                    for x, left in fiber:
                        for y, right in tuples[(g+1) % r]:
                            delta = np.subtract(right, left)
                            weight = (binary_monomial(point, c, r)
                                      * binary_monomial(point, (x << width)+y, r))
                            direct += weight*np.outer(delta, delta)
                assert np.array_equal(direct % r, observed)
                for row in observed.tolist():
                    full.add(row)
            assert full.rank == n
        print(f"E/F5 (a,b)=({a},{b}), r={r}: {len(raw)} distinct unguarded points; "
              f"ranks across all {r} targets {ranks}/{n}.")
        if r == 2:
            print("Both complete 256-point F2 grids: rank 6; direct pairs agree with "
                  "neighboring fibers at every point.")


def displacement_identity():
    r = 19
    cycle = curve_cycle(13, 4, 1)
    assert len(cycle) == r
    blocks = []
    for alpha, beta in ((1, 0), (2, 7)):
        block = [(p, (alpha*p[0]+beta) % 13+3*(p[1] % 2)) for p in cycle[1:]
                 if (alpha*p[0]+beta) % 13 < 3]
        blocks.append(block)
    tuples, width = full_tuples(blocks, cycle)
    side = moments(tuples, (1,)*width, r)
    matrices = [displacement(side, side, g, r) for g in range(r)]
    total = sum(matrices) % r
    exact = np.zeros_like(total)
    count = prod(len(block)**2 for block in blocks)
    start = 0
    for block in blocks:
        n = len(block)
        end = start+n
        exact[start:end, start:end] = (4*(count//n)**2 % r
                                       * (n*np.eye(n, dtype=np.int64)-1))
        start = end
    assert np.array_equal(exact % r, total)
    rank = rank_and_det(total.tolist(), r)[0]
    assert rank == sum(len(b)-1 for b in blocks)
    tested = 0
    for j in range(start):
        i = next(i for i in range(start) if total[i, j] % r)
        seq = [int(g[i, j]) for g in matrices]
        rotations = [seq[h:]+seq[:h] for h in range(r)]
        assert rank_and_det(rotations, r)[0] == r
        tested += 1
    print(f"Order-19 affine-mask control: aggregate shift rank {rank}; "
          f"{tested} nonconstant coordinate directions have cyclic translate rank 19.")


def order_at_one(coefficients, r):
    # Synthetic division by U-1, with ascending coefficients.
    coeff = [int(v) % r for v in coefficients]
    while coeff and not coeff[-1]:
        coeff.pop()
    if not coeff:
        return r
    order = 0
    while sum(coeff) % r == 0:
        quotient = [0]*(len(coeff)-1)
        carry = 0
        for j in range(len(coeff)-1, 0, -1):
            carry = (carry+coeff[j]) % r
            quotient[j-1] = carry
        order += 1
        coeff = quotient
    return order


def weighted_coefficients():
    p, r, threshold = 251, 257, 23
    cycle = curve_cycle(p, 20, 1)
    assert len(cycle) == r and pow(3, 128, r) == r-1
    blocks = []
    for beta in (0, 23):
        block = [(j, (x+beta) % p+threshold*(y % 2))
                 for j, (x, y) in enumerate(cycle[1:], 1) if (x+beta) % p < threshold]
        blocks.append(block)
    density = len({(a+b) % r for a, _ in blocks[0] for b, _ in blocks[1]})
    assert 8*density > 7*r and 4*sum(len(b)-1 for b in blocks) < r
    slots = (blocks[0], blocks[0], blocks[1], blocks[1])*2
    # p has eight bits. Each six-bit slot address occupies three base-four digits.
    def evaluate(t, h, pins):
        sides = []
        local_orders = []
        for side in range(2):
            distribution = np.zeros(r, dtype=np.int64)
            distribution[0] = 1
            for k in range(4*side, 4*side+4):
                shift = 3*(7-k)
                points = slots[k] if k not in pins else [slots[k][pins[k]]]
                factor = np.zeros(r, dtype=np.int64)
                updated = np.zeros_like(distribution)
                for log, address in points:
                    exponent = sum(((address >> (2*j)) & 3)*pow(t, shift+j, r-1)
                                   for j in range(3))
                    weight = pow(3, (h*exponent) % (r-1), r)
                    factor[log if side else -log] = weight
                    updated += weight*np.roll(distribution, log)
                local_orders.append(order_at_one(factor, r))
                distribution = updated % r
            sides.append(distribution)
        left, right = sides
        counts = np.array([int(left @ np.roll(right, -j)) % r for j in range(r)])
        actual = order_at_one(counts, r)
        assert actual == sum(local_orders) < r
        return counts, actual

    zero_at_generator = positive_order = None
    for t in range(1, 12):
        for h in range(55):
            counts, order = evaluate(t, h, {})
            if not counts[1] and zero_at_generator is None:
                zero_at_generator = (t, h, int(np.count_nonzero(counts)))
            if order and positive_order is None:
                positive_order = (t, h, order)
            if zero_at_generator and positive_order:
                break
        if zero_at_generator and positive_order:
            break
    assert zero_at_generator and positive_order
    for pins in ({1: 0}, {0: 1, 7: 2}):
        evaluate(positive_order[0], positive_order[1], pins)
    print(f"Order-257 affine-mask control: {density} half sums; zero coefficient at "
          f"the generator for (t,h,nonzero displacements)={zero_at_generator}; "
          f"positive order at one (t,h,order)={positive_order}.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without Python optimization.")
    finite_grids()
    displacement_identity()
    weighted_coefficients()
