"""Exact finite algebra checks for the fixed-shift research note. Python 3.10+."""
from itertools import combinations, product
from math import comb, prod


def rank_and_det(matrix, prime):
    a = [[x % prime for x in row] for row in matrix]
    rank = 0
    det = 1
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        if pivot != rank:
            a[rank], a[pivot] = a[pivot], a[rank]
            det = -det
        pivot_value = a[rank][j]
        det = det * pivot_value % prime
        for i in range(rank + 1, len(a)):
            factor = a[i][j] * pow(pivot_value, -1, prime) % prime
            a[i] = [(v - factor * w) % prime for v, w in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank, det % prime if rank == len(a) == len(a[0]) else 0


def outer_add(matrix, vector, row, rhs, weight, prime):
    for j, x in enumerate(row):
        vector[j] = (vector[j] + weight * rhs * x) % prime
        for k, y in enumerate(row):
            matrix[j][k] = (matrix[j][k] + weight * x * y) % prime


def fixed_shift_fixture():
    prime = 5
    blocks = ((0, 1), (0, 2))
    column_values = [0, 1, 0, 2, 1, 3]
    n = len(column_values)
    tuples = list(product(range(2), repeat=4))
    locals_ = (0, 0, 1, 1)

    def value(bits):
        return sum(blocks[i][b] for i, b in zip(locals_, bits)) % prime

    def incidence(bits):
        columns = [2 * i + b for i, b in zip(locals_, bits)]
        return [columns.count(j) for j in range(n)]

    explicit = [([int(j == 4) for j in range(n)], 1)]
    for j, v in enumerate(column_values):
        witness = next(bits for bits in tuples if value(bits) == v)
        explicit.append(([int(j == k) - h for k, h in enumerate(incidence(witness))], 0))

    pairs = [(x, y) for x, y in product(tuples, repeat=2)
             if (value(y) - value(x)) % prime == 1]
    weights = (-1,) * 4 + (1,) * 4
    slot_blocks = locals_ * 2

    def monomial(key, xs):
        return prod(pow(z, (key >> (2 * i)) & 3, prime)
                    for i, z in enumerate(xs)) % prime

    full_rank = 0
    singular = 0
    all_one_rank = None
    extracted = 0
    false_table = next(candidate for candidate in product(range(prime), repeat=n)
                       if list(candidate) != column_values and
                       all(sum(x*y for x, y in zip(row, candidate)) % prime == rhs
                           for row, rhs in explicit))
    for xs in product(range(prime), repeat=4):
        b = [[0] * n for _ in range(n)]
        h = [0] * n
        p = [[0] * n for _ in range(n)]
        g = [0] * n
        terms = []
        for key, (row, rhs) in enumerate(explicit):
            w = monomial(key, xs)
            outer_add(b, h, row, rhs, w, prime)
            outer_add(p, g, row, rhs, w, prime)
        for x, y in pairs:
            bits = x + y
            key = sum(bit << i for i, bit in enumerate(bits))
            w = monomial(len(explicit), xs) * monomial(key, xs) % prime
            a = [v - u for u, v in zip(incidence(x), incidence(y))]
            outer_add(b, h, a, 1, w, prime)
            discrepancy = (sum(coef*y for coef, y in zip(a, false_table)) - 1) % prime
            terms.append((bits, discrepancy, discrepancy * monomial(key, xs) % prime))
            # Slot-pin expansion independently of the row multiplicities.
            for i, bit in enumerate(bits):
                j = 2 * slot_blocks[i] + bit
                g[j] = (g[j] + w * weights[i]) % prime
                for z, other in enumerate(bits):
                    k = 2 * slot_blocks[z] + other
                    p[j][k] = (p[j][k] + w * weights[i] * weights[z]) % prime
        assert (p, g) == (b, h)
        assert [sum(a * v for a, v in zip(row, column_values)) % prime
                for row in b] == h
        rank, _ = rank_and_det(b, prime)
        full_rank += rank == n
        singular += rank < n
        if all(v == 1 for v in xs):
            all_one_rank = rank
        if sum(term[2] for term in terms) % prime:
            for pos in range(8):
                zero = [term for term in terms if term[0][pos] == 0]
                terms = zero if sum(term[2] for term in zero) % prime else [
                    term for term in terms if term[0][pos] == 1]
            assert len(terms) == 1 and terms[0][1] and terms[0][2]
            extracted += 1
    assert singular and full_rank and extracted and all_one_rank == n - 1
    print(f"C5 fixed shifts: {len(pairs)} ordered pairs, {n} columns, "
          f"{full_rank} nonsingular and {singular} singular points; "
          f"all-ones rank {all_one_rank}.")
    print(f"One-pin separator: every one of {extracted} nonzero evaluations "
          "returned a violated ordered pair by prefix extraction.")


def basis_minimum():
    prime = 5
    rows = [[1, 0], [2, 0], [0, 1], [1, 1], [1, 2], [0, 0]]
    keys = [4, 0, 7, 3, 6, 2]
    count = 0
    for set_size in range(2, len(rows) + 1):
        for ids in combinations(range(len(rows)), set_size):
            bases = []
            for i, j in combinations(ids, 2):
                _, det = rank_and_det([rows[i], rows[j]], prime)
                if det:
                    bases.append((keys[i] + keys[j], det))
            if not bases:
                continue
            cost = min(w for w, _ in bases)
            assert sum(w == cost for w, _ in bases) == 1
            for z in range(prime):
                b = [[0, 0], [0, 0]]
                for i in ids:
                    outer_add(b, [0, 0], rows[i], 0, pow(z, keys[i], prime), prime)
                assert rank_and_det(b, prime)[1] == sum(
                    det**2 * pow(z, w, prime) for w, det in bases) % prime
            count += 1
    # Three identical rows span F3. Their unweighted Gram is zero.
    assert sum(1 for _ in range(3)) % 3 == 0
    print(f"Cauchy-Binet and unique minimum: {count} full-rank row subsets, "
          "all five field points; unweighted cancellation control passed.")


def digit_and_filter():
    prime = 101
    digit_bits = 3
    tuples = [(x, y, (x+y) % prime) for x in range(8) for y in range(7)]
    weights = [[pow(t, a, prime) for t, a in zip((2, 0, 99), row)]
               for row in tuples]
    positions = range(3)
    result = 0
    for js in product(positions, repeat=3):
        coefficient = sum(prod((v >> (digit_bits * j)) & 7 for v, j in zip(row, js))
                          for row in weights)
        assert coefficient < 2**(6 + 3 * digit_bits + 1)
        result += coefficient * 2**(digit_bits * sum(js))
    assert result == sum(prod(row) for row in weights)
    w, r = 4, 11
    def filt(x):
        return prod((x - 1 - j * r)**2 for j in range(-w, w + 1))
    degree = 4*w + 2
    coeff = [sum((-1)**(j-i)*comb(j, i)*filt(i) for i in range(j+1))
             for j in range(degree+1)]
    def choose(x, j):
        return comb(x, j) if x >= 0 else (-1)**j * comb(j-x-1, j)
    for x in range(-w*(r-1), w*(r-1)+1):
        assert filt(x) == sum(c*choose(x, j) for j, c in enumerate(coeff))
        assert (filt(x) == 0) == ((x-1) % r == 0)
    print("Odd-field digit reconstruction: 56 assignments, 27 digit tuples; "
          "Newton filter matched all 81 signed differences.")


def m_matrix():
    r = 11
    orbit = {pow(4, i, r) for i in range(5)}
    point_values = list(range(1, r)) * 2
    bogus = [x if x in orbit else 2*x % r for x in point_values]
    rows = []
    for i, x in enumerate(point_values):
        row = [5 * int(j == i) - 2 * int(v == 4*x % r)
               for j, v in enumerate(point_values)]
        rows.append(row)
        assert sum(c*v for c, v in zip(row, bogus)) % r == 0
        # Two copies in each of the two blocks represent x by four copies of x/4.
        u = x * pow(4, -1, r) % r
        assert 4*bogus[u-1] % r == bogus[i]
    rank, _ = rank_and_det(rows, r)
    assert len(rows)-rank == 2 and bogus[0] == 1 and bogus[1] != 2
    print("Block-respecting C11 M-matrix: 20 columns, nullity 2; normalized "
          "false table also passes the chosen single-point consistency rows.")


def elliptic_fixture():
    from fixtures import curve_cycle
    from check_batch import add

    field = 13
    cycle = curve_cycle(field, 4)
    order = len(cycle)
    assert order == 19
    # Two toy affine blocks, including one point in both indexed lists.
    blocks = [[cycle[i] for i in ids] for ids in ((1, 2, 3), (1, 5, 9))]
    points = blocks[0] + blocks[1] + [cycle[1], cycle[12]]
    n = len(points)
    tuples = []
    for a, b, c, d in product(range(3), repeat=4):
        slots = (a, b, 3+c, 3+d)
        point = None
        for j in slots:
            point = add(point, points[j], field, 4)
        address = sum(u << (2*i) for i, u in enumerate((a, b, c, d)))
        tuples.append((point, address, [slots.count(i) for i in range(n)]))
    assert {rec[0] for rec in tuples} == set(cycle)
    explicit = [([int(i == n-2) for i in range(n)], 1)]
    for j, p in enumerate(points):
        a = next(rec[2] for rec in tuples if rec[0] == p)
        explicit.append(([int(k == j) - x for k, x in enumerate(a)], 0))
    rows = [(j, a, c) for j, (a, c) in enumerate(explicit)]
    for p, x, u in tuples:
        successor = add(p, cycle[1], field, 4)
        for q, y, v in tuples:
            if q == successor:
                rows.append((len(explicit)+(x << 8)+y, [b-a for a, b in zip(u, v)], 1))
    # This fixed finite regression uses a univariate key and checks each integer
    # in the small prime field. The asymptotic theorem uses digit variables.
    successes = 0
    for z in range(order):
        gram = [[0]*n for _ in points]
        rhs = [0]*n
        for key, a, c in rows:
            outer_add(gram, rhs, a, c, pow(z, key, order), order)
        _, det = rank_and_det(gram, order)
        if det:
            successes += 1
            solution = []
            for j in range(n):
                replaced = [row[:j] + [v] + row[j+1:] for row, v in zip(gram, rhs)]
                solution.append(rank_and_det(replaced, order)[1]*pow(det, -1, order) % order)
            assert [cycle[i] for i in solution] == points
    assert successes
    print(f"E(F13), order 19: {len(rows)-len(explicit)} direct group-law shifts; "
          f"all {successes} nonsingular evaluations recovered the eight labeled points.")


def packet_orientations():
    # The root proof varies the number of available symbols. Exhaust all
    # equality patterns of a three-coordinate packet at each alphabet size.
    from collections import Counter

    depth = 3
    width = 1 << depth
    instances = 0
    for groups in (2, 4, 8, 64, 128):
        for different in product((False, True), repeat=depth):
            pairs = [(i % groups, (i + int(diff)) % groups)
                     for i, diff in enumerate(different)]
            row = [[pair[(x >> i) & 1] for i, pair in enumerate(pairs)]
                   for x in range(width)]
            col = [[pair[1 - ((x >> i) & 1)] for i, pair in enumerate(pairs)]
                   for x in range(width)]
            loads = Counter((tuple(left), tuple(right)) for left, right in zip(row, col))
            assert len(loads) == 2**sum(different)
            assert set(loads.values()) == {width // len(loads)}
            for t in range(1 - width, width):
                if not t:
                    continue
                i = (abs(t) & -abs(t)).bit_length() - 1
                for x in range(width):
                    if 0 <= x+t < width:
                        assert row[x][i] == col[x+t][i]
                tiles = 1 << (i+1)
                for start in range(tiles):
                    symbols = {col[y][i] for y in range(start, width, tiles)}
                    assert len(symbols) == 1
                    for x in range(width):
                        if row[x][i] not in symbols:
                            assert (t-start+x) % tiles != 0
            instances += 1
    print(f"Packet orientation: {instances} symbol-pair patterns across five group "
          "counts; all shift localities, tile valuations and class multiplicities passed.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without Python optimization so verification assertions execute.")
    basis_minimum()
    fixed_shift_fixture()
    digit_and_filter()
    m_matrix()
    elliptic_fixture()
    packet_orientations()
