"""Exact finite-field checks for Lever A: O(1)-slot pinning barriers, Forest Theorem, and decoupled fiber evaluations.

Python 3.10+ standard library only. Run without -O so verification assertions execute:
    python3 -B research/deterministic_endgame_barriers/check_pinning_barriers.py
"""
from itertools import combinations, combinations_with_replacement, product
from math import comb
import random


def ec_add(left, right, p, a):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u and (y + v) % p == 0:
        return None
    if x == u:
        slope = (3 * x * x + a) * pow(2 * y, -1, p) % p
    else:
        slope = (v - y) * pow(u - x, -1, p) % p
    w = (slope * slope - x - u) % p
    return w, (slope * (x - w) - y) % p


def curve_cycle(p, a, b=1):
    assert (4 * a**3 + 27 * b**2) % p != 0
    roots = {}
    for y in range(p):
        roots.setdefault(y * y % p, []).append(y)
    affine = [(x, y) for x in range(p) for y in roots.get((x**3 + a * x + b) % p, ())]
    gen = affine[0]
    cycle = [None]
    pt = gen
    while pt is not None:
        assert pt not in cycle
        cycle.append(pt)
        pt = ec_add(pt, gen, p, a)
    assert set(cycle[1:]) == set(affine)
    return cycle


class RowSpace:
    def __init__(self, prime, columns):
        self.prime = prime
        self.columns = columns
        self.pivots = {}

    def add(self, row):
        p = self.prime
        v = [x % p for x in row]
        for j in range(self.columns):
            if not v[j]:
                continue
            pivot = self.pivots.get(j)
            if pivot is None:
                inv = pow(v[j], -1, p)
                self.pivots[j] = [x * inv % p for x in v]
                return True
            m = v[j]
            v = [(x - m * y) % p for x, y in zip(v, pivot)]
        return False

    def add_matrix(self, mat):
        for row in mat:
            self.add(row)

    @property
    def rank(self):
        return len(self.pivots)


def cyclic_conv(a, b, r):
    out = [0] * r
    nz_a = [(i, v) for i, v in enumerate(a) if v]
    nz_b = [(j, v) for j, v in enumerate(b) if v]
    for i, va in nz_a:
        for j, vb in nz_b:
            idx = i + j
            if idx >= r:
                idx -= r
            out[idx] = (out[idx] + va * vb) % r
    return out


def slot_factored_moments(blocks, cpb, slot_weights, r):
    """Exact fiber moments (mu, nu, H) in F_r[U]/(U^r - 1) for cpb*m slot-factored weights."""
    m = len(blocks)
    num_slots = cpb * m
    offsets = []
    pos = 0
    for blk in blocks:
        offsets.append(pos)
        pos += len(blk)
    n0 = pos

    P = []
    for l in range(num_slots):
        blk = blocks[l // cpb]
        w = slot_weights[l]
        poly = [0] * r
        for a_log, wt in zip(blk, w):
            if wt % r:
                poly[a_log % r] = (poly[a_log % r] + wt) % r
        P.append(poly)

    unit = [0] * r
    unit[0] = 1
    mu = unit
    for l in range(num_slots):
        mu = cyclic_conv(mu, P[l], r)

    pref = [unit]
    for l in range(num_slots):
        pref.append(cyclic_conv(pref[-1], P[l], r))
    suff = [unit] * (num_slots + 1)
    for l in range(num_slots - 1, -1, -1):
        suff[l] = cyclic_conv(P[l], suff[l + 1], r)
    L1 = [cyclic_conv(pref[l], suff[l + 1], r) for l in range(num_slots)]

    L2 = {}
    for l1 in range(num_slots):
        mid = unit
        for l2 in range(l1 + 1, num_slots):
            L2[(l1, l2)] = cyclic_conv(cyclic_conv(pref[l1], mid, r), suff[l2 + 1], r)
            mid = cyclic_conv(mid, P[l2], r)

    nu = [[0] * n0 for _ in range(r)]
    H = [[[0] * n0 for _ in range(n0)] for _ in range(r)]

    for i, blk in enumerate(blocks):
        off = offsets[i]
        slots_i = range(cpb * i, cpb * (i + 1))
        for idx, a_log in enumerate(blk):
            c = off + idx
            a = a_log % r
            for g in range(r):
                src = (g - a) % r
                val = 0
                for sl in slots_i:
                    wt = slot_weights[sl][idx] % r
                    if wt:
                        val = (val + wt * L1[sl][src]) % r
                nu[g][c] = val
                H[g][c][c] = val

    for sl1 in range(num_slots):
        i1 = sl1 // cpb
        blk1 = blocks[i1]
        off1 = offsets[i1]
        w1 = slot_weights[sl1]
        for sl2 in range(sl1 + 1, num_slots):
            i2 = sl2 // cpb
            blk2 = blocks[i2]
            off2 = offsets[i2]
            w2 = slot_weights[sl2]
            p12 = L2[(sl1, sl2)]
            for idx1, a1 in enumerate(blk1):
                wt1 = w1[idx1] % r
                if not wt1:
                    continue
                c1 = off1 + idx1
                for idx2, a2 in enumerate(blk2):
                    wt2 = w2[idx2] % r
                    if not wt2:
                        continue
                    c2 = off2 + idx2
                    coeff = (wt1 * wt2) % r
                    sh = (a1 + a2) % r
                    for g in range(r):
                        v = (coeff * p12[(g - sh) % r]) % r
                        if c1 == c2:
                            H[g][c1][c1] = (H[g][c1][c1] + 2 * v) % r
                        else:
                            H[g][c1][c2] = (H[g][c1][c2] + v) % r
                            H[g][c2][c1] = (H[g][c2][c1] + v) % r
    return mu, nu, H


def shift_gram_from_moments(mom_x, mom_y, r, step=1):
    mux, nux, hx = mom_x
    muy, nuy, hy = mom_y
    n0 = len(nux[0])
    mat = [[0] * n0 for _ in range(n0)]
    for g in range(r):
        gp = (g + step) % r
        mx, my = mux[g], muy[gp]
        vx, vy = nux[g], nuy[gp]
        hxg, hyg = hx[g], hy[gp]
        for j in range(n0):
            vxj, vyj = vx[j], vy[j]
            row = mat[j]
            hxgj, hygj = hxg[j], hyg[j]
            for k in range(n0):
                row[k] = (row[k] + mx * hygj[k] + my * hxgj[k] - vxj * vy[k] - vyj * vx[k]) % r
    return mat


def half_sumset_size(blocks, r):
    s = {0}
    for blk in blocks:
        s = {(x + y) % r for x in s for y in blk}
    return len(s)


def true_shift_rank(blocks, cpb, r):
    m = len(blocks)
    n0 = sum(len(b) for b in blocks)
    tgt = n0 - m
    offsets = []
    pos = 0
    for b in blocks:
        offsets.append(pos)
        pos += len(b)
    first = [None] * r
    ranges = [range(len(blocks[l // cpb])) for l in range(cpb * m)]
    for idx in product(*ranges):
        val = 0
        cnt = [0] * n0
        for l, sel in enumerate(idx):
            bi = l // cpb
            val = (val + blocks[bi][sel]) % r
            cnt[offsets[bi] + sel] += 1
        if first[val] is None:
            first[val] = cnt
    sp = RowSpace(r, n0)
    for g in range(r):
        if first[g] is not None and first[(g + 1) % r] is not None:
            sp.add([b - a for a, b in zip(first[g], first[(g + 1) % r])])
    if sp.rank == tgt:
        return tgt
    for idx in product(*ranges):
        val = 0
        cnt = [0] * n0
        for l, sel in enumerate(idx):
            bi = l // cpb
            val = (val + blocks[bi][sel]) % r
            cnt[offsets[bi] + sel] += 1
        if first[(val + 1) % r] is not None or first[(val - 1) % r] is not None:
            if cnt != first[val]:
                sp.add([b - a for a, b in zip(first[val], cnt)])
                if sp.rank == tgt:
                    return tgt
    return sp.rank


def order19_blocks():
    cycle = curve_cycle(13, 4, 1)
    assert len(cycle) == 19
    return [
        [i for i, (x, y) in enumerate(cycle[1:], 1) if (a * x + b) % 13 < 3]
        for a, b in ((1, 0), (4, 1))
    ]


def weighted_fibers(blocks, weights, r):
    """Enumerate doubled tuples in two labeled blocks, grouped by their sum."""
    left, right = blocks
    size = len(left) + len(right)
    fibers = [[] for _ in range(r)]
    for i, j, k, l in product(range(len(left)), range(len(left)), range(len(right)), range(len(right))):
        value = (left[i] + left[j] + right[k] + right[l]) % r
        weight = weights[0][i] * weights[1][j] * weights[2][k] * weights[3][l] % r
        counts = [0] * size
        for column in (i, j, len(left) + k, len(left) + l):
            counts[column] += 1
        fibers[value].append((counts, weight))
    return fibers


def dot(left, right):
    return sum(x * y for x, y in zip(left, right, strict=True))


def projected_fiber(fiber, u):
    """Return sums of w, wc, ws, wcs and ws² for the tuple score s = c·u."""
    mass = score_sum = square_sum = 0
    counts_sum = [0] * len(u)
    scored_counts = [0] * len(u)
    for counts, weight in fiber:
        score = dot(counts, u)
        mass += weight
        score_sum += weight * score
        square_sum += weight * score**2
        for j, count in enumerate(counts):
            counts_sum[j] += weight * count
            scored_counts[j] += weight * count * score
    return mass, counts_sum, score_sum, scored_counts, square_sum


def verify_decoupled_evaluator():
    rng = random.Random(42)
    for r, blocks in ((5, [[0, 1], [0, 2]]), (19, order19_blocks())):
        wx = [[rng.randrange(r) for _ in block] for block in blocks for _ in range(2)]
        wy = [[rng.randrange(r) for _ in block] for block in blocks for _ in range(2)]
        matrix = shift_gram_from_moments(
            slot_factored_moments(blocks, 2, wx, r), slot_factored_moments(blocks, 2, wy, r), r
        )
        xs, ys = weighted_fibers(blocks, wx, r), weighted_fibers(blocks, wy, r)
        size = sum(map(len, blocks))
        direct = [[0] * size for _ in range(size)]
        for g in range(r):
            for (cx, weight_x), (cy, weight_y) in product(xs[g], ys[(g + 1) % r]):
                row = [y - x for x, y in zip(cx, cy)]
                weight = weight_x * weight_y
                for j, v in enumerate(row):
                    for k, w in enumerate(row):
                        direct[j][k] += weight * v * w
        assert matrix == [[v % r for v in row] for row in direct]
        if r == 5:
            print("Decoupled fiber-moment evaluator: exact match against brute-force pairs on C_5 (decoupled random weights).")
            continue

        # Contract the enumerated fibers along u before multiplying the two sides.
        u = [rng.randrange(r) for _ in range(size)]
        action = [0] * size
        quadratic = 0
        for g in range(r):
            mx, nx, sx, gx, tx = projected_fiber(xs[g], u)
            my, ny, sy, gy, ty = projected_fiber(ys[(g + 1) % r], u)
            for j in range(size):
                action[j] += mx * gy[j] + my * gx[j] - sx * ny[j] - sy * nx[j]
            quadratic += mx * ty + my * tx - 2 * sx * sy
        matrix_action = [dot(row, u) % r for row in matrix]
        assert [v % r for v in action] == matrix_action
        assert quadratic % r == dot(u, matrix_action) % r
        print("Order-19 curve over F_13: weighted Gram matrix matches direct pairs; "
              "Bu and u^TBu match independent fiber sums.")


def verify_weighted_two_pin_recovery():
    r = 7
    blocks = [[1, 2, 3], [1, 2, 4], [1, 3, 5]]
    weights = [[3, 2, 3], [1, 5, 4], [5, 3, 5]]
    size = sum(map(len, blocks))
    witness = blocks[0].index(2)  # coordinate of e_(0,2)
    free = slot_factored_moments(blocks, 1, weights, r)
    pinned = []
    for i, row in enumerate(weights):
        for j in range(len(row)):
            restricted = list(weights)
            restricted[i] = [w if k == j else 0 for k, w in enumerate(row)]
            pinned.append(slot_factored_moments(blocks, 1, restricted, r))

    source, target, both = (RowSpace(r, size) for _ in range(3))
    for moment in pinned:
        for space, left, right in ((source, moment, free), (target, free, moment)):
            matrix = shift_gram_from_moments(left, right, r)
            assert all(row[witness] == 0 for row in matrix)
            space.add_matrix(matrix)
    detects_witness = False
    for left, right in product(pinned, repeat=2):
        matrix = shift_gram_from_moments(left, right, r)
        detects_witness |= any(row[witness] != 0 for row in matrix)
        both.add_matrix(matrix)
    assert source.rank == target.rank == 5
    assert both.rank == size - len(blocks) == 6 and detects_witness
    print("Weighted F_7 fixture: source-only and target-only pins each have rank 5/6 "
          "and kill e_(0,2); two pins detect it and have rank 6/6.")


def verify_theorem_a1_lucas_collapse():
    """Verify Theorem A1: m = r - 1, F_i = {1, 2}, doubled tuples => K_{A1..Ap} == 0 for p <= r - 7, deficient E."""
    for r in (11, 13):
        m = r - 1
        blks = [[1, 2]] * m
        assert half_sumset_size(blks, r) == r
        n0 = 2 * m
        ones = [[1, 1] for _ in range(2 * m)]
        mom_1 = slot_factored_moments(blks, 2, ones, r)
        b11 = shift_gram_from_moments(mom_1, mom_1, r)
        assert all(v == 0 for row in b11 for v in row)

        moms_pin = {}
        for i in range(m):
            for idx in range(2):
                w = [list(row) for row in ones]
                w[2 * i] = [1 if t == idx else 0 for t in range(2)]
                mx = slot_factored_moments(blks, 2, w, r)
                moms_pin[(i, idx)] = mx
                mat1 = shift_gram_from_moments(mx, mom_1, r)
                assert all(v == 0 for row in mat1 for v in row)

        for mx in moms_pin.values():
            for my in moms_pin.values():
                mat2 = shift_gram_from_moments(mx, my, r)
                assert all(v == 0 for row in mat2 for v in row)

        # Spot-check 3-pin and 4-pin orbits (all must vanish since p <= 4 <= r - 7)
        for p_x in [(0, 1), (0, 2)]:
            for vx in [(0, 0), (0, 1), (1, 1)]:
                wx = [list(row) for row in ones]
                wx[p_x[0]] = [1 if t == vx[0] else 0 for t in range(2)]
                wx[p_x[1]] = [1 if t == vx[1] else 0 for t in range(2)]
                mx = slot_factored_moments(blks, 2, wx, r)
                for p_y in [(0, 1), (2, 4)]:
                    for vy in [(0, 0), (0, 1), (1, 1)]:
                        wy = [list(row) for row in ones]
                        wy[p_y[0]] = [1 if t == vy[0] else 0 for t in range(2)]
                        wy[p_y[1]] = [1 if t == vy[1] else 0 for t in range(2)]
                        my = slot_factored_moments(blks, 2, wy, r)
                        mat4 = shift_gram_from_moments(mx, my, r)
                        assert all(v == 0 for row in mat4 for v in row)

        # Construct valid deficient explicit matrix E killing u_{(0,2)} = 1, u_{(1,2)} = -1
        n_full = n0 + 2
        labels = [1, 2] * m + [1, 2]
        u_ker = [0] * n_full
        u_ker[1] = 1
        u_ker[3] = r - 1
        zero_slots = [2, 2, 3, 3]
        explicit = [[int(j == n_full - 2) for j in range(n_full)]]
        for col_j, g in enumerate(labels):
            k = (g - 2 * m) % r
            twos = ([0] + zero_slots[:k - 1]) if col_j == 1 else (
                ([1] + zero_slots[:k - 1]) if col_j == 3 else zero_slots[:k]
            )
            cnt = [0] * n_full
            for bi in range(m):
                c2 = twos.count(bi)
                cnt[2 * bi] = 2 - c2
                cnt[2 * bi + 1] = c2
            assert sum(cnt[j] * labels[j] for j in range(n_full)) % r == g
            row = [int(j == col_j) - cnt[j] for j in range(n_full)]
            assert sum(a * b for a, b in zip(row, u_ker)) % r == 0
            explicit.append(row)
        sp_E = RowSpace(r, n_full)
        sp_E.add_matrix(explicit)
        assert sp_E.rank == n_full - 2
        print(f"Theorem A1 (Lucas collapse, r={r}, m={m}, F={{1,2}}): all1=K1=K2=K4=0 mod {r}; "
              f"rank(E)=rank(E+K2)={sp_E.rank}/{n_full}.")


def verify_theorem_a2_small_m_binomial():
    """Verify Theorem A2: strict 4m < r, |D_0| = r, positive binomial weights => K_{A1..Ap} == 0 for p <= m - 3."""
    for r, m, d in [(19, 3, 6), (29, 4, 7), (31, 5, 6), (37, 6, 6)]:
        assert m * d == r - 1 and 4 * m < r and d >= 5
        F = list(range(1, d + 2))
        L = d + 1
        w = [comb(d, j) % r for j in range(L)]
        assert all(v != 0 for v in w) and sum(w) % r != 0
        blks = [F] * m
        assert half_sumset_size(blks, r) == r
        n0 = m * L
        base_w = [list(w) for _ in range(2 * m)]
        mom_1 = slot_factored_moments(blks, 2, base_w, r)
        assert all(v == 0 for row in shift_gram_from_moments(mom_1, mom_1, r) for v in row)

        # After p external pins, a Gram entry fixes at most p+2 slots.
        # Sum the integer binomial coefficients in every possible residue class.
        for fixed in range(m):
            free_degree = (4 * m - fixed) * d
            for residue in range(r - 1 - fixed * d, r):
                count = sum(comb(free_degree, k) for k in range(residue, free_degree + 1, r))
                assert count % r == 0

        if m >= 4:  # One external pin is in range.
            pinned = []
            for i in range(m):
                for idx in range(L):
                    wp = [list(row) for row in base_w]
                    wp[2 * i] = [w[idx] if t == idx else 0 for t in range(L)]
                    mx = slot_factored_moments(blks, 2, wp, r)
                    pinned.append(mx)
                    assert all(v == 0 for row in shift_gram_from_moments(mx, mom_1, r) for v in row)
            if m >= 5:  # One external pin on each side is in range.
                for mx, my in product(pinned, repeat=2):
                    assert all(v == 0 for row in shift_gram_from_moments(mx, my, r) for v in row)

        n_full = n0 + 2
        labels = F * m + [1, 2]
        u_ker = [int(g == d + 1) for g in labels]
        normalization = [int(j == n_full - 2) for j in range(n_full)]
        assert any(u_ker) and u_ker[n_full - 2] == 0
        explicit = [normalization]
        tuples = {}
        for g in F:
            entries = [d + 1] if g == d + 1 else []
            reserved = len(entries)
            remaining = (g - 2 * m - reserved * d) % r
            for _ in range(2 * m - reserved):
                increment = min(d - 1, remaining)
                entries.append(1 + increment)
                remaining -= increment
            assert remaining == 0 and sum(entries) % r == g
            assert len(entries) == 2 * m and all(value in F for value in entries)
            assert entries.count(d + 1) == int(g == d + 1)
            tuples[g] = entries
        for col_j, g in enumerate(labels):
            row = [int(j == col_j) for j in range(n_full)]
            for slot, value in enumerate(tuples[g]):
                row[(slot // 2) * L + value - 1] -= 1
            assert sum(a * b for a, b in zip(row, u_ker)) % r == 0
            explicit.append(row)
        sp_E = RowSpace(r, n_full)
        sp_E.add_matrix(explicit)
        assert sp_E.rank < n_full
        print(f"Theorem A2 (r={r}, m={m}, d={d}): binomial residue sums vanish through p={m - 3}; "
              f"all {len(tuples)} target decompositions verified; rank(E)={sp_E.rank}/{n_full}.")

    # Also verify cpb=1 dense small-m binomial collapse (r=13, m=4, d=3 -> K2=5/12; r=19, m=6, d=3 -> K2=1/18)
    for r, m, d, expected_k2 in [(13, 4, 3, 5), (19, 6, 3, 1)]:
        F = list(range(1, d + 2))
        w = [comb(d, j) % r for j in range(d + 1)]
        blks = [F] * m
        n0 = m * (d + 1)
        base_w = [list(w) for _ in range(m)]
        moms_pin = {}
        for i in range(m):
            for idx in range(d + 1):
                wp = [list(row) for row in base_w]
                wp[i] = [w[idx] if t == idx else 0 for t in range(d + 1)]
                moms_pin[(i, idx)] = slot_factored_moments(blks, 1, wp, r)
        sp2 = RowSpace(r, n0)
        for mx in moms_pin.values():
            for my in moms_pin.values():
                sp2.add_matrix(shift_gram_from_moments(mx, my, r))
        assert sp2.rank == expected_k2
        print(f"Theorem A2 (cpb=1 variant, r={r}, m={m}, d={d}, 2m={2*m}<{r}, |D0|={r}/{r}): "
              f"K2 rank={sp2.rank}/{n0 - m} (nullity={n0 - m - sp2.rank}).")


def verify_proposition_a3_forest_theorem():
    """Verify Proposition A3: exact equality rank(K2(b0)) == true_shift_rank for m=2, cpb=1 across ALL pairs."""
    for r, expected_pairs in [(5, 66), (7, 1653)]:
        subsets = [list(c) for sz in range(2, r) for c in combinations(range(1, r), sz)]
        checked = 0
        for i, b1 in enumerate(subsets):
            for b2 in subsets[i:]:
                blks = [b1, b2]
                tr = true_shift_rank(blks, 1, r)
                n0 = len(b1) + len(b2)
                moms_b0 = []
                for idx in range(len(b1)):
                    w = [[1 if t == idx else 0 for t in range(len(b1))], [1] * len(b2)]
                    moms_b0.append(slot_factored_moments(blks, 1, w, r))
                sp = RowSpace(r, n0)
                for mx in moms_b0:
                    for my in moms_b0:
                        sp.add_matrix(shift_gram_from_moments(mx, my, r))
                assert sp.rank == tr
                checked += 1
        assert checked == expected_pairs
        print(f"Proposition A3 (Forest Theorem, m=2, cpb=1, r={r}): all {checked} subset pairs match true shift rank.")


def verify_unweighted_4m_lt_r_controls():
    """Verify unweighted 0/1 exhaustive controls at 2*cpb*m < r, elliptic curves, and Paley QR blocks."""
    # An unweighted interval can repeat a nonzero difference even well below sqrt(r).
    interval = range(1, 9)
    prime = 65537
    unit_differences = sum((a - b) % prime == 1 for a, b in product(interval, repeat=2))
    assert unit_differences == len(interval) - 1
    print(f"Difference multiplicity: an eight-element interval modulo {prime} has "
          f"{unit_differences} ordered pairs with difference one.")

    def canon_classes(r, min_sz=2):
        seen = {}
        for sz in range(min_sz, r):
            for comb_s in combinations(range(1, r), sz):
                key = min(tuple(sorted((x - c) % r for x in comb_s)) for c in comb_s)
                if key not in seen:
                    seen[key] = list(comb_s)
        return list(seen.values())

    # 1. Exhaustive all 969 canonical triples at r=7, m=3, cpb=1 (2m=6 < 7)
    r = 7
    reps7 = canon_classes(r)
    triples = list(combinations_with_replacement(reps7, 3))
    assert len(triples) == 969
    lam_drops = 0
    for b1, b2, b3 in triples:
        blks = [b1, b2, b3]
        n0 = len(b1) + len(b2) + len(b3)
        tgt = n0 - 3
        ones = [[1] * len(b) for b in blks]
        mom_1 = slot_factored_moments(blks, 1, ones, r)
        sp1 = RowSpace(r, n0)
        sp_lam = RowSpace(r, n0)
        lam = b1 + b2 + b3
        for i, blk in enumerate(blks):
            for idx in range(len(blk)):
                w = [list(row) for row in ones]
                w[i] = [1 if t == idx else 0 for t in range(len(blk))]
                mx = slot_factored_moments(blks, 1, w, r)
                mat = shift_gram_from_moments(mx, mom_1, r)
                sp1.add_matrix(mat)
                sp_lam.add([sum(lam[j] * mat[j][k] for j in range(n0)) % r for k in range(n0)])
        if sp1.rank < tgt:
            assert sp1.rank == true_shift_rank(blks, 1, r)
        if 8 * half_sumset_size(blks, r) > 7 * r and sp_lam.rank < tgt:
            lam_drops += 1
    assert lam_drops == 66
    print(f"Unweighted exhaustive control (r=7, m=3, cpb=1, 2m=6<7): 0/969 K1/K2 drops "
          f"(vs {lam_drops} dense drops of 1-row lambda^T K1).")

    # 2. Exhaustive all 185 identical pairs (F, F) at r=11, m=2, cpb=2 (4m=8 < 11)
    r = 11
    reps11 = canon_classes(r)
    assert len(reps11) == 185
    dense11 = 0
    for F in reps11:
        blks = [F, F]
        is_dense = (8 * half_sumset_size(blks, r) > 7 * r)
        if is_dense:
            dense11 += 1
        n0 = 2 * len(F)
        tgt = n0 - 2
        ones = [[1] * len(F) for _ in range(4)]
        mom_1 = slot_factored_moments(blks, 2, ones, r)
        sp1 = RowSpace(r, n0)
        for i in range(2):
            for idx in range(len(F)):
                w = [list(row) for row in ones]
                w[2 * i] = [1 if t == idx else 0 for t in range(len(F))]
                mx = slot_factored_moments(blks, 2, w, r)
                sp1.add_matrix(shift_gram_from_moments(mx, mom_1, r))
                if sp1.rank == tgt:
                    break
            if sp1.rank == tgt:
                break
        if sp1.rank < tgt:
            assert not is_dense
            assert sp1.rank == true_shift_rank(blks, 2, r)
    assert dense11 == 130
    print(f"Unweighted exhaustive control (r=11, m=2, cpb=2, 4m=8<11): 0/185 identical (F,F) K1/K2 drops ({dense11} dense).")

    # 3. Elliptic curve controls: E(F_13) r=19 (all-ones drop 7/8 -> K1=8/8), PR#4 order-37 (12/12), E(F_251) r=257 (50/50)
    blks19 = order19_blocks()
    ones19 = [[1] * len(blks19[l // 2]) for l in range(4)]
    m19 = slot_factored_moments(blks19, 2, ones19, 19)
    sp19_all1 = RowSpace(19, 10); sp19_all1.add_matrix(shift_gram_from_moments(m19, m19, 19))
    sp19_k1 = RowSpace(19, 10)
    for idx in range(len(blks19[0])):
        w = [list(row) for row in ones19]
        w[0] = [1 if t == idx else 0 for t in range(len(blks19[0]))]
        sp19_k1.add_matrix(shift_gram_from_moments(slot_factored_moments(blks19, 2, w, 19), m19, 19))
    assert sp19_all1.rank == 7 and sp19_k1.rank == 8

    blks37 = [[1, 3, 5, 10, 18, 19, 27, 32, 34, 36], [18, 19, 21, 31]]
    m37 = slot_factored_moments(blks37, 2, [[1] * len(blks37[l // 2]) for l in range(4)], 37)
    sp37 = RowSpace(37, 14); sp37.add_matrix(shift_gram_from_moments(m37, m37, 37))
    assert sp37.rank == 12

    cyc251 = curve_cycle(251, 20, 1)
    blks257 = [[j for j, (x, y) in enumerate(cyc251[1:], 1) if (x + beta) % 251 < 23] for beta in (0, 23)]
    m257 = slot_factored_moments(blks257, 2, [[1] * len(blks257[l // 2]) for l in range(4)], 257)
    sp257 = RowSpace(257, 52); sp257.add_matrix(shift_gram_from_moments(m257, m257, 257))
    assert sp257.rank == 50
    print("Elliptic curve controls: E/F13(r=19) all1=7/8 rescued to K1(b0)=8/8; PR#4 r=37 all1=12/12; E/F251(r=257) all1=50/50.")

    # 4. Paley QR blocks Q_r (r in {7, 11, 19}, where (Q_r + 1/2)^3 = 0 mod r)
    paley_checked = 0
    for r in (7, 11, 19):
        qr = sorted({(x * x) % r for x in range(1, r)})
        for cpb, m in [(1, 2), (1, 3), (2, 2)]:
            if 2 * cpb * m >= r:
                continue
            blks = [qr] * m
            n0 = m * len(qr); tgt = n0 - m
            ones = [[1] * len(qr) for _ in range(cpb * m)]
            m1 = slot_factored_moments(blks, cpb, ones, r)
            sp1 = RowSpace(r, n0)
            for i in range(m):
                for idx in range(len(qr)):
                    w = [list(row) for row in ones]
                    w[cpb * i] = [1 if t == idx else 0 for t in range(len(qr))]
                    sp1.add_matrix(shift_gram_from_moments(slot_factored_moments(blks, cpb, w, r), m1, r))
                    if sp1.rank == tgt:
                        break
                if sp1.rank == tgt:
                    break
            assert sp1.rank == tgt
            paley_checked += 1
    print(f"Unweighted Paley QR control (r in {{7,11,19}}): {paley_checked}/{paley_checked} instances reach full K1/K2 rank.")


def verify_fixed_regime_a1_adversary():
    """Verify fiber-respecting rank drops against fixed small Regime A-I lists at r=7 (coupled L=2 & 1-sided zero L=2)."""
    r, sz = 7, 4
    slot_bits, radix_bits, s_side = 2, 2, 4
    addrs = [list(range(sz)), list(range(sz))]

    def build_w(xv, yv):
        def side(vs):
            out = []
            for l in range(4):
                bs = l * slot_bits
                wts = []
                for a_addr in addrs[l // 2]:
                    w = 1
                    for b in range(slot_bits):
                        if (a_addr >> b) & 1:
                            pos = bs + b
                            w = (w * pow(vs[pos // radix_bits], 1 << (pos % radix_bits), r)) % r
                    wts.append(w)
                out.append(wts)
            return out
        return side(xv), side(yv)

    def stack_rank(blks, pts):
        sp = RowSpace(r, 2 * sz)
        for xv, yv in pts:
            wx, wy = build_w(xv, yv)
            sp.add_matrix(shift_gram_from_moments(
                slot_factored_moments(blks, 2, wx, r),
                slot_factored_moments(blks, 2, wy, r),
                r,
            ))
        return sp.rank

    # Coupled L=2 drops on [[1,2,3,4],[1,2,3,4]] (g=3 primitive root mod 7)
    coupled_l2 = [((1, 1, 1, 1), (1, 1, 1, 1)), ((3, 3, 3, 3), (3, 3, 3, 3))]
    rk_c = stack_rank([[1, 2, 3, 4], [1, 2, 3, 4]], coupled_l2)
    assert rk_c == 4

    # 1-sided zero-projection L=2 drops on [[1,3,5,6],[1,3,5,6]]
    zero_l2 = [((1, 1, 1, 1), (1, 1, 1, 1)), ((0, 1, 1, 1), (1, 1, 1, 1))]
    rk_z = stack_rank([[1, 3, 5, 6], [1, 3, 5, 6]], zero_l2)
    assert rk_z == 5

    # Decoupled asymmetric L=2 recovers full rank 6/6 on both adversarial pairs
    asym_l2 = [((3, 2, 6, 4), (6, 4, 2, 3)), ((2, 4, 1, 2), (2, 6, 5, 3))]
    assert stack_rank([[1, 2, 3, 4], [1, 2, 3, 4]], asym_l2) == 6
    assert stack_rank([[1, 3, 5, 6], [1, 3, 5, 6]], asym_l2) == 6
    print(f"Regime A-I fixed-list adversary (r=7, dense m=2): coupled L=2 drops to {rk_c}/6; "
          f"1-sided zero-proj L=2 drops to {rk_z}/6; decoupled asymmetric L=2 hits 6/6 on both.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without Python optimization so assertions execute.")
    verify_decoupled_evaluator()
    verify_weighted_two_pin_recovery()
    verify_theorem_a1_lucas_collapse()
    verify_theorem_a2_small_m_binomial()
    verify_proposition_a3_forest_theorem()
    verify_unweighted_4m_lt_r_controls()
    verify_fixed_regime_a1_adversary()
    print("CHECK_PINNING_BARRIERS: ALL CHECKS PASSED.")
