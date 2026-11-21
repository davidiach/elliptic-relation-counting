"""Finite checks of arithmetic joins and boundary ranks. Standard library only."""

from collections import Counter
from itertools import product
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hillclimb_001"))
import common
import arithmetic_alternating_chord as c38
import simplex_contraction as sc


def chord_pairs(p, a, u, w, y, odd, square_roots):
    sign = -1 if odd else 1
    g = (w - u) % p
    h = (u*u + u*w + w*w + a) % p
    if g:
        disc = (y*y - sign*g*h) % p
        slopes = {(root - y) * pow(g, -1, p) % p for root in square_roots[disc]}
    elif y:
        slopes = {-sign*h * pow(2*y, -1, p) % p}
    else:
        assert h    # The anchored point lies on a nonsingular curve.
        slopes = set()
    return {((s*s-u-w) % p, s) for s in slopes}


def check_chord_slice(p, a):
    points = c38.points(p, a, 1)
    square_roots = {value: [] for value in range(p)}
    for root in range(p):
        square_roots[root * root % p].append(root)
    counts = Counter()
    for odd in (False, True):
        slot = 3 if odd else 2
        for (anchor_x, anchor_y), other_x in product(points, range(p)):
            u, w = (anchor_x, other_x) if odd else (other_x, anchor_x)
            sign = -1 if odd else 1
            g = (w - u) % p
            h = (u*u + u*w + w*w + a) % p
            disc = (anchor_y**2 - sign*g*h) % p
            assert disc == (other_x**3 + a*other_x + 1) % p
            if not g and anchor_y:
                counts["linear"] += 1
            elif not g:
                counts["two_torsion_coincident"] += 1
            pairs = chord_pairs(p, a, u, w, anchor_y, odd, square_roots)
            # Read the two actual equations for this slot from the vendor module.
            constants = dict(a=a, Xm=w, Ym=0, p=p)
            for j in range(1, 4):
                constants[f"alpha{j}"], constants[f"beta{j}"] = 3, 1
            source = dict.fromkeys(c38.source_names(3), 0)
            source.update(x1=u, X2=u if odd else w, Y2=anchor_y)
            found = set()
            for v, slope in product(range(p), repeat=2):
                source[f"x{slot}"], source[f"l{slot}"] = v, slope
                outputs = c38.operations(source, constants, 3)
                start = 2 * (slot - 2)
                if outputs[start] % p == outputs[start + 1] % p == 0:
                    found.add((v, slope))
            assert found == pairs and len(found) <= 2
            counts["slices"] += 1
            counts["literal_pairs"] += p**2
            counts[f"roots_{len(pairs)}"] += 1
    return counts


def check_deformed_join():
    """Compare two analytically supported chord slices with a normalized packet join."""
    p, a = 13, 4
    roots = {d: [s for s in range(p) if s*s % p == d] for d in range(p)}
    options = []
    for odd in (False, True):
        for (x, y), other in product(c38.points(p, a, 1), range(p)):
            u, w = (x, other) if odd else (other, x)
            pairs = chord_pairs(p, a, u, w, y, odd, roots)
            legal = [(v, s) for v, s in sorted(pairs) if v != u and (3*v+1) % 16 < 8]
            if len(legal) == 2:
                words = tuple(tuple((n >> j) & 1 for n in pair for j in range(4))
                              for pair in legal)
                if min(map(sum, words)) in (1, 2):
                    options.append((odd, words))
    target = next(words for odd, words in options if not odd)
    source = next(words for odd, words in options if odd and words != target)
    target = dict(zip(target, (2, 5)))
    source = dict(zip(source, (3, 7)))
    ports = tuple(range(8))
    net = sc.network("conditioned-chord-pair", [(ports, lambda key: target.get(key, 0))],
                     [(ports, lambda key: source.get(key, 0))], sc.Alphabet(1, 1), p=8)
    labels = {"r0": (0, 0), "c0": (1, 0)}
    packets = sc.geom.packetize(net, labels, 4, 2)
    assert len(packets) == 1 and len(packets[0]["positions"]) == 8
    ring = sc.Ring(60, (("z", 5),))
    weights, _ = sc.geom.local_polynomials(4, 1, 5)
    W = {delta: ring.polynomial(terms, ring.variable("z")) for delta, terms in weights.items()}
    positions = packets[0]["positions"]
    direct = ring.zero
    for xbits, xweight in source.items():
        raw = {}
        for pos in positions:
            value = ring.zero
            for incoming, edge in positions.items():
                if xbits[edge]:
                    value = ring.add(value, W[incoming-pos])
            raw[pos] = value
        for ybits, yweight in target.items():
            value = ring.val(xweight*yweight)
            for pos, edge in positions.items():
                value = ring.mul(value, raw[pos] if ybits[edge] else ring.sub(ring.one, raw[pos]))
            direct = ring.add(direct, value)

    # Independent reference branch uses the native normalized crossing modes.
    _, modes, forms = sc.modes_and_forms(packets, 4, 1, ring, "z", 1)
    engine = sc.Engine(net.alphabet, groups=2)
    count = 0
    full = binary = ring.zero
    for beta in sc.geom.assignments(modes, tuple(range(len(modes))), ring.p, 4):
        count += 1
        term = ring.one
        for group, support in enumerate((target, source)):
            child = engine.child(net, labels, packets, 0, group, beta, modes, forms, ring)
            total = ring.zero
            packet = child.cols[0]
            for bits, weight in support.items():
                port_values = tuple(bits[e] for _, e in packet.ports)
                entry = packet.entry(port_values, ring)
                total = ring.add(total, ring.scale(entry, weight))
            term = ring.mul(term, total)
        full = ring.add(full, term)
        if max(beta) <= 1:
            binary = ring.add(binary, term)
    assert full == direct and full != binary and any(full[1:])
    print(f"Deformed mixed packet: {len(modes)} modes, {count} exponent vectors, "
          f"{len(target)*len(source)} weighted support pairs, jet={full}")


def negate(point, p):
    return None if point is None else (point[0], -point[1] % p)


def hash_slot(points, bits, alpha, beta, threshold):
    return tuple(point for point in points
                 if (alpha*point[0]+beta) % (1 << bits) < threshold)


def interval_oracle(bases, add):
    """Materialize and charge ordinary convolution tables for tiny test groups."""
    tables = {}
    additions = 0
    for left in range(len(bases)):
        row = Counter({None: 1})
        for right in range(left, len(bases)):
            updated = Counter()
            for state, value in row.items():
                for factor in bases[right]:
                    updated[add(state, factor)] += value
                    additions += 1
            row = updated
            tables[left, right] = row

    def coefficient(left, right, state):
        if left > right:
            return int(state is None)
        return tables[left, right][state]

    return coefficient, additions


def sparse_guard_counts(bases, group, add, p):
    interval, work = interval_oracle(bases, add)
    corrections = []
    calls = 0

    def at(right, state):
        nonlocal calls
        value = interval(0, right, state)
        calls += 1
        for index, row in enumerate(corrections, start=1):
            for location, delta in row.items():
                value += delta * interval(index+1, right, add(state, negate(location, p)))
                calls += 1
        return value

    supports = []
    for index in range(1, len(bases)):
        slot = bases[index]
        incoming = {point: at(index-1, point) for point in slot}
        row = Counter()
        for point in slot:
            row[None] -= incoming[negate(point, p)]
            twice = add(point, point)
            if twice is not None:
                row[twice] -= incoming[point]
        corrections.append({key: value for key, value in row.items() if value})
        supports.append(len(corrections[-1]))
    return {state: at(len(bases)-1, state) for state in group}, work, calls, supports


def literal_counts(bases, specs, add, p, a):
    """Count ordered tuples and verify all C38 arithmetic equations on surviving paths."""
    counts = Counter()
    for factors in product(*bases):
        partial = factors[0]
        m = len(factors)
        sources = {f"x{j}": point[0] for j, point in enumerate(factors, 1)}
        for j, point in enumerate(factors[1:], 2):
            if partial[0] == point[0]:
                break
            sources[f"l{j}"] = (point[1]-partial[1]) * pow(point[0]-partial[0], -1, p) % p
            partial = add(partial, point)
            assert partial is not None
            if j < m:
                sources[f"X{j}"] = partial[0]
                if j % 2 == 0:
                    sources[f"Y{j}"] = partial[1]
        else:
            constants = dict(a=a, Xm=partial[0], Ym=partial[1], p=p)
            for j, (alpha, beta, _) in enumerate(specs, 1):
                constants[f"alpha{j}"], constants[f"beta{j}"] = alpha, beta
            raw = c38.operations(sources, constants, m)
            assert all(value % p == 0 for value in raw[:c38.arithmetic_count(m)])
            lifted = dict(sources)
            for i, value in enumerate(raw[:c38.arithmetic_count(m)]):
                lifted[f"q{i}"] = value//p + 8*p
                assert 0 < lifted[f"q{i}"] < 16*p
            outputs = c38.operations(lifted, constants, m, mode="general")
            assert not any(outputs[:c38.arithmetic_count(m)])
            modulus = 1 << p.bit_length()
            for hashed, (_, _, threshold) in zip(outputs[c38.arithmetic_count(m):], specs):
                assert hashed % modulus < threshold
            counts[partial] += 1
    return counts


def rank(matrix, modulus):
    rows = [[x % modulus for x in row] for row in matrix]
    found = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(found, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[found], rows[pivot] = rows[pivot], rows[found]
        inv = pow(rows[found][column], -1, modulus)
        for i in range(found+1, len(rows)):
            scale = rows[i][column] * inv % modulus
            if scale:
                for j in range(column, len(rows[i])):
                    rows[i][j] = (rows[i][j] - scale * rows[found][j]) % modulus
        found += 1
    return found


def check_ranks(points, bases, add, p):
    n = len(points)
    composed = [[int(i == j) for j in range(n)] for i in range(n)]
    assert c38.prime(n+1)
    cumulative_size = 0
    results = []
    for slot in bases:
        if not slot:
            break
        matrix = []
        for u in points:
            possible = {add(u, factor) for factor in slot if u[0] != factor[0]}
            row = []
            for w in points:
                unrestricted = int(add(w, negate(u, p)) in slot)
                correction = int(u in slot and w == add(u, u))
                assert unrestricted - correction == int(w in possible)
                row.append(int(w in possible))
            matrix.append(row)
        s = len(slot)
        one_rank = rank(matrix, 1009)
        binary_rank = rank(matrix, 2)
        assert one_rank >= n-s
        assert binary_rank >= (n-s+s-1)//s
        composed = [[sum(x*y for x, y in zip(row, col)) for col in zip(*matrix)]
                    for row in composed]
        cumulative_size += s
        block_rank = rank(composed, 1009)
        assert block_rank >= n - cumulative_size
        results.append((s, one_rank, binary_rank, block_rank))
    return results


def check_chain(p, a, specs):
    points = c38.points(p, a, 1)
    group = (None, *points)
    def add(u, v):
        return c38.point_add(u, v, p, a)
    bases = [hash_slot(points, p.bit_length(), *spec) for spec in specs]
    count, work, calls, supports = sparse_guard_counts(bases, group, add, p)
    literal = literal_counts(bases, specs, add, p, a)
    assert count == {point: literal[point] for point in group}
    assert count[None] == 0
    ranks = check_ranks(points, bases, add, p) if c38.prime(len(group)) else []
    print(f"p={p}, group={len(group)}, slot sizes={list(map(len,bases))}, "
          f"valid tuples={sum(literal.values())}, correction supports={supports}, "
          f"oracle construction additions={work}, coefficient calls={calls}, "
          f"ranks (slot size, mod 1009, mod 2, product mod 1009)={ranks}")
    return count, supports


def check_coverage_scale():
    p, a, m, threshold = 211, 1, 4, 128
    points = c38.points(p, a, 1)
    r, modulus = len(points) + 1, 1 << p.bit_length()
    assert c38.prime(r) and r >= 30*m
    assert r*threshold >= 20*m*modulus
    assert (r*threshold)**m >= 2**m*r*modulus**m
    assert r*(threshold//2) < 20*m*modulus
    bases = [hash_slot(points, p.bit_length(), alpha, beta, threshold)
             for alpha, beta in ((7, 9), (11, 13), (59, 73), (19, 37))]
    def add(u, v):
        return c38.point_add(u, v, p, a)
    left = {add(u, v) for u in bases[0] for v in bases[1]}
    right = {add(u, v) for u in bases[2] for v in bases[3]}
    full = {add(u, v) for u in left for v in right}
    valid = set(bases[0])
    for slot in bases[1:]:
        valid = {add(u, v) for u in valid for v in slot if u[0] != v[0]}
        assert None not in valid
    assert valid <= full and len(valid) == r-1
    print(f"Coverage-scale fixture: p={p}, r={r}, m={m}, threshold={threshold}, "
          f"slot sizes={list(map(len, bases))}, half supports=({len(left)},{len(right)}), "
          f"full support={len(full)}, guarded affine support={len(valid)}")


if __name__ == "__main__":
    common.verify_sources()
    diagnostics = Counter()
    for p, a in ((7, 1), (11, 1), (13, 4)):
        diagnostics.update(check_chord_slice(p, a))
    assert all(diagnostics[key] for key in ("two_torsion_coincident", "linear",
                                            "roots_0", "roots_1", "roots_2"))
    print("Original E/D comparisons", dict(diagnostics), flush=True)
    check_deformed_join()
    cases = [
        (13, 4, [(1, 0, 4), (3, 2, 4), (5, 1, 4), (7, 0, 4), (1, 0, 4)]),
        (23, 5, [(1, 0, 4), (3, 0, 4), (5, 0, 4), (7, 0, 4)]),
        (13, 1, [(1, 0, 8), (3, 1, 8), (5, 2, 8), (1, 0, 8)]),
        (11, 1, [(1, 0, 4), (1, 0, 16), (3, 2, 4)]),
    ]
    for case in cases:
        counts, corrections = check_chain(*case)
        assert sum(counts.values()) and any(corrections)
    empty, _ = check_chain(13, 4, [(1, 0, 4), (1, 1, 1), (5, 2, 4)])
    assert not any(empty.values())
    check_coverage_scale()
