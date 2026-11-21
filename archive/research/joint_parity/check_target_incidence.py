"""Exact finite checks for the determinant and fresh-hash claims in target_incidence.tex."""

from collections import Counter
from itertools import combinations, product
from math import prod
from pathlib import Path
from random import Random
from runpy import run_path

fixture = run_path(str(Path(__file__).resolve().parents[1] / "slot_arity" / "verify.py"))
curve_cycle, window = fixture["curve_cycle"], fixture["window"]


def evaluation(point, dimension, p):
    x, y = point
    return ([pow(x, j, p) for j in range(dimension // 2 + 1)] +
            [y * pow(x, j, p) % p for j in range((dimension - 1) // 2)])


def determinant(rows, p):
    matrix = [list(row) for row in rows]
    value = 1
    for col in range(len(matrix)):
        pivot = next((i for i in range(col, len(matrix)) if matrix[i][col] % p), None)
        if pivot is None:
            return 0
        if pivot != col:
            matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
            value = -value
        diagonal = matrix[col][col] % p
        value = value * diagonal % p
        inverse = pow(diagonal, -1, p)
        for i in range(col + 1, len(matrix)):
            factor = matrix[i][col] * inverse % p
            for j in range(col + 1, len(matrix)):
                matrix[i][j] = (matrix[i][j] - factor * matrix[col][j]) % p
    return value


def check_divisors():
    p = 13
    curve = curve_cycle(p, 4)
    order = len(curve)
    checked = singular = 0
    for d in range(3, 7):
        rows = {k: evaluation(curve[k], d, p) for k in range(1, order)}
        for indices in combinations(range(1, (order + 1) // 2), d):
            for signs in product((1, -1), repeat=d):
                chosen = [(k * s) % order for k, s in zip(indices, signs)]
                det = determinant([rows[k] for k in chosen], p)
                assert (det == 0) == (sum(chosen) % order == 0)
                checked += 1
                singular += det == 0
    assert 0 < singular < checked
    print(f"Divisor test: {checked:,} signed configurations, {singular:,} singular.")


def gram(rows, weights, p):
    size = len(rows[0])
    return [[sum(w * row[i] * row[j] for row, w in zip(rows, weights)) % p
             for j in range(size)] for i in range(size)]


def mixed_square(slots, weights, fixed, p):
    # Subtract the fixed row from each chosen row and expand along the constant column.
    m = len(slots)
    grams = []
    for rows, slot_weights in zip(slots, weights):
        projected = [[(v - t) % p for v, t in zip(row[1:], fixed[1:])] for row in rows]
        grams.append(gram(projected, slot_weights, p))
    total = 0
    for flags in product((0, 1), repeat=m):
        matrix = [[sum(g[i][j] for flag, g in zip(flags, grams) if flag)
                   for j in range(m)] for i in range(m)]
        total += (-1)**(m - sum(flags)) * determinant(matrix, p)
    return total % p


def target_form(slots, weights, p):
    d = len(slots) + 1
    grams = [gram(rows, w, p) for rows, w in zip(slots, weights)]
    form = [[0] * d for _ in range(d)]
    for flags in product((0, 1), repeat=len(slots)):
        matrix = [[sum(g[i][j] for flag, g in zip(flags, grams) if flag)
                   for j in range(d)] for i in range(d)]
        sign = (-1)**(len(slots) - sum(flags))
        for i in range(d):
            for j in range(d):
                minor = [[v for col, v in enumerate(row) if col != i]
                         for r, row in enumerate(matrix) if r != j]
                form[i][j] += sign * (-1)**(i+j) * determinant(minor, p)
    return form


def check_fresh_seeds():
    p, a, modulus, width = 13, 4, 16, 2
    curve = curve_cycle(p, a)
    order, d = len(curve), 3
    masks = Counter()
    for alpha, beta in product(range(1, modulus, 2), range(modulus)):
        present = window(curve, modulus, width, alpha, beta)
        masks[tuple(i for i, hit in enumerate(present) if hit)] += 1
    rows = [None] + [evaluation(point, d, p) for point in curve[1:]]
    totals = Counter()
    ambiguous = {}
    witness = None
    for left, right in product(masks, repeat=2):
        mass = masks[left] * masks[right]
        zero_count = 0
        for t in range(1, order):
            fixed = rows[-t]
            pairs = list(product(left, right))
            minors = [determinant([rows[u], rows[v], fixed], p) for u, v in pairs]
            z2 = sum(det**2 for det in minors) % p
            weighted = sum((u+1)*(v+2)*det**2 for (u, v), det in zip(pairs, minors)) % p
            if left and right:
                selected = [[rows[u] for u in left], [rows[v] for v in right]]
                assert z2 == mixed_square(selected, [[1]*len(left), [1]*len(right)], fixed, p)
                assert weighted == mixed_square(selected, [[u+1 for u in left],
                                                            [v+2 for v in right]], fixed, p)
            else:
                assert z2 == weighted == 0
            h = sum((u+v) % order == t for u, v in pairs)
            good = (len({curve[j][0] for j in left} & {curve[j][0] for j in right}) == 0
                    and all(curve[j][0] != curve[t][0] for j in (*left, *right)))
            if good:
                assert sum(det == 0 for det in minors) == h
                assert (len(pairs) - sum(pow(det, p-1, p) for det in minors)) % p == h % p
                key = (t, len(left), len(right), z2)
                if key in ambiguous and ambiguous[key] != bool(h) and witness is None:
                    witness = key
                ambiguous[key] = bool(h)
            totals["trials"] += mass
            totals["regular"] += mass * good
            totals["positive"] += mass * bool(h)
            totals["zero_bit_wrong"] += mass * ((z2 == 0) != bool(h))
            zero_count += z2 == 0
        assert zero_count <= 2*d or zero_count == order-1
    assert witness is not None and totals["trials"] == 128**2 * (order - 1)
    print(f"Actual affine seeds: {totals['trials']:,} seed-pair/target cases, "
          f"{totals['regular']:,} regular, {totals['positive']:,} positive; "
          f"unweighted zero test wrong in {totals['zero_bit_wrong']:,}.")
    print("Both answers occur on regular inputs at the same "
          f"(target index, slot sizes, squared moment) = {witness}.")


def check_multi_slot_moments():
    p, a = 23, 5
    curve = curve_cycle(p, a)
    checked = 0
    for m in (3, 4, 5):
        d = m + 1
        slots = [[evaluation(curve[k], d, p) for k in (i, -i)] for i in range(1, m+1)]
        weights = [[(2*i+1) % p, (3*i+4) % p] for i in range(m)]
        form = target_form(slots, weights, p)
        for target in range(1, len(curve)):
            fixed = evaluation(curve[-target], d, p)
            brute = 0
            for choices in product((0, 1), repeat=m):
                selected = [rows[c] for rows, c in zip(slots, choices)]
                factor = prod(w[c] for w, c in zip(weights, choices))
                brute += factor * determinant([*selected, fixed], p)**2
            assert brute % p == mixed_square(slots, weights, fixed, p)
            assert brute % p == sum(fixed[i] * form[i][j] * fixed[j]
                                    for i in range(d) for j in range(d)) % p
            checked += 1
    print(f"Weighted mixed coefficient agrees with direct expansion at {checked} more targets.")


def check_sampled_norm():
    p, a = 101, 17
    curve = curve_cycle(p, a)
    target = (2, 89)
    xs = (0, 47, 71, 21, 3)
    signed_indices = [[i for i in range(1, len(curve)) if curve[i][0] == x] for x in xs]
    target_index = curve.index(target)
    assert all(len(pair) == 2 for pair in signed_indices)
    assert all(sum(choice) % len(curve) != target_index for choice in product(*signed_indices))

    def residual(x):
        u, v = x**3 + 82*x + 18, 82*x + 39
        return (u*u - (x**3 + a*x + 1)*v*v - (x-target[0])*prod(x-i for i in xs)) % p

    x, y = target
    assert (82*x + 39) % p == 1 and (x**3 + 82*x + 18) % p == y
    assert residual(1) == 0 and any(residual(z) for z in range(p))
    print("Sampled norm passes the target and one extra node; all 32 sign choices miss.")


def solve_many(matrix, right_hand_sides, p):
    n = len(matrix)
    rows = [list(row) + [rhs[i] for rhs in right_hand_sides] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j] % p)
        rows[j], rows[pivot] = rows[pivot], rows[j]
        inverse = pow(rows[j][j] % p, -1, p)
        rows[j] = [v * inverse % p for v in rows[j]]
        for i in range(n):
            if i != j and rows[i][j]:
                factor = rows[i][j]
                rows[i] = [(a - factor*b) % p for a, b in zip(rows[i], rows[j])]
    return [list(col) for col in zip(*(row[n:] for row in rows))]


def check_higher_moments():
    p, a, m, k = 101, 17, 3, 5
    curve = curve_cycle(p, a)
    d, n = m + 1, (m + 1) * k
    candidates = curve[1:n+2]
    evaluations = [evaluation(point, n, p) for point in candidates]
    basis_rows = fixture['basis'](evaluations, p)
    assert len(basis_rows) == n
    nodes = [candidates[evaluations.index(row)] for row in basis_rows]
    slot_indices = []
    for alpha, beta in ((1, 4), (7, 41), (13, 17)):
        present = window(curve, 128, 32, alpha, beta)
        slot_indices.append([i for i, included in enumerate(present) if included])
    slots = [[curve[i] for i in indices] for indices in slot_indices]
    weights = [[(j+2)*i % p for i in indices] for j, indices in enumerate(slot_indices)]
    rhs = []
    for slot, ws in zip(slots, weights):
        lifted = [evaluation(point, n, p) for point in slot]
        rhs.append([sum(w*v for w, v in zip(ws, col)) % p for col in zip(*lifted)])
    node_weights = solve_many(list(zip(*basis_rows)), rhs, p)
    fixed = evaluation(curve[-11], d, p)

    def moments(point_sets, coefficients):
        records = [[(evaluation(point, d, p), w) for point, w in zip(points, ws)]
                   for points, ws in zip(point_sets, coefficients)]
        result = [0] * (k+1)
        for choice in product(*records):
            rows, ws = zip(*choice)
            delta = determinant([*rows, fixed], p)
            term = prod(ws) % p
            for j in range(k+1):
                result[j] = (result[j] + term) % p
                term = term * delta % p
        return result

    projected = moments([nodes]*m, node_weights)
    assert projected == moments(slots, weights) and any(projected[1:])
    print(f"Interpolation: moments 0 through {k} agree on {n**m:,} public node tuples "
          f"and {prod(map(len, slots)):,} original tuples.")


def check_batched_chords():
    # Use the same five-abscissa negative instance as the signed norm check.
    vendor = Path(__file__).resolve().parents[1] / "hillclimb_001" / "vendor"
    chord = run_path(str(vendor / "arithmetic_alternating_chord.py"))
    p, m = 101, 5
    xs = (0, 47, 71, 21, 3)
    public = dict(Xm=2, Ym=89, a=17)
    words = {name: 1 for name in chord["source_names"](m)}
    for i, x in enumerate(xs, 1):
        words[f'x{i}'] = x
        public[f'alpha{i}'], public[f'beta{i}'] = 1, (-x) % 128
    e = chord["arithmetic_count"](m)

    def residual(y2, y4):
        assignment = dict(words, Y2=y2, Y4=y4)
        return chord["operations"](assignment, public, m)[:e]

    base = residual(0, 0)
    columns = [[(new-old) % p for new, old in zip(residual(*ys), base)]
               for ys in ((1, 0), (0, 1))]
    assert [[column[row] for column in columns] for row in (1, 5)] == [[2, 0], [0, 2]]
    random = Random(43270)
    accepted = 0
    for _ in range(250):
        batch = [[random.randrange(p) for _ in range(e)] for _ in range(2)]
        matrix = [[sum(a*b for a, b in zip(row, col)) % p for col in columns]
                  for row in batch]
        det = determinant(matrix, p)
        if det == 0:
            continue
        rhs = [-sum(a*b for a, b in zip(row, base)) % p for row in batch]
        u, v = matrix
        y2 = (rhs[0]*v[1] - u[1]*rhs[1]) * pow(det, -1, p) % p
        y4 = (u[0]*rhs[1] - rhs[0]*v[0]) * pow(det, -1, p) % p
        fabricated = residual(y2, y4)
        assert any(value % p for value in fabricated)
        assert all(sum(a*b for a, b in zip(row, fabricated)) % p == 0 for row in batch)
        accepted += 1
    assert accepted > 230
    print(f"Random chord checks: {accepted} of 250 two-row batches had full rank and "
          "accepted an explicit false witness.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_divisors()
    check_fresh_seeds()
    check_multi_slot_moments()
    check_higher_moments()
    check_sampled_norm()
    check_batched_chords()
