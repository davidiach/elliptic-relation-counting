"""Compare the product-free chord test with literal products modulo p squared."""
from collections import Counter, defaultdict
from itertools import product
from math import prod
from random import Random

from flint import fmpz_mat, nmod_poly

from check_nonlinear_branches import mixed_fixture
from check_section_deformation import lift
from fixtures import determinant, evaluation
from probe_signed_moment import Probe


def negative(point, p):
    return None if point is None else (point[0], -point[1] % p)


def add_with_slope(a, b, p, linear):
    if a is None or b is None:
        return (b if a is None else a), 0
    x, y = a
    u, v = b
    if x == u and (y+v) % p == 0:
        return None, 0
    slope = ((3*x*x+linear)*pow(2*y, -1, p) if x == u else (y-v)*pow(x-u, -1, p)) % p
    new_x = (slope*slope-x-u) % p
    return (new_x, (slope*(x-new_x)-y) % p), slope


def half_table(slots, p, linear):
    table = [(None, 0, 0)]
    for slot in slots:
        next_table = []
        for point, weight, accumulated in table:
            for item, speed in slot:
                result, slope = add_with_slope(point, item, p, linear)
                next_table.append((result, (weight+speed) % p, (accumulated+slope) % p))
        table = next_table
    return table


def rational_tree(entries, p):
    if len(entries) == 1:
        x, weight = entries[0]
        return nmod_poly([-x, 1], p), nmod_poly([weight], p)
    split = len(entries)//2
    b, n = rational_tree(entries[:split], p)
    d, r = rational_tree(entries[split:], p)
    return b*d, n*d+r*b


def evaluation_tree(xs, p):
    if len(xs) == 1:
        return nmod_poly([-xs[0], 1], p), xs[0], None
    split = len(xs)//2
    left = evaluation_tree(xs[:split], p)
    right = evaluation_tree(xs[split:], p)
    return left[0]*right[0], left, right


def values(poly, node, out):
    modulus, left, right = node
    remainder = poly % modulus
    if right is None:
        out[left] = int(remainder[0])
    else:
        values(remainder, left, out)
        values(remainder, right, out)


def finite_parts(weights, xs, tree, p):
    if not weights:
        return dict.fromkeys(xs, 0)
    den, num = rational_tree(list(weights.items()), p)
    tables = []
    for f in (den, num, den.derivative(), num.derivative(), den.derivative().derivative()):
        table = {}
        values(f, tree, table)
        tables.append(table)
    d, n, dp, np, dpp = tables
    out = {}
    for x in xs:
        if x in weights:
            out[x] = (2*np[x]-weights[x]*dpp[x])*pow(2*dp[x], -1, p) % p
        else:
            out[x] = n[x]*pow(d[x], -1, p) % p
    return out


def cross_slope(left, right, p, linear):
    """Sum (u+v)*slope(P,Q), with duplicate and equal-abscissa corrections."""
    xs = sorted({pt[0] for pt, u, c in left if pt is not None})
    if not xs:
        return 0
    inputs = [defaultdict(int) for _ in range(4)]
    counts, summed = Counter(), Counter()
    for pt, v, _ in right:
        if pt is None:
            continue
        x, y = pt
        counts[pt] += 1
        summed[pt] += v
        for table, coefficient in zip(inputs, (1, y, v, v*y)):
            table[x] = (table[x]+coefficient) % p
    tree = evaluation_tree(xs, p)
    f, fy, g, gy = [finite_parts(weights, xs, tree, p) for weights in inputs]
    result = 0
    for pt, u, _ in left:
        if pt is not None:
            x, y = pt
            result += u*(y*f[x]-fy[x])+y*g[x]-gy[x]
            if y and counts[pt]:
                result += (u*counts[pt]+summed[pt])*(3*x*x+linear)*pow(2*y, -1, p)
    return result % p


def normalized_difference(slots, p, linear):
    if not all(slots):
        return 0, 0
    # The last slot is the fixed endpoint, so split only the varying slots.
    split = (len(slots)-1)//2
    left = half_table(slots[:split], p, linear)
    right = half_table(slots[split:], p, linear)
    counts, weights = Counter(), Counter()
    for pt, v, _ in right:
        counts[pt] += 1
        weights[pt] += v
    hits, speed = 0, 0
    for pt, u, _ in left:
        other = negative(pt, p)
        hits += counts[other]
        speed += u*counts[other]+weights[other]
    if hits:
        return hits, 0 if hits >= 2 else speed % p
    def totals(half):
        return sum(u for _, u, c in half), sum(c for _, u, c in half), sum(u*c for _, u, c in half)
    lu, lc, luc = totals(left)
    ru, rc, ruc = totals(right)
    val = len(right)*luc+len(left)*ruc+lu*rc+ru*lc+cross_slope(left, right, p, linear)
    size = prod(map(len, slots))
    for i, a in enumerate(slots):
        first = [(pt, u, 0) for pt, u in a]
        for b in slots[:i]:
            second = [(negative(pt, p), -v, 0) for pt, v in b]
            val += (size//(len(a)*len(b)))*cross_slope(first, second, p, linear)
    return 0, val % p


def literal(slots, p, linear, offset):
    modulus = p*p
    lifted = [[evaluation(lift(pt, u, linear, p, offset), len(slots), modulus)
               for pt, u in s] for s in slots]
    value = 1
    for choice in product(*lifted):
        det = determinant(list(choice), modulus) if len(slots) <= 4 else int(fmpz_mat(list(choice)).det())
        value = value*det % modulus
    return value


def check_raw_resolvent(probe, slots, target, moving_slot, hits, same):
    p = probe.p
    x = probe.curve[target][0]
    poles = Counter()
    for chosen in product(*slots):
        point, digit = chosen[moving_slot]
        if digit >> 1:
            poles[probe.curve[point][0]] -= 1
            s = sum(i for i, _ in chosen) % probe.n
            if s:
                poles[probe.curve[s][0]] += 1
    if hits == 0:
        assert poles[x] == 0
        scaled = sum(k*pow(x-z, -1, p) for z, k in poles.items() if z != x)
        assert same == (scaled % p == 0)
    elif hits == 1:
        assert poles[x] in (0, 2) and same == (poles[x] == 0)


def run(p, linear):
    probe = Probe(p, linear, (0,))
    original = mixed_fixture(probe)
    used = {probe.curve[i][0] for slot in original for i, _ in slot}
    population = Counter(sum(i for i, a in c) % probe.n for c in product(*original))
    eligible = [t for t in range(1, probe.n) if probe.curve[t][0] not in used]
    test_targets = eligible[:22]
    for k in (0, 1, 2):
        example = next((v for v in eligible if min(2, population[v]) == k), None)
        if example is not None and example not in test_targets:
            test_targets.append(example)
    rng = Random(908)
    seen = Counter()
    for target in test_targets:
        data = [[(probe.curve[i], a) for i, a in slot] for slot in original]+[[(probe.curve[-target], 0)]]
        base = [[(pt, 0) for pt, _ in s] for s in data]
        offset = rng.randrange(p)
        reference = literal(base, p, linear, offset)
        endpoint_move = literal([*base[:-1], [(data[-1][0][0], 1)]], p, linear, offset)
        for k in range(10):
            if k < 6:
                moving = [[(pt, (a >> (k%2)) & 1 if j == k//2 else 0) for pt, a in s]
                          for j, s in enumerate(data)]
            elif k == 6:
                moving = base
            else:
                moving = [[(pt, rng.randrange(p)) for pt, _ in s] for s in data]
            hits, scalar = normalized_difference(moving, p, linear)
            moved = literal(moving, p, linear, offset)
            same = scalar == 0
            assert same == (moved == reference)
            if hits < 2:
                unit = reference if hits == 0 else (endpoint_move-reference)//p
                assert scalar == ((moved-reference)//p)*pow(unit, -1, p) % p
            if k in (1, 3, 5):
                check_raw_resolvent(probe, original, target, k//2, hits, same)
            seen[min(2, hits), same] += 1
    assert all(seen[k] for k in ((0, True), (0, False), (1, True), (1, False)))
    assert seen[2, False] == 0 and (seen[2, True] or all(population[t] < 2 for t in eligible))
    print(f"F_{p}: {sum(seen.values())} literal comparisons, equality/inequality counts "
          f"for zero hits {seen[0,True]}/{seen[0,False]}, one hit {seen[1,True]}/{seen[1,False]}, "
          f"multiple hits {seen[2,True]}/{seen[2,False]}.")
    print(f"  Verified the raw-bit reciprocal formula or singleton residue on "
          f"{3*sum(min(2,population[t]) < 2 for t in test_targets)} comparisons.")
    print(f"  Matched the full normalized field digit on all {sum(n for (h, _), n in seen.items() if h < 2)} "
          "empty and singleton comparisons.")


def exceptional_chords():
    p, linear = 17, 1
    points = [(x, y) for x, y in product(range(p), repeat=2) if (y*y-x**3-linear*x-1) % p == 0]
    points = [None, *points, points[0], points[1]]
    rng = Random(34)
    for _ in range(12):
        left = [(pt, rng.randrange(p), 0) for pt in points]
        right = [(pt, rng.randrange(p), 0) for pt in reversed(points)]
        actual = sum((u+v)*add_with_slope(a, b, p, linear)[1] for a, u, _ in left for b, v, _ in right)
        assert actual % p == cross_slope(left, right, p, linear)
    assert any(pt and not pt[1] for pt in points)
    print("F_17: 12 full chord transforms with identities, opposite and repeated points, tangents and rational two-torsion.")


def six_slots():
    p, linear = 1031, 9
    probe = Probe(p, linear, (0,))
    slots = []
    for x in range(0, p-1, 2):
        slot = [(i, 2*j+bit) for j in (0, 1) for i, bit in probe.by_x.get(x+j, ())]
        if len(slot) == 2:
            slots.append(slot)
        if len(slots) == 6:
            break
    assert len(slots) == 6
    used = {probe.curve[i][0] for slot in slots for i, _ in slot}
    population = Counter(sum(i for i, _ in chosen) % probe.n for chosen in product(*slots))
    by_count = {}
    for t in range(1, probe.n):
        if probe.curve[t][0] not in used:
            by_count.setdefault(min(2, population[t]), t)
    assert 0 in by_count and 1 in by_count
    rng = Random(337)
    checked = 0
    for h, t in sorted(by_count.items()):
        data = [[probe.curve[i] for i, _ in slot] for slot in slots]+[[probe.curve[-t]]]
        base = [[(pt, 0) for pt in slot] for slot in data]
        reference = literal(base, p, linear, 15)
        unit = reference if h == 0 else (literal([*base[:-1], [(data[-1][0], 1)]], p, linear, 15)-reference)//p
        for _ in range(10):
            moving = [[(pt, rng.randrange(p)) for pt in slot] for slot in data]
            hits, scalar = normalized_difference(moving, p, linear)
            moved = literal(moving, p, linear, 15)
            assert min(2, hits) == h and (scalar == 0) == (moved == reference)
            if h < 2:
                assert scalar == ((moved-reference)//p)*pow(unit, -1, p) % p
            checked += 1
    print(f"F_{p}, six width-two slots: {checked} comparisons with independent integer determinants, "
          f"hit-count classes {sorted(by_count)} (2 denotes two or more).")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run(101, 17)
    run(1031, 9)
    exceptional_chords()
    six_slots()
