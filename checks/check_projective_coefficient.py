"""Compare a truncated coefficient product with signed determinants over small fields."""
from collections import Counter
from itertools import product
from math import prod
from random import Random

from flint import nmod_poly

from fixtures import determinant
from probe_signed_moment import Probe


def multiply(left, right, p):
    answer = {}
    for a, x in left.items():
        for b, y in right.items():
            term = tuple(u+v for u, v in zip(a, b))
            if max(term) < p:
                answer[term] = (answer.get(term, 0)+x*y) % p
    return {term: value for term, value in answer.items() if value}


def power(row, p):
    m = len(row)
    linear = {tuple(int(k == j) for k in range(m)): a for j, a in enumerate(row) if a}
    value = {(0,)*m: 1}
    for _ in range(p-1):
        value = multiply(value, linear, p)
    return value


def coefficient(slots, p, m):
    value = {(0,)*m: 1}
    for slot in slots:
        factor = {}
        for row, weight in slot:
            for term, a in power(row, p).items():
                factor[term] = (factor.get(term, 0)+weight*a) % p
        value = multiply(value, factor, p)
    return value.get((p-1,)*m, 0)


def directions(p, m):
    return [(0,)*i+(1,)+tail for i in range(m) for tail in product(range(p), repeat=m-i-1)]


def check(slots, p, m):
    signed_full_rank = 0
    for choices in product(*slots):
        rows, weights = zip(*choices)
        signed_full_rank += prod(weights)*pow(determinant(rows, p), p-1, p)
    assert coefficient(slots, p, m) == signed_full_rank % p
    return (prod(sum(w for _, w in slot) for slot in slots)-signed_full_rank) % p


def section_points(coefficients, endpoint, p, curve_poly):
    constant = -sum(c*u for c, u in zip(coefficients, endpoint))
    row = [constant, *coefficients]
    split = len(row)//2+1
    a = nmod_poly(row[:split], p)
    b = nmod_poly(row[split:], p)
    norm = a*a-curve_poly*b*b
    assert norm and norm.degree() <= len(row)
    points = set()
    for root, _ in norm.roots():
        x = int(root)
        if b(x):
            points.add((x, -int(a(x))*pow(int(b(x)), -1, p) % p))
        else:
            assert a(x) == 0
            points.update((x, int(y)) for y, _ in nmod_poly([-int(curve_poly(x)), 0, 1], p).roots())
    return points


def hidden_seed_check(probe, target, tail, sections):
    visible_x = {probe.curve[i][0] for slot in tail for i, _ in slot}
    known = [{probe.curve[i] for i, _ in slot} for slot in tail]
    parents = Counter()
    hits = Counter()
    first_feedback = Counter()
    seeds = probe.modulus**2//2
    for alpha in range(1, probe.modulus, 2):
        for beta in range(probe.modulus):
            anchor = [(i, 2*u+sign) for x, points in probe.by_x.items()
                      if (u := (alpha*x+beta) % probe.modulus) < probe.width for i, sign in points]
            if any(probe.curve[i][0] in visible_x | {probe.curve[target][0]} for i, _ in anchor):
                continue
            address, _ = probe.singleton([anchor, *tail], target)
            if address is None:
                continue
            parents[address] += 1
            selected = {probe.curve[i] for i, _ in anchor}
            for node, points in sections:
                if selected & points:
                    first_feedback[node] += 1
                    if all(row & points for row in known):
                        hits[node] += 1
    marginal = seeds*probe.width//probe.modulus
    assert 0 < max(parents.values()) <= marginal
    assert sum(hits.values()) == sum(parents.values())
    assert 0 < max(hits.values()) <= marginal
    assert 0 < max(first_feedback.values()) <= 3*marginal
    print(f"All {seeds} hidden anchor seeds with visible tails fixed: {sum(parents.values())} regular guarded singletons; "
          f"largest address/section masses {max(parents.values())}/{max(hits.values())} seeds (bound {marginal}).")


def curve_checks():
    probe = Probe(13, 4, (0,))
    m, target = 3, 1
    rng = Random(494)
    cases = []
    for _ in range(2000):
        slots = probe.planted(target, rng)
        xs = [{probe.curve[i][0] for i, _ in slot} for slot in slots]
        if any(probe.curve[target][0] in a for a in xs):
            continue
        if any(xs[i] & xs[j] for i in range(m) for j in range(i)):
            continue
        witness, _ = probe.singleton(slots, target)
        if witness is not None:
            cases.append((slots, witness))
            break
    assert cases
    slots, witness = cases[0]
    endpoint = probe.vectors[-target][1:]
    curve_poly = nmod_poly([1, 4, 0, 1], probe.p)
    # Recover a valid description of this sampled fixture, then use only the three
    # seed words for the incidence tests after each section has been factored.
    anchor = next((alpha, beta) for alpha in range(1, probe.modulus, 2) for beta in range(probe.modulus)
                  if set(slots[0]) == {(i, 2*u+sign) for x, points in probe.by_x.items()
                                        if (u := (alpha*x+beta) % probe.modulus) < probe.width
                                        for i, sign in points})
    fixed = [probe.curve[slot[0][0]][0] & ~1 for slot in slots[1:]]
    supported, sections = [], []
    for v in directions(probe.p, m):
        points = section_points(v, endpoint, probe.p, curve_poly)
        sections.append((v, points))
        expected = {probe.curve[i] for i in range(1, probe.n)
                    if sum(c*(a-b) for c, a, b in zip(v, probe.vectors[i][1:], endpoint)) % probe.p == 0}
        assert points == expected
        memberships = [[] for _ in slots]
        for x, y in points:
            sign = int(y > (-y) % probe.p)
            u = (anchor[0]*x+anchor[1]) % probe.modulus
            if u < probe.width:
                memberships[0].append(2*u+sign)
            for j, prefix in enumerate(fixed, 1):
                if x & ~1 == prefix:
                    memberships[j].append(2*(x & 1)+sign)
        for row, slot in zip(memberships, slots):
            assert sorted(row) == sorted(addr for i, addr in slot if probe.curve[i] in expected)
        if all(memberships):
            assert all(len(row) == 1 for row in memberships)
            supported.append(sum(row[0]*probe.block**(m-1-i) for i, row in enumerate(memberships)))
    assert supported == [witness]
    # The quartic-slot basis includes xy. This section has A(x)=B(x)=0
    # at the target abscissa and must recover both rational lifts there.
    x0, y0 = probe.curve[-target]
    points = section_points((-y0, 0, -x0, 1), [x0, x0*x0, y0, x0*y0], probe.p, curve_poly)
    assert points == {pt for pt in probe.curve[1:] if pt[0] == x0 or pt[1] == y0}
    assert (x0, -y0 % probe.p) in points
    rows = [[(tuple((u-v) % probe.p for u, v in zip(probe.vectors[i][1:], endpoint)), addr)
             for i, addr in slot] for slot in slots]
    basis = {row: power(row, probe.p) for slot in rows for row, _ in slot}
    for challenge in range(probe.size):
        value = {(0,)*m: 1}
        mass = 1
        for i, slot in enumerate(rows):
            component = challenge // probe.block**(m-i-1) % probe.block
            factor = {}
            slot_sum = 0
            for row, address in slot:
                weight = (-1)**((address & component).bit_count())
                slot_sum += weight
                for term, a in basis[row].items():
                    factor[term] = (factor.get(term, 0)+weight*a) % probe.p
            mass *= slot_sum
            value = multiply(value, {t: a for t, a in factor.items() if a}, probe.p)
        assert (mass-value.get((probe.p-1,)*m, 0)) % probe.p == (-1)**((challenge & witness).bit_count()) % probe.p
    print(f"Regular sampled singleton over F_{probe.p}: all {probe.size} address parities agree with capped coefficient extraction.")
    print(f"Factoring all {len(directions(probe.p, m))} sections and testing the three seeds finds one witness section; shared-root control passed.")
    hidden_seed_check(probe, target, slots[1:], sections)


def main():
    rng = Random(2356)
    for p, m in ((3, 2), (3, 3), (5, 3), (7, 2)):
        for _ in range(40):
            slots = [[(tuple(rng.randrange(p) for _ in range(m)), rng.choice((-1, 1)))
                      for _ in range(rng.randrange(1, 4))] for _ in range(m)]
            check(slots, p, m)
        nodes = directions(p, m)
        for u in nodes:
            pivot = next(i for i, v in enumerate(u) if v)
            cov = tuple(int(j == pivot) for j in range(m))
            slots = [[(cov, 1)]]
            for k in range(m):
                if k == pivot:
                    continue
                row = tuple((int(j == k)-u[k]*cov[j]) % p for j in range(m))
                slots.append([(cov, 1), (row, -1)])
            assert coefficient(slots, p, m) == (-1)**(m+1) % p
            for v in nodes:
                value = prod(sum(w*pow(sum(a*b for a, b in zip(row, v)), p-1, p)
                                 for row, w in slot) for slot in slots) % p
                assert value == int(u == v)
        print(f"F_{p}, m={m}: 40 signed grids and all {len(nodes)} projective point indicators verified.")
    curve_checks()


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
