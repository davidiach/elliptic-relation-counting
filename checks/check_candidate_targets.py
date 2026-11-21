"""Check translated, sign-free chord trees and their fixed-target candidate counts."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb

from check_block_support import group_tables
from check_strided_isolation import finite_laws


def polynomial_allows(x, u, z, p, linear, constant):
    if x is None:
        return z == u
    if u is None:
        return z == x
    if z is None:
        return x == u
    value = ((x-u)**2*z*z - 2*((x+u)*(x*u+linear)+2*constant)*z
             + (x*u-linear)**2 - 4*constant*(x+u))
    return value % p == 0


def orbit_tables(p, linear, constant):
    points, add, neg = group_tables(p, linear, constant)
    labels = [pt[0] if pt else None for pt in points]
    representatives = dict(zip(labels, range(len(points))))
    allowed = {}
    for x, i in representatives.items():
        for u, j in representatives.items():
            roots = {labels[add[i][j]], labels[add[i][neg[j]]]}
            allowed[x, u] = roots
    return points, add, neg, labels, allowed


def check_roots(p, linear, constant):
    points, add, _, xs, allowed = orbit_tables(p, linear, constant)
    for (x, u), real_roots in allowed.items():
        assert real_roots == {z for z in [None, *range(p)] if polynomial_allows(
            x, u, z, p, linear, constant)}
    leaf_checks = 0
    for i, pt in enumerate(points[1:], 1):
        for j, offset in enumerate(points[1:], 1):
            x, y = pt
            r, h = offset
            if not h or r == x:
                continue
            v = xs[add[i][j]]
            numerator = y*y + h*h - (x-r)**2*(v+x+r)
            assert numerator * pow(2*h, -1, p) % p == y
            leaf_checks += 1
    print(f"F_{p}, N={len(points)}: {len(allowed)} projective input pairs checked against "
          f"every projective root, {leaf_checks} translated signed-leaf identities.")


def check_pair_moment():
    p, linear, constant = 13, 4, 1
    points, add, neg = group_tables(p, linear, constant)
    n = len(points)
    positions = (0,)
    width, arity = 2, 3
    affine, fixed = finite_laws(p, points, positions)
    all_targets = [(t, q) for t in range(1, n) if t != neg[t]
                   for q in range(n) if q not in (t, neg[t])]
    joint = Counter()
    mass = 0
    for a, b, c in product(affine, fixed, fixed):
        weight = affine[a]*fixed[b]*fixed[c]
        mass += weight
        h = Counter(add[add[i][j]][k] for (i, _), (j, _), (k, _) in product(a, b, c))
        for t, q in all_targets:
            joint[t, q] += weight*h[t]*h[q]
    rho = Fraction(width, 1 << p.bit_length())
    nu = n*rho
    bound = rho**2*(2*width*nu)**(arity-1)*(1 + 4*(arity-1)*2/nu)
    largest = Fraction(max(joint.values()), mass)
    assert largest <= bound
    print(f"F_{p}: {mass:,} weighted mixed seeds, {len(all_targets)} ordered target pairs; "
          f"largest cross moment {largest}, proved upper bound {bound}.")


def checked_offset(add, neg, arity):
    for t in range(1, len(add)):
        multiples = [0]
        for _ in range(2*arity+1):
            multiples.append(add[multiples[-1]][t])
        qs = [multiples[1+2*j] for j in range(1, arity+1)]
        if all(q not in (t, neg[t]) for q in qs):
            return t, qs
    raise AssertionError("no eligible offset")


def check_singletons():
    p, linear, constant, arity = 257, 1, 1, 3
    points, add, neg, xs, allowed = orbit_tables(p, linear, constant)
    target, shifts = checked_offset(add, neg, arity)
    target_shift = target
    for _ in range(arity):
        target_shift = add[target_shift][target]
    root = xs[target_shift]
    signed_translates = [(shifted, neg[shifted]) for shifted in (add[i][target] for i in range(len(points)))]
    by_x = defaultdict(list)
    for i in range(1, len(points)):
        by_x[xs[i]].append(i)
    hits = {q: Counter() for q in [target, *shifts]}
    guarded = set()
    for i, j in product(range(1, len(points)), repeat=2):
        subtotal = add[i][j]
        for q, counts in hits.items():
            k = add[neg[subtotal]][q]
            if not k:
                continue
            code = (xs[i], xs[j], xs[k])
            counts[code] += 1
            if q == target and subtotal and xs[i] != xs[j] and xs[subtotal] != xs[k]:
                guarded.add(code)
    singletons = false_addresses = surplus_assignments = surplus_signs = 0
    for code in guarded:
        if hits[target][code] != 1:
            continue
        singletons += 1
        streams = [by_x[x] for x in code]
        wrong_here = surplus_here = sign_hits = 0
        for i, j, k in product(*streams):
            u, v, w = (xs[add[index][target]] for index in (i, j, k))
            assignments = sum(root in allowed[mid, w] for mid in allowed[u, v])
            true = add[add[i][j]][k] == target
            assert not true or assignments
            wrong_here += bool(assignments) and not true
            surplus_here += assignments - int(true)
            sign_hits += sum(add[add[a][b]][c] == target_shift for a, b, c in
                             product(signed_translates[i], signed_translates[j], signed_translates[k]))
        sign_bound = sum(comb(arity, j)*hits[q][code] for j, q in enumerate(shifts, 1))
        assert sign_hits == sign_bound + 1
        assert 0 <= wrong_here <= surplus_here <= sign_bound
        false_addresses += wrong_here
        surplus_assignments += surplus_here
        surplus_signs += sign_bound
    assert 0 < false_addresses <= surplus_assignments <= surplus_signs
    assert 3*surplus_signs < singletons
    print(f"F_{p}, N={len(points)}, width one, three slots, fixed R=T: {singletons:,} "
          f"guarded singleton abscissa patterns; {false_addresses:,} wrong accepted addresses, "
          f"{surplus_assignments:,} surplus internal assignments, {surplus_signs:,} extra sign patterns.")
    print(f"  Conditional means {false_addresses/singletons:.6f}, "
          f"{surplus_assignments/singletons:.6f}, {surplus_signs/singletons:.6f}, respectively.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_roots(13, 4, 1)
    check_roots(19, -1, 0)
    check_pair_moment()
    check_singletons()
