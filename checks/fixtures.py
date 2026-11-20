"""Independent arithmetic for the small curve fixtures used in this paper."""
from collections import Counter
from functools import lru_cache
from itertools import product


@lru_cache
def curve_cycle(p, a):
    """Enumerate y²=x³+ax+1 as multiples of the lexicographically first point."""
    roots = {}
    for y in range(p):
        roots.setdefault(y*y % p, []).append(y)
    points = [(x, y) for x in range(p) for y in roots.get((x**3+a*x+1) % p, ())]
    origin = points[0]
    gx, gy = origin
    cycle = [None, origin]
    seen = {origin}
    while True:
        x, y = cycle[-1]
        if x == gx and (y+gy) % p == 0:
            break
        numerator, denominator = ((3*x*x+a, 2*y) if x == gx else (gy-y, gx-x))
        slope = numerator*pow(denominator, -1, p) % p
        new_x = (slope*slope-gx-x) % p
        point = (new_x, (slope*(x-new_x)-y) % p)
        assert point not in seen
        seen.add(point)
        cycle.append(point)
    assert seen == set(points)
    return cycle


def hash_families(curve, modulus, width):
    masses = Counter()
    for alpha, beta in product(range(1, modulus, 2), range(modulus)):
        kept = tuple(i for i, (x, _) in enumerate(curve[1:], 1)
                     if (alpha*x+beta) % modulus < width)
        masses[kept] += 1
    return list(masses), list(masses.values())


def evaluation(point, dimension, p):
    x, y = point
    return ([pow(x, j, p) for j in range(dimension//2+1)]
            + [y*pow(x, j, p) % p for j in range((dimension-1)//2)])


def determinant(rows, p):
    """Laplace expansion, independent of the signed-orbit equations."""
    if len(rows) == 1:
        return rows[0][0] % p
    total = 0
    for j, coefficient in enumerate(rows[0]):
        if coefficient:
            minor = [row[:j]+row[j+1:] for row in rows[1:]]
            total += (-1)**j * coefficient * determinant(minor, p)
    return total % p
