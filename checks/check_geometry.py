"""Verify the generic-degree examples and the chord quartic on fixed finite curves."""
from collections import Counter
from itertools import product

from fixtures import curve_cycle
from check_batch import add


def check_degree():
    points = curve_cycle(1031, 9)
    order = len(points)
    counts = []
    for arity in range(3, 9):
        projections = set()
        for signs in product((-1, 1), repeat=arity-1):
            index = (250 - sum(s*(1 << i) for i, s in enumerate(signs))) % order
            assert index
            projections.add(points[index][0])
        assert len(projections) == 2**(arity-1)
        counts.append(len(projections))
    print("Distinct final abscissas at arities 3 through 8:", counts)


def check_chords(p, a, c, zeros, nonzeros):
    points = [(x, y) for x, y in product(range(p), repeat=2) if (y*y-x**3-a*x-c) % p == 0]
    counts = Counter()
    for left, right in product(points, repeat=2):
        if left[0] != right[0]:
            counts[left[0], right[0], add(left, right, p, a)] += 1
    positive = negative = 0
    for target in points:
        w, t = target
        if not t:
            continue
        for u, v in product(range(p), repeat=2):
            if u == v:
                continue
            quartic = ((u-v)**2*w*w - 2*((u+v)*(u*v+a)+2*c)*w
                       + (u*v-a)**2 - 4*c*(u+v))
            is_zero = quartic % p == 0
            assert counts[u, v, target] == int(is_zero)
            positive += int(is_zero)
            negative += int(not is_zero)
    assert (positive, negative) == (zeros, nonzeros)
    print(f"Field {p}: {positive:,} quartic zeros and {negative:,} nonzeros checked.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_degree()
    for fixture in ((13, 4, 1, 288, 2520), (17, 2, 2, 288, 4608),
                    (19, 1, 0, 308, 5848), (37, 3, 1, 1124, 44164)):
        check_chords(*fixture)
