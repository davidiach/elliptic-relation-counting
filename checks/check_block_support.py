"""Check block-sum support and collision formulas under exact small mixed seed laws."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import prod

from fixtures import hash_families


def group_tables(p, a, c):
    """Index the full curve group, including curves with several points of order two."""
    assert (4*a**3 + 27*c*c) % p != 0
    points = [None] + [(x, y) for x, y in product(range(p), repeat=2)
                       if (y*y-x**3-a*x-c) % p == 0]
    indices = {point: i for i, point in enumerate(points)}

    def add(left, right):
        if left is None or right is None:
            return right if left is None else left
        x, y = left
        u, v = right
        if x == u and (y+v) % p == 0:
            return None
        num, denom = (3*x*x+a, 2*y) if left == right else (v-y, u-x)
        slope = num*pow(denom, -1, p) % p
        nx = (slope*slope-x-u) % p
        return nx, (slope*(x-nx)-y) % p

    table = [[indices[add(left, right)] for right in points] for left in points]
    negatives = [row.index(0) for row in table]
    return points, table, negatives


def fixed_bit_law(curve, bits, positions):
    free_mask = sum(1 << j for j in positions)
    fixed_mask = (1 << bits) - 1 - free_mask
    return Counter(tuple(i for i, point in enumerate(curve) if i and point[0] & fixed_mask == pattern)
                   for pattern in range(1 << bits) if not pattern & free_mask)


def convolution(left, right, add):
    result = Counter()
    for a, n in left.items():
        for b, m in right.items():
            result[add[a][b]] += n*m
    return result


def difference_mass(law, add, negatives):
    result = Counter()
    for mask, weight in law.items():
        for a, b in product(mask, repeat=2):
            result[add[negatives[a]][b]] += weight
    return result


def verify(p, a, c, slots, positions, cap):
    curve, add, negatives = group_tables(p, a, c)
    order = len(curve)
    modulus = 1 << p.bit_length()
    width = 1 << len(positions)
    runs = sum(j-1 not in positions for j in positions)
    masks, weights = hash_families(curve, modulus, width)
    affine = Counter(dict(zip(masks, weights)))
    fixed = fixed_bit_law(curve, p.bit_length(), positions)
    laws = [affine] + [fixed]*(slots-1)
    assert Fraction(1, 3) <= Fraction(order-1, modulus) <= 2
    diffs = [difference_mass(law, add, negatives) for law in laws]
    for mass, law in zip(diffs, laws):
        denominator = sum(law.values())
        assert Fraction(mass[0], denominator) == Fraction(width*(order-1), modulus)
        assert max(Fraction(mass[d], denominator) for d in range(1, order)) <= Fraction(6*(1 << runs)*width**2, modulus)

    @lru_cache(None)
    def sums(parts):
        if not parts:
            return Counter({0: 1})
        return convolution(sums(parts[:-1]), Counter(parts[-1]), add)

    subsets = [tuple(s) for r in range(1, slots+1) for s in combinations(range(slots), r)]
    all_slots = tuple(range(slots))
    # For two unlabeled blocks, specifying the block that contains slot zero suffices.
    partitions = [(u, tuple(i for i in all_slots if i not in u)) for u in subsets if 0 in u]
    masses = {u: [0, 0, 0] for u in subsets}
    total_weight = prod(sum(law.values()) for law in laws)
    eligible_weight = eligible_hits = 0
    for root in product(*laws):
        weight = prod(laws[i][mask] for i, mask in enumerate(root))
        sizes = {(): 1}
        for subset in subsets:
            counts = sums(tuple(root[i] for i in subset))
            sizes[subset] = len(counts)
            acc = masses[subset]
            acc[0] += weight*len(counts)
            acc[1] += weight*sum(counts.values())
            acc[2] += weight*sum(v*v for v in counts.values())
        eligible = any(max(sizes[left], sizes[right]) <= cap for left, right in partitions)
        if eligible:
            assert sizes[all_slots] <= cap**2
            eligible_weight += weight
            eligible_hits += weight*len(set(sums(root)) - {0})
    for subset, (support, tuples, energy) in masses.items():
        diff_convolution = Counter({0: 1})
        law_weight = 1
        for i in subset:
            diff_convolution = convolution(diff_convolution, diffs[i], add)
            law_weight *= sum(laws[i].values())
        assert Fraction(energy, total_weight) == Fraction(diff_convolution[0], law_weight)
        assert support*energy >= tuples**2
        r = len(subset)
        assert Fraction(tuples, total_weight) == Fraction(width*(order-1), modulus)**r
        bound = Fraction(2*r*(1 << runs)*(2*width)**(2*r), modulus)
        assert Fraction(energy-tuples, total_weight) <= bound
    assert 0 < eligible_hits <= eligible_weight*cap**2
    assert 0 < eligible_weight < total_weight
    print(f"p={p}, m={slots}, w={width}: {total_weight:,} weighted roots; "
          f"{len(subsets)} block collision identities; cap {cap}: "
          f"{eligible_weight:,} admissible roots, {eligible_hits:,} nonidentity root-target pairs.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    verify(13, 4, 1, 3, (0, 2), 4)
    verify(13, 4, 1, 4, (1,), 4)
    verify(19, -1, 0, 3, (0, 2), 4)
