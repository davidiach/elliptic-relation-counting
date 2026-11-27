"""Exact threshold-three specialization of the seed potential on the order-19 curve."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, product

from check_displacement_kernels import CLASSES, R, affine_block
from check_short_relations import curve13, rank, short_fixture


def profile(block):
    points = set(block)
    return (len(points), *(sum((x + s) % R in points for x in points) for s in CLASSES))


def scaled_energy(left, right):
    # Bilinear in profiles, so this also works for sums over seed completions.
    mass = left[0] * right[0]
    collisions = mass + 2 * sum(a*b for a, b in zip(left[1:], right[1:]))
    return R * 169 * collisions - 5832 * mass


def sum_profiles(records):
    return tuple(map(sum, zip(*(p for _, p in records))))


def prefix_choice(records):
    # Each seed consists of four unsigned four-bit digits (a0,b0,a1,b1).
    # Preserve distinct map choices with the same block when taking expectations.
    available = [records, records]
    for slot in range(2):
        other = available[1-slot]
        other_total = sum_profiles(other)
        for bit in reversed(range(16)):
            groups = [[], []]
            for key, p in available[slot]:
                groups[(key >> bit) & 1].append((key, p))
            def average(group):
                return Fraction(scaled_energy(sum_profiles(group), other_total),
                                len(group) * len(other))
            available[slot] = min((g for g in groups if g), key=average)
    return tuple(group[0][0] for group in available)


def main():
    cycle, add = curve13()
    maps = list(product(range(1, 13), range(13)))
    by_profile = defaultdict(set)
    blocks = {}
    records = []
    for (a, b), (c, d) in product(maps, repeat=2):
        key = (a << 12) | (b << 8) | (c << 4) | d
        block = affine_block(cycle, [(a, b), (c, d)])
        p = profile(block)
        records.append((key, p))
        blocks[key] = block
        by_profile[p].add(block)

    minima = []
    best = None
    for left, right in combinations_with_replacement(by_profile, 2):
        value = scaled_energy(left, right)
        if best is None or value < best:
            best = value
            minima.clear()
        if value == best:
            minima.append((left, right))

    coverage = Counter()
    for left, right in minima:
        if left == right:
            pairs = combinations_with_replacement(sorted(by_profile[left]), 2)
        else:
            pairs = product(by_profile[left], by_profile[right])
        for pair in pairs:
            f = short_fixture(R, pair, cycle, add)
            target = f['n'] - 3
            assert all(rank(f['by_t'][s], R, f['n']) == target for s in range(1, R))
            collisions = sum(k*k for k in f['d'].values())
            assert best == R*169*collisions - 5832*f['W']
            coverage[len(f['d'])] += 1
    assert sum(map(len, by_profile.values())) == 3600
    assert (len(by_profile), len(minima), best) == (700, 9, -45998)
    assert coverage == {R: 144}
    print('All 6,481,800 unordered parity-split block pairs: 144 global minimizers; '
          'all cover 19/19 values and attain the target rank at every nonzero shift.')

    seeds = prefix_choice(records)
    pair = tuple(blocks[key] for key in seeds)
    f = short_fixture(R, pair, cycle, add)
    ranks = [rank(f['by_t'][s], R, f['n']) for s in CLASSES]
    value = scaled_energy(*map(profile, pair))
    assert seeds == (0x6165, 0xc6cc)
    assert (value, len(f['d']), ranks, f['n']-3) == (-26471, 13, [5]*9, 5)
    assert scaled_energy(*map(profile, pair)) == R*169*sum(k*k for k in f['d'].values()) - 5832*f['W']
    print(f'Prefix minimum of conditional expectations: seeds {[hex(k) for k in seeds]}, '
          f'blocks {pair}, scaled objective {value}, support {len(f["d"])}/19, '
          f'ranks {ranks}, target {f["n"]-3}.')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions required')
    main()
