"""Exact checks of trimmed-fiber selectors on two abstract cyclic groups."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb

from verify import modular_basis


def unweighted_count(prime, blocks, target, prefix, pins):
    law = Counter({0: 1})
    for i, block in enumerate(blocks):
        allowed = [x for x in block if (i >= len(prefix) or prefix[i] == x)
                   and (i not in pins or pins[i] == x)]
        next_law = Counter()
        for g, count in law.items():
            for x in allowed:
                next_law[(g + x) % prime] += count
        law = next_law
    return law[target]


def signed_count(prime, blocks, excluded, target, prefix, order):
    total = 0
    for h in range(order + 1):
        for places in combinations(range(len(blocks)), h):
            for labels in product(*(excluded[i] for i in places)):
                count = unweighted_count(prime, blocks, target, prefix,
                                         dict(zip(places, labels)))
                total += (-1)**h * count
    return total


def check_trim(prime, half, excluded_half, order, exhaustive, paths=None):
    blocks = half * 3
    excluded = excluded_half * 3
    half_values = Counter(sum(x) % prime for x in product(*half))
    energy = Fraction(prime * sum(n*n for n in half_values.values()),
                      sum(half_values.values())**2) - 1
    assert energy <= Fraction(1, 7)
    direct = Counter()
    if exhaustive:
        for x in product(*blocks):
            h = sum(a in bad for a, bad in zip(x, excluded))
            score = sum((-1)**i * comb(h, i) for i in range(min(h, order)+1))
            direct[sum(x) % prime] += score
    root_masses = []
    for target in range(prime):
        prefix = ()
        mass = signed_count(prime, blocks, excluded, target, prefix, order)
        root_masses.append(mass)
        if exhaustive:
            assert mass == direct[target]
        assert mass > 0
        for block in blocks:
            children = [signed_count(prime, blocks, excluded, target, prefix+(a,), order)
                        for a in block]
            assert sum(children) == mass
            index = next(i for i, value in enumerate(children) if value > 0)
            mass = children[index]
            prefix += (block[index],)
        assert mass == 1 and sum(prefix) % prime == target
        assert all(a not in bad for a, bad in zip(prefix, excluded))
        if paths is not None:
            paths[target] = prefix
    return energy, min(root_masses)


def check_corank(paths):
    prime = 11
    half = (tuple(range(1, 10)),) * 2
    nonref = tuple(block[1:] for block in half)
    coords = [(i, a) for i, block in enumerate(half) for a in block[1:]]

    def incidence(x):
        return tuple(sum(x[j] == a for j in range(i, len(x), 2)) for i, a in coords)

    def make_row(src, dst):
        assert (sum(dst)-sum(src)) % prime == 1
        return tuple(a-b for a, b in zip(incidence(dst), incidence(src)))

    choices = list(product(*(nonref*2)))
    by_sum = {sum(x) % prime: x for x in choices}
    assert len(by_sum) == prime
    rows = [make_row(x, by_sum[(sum(x)+1) % prime]) for x in choices]
    assert len(modular_basis(rows, prime)) == len(coords)-2
    src = (1, 2, 2, 2)
    dst = next(y for y in product(nonref[0], (1,), *nonref)
               if (sum(y)-sum(src)) % prime == 1)
    rows.append(make_row(src, dst))
    assert all(sum(row) == 0 for row in rows)
    rank = len(modular_basis(rows, prime))
    assert rank == len(coords)-1 and 3*len(half) < prime
    augmented = modular_basis([row+(1,) for row in rows], prime)
    assert len(augmented) == rank
    origin = incidence(paths[0])
    for t, path in paths.items():
        delta = tuple(x-y for x, y in zip(incidence(path), origin))
        assert len(modular_basis(augmented+[delta+(t,)], prime)) == rank
    assert len(paths) == prime
    # Pad each genuine relation with identical entries in a third block.
    # Its eight nonreference labels remain completely free in the quotient.
    width = len(coords)+8
    padded = modular_basis([row+(0,)*8 for row in rows], prime)
    padded_rhs = modular_basis([row+(0,)*8+(1,) for row in rows], prime)
    fractions = []
    for block in range(3):
        classes = []
        for a in range(9):
            v = tuple(int(a > 0 and j == 8*block+a-1) for j in range(width))
            for cell in classes:
                delta = tuple(x-y for x, y in zip(v, cell[0]))
                if len(modular_basis(padded+[delta], prime)) == rank:
                    cell.append(v)
                    break
            else:
                classes.append([v])
        fractions.append(Fraction(max(map(len, classes)), 9))
    for target, path in paths.items():
        delta = tuple(x-y for x, y in zip(incidence(path), origin))
        assert len(modular_basis(padded_rhs+[delta+(0,)*8+(target,)], prime)) == rank
    assert fractions == [Fraction(8, 9), Fraction(8, 9), Fraction(1, 9)]
    retained = fractions[0] * fractions[1]
    print("C11 one omitted block, rank/dimension, retained/full product:",
          (rank, width, retained, retained*fractions[2]))
    return rank, len(coords), len(paths)


if __name__ == "__main__":
    half = (tuple(range(3)), tuple(range(4)))
    excluded = ((2,), (3,))
    print("C5 minimum mass by odd order:",
          [min(signed_count(5, half*3, excluded*3, t, (), j) for t in range(5))
           for j in (1, 3, 5)])
    print("C5 fifth-order paths for all targets, chi-square and min mass:",
          check_trim(5, half, excluded, 5, True))
    paths = {}
    print("C11 first-order paths for all targets, chi-square and min mass:",
          check_trim(11, (tuple(range(1, 10)),)*2, ((1,),)*2, 1, False, paths))
    print("C11 known-shift rank, dimension, certified target values:", check_corank(paths))
