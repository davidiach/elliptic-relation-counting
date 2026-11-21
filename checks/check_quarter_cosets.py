"""Check the two-quarter inequalities on complete small curve groups and seed laws."""
from collections import Counter
from itertools import product

from check_block_support import convolution, fixed_bit_law, group_tables
from fixtures import hash_families


def subgroups(add):
    # Every finite elliptic-curve subgroup has at most two generators.
    generated = []
    for generator in range(len(add)):
        subgroup = {0}
        value = generator
        while value not in subgroup:
            subgroup.add(value)
            value = add[value][generator]
        generated.append(subgroup)
    return {frozenset(add[a][b] for a in left for b in right)
            for left, right in product(generated, repeat=2)}


def quarter_law(first, second, add):
    law = Counter()
    for (left, wl), (right, wr) in product(first.items(), second.items()):
        state = convolution(Counter(left), Counter(right), add)
        law[tuple(sorted(state.items()))] += wl * wr
    return law


def project(counts, labels):
    out = Counter()
    for x, n in counts.items():
        out[labels[x]] += n
    return out


def excess(counts):
    return sum(n * (n - 1) for n in counts.values())


def check(p, a, c):
    curve, add, negatives = group_tables(p, a, c)
    masks, weights = hash_families(curve, 1 << p.bit_length(), 2)
    affine = Counter(dict(zip(masks, weights)))
    fixed = fixed_bit_law(curve, p.bit_length(), (0,))
    other = fixed_bit_law(curve, p.bit_length(), (1,))
    left = quarter_law(affine, fixed, add)
    right = quarter_law(other, fixed, add)
    groups = [g for g in subgroups(add) if len(g) > 1]
    strict_cases = 0
    checked = 0
    for subgroup in groups:
        labels = [min(add[x][h] for h in subgroup) for x in range(len(curve))]

        def states(law):
            ans = []
            for items, weight in law.items():
                sums = dict(items)
                projected = project(sums, labels)
                pairs = sum(nx * ny for x, nx in items for y, ny in items
                            if add[negatives[x]][y] in subgroup) - sum(sums.values())
                assert excess(projected) == pairs
                ans.append((sums, weight, pairs, sum(sums.values())))
            return ans

        l_states, r_states = states(left), states(right)
        sum_products = 0
        for (l, wl, jl, nl), (r, wr, jr, nr) in product(l_states, r_states):
            sums = convolution(l, r, add)
            max_mass = max(project(sums, labels).values(), default=0)
            diff = max(0, max_mass - nl - nr)
            assert diff * diff <= jl * jr
            strict_cases += int(diff > 0)

            support = project(dict.fromkeys(sums, 1), labels)
            for den in (1, 2, 4):
                covered = sum(t for t in support.values() if t >= 2 and den * t >= len(subgroup))
                assert len(subgroup) * covered <= 2 * den * excess(support)
            sum_products += wl * wr * jl * jr
            checked += 1
        assert sum_products == sum(w * j for _, w, j, _ in l_states) * sum(
            w * j for _, w, j, _ in r_states)
    assert strict_cases > 0
    print(f"F_{p}, group size {len(curve)}: {len(groups)} nontrivial subgroups, "
          f"{checked:,} seed-state/subgroup checks, {strict_cases:,} where the two quarter sizes "
          "alone did not bound a coset's tuple mass. Three density thresholds checked in each case.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check(13, 4, 1)
    check(19, -1, 0)
