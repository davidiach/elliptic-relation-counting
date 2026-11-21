"""Check one-sided group fingerprints and cut inequalities on a full mixed seed law."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import prod

from check_strided_isolation import curve, finite_laws


def minimum_cross_pairs(witnesses, rectangles):
    small, extra = divmod(witnesses, rectangles)
    return rectangles * small * (small - 1) + 2 * extra * small


def reduce_row(row, basis):
    while row:
        pivot = row.bit_length() - 1
        if pivot not in basis:
            break
        row ^= basis[pivot]
    return row


def check_correction():
    # One corrupted value in an exact function from Z/53 to each output group.
    prime = 53
    for order, scale in ((2, 0), (prime, 7)):
        table = [scale * x % order for x in range(prime)]
        table[9] = (table[9] + 1) % order
        errors = sum(table[(x+y) % prime] != (table[x] + table[y]) % order
                     for x, y in product(range(prime), repeat=2))
        defect = Fraction(errors, prime * prime)
        assert 0 < 6 * defect < 1
        corrected = []
        for v in range(prime):
            counts = Counter((table[(v+y) % prime] - table[y]) % order for y in range(prime))
            mode, tally = counts.most_common(1)[0]
            assert Fraction(tally, prime) >= 1 - 2 * defect
            corrected.append(mode)
        assert all(corrected[(x+y) % prime] == (corrected[x] + corrected[y]) % order
                   for x, y in product(range(prime), repeat=2))
        assert corrected == [scale * x % order for x in range(prime)]
    print("  Modal increments corrected all 53 positions of each perturbed additive table, "
          "for output groups of orders 2 and 53.")


def main():
    p, arity, buckets = 13, 3, 2
    points, add = curve(p, 4, 1)
    n = len(points)
    neg = [row.index(0) for row in add]
    affine, fixed = finite_laws(p, points, (0,))
    laws = [affine, fixed, fixed]
    # Index coefficients of all binary additive checks depending on the local legal mask.
    parameters = {("target", t): t for t in range(n)}
    for i, law in enumerate(laws):
        for local in law:
            for point, _ in local:
                parameters[i, local, point] = len(parameters)
    true_basis = {}
    wrong_rows = Counter()
    for law in laws:
        incidence = Counter()
        for local, weight in law.items():
            for point, _ in local:
                incidence[point] += weight
        assert len(incidence) == n - 1
        assert set(incidence.values()) == {Fraction(sum(law.values()), 8)}
    mask = [t for t in range(1, n) if add[t][t]]
    rho = Fraction(2, 1 << p.bit_length())
    tag = [(x.bit_count() & 1) for x in range(n)]

    def group_sum(items):
        result = 0
        for point in items:
            result = add[result][point]
        return result

    def accepted_in_rectangle(prefix_sums, outside_sums, target):
        left = Counter(tag[v] for v in prefix_sums)
        right = Counter(tag[add[target][neg[v]]] for v in outside_sums)
        return [left[j] * right[j] for j in range(buckets)]

    weighted_seeds = singleton_events = 0
    wrong = Counter()
    outside_sizes = Counter()
    private_trial = Counter()
    forced_wrong = Counter()
    for a, b, c in product(*laws):
        seed = (a, b, c)
        weight = prod(law[slot] for law, slot in zip(laws, seed))
        slots = [tuple(point for point, _ in slot) for slot in seed]
        weighted_seeds += weight
        for k in range(1, arity + 1):
            outside_sizes[k] += weight * prod(map(len, slots[k:]))
        xs = [{points[point][0] for point in slot} for slot in slots]
        if not all(slots) or any(xs[i] & xs[j] for i in range(arity) for j in range(i)):
            continue
        words = defaultdict(list)
        for u in product(*slots):
            words[group_sum(u)].append(u)
        used = set.union(*xs)
        witnesses = {}
        for t in mask:
            if points[t][0] in used or len(words[t]) != 1:
                continue
            u, = words[t]
            partial = add[u[0]][u[1]]
            if partial and points[partial][0] != points[u[2]][0]:
                witnesses[t] = u
        if not witnesses:
            continue
        row_without_target = {
            u: sum(1 << parameters[i, seed[i], pt] for i, pt in enumerate(u))
            for u in product(*slots)
        }
        for t, witness in witnesses.items():
            term = 1 << parameters["target", t]
            row = reduce_row(row_without_target[witness] | term, true_basis)
            if row:
                true_basis[row.bit_length() - 1] = row
            for u, partial in row_without_target.items():
                if u != witness:
                    wrong_rows[partial | term] += weight
        singleton_events += weight * len(witnesses)
        for k in range(1, arity + 1):
            prefixes = [group_sum(u) for u in product(*slots[:k])]
            suffixes = [group_sum(v) for v in product(*slots[k:])]
            by_outside = Counter(u[k:] for u in witnesses.values())
            minimal = sum(minimum_cross_pairs(tally, buckets) for tally in by_outside.values())
            forced_wrong[k] += weight * minimal
            false = 0
            for target, u in witnesses.items():
                areas = accepted_in_rectangle(prefixes, suffixes, target)
                witness_tag = tag[group_sum(u[:k])]
                assert areas[witness_tag] > 0
                count = sum(areas)
                brute_force = sum(tag[group_sum(v[:k])] == tag[add[target][neg[group_sum(v[k:])]]]
                                  for v in product(*slots))
                assert count == brute_force
                false += count - 1
                private_trial[k] += Fraction(weight, buckets * areas[witness_tag])
            assert false >= minimal
            wrong[k] += weight * false
    total = weighted_seeds * len(mask)
    assert weighted_seeds == 8192 and singleton_events
    print(f"F_{p}, group order {n}: {weighted_seeds:,} weighted mixed seeds, "
          f"{len(mask)} eligible targets, {singleton_events:,} regular guarded singletons.")
    for k in range(1, arity + 1):
        eps = Fraction(wrong[k], singleton_events)
        mean_out = Fraction(outside_sizes[k], weighted_seeds)
        assert mean_out == ((n - 1) * rho) ** (arity - k)
        expected_singleton_targets = Fraction(singleton_events, weighted_seeds)
        assert buckets * (1 + eps) >= expected_singleton_targets / mean_out
        if k < arity:
            trial = private_trial[k] / total
            q = Fraction(singleton_events, total)
            assert q / (buckets * (1 + eps)) <= trial <= (n * rho) ** min(k, arity - k) / n
        print(f"  Prefix {k}: two-tag test kept every witness; {wrong[k]:,} wrong "
              f"address/target pairs, including at least {forced_wrong[k]:,} forced by rectangles.")
    assert forced_wrong[arity] > 0
    planted_mass = sum(prod(len(s) * law[s] for s, law in zip(seed, laws))
                       for seed in product(*laws))
    assert Fraction(planted_mass, weighted_seeds) == ((n - 1) * rho) ** arity
    q_mask = Fraction(singleton_events, total)
    assert Fraction(singleton_events, planted_mass) == len(mask) * q_mask / (
        (n - 1) * rho) ** arity
    forced_additive = sum(count for row, count in wrong_rows.items()
                          if reduce_row(row, true_basis) == 0)
    print(f"  All complete binary additive tests with local-mask tags together leave "
          f"{forced_additive:,} of {sum(wrong_rows.values()):,} wrong pairs "
          f"(witness-constraint rank {len(true_basis)}).")

    # Test all binary projection rows on two distinct encoded group points.
    width = (n - 1).bit_length()
    rows = range(1 << width)
    for left in range(n):
        for right in range(left):
            matches = sum(((row & (left ^ right)).bit_count() & 1) == 0 for row in rows)
            assert 2 * matches == len(rows)
    print(f"  Every distinct group-label pair collides on exactly half of {len(rows)} "
          "binary rows; independent rows give the stated per-address expectation.")
    check_correction()


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
