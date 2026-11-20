"""Exact coverage of optimal short completion lists when one sampled slot is hidden."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import prod

from check_block_support import fixed_bit_law, group_tables
from fixtures import hash_families


def optimal_lists(law, xbit, key_count, caps):
    """Maximize coverage over all sublists of each possible complement set."""
    denominator = sum(law.values())
    x_masks = Counter()
    for points, weight in law.items():
        bits = 0
        for point in points:
            bits |= xbit[point]
        x_masks[bits] += weight
    hits = [sum(weight for mask, weight in x_masks.items() if mask & subset)
            for subset in range(1 << key_count)]
    optimum = {}
    for k in caps:
        best = [mass if subset.bit_count() <= k else 0 for subset, mass in enumerate(hits)]
        # Subset-max transform. After processing all bits this maximizes over
        # every sublist of size ≤ k, including when different candidates share a seed.
        for j in range(key_count):
            bit = 1 << j
            for subset in range(1 << key_count):
                if subset & bit:
                    best[subset] = max(best[subset], best[subset ^ bit])
        optimum[k] = best
    return hits, optimum, denominator


def valid_chain(choices, curve, add):
    partial = choices[0]
    for point in choices[1:]:
        if not partial or curve[partial][0] == curve[point][0]:
            return False
        partial = add[partial][point]
    return bool(partial)


def verify(p, a, c, slot_count, free):
    curve, add, neg = group_tables(p, a, c)
    xs = sorted({point[0] for point in curve[1:]})
    by_x = {x: 1 << i for i, x in enumerate(xs)}
    xbit = [0] + [by_x[point[0]] for point in curve[1:]]
    modulus, width = 1 << p.bit_length(), 1 << len(free)
    masks, weights = hash_families(curve, modulus, width)
    affine, fixed = Counter(dict(zip(masks, weights))), fixed_bit_law(curve, p.bit_length(), free)
    laws = [affine] + [fixed]*(slot_count-1)
    caps = (1, 2, 4)
    totals = []
    for hidden in range(slot_count):
        hidden_law = laws[hidden]
        hits, best, denominator = optimal_lists(hidden_law, xbit, len(xs), caps)
        assert all(hits[1 << i]*modulus == denominator*width for i in range(len(xs)))
        if hidden:
            fixed_mask = modulus-1-sum(1 << j for j in free)
            for subset in range(1 << len(xs)):
                patterns = {x & fixed_mask for i, x in enumerate(xs) if subset & (1 << i)}
                for k in caps:
                    assert best[k][subset]*modulus == min(k, len(patterns))*denominator*width
        outside = [law for j, law in enumerate(laws) if j != hidden]
        peak = [0]*len(caps)
        target_pairs = [0, 0]
        for roll in product(*outside):
            outside_weight = prod(law[mask] for law, mask in zip(outside, roll))
            ordinary, guarded = [0]*len(curve), [0]*len(curve)
            for choices in product(*roll):
                partial = 0
                for v in choices:
                    partial = add[partial][v]
                for target in range(1, len(curve)):
                    missing = add[target][neg[partial]]
                    if missing:
                        ordinary[target] |= xbit[missing]
                        assigned = (*choices[:hidden], missing, *choices[hidden:])
                        if valid_chain(assigned, curve, add):
                            guarded[target] |= xbit[missing]
            # Independent incidence calculation on all full masks: add the
            # actual signed tuples forwards and accumulate their target sets.
            direct = [[0]*len(curve), [0]*len(curve)]
            for mask, weight in hidden_law.items():
                targets, legal = set(), set()
                for assignment in product(*roll[:hidden], mask, *roll[hidden:]):
                    partial = 0
                    for v in assignment:
                        partial = add[partial][v]
                    if partial:
                        targets.add(partial)
                        if valid_chain(assignment, curve, add):
                            legal.add(partial)
                for counts, found in zip(direct, (targets, legal)):
                    for target in found:
                        counts[target] += weight
            for j, complements in enumerate((ordinary, guarded)):
                assert [hits[subset] for subset in complements] == direct[j]
                target_pairs[j] += outside_weight*sum(direct[j])
                for subset in complements[1:]:
                    for k_idx, k in enumerate(caps):
                        attained = best[k][subset]
                        assert attained*modulus <= k*width*denominator
                        assert attained <= hits[subset]
                        peak[k_idx] = max(peak[k_idx], attained)
        totals.append(target_pairs)
        kind = "affine" if hidden == 0 else "fixed"
        rates = ", ".join(str(Fraction(v, denominator)) for v in peak)
        print(f"  {kind} slot {hidden+1}: peak optimal K=1,2,4 coverage {rates}")
    assert all(row == totals[0] for row in totals)
    roots = prod(sum(law.values()) for law in laws)
    print(f"p={p}, m={slot_count}, w={width}, {roots:,} weighted roots: "
          f"{totals[0][0]:,} ordinary and {totals[0][1]:,} guarded nonidentity root-target hits.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    verify(13, 4, 1, 3, (0, 2))
    verify(13, 4, 1, 4, (1,))
    verify(19, -1, 0, 3, (0, 2))
