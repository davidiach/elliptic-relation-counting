"""Check one-sided cuts, bounded product quotients and restricted quotient lists."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import isqrt, prod

from check_block_support import convolution, fixed_bit_law, group_tables
from check_seed_cuts import valid_chain
from fixtures import hash_families


def cut_check(p, a, c):
    curve, add, neg = group_tables(p, a, c)
    modulus = 1 << p.bit_length()
    masks, weights = hash_families(curve, modulus, 2)
    hidden = Counter(dict(zip(masks, weights)))
    fixed = fixed_bit_law(curve, p.bit_length(), (0,))
    target = next(i for i in range(1, len(curve)) if add[i][i])
    root_weight = sum(hidden.values()) * sum(fixed.values())**2
    singletons = wrong = 0
    harmonic = representative = 0
    for tails in product(fixed, repeat=2):
        tail_weight = prod(fixed[slot] for slot in tails)
        outside = defaultdict(list)
        for a_tail, b_tail in product(*tails):
            required = add[target][neg[add[a_tail][b_tail]]]
            outside[required % 2].append((a_tail, b_tail, required))
        for slot, weight in hidden.items():
            matches = [(r, x, y) for group in outside.values() for x, y, r in group if r in slot]
            if len(matches) != 1 or not valid_chain(matches[0], curve, add):
                continue
            actual, _, _ = matches[0]
            # Every inside member splices with every outside address of its color.
            feasible = sum(len(outside[point % 2]) for point in slot)
            mass = weight * tail_weight
            wrong += mass * (feasible - 1)
            singletons += mass
            harmonic += Fraction(mass, feasible)
            representative += Fraction(mass, len(outside[actual % 2]))
    bound = Fraction(4, modulus)  # Two boundary labels, each forced point has marginal 2/M.
    assert 0 < singletons < root_weight
    assert harmonic <= representative <= root_weight * bound
    assert Fraction(wrong, singletons) >= Fraction(singletons, root_weight)/bound - 1
    assert 3 * wrong > singletons
    print(f"F_{p}, {root_weight:,} weighted three-slot roots: {singletons:,} guarded singletons "
          f"at one fixed target. The one-bit cut retains the witness and averages "
          f"{float(Fraction(wrong, singletons)):.4f} wrong addresses; "
          f"reciprocal-count mass {float(harmonic/root_weight):.5f} <= representative mass "
          f"{float(representative/root_weight):.5f} <= {float(bound):.5f}.")


def short_prime_check(p, width, left, right, dx, dy):
    value = (left + dx) * (right + dy) % p
    cases = [(left + i)*(right + j) - value for i, j in product(range(width), repeat=2)]
    lo, hi = min(cases)//p, max(cases)//p
    height = hi - lo + 1
    assert height <= 2*width + 1
    true = [case % p == 0 for case in cases]
    assert 0 < sum(true) < len(cases)
    for small in range(2, height + 1):
        assert all(any((case - p*z) % small == 0 for z in range(lo, hi + 1)) for case in cases)
    for candidate in range(max(2, height + 1), 100000):
        if any(candidate % d == 0 for d in range(2, isqrt(candidate) + 1)) or candidate == p:
            continue
        inverse = pow(p, -1, candidate)
        accepts = [(case * inverse - lo) % candidate < height for case in cases]
        if accepts == true:
            print(f"{p.bit_length()}-bit p, two width-{width} intervals: "
                  f"{height} quotient values; validated {candidate.bit_length()}-bit prime "
                  f"{candidate} on all {len(cases)} inputs, {sum(true)} true.")
            return
    raise AssertionError("no good prime in fixture search range")


def quotient_check(p, a, c):
    curve, add, neg = group_tables(p, a, c)
    modulus = 1 << p.bit_length()
    masks, weights = hash_families(curve, modulus, 2)
    affine = Counter(dict(zip(masks, weights)))
    fixed = fixed_bit_law(curve, p.bit_length(), (0,))
    target = next(i for i in range(1, len(curve)) if add[i][i])
    final_weight = sum(fixed.values())
    caps = (1, 2, 4)
    covered = Counter()
    total_mass = 0
    for (f, wf), (g, wg) in product(affine.items(), fixed.items()):
        partial = convolution(Counter(f), Counter(g), add)
        collisions = sum(n*(n-1) for n in partial.values())
        prefix_weight = wf * wg
        for leaf, wj in fixed.items():
            choices = []
            for row in product(f, g, leaf):
                s = add[row[0]][row[1]]
                missing = add[target][neg[add[s][row[2]]]]
                z = curve[s][0]*curve[row[2]][0]//p if s else None
                choices.append((missing, row, z))
            local = defaultdict(Counter)
            for last, wl in fixed.items():
                hits = [(r, row, z) for r, row, z in choices if r in last]
                if len(hits) != 1:
                    continue
                r, row, z = hits[0]
                if z is None or not valid_chain((*row, r), curve, add):
                    continue
                local[row[-1]][z] += wl
                total_mass += prefix_weight * wj * wl
            for cap in caps:
                optimum = sum(sum(n for _, n in hist.most_common(cap)) for hist in local.values())
                bound = 0
                for j in leaf:
                    allowed = {z for z, _ in local[j].most_common(cap)}
                    v = curve[j][0]
                    possible = {r for r in range(1, len(curve)) if curve[r][0]*v//p in allowed}
                    size_bound = len(curve) if v == 0 else 2 * cap * ((p + v-1)//v)
                    assert len(possible) <= size_bound
                    candidates = sum(partial[s] for s in possible)
                    assert candidates <= len(possible) + collisions
                    # Forcing the last leaf multiplies each surviving tuple by at most rho.
                    bound += candidates
                assert optimum * modulus <= bound * final_weight * 2
                covered[cap] += prefix_weight * wj * optimum
    assert total_mass > 0 and 0 < covered[1] <= covered[2] <= covered[4] <= total_mass
    assert covered[1] < total_mass
    roots = sum(affine.values()) * final_weight**3
    rates = ", ".join(f"{float(Fraction(covered[h], total_mass)):.4f}" for h in caps)
    print(f"F_{p}, {roots:,} weighted four-slot roots: {total_mass:,} guarded singletons. "
          f"Optimal prefix- and leaf-dependent quotient lists at caps 1,2,4 cover {rates} "
          "conditionally; forced-tail and collision bounds checked for every local list.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    cut_check(13, 4, 1)
    short_prime_check(1000000007, 8, 380427101, 889291300, 3, 2)
    short_prime_check(2**61-1, 16, 910891223740200133, 1505050049293009348, 8, 5)
    quotient_check(13, 4, 1)
