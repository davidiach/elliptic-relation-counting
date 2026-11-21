"""Exact tests of a zero-count readout on genuine pooled inputs; these are small finite curves."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from random import Random

from pooled_zero import pooled_certificate

from fixtures import curve_cycle, hash_families, evaluation, determinant

ORDERS = (8, 9, 10, 12)


def pooled_observation(sets, target, p):
    fixed = evaluation((target[0], -target[1] % p), 4, p)
    sums = [0]*12
    count = 0
    for points in product(*sets):
        value = determinant([*(evaluation(point, 4, p) for point in points), fixed], p)
        count += value == 0
        power = 1
        for j in range(12):
            power = power*value % p
            sums[j] += power
    return tuple(z % p for z in sums), count


def tally(stats, cache, sets, target, p, mass, true_hits=None):
    values, truth = pooled_observation(sets, target, p)
    if true_hits is not None:
        assert truth == true_hits
    if p == 13:
        # Each nonzero determinant contributes one to the raw twelfth moment.
        assert 0 <= truth <= 4 and (16 - values[11]) % p == truth
        stats["Fermat exact"] += mass
    if values not in cache:
        cache[values] = pooled_certificate(p, values, ORDERS)
    stats["total"] += mass
    stats["positive"] += mass*(truth > 0)
    for j, possible in cache[values].items():
        if possible is not None:
            assert truth in possible
        if possible is not None and len(possible) == 1:
            stats[j, "exact"] += mass
        if possible == frozenset({0}) or (possible and 0 not in possible):
            stats[j, "bit"] += mass
            stats[j, "positive"] += mass*(truth > 0)
        elif truth:
            stats[j, "miss"] += mass
    return cache[values]


def report(stats):
    print(f"  total {stats['total']:,}, positives {stats['positive']:,}")
    for j in ORDERS:
        print(f"  prefix {j}: certified counts {stats[j,'exact']:,}, bits {stats[j,'bit']:,}, "
              f"positive certificates {stats[j,'positive']:,}, "
              f"errors defaulting undecided to absent {stats[j,'miss']:,}.")


def check_p17():
    p, target = 17, (0, 6)
    fibers = {x: [(x, y) for y in range(p) if (y*y-x**3-2*x-2) % p == 0]
              for x in range(p)}
    xs = [x for x, pair in fibers.items() if x != target[0] and len(pair) == 2]
    stats, cache = Counter(), {}
    for a, b in combinations(xs, 2):
        for u, v in combinations([x for x in xs if x not in (a,b)], 2):
            tally(stats, cache, (fibers[a], fibers[b], fibers[u]+fibers[v]), target, p, 1)
    print("F_17 equally weighted abscissa patterns:")
    report(stats)
    assert stats["total"] == 420 and stats[12, "exact"] == 385


def check_actual_seeds():
    p, target = 13, 1
    curve = curve_cycle(p, 4)
    masks, masses = hash_families(curve, 16, 2)
    abscissas = [set(curve[i][0] for i in mask) for mask in masks]
    stats, cache, keys, reference = Counter(), {}, {}, Counter()
    for idx in product(range(len(masks)), repeat=3):
        if sorted(len(masks[i]) for i in idx) != [2,2,4]:
            continue
        union = set().union(*(abscissas[i] for i in idx))
        if len(union) != 4 or curve[target][0] in union:
            continue
        chosen = [masks[i] for i in idx]
        hits = sum(sum(row) % len(curve) == target for row in product(*chosen))
        sets = [[curve[i] for i in mask] for mask in chosen]
        mass = masses[idx[0]]*masses[idx[1]]*masses[idx[2]]
        results = tally(stats, cache, sets, curve[target], p, mass, hits)
        # Only unresolved bins need a trained decision. The lookup saves
        # repeating the same raw-moment computations on this finite fixture.
        possible = results[9]
        if possible is None or len(possible) > 1:
            keys[idx] = possible
            reference[possible, hits] += mass
    print("F_13 regular two-orbit inputs, exact fresh-seed multiplicities:")
    report(stats)
    assert stats["Fermat exact"] == stats["total"] == 168480
    print(f"  Raw twelfth moment alone: {stats['Fermat exact']:,} exact counts, zero misses.")
    assert stats[12, "exact"] > 0 and stats[12, "positive"] > 0
    check_planted_readout(curve, masks, keys, reference, stats)
    check_marginals(len(curve), target, masks, masses)


def check_marginals(order, target, masks, masses):
    relation_mass, marked_mass = Counter(), Counter()
    for i, j, k in product(range(len(masks)), repeat=3):
        a, b, c = masks[i], masks[j], masks[k]
        d = len(a)*len(b)*len(c)
        if not d:
            continue
        h = sum((target-u-v) % order in c for u in a for v in b)
        w = masses[i]*masses[j]*masses[k]
        for key in ((), (i,), (i,j)):
            relation_mass[key] += w*h
            marked_mass[key] += w*d
    for idx, mass in marked_mass.items():
        if not idx:
            continue
        d = len(masks[idx[0]]) if len(idx) == 1 else len(masks[idx[0]])*len(masks[idx[1]])
        h = sum(sum(part) % order == target for part in product(*(masks[i] for i in idx)))
        t = 3-len(idx)
        predicted = ((1 + (-1)**t*(Fraction(order*h,d)-1)/(order-1)**t)
                     / (1+Fraction(1,(order-1)**3)))
        observed = Fraction(relation_mass[idx]*marked_mass[()], mass*relation_mass[()])
        assert observed == predicted
    print(f"  Planted-to-matched marginal likelihood verified on {len(marked_mass)-1:,} "
          "one- and two-slot patterns.")


def check_planted_readout(curve, masks, keys, reference, stats):
    modulus, width, target = 16, 2, 1
    order = len(curve)
    indices = {mask: i for i, mask in enumerate(masks)}
    seed_masks = {}
    for alpha, beta in product(range(1, modulus, 2), range(modulus)):
        mask = tuple(i for i in range(1, order)
                     if (alpha*curve[i][0]+beta) % modulus < width)
        seed_masks[alpha, beta] = indices[mask]
    samples = 60000
    tables = []
    marked_counts = Counter()
    for plant, seed in ((False, 54321), (True, 83914)):
        rng = Random(seed)
        counts = Counter()
        for _ in range(samples):
            if plant:
                while True:
                    i, j = rng.randrange(1, order), rng.randrange(1, order)
                    k = (target-i-j) % order
                    if k:
                        break
                anchors = (i, j, k)
            idx = []
            for slot in range(3):
                alpha = rng.randrange(1, modulus, 2)
                beta = ((rng.randrange(width)-alpha*curve[anchors[slot]][0]) % modulus
                        if plant else rng.randrange(modulus))
                idx.append(seed_masks[alpha,beta])
            idx = tuple(idx)
            if idx in keys:
                key = keys[idx]
                counts[key] += 1
                if plant:
                    # Count completions of the four small-slot pairs on
                    # training plants only. All test decisions use moments.
                    left, mid, right = sorted((masks[i] for i in idx), key=len)
                    h = sum((target-a-b) % order in right for a in left for b in mid)
                    assert h >= 1
                    marked_counts[key,h] += 1
        tables.append(counts)
    fresh, planted = tables
    # mu = ((order-1)^3+1)/order * (width/modulus)^3.
    mu_num, mu_den = ((order-1)**3+1)*width**3, order*modulus**3
    positive_bins = {key for key in planted if (key is None or 0 in key)
                     and 2*mu_num*planted[key] > mu_den*fresh[key]}
    fp = sum(mass for (key, h), mass in reference.items() if key in positive_bins and not h)
    fn = sum(mass for (key, h), mass in reference.items()
             if h and (key is None or 0 in key) and key not in positive_bins)
    print(f"  Nine moments, {samples:,} independent plants and {samples:,} unlabeled fresh seeds: "
          f"{len(positive_bins)} bins predicted positive, exact errors {fp+fn:,} "
          f"({fp:,} false positives, {fn:,} false negatives) of {stats['total']:,}.")
    print(f"  Defaulting unresolved cases to absent makes {stats[9,'miss']:,} errors.")
    # Q(H=j,v) = j P(H=j,v)/mu. Scale all masses by a common
    # 12*mu_den to keep division by j=1,2,3,4 exact.
    def count_decision(key):
        possible = sorted(key) if key is not None else list(range(5))
        w = {j: (12//j)*mu_num*marked_counts[key,j] for j in range(1,5)}
        w[0] = 12*mu_den*fresh[key] - sum(w.values())
        return max(possible, key=w.get)
    errors = sum(mass for (key, h), mass in reference.items() if count_decision(key) != h)
    baseline = sum(mass for (key, h), mass in reference.items() if (min(key) if key else 0) != h)
    print(f"  Plant-calibrated integer count: {errors:,} exact-population mistakes; "
          f"smallest candidate default: {baseline:,}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_p17()
    check_actual_seeds()
