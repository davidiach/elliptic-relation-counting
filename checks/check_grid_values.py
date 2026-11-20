"""Check field-value rigidity on legal abscissa grids and the success-conditioned size bound."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import prod

from fixtures import curve_cycle, hash_families


def chord_in_last(u, v, a, c):
    return (u-v)**2, -2*((u+v)*(u*v+a)+2*c), (u*v-a)**2-4*c*(u+v)


def eliminant(xs, target_x, a, p):
    q, r, s = chord_in_last(xs[0], xs[1], a, 1)
    if len(xs) == 2:
        return (q*target_x**2 + r*target_x + s) % p
    t, u, v = chord_in_last(xs[2], target_x, a, 1)
    # The homogeneous resultant of the two chord quadratics in the intermediate abscissa.
    return ((q*v-s*t)**2 - (q*u-r*t)*(r*v-s*u)) % p


def leading_weights(nodes, p):
    return {x: pow(prod(x-y for y in nodes if y != x), -1, p) for x in nodes}


def value_grid():
    p, a = 101, 17
    curve = curve_cycle(p, a)
    size = len(curve)
    lifts = {}
    for label, (x, y) in enumerate(curve[1:], 1):
        lifts.setdefault(x, []).append(label)
    mask = sum(1 << bit for bit in (0, 1, 3, 4, 6))
    grids = [sorted(x for x in lifts if (11*x+7) % 128 < 32)]
    grids.extend(sorted(x for x in lifts if x & ~mask == code) for code in (0, 4))
    target_label = 1
    target_x = curve[target_label][0]
    for m in (2, 3):
        degree = 2**(m-1)
        assert min(map(len, grids[:m])) > degree
        nodes = [axis[:degree+1] for axis in grids[:m]]
        weights = [leading_weights(axis, p) for axis in nodes]
        leading = sum(eliminant(xs, target_x, a, p) *
                      prod(w[x] for w, x in zip(weights, xs)) for xs in product(*nodes)) % p
        assert leading == 1
        zeros = 0
        for xs in product(*grids[:m]):
            hit = any(sum(labels) % size == target_label
                      for labels in product(*(lifts[x] for x in xs)))
            assert (eliminant(xs, target_x, a, p) == 0) == hit
            zeros += int(hit)
        print(f"m={m}, field {p}, grid {list(map(len, grids[:m]))}, "
              f"top degree {m*degree}, top coefficient {leading}, {zeros} hitting tuples")


def size_bound():
    p, a, bits, width, m = 13, 4, 4, 4, 3
    modulus = 1 << bits
    curve = curve_cycle(p, a)
    size = len(curve)
    rho = Fraction(width, modulus)
    nu = size*rho
    intensity = size**(m-1)*rho**m
    anchors, masses = hash_families(curve, modulus, width)
    free = (1 << 0) + (1 << 2)
    patterns = [tuple(j for j, (x, y) in enumerate(curve[1:], 1) if x & ~free == code)
                for code in range(modulus) if code & free == 0]
    successes = Counter()
    bad = Counter()
    trials = 0
    for anchor, mass in zip(anchors, masses):
        for left, right in product(patterns, repeat=2):
            slots = (anchor, left, right)
            minimum = min(len({curve[j][0] for j in slot}) for slot in slots)
            trials += mass
            guarded = set()
            for i, j, k in product(*slots):
                if curve[i][0] == curve[j][0]:
                    continue
                mid = (i+j) % size
                if not mid or curve[mid][0] == curve[k][0]:
                    continue
                target = (mid+k) % size
                if target:
                    guarded.add(target)
            for target in guarded:
                successes[target] += mass
                for h in range(1, width+1):
                    if minimum <= h:
                        bad[target, h] += mass
    assert trials == 2048
    for target in range(1, size):
        kappa = Fraction(successes[target], trials)
        assert kappa
        for h in range(1, width+1):
            conditional = Fraction(bad[target, h], successes[target])
            assert conditional <= 2*m*h*intensity/(nu*kappa)
    print(f"Conditional size bound checked on {trials:,} weighted seeds, "
          f"{size-1} targets, {width} thresholds")


def dictionary_bound():
    p, modulus, width = 13, 16, 2
    curve = curve_cycle(p, 4)
    size = len(curve)
    anchor, mass = hash_families(curve, modulus, width)
    affine = list(zip(anchor, mass))
    fixed = [(tuple(j for j, (x, y) in enumerate(curve[1:], 1) if x & ~1 == code), 1)
             for code in range(0, modulus, 2)]
    rho = Fraction(width, modulus)

    def differences(law):
        result = Counter()
        for mask, weight in law:
            for x, y in product(mask, repeat=2):
                result[(y-x) % size] += weight
        return result

    for laws in ((affine, fixed), (fixed, fixed)):
        pair_mass = prod(sum(weight for _, weight in law) for law in laws)
        diffs = [differences(law) for law in laws]
        difference_energy = sum(diffs[0][d]*diffs[1][-d % size] for d in range(size))
        observed_energy = 0
        survives = Counter()
        for (a, wa), (b, wb) in product(*laws):
            counts = Counter((i+j) % size for i, j in product(a, b))
            observed_energy += wa*wb*sum(x*x for x in counts.values())
            most_likely = sorted(counts, key=lambda state: (-counts[state], state))
            for cap in (1, 2):
                selected = most_likely[:cap]
                for suffix, weight in fixed:
                    for target in {(s+k) % size for s, k in product(selected, suffix)}:
                        survives[cap, target] += wa*wb*weight
        assert observed_energy == difference_energy
        energy = Fraction(observed_energy, pair_mass)
        total_mass = pair_mass*sum(weight for _, weight in fixed)
        for cap in (1, 2):
            for target in range(size):
                prob = Fraction(survives[cap, target], total_mass)
                assert prob**2 <= rho**2*cap*energy
        kind = "affine/fixed" if laws[0] is affine else "fixed/fixed"
        print(f"{kind} block energy {energy}, adaptive dictionary bounds at {size} targets")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    value_grid()
    size_bound()
    dictionary_bound()
