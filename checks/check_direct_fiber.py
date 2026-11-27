"""Check direct signed group incidence, translate fibers and Fourier sampling costs."""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from random import Random

import numpy as np
from flint import arb, ctx

from check_address_parity import character, transform
from check_strided_isolation import finite_laws
from fixtures import curve_cycle, determinant
from probe_signed_moment import Probe


def unpaired(xs):
    return any(n % 2 for n in Counter(xs).values())


def energy(curve, mask):
    n = len(curve)
    counts = Counter()
    pairs = {}
    for left, right in product(mask, repeat=2):
        key = (left+right) % n
        counts[key] += 1
        pairs.setdefault(key, set()).add(tuple(sorted((curve[left][0], curve[right][0]))))
    bad = any(unpaired(first+second) for sets in pairs.values()
              for first, second in combinations_with_replacement(sets, 2))
    result = sum(v*v for v in counts.values())
    if not bad:
        assert result <= 12*len(mask)**2
    return result, bad


def check_translates():
    for p, linear in ((101, 17), (1031, 9)):
        curve = curve_cycle(p, linear)
        xs = {}
        for i, (x, _) in enumerate(curve[1:], 1):
            xs.setdefault(x, []).append(i)
        trials = maximum = 0
        for ds in combinations_with_replacement((0, 1, 3, 5, 9), 4):
            if not unpaired(ds):
                continue
            zeros = 0
            for shift in range(p):
                slots = [xs.get((shift+d) % p, ()) for d in ds]
                zeros += any(sum(row) % len(curve) == 0 for row in product(*slots))
            assert zeros <= 16
            trials += 1
            maximum = max(maximum, zeros)
        print(f"p={p}: {trials} unpaired offset quadruples, all {p} shifts each; "
              f"largest zero count {maximum} (bound 16).")


def check_masks():
    p = 101
    curve = curve_cycle(p, 17)
    affine, fixed = finite_laws(p, curve, (0, 2))
    for name, law in (("affine", affine), ("fixed-bit", fixed)):
        total = sum(law.values())
        bad = 0
        largest = Fraction()
        cache = {}
        for row, weight in law.items():
            mask = tuple(i for i, _ in row)
            if mask not in cache:
                cache[mask] = energy(curve, mask)
            quad, exceptional = cache[mask]
            if exceptional:
                bad += weight
            elif mask:
                largest = max(largest, Fraction(quad, len(mask)**2))
        assert largest <= 12
        print(f"p={p}: {name} width-four law: {bad}/{total} masks have an unpaired "
              f"four-point relation; maximum good energy/size^2 {largest}.")


def first_character_phase_order(rows, targets, n):
    # Each column gives the integer coefficients of a sum of n-th roots of unity.
    coefficients = np.zeros((n, rows.shape[1]), dtype=np.int64)
    coefficients[targets] = rows
    reflected = coefficients[(-np.arange(n)) % n]
    real = np.all(coefficients == reflected, axis=0)
    imaginary = np.all(coefficients == -reflected, axis=0)
    assert np.all(real ^ imaginary), "expected each fixture column on one axis"

    # Conjugation fixes a real column and negates an imaginary column. Certify
    # the remaining sign, then order the four rays from -pi/2 through +pi.
    with ctx.workprec(96):
        sines, cosines = zip(*[(arb(2*g)/n).sin_cos_pi() for g in targets])
        rays = []
        for is_real, column in zip(real, rows.T.tolist()):
            coordinates = cosines if is_real else sines
            value = sum(c * x for c, x in zip(column, coordinates) if c)
            assert value > 0 or value < 0, "phase sign is unresolved"
            if is_real:
                rays.append(1 if value > 0 else 3)
            else:
                rays.append(2 if value > 0 else 0)
    # The numeric challenge is the secondary key for every exact phase tie.
    return sorted(range(len(rays)), key=lambda challenge: (rays[challenge], challenge))


def check_fourier():
    probe = Probe(1031, 9, (0, 2))
    rng = Random(13041)
    target = 1
    for _ in range(1000):
        slots = probe.planted(target, rng)
        addr, legal = probe.singleton(slots, target)
        if addr is None or min(map(len, slots)) < 4:
            continue
        abscissas = [{probe.curve[i][0] for i, _ in slot} for slot in slots]
        if (sum(map(len, abscissas)) != len(set.union(*abscissas))
                or any(probe.curve[target][0] in group for group in abscissas)):
            continue
        if any(energy(probe.curve, [i for i, _ in slot])[1] for slot in slots):
            continue
        break
    else:
        raise AssertionError("No regular singleton fixture found.")
    n, width, grid = probe.n, probe.size, len(legal)
    address_rows = {}
    for choices in product(*slots):
        labels, addresses = zip(*choices)
        endpoint = sum(labels) % n
        address = (addresses[0]*probe.block+addresses[1])*probe.block+addresses[2]
        address_rows.setdefault(endpoint, []).append(address)
        delta = determinant([probe.vectors[i] for i in labels]+[probe.vectors[-target]], probe.p)
        assert (1-pow(delta, probe.p-1, probe.p)) == int(endpoint == target)
    targets = sorted(address_rows)
    rows = []
    for t in targets:
        indicator = [0]*width
        for v in address_rows[t]:
            indicator[v] = 1
        rows.append(transform(indicator))
    rows = np.array(rows, dtype=np.int64)
    multiplicities = [len(address_rows[t]) for t in targets]
    assert np.array_equal(rows @ rows.T, width*np.diag(multiplicities))
    assert np.array_equal(rows[targets.index(target)], [character(addr, r) for r in range(width)])

    def phasor(labels):
        indices = np.multiply.outer(labels, np.arange(n)) % n
        return np.exp((2j*np.pi/n)*indices)

    factors = []
    for slot in slots:
        labels, addresses = zip(*slot)
        signs = np.array([[character(a, r) for a in addresses] for r in range(probe.block)])
        values = signs @ phasor(labels)
        fourth = (abs(values)**4).mean(axis=0)
        assert np.allclose((abs(values)**2).mean(axis=0), len(slot))
        assert fourth.mean() <= 12*len(slot)**2 + 1e-7
        factors.append(values)
    all_terms = (factors[0][:, None, None, :]*factors[1][None, :, None, :]
                 *factors[2][None, None, :, :]).reshape(width, n)
    phased = all_terms*phasor([-target])[0]
    truth = rows[targets.index(target)]
    assert np.allclose(phased.mean(axis=1), truth)
    uniform_mse = (abs(phased-truth[:, None])**2).mean()
    optimal_mse = ((abs(phased).mean(axis=1))**2).mean()-1
    assert abs(uniform_mse-(grid-1)) < 1e-8
    assert grid/12**len(slots)-1 <= optimal_mse < uniform_mse
    print(f"p={probe.p}: {grid} legal tuples, {len(targets)} represented targets, "
          f"{width} challenges; exact signed-row Gram matrix checked.")
    print(f"  One-sample complex mean squared error: uniform {uniform_mse:.6f}, "
          f"optimal challenge-dependent importance proposal {optimal_mse:.6f}.")
    order = first_character_phase_order(rows, targets, n)
    singletons = rows[np.array(multiplicities) == 1][:, order]
    prefixes = np.cumsum(singletons, axis=1)
    prefixes = np.column_stack((np.zeros(len(singletons), dtype=int), prefixes))
    best = 2*(prefixes.max(axis=1)-prefixes.min(axis=1))/width
    best[singletons.sum(axis=1) != 0] = 1  # The constant address character.
    for power in range(width.bit_length()):
        length = 2**power
        blocks = singletons.reshape(len(singletons), -1, length).sum(axis=2)/width
        assert np.all((blocks**2).sum(axis=0) <= length/width + 1e-12)
    assert np.all(best <= 1)
    assert (best**2).sum() <= 16*width.bit_length()**2+2
    print(f"  {len(singletons)} singleton targets: sum of best squared sign "
          f"correlations using one interval in the first character's phase order "
          f"{(best**2).sum():.6f}; all dyadic block bounds verified.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_translates()
    check_masks()
    check_fourier()
