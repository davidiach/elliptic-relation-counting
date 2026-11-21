"""Out-of-sample tests of planted, full-moment decoders at balanced widths on toy curves."""

from collections import Counter
from fractions import Fraction
from math import log2, prod, sqrt
from pathlib import Path
from runpy import run_path

import numpy as np

source = run_path(str(Path(__file__).with_name("check_target_incidence.py")))
curve_cycle = source["curve_cycle"]


def study(p, a, widths, seed, n_plant=100000, n_train=25000, n_test=20000):
    curve = curve_cycle(p, a)
    order = len(curve)
    b, m = p.bit_length(), len(widths)
    assert m == 3
    modulus = 1 << b
    rho_product = Fraction(prod(widths), modulus**m)
    lam = order**(m-1)*rho_product
    assert 1/(2*log2(b)) < lam <= 1/log2(b)
    mu = Fraction((order-1)**m-(-1)**m, order)*rho_product
    target = 1
    x_t, y_t = curve[target]
    xs = np.array([-1]+[pt[0] for pt in curve[1:]])
    lifts = [[] for _ in range(modulus)]
    for i in range(1, order):
        lifts[curve[i][0]].append(i)
    columns = np.array([[0, 0, 0]] + [[(x-x_t) % p, (x*x-x_t*x_t) % p, (y+y_t) % p]
                                       for x, y in curve[1:]], dtype=np.int64)
    rng = np.random.default_rng(seed)
    projectors = rng.integers(0, p, (7, 2))

    def draw_seeds(planted, count):
        slopes = 2*rng.integers(0, modulus//2, size=(count, m)) + 1
        if not planted:
            return slopes, rng.integers(0, modulus, (count, m))
        left = rng.integers(1, order, size=(count, 2))
        last = (target-left.sum(axis=1)) % order
        invalid = last == 0
        while invalid.any():
            left[invalid] = rng.integers(1, order, size=(int(invalid.sum()), 2))
            last[invalid] = (target-left[invalid].sum(axis=1)) % order
            invalid = last == 0
        chosen = np.column_stack([left, last])
        addresses = rng.integers(np.zeros(m, dtype=int), widths, size=(count, m))
        return slopes, (addresses-slopes*xs[chosen]) % modulus

    inverses = {alpha: pow(alpha, -1, modulus) for alpha in range(1, modulus, 2)}

    def observe(slopes, offsets, want_truth):
        slots = []
        occupied = {x_t}
        regular = True
        for slope, beta, width in zip(slopes, offsets, widths):
            inv = inverses[slope]
            included = []
            for u in range(width):
                x = ((u-int(beta))*inv) % modulus
                if lifts[x]:
                    if x in occupied:
                        regular = False
                    occupied.add(x)
                    included.extend(lifts[x])
            slots.append(np.asarray(included, dtype=int))
        size = prod(map(len, slots))
        if want_truth and size:
            u, v, w = slots
            sums = (u[:, None, None] + v[None, :, None] + w[None, None, :]) % order
            truth = bool((sums == target).any())
        else:
            truth = False
        if not regular:
            return None, truth
        if size:
            left, middle, right = (columns[s] for s in slots)
            cross = np.cross(left[:, None, :], middle[None, :, :]).reshape(-1, 3)
            delta = (-cross @ right.T).ravel() % p
            if want_truth:
                assert bool((delta == 0).any()) == truth
            values = []
            power = delta.copy()
            for _ in range(6):
                values.append(int(power.sum() % p))
                power = power*delta % p
        else:
            values = [0]*6
        compressed = (np.array([1, *values]) @ projectors) % p
        return ((size, *values), (size,), (size, int(compressed[0])),
                (size, *map(int, compressed))), truth

    trained = [Counter() for _ in range(4)]
    fresh = [Counter() for _ in range(4)]
    for plant, count, tables in ((True, n_plant, trained), (False, n_train, fresh)):
        alphas, betas = draw_seeds(plant, count)
        for slopes, offsets in zip(alphas, betas):
            features, _ = observe(slopes, offsets, False)
            if features is not None:
                for table, key in zip(tables, features):
                    table[key] += 1
    # The sparse full-signature dictionary estimates support from positive examples.
    # For the two compressed features, minimize the empirical planted/unlabeled risk score.
    positive = [set(trained[0])]
    for positives, unlabelled in zip(trained[1:], fresh[1:]):
        positive.append({key for key, q_count in positives.items()
                         if 2*mu.numerator*q_count*n_train > mu.denominator*unlabelled[key]*n_plant})
    alphas, betas = draw_seeds(False, n_test)
    mistakes = [0]*4
    all_hits = regular = regular_hits = irregular_hits = 0
    for slopes, offsets in zip(alphas, betas):
        features, truth = observe(slopes, offsets, True)
        all_hits += truth
        if features is None:
            irregular_hits += truth
            for i in range(len(mistakes)):
                mistakes[i] += truth
            continue
        regular += 1
        regular_hits += truth
        for i, (table, key) in enumerate(zip(positive, features)):
            mistakes[i] += ((key in table) != truth)
    print(f"F_{p}, order {order}, widths {widths}, lambda {float(lam):.4f}: "
          f"{regular}/{n_test} regular, {all_hits} total positives ({regular_hits} regular).")
    print(f"  {n_plant:,} plants, {n_train:,} separate unlabeled training draws; "
          f"{len(positive[0]):,} distinct eligible planted full signatures.")
    for label, errors in zip(("predict absent", "full dictionary", "tuple count alone",
                              "one projection", "two projections"),
                             (all_hits, *mistakes)):
        rate = errors / n_test
        print(f"  {label}: {errors}/{n_test} errors ({100*rate:.2f}%, "
              f"test SE {100*sqrt(rate*(1-rate)/n_test):.2f} percentage points), "
              f"regular-only {(errors-irregular_hits)}/{regular}.")
    assert all_hits == regular_hits + irregular_hits and 0 < regular_hits < regular


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    study(101, 17, (4, 4, 4), 121617)
    study(1031, 9, (8, 16, 16), 121619)
