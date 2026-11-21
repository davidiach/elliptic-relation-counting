"""Finite experiments with full field-valued determinant moments, using exact seed masses."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb, prod
from pathlib import Path
from runpy import run_path

import numpy as np

source = run_path(str(Path(__file__).with_name("check_target_incidence.py")))
curve_cycle, window = source["curve_cycle"], source["window"]
evaluation, determinant = source["evaluation"], source["determinant"]


def hash_families(curve, modulus, width):
    masses = Counter()
    for alpha, beta in product(range(1, modulus, 2), range(modulus)):
        mask = tuple(i for i, keep in enumerate(window(curve, modulus, width, alpha, beta))
                     if keep)
        masses[mask] += 1
    return list(masses), list(masses.values())


def aggregate(array, indicator):
    # Contract one point axis at a time, moving the new mask axis to the back.
    for _ in range(array.ndim):
        array = np.tensordot(array, indicator.T, axes=(0, 0))
    return array


def sign_orbit_product(size, moments, p):
    """Product and zero multiplicity on a sign orbit, or their pooled predictions."""
    mean_powers = [1, *(int(z) * pow(size, -1, p) % p for z in moments)]
    center = mean_powers[1]

    def central(j):
        return sum(comb(j, i)*mean_powers[i]*pow(-center, j-i, p) for i in range(j+1)) % p

    s1, fourth, sixth = central(2), central(4), central(6)
    s2 = (3*s1*s1 - fourth) * pow(2, -1, p) % p
    s3 = (sixth - 15*fourth*s1 + 30*s1**3) * pow(16, -1, p) % p
    e2 = (s1*s1 - s2) * pow(2, -1, p) % p
    e3 = (s1**3 - 3*s1*s2 + 2*s3) * pow(6, -1, p) % p
    quartic = [s1*s1 - 4*e2, 0, -2*s1, 0, 1]
    poly = [0]*9
    for i, a in enumerate(quartic):
        for j, b in enumerate(quartic):
            poly[i+j] += a*b
    poly[2] -= 64*e3
    value = sum(c*pow(center, j, p) for j, c in enumerate(poly)) % p
    multiplicity = 0
    while len(poly) > 1:
        quotient = [poly[-1]]
        for coefficient in reversed(poly[1:-1]):
            quotient.append((coefficient + center*quotient[-1]) % p)
        if (poly[0] + center*quotient[-1]) % p:
            break
        multiplicity += 1
        poly = quotient[::-1]
    return value, multiplicity


def exact_decoder(p, a, arity, width, target, orders):
    curve = curve_cycle(p, a)
    order = len(curve)
    masks, masses = hash_families(curve, 1 << p.bit_length(), width)
    indicator = np.zeros((len(masks), order - 1), dtype=np.int64)
    for i, mask in enumerate(masks):
        if mask:
            indicator[i, np.array(mask) - 1] = 1
    abscissas = [{curve[i][0] for i in mask} for mask in masks]
    rows = [evaluation(point, arity+1, p) for point in curve[1:]]
    fixed = evaluation(curve[-target], arity+1, p)
    delta = np.empty((order-1,)*arity, dtype=np.int64)
    hit = np.empty_like(delta)
    for index in np.ndindex(delta.shape):
        delta[index] = determinant([*(rows[i] for i in index), fixed], p)
        hit[index] = (sum(index) + arity) % order == target
    true_count = aggregate(hit, indicator)
    zero_count = aggregate((delta == 0).astype(np.int64), indicator)
    moments = []
    power = np.ones_like(delta)
    for j in range(1, max(orders)+1):
        power = (power * delta) % p
        moments.append(aggregate(power, indicator) % p)

    tables = {k: defaultdict(lambda: [0, 0]) for k in orders}
    flag_tables = {k: defaultdict(lambda: [0, 0]) for k in orders}
    regular_mass = hits = 0
    orbit_mass = orbit_hits = pool_mass = pooled_errors = 0
    for index in np.ndindex(true_count.shape):
        sets = [abscissas[i] for i in index]
        if (any(curve[target][0] in xs for xs in sets)
                or sum(map(len, sets)) != len(set().union(*sets))):
            continue
        assert true_count[index] == zero_count[index]
        label = int(true_count[index] > 0)
        weight = prod(masses[i] for i in index)
        regular_mass += weight
        hits += weight * label
        sizes = tuple(len(masks[i]) for i in index)
        values = tuple(int(moment[index]) for moment in moments)
        size = prod(sizes)
        if arity == 3 and size % p:
            predicted_product, predicted_hits = sign_orbit_product(size, values, p)
            if sizes == (2, 2, 2):
                true_product = prod(int(delta[tuple(u-1 for u in choice)])
                                    for choice in product(*(masks[i] for i in index))) % p
                assert predicted_product == true_product and (true_product == 0) == bool(label)
                assert predicted_hits == true_count[index]
                orbit_mass += weight
                orbit_hits += weight*label
            else:
                pool_mass += weight
                pooled_errors += weight*((predicted_product == 0) != label)
        for k in orders:
            tables[k][sizes + values[:k]][label] += weight
            flags = tuple(v == 0 for v in values[:k])
            flag_tables[k][sizes + flags][label] += weight

    print(f"F_{p}, m={arity}, width={width}, target={target}: regular {regular_mass:,} "
          f"of {sum(masses)**arity:,} seeds, {hits:,} hits.")
    outcomes = []
    for k in orders:
        values = sum(min(counts) for counts in tables[k].values())
        flags = sum(min(counts) for counts in flag_tables[k].values())
        assert values <= flags and (not outcomes or values <= outcomes[-1])
        outcomes.append(values)
        print(f"  k={k:2}: full mistakes {values:>9,}, flags {flags:>9,}, "
              f"full value cells {len(tables[k]):,}")
    assert outcomes[0] > outcomes[-1]
    if arity == 3:
        print(f"  Sign-orbit formula exact on {orbit_mass:,} seed triples ({orbit_hits:,} hits); "
              f"pooling multiple nonempty orbits makes {pooled_errors:,}/{pool_mass:,} errors.")
    return regular_mass, hits, outcomes


def cross_target_lookup():
    p, a, width, k = 61, 4, 2, 6
    curve = curve_cycle(p, a)
    order = len(curve)
    masks, masses = hash_families(curve, 64, width)
    points = np.array(curve[1:])
    labels = np.arange(1, order)
    incidence = np.array([[int(i in mask) for i in labels] for mask in masks])
    sizes = incidence.sum(axis=1)
    degrees = np.outer(sizes, sizes)
    seed_masses = np.outer(masses, masses)
    disjoint = incidence @ incidence.T == 0

    def features(t):
        x, y = (points - curve[-t]).T
        # Determinant of rows (1,x,y) for P,Q,-T.
        delta = (np.outer(x, y) - np.outer(y, x)) % p
        h = incidence @ ((labels[:, None]+labels) % order == t).astype(int) @ incidence.T
        target_free = incidence[:, t-1] == 0
        keep = disjoint & target_free[:, None] & target_free & (degrees > k)
        collisions = incidence @ (delta == 0).astype(int) @ incidence.T
        assert np.array_equal(h[keep], collisions[keep])
        sketches = [degrees[keep].copy() for _ in range(3)]
        power = np.ones_like(delta)
        unit = pow(curve[t][1], -1, p)
        for j in range(1, k+1):
            power = power * delta % p
            value = (incidence @ power @ incidence.T)[keep] % p
            # Raw at t, raw at -t, and the sign-invariant normalization by y(t)^j.
            variants = (value, ((-1)**j * value) % p, value * pow(unit, j, p) % p)
            for s, v in enumerate(variants):
                sketches[s] = sketches[s] * p + v
        return sketches, h[keep], seed_masses[keep], int((h*seed_masses).sum()), keep

    test, counts_test, weights_test, unconditioned_hit_mass, eligible = features(1)
    labels_test = counts_test > 0
    hits = int(weights_test[labels_test].sum())
    print(f"F_61 held-out target, regular with >6 tuples: {weights_test.sum():,} seeds, {hits:,} hits.")
    keys = [np.unique(test[v]) for v in (0, 2)]
    training = [np.zeros((len(ks), 2), dtype=np.int64) for ks in keys]

    def match_indices(ks, observations):
        indices = np.searchsorted(ks, observations)
        valid = (indices < len(ks)) & (ks[np.minimum(indices, len(ks)-1)] == observations)
        return indices[valid], valid

    for t in range(2, (order+1)//2):
        observed, counts, weights, _, _ = features(t)
        truth = counts > 0
        # Omit both signs of target 1; give the raw dictionary both signs of every training target.
        for variant, key_index in ((0, 0), (1, 0), (2, 1)):
            indices, valid = match_indices(keys[key_index], observed[variant])
            np.add.at(training[key_index], (indices, truth[valid].astype(int)), weights[valid])
    for j, label in enumerate(("raw", "normalized")):
        observed = test[(0, 2)[j]]
        positions = np.searchsorted(keys[j], observed)
        counts = training[j][positions]
        seen = counts.sum(axis=1) > 0
        predictions = counts[:, 1] > counts[:, 0]  # An unseen key and a tie both default to zero.
        errors = int(weights_test[predictions != labels_test].sum())
        own = np.zeros((len(keys[j]), 2), dtype=np.int64)
        np.add.at(own, (positions, labels_test.astype(int)), weights_test)
        optimum = int(own.min(axis=1).sum())
        print(f"  {label}: {len(keys[j]):,} test signatures, own-target optimum {optimum:,}, "
              f"held-out errors {errors:,}, seen seed mass {weights_test[seen].sum():,}.")
        assert optimum == 0
        if j == 0:
            total = int(weights_test.sum())
            frequency = own.sum(axis=1) / total
            for sample_size in (1000, 10000, 100000):
                miss = np.dot(own[:, 1], np.exp(sample_size * np.log1p(-frequency))) / total
                print(f"  {sample_size:,} same-target labeled draws: expected unseen-positive "
                      f"error {100*miss:.4f}% (conditional, exact occupancy formula).")
            planted_mass = np.zeros(len(keys[j]), dtype=np.int64)
            np.add.at(planted_mass, positions, counts_test * weights_test)
            planted_frequency = planted_mass / unconditioned_hit_mass
            for sample_size in (1000, 10000, 100000):
                miss = np.dot(own[:, 1], np.exp(sample_size*np.log1p(-planted_frequency))) / total
                print(f"  {sample_size:,} unconditioned planted draws: expected dictionary "
                      f"error on test regime {100*miss:.4f}%.")
            # Decode the base-p packing to apply fixed projections chosen independently of labels.
            remaining = observed.copy()
            digits = []
            for _ in range(k):
                digits.append(remaining % p)
                remaining //= p
            values = np.array(digits[::-1]).T
            rng = np.random.default_rng(20261120)
            projectors = rng.integers(0, p, size=(k+1, 10))
            projections = (values @ projectors[1:] + projectors[0]) % p
            squares = set((v*v) % p for v in range(1, p))
            for mode, sizes_to_try in (("field", (1, 2, 3)), ("square", (3, 6, 10))):
                for size in sizes_to_try:
                    codes = remaining.copy()  # Include the exact tuple count in each feature.
                    for col in projections[:, :size].T:
                        if mode == "field":
                            codes = codes*p + col
                        else:
                            codes = codes*2 + np.isin(col, list(squares))
                    classes, group = np.unique(codes, return_inverse=True)
                    n = np.zeros(len(classes), dtype=np.int64)
                    positives = n.copy()
                    planted_counts = n.copy()
                    np.add.at(n, group, weights_test)
                    np.add.at(positives, group, weights_test*labels_test)
                    np.add.at(planted_counts, group, weights_test*counts_test)
                    prediction = 2*planted_counts > n
                    risk = int(np.where(prediction, n-positives, positives).sum())
                    ideal = int(np.minimum(positives, n-positives).sum())
                    print(f"  {size} {mode} projections: {len(classes):,} cells, "
                          f"optimal {ideal:,} errors, planted-score decoder {risk:,}.")

    # This trial learns using only planted draws and their moment signatures.
    # Cache feature evaluations on the toy mask family to avoid recomputing small determinants.
    signatures = np.full(eligible.shape, -1, dtype=np.int64)
    signatures[eligible] = test[0]
    mask_ids = {mask: j for j, mask in enumerate(masks)}
    address_to_mask = np.empty((32, 64), dtype=np.int64)
    for alpha, beta in product(range(1, 64, 2), range(64)):
        mask = tuple(i for i, present in enumerate(window(curve, 64, width, alpha, beta))
                     if present)
        address_to_mask[alpha//2, beta] = mask_ids[mask]
    rng = np.random.default_rng(601376)
    sample_size = 100000
    left = rng.choice(labels[labels != 1], size=sample_size)
    right = (1-left) % order
    training_ids = []
    for chosen in (left, right):
        half = rng.integers(0, 32, size=sample_size)
        u = rng.integers(0, width, size=sample_size)
        beta = (u - (2*half+1)*points[chosen-1, 0]) % 64
        training_ids.append(address_to_mask[half, beta])
    learned = np.unique(signatures[tuple(training_ids)])
    learned = learned[learned >= 0]
    predicted = np.isin(test[0], learned)
    false_positives = int(weights_test[predicted & ~labels_test].sum())
    errors = int(weights_test[predicted != labels_test].sum())
    assert false_positives == 0
    print(f"  Actual 100,000-plant trial: {len(learned):,} stored positive signatures, "
          f"{errors:,} test errors ({false_positives} false positives).")


def check_random_weight_rank():
    # Conditional linear map from one slot's unobserved weights to three field outputs.
    p, a, target, arity, k = 13, 4, 1, 2, 2
    curve = curve_cycle(p, a)
    masks, _ = hash_families(curve, 16, 4)
    d = arity+1
    fixed = evaluation(curve[-target], d, p)
    for left, right in product(masks, repeat=2):
        x_left = {curve[i][0] for i in left}
        x_right = {curve[i][0] for i in right}
        if (len(left) != 2 or len(right) < k+1 or x_left & x_right
                or curve[target][0] in x_left | x_right):
            continue
        values = [[determinant([evaluation(curve[u], d, p), evaluation(curve[v], d, p), fixed], p)
                   for v in right] for u in left]
        if len(set(values[0])) >= k+1:
            break
    else:
        raise AssertionError("Missing small Vandermonde fixture")
    bad = 0
    for weights in product(range(p), repeat=len(left)):
        columns = [[sum(w * pow(row[j], s, p) for row, w in zip(values, weights)) % p
                    for s in range(k+1)] for j in range(len(right))]
        bad += len(source["fixture"]["basis"](columns, p)) < k+1
    assert bad * p <= (k+1)*(arity-1) * p**len(left)
    print(f"Shared-weight moment vector: {bad} rank-deficient maps among {p**len(left)}.")


def check_moment_ambiguity():
    from itertools import combinations
    p, k, size = 29, 2, 4
    assert comb(p, size) > p**k
    seen = {}
    for values in combinations(range(p), size):
        moments = tuple(sum(pow(v, j, p) for v in values) % p for j in range(1, k+1))
        if moments not in seen:
            seen[moments] = set(values)
            continue
        first, second = seen[moments], set(values)
        shift = next(iter(first - second))
        shifted = [{(v-shift) % p for v in group} for group in (first, second)]
        assert (0 in shifted[0]) != (0 in shifted[1])
        assert all((sum(pow(v, j, p) for v in shifted[0]) -
                    sum(pow(v, j, p) for v in shifted[1])) % p == 0 for j in range(k+1))
        print("Two distinct ordinary four-element spectra have identical full moments 0..2 "
              "and opposite membership at zero.")
        return
    raise AssertionError("Missing finite-alphabet collision")


def check_spectra():
    p, a, width, m = 13, 4, 2, 3
    curve = curve_cycle(p, a)
    order, d, modulus = len(curve), m+1, 16
    masks, masses = hash_families(curve, modulus, width)
    incidence = np.array([[int(i in mask) for i in range(1, order)] for mask in masks])
    sizes = list(map(len, masks))
    supports = [{curve[i][0] for i in mask} for mask in masks]
    rows = [evaluation(pt, d, p) for pt in curve[1:]]
    spectra = []
    for target in (1, 2):
        fixed = evaluation(curve[-target], d, p)
        deltas = np.empty((order-1,)*m, dtype=np.int64)
        for idx in np.ndindex(deltas.shape):
            deltas[idx] = determinant([*(rows[i] for i in idx), fixed], p)
        spectra.append(np.array([aggregate((deltas == v).astype(np.int64), incidence)
                                 for v in range(p)]))
    collision = cross = 0
    checked = False
    for idx in product(range(len(masks)), repeat=m):
        sets = [supports[i] for i in idx]
        union = set().union(*sets)
        if len(union) != sum(map(len, sets)) or curve[1][0] in union:
            continue
        weight = prod(masses[i] for i in idx)
        population = prod(sizes[i] for i in idx)
        bins = spectra[0][(slice(None), *idx)]
        assert bins.sum() == population
        collision += weight * int(bins @ bins - population)
        if curve[2][0] not in union:
            cross += weight * int(bins @ spectra[1][(slice(None), *idx)])
        values = np.flatnonzero(bins[1:] % p) + 1
        if not checked and 0 < population < p and len(values) >= 4 and bins[0]:
            s = len(values)
            seq = [sum(int(bins[v])*pow(int(v), j, p) for v in values) % p
                   for j in range(1, 2*s+1)]
            hankel = [seq[i:i+s] for i in range(s)]
            assert determinant(hankel, p) != 0
            checked = True
    assert checked
    nu = Fraction(order*width, modulus)
    lam = nu**m / order
    proper = d * order * lam**2 * (1 + 2/nu)**m
    total = sum(masses)**m
    assert Fraction(collision, total) <= proper + m*d*(2**m-1)*lam
    assert Fraction(cross, total) <= proper + m*d*2**m*lam
    print(f"Checked collision and cross-target bounds on all {total:,} seed triples, "
          "and a nonsingular full-order moment Hankel matrix.")


def check_planting():
    p, a, width, m, target = 13, 4, 2, 2, 1
    curve = curve_cycle(p, a)
    modulus, order = 16, len(curve)
    masks, masses = hash_families(curve, modulus, width)
    accepts = [[i for i, mask in enumerate(masks) if point in mask] for point in range(order)]
    fresh_mass = np.outer(masses, masses)
    planted = np.zeros_like(fresh_mass)
    for left in range(1, order):
        right = (target-left) % order
        if right:
            indices = np.ix_(accepts[left], accepts[right])
            planted[indices] += fresh_mass[indices]
    h = np.array([[sum((u+v) % order == target for u, v in product(left, right))
                   for right in masks] for left in masks])
    assert np.array_equal(planted, fresh_mass*h)
    # Derive each slot's conditional seed law directly with alpha uniform and uniform retained u.
    seeds = sum(masses)
    rho = Fraction(width, modulus)
    for i in range(1, order):
        x = curve[i][0]
        drawn = Counter()
        for alpha, u in product(range(1, modulus, 2), range(width)):
            beta = (u - alpha*x) % modulus
            mask = tuple(j for j, hit in enumerate(window(curve, modulus, width, alpha, beta)) if hit)
            drawn[mask] += 1
        assert sum(drawn.values()) == seeds*rho
        assert all(drawn[mask] == (mass if i in mask else 0) for mask, mass in zip(masks, masses))
    mu = Fraction(int(planted.sum()), seeds**m)
    assert mu == Fraction((order-1)**m-(-1)**m, order) * rho**m
    for truth in (h > 0, h % 2 == 1):
        discrepancy = int((fresh_mass * (h-truth)).sum())
        for prediction in (h == 0, h % 3 == 1, fresh_mass % 3 == 0, planted % 5 == 0):
            score = int(planted.sum() + fresh_mass[prediction].sum() - 2*planted[prediction].sum())
            risk = int(fresh_mass[prediction != truth].sum())
            assert abs(score-risk) <= discrepancy
    print("Planted seed law and unlabeled risk correction agree with exact enumeration.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    exact_decoder(13, 4, 2, 2, 1, (0, 1, 2, 3, 4, 5, 6))
    exact_decoder(23, 5, 2, 4, 1, (0, 1, 2, 3, 4, 5, 6))
    exact_decoder(13, 4, 3, 2, 1, (0, 1, 2, 3, 4, 5, 6))
    cross_target_lookup()
    check_random_weight_rank()
    check_moment_ambiguity()
    check_spectra()
    check_planting()
