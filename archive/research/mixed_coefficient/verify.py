"""Check the mixed-coefficient note. Run with python -B and NumPy."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import prod, sqrt
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "slot_arity"))
from verify import curve_cycle, is_prime, point_add, points, window


def factorial_bound(order, torsion, modulus, widths):
    sizes = [Fraction(order * width, modulus) for width in widths]
    intensity = prod(sizes) / order
    bound = intensity**2 * prod(1 + 2 / size for size in sizes)
    bound += intensity * torsion * (2**len(widths) - 2) / order
    return intensity, bound


def families(curve, modulus, width):
    """Collapse seeds with identical masks, preserving their exact masses."""
    masks = Counter()
    for alpha in range(1, modulus, 2):
        for beta in range(modulus):
            mask = window(curve, modulus, width, alpha, beta)
            masks[tuple(i for i, present in enumerate(mask) if present)] += 1
    return list(masks), list(masks.values())


def check_two_abscissas(modulus):
    xs = np.arange(modulus)
    for exp in range(modulus.bit_length()):
        width = 1 << exp
        hits = np.array([((alpha * xs + beta) % modulus < width)
                         for alpha in range(1, modulus, 2)
                         for beta in range(modulus)], dtype=np.int64)
        observed = hits.T @ hits
        different = (xs[:, None] - xs) % width != 0
        expected = different * width**2 + np.eye(modulus, dtype=np.int64) * modulus * width
        assert np.array_equal(2 * observed, expected)


def table(curve, p, a):
    index = {point: i for i, point in enumerate(curve)}
    return [[index[point_add(left, right, p, a)] for right in curve] for left in curve]


def counts(selected, curve, addition):
    ordinary = np.zeros(len(curve), dtype=np.int64)
    guarded = np.zeros_like(ordinary)
    for first, *tail in product(*selected):
        total = first
        valid = True
        for point in tail:
            valid = valid and total != 0 and curve[total][0] != curve[point][0]
            total = addition[total][point]
        ordinary[total] += 1
        guarded[total] += valid and total != 0
    return ordinary, guarded


def exhaustive(p, a, curve, widths):
    order = len(curve)
    modulus = 1 << p.bit_length()
    addition = table(curve, p, a)
    torsion = sum(addition[j][j] == 0 for j in range(order))
    slots, masses = zip(*(families(curve, modulus, width) for width in widths))
    shape = tuple(map(len, slots))
    total_mass = prod(map(sum, masses))
    ordinary = np.zeros(shape + (order,), dtype=np.int64)
    guarded = np.zeros_like(ordinary)
    size_product = np.empty(shape, dtype=np.int64)
    weights = np.empty_like(size_product)
    for idx in np.ndindex(shape):
        selected = [slot[k] for slot, k in zip(slots, idx)]
        ordinary[idx], guarded[idx] = counts(selected, curve, addition)
        weights[idx] = prod(mass[k] for mass, k in zip(masses, idx))
        size_product[idx] = prod(map(len, selected))
    assert weights.sum() == total_mass
    intensity, bound = factorial_bound(order, torsion, modulus, widths)
    arity = len(widths)
    rho_product = Fraction(prod(widths), modulus**arity)
    mean = rho_product * (((order - 1)**arity - (-1)**arity) // order)
    probabilities = []
    for t in range(1, order):
        h, x = ordinary[..., t], guarded[..., t]
        assert Fraction(int((h * weights).sum()), total_mass) == mean
        if addition[t][t] == 0:
            continue
        factorial = Fraction(int((h * (h-1) * weights).sum()), total_mass)
        unique = Fraction(int(((h == 1) * (x == 1) * weights).sum()), total_mass)
        guarded_mean = Fraction(int((x * weights).sum()), total_mass)
        assert factorial <= bound and unique >= guarded_mean - factorial
        assert guarded_mean >= intensity * (1 - Fraction(3 * arity, order))
        probabilities.append(unique)
    assert probabilities and min(probabilities) > 0
    if torsion == 1:
        assert np.all(ordinary[..., 0] % 2 == 0)
    print(f"order {order}, widths {widths}, |G[2]|={torsion}, {total_mass:,} seed vectors, "
          f"{len(probabilities)} eligible targets, exact singleton probabilities "
          f"from {float(min(probabilities)):.6f} to {float(max(probabilities)):.6f}")
    return ordinary, size_product, slots, masses


def check_spectrum(curve, modulus, widths, outputs):
    ordinary, size_product, slots, masses = outputs
    r, n = len(curve), (len(curve) - 1) // 2
    axis = np.arange(1, n + 1)
    # The curve_cycle order indexes all points by multiples of the first point.
    W = 2 * np.cos(2 * np.pi * np.outer(axis, axis) / r)
    assert np.allclose(W @ W.T, r * np.eye(n) - 2, atol=1e-11)
    xs = np.array([curve[j][0] for j in axis])
    probabilities = [np.array(mass) / sum(mass) for mass in masses]
    weights = probabilities[0]
    for prob in probabilities[1:]:
        weights = np.multiply.outer(weights, prob)
    grams, centered = [], []
    for width, slot, prob in zip(widths, slots, probabilities):
        rho = width / modulus
        bits = np.array([[int(j) in mask for j in axis] for mask in slot], dtype=float)
        features = bits @ W.T
        A = (features.T * prob) @ features / (r * rho)
        incidence_pairs = (xs[:, None] - xs) % width == 0
        deficit = (2 - rho) / r + rho / r * (W @ incidence_pairs @ W.T)
        assert np.allclose(np.eye(n) - A, deficit, atol=1e-12)
        centered_features = features - prob @ features
        Ac = (centered_features.T * prob) @ centered_features / (r * rho)
        assert np.allclose(np.eye(n) - Ac, deficit + rho / r, atol=1e-12)
        grams.append(A)
        centered.append(Ac)
    intensity = float(factorial_bound(r, 1, modulus, widths)[0])
    scale = intensity / r
    deficit_bound = 4 * sum(1 + n * w / modulus for w in widths)
    seed_mean = np.sum(weights * size_product) / r
    seed_variance = np.sum(weights * (size_product / r - seed_mean)**2)
    assert np.isclose(seed_mean, intensity * (1 - 1/r)**len(widths))
    variance_cap = seed_mean**2 * (prod(1 + 2*modulus/((r-1)*w) for w in widths) - 1)
    assert seed_variance <= variance_cap + 1e-12
    splits = [subset for k in range(1, len(widths)) for subset in combinations(range(len(widths)), k)
              if 0 in subset]
    for t in range(1, r):
        residual = ordinary[..., t] - size_product / r
        full_interaction = residual.copy()
        for i, prob in enumerate(probabilities):
            mean = np.tensordot(full_interaction, prob, axes=(i, 0))
            full_interaction -= np.expand_dims(mean, i)
        D = 2 * np.cos(2 * np.pi * axis * t / r)
        for observed, factors in ((residual, grams), (full_interaction, centered)):
            energy = np.sum(weights * observed**2)
            predicted = scale * (D @ np.prod(factors, axis=0) @ D)
            assert abs(energy - predicted) < 1e-11
            assert scale * (r-2-deficit_bound) - 1e-12 <= energy <= scale * (r-2) + 1e-12
        for left in splits:
            axes = left + tuple(i for i in range(len(widths)) if i not in left)
            rows = prod(residual.shape[i] for i in left)
            matrix = (residual * np.sqrt(weights)).transpose(axes).reshape(rows, -1)
            sv = np.linalg.svd(matrix, compute_uv=False)
            assert sv[0]**2 <= 4 * scale + 1e-12
            for q in range(n + 1):
                assert sum(sv[q:]**2) + 1e-12 >= scale * (r-2-deficit_bound-4*q)
            # Give a predictor exact parity on two rows selected after seeing the truth.
            # Its output bit matrix still has real rank at most two.
            mass = weights.transpose(axes).reshape(rows, -1)
            parity = (ordinary[..., t] % 2).transpose(axes).reshape(rows, -1)
            best_rows = np.argsort((mass * parity).sum(axis=1))[-2:]
            predictor = np.zeros_like(parity)
            predictor[best_rows] = parity[best_rows]
            u = (mass * predictor).sum()
            diff = (mass * (predictor != parity)).sum() - (mass * parity).sum()
            floor = (1 - 2 * seed_mean) * u - 2 * (sqrt(seed_variance) + 2*sqrt(2*scale)) * sqrt(u)
            assert diff + 1e-12 >= floor
    print(f"order {r}, all {r-1} nonidentity targets: Fourier and centered energies, "
          f"singular-value and parity bounds for all {len(splits)} seed splits")


def check_restrictions(curve, p, a, modulus, seeds, widths):
    addition = table(curve, p, a)

    def mask(seed, width):
        return tuple(j for j, hit in enumerate(window(curve, modulus, width, *seed)) if hit)

    root = [mask(seed, width) for seed, width in zip(seeds, widths)]
    h, x = counts(root, curve, addition)
    witnessed = 0
    for target, count in enumerate(h):
        if count % 2 == 0:
            continue
        masks = root.copy()
        for i, ((alpha, beta), initial) in enumerate(zip(seeds, widths)):
            width = initial
            while width > 1:
                width //= 2
                low = mask((alpha, beta), width)
                high = mask((alpha, beta-width), width)
                assert not set(low).intersection(high) and set(masks[i]) == set(low).union(high)
                masks[i] = low
                low_h, _ = counts(masks, curve, addition)
                if low_h[target] % 2 == 0:
                    beta -= width
                    masks[i] = high
        assert all(len(slot) <= 2 for slot in masks)
        leaf_h, leaf_x = counts(masks, curve, addition)
        assert leaf_h[target] % 2 == 1
        assert h[target] != 1 or leaf_x[target] == x[target]
        witnessed += int(h[target] == x[target] == 1)
    return witnessed


def large_fixture():
    p, a, r, modulus = 65521, 44, 65419, 65536
    assert is_prime(p) and is_prime(r) and (4 * a**3 + 27) % p
    squares = {y*y % p for y in range(1, p)}
    rhs = ((x**3 + a*x + 1) % p for x in range(p))
    assert 1 + sum(2 if value in squares else int(value == 0) for value in rhs) == r
    widths = (32, 32, 32)
    intensity, factorial = factorial_bound(r, 1, modulus, widths)
    assert Fraction(1, 3) < intensity <= Fraction(2, 3)
    isolation = intensity * (1 - Fraction(3 * len(widths), r)) - factorial
    q = 100
    deficit = 4 * sum(1 + Fraction((r-1) * width, 2 * modulus) for width in widths)
    remaining = (r - 2 - deficit - 4*q) / r
    assert isolation > Fraction(19, 100) and remaining > Fraction(99, 100)
    print(f"certified order {r}, widths {widths}, lambda={float(intensity):.8f}, "
          f"guarded isolation >= {float(isolation):.8f}, "
          f"rank-{q} residual MSE >= {float(remaining):.8f} times lambda")
    widths = (16, 32, 32)
    intensity, factorial = factorial_bound(r, 1, modulus, widths)
    mean = intensity * Fraction(r-1, r)**len(widths)
    variance_cap = mean**2 * (prod(1 + Fraction(2*modulus, (r-1)*w) for w in widths) - 1)
    # Minimize the prediction inequality over all positive prediction rates.
    gain = (sqrt(variance_cap) + 2*sqrt(q*intensity/r))**2 / (1-2*mean)
    isolation = intensity * (1 - Fraction(3 * len(widths), r)) - factorial
    error_floor = float(isolation) - gain
    assert error_floor > 0.1
    print(f"same curve, widths {widths}: rank-{q} Boolean parity error >= {error_floor:.8f}")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    for modulus in (8, 16, 32):
        check_two_abscissas(modulus)
    for p, a, widths in ((13, 4, (1, 2, 2)), (13, 4, (1, 1, 1, 1)), (23, 5, (1, 1, 1))):
        curve = curve_cycle(p, a)
        outputs = exhaustive(p, a, curve, widths)
        check_spectrum(curve, 1 << p.bit_length(), widths, outputs)
    composite = [None] + points(7, -1, 0)
    exhaustive(7, -1, composite, (1, 2, 2))
    restricted = sum(check_restrictions(curve_cycle(13, 4), 13, 4, 16,
                                       ((3, b), (7, b+1), (1, 7-b)), (2, 2, 4))
                     for b in range(16))
    assert restricted > 0
    print(f"dyadic parity restriction recovered {restricted} isolated guarded target instances")
    large_fixture()
