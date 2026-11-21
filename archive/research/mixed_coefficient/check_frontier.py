"""Finite checks of the separated-approximation bounds in frontier.tex."""

from itertools import combinations
from fractions import Fraction
from math import pi, prod, sin
from pathlib import Path
from runpy import run_path

import numpy as np

# Reuse the project's audited curve arithmetic and exact hash seed masses.
fixture = run_path(str(Path(__file__).with_name("verify.py")))


def squared_tails(values):
    return np.r_[np.cumsum(values[::-1] ** 2)[::-1], 0.0]


def check_curve(prime, coefficient, widths):
    group = fixture["curve_cycle"](prime, coefficient)
    addition = fixture["table"](group, prime, coefficient)
    order = len(group)
    pairs = (order - 1) // 2
    modulus = 1 << prime.bit_length()
    intensity = prod(order * width / modulus for width in widths) / order
    scale = intensity / order
    families, masses = zip(*(fixture["families"](group, modulus, w) for w in widths))
    shape = tuple(map(len, families))
    seed_mass = np.ones(shape)
    slot_product = np.ones(shape)
    for axis, (masks, mass) in enumerate(zip(families, masses)):
        axis_shape = [1] * len(widths)
        axis_shape[axis] = -1
        seed_mass *= (np.array(mass) / sum(mass)).reshape(axis_shape)
        slot_product *= np.array(list(map(len, masks))).reshape(axis_shape)
    ordinary = np.empty((*shape, order), dtype=np.int64)
    guarded = np.empty_like(ordinary)
    for idx in np.ndindex(shape):
        slots = [family[i] for family, i in zip(families, idx)]
        ordinary[idx], guarded[idx] = fixture["counts"](slots, group, addition)

    assert np.isclose(seed_mass.sum(), 1.0)
    tests = 0
    ranks = np.arange(pairs + 1)
    offset = slot_product / order
    deficit_cap = 4 * sum(1 + pairs * w / modulus for w in widths)
    for target in range(1, order):
        h = ordinary[..., target]
        x = guarded[..., target]
        residual = h - offset
        diagonal = 2 * np.cos(2 * pi * target * np.arange(1, pairs + 1) / order)
        ideal = squared_tails(np.sort(abs(diagonal))[::-1])
        for bit in (h % 2, (x > 0).astype(int)):
            assert np.all((h - bit) ** 2 <= 2 * h * (h - 1) + ((h == 1) & (x == 0)))
        for length in range(1, len(shape)):
            for left in combinations(range(len(shape)), length):
                if 0 not in left:
                    continue
                axes = left + tuple(i for i in range(len(shape)) if i not in left)
                rows = prod(shape[i] for i in left)

                def error_by_rank(data):
                    matrix = (data * np.sqrt(seed_mass)).transpose(axes).reshape(rows, -1)
                    singular = np.linalg.svd(matrix, compute_uv=False)
                    tail = squared_tails(singular)
                    return tail[np.minimum(ranks, len(singular))]

                observed = error_by_rank(residual) / scale
                deficit = ideal[0] - observed[0]
                assert -1e-9 <= deficit <= deficit_cap + 1e-9
                assert np.all(observed <= ideal + 1e-9)
                assert np.all(observed >= ideal - deficit - 1e-9)
                for bit in (h % 2, (x > 0).astype(int)):
                    bit_error = np.sqrt(error_by_rank(bit))
                    distance = np.sqrt(np.sum(seed_mass * (residual - bit) ** 2))
                    assert np.all(abs(bit_error - np.sqrt(observed * scale)) <= distance + 1e-9)
                tests += len(ranks)
    print(f"Order {order}, widths {widths}: {tests} target, split and rank combinations passed.")


def check_integral():
    order = 65537
    pairs = (order - 1) // 2
    diagonal = 2 * np.cos(2 * pi * np.arange(1, pairs + 1) / order)
    tail = squared_tails(np.sort(abs(diagonal))[::-1]) / order
    fractions = np.arange(pairs + 1) / pairs
    limit = 1 - fractions - np.sin(pi * fractions) / pi
    assert np.max(abs(tail - limit)) < 5 / order
    half_energy = 0.26474189536615045
    assert abs(half_energy + sin(pi * half_energy) / pi - 0.5) < 1e-14
    print(f"Checked all {pairs + 1} cosine tails at order {order}.")


def check_smaller_width():
    curve = fixture["curve_cycle"](1031, 9)
    order, modulus = len(curve), 2048
    assert order == 1019
    widths = [1 << i for i in range(modulus.bit_length())]
    new = next(w for w in widths if order*w >= 40*modulus
               and order*w*w >= modulus*modulus)
    old = next(w for w in widths if order*w >= 40*modulus
               and order*w*w >= 4*modulus*modulus)
    assert (new, old) == (128, 256)
    # For target 1 in cyclic coordinates, the second point is determined.
    tuples = [(curve[i][0], curve[(1-i) % order][0]) for i in range(2, order)
              if curve[i][0] != curve[(1-i) % order][0]]
    first, second = np.asarray(tuples, dtype=np.int64).T

    def pair_weights(xs):
        delta = xs[:, None] - xs[None, :]
        return np.where(delta == 0, new*modulus, (delta % new != 0) * new**2)

    mean = len(tuples) * Fraction(new, modulus)**2
    moment = Fraction(int((pair_weights(first) * pair_weights(second)).sum()), modulus**4)
    population = Fraction(order * new, modulus)
    gamma = 1 - Fraction(6, order)
    cap = ((1 + 2/population)**2 - (2/population)**2 + order/population**2) / gamma
    cap += 2 / (population**2 * gamma**2)  # the group has one 2-torsion point
    assert moment / mean**2 <= cap < Fraction(5, 2)
    print(f"Order {order}, two-slot width {old} to {new}: exact moment coverage "
          f"lower bound {float(mean**2/moment):.6f}.")


def check_fixed_multipliers():
    curve = fixture["curve_cycle"](13, 4)
    order, modulus = len(curve), 16
    addition = fixture["table"](curve, 13, 4)
    for left_width, right_width in ((1, 1), (1, 2)):
        sides = []
        for alpha, width in ((3, left_width), (7, right_width)):
            sides.append([[i for i, pt in enumerate(curve) if pt is not None
                           and (alpha*pt[0] + beta) % modulus < width]
                          for beta in range(modulus)])
        counts = np.array([[fixture["counts"]((l, r), curve, addition)[1] for r in sides[1]]
                           for l in sides[0]])
        widths_product = left_width * right_width
        for target in range(1, order):
            a = counts[..., target]
            assert a.sum() >= (order - 6) * widths_product
            assert max(a.sum(axis=0)) <= 2 * widths_product
            assert max(a.sum(axis=1)) <= 2 * widths_product
            op = np.linalg.svd(a.astype(float), compute_uv=False)[0]
            assert np.sum(a*a) / op**2 >= (order - 6)/(4 * widths_product) - 1e-9
    print("Fixed multipliers: guarded row, column and stable-rank bounds passed.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_curve(13, 4, (1, 1, 1))
    check_curve(23, 5, (1, 2))
    check_integral()
    check_smaller_width()
    check_fixed_multipliers()
