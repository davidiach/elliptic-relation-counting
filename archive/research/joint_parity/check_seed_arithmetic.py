"""Finite checks of seed_arithmetic.tex using exact counts on small curves."""

from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np

fixture = run_path(str(Path(__file__).resolve().parents[1] / "mixed_coefficient" / "verify.py"))


def check_progression_overlaps():
    trials = 0
    for modulus in (16, 32, 64, 128, 256, 512):
        for width in (2, 4, 8, 16):
            if 2 * width**2 >= modulus:
                continue
            for slope in range(3, modulus - 1, 2):
                intersections = Counter((x-slope*y) % modulus
                                        for x, y in product(range(width), repeat=2))
                assert max(intersections.values()) <= width // 2
                if slope == 1 + modulus // 2:
                    assert max(intersections.values()) == width // 2
                trials += modulus
    print(f"Checked {trials:,} combinations of direction, displacement and width.")


def check_mask_metric():
    curve = fixture["curve_cycle"](23, 5)
    order, modulus = len(curve), 32
    xs = np.array([-1] + [pt[0] for pt in curve[1:]])
    affine = np.arange(order) != 0
    shifts = np.arange(modulus)[:, None]
    residues = (xs - shifts) % modulus
    shapes = ((0, 1), (0, 3), (4, 5), (1, 3))
    width = 2
    masks = [np.isin(residues, shape) & affine for shape in shapes]
    families, masses = fixture["families"](curve, modulus, 2)
    # Every representative below is weighted by its exact fresh-seed multiplicity.
    dist = np.array(masses)
    denomin = modulus * sum(masses)
    slot = np.array([[int(j in mask) for j in range(order)] for mask in families])
    other_width = Fraction(2, modulus)
    tests = 0
    for arity in (2, 3):
        if arity == 2:
            others, weights, total = slot, dist, denomin
        else:
            pairs = list(product(range(len(slot)), repeat=2))
            others = np.array([sum(slot[a, j] * np.roll(slot[b], j) for j in range(order))
                               for a, b in pairs])
            weights = np.array([dist[a]*dist[b] for a, b in pairs])
            total = denomin * sum(masses)
        for target in range(1, order):
            kernel = others[:, (target - np.arange(order)) % order].T
            values = [mask.astype(int) @ kernel for mask in masks]
            root_sum = int((values[0] * weights).sum())
            for i, j in product(range(len(shapes)), repeat=2):
                d = len(set(shapes[i]) ^ set(shapes[j]))
                h = (masks[i] ^ masks[j]).astype(int) @ kernel
                assert np.array_equal(h % 2, (values[i] ^ values[j]) % 2)
                mean = Fraction(int((h * weights).sum()), total)
                factorial = Fraction(int((h*(h-1)*weights).sum()), total)
                assert mean == Fraction(root_sum * d, total * width)
                lam = Fraction(order**(arity-1)*d, modulus) * other_width**(arity-1)
                cap = (Fraction(2*modulus, order) * (1+2/(order*other_width))**(arity-1)
                       * lam**2 + (8+2*(2**(arity-1)-2)) * lam*d/order)
                assert factorial <= cap
                tests += 1
    print(f"Exact arbitrary-shape means, parity distances and factorial bounds: {tests} cases.")


def check_rounding_and_runs():
    curve = fixture["curve_cycle"](13, 4)
    order, modulus, width = len(curve), 16, 2
    seeds = list(product(range(1, modulus, 2), range(modulus)))
    xs = np.array([-1] + [pt[0] for pt in curve[1:]])
    affine = np.arange(order) != 0
    alpha, beta = np.array(seeds).T
    raw = (alpha[:, None]*xs + beta[:, None]) % modulus
    masks = (raw < width) & affine
    rounded = ((alpha[:, None]*xs + (beta - beta % 2)[:, None]) % modulus < width) & affine
    assert np.array_equal(masks.sum(axis=0), rounded.sum(axis=0))
    assert np.all(np.count_nonzero(masks ^ rounded, axis=0) <= (modulus//2)*2)
    lam = Fraction(order * width**2, modulus**2)
    seed_pairs = len(seeds)**2
    errors = []
    run_checks = 0
    for target in range(1, order):
        partner = (target - np.arange(order)) % order
        h = masks.astype(int) @ masks[:, partner].T
        guard = affine & (partner != 0) & (xs != xs[partner])
        exact_guarded = (masks * guard).astype(int) @ masks[:, partner].T
        new = (rounded * guard).astype(int) @ rounded[:, partner].T
        wrong = np.count_nonzero((new != 0) != (h % 2))
        term = (h*(h-1)//2 + h-exact_guarded).sum()
        assert Fraction(wrong, seed_pairs) <= 4*lam + Fraction(int(term), seed_pairs)
        errors.append(wrong / seed_pairs)
        for ai in range(1, modulus, 2):
            counts = h[(ai//2)*modulus:(ai//2+1)*modulus, :]
            starts = (-ai*xs) % modulus
            for sj in range(len(seeds)):
                occurrences = starts[affine & masks[sj, partner]]
                # Each eligible complementary tuple contributes one arc, including duplicates.
                expanded = sum(((np.arange(modulus)-start) % modulus < width).astype(int)
                               for start in occurrences) if len(occurrences) else np.zeros(modulus)
                assert np.array_equal(expanded, counts[:, sj])
                true = counts[:, sj] % 2
                runs = int(np.sum((true == 1) & (np.roll(true, 1) == 0)))
                assert runs <= len(occurrences)
                run_checks += 1
    print(f"Order {order}: {seed_pairs:,} joint seeds per target, "
          f"rounded existence errors {min(errors):.4f} to {max(errors):.4f}; "
          f"{run_checks:,} offset traces checked.")


def check_incidence_rank():
    for q in (2, 4, 8, 16, 32, 64):
        pivots = {}
        for a, b in product(range(q), repeat=2):
            row = sum(1 << (u*q + (a*u+b) % q) for u in range(q))
            while row:
                column = row.bit_length()-1
                if column not in pivots:
                    pivots[column] = row
                    break
                row ^= pivots[column]
        valuation = 0
        expected = 0
        for j in range(q):
            valuation += j & -j
            expected += max(0, q - valuation)
        assert len(pivots) == expected
        print(f"Modulus {q}: binary incidence rank {expected} out of {q*q}.")


def check_slope_dictionary():
    # The recentered inverse mask reveals its first endpoint in the canonical orientation.
    curve = fixture["curve_cycle"](23, 5)
    order, modulus, width = len(curve), 32, 2
    xs = np.array([-1] + [pt[0] for pt in curve[1:]])
    centers = np.arange(modulus)[:, None]
    directions = list(range(1, modulus//2, 2))
    def translated_masks(direction):
        return (np.isin((xs-centers) % modulus, [(t*direction) % modulus for t in range(width)])
                & (xs != -1))
    masks = [translated_masks(direction) for direction in directions]
    for direction, mask in zip(directions, masks):
        reversed_masks = translated_masks(modulus-direction)
        assert np.array_equal(mask, reversed_masks[(np.arange(modulus)+(width-1)*direction) % modulus])
    # For each center, canonical direction and its negative give identical masks.
    other = np.concatenate(masks).astype(int)
    errors = []
    for target in range(1, order):
        partner = (target-np.arange(order)) % order
        truth = np.array([mask.astype(int) @ other[:, partner].T % 2 for mask in masks])
        # Merge consecutive pairs of canonical directions, keeping this same center.
        counts = truth[::2] + truth[1::2]
        errors.append(int(np.minimum(counts, 2-counts).sum()) / truth.size)
    assert min(errors) > 0
    print(f"Order {order}: exact sign dictionary of {len(directions)} directions; "
          f"merging pairs gives Bayes errors {min(errors):.4f} to {max(errors):.4f}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_progression_overlaps()
    check_mask_metric()
    check_rounding_and_runs()
    check_incidence_rank()
    check_slope_dictionary()
