"""Finite checks of fresh_windows.tex. Requires NumPy and enabled assertions."""

from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
from runpy import run_path

import numpy as np

fixtures = run_path(str(Path(__file__).resolve().parents[1] / "mixed_coefficient" / "verify.py"))


def check_zoom():
    cases = 0
    for bits in range(2, 9):
        for exponent in range(1, bits):
            width = 1 << exponent
            for depth in range(exponent):
                length = width >> depth
                half = length // 2
                multiplier = (2 << depth) - 1
                image = [half + multiplier*y + eta
                         for y in range(length) for eta in range(multiplier)]
                assert image == list(range(half, 2*width - half))
                assert 2*width <= 1 << bits
                assert all((value < width) == (i // multiplier < half)
                           for i, value in enumerate(image))
                cases += 1
    print(f"Zoom classification and uniform marked image: {cases} dyadic cases.")


def check_actual_hashes():
    curve = fixtures["curve_cycle"](13, 4)
    order, modulus = len(curve), 16
    width = 4
    seeds = list(product(range(1, modulus, 2), range(modulus)))
    where = {seed: i for i, seed in enumerate(seeds)}
    xs = np.array([-1] + [pt[0] for pt in curve[1:]])
    alpha, beta = np.array(seeds).T
    addresses = (alpha[:, None]*xs + beta[:, None]) % modulus
    affine = np.arange(order) != 0
    small = ((addresses < width) & affine).astype(int)
    big = ((addresses < 2*width) & affine).astype(int)
    trials = len(seeds)**2
    # Introduce a fixed set of systematic errors on fresh seeds.
    faults = (np.arange(len(seeds))[:, None] + 3*np.arange(len(seeds))) % 79 == 0
    examples = 0
    for target in range(1, order):
        partner = (target - np.arange(order)) % order
        right = small[:, partner]
        root = small @ right.T
        enlarged = big @ right.T
        predictor = (root % 2).astype(bool) ^ faults
        factorial = int(np.sum(enlarged*(enlarged - 1)))
        for depth in range(2):
            length = width >> depth
            half, q = length // 2, (2 << depth) - 1
            bad = 0
            for point in range(1, order):
                left_rows = np.flatnonzero(small[:, point])
                right_rows = np.flatnonzero(right[:, point])
                if not len(right_rows):
                    continue
                u = addresses[left_rows, point]
                j = (u // length) * length
                desired = (u - j < half)[:, None]
                for eta in range(q):
                    new_alpha = (q*alpha[left_rows]) % modulus
                    new_beta = (q*(beta[left_rows]-j) + half + eta) % modulus
                    new_rows = [where[pair] for pair in zip(new_alpha, new_beta)]
                    observed = predictor[np.ix_(new_rows, right_rows)]
                    bad += int(np.count_nonzero(observed != desired))
            assert bad <= q * (factorial + np.count_nonzero(faults))
            examples += 1
    print(f"Marked error inequality: {examples} target-depth checks, "
          f"all {trials:,} actual seed pairs per target.")

    # Two narrower root slots and fourfold enlargement in either one.
    w = 2
    root_mask = ((addresses < w) & affine).astype(int)
    clean_count = 0
    for target in range(1, order):
        partner = (target - np.arange(order)) % order
        right = root_mask[:, partner]
        h = root_mask @ right.T
        guarded_points = affine & (partner != 0) & (xs != xs[partner])
        guarded = (root_mask * guarded_points) @ right.T
        big_left = big @ right.T
        big_right = root_mask @ big[:, partner].T
        clean = (guarded == 1) & (big_left == 1) & (big_right == 1)
        assert np.all(clean.astype(int) >= guarded - h*(big_left + big_right - 2))
        clean_count += int(clean.sum())
        for residue in range(3):
            inv = pow(3, -1, modulus)
            moved = [where[(int(inv*a % modulus), int(inv*(b-residue) % modulus))]
                     for a, b in seeds]
            assert sorted(moved) == list(range(len(seeds)))
            original_residue = root_mask * (addresses % 3 == residue)
            left_answer = root_mask[moved] @ right.T
            right_answer = root_mask @ root_mask[moved][:, partner].T
            assert np.array_equal(left_answer[clean] % 2, (original_residue @ right.T)[clean])
            desired = root_mask @ original_residue[:, partner].T
            assert np.array_equal(right_answer[clean] % 2, desired[clean])
    assert clean_count > 0
    print(f"Same-width residue probes: all {clean_count:,} clean target/seed triples passed.")


def check_buckets():
    for k, a, c in ((4, 2, 3), (4, 4, 3), (5, 2, 3)):
        numerator = 0
        odd_outputs = 0
        for left in product(range(k), repeat=a):
            odd = {i for i in range(k) if left.count(i) % 2}
            numerator += k**c - (k - 2*len(odd))**c
            for right in product(range(k), repeat=c):
                odd_outputs += sum(i in odd for i in right) % 2
        total = k**(a+c)
        assert 2 * odd_outputs == numerator
        collision_mean = Fraction(a*c, k)
        assert (collision_mean - Fraction(a*c*(a*c-1), k*k)
                <= Fraction(odd_outputs, total) <= collision_mean)
    print("Exact one-pass bucket parity error matched on three disjoint-support fixtures.")


def check_repeated_sketch_law():
    def walk(k, steps):
        counts = [1] + [0] * (2**k-1)
        for _ in range(steps):
            counts = [sum(counts[v ^ (1 << j)] for j in range(k)) for v in range(2**k)]
        return counts

    for k in range(2, 7):
        for steps in range(1, 21):
            counts = walk(k, steps)
            total = k**steps
            diff = [Fraction(count, total) -
                    (Fraction(2, 2**k) if (v.bit_count()-steps) % 2 == 0 else 0)
                    for v, count in enumerate(counts)]
            spectrum = sum(comb(k, j)*Fraction(k-2*j, k)**(2*steps) for j in range(1, k))
            assert 2**k * sum(d*d for d in diff) == spectrum
            assert sum(abs(d) for d in diff)**2 <= spectrum
    # Both support sizes are four. Their parity-vector marginals are even
    # regardless of how many support elements are shared.
    k = 3
    for overlap in (0, 1, 2):
        counts = np.zeros((2**k, 2**k), dtype=int)
        for labels in product(range(k), repeat=8-overlap):
            v = w = 0
            for h in labels[:4]:
                v ^= 1 << h
            for h in labels[4-overlap:]:
                w ^= 1 << h
            counts[v, w] += 1
        unmatched = walk(k, 4-overlap)
        shared = walk(k, overlap)
        for v, w in np.ndindex(counts.shape):
            predicted = sum(n * unmatched[v ^ s] * unmatched[w ^ s] for s, n in enumerate(shared))
            assert counts[v, w] == predicted
            if v.bit_count() % 2 or w.bit_count() % 2:
                assert counts[v, w] == 0
    print("Bucket-walk spectral identity: 100 exact laws and three joint overlap cases.")


def check_copy_tensor():
    matrices = [np.array(v).reshape(2, 2) for v in product((0, 1), repeat=4)
                if (v[0]*v[3] - v[1]*v[2]) % 2]
    equality = np.zeros((2, 2, 2), dtype=int)
    equality[0, 0, 0] = equality[1, 1, 1] = 1
    for a, b, c in product(matrices, repeat=3):
        f = np.einsum("ia,jb,kc,abc->ijk", a, b, c, equality) % 2
        mixed = sum(f[0, j, k] * f[1, 1-j, 1-k] for j, k in product(range(2), repeat=2))
        assert mixed % 2 == 1
        assert {sum(idx) % 2 for idx in np.ndindex(f.shape) if f[idx]} == {0, 1}
    print(f"Ternary equality invariant retained in all {len(matrices)**3} binary bases.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_zoom()
    check_actual_hashes()
    check_buckets()
    check_repeated_sketch_law()
    check_copy_tensor()
