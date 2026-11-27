"""Finite checks of the cost expressions, tied-label expectations and balanced cuts.

Run with Python 3.10+ and assertions enabled.
"""
from fractions import Fraction
from itertools import product
import math


def verify_proposition_b1():
    # t is a power of two, so S = 8mb and beta = 4b also have exact integer logs.
    samples = (64, 128, 256)
    for t in samples:
        L = t**6
        m = t
        x = t**2  # log2(R)
        lam = L + 3 + (m.bit_length() - 1)
        log_beta = L + 2
        u = lam // (2**17 * x)
        a = 2 * u
        log_D = 2 * a * x
        log_H = lam - 0.5 * math.log2(lam)
        log_delta = 0.9 * lam

        # Evaluate the root packing and precision estimates in log space.
        assert u >= 1 and x >= math.log2(max(64, 4 * m))
        assert log_H > math.log2(a + 2) + 2 * log_D + log_delta
        assert log_H >= a + log_D and a + 1 <= log_delta
        assert lam + a + 1 <= 4 * lam
        assert u - math.log2(u) > 0.5 * math.log2(lam) + 2
        assert a * x * 2**17 <= 2 * lam and 2 * log_D * 2**17 <= 8 * lam
        assert 2 * (x + 1) + math.log2(lam) + 2 < log_beta

        # The three normalized core costs agree exactly at these parameters.
        assert Fraction(1, m) == m * Fraction(math.isqrt(x), math.isqrt(L))
        assert Fraction(1, m) == Fraction(1, math.isqrt(x))
        gap_core = min(x, (log_beta - math.log2(a + 2)) / 2) - 1
        gap_lookup = (4 * m).bit_length() - 1
        rho_core = max(1, math.floor(min(gap_core, log_beta / 10)))
        rho_with_lookups = max(1, min(rho_core, gap_lookup))
        assert rho_with_lookups == gap_lookup < rho_core

    # Fix total log fan-out and enumerate the positive parts of two- and three-level splits.
    for total in (10, 20, 40):
        for depth in (2, 3):
            for prefix in product(range(1, total), repeat=depth - 1):
                last = total - sum(prefix)
                if last > 0:
                    cost = sum(2**j * math.sqrt(x) for j, x in enumerate((*prefix, last)))
                    assert cost > math.sqrt(total)
    print(f"Proposition B1: core budget terms agree and sampled guards pass for m={samples}; "
          "the lookup degree caps rho at log2(4m) in each sample.")
    print("Proposition B1: one level minimizes the comparison expression in all tested "
          "two- and three-level splits.")


def verify_proposition_b2():
    # Count class collisions directly under full and first-coordinate ties.
    for R, a in ((2, 1), (2, 2), (3, 1), (3, 2), (4, 1)):
        labels = list(product(range(R), repeat=a))
        q = len(labels)
        D = q**2
        for tied_coordinates in sorted({1, a}):
            same_position = 0
            opposite_positions = 0
            samples = 0
            for h1, h2 in product(labels, repeat=2):
                if h1[:tied_coordinates] != h2[:tied_coordinates]:
                    continue
                for other1, other2 in product(labels, repeat=2):
                    edge1 = (h1, other1)
                    same_position += edge1 == (h2, other2)
                    opposite_positions += edge1 == (other2, h2)
                    samples += 1
            independent = Fraction(1, D)
            assert Fraction(same_position, samples) - independent == Fraction(R**tied_coordinates - 1, D)
            assert Fraction(opposite_positions, samples) - independent == 0
    print("Proposition B2: same-position covariance is 1/q-1/D for a full tie and (R-1)/D "
          "for one coordinate; opposite-position covariance is zero in both models.")

    # Two source clusters, two edges per cluster, and four independently labeled leaves.
    # With two symbols there are four classes; a perfectly balanced load is one per class.
    potentials = []
    for h1, h2, v1, v2, v3, v4 in product(range(2), repeat=6):
        edges = [(h1, v1), (h1, v2), (h2, v3), (h2, v4)]
        potentials.append(sum((edges.count(pair) - 1)**2 for pair in product(range(2), repeat=2)))
    mean = Fraction(sum(potentials), len(potentials))
    assert mean == 4 and min(potentials) == 0
    print(f"Proposition B2: across {len(potentials)} four-edge labelings, the mean packing "
          f"potential is {mean} and the minimum is {min(potentials)}.")


def verify_theorem_b3_bipartite_cutwidth():
    fractions = []
    for left, right in ((64, 64), (128, 64), (256, 128), (64, 256)):
        edges = left * right
        cuts = []
        for a in range(left + 1):
            for c in range(right + 1):
                weight = a * right + c * left
                if 2 * edges <= 3 * weight <= 4 * edges:
                    cuts.append(a * (right - c) + c * (left - a))
        fraction = Fraction(min(cuts), edges)
        assert fraction >= Fraction(4, 9)
        fractions.append(fraction)
    print(f"Theorem B3: minimum balanced cut fraction over four tested aspect ratios is "
          f"{min(fractions)} >= 4/9.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without Python optimization so assertions execute.")
    verify_proposition_b1()
    verify_proposition_b2()
    verify_theorem_b3_bipartite_cutwidth()
    print("CHECK_RESIDUAL_BARRIERS: ALL CHECKS PASSED.")
