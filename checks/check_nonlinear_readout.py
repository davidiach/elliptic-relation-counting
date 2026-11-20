"""Exact finite-population comparison of all two-character sign tables and a collision test."""
from fractions import Fraction
from itertools import product

import numpy as np

from check_address_parity import character
from check_nonlinear_branches import mixed_fixture, target_moments
from check_strided_isolation import finite_laws
from fixtures import determinant
from probe_signed_moment import Probe


def sign_tables(bin_count):
    # Opposite tables have identical absolute Walsh coefficients.
    return np.array([(1, *row) for row in product((-1, 1), repeat=bin_count-1)],
                    dtype=np.int64)


def check_target(target):
    probe = Probe(13, 4, (0,))
    affine, fixed = finite_laws(probe.p, probe.curve, probe.positions)
    size = probe.size
    walsh = np.array([[character(a, r) for r in range(size)] for a in range(size)],
                     dtype=np.int64)
    assert np.array_equal(walsh @ walsh, size*np.eye(size, dtype=np.int64))
    chars = np.array(probe.legendre)
    tables = sign_tables(9)
    # Use one common integer denominator for uniform choice among the wrong legal tuples.
    # Here each slot on a singleton has either two or four points.
    common = int(np.lcm.reduce([7, 15, 31, 63]))
    witness_sum = np.zeros(len(tables), dtype=np.int64)
    wrong_sum = np.zeros(len(tables), dtype=np.int64)
    collision = np.zeros((2, 2), dtype=np.int64)
    total = roots = 0
    checked_collision = False
    for slots in product(affine, fixed, fixed):
        mass = affine[slots[0]] * fixed[slots[1]] * fixed[slots[2]]
        total += mass
        marked, legal = probe.singleton(slots, target)
        if marked is None:
            continue
        wrong = [a for a in legal if a != marked]
        assert wrong and common % len(wrong) == 0
        roots += mass
        scale = common // len(wrong)
        values = np.zeros(size, dtype=np.int64)
        for (i, a), (j, b), (k, c) in product(*slots):
            address = (a*probe.block+b)*probe.block+c
            rows = [probe.vectors[v] for v in (i, j, k, -target)]
            values[address] = determinant(rows, probe.p)
        assert values[marked] == 0
        moments = [(walsh @ (values**j)) % probe.p for j in (1, 2)]
        raw_bins = probe.p*moments[0] + moments[1]
        coarse_bins = 3*(chars[moments[0]]+1) + chars[moments[1]]+1
        for index, bins in enumerate((coarse_bins, raw_bins)):
            labels, inverse = np.unique(bins, return_inverse=True)
            fiber = (inverse[:, None] == np.arange(len(labels))).astype(np.int64)
            coefficients = walsh @ fiber
            energy = (coefficients**2).sum(axis=1)
            assert energy.sum() == size**2
            collision[index, 0] += mass*energy[marked]
            collision[index, 1] += mass*scale*energy[wrong].sum()
            if not checked_collision:
                pair_sum = sum(character(marked, r ^ s) for r in range(size)
                               for s in range(size) if bins[r] == bins[s])
                assert pair_sum == energy[marked]
            if index == 0:
                amplitudes = np.abs(tables[:, labels] @ coefficients.T)
                witness_sum += mass*amplitudes[:, marked]
                wrong_sum += mass*scale*amplitudes[:, wrong].sum(axis=1)
        checked_collision = True
    assert total == 8192
    assert roots == {1: 940, 2: 1054}[target]
    tables_shaped = tables.reshape(-1, 3, 3)
    first_only = (tables_shaped == tables_shaped[:, :, :1]).all(axis=(1, 2))
    second_only = (tables_shaped == tables_shaped[:, :1, :]).all(axis=(1, 2))
    best = int(witness_sum.argmax())
    # This separate contrast objective is maximized over the same sign tables.
    gaps = common*witness_sum - wrong_sum
    contrast = int(gaps.argmax())
    den = size*roots
    print(f"p=13, T=[{target}]G: {roots}/{total} guarded singletons, "
          f"{len(tables)} joint binary tables up to global sign.")
    print(f"  Best absolute witness correlation: first only {witness_sum[first_only].max()/den:.6f}, "
          f"second only {witness_sum[second_only].max()/den:.6f}, "
          f"joint {witness_sum[best]/den:.6f} (same-table wrong {wrong_sum[best]/(common*den):.6f}).")
    print(f"  Largest witness-minus-wrong gap among binary tables {gaps[contrast]/(common*den):.6f} "
          f"(witness {witness_sum[contrast]/den:.6f}, wrong {wrong_sum[contrast]/(common*den):.6f}).")
    for label, (hit, miss) in zip(("character pair", "field pair"), collision):
        left = Fraction(int(hit), size*den)
        right = Fraction(int(miss), size*den*common)
        print(f"  Two-challenge energy, {label}: witness {float(left):.6f}, wrong {float(right):.6f}.")


def check_target_average():
    probe = Probe(101, 17, (0,))
    slots = mixed_fixture(probe)
    targets = range(1, probe.n)
    size = probe.size
    walsh = np.array([[character(a, r) for r in range(size)] for a in range(size)])
    # Pack M1 and M2 without identifying any distinct pairs of field values.
    moments = np.array([target_moments(probe, slots, t) for t in targets])
    first, second = moments[:, :, 1], moments[:, :, 2]
    values = first*probe.p + second
    pole_bound = 2*4
    signature = values[:pole_bound+1].T
    _, generic_bin = np.unique(signature, axis=0, return_inverse=True)
    generic_equal = generic_bin[:, None] == generic_bin
    actual_equal = values[:, :, None] == values[:, None, :]
    assert (actual_equal[:, generic_equal]).all()
    extra = actual_equal & ~generic_equal
    max_extra = int(extra.sum(axis=0).max())
    assert 0 < max_extra <= pole_bound
    generic_energy = ((walsh @ generic_equal)*walsh).sum(axis=1)
    assert generic_energy.min() >= 0 and generic_energy.sum() == size**2
    generic_sum = actual_sum = count = 0
    for k, t in enumerate(targets):
        marked, _ = probe.singleton(slots, t)
        if marked is None:
            continue
        count += 1
        signs = walsh[marked]
        actual = int(signs @ actual_equal[k] @ signs)
        assert abs(actual - generic_energy[marked]) <= extra[k].sum()
        generic_sum += int(generic_energy[marked])
        actual_sum += actual
    assert generic_sum <= size**2
    assert actual_sum <= size**2+extra.sum() <= (1+pole_bound)*size**2
    print(f"p=101: {len(targets)} targets, {count} singletons, {size} challenges, "
          f"{len(set(generic_bin))} distinct symbolic moment pairs.")
    print(f"  Most extra target collisions for a distinct pair {max_extra} (bound {pole_bound}); "
          f"sum of singleton energies {actual_sum}/{size**2} (bound {1+pole_bound}).")
    check_character_mixtures(probe, slots, walsh, first, second)


def check_character_mixtures(probe, slots, walsh, first, second):
    small_p = 13
    small_chi = np.array(Probe(small_p, 4, (0,)).legendre)
    translates = small_chi[np.add.outer(np.arange(small_p), np.arange(small_p)) % small_p]
    assert np.array_equal(translates @ translates.T, small_p*np.eye(small_p, dtype=int)-1)
    p = probe.p
    parameters = np.arange(p)[:, None]
    chi = np.array(probe.legendre)
    totals = [0, 0]
    for t in range(1, probe.n):
        marked, _ = probe.singleton(slots, t)
        if marked is None:
            continue
        signs = walsh[marked]
        a, b = second[t-1], first[t-1]**2 % p
        shift = (a+b) % p
        active = b != 0
        ratios = np.array([int(x)*pow(int(y), -1, p) % p for x, y in zip(a[active], b[active])])
        for idx, (rows, labels, mask) in enumerate((
                (chi[(parameters+shift) % p], shift, np.ones_like(active)),
                (chi[(a+parameters*b) % p]*active, ratios, active))):
            coefficients = rows @ signs
            square_sum = int(coefficients @ coefficients)
            same = labels[:, None] == labels
            restricted = signs[mask]
            collision = int(restricted @ same @ restricted)
            assert square_sum == p*collision - int(restricted.sum())**2
            totals[idx] += square_sum
    den = p*probe.size**2
    assert totals[0] <= (1+2*4)*den and totals[1] <= (1+8*4)*den
    print(f"  All {p} character parameters: singleton squared coefficient sums "
          f"{totals[0]}/{den} (shift) and {totals[1]}/{den} (active relative coefficient).")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_target(1)
    check_target(2)
    check_target_average()
