"""Finite target-averaged checks for adaptive equality readouts of signed moments."""
from itertools import product

from check_address_parity import transform
from fixtures import determinant
from probe_signed_moment import Probe


def mixed_fixture(probe):
    # Each mask is a valid width-two fixed-bit sample. The first mask is
    # also obtained with odd affine multiplier one and the indicated shift.
    slots = []
    for offset in range(0, probe.p-1, 2):
        mask = tuple((point, 2*a + sign) for a in (0, 1)
                     for point, sign in probe.by_x.get(offset+a, ()))
        if len(mask) == 4:
            slots.append(mask)
        if len(slots) == 3:
            break
    assert len(slots) == 3
    return slots


def target_moments(probe, slots, target):
    powers = [[0]*probe.size for _ in range(4)]
    for choices in product(*slots):
        labels, addresses = zip(*choices)
        v = (addresses[0]*probe.block+addresses[1])*probe.block+addresses[2]
        value = determinant([probe.vectors[i] for i in labels]+[probe.vectors[-target]], probe.p)
        for j, table in enumerate(powers):
            table[v] = pow(value, j, probe.p)
    # Walsh transformation on addresses gives all signed moment challenges.
    return list(zip(*[[a % probe.p for a in transform(table)] for table in powers]))


def branches(z, label, r, p):
    z0, z1, z2, z3 = z
    return ((z0*z2-z1*z1-label) % p,
            (z1*z3-z2*z2-label) % p,
            (z3-z1*z2-(r % 3)-label) % p)


def leaf(first, second, r, label):
    return 1 if (first + second + (r >> 1).bit_count() + label) % 2 else -1


def verify(p, linear):
    probe = Probe(p, linear, (0,))
    slots = mixed_fixture(probe)
    targets = list(range(1, probe.n))
    n, width = len(targets), probe.size
    assert n > 16
    scores = [[] for _ in targets]
    default = [[], []]
    advice = [probe.curve[t][0] % 2 for t in targets]
    spectra = [target_moments(probe, slots, t) for t in targets]
    degree = (8, 16, 12)
    exceptional = departures = zero_paths = total_charge = 0
    for r in range(width):
        for label in (0, 1):
            values = [branches(spec[r], label, r, p) for spec in spectra]
            # n>16. Agreement with zero on every target detects an identity
            # for these functions because every tested pole bound is <=16.
            identically_zero = [all(z[j] == 0 for z in values) for j in range(3)]
            for j, is_zero in enumerate(identically_zero):
                if not is_zero:
                    assert sum(z[j] == 0 for z in values) <= degree[j]
            generic_first = identically_zero[0]
            j = 1 if generic_first else 2
            generic_second = identically_zero[j]
            generic_leaf = leaf(generic_first, generic_second, r, label)
            charge = (degree[0] if not generic_first else 0)+(degree[j] if not generic_second else 0)
            total_charge += charge
            zero_paths += int(generic_first)
            default[label].append(generic_leaf)
            changed = 0
            for i, z in enumerate(values):
                if advice[i] != label:
                    continue
                first = z[0] == 0
                second = z[1 if first else 2] == 0
                changed += int((first, second) != (generic_first, generic_second))
                actual = leaf(first, second, r, label)
                departures += int(actual != generic_leaf)
                scores[i].append(actual)
            assert changed <= charge
            exceptional += changed
    assert 0 < departures <= exceptional <= total_charge and 0 < zero_paths < 2*width
    generic_coeffs = [transform(row) for row in default]
    singleton_count = actual_sum = base_sum = base_energy = 0
    for i, t in enumerate(targets):
        address, _ = probe.singleton(slots, t)
        if address is None:
            continue
        observed = transform(scores[i])[address]
        baseline = generic_coeffs[advice[i]][address]
        singleton_count += 1
        actual_sum += abs(observed)
        base_sum += abs(baseline)
        base_energy += baseline*baseline
    assert singleton_count and base_energy <= 2*width**2
    assert base_sum**2 <= 2*singleton_count*width**2
    assert actual_sum <= base_sum + 2*departures
    print(f"p={p}: {n} targets x {width} challenges, {singleton_count} singleton targets, "
          f"{exceptional} departures from generic transcripts ({departures} changed scores).")
    print(f"  Generic/actual singleton coefficient sums in units of 1/{width}: "
          f"{base_sum}/{actual_sum}; checked Parseval and all generic path degree budgets.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    verify(101, 17)
