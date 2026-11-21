"""Finite reference contractions and rational-height diagnostics for scalar projection."""

from collections import Counter
from math import gcd, isqrt, prod
from random import Random

from block_unit import PowerBlocks, combine_components
from exception_product import rough_atoms
from prime_projection import (auxiliary_component, auxiliary_product, height_bits,
                              reconstruct_positive, scalar_product, stream_joint_count, unit_log_table)
from target import lift_windows, point_add
from verify_block_unit import GuardedOracle, literal_unit
from verify_exception_product import oriented_weights, reference


def image(atoms, ell, logs):
    def strip(v):
        k = 0
        while v % ell == 0:
            v //= ell
            k += 1
        return k, logs[v % ell]
    delta = character = 0
    for _, a, d in atoms:
        ka, la = strip(a)
        kd, ld = strip(d)
        delta += ka-kd
        character += la-ld
    return delta, character


def atom_summary(atoms):
    # Histograms only in the finite oracle; the graph uses three ordered-pair counts.
    net = Counter(a for _, a, _ in atoms)
    net.subtract(d for _, _, d in atoms)
    return sum(s for s, _, _ in atoms), sum(x*x for x in net.values()), net


def small():
    checks = 0
    for p, a, c in ((211, 1, 1), (167, 1, 3)):
        guard = GuardedOracle(p, a, c)
        P, Q = guard.affine[0], guard.affine[8]
        T = point_add(P, point_add(P, P, p, a), p, a)
        first = (1, -P[0], 1)
        modulus = 1 << p.bit_length()
        base, _ = unit_log_table(p)
        power = PowerBlocks.prepare(p, base, 3)
        classes = {}
        def witness(r):
            if r not in classes:
                matches = [e for e in range(p-1) if power.accepts(e, r)]
                assert len(matches) == 1
                classes[r] = sum(z < matches[0] for z in range(p-1))
            return classes[r]
        cases = [(first,), (first,)*3, (first, (1, 0, modulus)),
                 ((first[0], first[1], 4), (15, -15*Q[0], 4), (23, 0, 4)),
                 (first, (1, -(modulus-1), 1))]
        if p == 211:
            cases += [((1, 0, 2), (1, 205, 8)), ((1, 0, 2), (11, 185, 8))]
        for seeds in cases:
            slots = lift_windows(p, a, c, modulus, seeds)
            full = scalar_product(power, seeds, T, None,
                                  lambda box, target, F, cert: stream_joint_count(
                                      p, a, c, box, target, F, witness))
            assert full == literal_unit(p, a, c, T, slots)
            for factors in ([], [(2, 2), (3, 1)]):
                F = prod(q**e for q, e in factors)
                def expected(part, prefix):
                    return reference(oriented_weights(p, a, T, slots[:len(part)], prefix), p, F)[0]
                want = combine_components(p, seeds, list(map(len, slots)), expected)
                calls = []
                def joint(box, target, cover, cert):
                    assert box == tuple(seeds) and cert == power
                    calls.append(stream_joint_count(p, a, c, box, target, cover, witness))
                    return calls[-1]
                assert scalar_product(power, seeds, T, factors, joint) == want and len(calls) == 1
                def atoms(part, target, prefix, cover):
                    return [rough_atoms(w, p, cover)
                            for w in oriented_weights(p, a, target, slots[:len(part)], prefix)]
                def summary(part, target, prefix, cover):
                    return atom_summary(atoms(part, target, prefix, cover))[:2]
                def projection(part, target, prefix, cover, ell, logs):
                    return image(atoms(part, target, prefix, cover), ell, logs)
                assert auxiliary_product(p, a, c, seeds, T, factors, summary, projection) == want
                if p == 211 and list(map(len, slots)) == [1, 111] and F == 12:
                    print(f"Dense fixture: base {base}, joint count {calls[0]}, "
                          f"complete correction {want}, full unit {full}.")
                checks += 1
    print(f"{checks} seed/cover cases matched the literal correction by both projections.")


def reconstruction():
    cases = 0
    for h in range(1, 34):
        modulus = 1 << (2*h.bit_length()+1)
        for n in range(1, h+1):
            for d in range(1, h+1, 2):
                if gcd(n, d) != 1 or gcd(n, modulus) != 1:
                    continue
                residue = n*pow(d, -1, modulus) % modulus
                assert reconstruct_positive(residue, modulus, h) == (n, d)
                cases += 1
    atoms = [(False, 6, 7), (True, 5, 3)]
    negatives, energy, _ = atom_summary(atoms)
    seen = {}
    def counts(ell, logs):
        v, e = image(atoms, ell, logs)
        seen[ell] = v
        return v, e
    assert auxiliary_component(211, 2, negatives, energy, counts) == -10*pow(7, -1, 211) % 211
    assert seen[3] == 0 and seen[5] == 1 and seen[7] == -1
    print(f"{cases} positive rational reconstructions; cancelled, numerator and denominator primes checked.")


def dense_probe():
    rng = Random(18112026)
    p, a, c, t = 65521, 44, 1, (0, 1)
    g, logs = unit_log_table(p)  # finite oracle only, explicitly field-sized
    power = PowerBlocks.prepare(p, g, 7)
    heights = []
    certificates = []
    for _ in range(12):
        seeds = [(rng.randrange(1, 65536, 2), rng.randrange(65536), w) for w in (8, 16, 16, 16)]
        slots = lift_windows(p, a, c, 65536, seeds)
        def component(part, prefix):
            return reference(oriented_weights(p, a, t, slots[:len(part)], prefix), p, 12)[0]
        target = combine_components(p, seeds, list(map(len, slots)), component)
        assert scalar_product(power, seeds, t, [(2, 2), (3, 1)],
                              lambda box, target, F, cert: stream_joint_count(
                                  p, a, c, box, target, F, logs.__getitem__)) == target
        atoms = [rough_atoms(w, p, 12) for w in oriented_weights(p, a, t, slots, False)]
        _, energy, net = atom_summary(atoms)
        n = prod(k**v for k, v in net.items() if v > 0)
        d = prod(k**(-v) for k, v in net.items() if v < 0)
        reduced = max(n, d)//gcd(n, d)
        exact_h = (reduced-1).bit_length()
        certified_h = height_bits(isqrt(p), len(atoms), energy)
        assert exact_h <= certified_h
        heights.append(exact_h)
        certificates.append(certified_h)
    print(f"12 independent four-slot joint counts matched; terminal reduced rational heights "
          f"{min(heights)}..{max(heights)} bits, certified upper bounds "
          f"{min(certificates)}..{max(certificates)} bits.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    reconstruction()
    small()
    dense_probe()
