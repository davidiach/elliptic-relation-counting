"""Independent finite oracles for joint exceptional-weight products."""

from bisect import bisect_left
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod
from random import Random

from block_unit import combine_components
from exception_product import (Bracket, atom_pins, balanced_component, exception_product,
                               exception_product_balanced, exception_product_weighted,
                               farey_bracket, farey_component, prime_columns, ratio_atoms, rough_atoms)
from target import chord_scale, lift_windows, point_add
from verify_block_unit import GuardedOracle, literal_unit, negate


def check_coordinates():
    cases = 0
    for p in (3, 5, 7, 13, 31, 43):
        for cap in range(1, p):
            fractions = sorted({Fraction(a, d) for d in range(1, cap+1) for a in range(d+1)})
            for x in range(1, p):
                idx = bisect_left(fractions, Fraction(x, p))
                l, u = fractions[idx-1:idx+1]
                expected = Bracket(l.numerator, l.denominator, u.numerator, u.denominator)
                assert farey_bracket(x, p, cap) == expected and expected.valid(x, p, cap)
                assert not Bracket(l.numerator, l.denominator, l.numerator, l.denominator).valid(x, p, cap)
                cases += 1
    rng = Random(18112026)
    for p in (167, 211, 65521, (1 << 127)-1):
        residues = range(1, p) if p < 300 else [1, 2, p-1, *(rng.randrange(1, p) for _ in range(256))]
        for r in residues:
            negative, a, d = ratio_atoms(r, p)
            assert 1 <= min(a, d) <= max(a, d) <= isqrt(p) and gcd(a, d) == 1
            assert r*d % p == (-a if negative else a) % p
            if r != pow(r, -1, p):
                assert ratio_atoms(pow(r, -1, p), p) == (negative, d, a)
    print(f"{cases} brackets agree with exhaustive Farey neighbours; ratio checks through 127 bits.")


def oriented_weights(p, a, target, slots, prefix):
    choices = []
    for j, slot in enumerate(slots):
        if j == 0 or (prefix and j == len(slots)-1):
            choices.append(slot)
        else:
            choices.append([v for q in slot for v in (q, negate(q, p))])
    for row in product(*choices):
        total = None
        for v in (row[:-1] if prefix else row):
            total = point_add(total, v, p, a)
        if prefix:
            yield chord_scale(total, row[-1], p)
        elif total is None:
            yield p-1
        else:
            yield (target[0]-total[0]) % p or 1


def reference(weights, p, cover):
    correction = covered = 1
    for w in weights:
        h = min(w, p-w)
        d = gcd(h, cover)
        correction = correction*(h//d) % p
        covered = covered*(d if w == h else -d) % p
    return correction, covered


def raw_atom_oracle(guard):
    def query(box, target, prefix, cover, pins):
        counts = [0]*len(pins)
        for weight, body in guard.weight_paths(box, target, prefix):
            atoms = rough_atoms(weight, guard.p, cover)
            matched = [j for j, pin in enumerate(pins) if pin.matches(atoms)]
            if matched:
                valid = guard.good_offsets(body)
                for j in matched:
                    counts[j] += valid
        return counts
    return query


def collision_counts(atoms, lo, hi):
    """Literal ordered pairs for the verifier. The reduction stores no atom support."""
    left = [a for _, a, _ in atoms if lo <= a <= hi]
    right = [d for _, _, d in atoms if lo <= d <= hi]
    def matches(xs, ys):
        return sum(x == y for x in xs for y in ys)
    return len(left), len(right), matches(left, left), matches(left, right), matches(right, right)


def pair_check(p, atoms):
    calls = 0
    def interval(lo, hi):
        nonlocal calls
        calls += 1
        return collision_counts(atoms, lo, hi)
    sign = sum(negative for negative, _, _ in atoms)
    return balanced_component(p, sign, interval), calls


def verify_seed_cases():
    for p, a, c in [(211, 1, 1), (167, 1, 3)]:
        guard = GuardedOracle(p, a, c)
        P, Q = guard.affine[0], guard.affine[8]
        T = point_add(P, point_add(P, P, p, a), p, a)
        modulus = 1 << p.bit_length()
        first = (1, -P[0], 1)
        slots3 = tuple((alpha, -alpha*v[0], 4)
                       for alpha, v in zip((5, 15, 23), (P, Q, point_add(P, Q, p, a))))
        cases = [(slots3, P), ((first,)*3, P), ((first,), P),
                 (((1, -(modulus-1), 1),), P)]
        if p == 211:
            cases.extend([(((1, 0, 2), (1, 205, 8)), T), (((1, 0, 2), (11, 185, 8)), T),
                          ((first, (1, 0, modulus)), T)])
        for seeds, target in cases:
            lifted = lift_windows(p, a, c, modulus, seeds)
            sizes = list(map(len, lifted))
            expected_unit = literal_unit(p, a, c, target, lifted)
            for factors in ([], [(2, 2), (3, 1)]):
                cover = prod(q**e for q, e in factors)
                raw = raw_atom_oracle(guard)
                def expected(part, prefix):
                    weights = oriented_weights(p, a, target, lifted[:len(part)], prefix)
                    return reference(weights, p, cover)
                want = combine_components(p, seeds, sizes, lambda part, prefix: expected(part, prefix)[0])
                covered = combine_components(p, seeds, sizes, lambda part, prefix: expected(part, prefix)[1])
                assert want*covered % p == expected_unit
                got = exception_product(p, a, c, len(guard.affine)+1, seeds, target, factors, raw)
                assert got == want
                def intervals(part, target, prefix, cover, lo, hi):
                    weights = oriented_weights(p, a, target, lifted[:len(part)], prefix)
                    return collision_counts([rough_atoms(w, p, cover) for w in weights], lo, hi)
                assert exception_product_balanced(p, a, c, len(guard.affine)+1, seeds, target,
                                                  factors, raw, intervals) == want
                def weighted(part, target, prefix, cover, powers):
                    count = 0
                    for w in oriented_weights(p, a, target, lifted[:len(part)], prefix):
                        negative, n, d = rough_atoms(w, p, cover)
                        if powers:
                            count += sum((n % power == 0) + (d % power != 0) for power in powers)
                        else:
                            count += negative
                    return count
                assert exception_product_weighted(p, a, c, seeds, target, factors, weighted) == want
                if all(sizes):
                    for k, prefix in [(len(seeds), False), *[(j, True) for j in range(2, len(seeds)+1)]]:
                        weights = list(oriented_weights(p, a, target, lifted[:k], prefix))
                        atoms = [rough_atoms(w, p, cover) for w in weights]
                        paired, calls = pair_check(p, atoms)
                        assert paired == reference(weights, p, cover)[0]
                        if p == 211 and sizes == [1, 111] and cover == 12 and not prefix:
                            exc = {min(w, p-w) for w in weights if cover % min(w, p-w)}
                            direct = farey_component(p, 223, seeds, target, False, cover, raw, batch_size=3)
                            assert len(exc) == 78 and direct == paired
                            print(f"Dense terminal: {len(exc)} distinct exceptions, correction {paired}, "
                                  f"{len(list(atom_pins(isqrt(p))))} guarded pin types, "
                                  f"{1+len(list(prime_columns(isqrt(p))))} oriented weighted counts, "
                                  f"{calls} balance intervals. "
                                  f"Combined prefix and terminal correction {got}, complete unit {expected_unit}.")
                if p == 211 and sizes == [1, 5] and cover == 12:
                    print(f"Collision fixture: correction {got}, covered factor {covered}, full unit {expected_unit}.")
        print(f"F_{p}: {2*len(cases)} seed/cover cases matched the literal product.")


def check_balance_extremes():
    p = 211
    residues = [r for r in range(2, p//2+1) if r < pow(r, -1, p) <= p//2]
    weights = [v for r in residues for v in (r, pow(r, -1, p))]
    answer, calls = pair_check(p, [ratio_atoms(w, p) for w in weights])
    assert answer == 1 and calls == 1
    print(f"{len(weights)} distinct reciprocal-paired rough values cancel in one balance interval.")
    p, a, c, T = 65521, 44, 1, (0, 1)
    window = lift_windows(p, a, c, 65536, ((1, 0, 64),))
    weights = list(oriented_weights(p, a, T, window, False))
    magnitudes = [min(w, p-w) for w in weights]
    atoms = [ratio_atoms(h, p) for h in magnitudes]
    assert all(atom == (False, h, 1) for atom, h in zip(atoms, magnitudes))
    k = sum(h > 1 for h in set(magnitudes)) + 1
    answer, calls = pair_check(p, atoms)
    assert k > 20 and calls >= k and answer == prod(magnitudes) % p
    print(f"One width-64 seed over F_65521: {k} imbalanced labels, {calls} balance intervals.")


def independent_seed_probe():
    """Finite reference histograms to measure cancellation in independently drawn slots."""
    p, a, c, target = 65521, 44, 1, (0, 1)
    rng = Random(18112026)
    observed = []
    columns = list(prime_columns(isqrt(p)))
    # Public small-integer divisibility table for the finite reference calculation.
    factors = [[] for _ in range(isqrt(p)+1)]
    for i, (_, powers) in enumerate(columns):
        for power in powers:
            for value in range(power, len(factors), power):
                factors[value].append(i)
    for _ in range(12):
        seeds = [(rng.randrange(1, 65536, 2), rng.randrange(65536), w) for w in (8, 16, 16, 16)]
        slots = lift_windows(p, a, c, 65536, seeds)
        if not all(slots):
            continue
        exceptional = set()
        net = Counter()
        sign = 0
        direct = 1
        tallies = [0]*len(columns)
        for w in oriented_weights(p, a, target, slots, False):
            h = min(w, p-w)
            rough = h//gcd(h, 12)
            if rough != 1:
                exceptional.add(h)
            negative, n, d = ratio_atoms(rough, p)
            net[n] += 1
            net[d] -= 1
            sign += negative
            direct = direct*rough % p
            for i in factors[n]:
                tallies[i] += 1
            for i in factors[d]:
                tallies[i] -= 1
        got = -1 if sign % 2 else 1
        for value, mass in net.items():
            got = got*pow(value, mass, p) % p
        assert got == direct
        prime_answer = -1 if sign % 2 else 1
        for (q, _), tally in zip(columns, tallies):
            prime_answer = prime_answer*pow(q, tally, p) % p
        assert prime_answer == direct
        observed.append((len(exceptional), sum(v != 0 for v in net.values())))
    assert observed
    ds, ks = zip(*observed)
    print(f"Independent four-slot terminal probe: {len(observed)}/12 nonempty; exceptions "
          f"{min(ds)}..{max(ds)}, imbalanced labels {min(ks)}..{max(ks)} out of {isqrt(p)}. "
          f"All prime-coordinate products matched, with {len(columns)+1} coordinates including sign.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_coordinates()
    verify_seed_cases()
    check_balance_extremes()
    independent_seed_probe()
