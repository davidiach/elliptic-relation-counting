"""Exact finite checks for the block certificate. No contraction timings are inferred."""

from itertools import product
from random import Random

from block_unit import PowerBlocks, unit_from_counts
from divisor_unit import public_remainder, unit_from_divisors
from target import chord_scale, lift_windows, point_add
from scalar import form, polynomial
from arithmetic_alternating_chord import points


def check_power_gates():
    cases = [(3, 2), (17, 3), (167, 5), (211, 2), (1019, 2)]
    comparisons = 0
    for p, g in cases:
        for width in {1, 3, 5, p.bit_length()+1}:
            cert = PowerBlocks.prepare(p, g, width)
            assert all(pow(g, i, p) != 1 for i in range(1, p-1))
            residue = 1
            for exponent in range(p-1):
                assert cert.accepts(exponent, residue)
                assert not cert.accepts(exponent, residue % (p-1)+1)
                assert not cert.accepts(exponent, 0)
                comparisons += 1
                residue = residue*g % p
            assert residue == 1
            assert not cert.accepts(p-1, 1) and not cert.accepts(-1, 1)
            for block, bits in enumerate(cert.widths()):
                for digit in range(1 << bits):
                    value = pow(g, digit << (width*block), p)
                    assert cert.matches_slices(block, digit, value, 3)
                    for k in range(p.bit_length()):
                        assert not cert.matches_slices(block, digit, value ^ (1 << k), 3)
                assert not cert.matches_slices(block, 1 << bits, 1, 3)
    print(f"{comparisons} bounded exponent cases across five fields; corrupted output bits rejected.")


def check_wide_words():
    rng = Random(20261118)
    p = (1 << 64)-(1 << 32)+1
    assert p-1 == (1 << 32)*3*5*17*257*65537
    assert all(pow(7, (p-1)//q, p) != 1 for q in (2, 3, 5, 17, 257, 65537))
    for bits in (5, 16, 27):
        cert = PowerBlocks.prepare(p, 7, bits)
        for e in [0, p-2, *(rng.randrange(p-1) for _ in range(48))]:
            assert cert.accepts(e, pow(7, e, p))
            for i, digit in enumerate(cert.digits(e)):
                value = pow(7, digit << (i*bits), p)
                assert cert.matches_slices(i, digit, value, 7)
    for width in (3, 8, 31, 64):
        for _ in range(128):
            cover = rng.getrandbits(rng.randrange(1, 1025))
            divisor = rng.randrange(1, 1 << width)
            assert public_remainder(cover, divisor, width) == cover % divisor
    print("150 block checks over a 64-bit field and 512 multi-limb remainder checks passed.")


def negate(P, p):
    return P[0], -P[1] % p


def literal_unit(p, a, c, target, slots):
    z_minus_t = polynomial([-target[0], 1], p)
    unit = 1
    for row in product(*slots):
        f = form([P[0] for P in row], p, a, c)
        while not f(target[0]):
            f = f//z_minus_t
        unit = unit*int(f(target[0])) % p
    return unit


class GuardedOracle:
    """Stream signed test tuples and offsets; accumulate only bit counts."""

    def __init__(self, p, a, c):
        self.p, self.a, self.c = p, a, c
        self.affine = points(p, a, c)
        self.anchor = self.affine[0]
        self.box_calls = 0
        self.divisor_boxes = 0
        self.nonempty_exceptions = 0
        self.branches = set()

    def good_offsets(self, path):
        count = 0
        for offset in self.affine:
            total = offset
            for Q in (*path, negate(offset, self.p)):
                if total[0] == Q[0]:
                    break
                total = point_add(total, Q, self.p, self.a)
            else:
                assert total == self.anchor
                count += 1
        return count

    def weight_paths(self, box, target, prefix):
        p = self.p
        slots = lift_windows(p, self.a, self.c, 1 << p.bit_length(), box)
        signed = [[Q for P in slot for Q in (P, negate(P, p))] for slot in slots]
        for row in product(*signed):
            path = row[:-1] if prefix else row
            total = None
            for Q in path:
                total = point_add(total, Q, p, self.a)
            if total is None:
                weight = 1 if prefix else p-1
                body = (*path, self.anchor)
                self.branches.add("prefix infinity" if prefix else "terminal infinity")
            else:
                if prefix:
                    weight = chord_scale(total, row[-1], p)
                    self.branches.add("tangent" if total[0] == row[-1][0] else "chord")
                else:
                    weight = (target[0]-total[0]) % p or 1
                    self.branches.add("target" if target[0] == total[0] else "ordinary")
                body = (*path, negate(total, p), self.anchor)
            yield weight, body

    def __call__(self, box, target, prefix, cert):
        self.box_calls += 1
        counts = [0]*cert.bit_count
        for weight, body in self.weight_paths(box, target, prefix):
            accepted = [e for e in range(self.p-1) if cert.accepts(e, weight)]
            assert len(accepted) == 1
            valid = self.good_offsets(body)
            for k in range(cert.bit_count):
                if accepted[0] & (1 << k):
                    counts[k] += valid
        return counts

    def divisors(self, box, target, prefix, cover, pins):
        self.divisor_boxes += 1
        counts = [0]*len(pins)
        for weight, body in self.weight_paths(box, target, prefix):
            matched = [k for k, pin in enumerate(pins) if pin.matches(weight, self.p, cover)]
            if matched:
                valid = self.good_offsets(body)
                for k in matched:
                    counts[k] += valid
        self.nonempty_exceptions += sum(pin.kind == "exception" and value > 0
                                        for pin, value in zip(pins, counts))
        return counts


def check_seeds():
    p, a, c = 211, 1, 1
    oracle = GuardedOracle(p, a, c)
    P, Q = oracle.affine[0], oracle.affine[8]
    S = point_add(P, Q, p, a)
    T = point_add(P, point_add(P, P, p, a), p, a)
    atoms = tuple((alpha, -alpha*point[0] % 256, 4)
                  for alpha, point in zip((5, 15, 23), (P, Q, S)))
    cases = [(atoms[:2], S, 59), (atoms, P, 53), ((atoms[0],)*3, P, 61),
             ((atoms[0],)*2, point_add(P, P, p, a), 84),
             (((1, 0, 2), (1, 205, 8)), T, 104),
             (((1, 0, 2), (11, 185, 8)), T, 200),
             (((1, 0, 1),), P, 1), (((1, -254, 1),), P, 1)]

    def check(oracle, base, cases):
        p, a, c = oracle.p, oracle.a, oracle.c
        answers = []
        # This public control covers every centered magnitude on the tiny field.
        full = [(q, p.bit_length()) for q in range(2, p//2+1)
                if all(q % d for d in range(2, int(q**0.5)+1))]
        for seeds, target, reference in cases:
            widths = (2, 5)
            slots = lift_windows(p, a, c, 1 << p.bit_length(), seeds)
            expected = literal_unit(p, a, c, target, slots)
            if reference is not None:
                assert expected == reference, (seeds, target, expected, reference)
            for width in widths:
                cert = PowerBlocks.prepare(p, base, width)
                result = unit_from_counts(cert, a, c, len(oracle.affine)+1, seeds, target, oracle)
                assert result == expected, (seeds, target, width, result, expected)
            for factors in [full, [(2, 2), (3, 1)]]:
                before = oracle.nonempty_exceptions
                result = unit_from_divisors(p, a, c, len(oracle.affine)+1, seeds, target,
                                            factors, oracle.divisors)
                assert result == expected, (seeds, factors, result, expected)
                if factors == full:
                    assert before == oracle.nonempty_exceptions
            answers.append(expected)
        assert oracle.nonempty_exceptions > 0
        print(f"F_{p}: exact units {answers}; {oracle.box_calls} block box vectors and "
              f"{oracle.divisor_boxes} divisor box vectors, two widths and two covers.")
        return answers

    assert len(oracle.affine)+1 == 223
    check(oracle, 2, cases)
    assert oracle.branches == {"prefix infinity", "terminal infinity", "tangent",
                               "chord", "target", "ordinary"}
    # A full-width last slot realizes every affine V. Fix the first slot at P.
    # Count distinct exceptional magnitudes without storing a weight distribution.
    def occurs(magnitude):
        for V in oracle.affine:
            total = point_add(P, V, p, a)
            if total is None or total[0] == T[0]:
                continue
            w = (T[0]-total[0]) % p
            if min(w, p-w) == magnitude:
                return True
        return False

    distinct = sum(12 % h != 0 and occurs(h) for h in range(1, p//2+1))
    lower = len(oracle.affine)//4 - 6  # ceil((2n-3)/4) = floor(2n/4), tau(12)=6
    assert distinct >= lower > 0
    print(f"Full last slot over F_211, cover 12: {distinct} exceptional magnitudes, lower bound {lower}.")
    second = GuardedOracle(167, 1, 3)
    assert len(second.affine)+1 == 181
    P, Q = second.affine[0], second.affine[9]
    seeds = ((7, -7*P[0], 4), (11, -11*Q[0], 4))
    values = check(second, 5, [(seeds, P, None), (seeds, point_add(P, Q, 167, 1), None)])
    assert any(pow(x, 2, 167) != 1 for x in values)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_power_gates()
    check_wide_words()
    check_seeds()
