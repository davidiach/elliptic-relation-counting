"""Check distributed sparse-radix division, without a full quotient transform."""

from collections import Counter
from fractions import Fraction
from math import ceil, log, log2, pi, sqrt
from random import Random


def owner_partition(owners, offsets, radius):
    """Carve balls of small relative boundary in the remaining owner graph."""
    steps = {d % owners for d in (*offsets, 1, *(-v for v in offsets), -1)}
    remaining = set(range(owners))
    groups = []
    while remaining:
        ball = {min(remaining)}
        shell = set(ball)
        best = set(ball)
        best_ratio = None
        for depth in range(radius):
            next_shell = {(u + d) % owners for u in shell for d in steps}
            next_shell &= remaining - ball
            ratio = Fraction(len(next_shell), len(ball))
            if best_ratio is None or ratio < best_ratio:
                best, best_ratio = set(ball), ratio
            if not next_shell:
                break
            shell = next_shell
            ball.update(shell)
        groups.append(best)
        remaining -= best
    return {owner: group for group, members in enumerate(groups) for owner in members}


def trace(raw, digits, radix, nq):
    """Generate the unique radix digits and carries, flagging a failed equation."""
    q, carry = [], [0]
    inverse = pow(digits[0], -1, radix)
    valid = True
    for i, value in enumerate(raw):
        numerator = value + carry[-1] - sum(
            weight * q[i - d] for d, weight in digits.items() if 0 < d <= i)
        digit = numerator * inverse % radix if i < nq else 0
        next_carry, remainder = divmod(numerator - digits[0] * digit, radix)
        valid &= remainder == 0
        q.append(digit)
        carry.append(next_carry)
    return valid and carry[-1] == 0, q, carry


def signed_presentation(p, radix):
    """Choose a shortest-support signed radix expansion using two carry states."""
    ordinary = []
    while p:
        p, digit = divmod(p, radix)
        ordinary.append(digit)
    paths = {0: ()}
    for value in ordinary + [0]:
        new_paths = {}
        for incoming, prefix in paths.items():
            for outgoing in (0, 1):
                digit = value + incoming - radix*outgoing
                if abs(digit) >= radix:
                    continue
                path = (*prefix, digit)
                if (outgoing not in new_paths or
                        sum(bool(x) for x in path) < sum(bool(x) for x in new_paths[outgoing])):
                    new_paths[outgoing] = path
        paths = new_paths
    return {d: v for d, v in enumerate(paths[0]) if v}


class SparseDivision:
    """Only crossing quotient digits and carries have public graph port IDs."""

    def __init__(self, digits, radix, nq, length, groups, raw_bound):
        self.digits, self.radix, self.nq = digits, radix, nq
        self.owners = [groups[i % len(groups)] for i in range(length)]
        self.by_group = {}
        for i, group in enumerate(self.owners):
            self.by_group.setdefault(group, []).append(i)
        self.raw_bound = raw_bound
        self.carry_bound = 1 + (raw_bound + (radix - 1) * sum(map(abs, digits.values()))) // (radix - 1)
        self.ports = set()
        for i, group in enumerate(self.owners):
            if i and group != self.owners[i-1]:
                self.ports.add(("c", i, group))
            for d in digits:
                j = i - d
                if d and 0 <= j < nq and group != self.owners[j]:
                    self.ports.add(("q", j, group))

    def messages(self, q, carry):
        return {port: (q if port[0] == "q" else carry)[port[1]] for port in self.ports}

    def accepts(self, raw, messages):
        if len(raw) != len(self.owners) or set(messages) != self.ports:
            return False
        if any(abs(v) > self.raw_bound for v in raw):
            return False
        for (kind, _, _), value in messages.items():
            if type(value) is not int:
                return False
            if kind == "q" and not 0 <= value < self.radix:
                return False
            if kind == "c" and abs(value) > self.carry_bound:
                return False
        inverse = pow(self.digits[0], -1, self.radix)
        for group, indices in self.by_group.items():
            q, carries = {}, {0: 0}
            for i in indices:
                incoming = (carries[i] if i == 0 or group == self.owners[i-1]
                            else messages["c", i, group])
                numerator = raw[i] + incoming
                for d, weight in self.digits.items():
                    j = i - d
                    if d == 0 or not 0 <= j < self.nq:
                        continue
                    digit = q[j] if group == self.owners[j] else messages["q", j, group]
                    numerator -= weight * digit
                q[i] = numerator * inverse % self.radix if i < self.nq else 0
                carries[i+1], remainder = divmod(numerator - self.digits[0] * q[i], self.radix)
                if remainder or abs(carries[i+1]) > self.carry_bound:
                    return False
                if i == len(raw)-1 and carries[i+1] != 0:
                    return False
            for (kind, j, _), sent in messages.items():
                if kind == "q" and self.owners[j] == group and sent != q[j]:
                    return False
                if kind == "c" and self.owners[j-1] == group and sent != carries[j]:
                    return False
        return True


def residual_coefficients(value, radix, length, rng):
    raw = [(value // radix**i) % radix for i in range(length)]
    # Redistribute carries so that some valid Fourier coefficients are negative.
    for i in range(length-1):
        c = rng.randrange(-2*radix, 2*radix)
        raw[i] += radix*c
        raw[i+1] -= c
    assert sum(d*radix**i for i, d in enumerate(raw)) == value
    return raw


def check_equations():
    rng = Random(60523)
    examples = [(13, 2, {0: -3, 2: 1}), (29, 2, {0: -3, 2: 2}),
                (101, 3, {0: -3, 1: 5, 2: 1}), (257, 2, {0: 1, 4: 1}),
                (65521, 4, {0: 1, 1: -1, 4: 1})]
    total, valid_count, corruptions = 0, 0, Counter()
    for p, k, digits in examples:
        radix = 1 << k
        assert sum(v * radix**d for d, v in digits.items()) == p
        alternative = signed_presentation(p, radix)
        assert sum(v * radix**d for d, v in alternative.items()) == p and len(alternative) <= len(digits)
        nq = ceil((p.bit_length()+4) / k)
        length = 2*nq + max(digits) + 5
        groups = owner_partition(11, tuple(d for d in digits if d), radius=2)
        checker = SparseDivision(digits, radix, nq, length, groups, 4*radix**2)
        for case in range(300):
            wanted = (0, 16*p-1)[case] if case < 2 else rng.randrange(16*p)
            error = rng.randrange(1,p) if case > 1 and case % 2 else 0
            raw = residual_coefficients(p*wanted+error, radix, length, rng)
            valid, q, carry = trace(raw, digits, radix, nq)
            messages = checker.messages(q, carry)
            assert valid == (error == 0) == checker.accepts(raw, messages)
            if valid:
                assert sum(d * radix**i for i, d in enumerate(q)) == wanted
                valid_count += 1
                if case == 0:
                    for port in sorted(checker.ports):
                        changed = dict(messages)
                        changed[port] = changed[port] ^ 1
                        assert not checker.accepts(raw, changed)
                        corruptions[port[0]] += 1
            total += 1
    assert corruptions["q"] and corruptions["c"]
    print(f"{total} distributed modular checks, {valid_count} exact quotient reconstructions, "
          f"{corruptions['q']} quotient and {corruptions['c']} carry message corruptions rejected")


def check_partitions():
    for owners, offsets, radius in [(64,(5,),6), (257,(8,21),7), (1024,(37,),12)]:
        groups = owner_partition(owners, offsets, radius)
        generators = set((1, *offsets))
        sizes = Counter(groups.values())
        boundary = sum(groups[i] != groups[(i+d) % owners] for i in groups for d in generators)
        g = len(generators)
        assert max(sizes.values()) <= (2*radius + 1)**g
        expansion = (2*radius + 1)**(g/radius)-1
        assert boundary <= 2*g*expansion*owners
        assert set(groups) == set(range(owners))
    print("Three distinct owner graphs met the ball-size and boundary estimates")


def check_compiler_outputs():

    from c38_network import arith

    rng = Random(143231)
    accepted = rejected = 0
    for b, length, offset, m in ((4,8,3,3),(7,16,27,5)):
        data, original, _ = arith.make_instance(b, length, offset, m, rng, "general")
        par = data["par"]
        p, radix = data["constants"]["p"], par["B"]
        digits = signed_presentation(p, radix)
        offsets = tuple(d for d in digits if d)
        checker = SparseDivision(digits, radix, ceil((b+4)/par["k"]), length,
                                 owner_partition(par["L2"], offsets, radius=1), 32*length*radix**2)
        for delta in (0,1,2):
            values = dict(original)
            values["l2"] = (values["l2"] + delta) % p
            for eq in range(arith.arithmetic_count(m)):
                values[f"q{eq}"] = 0
            hs = {key: arith.pack_boundary(arith.H(arith.limbs(value, par), par), par)
                  for key, value in values.items()}
            outputs = arith.middle_outputs(data, hs)
            reference = arith.operations(values, data["constants"], m, mode="general")
            for eq in range(arith.arithmetic_count(m)):
                restored = arith.H_inverse(arith.unpack_boundary(outputs[eq], par), par)
                raw = [v if 2*v < par["modulus"] else v-par["modulus"] for v in restored]
                assert sum(v*radix**i for i,v in enumerate(raw)) == reference[eq]
                valid, q, carry = trace(raw, digits, radix, checker.nq)
                assert checker.accepts(raw, checker.messages(q, carry)) == valid == (reference[eq] % p == 0)
                if delta == 0:
                    assert valid and sum(v*radix**i for i,v in enumerate(q)) == original[f"q{eq}"]
                accepted += valid
                rejected += not valid
    assert accepted and rejected
    print(f"Two C38 Fourier fixtures produced {accepted} divisible raw outputs and "
          f"{rejected} rejected outputs after changing slopes")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled")
    check_equations()
    check_partitions()
    check_compiler_outputs()
    K = 2*pi*sqrt(14*log2(90))/(9*log(2))
    print(f"Bounded-support conditional coefficient {2*sqrt(10*K):.12f}")
