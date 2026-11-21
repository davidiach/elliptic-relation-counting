"""Exact folded products, conjugate-paired Fourier boundaries, and wrapped carries."""
from math import isqrt
from random import Random

import arithmetic_alternating_chord as chord


def prime(value):
    return value > 1 and all(value % d for d in range(2, isqrt(value)+1))


class QuadraticRing:
    """Unramified quadratic ring (Z/q^e)[i] with q = 3 mod 4."""
    def __init__(self, q, exponent):
        self.q = q
        self.modulus = q**exponent

    def scale(self, a, scalar):
        return (a[0]*scalar % self.modulus, a[1]*scalar % self.modulus)

    def times(self, a, b):
        return ((a[0]*b[0]-a[1]*b[1]) % self.modulus,
                (a[0]*b[1]+a[1]*b[0]) % self.modulus)

    def conjugate(self, a):
        return a[0], -a[1] % self.modulus

    def inverse(self, a):
        norm = (a[0]**2+a[1]**2) % self.modulus
        return self.scale(self.conjugate(a), pow(norm, -1, self.modulus))

    def power(self, a, exponent):
        if exponent < 0:
            return self.power(self.inverse(a), -exponent)
        result = (1, 0)
        while exponent:
            if exponent & 1:
                result = self.times(result, a)
            a = self.times(a, a)
            exponent //= 2
        return result

    def total(self, terms):
        real, imag = 0, 0
        for x, y in terms:
            real += x
            imag += y
        return real % self.modulus, imag % self.modulus

    def principal_root(self, n):
        torus_order = self.modulus//self.q*(self.q+1)
        for r in range(self.q):
            a = (r, 1)
            unit = self.times(a, self.inverse(self.conjugate(a)))
            root = self.power(unit, torus_order//n)
            if self.power(root, n//2) == (self.modulus-1, 0):
                assert self.conjugate(root) == self.inverse(root)
                for t in range(1, n):
                    x, y = self.power(root, t)
                    assert ((x-1)**2+y**2) % self.q
                return root
        raise ValueError("No principal root of the requested order.")

    def fourier(self, values, root, inverse=False):
        if inverse:
            root = self.inverse(root)
        result = [self.total(self.times(a, self.power(root, i*t)) for i, a in enumerate(values))
                  for t in range(len(values))]
        if inverse:
            factor = pow(len(values), -1, self.modulus)
            result = [self.scale(a, factor) for a in result]
        return result


class WeightedTransform:
    def __init__(self, b, c, n, rows):
        self.b, self.c, self.n, self.rows = b, c, n, rows
        self.first = n//rows
        self.p = 2**b-c
        assert prime(self.p) and 0 < c < self.p and c % 2 and b % n
        progression = 8*n*c
        q = progression-1
        while not prime(q):
            q += progression
        self.positions = [(i*b+n-1)//n for i in range(n+1)]
        self.radices = [1 << (end-start) for start, end in zip(self.positions, self.positions[1:])]
        max_digit = max(self.radices)-1
        self.envelope = 2*c*n*max_digit**2 + max_digit
        self.exponent = 2
        while q**self.exponent <= 32*self.envelope:
            self.exponent += 1
        self.ring = QuadraticRing(q, self.exponent)
        mod = self.ring.modulus
        self.root = self.ring.principal_root(n)
        odd_order = (q-1)*(mod//q)//2
        base = c*pow(2, -b, mod) % mod
        assert pow(base, odd_order, mod) == 1
        beta = pow(base, pow(n, -1, odd_order), mod)
        assert pow(beta, n, mod) == base
        self.weights = [pow(2, pos, mod)*pow(beta, i, mod) % mod
                        for i, pos in enumerate(self.positions[:-1])]

    def digits(self, value):
        return [(value >> i) % radix for i, radix in zip(self.positions, self.radices)]

    def pack(self, untwisted):
        """Pack the first-frequency orbit representatives after removing row twiddles."""
        last = self.first//2
        assert untwisted[0][1] == untwisted[last][1] == 0
        assert all(untwisted[self.first-t] == self.ring.conjugate(untwisted[t])
                   for t in range(1, last))
        return [untwisted[0][0], *(word for t in range(1, last) for word in untwisted[t]),
                untwisted[last][0]]

    def unpack(self, packet):
        assert len(packet) == self.first
        middle = [(packet[i], packet[i+1]) for i in range(1, self.first-1, 2)]
        return [(packet[0], 0), *middle, (packet[-1], 0),
                *(self.ring.conjugate(a) for a in reversed(middle))]

    def row_twiddle(self, row, t):
        return self.ring.power(self.root, row*t)

    def send(self, digits):
        R = self.ring
        first_root = R.power(self.root, self.rows)
        packets = []
        weighted = [(d*w % R.modulus, 0) for d, w in zip(digits, self.weights)]
        for row in range(self.rows):
            packets.append(self.pack(R.fourier(weighted[row::self.rows], first_root)))
        assert sum(map(len, packets)) == self.n
        return packets

    def second_stage(self, packets):
        R = self.ring
        rows = [self.unpack(pkt) for pkt in packets]
        middle = []
        for t in range(self.first):
            column = [R.times(r[t], self.row_twiddle(j, t)) for j, r in enumerate(rows)]
            middle.append(R.fourier(column, R.power(self.root, self.first)))
        return middle

    def receive(self, middle):
        R = self.ring
        outputs = [[] for _ in range(self.rows)]
        for t, column in enumerate(middle):
            recovered = R.fourier(column, R.power(self.root, self.first), inverse=True)
            for j, value in enumerate(recovered):
                outputs[j].append(R.times(value, R.inverse(self.row_twiddle(j, t))))
        digits = [None]*self.n
        for row, values in enumerate(outputs):
            packet = self.pack(values)
            assert len(packet) == self.first
            recover = R.fourier(self.unpack(packet), R.power(self.root, self.rows), inverse=True)
            for j, (real, imag) in enumerate(recover):
                i = j*self.rows+row
                assert imag == 0
                reduced = real*pow(self.weights[i], -1, R.modulus) % R.modulus
                digits[i] = reduced if 2*reduced < R.modulus else reduced-R.modulus
        return digits

    def polynomial(self, a, b, subtraction, sign):
        """Evaluate sign*(a*b-subtraction) using only the two-stage transforms."""
        R = self.ring
        A, B, C = (self.second_stage(self.send(d)) for d in (a, b, subtraction))
        middle = [[R.scale(R.total([R.times(x,y), R.scale(z,-1)]),sign)
                   for x, y, z in zip(aa, bb, cc)] for aa, bb, cc in zip(A, B, C)]
        return self.receive(middle)

    def reference(self, a, b, subtraction, sign):
        result = [-x for x in subtraction]
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                wrap, dest = divmod(i+j, self.n)
                bit = self.positions[i]+self.positions[j]-self.positions[dest]-wrap*self.b
                assert bit in (0, 1)
                result[dest] += (x*y << bit) * self.c**wrap
        return [sign*x for x in result]

    def trace(self, raw, terminal):
        carry = self.c*terminal
        for digit, radix in zip(raw, self.radices):
            carry, remainder = divmod(digit+carry, radix)
            if remainder:
                return False
        return carry == terminal


def test(b, c, n, rows):
    transform = WeightedTransform(b, c, n, rows)
    zero = [0]*n
    basis = [[int(i==j) for j in range(n)] for i in range(n)]
    for a in basis:
        for v in basis:
            assert transform.polynomial(a, v, zero, 1) == transform.reference(a, v, zero, 1)
    rng = Random(120331+n)
    positive, negative, rejected = 0, 0, 0
    for _ in range(12):
        x, y = rng.randrange(1, transform.p), rng.randrange(1, transform.p)
        for perturb in (0, 1):
            z = (x*y+perturb) % transform.p
            a, v, sub = map(transform.digits, (x, y, z))
            for sign in (-1, 1):
                raw = transform.polynomial(a, v, sub, sign)
                assert raw == transform.reference(a, v, sub, sign)
                assert max(map(abs, raw)) <= transform.envelope
                folded = sum(d << i for d, i in zip(raw, transform.positions))
                assert (folded-sign*(x*y-z)) % transform.p == 0
                terminal, rem = divmod(folded, transform.p)
                assert transform.trace(raw, terminal) == (perturb == 0) == (rem == 0)
                assert not transform.trace(raw, terminal-1) and not transform.trace(raw, terminal+1)
                if perturb:
                    rejected += 1
                elif terminal > 0:
                    positive += 1
                elif terminal < 0:
                    negative += 1
    assert min(positive, negative, rejected) > 0
    print(f"p=2^{b}-{c}={transform.p}, n={n}, FFT ring q^e={transform.ring.q}^{transform.exponent}: "
          f"{n*n} basis products, {n} base words per paired boundary, "
          f"{positive} positive and {negative} negative quotient checks, {rejected} rejections")


class RingValue:
    """Let the unchanged vendor quadratic formulas operate at a single frequency."""
    def __init__(self, ring, value):
        self.ring, self.value = ring, value

    def __add__(self, other):
        extra = other.value if isinstance(other, RingValue) else (other, 0)
        return RingValue(self.ring, self.ring.total([self.value, extra]))

    def __mul__(self, other):
        if isinstance(other, RingValue):
            return RingValue(self.ring, self.ring.times(self.value, other.value))
        return RingValue(self.ring, self.ring.scale(self.value, other))

    __radd__ = __add__
    __rmul__ = __mul__


def check_chord_arithmetic(b, c, source_length):
    transform = WeightedTransform(b, c, 4, 2)
    rng = Random(272001+b)
    valid_equations = 0
    rejected_assignments = 0
    for m in (3, 4, 5):
        data, values, _ = chord.make_instance(b, source_length, c, m, rng, "general")
        assert chord.accept(data, chord.witness(data, values))
        constants = data["constants"]
        sources = {name: values[name] for name in chord.source_names(m, "special")}
        count = chord.arithmetic_count(m)
        for changed in (None, "x2", "X2", "Y2"):
            inputs = sources.copy()
            if changed is not None:
                inputs[changed] = (inputs[changed]+1) % transform.p
            truth = chord.operations(inputs, constants, m, mode="special")[:count]
            def spectra(mapping):
                return {key: transform.second_stage(transform.send(transform.digits(value)))
                        for key, value in mapping.items()}
            src, pub = spectra(inputs), spectra(constants)
            outputs = [[[] for _ in range(transform.first)] for _ in range(count)]
            for t in range(transform.first):
                for j in range(transform.rows):
                    def at_frequency(mapping):
                        return {key: RingValue(transform.ring, value[t][j])
                                for key, value in mapping.items()}
                    raw = chord.operations(at_frequency(src), at_frequency(pub), m, mode="special")
                    for i in range(count):
                        outputs[i][t].append(raw[i].value)
            passed = []
            for actual, expected in zip(outputs, truth):
                digits = transform.receive(actual)
                assert max(map(abs, digits)) <= 8*transform.envelope
                folded = sum(x << i for x, i in zip(digits, transform.positions))
                terminal = folded//transform.p
                ok = transform.trace(digits, terminal)
                assert (folded-expected) % transform.p == 0
                assert ok == (expected % transform.p == 0)
                passed.append(ok)
            if changed is None:
                assert all(passed)
                valid_equations += len(passed)
            else:
                assert not all(passed)
                rejected_assignments += 1
    print(f"p={transform.p}, actual C38 arithmetic, m=3,4,5: {valid_equations} "
          f"valid equation checks, {rejected_assignments} perturbed witnesses rejected")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    test(13, 13, 8, 2)
    test(19, 19, 16, 4)
    check_chord_arithmetic(7, 27, 16)
    check_chord_arithmetic(5, 3, 8)
