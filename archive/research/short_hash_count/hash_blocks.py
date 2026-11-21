"""Exact transform interfaces for a short dyadic preimage and a long public multiplier."""

from dataclasses import dataclass
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hillclimb_001"))
import common
import arithmetic_alternating_chord as arith


def ceildiv(a, b):
    return (a + b - 1) // b


@dataclass(frozen=True)
class Layout:
    b: int
    k: int
    ell: int
    nfft: int
    owners: int
    run: int

    def __post_init__(self):
        assert 0 < self.ell <= self.b and self.k > 0
        for n in (self.nfft, self.owners, self.run):
            assert n > 0 and n & (n-1) == 0
        assert self.nfft % self.owners == self.owners % self.run == 0
        assert self.payload >= self.short

    @property
    def long(self):
        return ceildiv(self.b, self.k)

    @property
    def short(self):
        return ceildiv(self.ell, self.k)

    @property
    def payload(self):
        return ((self.nfft - self.short + 1) // self.owners) * self.owners

    @property
    def blocks(self):
        return ceildiv(self.long, self.payload)


class InverseHash:
    def __init__(self, main, layout, alpha, beta):
        self.layout = l = layout
        self.radix = 1 << l.k
        self.modulus = main["modulus"]
        assert main["L"] % l.nfft == 0 and main["L2"] == l.owners
        assert main["k"] == l.k and main["b"] == l.b
        assert self.modulus > l.short * (self.radix-1)**2
        self.root = pow(main["omega"], main["L"]//l.nfft, self.modulus)
        self.v = l.owners
        self.u = l.nfft//self.v
        self.par = dict(L1=self.u, L2=self.v, L=l.nfft, omega=self.root, modulus=self.modulus)
        dyadic = 1 << l.b
        inverse = pow(alpha, -1, dyadic)
        self.constant = arith.encode((-inverse * beta) % dyadic, l.long, self.radix)
        a = arith.encode(inverse, l.long, self.radix)
        self.public_blocks = []
        for start in range(0, l.long, l.payload):
            block = a[start:start+l.payload]
            self.public_blocks.append(block + [0] * (l.nfft-len(block)))
        self.spectra = [arith.dft(block, self.root, self.modulus) for block in self.public_blocks]
        self.carry_max = l.short * (self.radix-1) + 1

    def canonical(self, word):
        return type(word) is int and 0 <= word < self.modulus

    def middle(self, k1, input_column):
        assert len(input_column) == self.v and all(map(self.canonical, input_column))
        g = pow(self.root, self.u, self.modulus)
        freq = arith.dft(input_column, g, self.modulus)
        return [arith.idft([x * spectrum[k1+self.u*k2] % self.modulus
                            for k2, x in enumerate(freq)], g, self.modulus)
                for spectrum in self.spectra]

    def raw_rows(self, row, source_row, output_rows):
        """Only this owner's H and inverse-H data are read here."""
        l = self.layout
        assert len(source_row) == self.u and all(map(self.canonical, source_row))
        recovered = arith.H_inverse_row(source_row, row, self.par)
        for j, value in enumerate(recovered):
            pos = j*self.v + row
            valid_bits = max(0, min(l.k, l.ell - pos*l.k))
            assert value < 1 << valid_bits
        assert len(output_rows) == l.blocks
        accumulated = {}
        for block, words in enumerate(output_rows):
            assert len(words) == self.u and all(map(self.canonical, words))
            for j, digit in enumerate(arith.H_inverse_row(words, row, self.par)):
                # A valid output coefficient is an ordinary nonnegative integer.
                assert digit <= l.short * (self.radix-1)**2
                pos = block*l.payload + j*self.v + row
                assert pos % self.v == row
                if pos < l.long:
                    accumulated[pos] = accumulated.get(pos, 0) + digit
        return accumulated

    def validate_owner(self, group, source_rows, output_rows, x_rows, carries):
        """x_rows are the already decoded main-C38 limbs for these residues."""
        l = self.layout
        raw = {}
        for row in range(group*l.run, (group+1)*l.run):
            raw.update(self.raw_rows(row, source_rows[row], [out[row] for out in output_rows]))
        for start in range(group*l.run, l.long, self.v):
            end = min(start+l.run, l.long)
            carry = carries[start] if start else 0
            assert type(carry) is int and 0 <= carry <= self.carry_max
            for pos in range(start, end):
                carry, digit = divmod(raw[pos]+self.constant[pos]+carry, self.radix)
                bits = min(l.k, l.b-pos*l.k)
                assert digit % (1 << bits) == x_rows[pos]
                assert 0 <= carry <= self.carry_max
            if end < l.long:
                assert carry == carries[end]

    def witness(self, short_index):
        """Construct test data. Carries here use independent schoolbook convolution."""
        l = self.layout
        digits = arith.encode(short_index, l.short, self.radix)
        source = arith.H(digits + [0] * (l.nfft-l.short), self.par)
        outputs = [[[0]*self.u for _ in range(self.v)] for _ in self.public_blocks]
        for k1 in range(self.u):
            cols = self.middle(k1, [row[k1] for row in source])
            for destination, column in zip(outputs, cols):
                for row, value in enumerate(column):
                    destination[row][k1] = value
        carry, carries = 0, {}
        for pos in range(l.long):
            if pos and pos % l.run == 0:
                carries[pos] = carry
            total = self.constant[pos] + carry
            for j, digit in enumerate(digits):
                if pos >= j:
                    block, within = divmod(pos-j, l.payload)
                    total += self.public_blocks[block][within]*digit
            carry = total // self.radix
        return source, outputs, carries

    def check(self, ports, x_rows):
        source, outputs, carries = ports
        for k1 in range(self.u):
            expected = self.middle(k1, [row[k1] for row in source])
            actual = [[row[k1] for row in out] for out in outputs]
            assert actual == expected
        for group in range(self.v//self.layout.run):
            self.validate_owner(group, source, outputs, x_rows, carries)
