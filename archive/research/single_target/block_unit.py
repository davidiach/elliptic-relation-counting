"""Windowed certificates and a guarded-count reduction for a single field unit."""

from dataclasses import dataclass
from math import prod
from typing import Callable


@dataclass(frozen=True)
class PowerBlocks:
    p: int
    digit_bits: int
    bases: tuple[int, ...]

    @classmethod
    def prepare(cls, p, base, digit_bits):
        """The caller supplies an odd prime and a base of certified order p-1."""
        if p < 3 or p % 2 == 0 or not 1 <= base < p or digit_bits < 1:
            raise ValueError("Invalid field, base or digit width.")
        width = (p-2).bit_length()
        bases = []
        for position in range(0, width, digit_bits):
            bases.append(base)
            if position + digit_bits < width:
                base = pow(base, 1 << digit_bits, p)
        return cls(p, digit_bits, tuple(bases))

    @property
    def bit_count(self):
        return (self.p-2).bit_length()

    def widths(self):
        return tuple(min(self.digit_bits, self.bit_count-start)
                     for start in range(0, self.bit_count, self.digit_bits))

    def digits(self, exponent):
        return tuple((exponent >> (i*self.digit_bits)) & ((1 << width)-1)
                     for i, width in enumerate(self.widths()))

    def matches_slices(self, block, digit, residue, slice_bits):
        """Each output slice is a separate local predicate on a copy of the digit."""
        if slice_bits < 1 or not 0 <= residue < self.p:
            return False
        if not 0 <= digit < 1 << self.widths()[block]:
            return False
        for start in range(0, self.p.bit_length(), slice_bits):
            mask = (1 << min(slice_bits, self.p.bit_length()-start))-1
            expected = (pow(self.bases[block], digit, self.p) >> start) & mask
            if ((residue >> start) & mask) != expected:
                return False
        return True

    def accepts(self, exponent, weight):
        """Scalar simulation with all uniquely forced internal words reconstructed."""
        if not 0 <= exponent < self.p-1 or not 0 < weight < self.p:
            return False
        residue = 1
        for base, digit in zip(self.bases, self.digits(exponent)):
            residue = residue*pow(base, digit, self.p) % self.p
        return residue == weight


def split_for_rounding(seeds, group_order, slack):
    """Iteratively partition a dyadic box until its guarded counts round exactly."""
    if slack*(1 << len(seeds)) >= group_order:
        raise ValueError("The finite offset rounding bound fails even at width one.")
    pending = [tuple(seeds)]
    while pending:
        box = pending.pop()
        if slack*prod(2*width for _, _, width in box) < group_order:
            yield box
            continue
        i = max(range(len(box)), key=lambda j: box[j][2])
        alpha, beta, width = box[i]
        half = width//2
        for offset in (half, 0):
            pending.append((*box[:i], (alpha, beta-offset, half), *box[i+1:]))


def slot_size(p, a, c, modulus, seed):
    alpha, beta, width = seed
    step = pow(alpha, -1, modulus)
    count = 0
    for i in range(width):
        x = (step*(i-beta)) % modulus
        if x < p:
            rhs = (x**3 + a*x + c) % p
            if rhs == 0:
                raise ValueError("The paired-slot model excludes rational two-torsion.")
            count += pow(rhs, (p-1)//2, p) == 1
    return count


def slot_sizes(p, a, c, seeds, target):
    modulus = 1 << p.bit_length()
    if not seeds or any(not (1 <= w <= modulus) or w & (w-1) for _, _, w in seeds):
        raise ValueError("Supply a nonempty sequence of dyadic windows.")
    x, y = target
    if not (0 <= x < p and 0 < y < p and (y*y-x**3-a*x-c) % p == 0):
        raise ValueError("Supply a canonical finite non-two-torsion target.")
    return [slot_size(p, a, c, modulus, seed) for seed in seeds]


def ordinary_counts(seeds, group_order, prefix, number, raw_counts):
    counts = [0]*number
    for box in split_for_rounding(seeds, group_order, 2*len(seeds)+11):
        raw = tuple(raw_counts(box))
        if len(raw) != number or any(z < 0 for z in raw):
            raise ValueError("Invalid guarded count vector.")
        for k, z in enumerate(raw):
            counts[k] += (z + group_order-1)//group_order
    symmetry = 4 if prefix else 2
    if any(n % symmetry for n in counts):
        raise ArithmeticError("The exact component counts must respect sign symmetry.")
    return [n//symmetry for n in counts]


def combine_components(p, seeds, sizes, component):
    if 0 in sizes:
        return 1
    unit = component(seeds, False)
    tail = 1
    for j in range(len(seeds), 1, -1):
        unit = unit*pow(component(seeds[:j], True), tail, p) % p
        tail = tail*2*sizes[j-1] % (p-1)
    return unit


def unit_from_counts(power: PowerBlocks, a, c, group_order, seeds, target,
                     raw_bit_counts: Callable):
    """Reduce the unit to exact raw guarded counts, one bit pin per component.

    raw_bit_counts(box, target, prefix, power) returns one raw count per bit,
    summed over the separate affine and infinity paths. Each raw count includes
    the allowed offsets and both signs of each slot. The callback is the
    contractor interface. The small-field verifier supplies an exhaustive oracle.
    """
    p = power.p
    seeds = tuple(seeds)
    sizes = slot_sizes(p, a, c, seeds, target)

    def component(box_seeds, prefix):
        counts = ordinary_counts(box_seeds, group_order, prefix, power.bit_count,
                                 lambda box: raw_bit_counts(box, target, prefix, power))
        exponent = sum(n << k for k, n in enumerate(counts)) % (p-1)
        return pow(power.bases[0], exponent, p)

    return combine_components(p, seeds, sizes, component)
