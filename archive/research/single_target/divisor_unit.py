"""Recover a field unit with integer divisor counts and streamed exceptional values."""

from dataclasses import dataclass
from math import gcd

from flint import fmpz

from block_unit import combine_components, ordinary_counts, slot_sizes


def public_remainder(value, divisor, word_bits):
    """Simulate the unique quotient/remainder ports of the public-limb check."""
    remainder = 0
    mask = (1 << word_bits)-1
    start = ((value.bit_length()-1)//word_bits)*word_bits
    for shift in range(start, -1, -word_bits):
        word = (value >> shift) & mask
        quotient, remainder = divmod((remainder << word_bits) + word, divisor)
        assert 0 <= quotient <= mask and 0 <= remainder < divisor
    return remainder


@dataclass(frozen=True)
class DivisorPin:
    kind: str
    low: int = 0
    high: int = 0

    def matches(self, weight, p, cover):
        magnitude = min(weight, p-weight)
        if self.kind == "negative":
            return weight > p//2
        if self.kind == "multiple":
            return magnitude % self.low == 0
        if self.kind == "exception":
            return (self.low <= magnitude <= self.high
                    and public_remainder(cover, magnitude, p.bit_length()) != 0)
        raise ValueError("Unknown divisor predicate.")


def factored_cover(factors, maximum):
    """Ignore prime powers above the largest possible centered magnitude."""
    seen = set()
    cover = 1
    powers = []
    for prime, exponent in factors:
        if prime in seen or exponent < 1 or not fmpz(prime).is_prime():
            raise ValueError("Supply a distinct prime with a positive exponent in each pair.")
        seen.add(prime)
        power = 1
        for _ in range(exponent):
            power *= prime
            if power > maximum:
                break
            cover *= prime
            powers.append((prime, power))
    return cover, powers


def unit_from_divisors(p, a, c, group_order, seeds, target, factors, raw_divisor_counts):
    """The raw-count callback has the same guarded semantics as the block reducer.

    Its last two arguments are a public integer and a sequence of DivisorPin
    predicates. Each returned entry represents one ordinary graph query before
    offset rounding. Only a single search path and its pending siblings are kept.
    """
    seeds = tuple(seeds)
    sizes = slot_sizes(p, a, c, seeds, target)
    maximum = p//2
    cover, powers = factored_cover(factors, maximum)

    def component(box_seeds, prefix):
        def counts(pins):
            return ordinary_counts(box_seeds, group_order, prefix, len(pins),
                                   lambda box: raw_divisor_counts(box, target, prefix, cover, pins))

        initial = [DivisorPin("negative"), DivisorPin("exception", 1, maximum)]
        initial.extend(DivisorPin("multiple", power) for _, power in powers)
        negatives, outside, *valuations = counts(initial)
        value = -1 if negatives % 2 else 1
        for (prime, _), tally in zip(powers, valuations):
            value = value*pow(prime, tally, p) % p

        pending = [(1, maximum, outside)] if outside else []
        while pending:
            low, high, count = pending.pop()
            if low == high:
                rough = low//gcd(low, cover)
                value = value*pow(rough, count, p) % p
                continue
            mid = (low+high)//2
            left, = counts([DivisorPin("exception", low, mid)])
            if not 0 <= left <= count:
                raise ArithmeticError("Child multiplicity exceeds its parent.")
            if left < count:
                pending.append((mid+1, high, count-left))
            if left:
                pending.append((low, mid, left))
        return value % p

    return combine_components(p, seeds, sizes, component)
