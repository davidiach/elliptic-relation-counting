"""Joint products of the exceptional cofactors in the paired-slot count model."""

from dataclasses import dataclass
from itertools import islice
from math import gcd, isqrt, prod

from flint import fmpz

from block_unit import combine_components, ordinary_counts, slot_sizes
from divisor_unit import factored_cover


@dataclass(frozen=True)
class Bracket:
    left: int
    left_den: int
    right: int
    right_den: int

    def valid(self, x, p, bound):
        return (0 <= self.left <= self.left_den <= bound
                and 0 <= self.right <= self.right_den <= bound
                and min(self.left_den, self.right_den) >= 1
                and self.left_den + self.right_den > bound
                and self.right*self.left_den - self.left*self.right_den == 1
                and p*self.left < x*self.left_den
                and x*self.right_den < p*self.right)


def farey_bracket(x, p, bound):
    """The order-bound Farey neighbours of x/p, by batched mediant steps."""
    if not (0 < x < p and 1 <= bound < p and gcd(x, p) == 1):
        raise ValueError("Expected a reduced interior fraction and a smaller order.")
    left, left_den, right, right_den = 0, 1, 1, 1
    left_error, right_error = x, p-x
    while left_den + right_den <= bound:
        if left_error > right_error:
            steps = min((left_error-1)//right_error, (bound-left_den)//right_den)
            left += steps*right
            left_den += steps*right_den
            left_error -= steps*right_error
        else:
            # Equality would put x/p at a mediant with denominator < p.
            steps = min((right_error-1)//left_error, (bound-right_den)//left_den)
            right += steps*left
            right_den += steps*left_den
            right_error -= steps*left_error
    return Bracket(left, left_den, right, right_den)


def ratio_atoms(residue, p):
    """Return (negative, A, D) with residue = sign*A/D and A,D <= floor(sqrt(p))."""
    if not 0 < residue < p:
        raise ValueError("Supply a canonical nonzero field element.")
    inverse = pow(residue, -1, p)
    x = min(residue, inverse)
    bracket = farey_bracket(x, p, isqrt(p))
    lower = x*bracket.left_den - p*bracket.left
    upper = p*bracket.right - x*bracket.right_den
    negative = upper < lower
    a, d = (upper, bracket.right_den) if negative else (lower, bracket.left_den)
    if residue != x:
        a, d = d, a
    return negative, a, d


def rough_atoms(weight, p, cover):
    magnitude = min(weight, p-weight)
    return ratio_atoms(magnitude//gcd(magnitude, cover), p)


@dataclass(frozen=True)
class AtomPin:
    side: str
    prime: int = 1
    power: int = 1

    def matches(self, atoms):
        negative, a, d = atoms
        if self.side == "sign":
            return negative
        if self.side == "numerator":
            return a % self.power == 0
        if self.side == "denominator":
            return d % self.power == 0
        raise ValueError("Unknown atom predicate.")


def prime_columns(bound):
    """Stream one public prime and all its admissible powers at a time."""
    for q in range(2, bound+1):
        if not fmpz(q).is_prime():
            continue
        power = q
        powers = []
        while power <= bound:
            powers.append(power)
            power *= q
        yield q, tuple(powers)


def atom_pins(bound):
    """Stream the public coordinate predicates with polynomial working storage."""
    yield AtomPin("sign")
    for q, powers in prime_columns(bound):
        for power in powers:
            yield AtomPin("numerator", q, power)
            yield AtomPin("denominator", q, power)


def farey_component(p, group_order, seeds, target, prefix, cover, raw_counts, batch_size=32):
    """Each vector entry represents one nonnegative guarded count before rounding."""
    if batch_size < 1:
        raise ValueError("The count batch must be nonempty.")
    queries = atom_pins(isqrt(p))
    value = 1
    while pins := tuple(islice(queries, batch_size)):
        counts = ordinary_counts(seeds, group_order, prefix, len(pins),
                                 lambda box: raw_counts(box, target, prefix, cover, pins))
        for pin, count in zip(pins, counts):
            if pin.side == "sign":
                if count % 2:
                    value = -value
            else:
                exponent = count if pin.side == "numerator" else -count
                value = value*pow(pin.prime, exponent, p) % p
    return value % p


def balanced_component(p, negative_count, interval_counts):
    """interval_counts(lo,hi) returns the exact oriented totals A,D,AA,AD,DD."""
    result = -1 if negative_count % 2 else 1
    pending = [(1, isqrt(p))]
    while pending:
        low, high = pending.pop()
        na, nd, aa, ad, dd = interval_counts(low, high)
        mass = na-nd
        energy = aa-2*ad+dd
        width = high-low+1
        if min(na, nd, aa, ad, dd) < 0 or width*energy < mass*mass:
            raise ArithmeticError("Inconsistent single and pair counts.")
        if not energy:
            continue
        if low == high:
            if energy != mass*mass:
                raise ArithmeticError("Inconsistent singleton collision count.")
            result = result*pow(low, mass, p) % p
            continue
        mid = (low+high)//2
        pending.extend([(mid+1, high), (low, mid)])
    return result % p


def exception_product(p, a, c, group_order, seeds, target, factors, raw_counts):
    """The exceptional correction alone, with all prefix continuation exponents."""
    seeds = tuple(seeds)
    sizes = slot_sizes(p, a, c, seeds, target)
    cover, _ = factored_cover(factors, p//2)
    return combine_components(p, seeds, sizes,
                              lambda part, prefix: farey_component(
                                  p, group_order, part, target, prefix, cover, raw_counts))


def exception_product_balanced(p, a, c, group_order, seeds, target, factors,
                               raw_counts, interval_counts):
    """Use a guarded sign count and oriented, unguarded single/pair collision counts."""
    seeds = tuple(seeds)
    sizes = slot_sizes(p, a, c, seeds, target)
    cover, _ = factored_cover(factors, p//2)
    def component(part, prefix):
        signs, = ordinary_counts(part, group_order, prefix, 1,
                                 lambda box: raw_counts(box, target, prefix, cover, (AtomPin("sign"),)))
        return balanced_component(
            p, signs, lambda low, high: interval_counts(part, target, prefix, cover, low, high))
    return combine_components(p, seeds, sizes, component)


def exception_product_weighted(p, a, c, seeds, target, factors, oriented_query):
    """One oriented sign query and one nonnegative witness count for each short prime.

    With empty powers the callback counts negative ratios. Otherwise one witness
    chooses a public power and either numerator divisibility or denominator
    nondivisibility. The callback returns the complete count before centering.
    """
    seeds = tuple(seeds)
    sizes = slot_sizes(p, a, c, seeds, target)
    cover, _ = factored_cover(factors, p//2)
    def component(part, prefix):
        total = prod(2*n for n in sizes[:len(part)])//(4 if prefix else 2)
        negatives = oriented_query(part, target, prefix, cover, ())
        if not 0 <= negatives <= total:
            raise ArithmeticError("Invalid oriented sign count.")
        value = -1 if negatives % 2 else 1
        for q, powers in prime_columns(isqrt(p)):
            count = oriented_query(part, target, prefix, cover, powers)
            offset = len(powers)*total
            if not 0 <= count <= 2*offset:
                raise ArithmeticError("Invalid prime valuation witness count.")
            value = value*pow(q, count-offset, p) % p
        return value % p
    return combine_components(p, seeds, sizes, component)
