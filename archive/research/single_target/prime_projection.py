"""Prime-free projections of the exceptional multiplicative correction."""

from itertools import product
from math import gcd, isqrt, prod

from block_unit import combine_components, slot_sizes
from divisor_unit import factored_cover
from target import chord_scale, lift_windows, point_add


def stream_joint_count(p, a, c, seeds, target, cover, exponent_witness):
    """Finite evaluation of the full-box tagged count, charging each supplied log evaluation.

    A canonical exponent e gives exactly e accepted rank choices. This evaluator
    streams the oriented tuples. The graph reduction counts their assignments.
    """
    if not seeds:
        raise ValueError("Supply at least one slot.")
    slots = lift_windows(p, a, c, 1 << p.bit_length(), seeds)
    choices = [slots[0]]
    for slot in slots[1:]:
        choices.append([signed for x, y in slot for signed in ((x, y), (x, p-y))])

    def rank_count(weight):
        if cover is None:
            return exponent_witness(weight)
        magnitude = min(weight, p-weight)
        return exponent_witness(magnitude//gcd(magnitude, cover))

    result = 0
    for row in product(*choices):
        partial = row[0]
        for point in row[1:]:
            if point[1] < p-point[1]:
                result += rank_count(chord_scale(partial, point, p))
            partial = point_add(partial, point, p, a)
        weight = p-1 if partial is None else (target[0]-partial[0]) % p or 1
        result += rank_count(weight)
    return result


def scalar_product(power, seeds, target, factors, ordinary_joint_count):
    """One ordinary count over the full hash box, component tag, exponent and rank.

    With factors=None use the full local field weight, including its sign.
    With a factored integer use the exceptional cofactor. The callback uses
    the supplied primitive base. The disjoint terminal and
    chord branches and all sign conventions are specified in prime_projection.tex.
    """
    seeds = tuple(seeds)
    if not seeds:
        raise ValueError("Supply at least one slot.")
    cover = None if factors is None else factored_cover(factors, power.p//2)[0]
    total = ordinary_joint_count(seeds, target, cover, power)
    if not isinstance(total, int) or total < 0:
        raise ValueError("Expected a nonnegative exact integer count.")
    return pow(power.bases[0], total % (power.p-1), power.p)


def height_bits(bound, population, energy):
    """Upper bound on both log heights after numerator/denominator cancellations."""
    if bound < 1 or population < 0 or energy < 0:
        raise ValueError("Invalid atom bound, population or collision energy.")
    # The positive and negative masses are equal since each occurrence supplies both.
    mass = min(population, energy//2, isqrt(bound*energy)//2)
    return mass*(bound-1).bit_length()


def reconstruct_positive(residue, modulus, height):
    """Unique positive reduced fraction with both entries at most height."""
    if modulus <= 2*height*height or not 0 < residue < modulus:
        raise ValueError("Insufficient precision or noncanonical unit.")
    prev, rem = modulus, residue
    prev_den, den = 0, 1
    while rem > height:
        quotient, next_rem = divmod(prev, rem)
        prev, rem = rem, next_rem
        prev_den, den = den, prev_den-quotient*den
    if not (0 < rem <= height and 0 < den <= height and gcd(rem, den) == 1
            and (rem-residue*den) % modulus == 0):
        raise ArithmeticError("The counts contradict the supplied height bound.")
    return rem, den


def auxiliary_primes():
    """Small odd moduli for reconstruction, unrelated to the atom factor base."""
    from flint import fmpz
    candidate = 3
    while True:
        if fmpz(candidate).is_prime():
            yield candidate
        candidate += 2


def unit_log_table(prime):
    """Return a generator and its full small-field logarithm table, with cost charged."""
    from flint import fmpz
    order_factors = [int(q) for q, _ in fmpz(prime-1).factor()]
    base = next(g for g in range(2, prime)
                if all(pow(g, (prime-1)//q, prime) != 1 for q in order_factors))
    logs = [0]*prime
    value = 1
    for e in range(prime-1):
        logs[value] = e
        value = value*base % prime
    return base, logs


def auxiliary_component(p, population, negative_count, energy, image_counts):
    """Recover the positive rational atom product and then apply the Farey sign.

    image_counts(ell, logs) returns signed valuation and signed stripped-unit
    logarithm sums. Each difference is obtained from ordinary additive counts.
    """
    if not 0 <= negative_count <= population:
        raise ValueError("Invalid oriented sign count.")
    h = height_bits(isqrt(p), population, energy)
    sign = -1 if negative_count % 2 else 1
    if h == 0:
        return sign % p
    height = 1 << h
    threshold = 2*height*height
    modulus, residue = 1, 0
    for ell in auxiliary_primes():
        if ell == p:  # keep the reconstruction independent of the output field
            continue
        base, logs = unit_log_table(ell)
        valuation, log_sum = image_counts(ell, logs)
        if valuation:
            continue
        image = pow(base, log_sum % (ell-1), ell)
        lift = (image-residue)*pow(modulus, -1, ell) % ell
        residue += modulus*lift
        modulus *= ell
        if modulus > threshold:
            num, den = reconstruct_positive(residue, modulus, height)
            return sign*num*pow(den, -1, p) % p


def auxiliary_product(p, a, c, seeds, target, factors, summary_counts, image_counts):
    """summary_counts returns the oriented negative count and integer collision energy."""
    seeds = tuple(seeds)
    sizes = slot_sizes(p, a, c, seeds, target)
    cover, _ = factored_cover(factors, p//2)
    def component(part, prefix):
        population = prod(2*n for n in sizes[:len(part)])//(4 if prefix else 2)
        negatives, energy = summary_counts(part, target, prefix, cover)
        return auxiliary_component(
            p, population, negatives, energy,
            lambda ell, logs: image_counts(part, target, prefix, cover, ell, logs))
    return combine_components(p, seeds, sizes, component)
