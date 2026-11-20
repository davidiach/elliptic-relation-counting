"""Independent finite-field rank checks using the binary packet kernels."""

from fractions import Fraction
from itertools import product
from math import isqrt


PRIME = 1_000_000_007
import packet_space


def normalized_product(rows, sources, columns):
    """Top Taylor terms after stripping the common grid polynomial at infinity."""
    values = [1] + [0] * (columns - 1)
    for x in rows:
        summed = [0] * columns
        for y in sources:
            # Divide by 1 - (y-x) z as a truncated formal series.
            previous = 0
            for j, coefficient in enumerate(values):
                previous = (coefficient + (y - x) * previous) % PRIME
                summed[j] += previous
        values = [value % PRIME for value in summed]
    return values


def marginal_inputs(width, order):
    # Every interval of order consecutive grid differences occurs.
    for shift in range(order - width, width):
        start = max(0, -shift)
        yield range(start, start + order), (max(0, shift),)
    if order == 1:
        return
    slack = int(order < width)
    if order % 2 and slack:
        yield [0, *range(2, order + 1)], (0, order)
    for gap in reversed(range(1, order, 2)):
        for start in range(slack + 1):
            for y in range(order + slack - gap):
                yield range(start, start + order), (y, y + gap)


def check_layer(width, order):
    columns = 2 * width - 1
    basis = {}
    examined = 0
    for rows, sources in marginal_inputs(width, order):
        assert len(rows) == order and all(0 <= x < width for x in rows)
        assert all(0 <= y < width for y in sources)
        examined += 1
        row = normalized_product(rows, sources, columns)
        for i in range(columns):
            coefficient = row[i]
            if coefficient == 0:
                continue
            if i not in basis:
                inverse = pow(coefficient, -1, PRIME)
                basis[i] = [value * inverse % PRIME for value in row[i:]]
                break
            for j, entry in enumerate(basis[i], i):
                row[j] = (row[j] - coefficient * entry) % PRIME
        if len(basis) == columns:
            return examined
    raise AssertionError((width, order, len(basis), columns))


def check_literal_marginals():
    """Compare to a direct sum of the literal indicator polynomials from checkpoint 009."""
    tested = 0
    for exponent in (1, 2):
        width = 1 << exponent
        kernels, derivative = packet_space.kernels(exponent)
        grid = [1]
        for node in range(-width, width):
            grid = packet_space.multiply(grid, [-node, 1])
        for bits in product((0, 1), repeat=width):
            sources = [y for y, bit in enumerate(bits) if bit]
            if not sources:
                continue
            indicators = packet_space.family(exponent, bits, kernels, derivative)
            for selection in product((0, 1), repeat=width):
                rows = [x for x, bit in enumerate(selection) if bit]
                if not rows:
                    continue
                marginal = [0]
                for label, polynomial in indicators.items():
                    if all(label[x] for x in rows):
                        marginal = packet_space.add(marginal, polynomial)
                order = len(rows)
                degree = (2 * width - 1) * order
                common = [derivative ** (width - order)]
                for _ in rows:
                    common = packet_space.multiply(common, grid)
                leading = normalized_product(rows, sources, 2 * width - 1)
                for j in range(len(leading)):
                    reconstructed = sum(
                        leading[t] * common[-1 - (j - t)] for t in range(j + 1)
                    )
                    assert (reconstructed - marginal[degree - j]) % PRIME == 0
                tested += 1
    return tested


def endpoint_weights(width, overlap):
    weights = [Fraction(1)]
    for j in range(overlap - 1):
        numerator = (width - overlap + 2*j + 2) * (width - 1 - j) * (overlap - 1 - j)
        denominator = (width + overlap - 2*j - 2) * (j + 1) * (width - overlap + 1 + j)
        weights.append(weights[-1] * Fraction(numerator, denominator))
    return weights


def check_endpoint_identity(width):
    a = endpoint_weights(width, 3)[1]
    b, c = endpoint_weights(width, 5)[1:3]
    d, e, f = endpoint_weights(width, 7)[1:4]
    multiplier = (d - b) * (a*b - c + 1) - (e - c - 1)*a + f - d - b
    reduced = Fraction(
        384 * width * (width - 5) * (width - 1),
        (width - 4) * (width - 2) * (width + 1)**2 * (width + 3)**2 * (width + 5),
    )
    assert multiplier == reduced and multiplier > 0


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    assert all(PRIME % p for p in range(2, isqrt(PRIME) + 1))
    print(f"Literal indicator marginal comparisons: {check_literal_marginals()}", flush=True)
    for width in range(8, 130, 2):
        check_endpoint_identity(width)
    print("Endpoint rational identity checked at 61 even widths", flush=True)
    for width in (2, 4, 8, 16, 32, 64):
        counts = [check_layer(width, k) for k in range(1, width + 1)]
        print(f"L={width}: full dimension {1 + width*(2*width-1)}, "
              f"{sum(counts)} marginal rows examined", flush=True)
    count = check_layer(128, 128)
    print(f"L=128: highest layer rank 255, {count} marginal rows examined", flush=True)
