"""Deterministic base-field points and streamed row-rank certificates.

The caller supplies a prime field and, for sparse tests, a certified primitive
root. The large-input construction of that root is LZR Proposition 3.3.
"""
from math import comb


def explicit_row_points(prime, radix, coordinates, row_count):
    if not 0 < row_count <= prime or radix < 2 or coordinates < 1:
        raise ValueError("invalid interpolation parameters")
    for z in range(row_count):
        point = [z]
        for _ in range(1, coordinates):
            point.append(pow(point[-1], radix, prime))
        yield tuple(point)


def sparse_witness_points(prime, generator, coordinates, terms, max_degree):
    if min(coordinates, terms) < 1 or max_degree < 0:
        raise ValueError("invalid sparsity parameters")
    slopes = 1 + (coordinates - 1) * comb(terms, 2)
    degrees = max_degree * sum(slopes**j for j in range(coordinates))
    if degrees >= prime - 1:
        raise ValueError("projected exponents may collide modulo the group order")
    for t in range(1, slopes + 1):
        step = [pow(generator, t**j, prime) for j in range(coordinates)]
        point = [1] * coordinates
        for _ in range(terms):
            yield tuple(point)
            point = [v * factor % prime for v, factor in zip(point, step)]


class RowSpace:
    def __init__(self, prime, columns):
        self.prime = prime
        self.columns = columns
        self.pivots = {}

    def add(self, row):
        if len(row) != self.columns:
            raise ValueError("row width differs from certificate width")
        v = [x % self.prime for x in row]
        for j in range(self.columns):
            if not v[j]:
                continue
            pivot = self.pivots.get(j)
            if pivot is None:
                inv = pow(v[j], -1, self.prime)
                self.pivots[j] = [x * inv % self.prime for x in v]
                return True
            multiplier = v[j]
            v = [(x - multiplier*y) % self.prime for x, y in zip(v, pivot)]
        return False

    @property
    def rank(self):
        return len(self.pivots)


def compress_stack(points, matrix_at, prime, columns):
    """Return attaining points (at most columns) and the certified row space."""
    selected = []
    space = RowSpace(prime, columns)
    for x in points:
        previous = space.rank
        for row in matrix_at(x):
            space.add(row)
        if space.rank > previous:
            selected.append(x)
        if space.rank == columns:
            break
    return selected, space
