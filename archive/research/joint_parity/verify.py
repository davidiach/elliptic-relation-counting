"""Finite joint norm computations on the project's small elliptic curves."""

from itertools import product
from math import comb, prod
from pathlib import Path
from random import Random
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hillclimb_001" / "vendor"))
from arithmetic_alternating_chord import point_add, points


def mat_power(matrix, exponent, p):
    result = np.eye(len(matrix), dtype=np.int64)
    while exponent:
        if exponent & 1:
            result = result @ matrix % p
        matrix = matrix @ matrix % p
        exponent //= 2
    return result


def determinant(matrix, p):
    rows = matrix.copy() % p
    value = 1
    for j in range(len(rows)):
        pivot = next((k for k in range(j, len(rows)) if rows[k, j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            rows[[j, pivot]] = rows[[pivot, j]]
            value = -value
        diagonal = int(rows[j, j])
        value = value * diagonal % p
        multiples = rows[j+1:, j] * pow(diagonal, -1, p) % p
        rows[j+1:] = (rows[j+1:] - multiples[:, None] * rows[j]) % p
    return value


def interpolation_weights(degree, p):
    """Column j contains the coefficients of the j-th Lagrange basis polynomial."""
    assert degree < p
    result = []
    for j in range(degree + 1):
        coeff = [1]
        denominator = 1
        for k in range(degree + 1):
            if k != j:
                coeff = (np.convolve(coeff, [-k, 1]) % p).tolist()
                denominator = denominator * (j-k) % p
        result.append(np.array(coeff) * pow(denominator, -1, p) % p)
    return np.array(result).T


def companion(xs, p):
    polynomial = [1]
    for x in xs:
        polynomial = (np.convolve(polynomial, [-x, 1]) % p).tolist()
    matrix = np.eye(len(xs), k=-1, dtype=np.int64)
    matrix[:, -1] = -np.array(polynomial[:-1]) % p
    return matrix


def quadratic_resultant(coeff, A, B, C, identity, p):
    """Homogeneous resultant with A*w^2 + B*w + C, including A=0."""
    degree = len(coeff) - 1
    powers_a, powers_c = [identity], [identity]
    lucas = [2 * identity, -B % p]
    ac = A @ C % p
    for k in range(1, degree + 1):
        powers_a.append(powers_a[-1] @ A % p)
        powers_c.append(powers_c[-1] @ C % p)
        if k >= 2:
            lucas.append((-B @ lucas[-1] - ac @ lucas[-2]) % p)
    result = np.zeros_like(identity)
    for i, coefficient in enumerate(coeff):
        mixed = coefficient @ powers_c[i] % p
        for j in range(i):
            mixed = (mixed + (coeff[j] @ powers_c[j] % p) @ lucas[i-j]) % p
        result = (result + (powers_a[degree-i] @ coefficient % p) @ mixed) % p
    return result


class SlotAlgebra:
    """Multiplication matrices in monomial coordinates, never a table of partial sums."""

    def __init__(self, abscissas, p, a, c):
        dims = list(map(len, abscissas))
        size = prod(dims)
        assert dims and min(dims) > 0 and p**2 * size < 2**60
        self.p = p
        self.identity = identity = np.eye(size, dtype=np.int64)
        self.coordinates = []
        for i, xs in enumerate(abscissas):
            before = np.eye(prod(dims[:i]), dtype=np.int64)
            after = np.eye(prod(dims[i+1:]), dtype=np.int64)
            self.coordinates.append(np.kron(np.kron(before, companion(xs, p)), after))
        coeff = [-self.coordinates[0] % p, identity]
        for x in self.coordinates[1:]:
            degree = 2 * (len(coeff) - 1)
            values = []
            for z in range(degree + 1):
                diff = x - z * identity
                A = diff @ diff % p
                B = -2 * ((x + z*identity) @ ((z*x + a*identity) % p) + 2*c*identity) % p
                cross = (z*x - a*identity) % p
                C = (cross @ cross - 4*c*(x + z*identity)) % p
                values.append(quadratic_resultant(coeff, A, B, C, identity, p))
            coeff = [sum(w * v for w, v in zip(row, values)) % p
                     for row in interpolation_weights(degree, p)]
        self.coeff = coeff

    def projectors(self, target_x):
        p, identity = self.p, self.identity
        prefix = identity
        odd = np.zeros_like(identity)
        for k in range(len(self.coeff)-1):
            hasse = sum((comb(i, k) * pow(target_x, i-k, p) % p) * term
                        for i, term in enumerate(self.coeff) if i >= k) % p
            prefix = prefix @ (identity - mat_power(hasse, p-1, p)) % p
            if k == 0:
                simple = prefix.copy()
            odd = (odd + (-1)**k * prefix) % p
            if not np.any(prefix):
                break
        for matrix in (odd, simple):
            assert np.array_equal(matrix @ matrix % p, matrix)
        return odd, simple

    def order_from_kernels(self, target_x):
        """Integer order of the scalar norm, from the common kernels of its jets."""
        p = self.p
        conditions = []
        order = 0
        for k in range(len(self.coeff)-1):
            hasse = sum((comb(i, k) * pow(target_x, i-k, p) % p) * term
                        for i, term in enumerate(self.coeff) if i >= k) % p
            conditions.append(hasse)
            nullity = len(self.identity) - field_rank(np.vstack(conditions), p)
            if nullity == 0:
                break
            order += nullity
        return order

    def bit(self, projector):
        norm = determinant(self.identity - 2 * projector, self.p)
        assert norm in (1, self.p-1)
        return int(norm != 1)

    def simple_norm_readout(self, target_x):
        """Use only the value and first derivatives of the scalar norm."""
        p = self.p
        value = sum(pow(target_x, i, p) * v for i, v in enumerate(self.coeff)) % p
        if determinant(value, p):
            return "empty", ()
        dz = sum((i*pow(target_x, i-1, p) % p) * self.coeff[i]
                 for i in range(1, len(self.coeff))) % p
        # Cofactor definition of a determinant derivative, valid at singular matrices.
        # The fixtures are small. Each minor here is shared by all directions.
        cofactor = np.zeros_like(value)
        for i, j in np.ndindex(value.shape):
            minor = np.delete(np.delete(value, i, axis=0), j, axis=1)
            cofactor[i, j] = (-1)**(i+j) * determinant(minor, p) % p
        if np.sum(cofactor * dz) % p == 0:
            return "multiple", ()
        unit = int(np.trace(cofactor)) % p
        assert unit
        coords = tuple(int(np.sum(cofactor * x)) * pow(unit, -1, p) % p for x in self.coordinates)
        return "singleton", coords


def children(seed):
    alpha, beta, width = seed
    yield seed
    if width > 1:
        yield from children((alpha, beta, width//2))
        yield from children((alpha, beta - width//2, width//2))


def hash_slots(curve, modulus, seeds):
    return [tuple(pt for pt in curve if (alpha*pt[0]+beta) % modulus < width)
            for alpha, beta, width in seeds]


def counts(slots, p, a):
    """Separate point-law reference, including identity and cancellation."""
    result = {}
    for row in product(*slots):
        total = None
        for point in row:
            total = point_add(total, point, p, a)
        abscissas = tuple(point[0] for point in row)
        by_abscissa = result.setdefault(total, {})
        by_abscissa[abscissas] = by_abscissa.get(abscissas, 0) + 1
    return result


def check_regular(slots, curve, p, a, c):
    expected = counts(slots, p, a)
    xs = [sorted({x for x, y in slot}) for slot in slots]
    if any(not slot for slot in xs):
        assert not expected
        return 0, 0, 0, 0
    algebra = SlotAlgebra(xs, p, a, c)
    checks = corrections = recovered = largest = 0
    for u in sorted({x for x, y in curve if y}):
        target = next(pt for pt in curve if pt[0] == u)
        rows = expected.get(target, {})
        total = sum(rows.values())
        assert algebra.order_from_kernels(u) == total
        odd, simple = algebra.projectors(u)
        assert algebra.bit(odd) == total % 2
        assert algebra.bit(simple) == len(rows) % 2
        assert int(np.trace(odd)) % p == sum(v % 2 for v in rows.values()) % p
        status, coordinates = algebra.simple_norm_readout(u)
        assert status == ("empty" if total == 0 else "singleton" if total == 1 else "multiple")
        if total == 1:
            witness, = rows
            readout = tuple(int(np.trace(x @ odd)) % p for x in algebra.coordinates)
            assert coordinates == readout == witness
            recovered += 1
        corrections += int(algebra.bit(odd) != algebra.bit(simple))
        checks += 1
        largest = max(largest, total)
    return checks, corrections, recovered, largest


def check_torsion():
    p, a, c = 7, -1, 0
    curve = points(p, a, c)
    slots = hash_slots(curve, 8, ((1, 0, 4), (3, 1, 4), (5, 3, 4)))
    options = []
    for slot in slots:
        branches = [(True, point) for point in slot if point[1] == 0]
        ordinary = sorted({x for x, y in slot if y})
        if ordinary:
            branches.append((False, ordinary))
        options.append(branches)
    expected = counts(slots, p, a)
    values = {t: 0 for t in curve if t[1]}
    for pattern in product(*options):
        offset, regular = None, []
        for fixed, value in pattern:
            if fixed:
                offset = point_add(offset, value, p, a)
            else:
                regular.append(value)
        if regular:
            algebra = SlotAlgebra(regular, p, a, c)
            for target in values:
                shifted = point_add(target, offset, p, a)
                odd, _ = algebra.projectors(shifted[0])
                values[target] ^= algebra.bit(odd)
    assert values and all(bit == sum(expected.get(t, {}).values()) % 2 for t, bit in values.items())
    return len(values)


def restriction_checks(root, curve, p, a, c):
    slots = hash_slots(curve, 16, root)
    xs = [sorted({x for x, _ in slot}) for slot in slots]
    algebra = SlotAlgebra(xs, p, a, c)
    targets = [pt for pt in curve if pt[1] < p-pt[1]]
    odd = [algebra.projectors(pt[0])[0] for pt in targets]
    checks = 0
    for box in product(*(children(s) for s in root)):
        selected = hash_slots(curve, 16, box)
        keep = algebra.identity
        for original, subset, coordinate in zip(xs, selected, algebra.coordinates):
            accepted = {x for x, _ in subset}
            local = np.zeros_like(keep)
            for x in accepted:
                lagrange = algebra.identity
                for y in original:
                    if x != y:
                        lagrange = lagrange @ (coordinate - y*algebra.identity) % p
                        lagrange = lagrange * pow(x-y, -1, p) % p
                local = (local + lagrange) % p
            keep = keep @ local % p
        expected = counts(selected, p, a)
        for target, projector in zip(targets, odd):
            assert algebra.bit(projector @ keep % p) == sum(expected.get(target, {}).values()) % 2
            checks += 1
    return checks


def field_rank(matrix, p):
    rows = matrix.copy()
    rank = 0
    for j in range(rows.shape[1]):
        indices = np.flatnonzero(rows[rank:, j])
        if indices.size:
            pivot = rank + indices[0]
            rows[[rank, pivot]] = rows[[pivot, rank]]
            multipliers = rows[rank+1:, j] * pow(int(rows[rank, j]), -1, p) % p
            rows[rank+1:] = (rows[rank+1:] - multipliers[:,None] * rows[rank]) % p
            rank += 1
    return rank


def power_ranks():
    """Inspect one powering route in evaluation coordinates, only on a finite grid."""
    p, a, c, mod = 65521, 44, 1, 65536
    rng = Random(18112026)
    squares = {y*y % p: y for y in range(1, p)}
    rhs = ((x**3 + a*x + c) % p for x in range(p))
    assert 1 + sum(2 if v in squares else int(v == 0) for v in rhs) == 65419
    widths = (8, 16, 16, 16)

    def advance(coeff, x, z):
        A = (x-z)**2 % p
        B = -2 * ((x+z) * ((x*z+a) % p) + 2*c) % p
        C = (((x*z-a) % p)**2 - 4*c*(x+z)) % p
        d = len(coeff)-1
        aa, cc, ll = [1], [1], [2, -B % p]
        for k in range(1, d+1):
            aa.append(aa[-1] * A % p)
            cc.append(cc[-1] * C % p)
            if k >= 2:
                ll.append((-B * ll[-1] - (A*C % p) * ll[-2]) % p)
        acc = 0
        for i, value in enumerate(coeff):
            mixed = value * cc[i] % p
            for j in range(i):
                mixed = (mixed + (coeff[j] * cc[j] % p) * ll[i-j]) % p
            acc = (acc + (aa[d-i]*value % p) * mixed) % p
        return acc

    for draw in range(1, 257):
        roots = []
        for width in widths:
            alpha, beta = rng.randrange(1, mod, 2), rng.randrange(mod)
            xs = [(pow(alpha, -1, mod)*(j-beta)) % mod for j in range(width)]
            roots.append([x for x in xs if x < p and (x**3+a*x+c) % p in squares])
        dims = tuple(map(len, roots))
        if min(dims) < 2:
            continue
        xx = np.ix_(*roots)
        coeff = [-xx[0], 1]
        for x in xx[1:-1]:
            degree = 2*(len(coeff)-1)
            values = [advance(coeff, x, z) for z in range(degree+1)]
            coeff = [sum(w*v for w, v in zip(row, values)) % p
                     for row in interpolation_weights(degree, p)]
        values = advance(coeff, xx[-1], 0)
        zeros = np.argwhere(values == 0)
        if len(zeros) != 1:
            continue
        candidates = [slot[j] for slot, j in zip(roots, zeros[0])]
        signed = []
        for x in candidates:
            y = squares[(x**3+a*x+c) % p]
            signed.append([(x,y), (x,-y % p)])
        if sum(counts(signed, p, a).get((0,1), {}).values()) != 1:
            continue
        flat = lambda v: v.reshape(prod(dims[:2]), -1)
        ranks = []
        power = values
        for _ in range((p-1).bit_length()):
            ranks.append(field_rank(flat(power), p))
            power = power * power % p
        # Verify the final power independently by ordinary modular exponentiation.
        final = np.array([pow(int(v), p-1, p) for v in values.flat]).reshape(dims)
        endpoint_rank = field_rank(flat(1-final), p)
        cap = min(prod(dims[:2]), prod(dims[2:]))
        assert endpoint_rank == 1 and max(ranks) == cap
        print(f"Order 65419 fixture, deterministic draw {draw}, grid {dims}, target (0,1): "
              f"ranks of successive squares {ranks}, odd projector rank {endpoint_rank}.")
        return
    raise AssertionError("No fixture with a verified singleton")


def run():
    p, a, c = 13, 4, 1
    curve = points(p, a, c)
    root = ((3, 4, 2), (7, 9, 2), (1, 8, 4))
    tests = corrections = recovered = empty = 0
    for box in product(*(children(s) for s in root)):
        slots = hash_slots(curve, 16, box)
        checks, changes, reads, _ = check_regular(slots, curve, p, a, c)
        tests += checks
        corrections += changes
        recovered += reads
        empty += int(checks == 0)
    print(f"All 63 dyadic boxes: {tests} nonempty box/target checks, {empty} empty boxes, "
          f"{recovered} singleton trace readouts, {corrections} zero-test parity corrections.")
    print(f"Root projectors reused for {restriction_checks(root, curve, p, a, c)} child/target checks.")
    # Repeated masks deliberately create multiple signed completions and counts exceeding p.
    seed = next((1, b, 2) for b in range(16) if len(hash_slots(curve, 16, [(1,b,2)])[0]) == 4)
    repeated = hash_slots(curve, 16, [seed]*4)
    checks, corrections, _, largest = check_regular(repeated, curve, p, a, c)
    assert largest >= p and corrections > 0
    print(f"Four repeated slots: {checks} target abscissas, maximum count {largest} over F_{p}, "
          f"{corrections} zero-test parity corrections.")
    print(f"Full rational two-torsion curve: {check_torsion()} non-2-torsion signed targets.")
    power_ranks()


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run()
