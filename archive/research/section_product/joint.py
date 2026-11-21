"""Exact section products on odd cyclic finite test groups, using integer convolutions."""

from math import gcd, prod
from pathlib import Path
import sys

from flint import fmpz_poly

research = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(research / "implicit_resultant"))
sys.path.insert(0, str(research / "scalar_norm"))
from hash_probe import batch_remainders
from scalar import form, polynomial, product_tree


def chord_scale(P, Q, p):
    if P is None or Q is None:
        return 1
    if P[0] == Q[0]:
        return 4*P[1]**2 % p
    return (P[0]-Q[0])**2 % p


class CurveTable:
    """The whole finite test group is charged, including its point index."""

    def __init__(self, p, a, c, r, generator, add):
        if r < 3 or r % 2 == 0 or gcd(r, p-1) != 1:
            raise ValueError("The group order must be odd and invertible modulo p-1.")
        self.p, self.a, self.c, self.r = p, a, c, r
        self.points = [None]
        current = generator
        for _ in range(1, r):
            if current is None:
                raise ValueError("Generator has smaller order than supplied.")
            self.points.append(current)
            current = add(current, generator, p, a)
        if current is not None:
            raise ValueError("Generator does not have the supplied order.")
        self.index = {self.points[k][0]: k for k in range(1, (r+1)//2)}

        # Solve the adjacent-scale recurrence, enforcing gamma(rP)=gamma(O)=1.
        unscaled = [1, 1]
        for point in self.points[1:]:
            denom = unscaled[-2]*chord_scale(point, generator, p)
            unscaled.append(unscaled[-1]**2 * pow(denom, -1, p) % p)
        first = pow(unscaled[r], -pow(r*r, -1, p-1), p)
        self.gamma = [v*pow(first, k*k % (p-1), p) % p for k, v in enumerate(unscaled[:r])]

    def windows(self, modulus, seeds):
        slots = []
        for alpha, beta, width in seeds:
            if modulus < self.p or modulus & (modulus-1) or width < 1 or width > modulus or width & (width-1):
                raise ValueError("Expected a dyadic window in a power-of-two hash modulus.")
            step = pow(alpha, -1, modulus)
            xs = (step*(j-beta) % modulus for j in range(width))
            slots.append([x for x in xs if x in self.index])
        return slots


def cyclic(f, r):
    # Each multiplication has degree at most 2r-2.
    return f.truncate(r) + f.right_shift(r)


def slot_moments(table, slots):
    """Count and first two exact integer moments of the encoded abscissa tuple."""
    r, p = table.r, table.p
    count, first, second = fmpz_poly([1]), fmpz_poly(), fmpz_poly()
    place = 1
    for i, slot in enumerate(slots):
        arrays = [[0]*r for _ in range(3)]
        for x in slot:
            j = table.index[x]
            digit = place*x
            for k in ([j, r-j] if i else [j]):
                arrays[0][k], arrays[1][k], arrays[2][k] = 1, digit, digit*digit
        base, linear, square = map(fmpz_poly, arrays)
        second = cyclic(second*base + 2*first*linear + count*square, r)
        first = cyclic(first*base + count*linear, r)
        count = cyclic(count*base, r)
        place *= p
    return count, first, second


def evaluations(f, moduli):
    return [int(g[0]) for g in batch_remainders(f, moduli)]


class SeedProduct:
    """Normalized univariate norm as one scalar times a divisor of elliptic sums."""

    def __init__(self, table, slots):
        if not slots:
            raise ValueError("Supply at least one slot.")
        self.table, self.arity = table, len(slots)
        p, r = table.p, table.r
        self.moments = {}
        self.scale = 1
        if not all(slots):
            return
        volume = prod(map(len, slots))
        for slot in slots:
            exponent = (1 << (self.arity-1))*(volume//len(slot))
            local = 1
            for x in slot:
                local = local*table.gamma[table.index[x]] % p
            self.scale = self.scale*pow(local, exponent, p) % p
        count, first, second = slot_moments(table, slots)
        if int(count[0]) % 2:
            self.scale = -self.scale % p
        for x, j in table.index.items():
            h = int(count[j]+count[r-j])
            if h:
                self.moments[x] = h, int(first[j]+first[r-j]), int(second[j]+second[r-j])
                self.scale = self.scale*pow(table.gamma[j], -h, p) % p

    def target_units(self, targets):
        """For each t return (vanishing order, first nonzero Taylor coefficient, N_z)."""
        targets = list(targets)
        if not targets:
            return []
        p = self.table.p
        linear = polynomial([0, 1], p)
        moduli = [linear-t for t in targets]
        orders = [self.moments.get(t, (0, 0, 0))[0] for t in targets]
        units = [self.scale]*len(targets)
        log_derivatives = [0]*len(targets)
        max_order = max((h for h, _, _ in self.moments.values()), default=0)
        for bit in range(max_order.bit_length()):
            weight = 1 << bit
            leaves = [(linear-x, []) for x, (h, _, _) in self.moments.items() if h & weight]
            if not leaves:
                continue
            f, _ = product_tree(leaves)
            values = evaluations(f, moduli)
            derivatives = evaluations(f.derivative(), moduli)
            for j, (v, d) in enumerate(zip(values, derivatives)):
                units[j] = units[j]*pow(v or d, weight, p) % p
                if orders[j] == 0:
                    log_derivatives[j] += weight*d*pow(v, -1, p)
        return [(h, unit, unit*ld % p if h == 0 else unit if h == 1 else 0)
                for h, unit, ld in zip(orders, units, log_derivatives)]

    def zero_gradient(self, target, order, unit, dz):
        """All additive derivatives at a zero, including one tuple with repeated completions."""
        if not order:
            raise ValueError("The target must be a zero.")
        h, first, second = self.moments[target]
        if h != order:
            raise ValueError("Target order mismatch.")
        return gradient_from_moments(self.table.p, self.table.a, self.table.c, self.arity,
                                     target, h, first, second, unit, dz)


def gradient_from_moments(p, a, c, arity, target, h, first, second, unit, dz):
    """Read all additive derivatives at a zero from its local unit and exact moments."""
    if h <= 0:
        raise ValueError("The target must be a zero.")
    if h*second != first*first:
        return [0]*(arity+2)
    code = first//h
    row = []
    for _ in range(arity):
        code, x = divmod(code, p)
        row.append(x)
    f = form(row, p, a, c)
    linear = polynomial([-target, 1], p)
    quotient, remainder = divmod(f, linear**h)
    if remainder or not quotient(target):
        raise ArithmeticError("The isolated tuple has the wrong root multiplicity.")
    cofactor = unit*pow(int(quotient(target)), -1, p) % p
    return [dz, cofactor, *(x*cofactor % p for x in row)]
