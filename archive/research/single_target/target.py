"""Single-target local unit and integer moments on paired abscissa slots."""

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import sys

from flint import fmpz

research = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(research / "section_product"))
sys.path.insert(0, str(research / "hillclimb_001" / "vendor"))
from arithmetic_alternating_chord import point_add
from joint import batch_remainders, chord_scale, gradient_from_moments, polynomial, product_tree


@dataclass(frozen=True)
class LocalValue:
    order: int
    unit: int
    target_derivative: int
    first_moment: int
    second_moment: int

    def zero_gradient(self, p, a, c, arity, t):
        return gradient_from_moments(p, a, c, arity, t, self.order, self.first_moment,
                                     self.second_moment, self.unit, self.target_derivative)


def lift_windows(p, a, c, modulus, seeds):
    """One rational lift per accepted abscissa, directly from its inverse window."""
    if modulus < p or modulus & (modulus-1):
        raise ValueError("Expected a power-of-two modulus covering the base field.")
    slots = []
    for alpha, beta, width in seeds:
        if width < 1 or width > modulus or width & (width-1):
            raise ValueError("Expected a positive dyadic width within the hash modulus.")
        step = pow(alpha, -1, modulus)
        slot = []
        for j in range(width):
            x = step*(j-beta) % modulus
            if x >= p:
                continue
            rhs = (x**3+a*x+c) % p
            if not rhs:
                raise ValueError("This evaluator uses paired slots without two-torsion.")
            if pow(rhs, (p-1)//2, p) == 1:
                y = int(fmpz(rhs).sqrtmod(p))
                slot.append((x, min(y, p-y)))
        slots.append(slot)
    return slots


def extend_half(records, slots, positions, p, a, target_y=None):
    """Propagate signed roots and the common projective scalar of their form."""
    scale, log_slope = 1, 0
    inverse_y = pow(target_y, -1, p) if target_y is not None else None
    for i in positions:
        slot = slots[i]
        place = p**i
        extended = []
        scale = pow(scale, 2*len(slot), p)
        log_slope *= 2*len(slot)
        for point, code in records:
            for other in slot:
                x, y = other
                scale = scale*chord_scale(point, other, p) % p
                if inverse_y is not None and point is not None:
                    u, v = point
                    if u == x:
                        log_slope += (3*u*u+a)*inverse_y*pow(2*v, -1, p)
                    else:
                        log_slope += 2*v*inverse_y*pow(u-x, -1, p)
                next_code = code+place*x
                extended.append((point_add(point, other, p, a), next_code))
                extended.append((point_add(point, (x, -y % p), p, a), next_code))
        records = extended
        log_slope %= p
    return records, scale, log_slope


def local_value(p, a, target, slots):
    """Return the unit, order, target derivative and exact moments for one target."""
    if not slots or target is None or not target[1]:
        raise ValueError("Supply at least one slot and a finite non-two-torsion target.")
    if not all(slots):
        return LocalValue(0, 1, 0, 0, 0)
    cut = max(1, len(slots)//2)
    left, left_scale, _ = extend_half([(P, P[0]) for P in slots[0]], slots,
                                      range(1, cut), p, a)
    right, right_scale, log_scale = extend_half([(target, 0)], slots,
                                                reversed(range(cut, len(slots))), p, a, target[1])
    buckets = defaultdict(lambda: [0, 0, 0])
    for point, code in left:
        stats = buckets[None if point is None else point[0]]
        stats[0] += 1
        stats[1] += code
        stats[2] += code*code
    empty = (0, 0, 0)
    infinity = buckets.get(None, empty)[0]
    x = polynomial([0, 1], p)
    leaves = [((x-u)**s[0], []) for u, s in buckets.items() if u is not None]
    monic = product_tree(leaves)[0] if leaves else polynomial([1], p)
    finite = sorted({point[0] for point, _ in right if point is not None})
    moduli = [(x-u)**max(2, buckets.get(u, empty)[0]+1) for u in finite]
    local = {}
    if moduli:
        for u, remainder in zip(finite, batch_remainders(monic, moduli)):
            k = buckets.get(u, empty)[0]
            local[u] = (int(remainder[k] if k else remainder[0]+u*remainder[1]) % p,
                        int(remainder[1]) if not k else 0)

    # The affine leading coefficient includes one minus sign per infinite left root.
    unit = pow(right_scale, len(left), p)*pow((-1)**(infinity % 2)*left_scale, len(right), p) % p
    log_slope = len(left)*log_scale
    inv_y = pow(target[1], -1, p)
    infinity_lead = pow(inv_y*pow(2, -1, p), 2, p)
    order = first = second = 0
    for point, right_code in right:
        u = None if point is None else point[0]
        k, first_left, second_left = buckets.get(u, empty)
        copies = 2 if point is None else 1
        order += copies*k
        first += copies*(first_left+k*right_code)
        second += copies*(second_left+2*right_code*first_left+k*right_code**2)
        if point is None:
            factor = pow(infinity_lead, k, p)
        else:
            value, derivative = local[u]
            velocity = point[1]*inv_y % p
            factor = value*pow(velocity, k, p)
            if not k:
                log_slope += derivative*velocity*pow(value, -1, p)
        unit = unit*factor % p
    dz = unit*log_slope % p if order == 0 else unit if order == 1 else 0
    return LocalValue(order, unit, dz, first, second)


def from_seeds(p, a, c, modulus, seeds, target):
    if target is None or not (0 <= target[0] < p and 0 < target[1] < p):
        raise ValueError("Expected a canonical finite non-two-torsion target.")
    x, y = target
    if (y*y-x**3-a*x-c) % p:
        raise ValueError("Target is outside the supplied curve.")
    return local_value(p, a, target, lift_windows(p, a, c, modulus, seeds))
