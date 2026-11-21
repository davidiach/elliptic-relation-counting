"""Three-block evaluation of the projective slot norm on small finite fields."""

from collections import Counter
from dataclasses import dataclass
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scalar_norm"))
from scalar import curve_jet, fixed_resultant, form, full_jet, polynomial, product_tree, split_norm


@dataclass
class Block:
    value: object
    degree: int
    roots: Counter
    derivatives: list
    keys: list


def make_block(indices, slot_values, target_index, p, a, c, gradients):
    """Only form factors belonging to one of the three small blocks."""
    leaves = []
    roots = Counter()
    keys = []
    if gradients:
        for i in indices:
            keys.append((i, "shift"))
            if i != target_index:
                keys.append((i, "scale"))
    single_degree = 2**(len(indices)-1)
    for row in product(*(slot_values[i] for i in indices)):
        factor = form(row, p, a, c)
        for point, count in factor.roots():
            roots[int(point)] += int(count)
        if factor.degree() < single_degree:
            roots[None] += single_degree-factor.degree()
        changes = []
        if gradients:
            for j, i in enumerate(indices):
                derivative = form(row, p, a, c, j)
                changes.append(derivative)
                if i != target_index:
                    changes.append(row[j]*derivative)
        leaves.append((factor, changes))
    degree = len(leaves)*single_degree
    if sum(roots.values()) != degree:
        raise ValueError("The paired slots and target must lift to rational curve points.")
    value, derivatives = product_tree(leaves)
    return Block(value, degree, roots, derivatives, keys)


def chord_section(x, y, p, a, c):
    """Homogeneous degree-two section in the final variable, in its affine chart."""
    if x is None and y is None:
        return polynomial([1], p)
    if x is None or y is None:
        finite = x if y is None else y
        return polynomial([-finite, 1], p)**2
    return polynomial([(x*y-a)**2 - 4*c*(x+y),
                       -2*((x+y)*(x*y+a)+2*c), (x-y)**2], p)


def vary_block(outer, middle, moving, p, a, c):
    """Differentiate the moving block; the other two stay at their base values."""
    outer_scale = pow(int(outer.value.leading_coefficient()), 2*middle.degree*moving.degree, p)
    middle_scale = pow(int(middle.value.leading_coefficient()), 2*outer.degree*moving.degree, p)
    value = outer_scale*middle_scale % p
    gradient = [0]*len(moving.derivatives)
    for x, multiplicity in outer.roots.items():
        leaves = [(chord_section(x, y, p, a, c)**count, []) for y, count in middle.roots.items()]
        section, _ = product_tree(leaves)
        v, deriv = fixed_resultant(section, 2*middle.degree, moving.value,
                                   moving.degree, moving.derivatives, p)
        power = pow(v, multiplicity, p)
        slope = multiplicity*pow(v, multiplicity-1, p) % p
        gradient = [(g*power + value*slope*d) % p for g, d in zip(gradient, deriv)]
        value = value*power % p
    return value, gradient


def three_blocks(slots, target, p, a, c, gradients=True):
    """The old left and right products are not inputs to any operation here."""
    m = len(slots)
    if not all(slots):
        return 1, 0, [(0, 0)]*m
    if m == 1:
        return split_norm(slots, target, p, a, c, gradients)
    values = [*slots, [target]]
    blocks = [make_block(list(range(j, m+1, 3)), values, m, p, a, c, gradients)
              for j in range(3)]
    first_value, all_derivatives = None, {}
    for j, moving in enumerate(blocks):
        # The smaller stationary block determines the degree of the temporary section.
        middle, outer = sorted((b for i, b in enumerate(blocks) if i != j), key=lambda b: b.degree)
        value, derivative = vary_block(outer, middle, moving, p, a, c)
        if not gradients:
            return value, None, []
        if first_value is None:
            first_value = value
        else:
            assert value == first_value
        all_derivatives.update(zip(moving.keys, derivative))
    return first_value, all_derivatives[m, "shift"], [
        (all_derivatives[i, "shift"], all_derivatives[i, "scale"]) for i in range(m)]


def implicit_jet(slots, target, p, a, c):
    return full_jet(slots, target, p, a, c, evaluate=three_blocks)


def implicit_curve_jet(point_slots, target, p, a, c, add):
    return curve_jet(point_slots, target, p, a, c, add, evaluate=three_blocks)
