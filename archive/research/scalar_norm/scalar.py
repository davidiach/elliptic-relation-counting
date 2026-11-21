"""Scalar and first-order finite-field calculations for the research note.

Polynomial arithmetic uses python-flint. Inputs are lists of distinct
abscissas in the paired parts of the original slots.
"""

from itertools import product

from flint import nmod_poly


def polynomial(coeffs, p):
    return nmod_poly([int(v) % p for v in coeffs], p)


class Dual:
    """First-order arithmetic for a polynomial and its directional derivative."""

    def __init__(self, value, tangent=0):
        self.value = value
        self.tangent = tangent

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.value + other.value, self.tangent + other.tangent)
        return Dual(self.value + other, self.tangent)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, -self.tangent)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if isinstance(other, Dual):
            return Dual(self.value*other.value, self.value*other.tangent + self.tangent*other.value)
        return Dual(self.value*other, self.tangent*other)

    __rmul__ = __mul__

    def __getitem__(self, i):
        return Dual(int(self.value[i]), int(self.tangent[i]))


def form(xs, p, a, c, direction=None):
    """The normalized projective summation form in an affine output coordinate."""
    z = polynomial([0, 1], p)
    zero, one = polynomial([], p), polynomial([1], p)
    if direction is not None:
        xs = [Dual(x, int(j == direction)) for j, x in enumerate(xs)]
        z, zero, one = Dual(z, zero), Dual(zero, zero), Dual(one, zero)
    result = z - xs[0]
    degree = 1
    for x in xs[1:]:
        A = (z-x)*(z-x)
        B = -2*((x+z)*(x*z+a)+2*c)
        C = (x*z-a)*(x*z-a)-4*c*(x+z)
        ap, cp, lucas = [one], [one], [2, -B]
        for k in range(1, degree+1):
            ap.append(ap[-1]*A)
            cp.append(cp[-1]*C)
            if k >= 2:
                lucas.append(-B*lucas[-1] - (A*C)*lucas[-2])
        total = zero
        for i in range(degree+1):
            mixed = cp[i]*result[i]
            for j in range(i):
                mixed += cp[j]*result[j]*lucas[i-j]
            total += ap[degree-i]*result[i]*mixed
        result = total
        degree *= 2
    return result if direction is None else result.tangent


def multiply_jets(left, right):
    f, df = left
    g, dg = right
    return f*g, [f*y + g*x for x, y in zip(df, dg)]


def product_tree(items):
    layer = list(items)
    if not layer:
        raise ValueError("The polynomial product tree needs a nonempty list.")
    while len(layer) > 1:
        parents = [multiply_jets(layer[i], layer[i+1]) for i in range(0, len(layer)-1, 2)]
        if len(layer) % 2:
            parents.append(layer[-1])
        layer = parents
    return layer[0]


def half_forms(slots, target, p, a, c, gradients):
    """One entry per half tuple. The backwards half starts with the target."""
    items = []
    for row in product(*slots):
        xs = row if target is None else (target, *reversed(row))
        base = form(xs, p, a, c)
        tangent = []
        if gradients:
            if target is not None:
                tangent.append(form(xs, p, a, c, 0))
            for j, x in enumerate(row):
                index = j if target is None else len(row)-j
                shift = form(xs, p, a, c, index)
                tangent.extend((shift, x*shift))
        items.append((base, tangent))
    degree = len(items) * 2**(len(xs)-1)
    return (*product_tree(items), degree)


def fixed_resultant(A, degree_a, B, degree_b, variations, p):
    """Homogeneous resultant and its derivatives in B with A held fixed."""
    n = A.degree()
    infinity = degree_a - n
    lead = int(A.leading_coefficient())
    monic = A*pow(lead, -1, p)
    finite = int(monic.resultant(B)) if n else 1
    changes = [0] * len(variations)
    if n and variations:
        if finite:
            gcd, inv, _ = B.xgcd(monic)
            assert gcd.is_one()
            trace_kernel = inv*monic.derivative() % monic
            changes = [finite*int((v*trace_kernel % monic)[n-1]) % p for v in variations]
        else:
            common = monic.gcd(B)
            if common.degree() == 1:
                root = -int(common[0]) % p
                rest = monic // common
                # A repeated root of the fixed form contributes zero to this derivative.
                if rest(root):
                    surviving_product = int(rest.resultant(B)) if rest.degree() else 1
                    changes = [surviving_product*int(v(root)) % p for v in variations]
    top = int(B[degree_b])
    base_at_inf = pow(top, infinity, p)
    deriv_at_inf = infinity*pow(top, infinity-1, p) % p if infinity else 0
    factor = (-1 if (infinity*degree_b) % 2 else 1)*pow(lead, degree_b, p)
    changes = [factor*(base_at_inf*d + finite*deriv_at_inf*int(v[degree_b])) % p
               for d, v in zip(changes, variations)]
    return factor*base_at_inf*finite % p, changes


def split_norm(slots, target, p, a, c, gradients=True):
    """Value, target derivative, and (translation, dilation) derivative per slot."""
    if not all(slots):
        return 1, 0, [(0, 0)]*len(slots)
    if len(slots) == 1:
        value, deriv, scale = 1, 0, 0
        for x in slots[0]:
            factor = (target-x) % p
            deriv, scale = (deriv*factor + value) % p, (scale*factor - x*value) % p
            value = value*factor % p
        return value, deriv, [(-deriv % p, scale)]
    cut = len(slots)//2
    A, dA, da = half_forms(slots[:cut], None, p, a, c, gradients)
    B, dB, db = half_forms(slots[cut:], target, p, a, c, gradients)
    value, right = fixed_resultant(A, da, B, db, dB, p)
    if not gradients:
        return value, None, []
    swapped, left = fixed_resultant(B, db, A, da, dA, p)
    sign = -1 if da*db % 2 else 1
    assert value == sign*swapped % p
    shifts = [sign*v % p for v in left] + right[1:]
    return value, right[0], list(zip(shifts[::2], shifts[1::2]))


def two_slot_jet(slots, target, p, a, c):
    """Full additive jet, including at nonzeros. The last two slots are batched."""
    m = len(slots)
    if not all(slots):
        return 1, [0]*(m+2)
    if m == 1:
        value, dz, changes = split_norm(slots, target, p, a, c)
        return value, [dz, dz, -changes[0][1] % p]
    degree = 2**(m-1)
    z = polynomial([0, 1], p)
    one = polynomial([1], p)
    roots = product_tree([(z-x, []) for x in slots[-2]])[0]
    total, gradient = 1, [0]*(m+2)
    for prefix in product(*slots[:-2]):
        leaves = []
        for v in slots[-1]:
            # Symmetry puts the penultimate slot in the output variable.
            xs = (*prefix, v, target)
            leaf = form(xs, p, a, c)
            dz = form(xs, p, a, c, len(xs)-1)
            leaves.append((leaf, [dz, one, *(one*x for x in prefix), z, one*v]))
        joint, variations = product_tree(leaves)
        # The fixed slot form is monic with no infinite roots.
        value, changes = fixed_resultant(roots, roots.degree(), joint,
                                         degree*len(slots[-1]), variations, p)
        gradient = [(old*value + total*change) % p for old, change in zip(gradient, changes)]
        total = total*value % p
    return total, gradient


def full_jet(slots, target, p, a, c, *, evaluate=split_norm):
    """Exact additive derivatives (target, constant marker, coordinate markers)."""
    value, dz, geometric = evaluate(slots, target, p, a, c)
    if not all(slots):
        return value, [0]*(len(slots)+2)
    if value:
        return two_slot_jet(slots, target, p, a, c)
    if dz:
        xs = [dilate*pow(shift, -1, p) % p for shift, dilate in geometric]
        own = int(form(xs, p, a, c).derivative()(target))
        cofactor = dz*pow(own, -1, p) % p
    else:
        # A single abscissa tuple can have several sign completions.
        restricted = [list(u) for u in slots]
        for i in range(len(slots)):
            while len(restricted[i]) > 1:
                mid = len(restricted[i])//2
                low, high = restricted[i][:mid], restricted[i][mid:]
                restricted[i] = low
                if evaluate(restricted, target, p, a, c, False)[0]:
                    restricted[i] = high
        xs = [u[0] for u in restricted]
        cofactor = 1
        for i, x in enumerate(xs):
            box = [[q] for q in xs[:i]] + [[q for q in slots[i] if q != x]] + list(slots[i+1:])
            cofactor = cofactor*evaluate(box, target, p, a, c, False)[0] % p
            if not cofactor:
                break
    return 0, [dz, cofactor, *(x*cofactor % p for x in xs)]


def curve_jet(point_slots, target, p, a, c, add, *, evaluate=split_norm):
    """Include rational two-torsion, using the common local parameter of the project note."""
    if target is None or target[1] == 0:
        raise ValueError("The target must satisfy 2T != O.")
    branches = []
    for slot in point_slots:
        fixed = [(point, None) for point in slot if point[1] == 0]
        paired = sorted({x for x, y in slot if y})
        branches.append(fixed + ([(None, paired)] if paired else []))
    value, gradient = 1, [0]*(len(point_slots)+2)
    for pattern in product(*branches):
        translated = target
        regular = []
        for fixed, paired in pattern:
            if paired is None:
                translated = add(translated, fixed, p, a)
            else:
                regular.append(paired)
        if not regular:
            continue
        v, derivatives = full_jet(regular, translated[0], p, a, c, evaluate=evaluate)
        dz, ds, *regular_marks = derivatives
        marks = iter(regular_marks)
        changes = [dz, ds, *(fixed[0]*ds if paired is None else next(marks)
                             for fixed, paired in pattern)]
        gradient = [(old*v+value*new) % p for old, new in zip(gradient, changes)]
        value = value*v % p
    return value, gradient
