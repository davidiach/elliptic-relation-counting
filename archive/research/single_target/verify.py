"""Independent literal-product, moment and guarded-character checks on finite curves."""

from collections import Counter
from itertools import product
from random import Random
from target import LocalValue, from_seeds, lift_windows, local_value, point_add, polynomial
from scalar import form, full_jet
from arithmetic_alternating_chord import points


def children(seed):
    alpha, beta, width = seed
    yield seed
    if width > 1:
        yield from children((alpha, beta, width//2))
        yield from children((alpha, beta-width//2, width//2))


def neg(P, p):
    return P[0], -P[1] % p


def literal(slots, T, p, a, c):
    t = T[0]
    h = first = second = dz = 0
    unit = value = 1
    linear = polynomial([-t, 1], p)
    for row in product(*slots):
        xs = [P[0] for P in row]
        f = form(xs, p, a, c)
        v = int(f(t))
        dz, value = (dz*v+value*int(f.derivative()(t))) % p, value*v % p
        while not f(t):
            f = f//linear
            h += 1
        unit = unit*int(f(t)) % p
        code = sum(x*p**i for i, x in enumerate(xs))
        for signed in product(*[(P, neg(P, p)) for P in row]):
            total = None
            for P in signed:
                total = point_add(total, P, p, a)
            if total == T:
                first += code
                second += code*code
    return LocalValue(h, unit, dz, first, second)


def check(slots, T, p, a, c, stats):
    result = local_value(p, a, T, slots)
    assert result == literal(slots, T, p, a, c), (slots, T, result)
    h = result.order
    stats['orders'][h] += 1
    if h:
        xs = [[P[0] for P in slot] for slot in slots]
        grad = result.zero_gradient(p, a, c, len(slots), T[0])
        assert (0, grad) == full_jet(xs, T[0], p, a, c)
        stats['gradients'] += 1
    return result


def small():
    p, a, c = 13, 4, 1
    group = points(p, a, c)
    targets = [P for P in group if 0 < P[1] < p-P[1]]
    stats = {'orders': Counter(), 'gradients': 0}
    parent = ((3, 4, 2), (7, 9, 2), (1, 8, 4))
    for box in product(*(children(seed) for seed in parent)):
        slots = lift_windows(p, a, c, 16, box)
        for slot, (alpha, beta, width) in zip(slots, box):
            assert {P[0] for P in slot} == {x for x, y in group if (alpha*x+beta) % 16 < width}
        for T in targets:
            check(slots, T, p, a, c, stats)
    point = (0, 1)
    for slots in ([[point]], [[point]]*6, [[point, targets[1]]]*4, [[point]]*7):
        for T in targets:
            check(slots, T, p, a, c, stats)
    assert max(stats['orders']) >= p
    print(f"{sum(stats['orders'].values())} literal unit, order, derivative and integer-moment checks, "
          f"{stats['gradients']} gradient checks at zeros, maximum order {max(stats['orders'])}.")


def large():
    p, a, c, T = 65521, 44, 1, (0, 1)
    rng = Random(18112026)
    for _ in range(2):
        seeds = [(rng.randrange(1, 65536, 2), rng.randrange(65536), w) for w in (8, 16, 16, 16)]
    slots = lift_windows(p, a, c, 65536, seeds)
    assert list(map(len, slots)) == [4, 7, 10, 10]
    value = from_seeds(p, a, c, 65536, seeds, T)
    assert value == literal(slots, T, p, a, c)
    assert (value.order, value.unit, value.target_derivative) == (1, 43446, 43446)
    assert value.second_moment == value.first_moment**2
    print(f"2,800-tuple fixture from seeds: {value}.")


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    small()
    large()
