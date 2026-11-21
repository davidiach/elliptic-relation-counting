"""Exhaustive backward-fiber checks and sampled linear batches on two negative fibers."""
from collections import Counter
from itertools import product

import numpy as np


def add(left, right, p, a):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u and (y + v) % p == 0:
        return None
    if x == u:
        slope = (3*x*x + a) * pow(2*y, -1, p) % p
    else:
        slope = (v-y) * pow(u-x, -1, p) % p
    w = (slope*slope - x - u) % p
    return w, (slope*(x-w) - y) % p


def signed_sums(leaves, p, a, b):
    sums = {None}
    for x in leaves:
        lifts = [(x, y) for y in range(p) if (y*y - x**3 - a*x - b) % p == 0]
        assert lifts
        sums = {add(partial, point, p, a) for partial in sums for point in lifts}
    return sums


def residuals(values, leaves, target, p, a):
    m = len(leaves)
    x = {1: leaves[0], m: target[0], **dict(zip(range(2, m), values[:m-2]))}
    slope = dict(zip(range(2, m+1), values[m-2:2*m-3]))
    y = {m: target[1], **dict(zip(range(2, m, 2), values[2*m-3:]))}
    result = []
    for j in range(2, m+1):
        u, v, w, ell = x[j-1], leaves[j-1], x[j], slope[j]
        e = (ell*ell - u - v - w) % p
        leading = w - v - (2*e if j % 2 else 0)
        if u == v or leading % p == 0:
            return None
        if j % 2:
            d = 2*ell*y[j-1] - 2*u*u - u*v - u*w + v*w - a
        else:
            d = 2*ell*y[j] + 2*w*w + u*w + v*w - u*v + a
        result.extend((e, d % p))
    for j in range(3, m+1, 2):
        left_chord = slope[j] * (x[j-1] - x[j])
        if j == m:
            link = left_chord - y[j-1] - y[m]
        else:
            right_chord = slope[j+1] * (x[j] - x[j+1])
            link = left_chord - right_chord - y[j-1] + y[j+1]
        result.append(link % p)
    return tuple(result)


def inverse_count(raw, leaves, target, p, a):
    """Invert each prescribed residual in degree at most two, starting at the target."""
    m = len(leaves)
    links = dict(zip(range(3, m+1, 2), raw[2*(m-1):]))
    current = [target]
    for j in range(m, 1, -1):
        previous = []
        e, d = raw[2*(j-2):2*(j-2)+2]
        v = leaves[j-1]
        for w, y in current:
            leading = (w - v - (2*e if j % 2 else 0)) % p
            if not leading:
                continue
            # R = (w-v)*ell^2 + 2*y*ell + (v-w)*(v+w+e) + 2*w^2+v*w+a.
            constant = (v-w)*(v+w+e) + 2*w*w + v*w + a
            if j % 2:
                linear = 2*(y + links[j])
                constant += d + 2*e*(v+2*w+e)
            else:
                linear = 2*y
                constant -= d
            roots = [ell for ell in range(p)
                     if (leading*ell*ell + linear*ell + constant) % p == 0]
            assert len(roots) <= 2
            for ell in roots:
                u = (ell*ell - v - w - e) % p
                if u != v:
                    prior_y = (ell*(u-w) - y - links.get(j, 0)) % p
                    previous.append((u, prior_y))
        current = previous
    return sum(u == leaves[0] for u, _ in current)


def fixture(p, a, b, leaves, target, expected, expected_hits):
    m = len(leaves)
    variables = 2*m - 3 + (m-1)//2
    equations = variables + 1
    assert (target[1]**2 - target[0]**3 - a*target[0] - b) % p == 0
    assert target not in signed_sums(leaves, p, a, b)
    fibers = Counter()
    for values in product(range(p), repeat=variables):
        raw = residuals(values, leaves, target, p, a)
        if raw is not None:
            fibers[raw] += 1
    zero = (0,)*equations
    assert zero not in fibers and inverse_count(zero, leaves, target, p, a) == 0
    assert max(fibers.values()) <= 2**(m-1)
    for index, (raw, multiplicity) in enumerate(fibers.items()):
        if index % max(1, len(fibers)//64) == 0:
            assert inverse_count(raw, leaves, target, p, a) == multiplicity
    directions = set()
    for raw in fibers:
        unit = pow(next(entry for entry in raw if entry), -1, p)
        directions.add(tuple(unit*entry % p for entry in raw))
    counts = sum(fibers.values()), len(fibers), max(fibers.values()), len(directions)
    assert counts == expected
    print(f"p={p}, m={m}, negative leaves={leaves}, target={target}: "
          f"assignments={counts[0]}, residuals={counts[1]}, max fiber={counts[2]}, "
          f"projective directions={counts[3]}")
    rays = np.array(list(directions), dtype=np.int32)
    rng = np.random.default_rng(23272)
    for rows, expected_yes in zip(range(equations-3, equations), expected_hits):
        yes = 0
        for _ in range(64):
            matrix = rng.integers(0, p, (equations, rows), dtype=np.int32)
            yes += bool(np.any(np.all((rays @ matrix) % p == 0, axis=1)))
        assert yes == expected_yes
        print(f"  {rows} of {equations} random equations: {yes}/64 false acceptances")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    fixture(13, 1, 1, (0, 1, 4), (0, 1), (22308, 20201, 4, 15221), (64, 64, 33))
    fixture(7, 0, 1, (0, 1, 0, 0), (2, 3), (61397, 51006, 6, 43422), (64, 62, 18))
