"""Exact finite-law checks of the mixed-slot factorial moment and signed-address advice."""
from collections import Counter
from fractions import Fraction
from itertools import product


def curve(p, a, c):
    points = [None] + [(x, y) for x, y in product(range(p), repeat=2)
                       if (y*y-x**3-a*x-c) % p == 0]
    indices = {point: i for i, point in enumerate(points)}

    def add(left, right):
        if left is None or right is None:
            return right if left is None else left
        x, y = left
        u, v = right
        if x == u and (y+v) % p == 0:
            return None
        numerator, denominator = (3*x*x+a, 2*y) if x == u else (v-y, u-x)
        slope = numerator*pow(denominator, -1, p) % p
        new_x = (slope*slope-x-u) % p
        return new_x, (slope*(x-new_x)-y) % p

    sums = [[indices[add(left, right)] for right in points] for left in points]
    return points, sums


def finite_laws(p, points, positions):
    modulus = 1 << p.bit_length()
    width = 1 << len(positions)
    free = sum(1 << j for j in positions)
    affine = Counter()
    boxes = Counter()
    for alpha in range(1, modulus, 2):
        for beta in range(modulus):
            mask = []
            for i, (x, y) in enumerate(points[1:], 1):
                address = (alpha*x+beta) % modulus
                if address < width:
                    mask.append((i, 2*address + int(y > (-y) % p)))
            affine[tuple(mask)] += 1
    for pattern in range(modulus):
        if pattern & free:
            continue
        mask = []
        for i, (x, y) in enumerate(points[1:], 1):
            if x & ~free == pattern:
                addr = sum(((x >> bit) & 1) << j for j, bit in enumerate(positions))
                mask.append((i, 2*addr + int(y > (-y) % p)))
        boxes[tuple(mask)] += 1
    return affine, boxes


def numbers(p, n, m, positions):
    modulus = 1 << p.bit_length()
    width = 1 << len(positions)
    rho = Fraction(width, modulus)
    lam = n**(m-1) * rho**m
    runs = sum(bit-1 not in positions for bit in positions)
    factorial_bound = rho**2 * (2*width)**(2*m-2) * (
        1 + Fraction(2*(m-1)*2**runs, width))
    return lam, factorial_bound


def verify(data, total, p, points, m, positions):
    lam, bound = numbers(p, len(points), m, positions)
    first, second, guarded, unique, absent, multiple = [
        Fraction(count, total) for count in data]
    assert second <= bound
    assert guarded <= first <= lam
    assert guarded >= (1-Fraction(3*m, len(points)))*lam
    assert unique >= guarded-second
    assert multiple <= second/2 and absent+unique+multiple <= 1
    return second, unique, multiple


def full_law():
    p, positions, m = 13, (0,), 3
    points, sums = curve(p, 4, 1)
    n = len(points)
    affine, boxes = finite_laws(p, points, positions)
    counts = [[0]*6 for _ in points]
    total = 0
    address_bits = 1+len(positions)
    for first, second, third in product(affine, boxes, boxes):
        mass = affine[first]*boxes[second]*boxes[third]
        total += mass
        raw = Counter()
        guarded_addresses = [set() for _ in points]
        for (i, ai), (j, aj), (k, ak) in product(first, second, third):
            subtotal = sums[i][j]
            target = sums[subtotal][k]
            raw[target] += 1
            if subtotal and points[i][0] != points[j][0] and points[subtotal][0] != points[k][0]:
                key = (ai << (2*address_bits)) + (aj << address_bits) + ak
                assert key not in guarded_addresses[target]
                guarded_addresses[target].add(key)
        for target in range(1, n):
            h = raw[target]
            x = len(guarded_addresses[target])
            singleton = h == x == 1
            values = (h, h*(h-1), x, int(singleton), int(h == 0), int(h >= 2))
            for w, v in enumerate(values):
                counts[target][w] += mass*v
            if singleton:
                key, = guarded_addresses[target]
                for bit in range(m*address_bits):
                    for chosen in (0, 1):
                        actual = sum((value >> bit) & 1 == chosen for value in guarded_addresses[target])
                        assert actual == int((key >> bit) & 1 == chosen)

    selected = [verify(counts[t], total, p, points, m, positions)
                for t in range(1, n) if sums[t][t]]
    assert len(selected) == n-1 and all(entry[1] and entry[2] for entry in selected)
    print(f"Field {p}: {total:,} mixed seeds, {len(selected)} targets, "
          f"per-target guarded-singleton counts {min(v[3] for v in counts[1:])}.."
          f"{max(v[3] for v in counts[1:])}.")


def width_one():
    p, m = 127, 3
    points, sums = curve(p, 1, 1)
    n = len(points)
    neg = {i: j for i, row in enumerate(sums) for j, value in enumerate(row) if value == 0}
    target = next(i for i in range(1, n) if sums[i][i])
    raw, guarded = Counter(), Counter()
    # At width one each slot chooses one uniform abscissa in the full binary universe.
    # Count all marked ordered pairs and force the third point. Empty roots are the complement.
    for i, j in product(range(1, n), repeat=2):
        subtotal = sums[i][j]
        k = sums[target][neg[subtotal]]
        if not k:
            continue
        code = tuple(points[v][0] for v in (i, j, k))
        raw[code] += 1
        if subtotal and points[i][0] != points[j][0] and points[subtotal][0] != points[k][0]:
            guarded[code] += 1
    total = (1 << p.bit_length())**m
    unique = sum(h == guarded[code] == 1 for code, h in raw.items())
    many = sum(h > 1 for h in raw.values())
    data = (sum(raw.values()), sum(h*(h-1) for h in raw.values()), sum(guarded.values()),
            unique, total-len(raw), many)
    factorial, _, _ = verify(data, total, p, points, m, ())
    lam, bound = numbers(p, n, m, ())
    assert (1-Fraction(3*m, n))*lam-bound > 0 and unique and many and factorial > 0
    print(f"Field {p}: {total:,} width-one abscissa patterns, {unique:,} guarded singletons, "
          f"{many:,} ambiguous, {total-len(raw):,} absent. The proven isolation lower bound is positive.")
    single = next(code for code, h in raw.items() if h == guarded[code] == 1)
    check_batch(p, points, sums, target, single)


def check_batch(p, points, sums, target, single):
    """Exact probability over all public two-row challenges for one free field scalar."""
    def projective(vector):
        first = next(x for x in vector if x % p)
        return tuple(x*pow(first, -1, p) % p for x in vector)

    slots = [[i for i, point in enumerate(points) if i and point[0] == x] for x in single]
    tx, ty = points[target]
    normals = ([(1, y, z) for y, z in product(range(p), repeat=2)]
               + [(0, 1, z) for z in range(p)] + [(0, 0, 1)])
    bad_rays = set()
    solutions = 0
    for i, j, k in product(*slots):
        subtotal = sums[i][j]
        if not subtotal or points[i][0] == points[j][0] or points[subtotal][0] == points[k][0]:
            continue
        result = points[sums[subtotal][k]]
        r1, r2 = (1, 0) if result is None else ((result[0]-tx) % p, (result[1]-ty) % p)
        checksum = (points[i][1]+points[j][1]+points[k][1]) % p
        for z in range(p):
            residual = (r1, r2, (z-checksum) % p)
            if any(residual):
                bad_rays.add(projective(residual))
            else:
                solutions += 1
    assert solutions == 1
    bad_planes = sum(any((nx*x+ny*y+nz*z) % p == 0 for x, y, z in bad_rays)
                     for nx, ny, nz in normals)
    # Rank-2 matrices have a uniformly random kernel ray, rank-1 matrices a kernel plane.
    # The all-zero matrix accepts every assignment.
    rank2 = (p**3-1)*(p**3-p)
    rank1 = len(normals)*(p**2-1)
    assert rank1+rank2+1 == p**6
    bad_matrices = 1 + rank1*bad_planes//len(normals) + rank2*len(bad_rays)//len(normals)
    failure = Fraction(bad_matrices, p**6)
    assert failure <= Fraction(8, p)
    print(f"Field {p}: exact two-check failure probability {float(failure):.5f} "
          f"over {p**6:,} challenge matrices on one guarded singleton (bound 8/{p}).")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    full_law()
    width_one()
