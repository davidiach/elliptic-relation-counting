"""Check ordinary quadratic-field incidences against products over Z/p²."""
from collections import Counter
from itertools import product

from check_section_deformation import lift
from fixtures import determinant, evaluation


class QuadraticCurve:
    """Elements x+p*y encode x+y*sqrt(d). The point at infinity is None."""
    def __init__(self, p, a):
        self.p, self.a = p, a
        self.d = next(d for d in range(2, p) if pow(d, (p-1)//2, p) == p-1)

    def add(self, x, y):
        p = self.p
        return (x%p+y%p)%p+p*((x//p+y//p)%p)

    def neg(self, x):
        p = self.p
        return (-x%p)+p*(-(x//p)%p)

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def mul(self, x, y):
        p = self.p
        return ((x%p)*(y%p)+self.d*(x//p)*(y//p))%p+p*(((x%p)*(y//p)+(x//p)*(y%p))%p)

    def inv(self, x):
        p = self.p
        a, b = x%p, x//p
        denominator = pow((a*a-self.d*b*b)%p, -1, p)
        return (a*denominator)%p+p*((-b*denominator)%p)

    def minus(self, point):
        return None if point is None else (point[0], self.neg(point[1]))

    def sum(self, point, other):
        if point is None or other is None:
            return other if point is None else point
        x, y = point
        u, v = other
        if x == u and y == self.neg(v):
            return None
        if x == u:
            slope = self.mul(self.add(self.mul(3, self.mul(x, x)), self.a), self.inv(self.add(y, y)))
        else:
            slope = self.mul(self.sub(v, y), self.inv(self.sub(u, x)))
        new_x = self.sub(self.sub(self.mul(slope, slope), x), u)
        return new_x, self.sub(self.mul(slope, self.sub(x, new_x)), y)

    def row(self, point, size):
        x, y = point
        powers = [1]
        for _ in range(size//2):
            powers.append(self.mul(powers[-1], x))
        return powers+[self.mul(y, v) for v in powers[:(size-1)//2]]

    def singular(self, rows):
        rows = [list(row) for row in rows]
        for j in range(len(rows)):
            pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
            if pivot is None:
                return True
            rows[j], rows[pivot] = rows[pivot], rows[j]
            inverse = self.inv(rows[j][j])
            for i in range(j+1, len(rows)):
                ratio = self.mul(rows[i][j], inverse)
                rows[i] = [self.sub(x, self.mul(ratio, y)) for x, y in zip(rows[i], rows[j])]
        return False


def make_fixture(p, a):
    curve = QuadraticCurve(p, a)
    affine = [(x, y) for x in range(p) for y in range(p) if (y*y-x**3-a*x-1)%p == 0]
    xs = sorted({x for x, y in affine if y})
    slots = [[(point, s) for s, point in enumerate(z for z in affine if z[0] == x)] for x in xs[:3]]
    words = list(product(*slots))
    counts = Counter(curve.sum(curve.sum(w[0][0], w[1][0]), w[2][0]) for w in words)
    targets = {}
    for target in affine:
        if target[0] not in xs[:3] and target[1]:
            targets.setdefault(counts[target], target)
    assert 0 in targets and 1 in targets
    for x in range(p):
        for y in range(1, p):
            if (curve.d*y*y-x**3-a*x-1)%p:
                continue
            marker = (x, p*y)
            multiples = [None, marker]
            for _ in range(7):
                multiples.append(curve.sum(multiples[-1], marker))
            if all(P is not None and (P[0] >= p or P[1] >= p) for P in multiples[1:]):
                return curve, affine, slots, targets, multiples
    raise AssertionError("fixture needs a twist marker with eight nonrational multiples")


def precision_product(p, a, slots, target, masks, endpoint_speed):
    size = len(slots)+1
    rows = [[evaluation(lift(P, (digit & mask).bit_count()%2, a, p, 3), size, p*p)
             for P, digit in slot] for mask, slot in zip(masks, slots)]
    end = evaluation(lift((target[0], -target[1]%p), endpoint_speed, a, p, 3), size, p*p)
    result = 1
    for chosen in product(*rows):
        result = result*determinant([*chosen, end], p*p) % (p*p)
    return result


def check_fixture(p, a):
    curve, affine, slots, targets, multiples = make_fixture(p, a)
    comparisons = determinant_tests = equalities = 0
    for count, target in sorted(targets.items()):
        references = [precision_product(p, a, slots, target, (0, 0, 0), j) for j in range(4)]
        for masks in product(range(2), repeat=3):
            weights = []
            shifted = [[(curve.sum(P, multiples[(mask & digit).bit_count()%2]), digit)
                        for P, digit in slot] for mask, slot in zip(masks, slots)]
            tags = (0, 2)
            enlarged = [curve.sum(P, curve.minus(multiples[j])) for P, _ in shifted[-1] for j in tags]
            assert len(enlarged) == len(set(enlarged)) and None not in enlarged
            sets = [{P[0] for P, _ in slot} for slot in shifted[:-1]]+[{P[0] for P in enlarged}, {target[0]}]
            assert all(not sets[i] & sets[j] for i in range(len(sets)) for j in range(i))
            for w in product(*slots):
                if curve.sum(curve.sum(w[0][0], w[1][0]), w[2][0]) == target:
                    weights.append(sum((digit & mask).bit_count()%2 for (_, digit), mask in zip(w, masks)))
            assert len(weights) == count
            observations = Counter()
            for w in product(*shifted):
                pts = [point for point, _ in w]
                shifted_total = curve.sum(curve.sum(pts[0], pts[1]), pts[2])
                for j in range(4):
                    observations[j] += shifted_total == curve.sum(target, multiples[j])
                for j in tags:
                    last = curve.sum(pts[-1], curve.minus(multiples[j]))
                    expected = curve.sum(shifted_total, curve.minus(multiples[j])) == target
                    determinant_rows = [curve.row(P, 4) for P in [*pts[:-1], last, curve.minus(target)]]
                    assert curve.singular(determinant_rows) == expected
                    determinant_tests += 1
            assert all(observations[j] == weights.count(j) for j in range(4))
            assert bool(sum(observations[j] for j in (0, 2))) == any(k%2 == 0 for k in weights)
            comparisons += 1
            if count == 1:
                moved = precision_product(p, a, slots, target, masks, 0)
                for j in range(4):
                    assert (moved == references[j]) == bool(observations[j])
                    equalities += 1
    torsion_count = sum(y == 0 for _, y in affine)
    print(f"p={p}, A={a}, {torsion_count} rational nonidentity two-torsion points, "
          f"hit counts {sorted(targets)}, {comparisons} complete challenge comparisons.")
    print(f"  Checked {determinant_tests} ordinary determinants with folded tags over F_(p^2), "
          f"{equalities} singleton comparisons against literal products modulo {p*p}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_fixture(41, 8)
    check_fixture(101, 17)
