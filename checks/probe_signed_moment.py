"""A joint first-moment arithmetic score tested against same-seed Walsh baselines."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_address_parity import character, transform
from check_strided_isolation import finite_laws
from fixtures import curve_cycle, determinant, evaluation


def exterior(a, b):
    return tuple(a[i]*b[j]-a[j]*b[i] for i, j in ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3)))


class Probe:
    def __init__(self, p, linear, positions):
        self.p = p
        self.curve = curve_cycle(p, linear)
        self.n = len(self.curve)
        self.modulus = 1 << p.bit_length()
        self.positions = positions
        self.width = 1 << len(positions)
        self.block = 2*self.width
        self.size = self.block**3
        self.vectors = [None] + [evaluation(point, 4, p) for point in self.curve[1:]]
        self.by_x = defaultdict(list)
        for i, (x, y) in enumerate(self.curve[1:], 1):
            self.by_x[x].append((i, int(y > (-y) % p)))
        squares = {j*j % p for j in range(1, p)}
        self.legendre = [0 if j == 0 else 1 if j in squares else -1 for j in range(p)]

    @lru_cache(maxsize=1024)
    def signed_rows(self, label_mask):
        rows = []
        for mask in range(self.block):
            row = [0]*4
            for point, address in label_mask:
                sign = character(address, mask)
                for k, value in enumerate(self.vectors[point]):
                    row[k] += sign*value
            rows.append(tuple(x % self.p for x in row))
        return rows

    def spectrum(self, slots, target, compare):
        rows1, rows2, rows3 = map(self.signed_rows, slots)
        endpoint = self.vectors[-target]
        scores = []
        right = []
        for c in rows3:
            x01,x02,x03,x12,x13,x23 = exterior(c, endpoint)
            right.append((x23, -x13, x12, x03, -x02, x01))
        checked = 0
        for i, j in product(range(self.block), repeat=2):
            left = exterior(rows1[i], rows2[j])
            for k, cofactors in enumerate(right):
                z = sum(a*b for a, b in zip(left, cofactors)) % self.p
                if compare and (i, j, k) in ((0, 0, 0), (0, 1, 2), (1, 2, 0), (2, 0, 1)):
                    literal = 0
                    for (a, x), (b, y), (c, w) in product(*slots):
                        sign = character(x, i)*character(y, j)*character(w, k)
                        literal += sign*determinant([self.vectors[v] for v in (a,b,c)]+[endpoint], self.p)
                    assert z == literal % self.p
                    checked += 1
                scores.append(self.legendre[z])
        return transform(scores), checked

    def singleton(self, slots, target):
        valid, hits, marked = [], 0, None
        for (i, a), (j, b), (k, c) in product(*slots):
            address = (a*self.block+b)*self.block+c
            valid.append(address)
            partial = (i+j) % self.n
            if (partial+k) % self.n != target:
                continue
            hits += 1
            if (partial and self.curve[i][0] != self.curve[j][0]
                    and self.curve[partial][0] != self.curve[k][0]):
                marked = address
        assert len(valid) == len(set(valid))
        return (marked if hits == 1 else None), valid

    def planted(self, target, rng):
        while True:
            a, b = rng.randrange(1, self.n), rng.randrange(1, self.n)
            c = (target-a-b) % self.n
            if c:
                break
        x1, x2, x3 = (self.curve[i][0] for i in (a,b,c))
        mul = 1+2*rng.randrange(self.modulus//2)
        offset = (rng.randrange(self.width) - mul*x1) % self.modulus
        inverse = pow(mul, -1, self.modulus)
        affine = []
        for addr in range(self.width):
            x = inverse*(addr-offset) % self.modulus
            affine.extend((i, 2*addr+sign) for i, sign in self.by_x.get(x, ()))
        free_mask = sum(1 << j for j in self.positions)
        fixed = []
        for marked in (x2, x3):
            box = []
            for addr in range(self.width):
                bits = sum(((addr >> j) & 1) << bit for j, bit in enumerate(self.positions))
                x = (marked & ~free_mask) | bits
                box.extend((i, 2*addr+sign) for i, sign in self.by_x.get(x, ()))
            fixed.append(tuple(sorted(box)))
        return tuple(sorted(affine)), *fixed


class Summary:
    def __init__(self, probe, target):
        self.probe, self.target = probe, target
        self.accepted = self.real = self.uniform = self.checks = 0
        self.legal = Fraction()

    def add(self, slots, mass):
        address, valid = self.probe.singleton(slots, self.target)
        if address is None:
            return False
        coeff, checked = self.probe.spectrum(slots, self.target, self.checks < 20)
        self.real += mass*abs(coeff[address])
        self.uniform += mass*sum(map(abs, coeff))
        self.legal += Fraction(mass*sum(abs(coeff[z]) for z in valid), len(valid))
        self.accepted += mass
        self.checks += checked
        return True

    def report(self, description):
        size, m = self.probe.size, self.accepted
        assert m and self.checks >= 20
        values = (Fraction(self.real, size*m), Fraction(self.uniform, size**2*m), self.legal/(size*m))
        print(f"p={self.probe.p}, T=[{self.target}]G, w={self.probe.width}, {description}, "
              f"singletons={m}: witness {float(values[0]):.6f}, uniform {float(values[1]):.6f}, "
              f"legal {float(values[2]):.6f}, {self.checks} explicit determinant comparisons")
        return values


def exact():
    probe = Probe(13, 4, (0,))
    affine, boxes = finite_laws(probe.p, probe.curve, probe.positions)
    totals = [Summary(probe, target) for target in (1, 2)]
    for first, second, third in product(affine, boxes, boxes):
        mass = affine[first]*boxes[second]*boxes[third]
        for summary in totals:
            summary.add((first, second, third), mass)
    for summary in totals:
        summary.report("exact law, 8192 seeds")


def sampled(p, linear, positions, seed, count):
    probe = Probe(p, linear, positions)
    summary = Summary(probe, 1)
    rng = Random(seed)
    draws = 0
    # Planting size-biases a seed by H. Conditional on H=X=1, this exactly gives
    # the intended fresh-seed distribution conditional on the same singleton event.
    while summary.accepted < count:
        summary.add(probe.planted(summary.target, rng), 1)
        draws += 1
    summary.report(f"sampled conditional law, {draws} plant draws")


def elimination(rows, p):
    a = [list(row) for row in rows]
    out = 1
    for j in range(len(a)):
        i = next((i for i in range(j, len(a)) if a[i][j] % p), None)
        if i is None:
            return 0
        if i != j:
            a[i], a[j] = a[j], a[i]
            out = -out
        pivot = a[j][j] % p
        out = out*pivot % p
        for i in range(j+1, len(a)):
            multiple = a[i][j]*pow(pivot, -1, p) % p
            for col in range(j+1, len(a)):
                a[i][col] = (a[i][col]-multiple*a[j][col]) % p
    return out % p


def check_slices(m):
    # Use a field where the quadratic character detects a minus sign.
    probe = Probe(1031, 9, (0,))
    p, dim = probe.p, m+1
    blocks = []
    used = set()
    for pattern in range(0, probe.modulus, 2):
        mask = []
        for k in range(2):
            # First slot has affine hash x-pattern; the others have free least bit.
            mask.extend((v, 2*k+s) for v, s in probe.by_x.get(pattern+k, ()))
        if len(mask) == 4 and not used.intersection(v for v, _ in mask):
            blocks.append(mask)
            used.update(v for v, _ in mask)
        if len(blocks) == m:
            break
    assert len(blocks) == m
    endpoint = evaluation(probe.curve[-1], dim, p)
    even, odd = dim//2+1, (dim-1)//2
    rows = []
    for mask in blocks:
        table = []
        for r in range(probe.block):
            row = [0]*dim
            for v, address in mask:
                sign = character(r, address)
                row = [x+sign*y for x, y in zip(row, evaluation(probe.curve[v], dim, p))]
            row = [x % p for x in row]
            assert not any(row[:even] if r & 1 else row[even:])
            table.append(row)
        rows.append(table)
    allowed = nonzero = 0
    for r in product(range(probe.block), repeat=m):
        inputs = [table[v] for table, v in zip(rows, r)] + [endpoint]
        raw = elimination(inputs, p)
        ej = [j for j in range(m) if not r[j] & 1]
        oj = [j for j in range(m) if r[j] & 1]
        if len(oj) not in (odd-1, odd):
            assert raw == 0
            continue
        allowed += 1
        if len(oj) == odd:
            ej.append(m)
        else:
            oj.append(m)
        permutation = ej+oj
        inversions = sum(i > j for h, i in enumerate(permutation) for j in permutation[h+1:])
        sign = (-1)**inversions
        e_det = elimination([inputs[j][:even] for j in ej], p)
        o_det = elimination([inputs[j][even:] for j in oj], p)
        assert raw == (sign*e_det*o_det) % p
        assert probe.legendre[raw] == probe.legendre[sign % p]*probe.legendre[e_det]*probe.legendre[o_det]
        nonzero += int(raw != 0)
    assert 0 < nonzero <= allowed < probe.block**m
    print(f"p={p}, m={m}: all {probe.block**m} Walsh masks agree with the sign-slice product; "
          f"{allowed} masks on permitted slices, {nonzero} nonzero determinants.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    exact()
    sampled(101, 17, (0,), 43501, 1500)
    sampled(1031, 9, (0,3), 431031, 1000)
    check_slices(4)
    check_slices(5)
