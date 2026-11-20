"""Check signed determinant spectra and power-sum recovery on isolated curve fibers."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from random import Random

from flint import nmod_mat

from check_address_parity import character
from check_strided_isolation import finite_laws
from fixtures import determinant
from probe_signed_moment import Probe


class Spectrum:
    def __init__(self, p, linear, positions):
        self.probe = Probe(p, linear, positions)
        self.seeds = 0
        self.recurrences = 0
        self.signed_support = 0
        self.distinct = 0
        self.max_count = 0
        self.regular_tuples = 0

    @lru_cache(maxsize=32768)
    def value(self, i, j, k):
        rows = [self.probe.vectors[t] for t in (i, j, k, -1)]
        return determinant(rows, self.probe.p)

    def regular(self, slots):
        used = {self.probe.curve[1][0]}
        for slot in slots:
            xs = {self.probe.curve[v][0] for v, _ in slot}
            if xs & used:
                return False
            used.update(xs)
        return True

    def inspect(self, slots, mass, check_recurrences):
        if not self.regular(slots):
            return 0
        rows = []
        block = self.probe.block
        for (i, a), (j, b), (k, c) in product(*slots):
            rows.append(((a*block+b)*block+c, self.value(i, j, k)))
        fibers = Counter(z for _, z in rows)
        self.regular_tuples += mass*len(rows)
        energy = sum(n*n for n in fibers.values())
        address, _ = self.probe.singleton(slots, 1)
        if address is None:
            return energy
        assert fibers[0] == 1 and max(fibers.values()) < self.probe.p
        nonzero = {z: n for z, n in fibers.items() if z}
        supports = 0
        parseval = Counter()
        for r in range(self.probe.size):
            signed = Counter()
            for a, z in rows:
                signed[z] += character(r, a)
            signed = {z: n for z, n in signed.items() if n}
            assert signed[0] == character(r, address)
            for z, n in signed.items():
                parseval[z] += n*n
            actual = sum(z != 0 for z in signed)
            assert actual >= sum(n == 1 for n in nonzero.values())
            supports += actual
            if check_recurrences:
                self.check_recovery(signed, actual, len(nonzero))
                self.recurrences += 1
        assert all(parseval[z] == self.probe.size*n for z, n in fibers.items())
        assert Fraction(supports, self.probe.size) >= sum(Fraction(1, n) for n in nonzero.values())
        self.signed_support += mass*supports
        self.distinct += mass*len(nonzero)
        self.seeds += mass
        self.max_count = max(self.max_count, max(fibers.values()))
        return energy

    def check_recovery(self, signed, expected_order, upper):
        p = self.probe.p
        def moment(j):
            return sum(n*pow(z, j, p) for z, n in signed.items()) % p
        moments = [moment(j) for j in range(2*upper+1)]
        assert (moments[0]-moment(p-1)) % p == signed[0] % p
        def hankel(row_count):
            return nmod_mat([[moments[i+j+1] for j in range(row_count)]
                             for i in range(row_count)], p)
        order = hankel(upper).rank() if upper else 0
        assert order == expected_order
        if order:
            rhs = nmod_mat([[(-moments[i+order+1]) % p] for i in range(order)], p)
            sol = hankel(order).solve(rhs)
            annihilator = [int(sol[i, 0]) for i in range(order)] + [1]
        else:
            annihilator = [1]
        recovered = sum(c*moments[i] for i, c in enumerate(annihilator))
        assert recovered*pow(annihilator[0], -1, p) % p == signed[0] % p

    def report(self, ensemble):
        print(f"p={self.probe.p}, {ensemble}: {self.seeds:,} regular guarded singletons, "
              f"mean nonzero bins {self.distinct/self.seeds:.4f}, mean signed recurrence order "
              f"{self.signed_support/(self.seeds*self.probe.size):.4f}, "
              f"{self.recurrences:,} independent Hankel and Fermat recoveries.")


def exact():
    test = Spectrum(13, 4, ())
    probe = test.probe
    affine, boxes = finite_laws(probe.p, probe.curve, ())
    energy = total = 0
    seen_tails = set()
    for left, middle, right in product(affine, boxes, boxes):
        mass = affine[left]*boxes[middle]*boxes[right]
        energy += mass*test.inspect((left, middle, right), mass, True)
        total += mass
        if test.regular((middle, right)):
            seen_tails.update((j,k) for j, _ in middle for k, _ in right)
    n = probe.n
    rho = Fraction(1, probe.modulus)
    nu = n*rho
    bound = 4*nu**3 + 4*(n*rho**2+2*rho)*2**4
    assert Fraction(energy, total) <= bound and total == 32768
    offdiagonal = Fraction(energy-test.regular_tuples, total)
    refined = 4*(n*rho**2+2*rho)*2**4 + nu*(2*4*rho + n*rho**2)
    assert 0 <= offdiagonal <= refined
    # A rational function difference here has pole order at most four.
    # Comparing all affine field points distinguishes equality from at most four roots.
    section_pairs = 0
    for middle, right in product(boxes, repeat=2):
        if not test.regular((middle, right)):
            continue
        tails = [(j,k) for j, _ in middle for k, _ in right]
        for a in tails:
            for sign in (1, -1):
                alike = 0
                for other in tails:
                    roots = sum(test.value(v, *a) == test.value(sign*v % n, *other)
                                for v in range(1, n))
                    assert roots == n-1 or roots <= 4
                    if roots == n-1 and sign == -1:
                        assert sum(a) % n == 0 and all((x+y) % n == 0 for x,y in zip(a, other))
                    if roots == n-1 and sign == 1 and a != other:
                        extra, = [y for x, y in zip(a, other) if x != y]
                        assert (sum(a)+extra) % n == 1
                    alike += roots == n-1
                    section_pairs += 1
                assert alike <= 2
    assert seen_tails and test.seeds
    test.report("exact width-one mixed law")
    print(f"  {total:,} weighted seeds, regular collision energy {Fraction(energy, total)}, "
          f"bound {bound}; {section_pairs:,} compared pairs of tail sections.")
    print(f"  Off-diagonal collision energy {offdiagonal}, refined bound {refined}.")


def sampled():
    test = Spectrum(101, 17, (0,))
    rng = Random(19007)
    draws = 0
    while test.seeds < 256:
        slots = test.probe.planted(1, rng)
        test.inspect(slots, 1, test.seeds < 8)
        draws += 1
    test.report(f"width-two conditional sample ({draws} planted draws)")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    exact()
    sampled()
