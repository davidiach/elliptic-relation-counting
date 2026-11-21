"""Direct projected-section parity and the exact first graded kernel, on small curves."""
from functools import lru_cache
from itertools import product
from math import prod

from flint import nmod_mat

from check_nonlinear_branches import mixed_fixture
from fixtures import determinant, evaluation
from probe_signed_moment import Probe


@lru_cache(maxsize=None)
def monomials(n, degree):
    if degree < 0:
        return ()
    if n == 1:
        return ((degree,),)
    return tuple((j, *tail) for j in range(degree+1)
                 for tail in monomials(n-1, degree-j))


def multiply(f, g, p):
    result = {}
    for a, x in f.items():
        for b, y in g.items():
            term = tuple(i+j for i, j in zip(a, b))
            result[term] = (result.get(term, 0)+x*y) % p
    return {c: z for c, z in result.items() if z}


class Slice:
    def __init__(self, n, degree, polynomials, p):
        self.p = p
        self.terms = monomials(n, degree)
        relations = []
        for f in polynomials:
            fdegree = sum(next(iter(f)))
            for term in monomials(n, degree-fdegree):
                relations.append(self.vector(multiply({term: 1}, f, p)))
        self.rules = []
        if relations:
            reduced, rank = nmod_mat(relations, p).rref()
            for j in range(rank):
                row = [int(reduced[j, i]) for i in range(len(self.terms))]
                self.rules.append((next(i for i, x in enumerate(row) if x), row))
        pivots = {i for i, _ in self.rules}
        self.basis = [i for i in range(len(self.terms)) if i not in pivots]

    def vector(self, f):
        return [f.get(t, 0) for t in self.terms]

    def coordinates(self, f):
        row = self.vector(f)
        for i, rule in self.rules:
            c = row[i]
            if c:
                row = [(a-c*b) % self.p for a, b in zip(row, rule)]
        return [row[i] for i in self.basis]


def verify(slots, probe, target, tag):
    n, p = len(slots), probe.p
    endpoint = evaluation(probe.curve[-target], n+1, p)
    forms = []
    polynomials = []
    for slot in slots:
        linear = {}
        for point, addr in slot:
            vec = evaluation(probe.curve[point], n+1, p)
            f = {tuple(int(i == j) for i in range(n)): (vec[j+1]-endpoint[j+1]) % p
                 for j in range(n)}
            linear[addr] = {a: c for a, c in f.items() if c}
        polynomial = {(0,)*n: 1}
        for f in linear.values():
            polynomial = multiply(polynomial, f, p)
        polynomials.append(polynomial)
        forms.append(linear)

    hits, grid = [], []
    for choices in product(*slots):
        points, addresses = zip(*choices)
        det = determinant([evaluation(probe.curve[v], n+1, p) for v in points]+[endpoint], p)
        hit = sum(points) % probe.n == target
        assert hit == (det == 0)
        address = 0
        for a in addresses:
            address = address*probe.block+a
        grid.append((address, det))
        if hit:
            hits.append(addresses)
    assert len(hits) <= 1

    @lru_cache(maxsize=None)
    def quotient(degree):
        return Slice(n, degree, polynomials[:-1], p)

    degrees = [len(slot) for slot in slots]
    critical = sum(a-1 for a in degrees[:-1])
    length = prod(degrees[:-1])
    found = None
    for k in range(critical+2):
        source, dest = quotient(k), quotient(k+degrees[-1])
        rows = [dest.coordinates(multiply({source.terms[i]: 1}, polynomials[-1], p))
                for i in source.basis]
        kernel, dimension = nmod_mat(rows, p).transpose().nullspace()
        assert dimension == (len(hits) if k >= critical else 0)
        if k >= critical:
            assert len(source.basis) == length
        if k == critical and dimension:
            separator = {source.terms[j]: int(kernel[i, 0]) for i, j in enumerate(source.basis)
                         if kernel[i, 0]}
            found = []
            after = quotient(k+1)
            for slot in forms:
                zeros = [a for a, f in slot.items()
                         if not any(after.coordinates(multiply(f, separator, p)))]
                assert len(zeros) == 1
                found.append(zeros[0])
            assert tuple(found) == hits[0]
    marked = 0
    if found:
        for a in found:
            marked = marked*probe.block+a
    challenges = probe.block**n
    for r in range(challenges):
        fermat = sum((-1)**((a & r).bit_count())*(1-pow(z, p-1, p)) for a,z in grid) % p
        direct = 0 if found is None else (-1)**((marked & r).bit_count())
        assert fermat == direct % p
    print(f"p={p}, m={n}, {tag}, T=[{target}]G: hit count {len(hits)}; "
          f"critical source degree {critical}, quotient dimension {length}, "
          f"{challenges} exact Fermat-parity comparisons.")


def projective_count():
    probe = Probe(13, 4, (0,))
    p = probe.p
    # All three masks contain the target's abscissa. The affine anchor and
    # the other two boxes use free positions 0, 2 and 3 respectively.
    slots = [[(point, 2*a+sign) for a, x in enumerate(xs)
              for point, sign in probe.by_x.get(x, ())] for xs in ((0, 1), (0, 4), (0, 8))]
    endpoint = probe.vectors[-1]
    sections = [(1, a, b) for a,b in product(range(p), repeat=2)]
    sections += [(0, 1, a) for a in range(p)]+[(0, 0, 1)]
    point_sections = {}
    for slot in slots:
        for point, _ in slot:
            row = [x-y for x, y in zip(probe.vectors[point][1:], endpoint[1:])]
            point_sections[point] = {s for s in sections if sum(x*y for x,y in zip(s, row)) % p == 0}
    sizes, terms = set(), []
    for (i,a), (j,b), (k,c) in product(*slots):
        count = len(point_sections[i] & point_sections[j] & point_sections[k])
        det = determinant([probe.vectors[v] for v in (i,j,k,-1)], p)
        assert count % p == (1-pow(det, p-1, p)) % p
        sizes.add(count)
        terms.append(((a*probe.block+b)*probe.block+c, det))
    assert sizes == {0, 1, p+1, p*p+p+1}
    for r in range(probe.size):
        weights = []
        for i, slot in enumerate(slots):
            digit = r//probe.block**(2-i) % probe.block
            weights.append([sum((-1)**((digit & a).bit_count()) for v,a in slot
                                if s in point_sections[v]) for s in sections])
        actual = sum(prod(col) for col in zip(*weights))
        fermat = sum((-1)**((a & r).bit_count())*(1-pow(d, p-1, p)) for a,d in terms)
        assert (actual-fermat) % p == 0
    print(f"p={p}, overlapping slots: all four nullities and {probe.size} signed "
          f"projective-incidence comparisons over {len(sections)} sections.")


def main():
    projective_count()
    probe = Probe(101, 17, (0,))
    slots = mixed_fixture(probe)
    for target, expected in ((12, 15), (13, 16)):
        # For two selected slots the residual third zero is T-i-j. Two
        # equal unordered triples give the same projective section.
        distinct_sections = {tuple(sorted((i, j, (target-i-j) % probe.n)))
                             for i, _ in slots[0] for j, _ in slots[1]}
        assert len(distinct_sections) == expected
        verify(slots, probe, target, f"{expected} distinct partial sections")
    verify(slots, probe, 2, "empty control")

    # Use disjoint width-two masks with exactly one liftable abscissa each.
    slots = []
    for offset in range(0, probe.p-1, 2):
        slot = [(i, 2*j+sign) for j in (0, 1) for i, sign in probe.by_x.get(offset+j, ())]
        if len(slot) == 2:
            slots.append(slot)
        if len(slots) == 4:
            break
    used = {probe.curve[i][0] for slot in slots for i, _ in slot}
    sums = [sum(v[0] for v in choices) % probe.n for choices in product(*slots)]
    target = next(t for t in range(1, probe.n) if sums.count(t) == 1 and probe.curve[t][0] not in used)
    verify(slots, probe, target, "four-slot singleton")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
