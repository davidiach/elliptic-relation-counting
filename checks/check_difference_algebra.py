"""Top graded membership versus the actual lifted products on regular curve grids."""
from collections import Counter
from itertools import product
from random import Random

from flint import nmod_mat

from check_direct_algebra import Slice, multiply
from check_nonlinear_branches import mixed_fixture
from check_section_deformation import lift
from fixtures import determinant, evaluation
from probe_signed_moment import Probe


def add(f, g, c, p):
    result = dict(f)
    for term, value in g.items():
        result[term] = (result.get(term, 0)+c*value) % p
    return {t: v for t, v in result.items() if v}


def check(p, linear):
    probe = Probe(p, linear, (0,))
    slots = mixed_fixture(probe)
    m, digit_base = 3, probe.block
    unit = {(0,)*m: 1}
    variables = [{tuple(int(i == j) for i in range(m)): 1} for j in range(m)]
    def polynomial(v):
        return {next(iter(var)): c % p for c, var in zip(v, variables) if c % p}
    def mul(fs):
        acc = unit
        for f in fs:
            acc = multiply(acc, f, p)
        return acc

    counts = Counter(sum(i for i, _ in choices) % probe.n for choices in product(*slots))
    used = {probe.curve[i][0] for slot in slots for i, _ in slot}
    targets = {}
    for t in range(1, probe.n):
        if probe.curve[t][0] not in used:
            targets.setdefault(counts[t], t)
    total_probes, comparisons = 0, 0
    rng = Random(p)
    for h, target in sorted(targets.items()):
        endpoint = probe.vectors[-target]
        forms = [[polynomial([a-b for a,b in zip(probe.vectors[i][1:], endpoint[1:])])
                  for i, _ in slot] for slot in slots]
        degrees = [len(f) for f in forms]
        D = sum(k-1 for k in degrees)
        generators = [mul(f) for f in forms]
        top = Slice(m, D, generators, p)
        assert len(top.basis) == max(h, 1)
        if h > 1:
            continue

        prev = Slice(m, D-1, generators, p)
        def coord(f):
            return top.coordinates(f)[0]
        def probe_vector(u):
            u_power = mul([u]*(D-1))
            return [coord(mul((var, u_power))) for var in variables]
        pairing = [[coord(mul((var, {prev.terms[t]: 1}))) for t in prev.basis] for var in variables]
        assert nmod_mat(pairing, p).rank() == (m if h == 0 else 1)

        nonzero_samples = 0
        for _ in range(24):
            u, v, r, s = [[rng.randrange(p) for _ in range(m)] for _ in range(4)]
            col1, col2 = [probe_vector(polynomial(row)) for row in (u, v)]
            minor = determinant([[sum(a*b for a,b in zip(row, col)) for col in (col1, col2)]
                                 for row in (r, s)], p)
            nonzero_samples += int(minor != 0)
            total_probes += 1
        assert bool(nonzero_samples) == (h == 0)
        if h == 0:
            continue

        while True:
            found = probe_vector(polynomial([rng.randrange(p) for _ in range(m)]))
            if any(found):
                break
        digits = []
        for slot, f in zip(slots, forms):
            selected = [digit for (_, digit), form in zip(slot, f)
                        if sum(form.get(next(iter(var)), 0)*v for var, v in zip(variables, found)) % p == 0]
            assert len(selected) == 1
            digits.append(selected[0])
        witnesses = [tuple(d for _, d in choices) for choices in product(*slots)
                     if sum(i for i, _ in choices) % probe.n == target]
        assert [tuple(digits)] == witnesses

        A, H = [], []
        for j, (slot, factors) in enumerate(zip(slots, forms)):
            a, denom = {}, {}
            for i, (_, digit) in enumerate(slot):
                other = mul(factors[:i]+factors[i+1:])
                a = add(a, other, digit*digit_base**j, p)
                denom = add(denom, other, 1, p)
            A.append(a)
            H.append(denom)
        G = {}
        for i, a in enumerate(A):
            G = add(G, mul([a]+H[:i]+H[i+1:]), 1, p)
        denominator = coord(mul(H))
        assert denominator
        decoded = coord(G)*pow(denominator, -1, p) % p
        assert decoded == sum(d*digit_base**j for j,d in enumerate(digits))

        modulus = p*p
        def lifted_product(tagged, speed):
            rows = [[evaluation(lift(probe.curve[i], digit*digit_base**j if tagged else 0,
                                     linear, p, 7), m+1, modulus) for i, digit in slot]
                    for j, slot in enumerate(slots)]
            end = evaluation(lift(probe.curve[-target], speed, linear, p, 7), m+1, modulus)
            result = 1
            for pts in product(*rows):
                result = result*determinant([*pts, end], modulus) % modulus
            return result
        Vc = lifted_product(True, 0)
        vals = set(range(min(p, digit_base**m))) | {p-1}
        for j in vals:
            member = coord(add(G, mul(H), -j, p)) == 0
            assert member == (Vc == lifted_product(False, j)) == (j == decoded)
            comparisons += 1
    print(f'p={p}, slot sizes {[len(s) for s in slots]}, hit counts {sorted(targets)}; '
          f'graded dimensions and empty/singleton pairing ranks verified, '
          f'{total_probes} compact minor tests, {comparisons} true precision-two comparisons.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    check(101, 17)
    check(1031, 9)
