"""Finite checks for the slot-arity note. Run with python -B and assertions on."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, factorial, isqrt, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hillclimb_001" / "vendor"))
from arithmetic_alternating_chord import point_add, points


def convolution(f, g):
    order = len(f)
    out = [0] * order
    for i, x in enumerate(f):
        if x:
            for j, y in enumerate(g):
                if y:
                    out[(i+j) % order] += x*y
    return out


def ordinary_product(factors, order):
    out = [1] + [0] * (order-1)
    for f in factors:
        out = convolution(out, f)
    return out


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def curve_cycle(p, a):
    affine = points(p, a, 1)
    r = len(affine)+1
    assert (4*a**3+27) % p and is_prime(r)
    cycle = [None]
    for _ in affine:
        cycle.append(point_add(cycle[-1], affine[0], p, a))
    assert set(cycle[1:]) == set(affine) and point_add(cycle[-1], affine[0], p, a) is None
    return cycle


def window(cycle, modulus, tau, alpha, beta):
    return [int(pt is not None and (alpha*pt[0]+beta) % modulus < tau) for pt in cycle]


def complete_classes(xs, modulus, tau):
    sizes = Counter(x % tau for x in xs)
    return sum(size == modulus//tau for size in sizes.values())


def basis(rows, prime):
    """Return independent original rows, using Q when prime is zero."""
    pivots = {}
    originals = []
    for original in rows:
        v = [x % prime for x in original] if prime else list(map(Fraction, original))
        for j in sorted(pivots):
            if v[j]:
                v = [x-v[j]*y for x, y in zip(v, pivots[j])]
                if prime:
                    v = [x % prime for x in v]
        pivot = next((i for i, x in enumerate(v) if x), None)
        if pivot is None:
            continue
        if prime:
            inv = pow(v[pivot], -1, prime)
            pivots[pivot] = [x*inv % prime for x in v]
        else:
            pivots[pivot] = [x/v[pivot] for x in v]
        originals.append(original)
    return originals


def one_target(f, g, t):
    return sum(x*g[(t-i) % len(g)] for i, x in enumerate(f))


def check_ranks(p, a):
    cycle = curve_cycle(p, a)
    r, modulus = len(cycle), 1 << p.bit_length()
    xs = {pt[0] for pt in cycle[1:]}
    n = len(xs)
    assert 2*n+1 == r
    tests = 0
    for tau in (1 << e for e in range(p.bit_length())):
        c = complete_classes(xs, modulus, tau)
        slot_rank = n-max(0, c-1)
        for alpha in (1, 3, modulus-1):
            windows = [window(cycle, modulus, tau, alpha, t) for t in range(modulus)]
            for prime in (0, 2):
                assert len(basis(windows, prime)) == slot_rank
                tests += 1
            fixed = next(f for f in windows if any(f))
            outputs = [convolution(f, fixed) for f in windows]
            for target in (0, 1, n):
                matrix = ([one_target(f, g, target) for g in windows] for f in outputs)
                assert len(basis(matrix, 0)) >= max(0, 2*slot_rank-(n+1))
                tests += 1
            # Contract a third shift axis by an actual adjacent-window difference.
            beta = next(i for i, h in enumerate(windows) if h != windows[i-1])
            difference = [x-y for x, y in zip(windows[beta], windows[beta-1])]
            shifted = [convolution(f, difference) for f in windows]
            for target in (1, n):
                lower = 2*slot_rank-n
                matrix = ([one_target(f, g, target) for g in windows] for f in shifted)
                assert len(basis(matrix, 2)) >= lower
                ranks = []
                for background in (windows[beta], windows[beta-1]):
                    rows = [convolution(f, background) for f in windows]
                    matrix = ([one_target(f, g, target) for g in windows] for f in rows)
                    ranks.append(len(basis(matrix, 2)))
                assert 2*max(ranks) >= lower
                tests += 2
    # A 2+2 slot split. Convolve bases of actual shifted windows, then
    # evaluate the fixed-target pairing on independent rows of each side.
    tau = 2 if p == 13 else 4
    assert complete_classes(xs, modulus, tau) == 0
    for prime, expected in ((1009, n+1), (2, n)):
        windows = [window(cycle, modulus, tau, 3, t) for t in range(modulus)]
        slots = basis(windows, prime)
        halves = basis((convolution(f, g) for f, g in product(slots, repeat=2)), prime)
        assert len(halves) == expected
        for target in range(r):
            actual = len(basis(([one_target(f, g, target) for g in halves] for f in halves), prime))
            assert actual == (0 if prime == 2 and target == 0 else expected)
            tests += 1
    print(f"Order {r}: {tests} rank checks, balanced ranks {n+1} over Q and {n} mod 2 at nonidentity targets.")
    return cycle


def large_window_span():
    p, a, modulus, m = 1031, 9, 2048, 8
    squares = {y*y % p for y in range(1, p)}
    values = [(x**3+a*x+1) % p for x in range(p)]
    assert (4*a**3+27) % p and 0 not in values
    xs = [x for x, v in enumerate(values) if v in squares]
    r = 2*len(xs)+1
    assert r == 1019 and is_prime(r) and r >= 30*m
    tau = next(1 << e for e in range(p.bit_length()+1)
               if r*(1 << e) >= 20*m*modulus and (r*(1 << e))**m >= r*(2*modulus)**m)
    assert tau == 512 and complete_classes(xs, modulus, tau) == 0
    # Bit-packed elimination on all 2,048 shifts, independent of the class formula.
    pivots = {}
    for beta in range(modulus):
        row = sum(1 << i for i, x in enumerate(xs) if (7*x+beta) % modulus < tau)
        while row:
            j = row.bit_length()-1
            if j not in pivots:
                pivots[j] = row
                break
            row ^= pivots[j]
    assert len(pivots) == len(xs) == 509
    print(f"Coverage-threshold fixture: order {r}, arity {m}, width {tau}, shift span {len(pivots)}.")


def repeated_hashes(cycle):
    r = len(cycle)
    geometry = (16, 2) if r == 19 else (32, 4)
    modulus, tau = geometry
    windows = [window(cycle, modulus, tau, 1, beta) for beta in range(modulus)]
    visible_shifts = [0] * r
    for h in windows:
        power = ordinary_product([h]*8, r)
        for t, coefficient in enumerate(power):
            bit = coefficient % 2
            assert bit == h[(t*pow(8, -1, r)) % r]
            visible_shifts[t] += bit
    assert visible_shifts == [0]+[tau]*(r-1)
    f = max(windows, key=sum)
    candidates = [window(cycle, modulus, tau, 3, beta) for beta in range(modulus)]
    g = max((h for h in candidates if h != f), key=sum)
    assert any(f) and any(g)
    configs = (([f], [8]), ([f, g], [3, 4]), ([f, g], [5, 2]))
    for masks, multiplicities in configs:
        ordinary = ordinary_product([h for h, m in zip(masks, multiplicities) for _ in range(m)], r)
        parity_factors = []
        labelled_tokens = []
        for j, (h, m) in enumerate(zip(masks, multiplicities)):
            for e in range(m.bit_length()):
                weighted_points = [(i << e) % r for i in range(r) if h[i]]
                labelled_tokens.extend((j, 1 << e, u) for u in weighted_points)
                if (m >> e) & 1:
                    parity_factors.append([int(i in weighted_points) for i in range(r)])
        parity = ordinary_product(parity_factors, r)
        assert all((a-b) % 2 == 0 for a, b in zip(ordinary, parity))
        # Enumerate sums of binary chunks, enforcing the multiplicity of each
        # mask. Check residues of literal counts against the allowed support.
        start = (tuple(0 for _ in masks), 0)
        states = {start}
        cumulative = set()
        support_at_precision = []
        for d in range(1, sum(multiplicities)+1):
            next_states = set()
            for used, target in states:
                for j, size, u in labelled_tokens:
                    if used[j]+size <= multiplicities[j]:
                        nxt = list(used)
                        nxt[j] += size
                        next_states.add((tuple(nxt), (target+u) % r))
            states = next_states
            cumulative.update(t for used, t in states if list(used) == multiplicities)
            q = d-len(parity_factors)+1
            if q > 0:
                visible = {i for i, x in enumerate(ordinary) if x % (1 << q)}
                assert visible <= cumulative
                support_at_precision.append(len(visible))
        if multiplicities == [8]:
            print(f"Order {r}, eight equal slots: slot size {sum(f)}, ordinary support "
                  f"{sum(x > 0 for x in ordinary)}, support at precisions 1..8 {support_at_precision}.")


def polarization(cycle):
    r = len(cycle)
    modulus, tau = (16 if r == 19 else 32), 1
    masks = []
    for j in range(8):
        candidates = [window(cycle, modulus, tau, 2*j+1, beta) for beta in range(modulus)]
        masks.append(next(f for f in candidates if sum(f) == 2 and f not in masks))
    for m in (3, 4, 8):
        probed = [0] * r
        for subset in range(1, 1 << m):
            selected = [masks[j] for j in range(m) if (subset >> j) & 1]
            mix = list(map(sum, zip(*selected)))
            power = ordinary_product([mix]*m, r)
            sign = (-1)**(m-len(selected))
            probed = [x+sign*y for x, y in zip(probed, power)]
        literal = ordinary_product(masks[:m], r)
        assert any(x % 2 for x in literal)
        assert probed == [factorial(m)*x for x in literal]
        v = m-m.bit_count()
        assert all((z >> v) % 2 == n % 2 for z, n in zip(probed, literal))
        assert all(z % (1 << v) == 0 for z in probed)
        residue_modulus = 1 << max(literal).bit_length()
        inverse_odd = pow(factorial(m) >> v, -1, residue_modulus)
        recovered = [((z % (residue_modulus << v)) >> v)*inverse_odd % residue_modulus
                     for z in probed]
        assert recovered == literal
    print(f"Order {r}: polarization and the first informative bit checked at arities 3, 4 and 8.")


def prime_square_coefficient(h, prime, target):
    # Use literal powers as the reference implementation on tiny groups.
    return ordinary_product([h]*prime, len(h))[target] % (prime*prime)


def prime_square_reconstruction(cycle):
    order, m = len(cycle), 4
    modulus = 16 if order == 19 else 32
    tau = modulus//4
    masks = []
    for alpha in (1, 3, 5, 7):
        candidates = [window(cycle, modulus, tau, alpha, beta) for beta in range(modulus)]
        masks.append(max((f for f in candidates if f not in masks), key=sum))
    literal = ordinary_product(masks, order)
    upper = min(prod(sum(f) for j, f in enumerate(masks) if i != j) for i in range(m))
    primes = []
    prime_product = 1
    for p in range(m, upper+2):
        if is_prime(p):
            primes.append(p)
            prime_product *= p
            if prime_product > upper:
                break
    assert len(primes) > 1 and prime_product > upper >= max(literal) > primes[0]
    # Include the identity and a target whose count genuinely needs several residues.
    targets = sorted({0, 1, literal.index(max(literal))})
    queries = 0
    for target in targets:
        answer, covered_modulus = 0, 1
        for p in primes:
            digit_sum = 0
            raw_sum = 0
            for bits in range(1 << m):
                mix = [0]*order
                for j, f in enumerate(masks):
                    if bits & (1 << j):
                        mix = [x+y for x, y in zip(mix, f)]
                for padding in range(p-m+1):
                    mix[0] = padding
                    weight = (-1)**(p-bits.bit_count()-padding)*comb(p-m, padding)
                    raw = prime_square_coefficient(mix, p, target)
                    raw_sum += weight*raw
                    pushed = sum(x for i, x in enumerate(mix) if (p*i) % order == target)
                    excess = (raw-pushed) % (p*p)
                    assert excess % p == 0
                    digit_sum += weight*(excess//p)
                    queries += 1
            canonical = raw_sum % (p*p)
            assert canonical % p == 0 and canonical//p == digit_sum % p
            residue = digit_sum*pow(factorial(p-1), -1, p) % p
            assert residue == literal[target] % p
            step = ((residue-answer)*pow(covered_modulus, -1, p)) % p
            answer += covered_modulus*step
            covered_modulus *= p
        assert answer == literal[target]
    print(f"Order {order}: prime-square reconstruction from {primes}, "
          f"{queries} one-coefficient queries, largest recovered count {max(literal)}.")


def shift_expectation(cycle):
    order, modulus, tau = len(cycle), 16, 2
    masks = [[window(cycle, modulus, tau, a, b) for b in range(modulus)] for a in (3, 5, 7)]
    total = [0] * order
    for left, middle in product(*masks[:2]):
        pair = convolution(left, middle)
        for right in masks[2]:
            for t, count in enumerate(convolution(pair, right)):
                total[t] += count
    all_affine = [0]+[1]*(order-1)
    unconstrained = ordinary_product([all_affine]*3, order)
    assert total == [tau**3*count for count in unconstrained]
    for arity in (2, 4, 8, 32, 128):
        prime_product = prod(p for p in range(arity+1, 16*arity+1) if is_prime(p))
        assert prime_product > 10*4**arity
    print("Three-slot first moment checked for every target across 4096 independent shifts.")


def noninvertible_pushforward_check():
    # The reconstruction uses raw powers correctly even for a prime
    # dividing the group order, when its point pushforward has no inverse.
    order, prime = 10, 5
    factors = [[int(i in pair) for i in range(order)] for pair in ((1, 9), (2, 8), (3, 7))]
    literal = ordinary_product(factors, order)
    total = [0]*order
    for subset in range(1 << len(factors)):
        mix = [sum(f[i] for j, f in enumerate(factors) if subset & (1 << j)) for i in range(order)]
        for t in range(3):
            mix[0] = t
            weight = (-1)**(prime-subset.bit_count()-t)*comb(2, t)
            total = [a+weight*b for a, b in zip(total, ordinary_product([mix]*prime, order))]
    for raw, count in zip(total, literal):
        assert raw % (prime**2) == factorial(prime)*count % (prime**2)
    print("Raw prime-square cancellation checked at all targets on the order-10 group, prime 5.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    for p, a in ((13, 4), (23, 5)):
        toy = check_ranks(p, a)
        repeated_hashes(toy)
        polarization(toy)
        prime_square_reconstruction(toy)
        if p == 13:
            shift_expectation(toy)
    noninvertible_pushforward_check()
    large_window_span()
