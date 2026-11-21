"""Finite checks of the run formula and the bounds for actual dyadic hashes."""

from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import comb
from pathlib import Path
from random import Random
import sys

from digit import target_digit

research = Path(__file__).resolve().parents[1]
# Both directories contain a verify.py, so load the prior fixture by path.
spec = spec_from_file_location("slot_arity_checks", research / "slot_arity" / "verify.py")
slots = module_from_spec(spec)
spec.loader.exec_module(slots)
from arithmetic_alternating_chord import point_add


def reference(h, prime):
    power = slots.ordinary_product([h]*prime, len(h))
    pushed = [0]*len(h)
    for j, value in enumerate(h):
        pushed[(prime*j) % len(h)] += value
    assert all((a-b) % prime == 0 for a, b in zip(power, pushed))
    return [(a-b)//prime % prime for a, b in zip(power, pushed)]


def run_checks():
    rng = Random(254919)
    checked = 0
    for order, step in ((1, 0), (8, 2), (10, 1), (13, 4)):
        def add(a, b):
            return (a+b) % order
        for prime in (2, 3, 5, 7):
            fixtures = [[], [(0, 1, prime)], [(0, 2, -1), (0, 2, 1)],
                        [(0, order+2, -prime), (0, order+1, prime+1)]]
            fixtures.extend([[(rng.randrange(order), rng.randrange(order+4),
                               rng.randrange(-prime-1, prime+2)) for _ in range(2)]
                             for _ in range(3)])
            for runs in fixtures:
                direct = [0]*order
                for start, length, weight in runs:
                    for j in range(length):
                        direct[(start+j*step) % order] += weight
                expected = reference(direct, prime)
                for target in range(order):
                    assert target_digit(runs, step, target, prime, add, 0) == expected[target]
                    checked += 1
    # Exercise the implementation with curve points and the group law as its
    # only access to points. Cyclic indices are confined to the independent reference.
    p, a = 13, 4
    cycle = slots.curve_cycle(p, a)
    runs = [(cycle[2], 5, 2), (cycle[9], 4, -1), (None, 1, 3)]
    h = [0]*len(cycle)
    for start, length, weight in runs:
        index = cycle.index(start)
        for j in range(length):
            h[(index+3*j) % len(cycle)] += weight
    for prime in (2, 3, 5):
        expected = reference(h, prime)
        for target, answer in zip(cycle, expected):
            assert target_digit(runs, cycle[3], target, prime, lambda x,y: point_add(x,y,p,a), None) == answer
            checked += 1
    print(f"Direct digit vs ordinary powers: {checked} target checks.")


def long_run():
    order, length, prime, target = 1019, 64, 11, 346
    # Maximum sum is 693, so an integer coefficient is an independent reference.
    ways = sum((-1)**j*comb(prime, j)*comb(target-j*length+prime-1, prime-1)
               for j in range(target//length+1))
    assert target % prime and ways % prime == 0
    orbits = ways//prime
    value = target_digit([(0, length, 1)], 1, target, prime, lambda x,y: (x+y) % order, 0)
    assert value == orbits % prime
    print(f"Single run of {length} points, prime {prime}: {orbits} nonconstant orbits, digit {value}.")


def difference_histogram(left, right, order):
    return Counter((v-u) % order for u in left for v in right)


def hash_moments():
    cycle = slots.curve_cycle(13, 4)
    order, modulus, width = len(cycle), 16, 2
    rho = Fraction(width, modulus)
    all_masks = Counter(tuple(i for i, bit in enumerate(slots.window(cycle, modulus, width, a, b)) if bit)
                        for a in range(1, modulus, 2) for b in range(modulus))
    seeds = sum(all_masks.values())
    selfs = {f: difference_histogram(f, f, order) for f in all_masks}
    for q in range(1, order):
        expected = Fraction(sum(n*selfs[f][q] for f, n in all_masks.items()), seeds)
        assert expected <= order*rho**2+rho
    moment = good = 0
    for (f, nf), (g, ng) in product(all_masks.items(), repeat=2):
        cross = difference_histogram(f, g, order)
        collisions = sum(c*(c-1) for c in cross.values())
        assert collisions == sum(selfs[f][q]*selfs[g][q] for q in range(1, order))
        moment += nf*ng*collisions
        if collisions == 0:
            good += nf*ng
        union = set(f) | set(g)
        near = difference_histogram(union, union, order)
        for q in range(1, order):
            lower = min(len(f), len(g))-cross[0]-min(selfs[f][q], selfs[g][q])-cross[q]-cross[-q % order]
            assert len(union)-near[q] >= lower
            if collisions == 0:
                assert len(union)-near[q] >= min(len(f), len(g))-3
    expected = Fraction(moment, seeds**2)
    bound = order**3*rho**4+2*order**2*rho**3+order*rho**2
    assert expected <= bound
    assert seeds**2-good <= seeds**2*(order-1)*(order*rho**2+rho)**2
    print(f"All {seeds**2} seed pairs: cross factorial moment {expected}, upper bound {bound}.")


def padding_checks():
    order = 9
    masks = [[int(i in pair) for i in range(order)] for pair in ((1, 8), (2, 7), (1, 8))]
    checked = 0
    for m, prime in ((2, 3), (2, 7), (3, 5), (3, 7)):
        d = prime-m
        projected = [0]*order
        for subset in range(1 << m):
            h = [sum(masks[j][i] for j in range(m) if subset & (1 << j)) for i in range(order)]
            powers = [[1]+[0]*(order-1)]
            for _ in range(m):
                powers.append(slots.convolution(powers[-1], h))
            stencil = [0]*order
            low_degree = [0]*order
            for t in range(d+1):
                weight = (-1)**(d-t)*comb(d,t)
                padded = h.copy()
                padded[0] += t
                exact = reference(padded, prime)
                stencil = [a+weight*b for a,b in zip(stencil, exact)]
                low_degree[0] += weight*((t**prime-t)//prime)
                for j in range(1, m+1):
                    factor = weight*(comb(prime,j)//prime)*t**(prime-j)
                    low_degree = [a+factor*b for a,b in zip(low_degree, powers[j])]
            assert all((x-y) % prime == 0 for x,y in zip(stencil,low_degree))
            sign = (-1)**(m-subset.bit_count())
            projected = [a+sign*b for a,b in zip(projected,low_degree)]
            checked += order
        literal = slots.ordinary_product(masks[:m],order)
        assert all((x+y) % prime == 0 for x,y in zip(projected,literal))  # (prime-1)! = -1.
    print(f"Padding cancellation and mixed coefficient: {checked} target checks.")


def digit_rank_check():
    cycle = slots.curve_cycle(13, 4)
    order, modulus, width, prime, m, target = len(cycle), 16, 2, 5, 3, 1
    windows = [[slots.window(cycle, modulus, width, a, b) for b in range(modulus)] for a in (1, 3, 5)]
    extra = windows[2]
    b = next(i for i, w in enumerate(extra) if w != extra[i-1])
    difference = [x-y for x,y in zip(extra[b],extra[b-1])]
    centered = []
    for left in windows[0]:
        f = [x-y for x,y in zip(left,windows[0][0])]
        k = slots.convolution(f, difference)
        row = []
        for right in windows[1]:
            g = [x-y for x,y in zip(right,windows[1][0])]
            row.append(-slots.one_target(k,g,target) % prime)
        centered.append(row)
    combined = [[0]*modulus for _ in range(modulus)]
    ranks = []
    d = prime-m
    for t in range(d+1):
        stacked = []
        for background in (extra[b-1], extra[b]):
            for left in windows[0]:
                row = []
                for right in windows[1]:
                    mix = [x+y+z for x,y,z in zip(left,right,background)]
                    mix[0] += t
                    row.append(reference(mix,prime)[target])
                stacked.append(row)
        ranks.append(len(slots.basis(stacked,prime)))
        diff = [[v-u for u,v in zip(stacked[i],stacked[modulus+i])] for i in range(modulus)]
        for i,j in product(range(modulus),repeat=2):
            combined[i][j] += (-1)**(d-t)*comb(d,t)*(diff[i][j]-diff[0][j]-diff[i][0]+diff[0][0])
    assert [[v % prime for v in row] for row in combined] == centered
    rank = len(slots.basis(centered,prime))
    assert sum(ranks) >= rank >= (order-1)//2-4
    print(f"Single-target digit mod {prime}: ranks {ranks} on two-background slices, mixed rank {rank}.")


def hashed_digit_reconstruction():
    p, a, prime, m = 13, 4, 5, 3
    cycle = slots.curve_cycle(p,a)
    masks = []
    for alpha in (1, 3, 5):
        candidates = [slots.window(cycle,16,2,alpha,b) for b in range(16)]
        masks.append(next(f for f in candidates if sum(f) == 2 and f not in masks))
    literal = slots.ordinary_product(masks,len(cycle))
    queries = 0
    for target in (0, literal.index(max(literal))):
        total = 0
        for subset in range(1 << m):
            runs = [(point,1,1) for j,f in enumerate(masks) if subset & (1 << j)
                    for point,selected in zip(cycle,f) if selected]
            for t in range(prime-m+1):
                weight = (-1)**(prime-subset.bit_count()-t)*comb(prime-m,t)
                value = target_digit(runs+[(None,1,t)],cycle[1],cycle[target],prime,
                                     lambda x,y: point_add(x,y,p,a),None)
                total += weight*value
                queries += 1
        assert -total % prime == literal[target] % prime
    print(f"Actual three-slot masks: {queries} direct digit queries reconstruct the two target residues.")


def actual_threshold():
    p, a, order, arity, modulus = 65521, 44, 65419, 4, 65536
    squares = {y*y % p for y in range(1, p)}
    rhs = [(x**3+a*x+1) % p for x in range(p)]
    assert 1+sum(2 if v in squares else int(v == 0) for v in rhs) == order
    assert slots.is_prime(p) and slots.is_prime(order) and (4*a**3+27) % p
    width = next(1 << e for e in range(17)
                 if order*(1 << e) >= 20*arity*modulus
                 and (order*(1 << e))**arity >= order*(2*modulus)**arity)
    base, point = (0, 1), None
    rng = Random(7341)
    seeds = [(rng.randrange(1, modulus, 2), rng.randrange(modulus)) for _ in range(arity)]
    masks = [set() for _ in seeds]
    for j in range(1, order):
        point = point_add(point, base, p, a)
        for f, (alpha, beta) in zip(masks, seeds):
            if (alpha*point[0]+beta) % modulus < width:
                f.add(j)
    assert point_add(point, base, p, a) is None
    for parts in (masks[:2], masks):
        union = set.union(*parts)
        hist = difference_histogram(union, union, order)
        best_runs = len(union)-max(hist[q] for q in range(1, order))
        assert best_runs > len(union)//2
        print(f"Order {order}, width {width}, {len(parts)} masks, {len(union)} points: "
              f"minimum {best_runs} support runs over all {order-1} steps.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run_checks()
    long_run()
    hash_moments()
    padding_checks()
    digit_rank_check()
    hashed_digit_reconstruction()
    actual_threshold()
