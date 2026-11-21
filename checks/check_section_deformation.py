"""Verify the tangent determinant product and its quadratic-jet character."""
from collections import Counter
from itertools import product

from fixtures import determinant, evaluation
from probe_signed_moment import Probe


def tangent_row(point, dimension, linear, p):
    x, y = point
    dx, dy = 2*y, 3*x*x+linear
    def power_derivative(j):
        return 0 if j == 0 else j*pow(x, j-1, p)*dx
    return ([power_derivative(j) % p for j in range(dimension//2+1)]
            + [(dy*pow(x, j, p)+y*power_derivative(j)) % p
               for j in range((dimension-1)//2)])


def quadratic_product(x, y, d, p):
    return ((x[0]*y[0]+d*x[1]*y[1]) % p,
            (x[0]*y[1]+x[1]*y[0]) % p)


def lift(point, speed, linear, p, offset):
    x, y = point
    error = (y*y-x**3-linear*x-1)//p
    free = (offset+x+3*y) % p
    if y % p:
        dx = free
        dy = ((3*x*x+linear)*dx-error)*pow(2*y, -1, p) % p
    else:
        dy = free
        dx = (error+2*y*dy)*pow(3*x*x+linear, -1, p) % p
    a = (x+p*(dx+speed*2*y)) % (p*p)
    b = (y+p*(dy+speed*(3*x*x+linear))) % (p*p)
    assert (b*b-a**3-linear*a-1) % (p*p) == 0
    return a, b


def precision_two(probe, linear, slots, target, hits, offset):
    p, m = probe.p, len(slots)
    modulus = p*p
    def full_product(move_slots, move_endpoint):
        rows = [[evaluation(lift(probe.curve[i], (digit*4**j if move_slots else 0),
                                 linear, p, offset), m+1, modulus) for i, digit in slot]
                for j, slot in enumerate(slots)]
        end = evaluation(lift(probe.curve[-target], int(move_endpoint), linear, p, offset), m+1, modulus)
        value = 1
        for chosen in product(*rows):
            value = value*determinant([*chosen, end], modulus) % modulus
        return value
    base = full_product(False, False)
    reference = full_product(False, True)
    moved = full_product(True, False)
    observed = 0 if base % p else 1 if base != reference else 2
    assert observed == min(2, len(hits))
    if observed == 1:
        k = ((moved-base)//p)*pow((reference-base)//p, -1, p) % p
        assert k == sum(digit*4**j for j, digit in enumerate(hits[0]))


def check(p, linear, m):
    probe = Probe(p, linear, (0,))
    slots = []
    for offset in range(0, p-1, 2):
        slot = [(point, 2*j+sign) for j in range(2)
                for point, sign in probe.by_x.get(offset+j, ())]
        if len(slot) == 4:
            slots.append(slot)
        if len(slots) == m:
            break
    assert len(slots) == m and 4**m <= p
    grid = list(product(*slots))
    used = {probe.curve[i][0] for slot in slots for i, _ in slot}
    counts = Counter(sum(v for v, _ in choices) % probe.n for choices in grid)
    targets = {}
    for t in range(1, probe.n):
        if probe.curve[t][0] not in used:
            targets.setdefault(counts[t], t)
    assert {0, 1, 2}.issubset(targets)

    d = next(i for i in range(2, p) if probe.legendre[i] == -1)
    shifts = [a for a in range(p) if all(probe.legendre[((a-j)**2-d) % p] == (-1)**j
                                          for j in range(m+1))]
    assert shifts
    a = shifts[0]
    total = 0
    for count, t in sorted(targets.items()):
        nonzero = 1
        derivatives, hits = [], []
        for choices in grid:
            pts = [probe.curve[i] for i, _ in choices]+[probe.curve[-t]]
            rows = [evaluation(pt, m+1, p) for pt in pts]
            value = determinant(rows, p)
            if value:
                nonzero = nonzero*value % p
            else:
                partials = [determinant(rows[:j]+[tangent_row(pt, m+1, linear, p)]+rows[j+1:], p)
                            for j, pt in enumerate(pts)]
                assert partials[0] and len(set(partials)) == 1
                derivatives.append((tuple(digit for _, digit in choices), partials))
            if sum(i for i, _ in choices) % probe.n == t:
                hits.append(tuple(digit for _, digit in choices))
        assert len(derivatives) == len(hits) == count
        for offset in (0, 7):
            precision_two(probe, linear, slots, t, hits, offset)

        def leading(speeds):
            value = (nonzero, 0)
            for digits, partials in derivatives:
                base = sum(partials[j]*x for j, x in enumerate(speeds(digits)))
                coeff = (base-a*partials[-1], -partials[-1])
                value = quadratic_product(value, coeff, d, p)
            return value

        for masks in product(range(4), repeat=m):
            def speeds(digits):
                return [((digit & mask).bit_count() % 2) for digit, mask in zip(digits, masks)]
            x, y = leading(speeds)
            expected = (-1)**sum(sum(speeds(digits)) for digits in hits)
            assert probe.legendre[(x*x-d*y*y) % p] == expected
            total += 1
        marker_x, marker_y = leading(lambda _: [1]+[0]*(m-1))
        assert probe.legendre[(marker_x**2-d*marker_y**2) % p] == (-1)**count
        if count == 1:
            x, y = leading(lambda digits: [v*4**j for j, v in enumerate(digits)])
            recovered = (a-x*pow(y, -1, p)) % p
            assert recovered == sum(v*4**j for j, v in enumerate(hits[0]))
    print(f'p={p}, m={m}, {len(grid)} signed tuples, {len(shifts)} parity-encoding shifts '
          f'for d={d}; hit counts {sorted(targets)}, {total} jet-character checks, '
          'all count markers and packed singleton addresses recovered.')
    print(f'  Three ordinary products modulo {p*p}, two choices of lifts: '
          'empty/single/multiple classification and singleton labels agree with enumeration.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    check(101, 17, 3)
    check(1031, 9, 4)
