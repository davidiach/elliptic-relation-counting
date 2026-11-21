"""Finite tests of the unary chord potential, its isogeny trace and the marked digit."""
from collections import Counter
from itertools import product
from math import prod

from check_chord_difference import add_with_slope, negative, normalized_difference
from check_nonlinear_branches import mixed_fixture
from probe_signed_moment import Probe


def potential(point, order, p, linear):
    multiple, chord_sum = None, 0
    for bit in bin(order)[2:]:
        multiple, slope = add_with_slope(multiple, multiple, p, linear)
        chord_sum = (2*chord_sum+slope) % p
        if bit == '1':
            multiple, slope = add_with_slope(multiple, point, p, linear)
            chord_sum = (chord_sum+slope) % p
    assert multiple is None
    return chord_sum*pow(order, -1, p) % p


def points(p, linear, constant):
    return [None]+[(x,y) for x,y in product(range(p), repeat=2)
                   if (y*y-x**3-linear*x-constant) % p == 0]


def check_isogenies(p, linear, constant):
    group = points(p, linear, constant)
    order = len(group)
    g = {P: potential(P, order, p, linear) for P in group}
    families = set()
    for P in group:
        assert g[P] == sum(add_with_slope(P,Q,p,linear)[1] for Q in group)*pow(order,-1,p) % p
        kernel, Q = [], None
        while Q not in kernel:
            kernel.append(Q)
            Q = add_with_slope(Q,P,p,linear)[0]
        if 1 < len(kernel) <= 6 and order//len(kernel) >= 3:
            families.add(frozenset(kernel))
    for P,Q in product(group, repeat=2):
        S, slope = add_with_slope(P,Q,p,linear)
        assert (g[P]+g[Q]-g[S]) % p == slope
    two = frozenset(P for P in group if add_with_slope(P,P,p,linear)[0] is None)
    if len(two) == 4:
        for P in group:
            assert sum(g[add_with_slope(P,h,p,linear)[0]] for h in two) % p == 2*g[add_with_slope(P,P,p,linear)[0]] % p
        families.add(two)

    checked = 0
    for H in families:
        tail = [pt for pt in H if pt is not None]
        a = (linear-5*sum(3*x*x+linear for x,y in tail)) % p
        b = (constant-7*sum(5*x*x*x+3*linear*x+2*constant for x,y in tail)) % p
        # Vélu's normalized quotient. Computing the map is part of this finite test.
        image = {}
        for P in group:
            if P in H:
                image[P] = None
                continue
            x,y = P
            for h in tail:
                u,v = add_with_slope(P,h,p,linear)[0]
                x += u-h[0]
                y += v-h[1]
            image[P] = x % p, y % p
            assert (y*y-x*x*x-a*x-b) % p == 0
        image_order = len(set(image.values()))
        assert image_order*len(H) == order
        for P in group:
            traced = sum(g[add_with_slope(P,h,p,linear)[0]] for h in H) % p
            assert traced == potential(image[P], image_order, p, a)
            checked += 1
    print(f'p={p}, order {order}: {order**2} chord identities, {len(families)} normalized '
          f'isogenies of degrees {sorted(set(map(len,families)))}, {checked} coset traces; '
          f'{order if len(two)==4 else 0} full doubling traces.')


def check_mixed():
    p, linear = 101, 17
    probe = Probe(p, linear, (0,))
    masks = mixed_fixture(probe)
    g = {pt: potential(pt, probe.n, p, linear) for pt in probe.curve}
    used = {probe.curve[i][0] for slot in masks for i,d in slot}
    checked = Counter()
    for target in range(1, probe.n):
        if probe.curve[target][0] in used:
            continue
        for mark in (1,3):
            slots = [[(probe.curve[i], (j+1)*((digit & mark).bit_count() % 2))
                      for i,digit in slot] for j,slot in enumerate(masks)]
            slots.append([(probe.curve[-target], mark-1)])
            # Full group convolution in the small fixture. Indices of the enumerated
            # curve are unnecessary here; we use the actual curve points as keys.
            states = {None: (1,0)}
            for slot in slots:
                after = {}
                for P,(n,u) in states.items():
                    for Q,v in slot:
                        S = add_with_slope(P,Q,p,linear)[0]
                        k,w = after.get(S,(0,0))
                        after[S] = k+n,(w+u+n*v) % p
                states = after
            h, speed = states.get(None,(0,0))
            value = 0 if h > 1 else speed
            if h == 0:
                value = -sum(w*g[S] for S,(n,w) in states.items())
                volume = prod(map(len,slots))
                for i,F in enumerate(slots):
                    value += volume//len(F)*sum(u*g[P] for P,u in F)
                    for j,K in enumerate(slots):
                        if i != j:
                            value += volume//(len(F)*len(K))*sum(u for P,u in F)*sum(g[Q] for Q,v in K)
                        if i < j:
                            pair = sum((u-v)*add_with_slope(P,negative(Q,p),p,linear)[1]
                                       for P,u in F for Q,v in K)
                            value += volume//(len(F)*len(K))*pair
            assert normalized_difference(slots,p,linear) == (h, value % p)
            checked[min(h,2)] += 1
    print(f'p={p}, mixed-slot unary formula versus balanced chords: {sum(checked.values())} '
          f'field-digit checks, including {checked[0]} empty, {checked[1]} singleton and '
          f'{checked[2]} multiple-witness cases.')


def check_density():
    p, linear, modulus = 19, -1, 32
    group = points(p, linear, 0)
    add = {(P,Q): add_with_slope(P,Q,p,linear)[0] for P,Q in product(group, repeat=2)}
    cycles = set()
    for P in group:
        entries, Q = set(), None
        while Q not in entries:
            entries.add(Q)
            Q = add[P,Q]
        cycles.add(frozenset(entries))
    kernels = {frozenset(add[P,Q] for P,Q in product(C,D)) for C,D in product(cycles, repeat=2)}
    cosets = {H: {frozenset(add[P,h] for h in H) for P in group} for H in kernels if len(H)>1}
    aff = Counter(tuple(P for P in group[1:] if (a*P[0]+t) % modulus < 2)
                  for a in range(1,modulus,2) for t in range(modulus))
    fixed = Counter(tuple(P for P in group[1:] if P[0]//2 == t) for t in range(modulus//2))
    actual, cover = Counter(), Counter()
    seeds = 0
    for (left, ml), (right, mr) in product(aff.items(),fixed.items()):
        mass = ml*mr
        states = Counter(add[P,Q] for P,Q in product(left,right))
        counts = Counter()
        for P,i in states.items():
            for Q,j in states.items():
                if P != Q:
                    counts[add[negative(P,p),Q]] += i*j
        actual.update({W: mass*k for W,k in counts.items()})
        for H, cells in cosets.items():
            sizes = [len(states.keys() & cell) for cell in cells]
            covered = sum(t for t in sizes if t>=2 and 2*t>=len(H))
            assert covered*len(H) <= 4*sum(counts[W] for W in H if W is not None)
            cover[H] += mass*covered
        seeds += mass
    max_difference = max(actual.values())
    assert all(total*len(H) <= 4*(len(H)-1)*max_difference for H,total in cover.items())
    print(f'p={p}, {seeds} weighted two-slot mixed seeds: collision accounting at density '
          f'1/2 for all {len(cosets)} nontrivial rational subgroups; largest expected number '
          f'of covered distinct half sums for a fixed subgroup {max(cover.values())/seeds:.4f}.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    check_isogenies(17, -1, 0)
    check_isogenies(19, -1, 0)
    check_isogenies(37, 2, 1)
    check_mixed()
    check_density()
