"""Check the fixed-bit transfer and signed group-character Gram identities."""
from cmath import exp, pi
from itertools import product
from math import prod, sqrt

from check_address_parity import character
from check_nonlinear_branches import mixed_fixture
from probe_signed_moment import Probe


def roots_of_unity(n):
    return [exp(2j*pi*x/n) for x in range(n)]


def close(first, second):
    assert abs(first-second) <= 2e-8*(1+max(abs(first), abs(second))), (first, second)


def check_fixed_bits(probe):
    p, n, modulus = probe.p, probe.n, probe.modulus
    roots, field_roots = roots_of_unity(n), roots_of_unity(p)
    free_positions = (0, 2)
    mask = sum(1 << j for j in free_positions)
    offsets = [a for a in range(modulus) if a & ~mask == 0]
    bases = [z for z in range(modulus) if z & mask == 0]
    width = len(offsets)
    assert width*len(bases) == modulus
    maxima = [0.0]*3
    tested_completions = 0
    for mode in range(1, n):
        trace = [sum(roots[(mode*point) % n] for point, _ in probe.by_x[x])
                 for x in range(p)]
        assert all(abs(value.imag) < 1e-10 for value in trace)
        sums = []
        for z in bases:
            retained = [point for point in range(1, n) if probe.curve[point][0] & ~mask == z]
            actual = sum(roots[(mode*point) % n] for point in retained)
            integer_translate = sum(trace[z+a] for a in offsets if z+a < p)
            close(actual, integer_translate)
            sums.append(actual)
        for h in (1, 2, 3):
            moment = sum(abs(value)**(2*h) for value in sums)/len(bases)
            maxima[h-1] = max(maxima[h-1], moment/width**h)
        if mode in (1, 2, n-1):
            # Four distinct offsets, the highest incomplete field block, and every
            # fixed-bit run enter the direct and Fourier sums independently.
            products = [prod(trace[(x+a) % p] for a in offsets) for x in range(p)]
            support = [z for z in range(p-max(offsets)) if z & mask == 0]
            direct = sum(products[z] for z in support)
            completed = 0j
            for k in range(p):
                window = sum(field_roots[(-k*z) % p] for z in support)/p
                hybrid = sum(field_roots[(k*x) % p]*value for x, value in enumerate(products))
                completed += window*hybrid
            close(direct, completed)
            tested_completions += 1
    print(f"p={p}, N={n}, width {width}, two free-bit runs: checked all "
          f"{len(bases)} fixed patterns and {n-1} nontrivial characters.")
    print("  Maximum single-slot normalized even moments for h=1,2,3: "
          + ", ".join(f"{v:.6f}" for v in maxima) + ".")
    print(f"  Checked {tested_completions} incomplete hybrid sums by direct Fourier completion.")


def check_gram(probe):
    n, block, size = probe.n, probe.block, probe.size
    roots = roots_of_unity(n)
    slots = mixed_fixture(probe)
    grid = {}
    for choices in product(*slots):
        (i, a), (j, b), (k, c) = choices
        address = (a*block+b)*block+c
        grid[address] = (i+j+k) % n
    legal = len(grid)
    assert legal == prod(map(len, slots))
    target, witness = next((t, w) for t in range(1, n)
                           if (w := probe.singleton(slots, t)[0]) is not None)
    represented = sorted(set(grid.values()))
    polynomial = [1+0j]
    for other in represented:
        if other == target:
            continue
        denominator = roots[target]-roots[other]
        expanded = [0j]*(len(polynomial)+1)
        for exponent, term in enumerate(polynomial):
            expanded[exponent] -= term*roots[other]/denominator
            expanded[exponent+1] += term/denominator
        polynomial = expanded
    assert len(polynomial) == len(represented)
    for total in represented:
        value = sum(c*roots[(d*total) % n] for d, c in enumerate(polynomial))
        close(value, int(total == target))
    print(f"  Interpolated the singleton indicator on {len(represented)} represented "
          f"group sums with {len(polynomial)} consecutive modes.")
    # A low-coherence selection enables a nonvacuous Gershgorin bound at toy size.
    coefficients = [1+2j, 2-1j, -1+1j, -2+0j]
    ratios = [prod(sum(roots[(d*point) % n] for point, _ in slot) for slot in slots)
              for d in range(n)]
    modes = [0]
    while len(modes) < len(coefficients):
        modes.append(min((d for d in range(n) if d not in modes),
                         key=lambda d: max(abs(ratios[(d-e) % n]) for e in modes)))
    features = []
    for d in modes:
        row = []
        for r in range(size):
            blocks = (r//block**2, (r//block) % block, r % block)
            from_slots = prod(sum(character(a, mask)*roots[(d*point) % n]
                                  for point, a in slot) for slot, mask in zip(slots, blocks))
            from_grid = sum(character(address, r)*roots[(d*total) % n]
                            for address, total in grid.items())
            close(from_slots, from_grid)
            row.append(from_slots)
        features.append(row)
    for d, row in zip(modes, features):
        for e, other in zip(modes, features):
            close(sum(x*y.conjugate() for x, y in zip(row, other))/size, ratios[(d-e) % n])
    scores = [sum(coefficient*row[r] for coefficient, row in zip(coefficients, features))
              for r in range(size)]
    norm = sqrt(sum(abs(value)**2 for value in scores)/size)
    coefficient = sum(x*character(witness, r) for r, x in enumerate(scores))/(norm*size)
    close(coefficient, sum(c*roots[(d*target) % n] for c, d in zip(coefficients, modes))/norm)
    off_diagonal = max(abs(ratios[(d-e) % n]) for d in modes for e in modes if d != e)
    lower_eigenvalue = legal-(len(modes)-1)*off_diagonal
    assert lower_eigenvalue > 0
    bound = sqrt(len(modes)/lower_eigenvalue)
    assert abs(coefficient) <= bound+1e-10 < 1
    print(f"  Checked all {size} challenges for four selected modes on {legal} legal "
          f"addresses and all {len(modes)**2} Gram entries.")
    print(f"  Guarded singleton T=[{target}]G, normalized witness coefficient "
          f"{abs(coefficient):.6f}, Gram upper bound {bound:.6f}.")


def check_composite_projection():
    probe = Probe(101, 6, (0,))
    n, block = probe.n, probe.block
    assert n == 114
    subgroup_size = 6
    quotient_step = n//subgroup_size
    low_modes = range(0, n, quotient_step)
    roots = roots_of_unity(n)
    slots = mixed_fixture(probe)
    target = next(t for t in range(1, n) if probe.singleton(slots, t)[0] is not None)
    grid = {}
    for (i, a), (j, b), (k, c) in product(*slots):
        grid[(a*block+b)*block+c] = (i+j+k) % n
    retained = {v: total for v, total in grid.items() if (total-target) % subgroup_size == 0}
    assert 1 < len(retained) < len(grid)
    products = [prod(sum(roots[(mode*i) % n] for i, _ in slot) for slot in slots)
                for mode in range(n)]
    for d, e in product(range(quotient_step), repeat=2):
        actual = sum(roots[((d-e)*total) % n] for total in retained.values())
        predicted = sum(roots[(-h*target) % n]*products[(d-e+h) % n]
                        for h in low_modes)/subgroup_size
        close(actual, predicted)
        if d == e:
            close(predicted, len(retained))
    weighted = {3: 1j, 3+quotient_step: 2, 9: -1, 9+2*quotient_step: .5j,
                9+3*quotient_step: 1}
    collapsed = {}
    for mode, weight in weighted.items():
        representative = mode % quotient_step
        coefficient = weight*roots[((mode-representative)*target) % n]
        collapsed[representative] = collapsed.get(representative, 0j)+coefficient
    assert len(collapsed) == 2

    def score(r, tuples, weights):
        return sum(character(v, r)*sum(c*roots[(d*t) % n] for d, c in weights.items())
                   for v, t in tuples.items())

    full_norm = restricted_norm = 0
    for r in range(probe.size):
        full_norm += abs(score(r, grid, weighted))**2
        actual = score(r, retained, weighted)
        close(actual, score(r, retained, collapsed))
        restricted_norm += abs(actual)**2
    assert 0 < restricted_norm < full_norm
    print(f"p={probe.p}, composite N={n}: checked {quotient_step**2} restricted Gram "
          f"entries on {len(retained)} of {len(grid)} tuples with six quotient classes.")
    print(f"  Five modes collapsed to two representatives on all {probe.size} challenges; "
          "verified norm contraction.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    fixture = Probe(101, 17, (0,))
    check_fixed_bits(fixture)
    check_gram(fixture)
    check_composite_projection()
