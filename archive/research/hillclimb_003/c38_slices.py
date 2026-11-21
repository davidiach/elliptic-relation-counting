"""Actual C38 first-tensor restrictions, built without a planted relation."""
from itertools import product
from pathlib import Path
import sys

PREVIOUS = Path(__file__).resolve().parent.parent/'hillclimb_002'
sys.path.insert(0, str(PREVIOUS))
from c38_network import Compiled, Pinned, arith, sc, common


def public_instance(b, L, p, m=2):
    assert arith.prime(p) and p.bit_length() == b and 31 % p
    constants = dict(p=p, a=1, Xm=0, Ym=1)
    for j in range(1, m+1):
        constants['alpha'+str(j)] = 1
        constants['beta'+str(j)] = 0
    data = dict(m=m, par=arith.parameters(b, L, (1 << b)-p, 'general'),
                constants=constants, hash_ell=b)
    return Compiled(data, 1)


def source_rows(compiled, overrides=None, row=0):
    """Provide locally valid source values and their actual range traces."""
    data = compiled.data
    p, m = data['par'], data['m']
    raw = {}
    for name in arith.source_names(m, 'general'):
        value = (1 if name == 'x1' else 2) if name.startswith('x') else 0
        raw[name] = arith.limbs(value, p)
    raw.update(overrides or {})
    hs = {name: arith.H(digits, p)[row] for name, digits in raw.items()}
    states = {}
    for name, digits in raw.items():
        bound = arith.limbs((16 if name.startswith('q') else 1)*data['constants']['p'], p)
        trace = [0]
        for digit, threshold in zip(digits[:p['n']], bound[:p['n']]):
            trace.append(int(digit-threshold-trace[-1] < 0))
        states[('range', name)] = trace
    for j in range(2, m+1):
        prev = 'x1' if j == 2 else 'X'+str(j-1)
        trace = [0]
        for a, b in zip(raw['x'+str(j)][:p['n']], raw[prev][:p['n']]):
            trace.append(int(trace[-1] or a != b))
        states[('nonzero', j)] = trace
    tensor = next(t for t in compiled.net.tensors if t.uid == 'F'+str(row))
    symbols = {}
    for binding in tensor.bindings:
        path = binding.path
        if path[0] == 'hs':
            _, name, actual_row, frequency, digit = path
            assert actual_row == row
            symbols[binding.edge] = (hs[name][frequency] // compiled.codec.Q**digit) % compiled.codec.Q
        else:
            field, name, position = path
            symbols[binding.edge] = states[field, name][position]
    assert tensor.entry(tuple(symbols[e] for e in tensor.ports), sc.Ring(1)) == (1,)
    return tensor, symbols, states


def independent_digits(compiled):
    """Find consecutive local raw coefficients with common range-state pins.

    At each free digit use the larger borrow-transition interval; all other
    digits are zero. Stop below the highest nonzero bound digit.
    """
    p = compiled.data['par']
    bound = arith.limbs(16*compiled.data['constants']['p'], p)
    highest = max(j for j, digit in enumerate(bound) if digit)
    positions = [j*p['L2'] for j in range(p['L1']//2) if j*p['L2'] < highest]
    choices = {}
    borrow = 0
    for t in range(p['n']):
        if t in positions:
            split = bound[t]+borrow
            lower = tuple(range(min(split, p['B'])))
            upper = tuple(range(min(split, p['B']), p['B']))
            selected = lower if len(lower) >= len(upper) else upper
            assert len(selected) >= p['B']//2
            assert len(selected) >= 2
            choices[t] = selected
            digit = selected[0]
        else:
            digit = 0
        borrow = int(digit-bound[t]-borrow < 0)
    assert borrow == 1 and choices
    return choices


def fourier_minor(compiled, source_count=1):
    """Return an identity minor of the actual first predicate's Fourier cut."""
    p = compiled.data['par']
    choices = independent_digits(compiled)
    positions = list(choices)
    assert positions == [j*p['L2'] for j in range(len(positions))]
    d = len(positions)
    names = ['q'+str(j) for j in range(source_count)]
    variants = []
    for bits in product((0, 1), repeat=d*source_count):
        overrides = {}
        for j, name in enumerate(names):
            digits = [0]*p['L']
            for k, position in enumerate(positions):
                digits[position] = choices[position][bits[j*d+k]]
            overrides[name] = digits
        variants.append(source_rows(compiled, overrides))
    tensor = variants[0][0]
    selected = {binding.edge:binding.path for binding in tensor.bindings
                if binding.path[0] == 'hs' and binding.path[1] in names}
    left = {edge for edge, path in selected.items() if path[3] < d}
    right = set(selected)-left
    pinned = set(tensor.ports)-set(selected)
    baseline = variants[0][1]
    assert all(all(symbols[e] == baseline[e] for e in pinned) for _, symbols, _ in variants)
    # Verify the two Vandermonde determinants are units modulo the Fourier modulus.
    omega = pow(p['omega'], p['L2'], p['modulus'])
    from math import gcd
    for frequencies in (range(d), range(d, 2*d)):
        points = [pow(omega, i, p['modulus']) for i in frequencies]
        assert all(gcd(points[i]-points[j], p['modulus']) == 1
                   for i in range(d) for j in range(i))
    matrix = []
    for i, (_, a, _) in enumerate(variants):
        row = []
        for j, (_, b, _) in enumerate(variants):
            combined = dict(baseline)
            combined.update({e:a[e] for e in left})
            combined.update({e:b[e] for e in right})
            value = tensor.entry(tuple(combined[e] for e in tensor.ports), sc.Ring(1))[0]
            assert value == int(i == j)
            row.append(value)
        matrix.append(row)
    return dict(b=p['b'], L=p['L'], limb_bits=p['k'], degree=d, sources=source_count,
                minor_dimension=len(matrix), evaluated_entries=len(matrix)**2,
                full_digit_set_sizes=[len(choices[position]) for position in positions],
                exact_larger_minor_dimension=__import__('math').prod(len(v) for v in choices.values())**source_count,
                left_symbol_ports=len(left), right_symbol_ports=len(right), fixed_ports=len(pinned))


class Ones:
    pc = 1
    axes = ()
    def __init__(self, uid, ports, alphabet):
        self.uid, self.ports, self.alphabet = uid, tuple(ports), alphabet
    def entry(self, symbols, ring):
        return ring.one


def copy_restriction(r):
    """r independent Boolean copies inside a genuine finite C38 F0 tensor.

    Uses L=8, four-bit limbs, and unrestricted lower bits of quotient limbs.
    This finite geometry is not the asymptotic C38 packing schedule.
    """
    m = max(2, r//6+2)
    compiled = public_instance(4, 8, 13, m)
    p = compiled.data['par']
    assert p['L1'] == 2 and p['k'] == 4
    count = (r+2)//3
    overrides = {('q'+str(j)): [8]+[0]*(p['L']-1) for j in range(count)}
    tensor, baseline, _ = source_rows(compiled, overrides)
    lookup = {binding.path:binding.edge for binding in tensor.bindings}
    pairs = []
    for j in range(count):
        for digit in range(3):
            if len(pairs) < r:
                pairs.append((lookup['hs', 'q'+str(j), 0, 0, digit],
                              lookup['hs', 'q'+str(j), 0, 1, digit]))
    free = {e for pair in pairs for e in pair}
    fixed = {e:x for e,x in baseline.items() if e not in free}
    restricted = Pinned(tensor, fixed)
    left = Ones('left', [a for a,b in pairs], compiled.codec.alphabet)
    right = Ones('right', [b for a,b in pairs], compiled.codec.alphabet)
    net = sc.Network('actual-C38-copy-restriction', [restricted], [left, right], sorted(free), compiled.codec.alphabet)
    accepted = []
    for bits in product((0, 1), repeat=len(free)):
        trial = dict(zip(net.edges, bits))
        actual = restricted.entry(tuple(trial[e] for e in restricted.ports), sc.Ring(1))[0]
        assert actual == int(all(trial[a] == trial[b] for a,b in pairs))
        if actual:
            accepted.append(trial)
    labels = {restricted.uid:(0,0), 'left':(1,1), 'right':(2,2)}
    return net, labels, pairs, accepted
