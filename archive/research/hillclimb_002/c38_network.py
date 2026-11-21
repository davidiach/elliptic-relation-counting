"""Executable C38 graph descriptors; compilation uses public data only.

General-prime branch only. Every state relay has two distinct edge IDs.
Original tensors receive only their incident symbols. No planted witness
or full assignment is captured by a tensor entry program.
"""
from collections import Counter, defaultdict
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
import sys

PREVIOUS = Path(__file__).resolve().parent.parent / 'hillclimb_001'
sys.path.insert(0, str(PREVIOUS))
import common
import c38_simplex as bridge
import simplex_contraction as sc

arith = bridge.arithmetic


class NonincidentRead(RuntimeError):
    pass


class LocalMap:
    """An assignment interface that makes out-of-incidence reads fatal."""
    def __init__(self, entries):
        self.entries = entries

    def __getitem__(self, key):
        if key not in self.entries:
            raise NonincidentRead(repr(key))
        return self.entries[key]


class LocalSequence:
    def __init__(self, length, entries):
        self.length, self.entries = length, entries

    def __len__(self):
        return self.length

    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self[i] for i in range(*key.indices(self.length))]
        if key < 0:
            key += self.length
        if not 0 <= key < self.length:
            raise IndexError(key)
        if key not in self.entries:
            raise NonincidentRead(str(key))
        return self.entries[key]

    def __iter__(self):
        return (self[i] for i in range(self.length))


@dataclass(frozen=True)
class Binding:
    edge: int
    path: tuple


class Predicate:
    pc = 1
    axes = ()

    def __init__(self, uid, data, codec, kind, index, bindings):
        self.uid, self.data, self.codec = uid, data, codec
        self.alphabet, self.kind, self.index = codec.alphabet, kind, index
        self.bindings = tuple(bindings)
        self.ports = tuple(binding.edge for binding in bindings)

    def local_assignment(self, symbols):
        d, c = self.data, self.codec
        p, m = d['par'], d['m']
        leaves = {b.path: symbol for b, symbol in zip(self.bindings, symbols)}
        def read(path):
            if path not in leaves:
                raise NonincidentRead(repr(path))
            return leaves[path]
        def bit(path):
            value = read(path)
            if value not in (0, 1):
                raise ValueError('nonbinary auxiliary symbol')
            return value
        def word(path):
            value = sum(read(path + (k,))*c.Q**k for k in range(c.width))
            if value >= c.modulus:
                raise ValueError('noncanonical Fourier word')
            return arith.word(value, p)
        ass = {}
        if self.kind in ('first', 'middle'):
            hs = {}
            for name in arith.source_names(m, 'general'):
                js = (self.index,) if self.kind == 'first' else range(p['L2'])
                ts = range(p['L1']) if self.kind == 'first' else (self.index,)
                hs[name] = LocalSequence(p['L2'], {
                    j: LocalSequence(p['L1'], {t: word(('hs', name, j, t)) for t in ts})
                    for j in js})
            ass['hs'] = LocalMap(hs)
        if self.kind in ('middle', 'closure'):
            js = range(p['L2']) if self.kind == 'middle' else range(self.index*p['R'], (self.index+1)*p['R'])
            ts = (self.index,) if self.kind == 'middle' else range(p['L1'])
            ass['outs'] = [LocalSequence(p['L2'], {
                j: LocalSequence(p['L1'], {t: word(('outs', eq, j, t)) for t in ts})
                for j in js}) for eq in range(arith.arithmetic_count(m)+m)]
        if self.kind == 'first':
            for field, keys in (('range', arith.source_names(m, 'general')),
                                ('nonzero', range(2, m+1))):
                states = {}
                for key in keys:
                    needed = set()
                    for j1 in range(p['L1']):
                        t = j1*p['L2']+self.index
                        if t < p['n']:
                            needed.update((t, t+1))
                    values = {t: (0 if t == 0 else 1 if t == p['n'] else bit((field, key, t)))
                              for t in needed}
                    states[key] = LocalSequence(p['n']+1, values)
                ass[field] = LocalMap(states)
        if self.kind == 'closure':
            carries = []
            for eq in range(arith.arithmetic_count(m)+m):
                stop = p['L'] if eq < arith.arithmetic_count(m) else (p['b']+p['k']-1)//p['k']
                needed = set()
                for start in range(self.index*p['R'], stop, p['L2']):
                    end = min(start+p['R'], stop)
                    if start:
                        needed.add(start)
                    if end < stop:
                        needed.add(end)
                carries.append(LocalMap({t: [bit(('carry_ports', eq, t, k))
                                             for k in range(p['carry_bits'])] for t in needed}))
            ass['carry_ports'] = carries
        return LocalMap(ass)

    def entry(self, symbols, ring):
        if len(symbols) != len(self.ports):
            raise ValueError('incorrect arity')
        if any(type(x) is not int or not 0 <= x < self.alphabet.size for x in symbols):
            return ring.zero
        try:
            assignment = self.local_assignment(symbols)
            getattr(arith, self.kind+'_predicate')(self.data, assignment, self.index)
            return ring.one
        except (AssertionError, ValueError):
            return ring.zero
        # NonincidentRead is deliberately not caught.


class Relay:
    pc = 1
    axes = ()

    def __init__(self, uid, ports, alphabet):
        self.uid, self.ports, self.alphabet = uid, tuple(ports), alphabet

    def entry(self, symbols, ring):
        return ring.val(int(len(symbols) == 2 and symbols[0] in (0, 1) and symbols[0] == symbols[1]))


class Pinned:
    def __init__(self, base, pins):
        self.base, self.uid = base, base.uid
        self.alphabet, self.pc, self.axes = base.alphabet, base.pc, base.axes
        self.fixed = {e: pins[e] for e in base.ports if e in pins}
        self.ports = tuple(e for e in base.ports if e not in pins)

    def entry(self, symbols, ring):
        values = dict(self.fixed)
        values.update(zip(self.ports, symbols))
        return self.base.entry(tuple(values[e] for e in self.base.ports), ring)


class Compiled:
    def __init__(self, data, t=2):
        assert data['par']['mode'] == 'general'
        assert set(data) == {'m', 'par', 'constants', 'hash_ell'}
        self.data = deepcopy(data)
        self.codec = bridge.Codec(t, data['par']['modulus'])
        p, m = self.data['par'], self.data['m']
        self.bindings, self.endpoints, self.paths = defaultdict(list), {}, {}
        self.categories = Counter()
        relays = []
        def edge(a, b, path, category):
            e = len(self.endpoints)
            self.endpoints[e] = (a, b)
            self.paths[e] = path
            self.bindings[a].append(Binding(e, path))
            self.bindings[b].append(Binding(e, path))
            self.categories[category] += 1
            return e
        def relay(a, b, path, category):
            name = 'relay:'+str(len(relays))
            ids = (edge(a, name, path, category), edge(b, name, path, category))
            relays.append(Relay(name, ids, self.codec.alphabet))
        for name in arith.source_names(m, 'general'):
            for j in range(p['L2']):
                for k in range(p['L1']):
                    for digit in range(self.codec.width):
                        edge('F'+str(j), 'M'+str(k), ('hs', name, j, k, digit), 'simplex_fourier')
        for eq in range(arith.arithmetic_count(m)+m):
            for j in range(p['L2']):
                for k in range(p['L1']):
                    for digit in range(self.codec.width):
                        edge('E'+str(j//p['R']), 'M'+str(k), ('outs', eq, j, k, digit), 'simplex_fourier')
            stop = p['L'] if eq < arith.arithmetic_count(m) else (p['b']+p['k']-1)//p['k']
            for boundary in range(p['R'], stop, p['R']):
                a, b = 'E'+str(((boundary-1) % p['L2'])//p['R']), 'E'+str((boundary % p['L2'])//p['R'])
                for digit in range(p['carry_bits']):
                    relay(a, b, ('carry_ports', eq, boundary, digit), 'binary_carry_routes')
        for field, keys in (('range', arith.source_names(m, 'general')), ('nonzero', range(2, m+1))):
            for key in keys:
                for boundary in range(1, p['n']):
                    relay('F'+str((boundary-1) % p['L2']), 'F'+str(boundary % p['L2']),
                          (field, key, boundary), field)
        predicates = []
        for prefix, kind, size in (('F', 'first', p['L2']), ('M', 'middle', p['L1']),
                                   ('E', 'closure', p['L2']//p['R'])):
            for j in range(size):
                uid = prefix+str(j)
                predicates.append(Predicate(uid, self.data, self.codec, kind, j, self.bindings[uid]))
        rows = [x for x in predicates if x.kind != 'middle']
        cols = [x for x in predicates if x.kind == 'middle'] + relays
        self.net = sc.Network('C38-general', rows, cols, tuple(self.endpoints), self.codec.alphabet)
        expected = bridge.ledger(data, self.codec)
        assert len(self.net.edges) == expected['total_symbol_edges']
        original = arith.graph_incidence(data)
        assert len(self.net.tensors) == original['vertices']
        for key in ('binary_carry_routes', 'range', 'nonzero'):
            assert self.categories[key] == original['categories'][key]

    def encode_witness(self, original):
        """A test/serialization helper; never used by compilation or entries."""
        encoded = bridge.encode_assignment(self.data, original, self.codec)
        values = {}
        for edge, path in self.paths.items():
            value = encoded
            for key in path:
                value = value[key]
            values[edge] = self.codec.rank(value)
        return values

    def decode_witness(self, assignment):
        """Recover the old assignment from an accepted complete graph labeling."""
        if not self.accepts(assignment):
            raise ValueError('rejected graph labeling')
        p, m = self.data['par'], self.data['m']
        paths = {}
        for e, path in self.paths.items():
            if path in paths and paths[path] != assignment[e]:
                raise ValueError('inconsistent relay')
            paths[path] = assignment[e]
        def word(path):
            value = sum(paths[path+(j,)]*self.codec.Q**j for j in range(self.codec.width))
            return arith.word(value, p)
        return dict(
            hs={key: [[word(('hs', key, j, k)) for k in range(p['L1'])] for j in range(p['L2'])]
                for key in arith.source_names(m, 'general')},
            outs=[[[word(('outs', eq, j, k)) for k in range(p['L1'])] for j in range(p['L2'])]
                  for eq in range(arith.arithmetic_count(m)+m)],
            quotients=[],
            carry_ports=[{j: [paths[('carry_ports', eq, j, k)] for k in range(p['carry_bits'])]
                          for j in range(p['R'], p['L'] if eq < arith.arithmetic_count(m)
                                         else (p['b']+p['k']-1)//p['k'], p['R'])}
                         for eq in range(arith.arithmetic_count(m)+m)],
            range={key: [0]+[paths[('range', key, j)] for j in range(1, p['n'])]+[1]
                   for key in arith.source_names(m, 'general')},
            nonzero={key: [0]+[paths[('nonzero', key, j)] for j in range(1, p['n'])]+[1]
                     for key in range(2, m+1)})

    def accepts(self, assignment):
        if set(assignment) != set(self.net.edges):
            return False
        ring = sc.Ring(1)
        return all(t.entry(tuple(assignment[e] for e in t.ports), ring) == ring.one
                   for t in self.net.tensors)

    def pin(self, pins):
        if not set(pins) <= set(self.net.edges):
            raise ValueError('unknown edge')
        if any(type(x) is not int or not 0 <= x < self.codec.Q for x in pins.values()):
            raise ValueError('invalid pinned symbol')
        rows, cols = [], []
        ring = sc.Ring(1)
        zero = False
        for original, target in ((self.net.rows, rows), (self.net.cols, cols)):
            for t in original:
                fixed = Pinned(t, pins)
                if fixed.ports:
                    target.append(fixed)
                elif fixed.entry((), ring) != ring.one:
                    zero = True
        if zero:
            # Preserve the free variables but multiply by an explicit zero scalar.
            rows.append(sc.Table('zero-pinned-factor', (), {(): (0,)}, 1, self.codec.alphabet))
        return sc.Network('C38-pinned', rows, cols, [e for e in self.net.edges if e not in pins], self.codec.alphabet)


def recover_by_counts(compiled, pins, count):
    """Recover every unpinned symbol by exact conditional counts."""
    result = dict(pins)
    if count(compiled.pin(result)) == 0:
        return None
    for edge in compiled.net.edges:
        if edge in result:
            continue
        for symbol in range(compiled.codec.Q):
            result[edge] = symbol
            if count(compiled.pin(result)) > 0:
                break
        else:
            raise AssertionError('positive count without extension')
    assert compiled.accepts(result)
    return result
