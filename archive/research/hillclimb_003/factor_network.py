"""Exact C38 refactoring into bounded-source arithmetic factors.

The first and closure stages split without new wires. The middle stage
duplicates a source word only when multiple equations consume it, with
explicit equality tensors and every added wire charged.
"""
from collections import Counter, defaultdict
from math import ceil
from c38_slices import arith, sc, source_rows
from c38_network import Binding, Pinned, Predicate, LocalMap


class Dependencies:
    def __init__(self, names=()):
        self.names = frozenset(names)
    def __add__(self, other):
        return Dependencies(self.names | getattr(other, 'names', frozenset()))
    __radd__ = __add__
    __mul__ = __add__
    __rmul__ = __add__
    def __neg__(self):
        return self


def equation_dependencies(m):
    constants = {key:Dependencies() for key in ('p','a','Xm','Ym')}
    constants.update({key+str(j):Dependencies() for key in ('alpha','beta') for j in range(1,m+1)})
    values = {key:Dependencies([key]) for key in arith.source_names(m,'general')}
    result = [value.names for value in arith.operations(values,constants,m,mode='general')]
    assert max(map(len,result)) <= 8
    assert max(sum(name in d for d in result) for name in values) <= 6
    return result


class Equality:
    pc = 1
    axes = ()
    def __init__(self, uid, ports, alphabet):
        self.uid, self.ports, self.alphabet = uid, tuple(ports), alphabet
    def entry(self, symbols, ring):
        return ring.val(int(bool(symbols) and all(type(x) is int and 0 <= x < self.alphabet.size
                                                and x == symbols[0] for x in symbols)))


class EquationMiddle(Predicate):
    def __init__(self, base, eq, dependencies, bindings):
        super().__init__(base.uid+':eq'+str(eq), base.data, base.codec, 'middle', base.index, bindings)
        self.eq, self.dependencies = eq, dependencies

    def entry(self, symbols, ring):
        if len(symbols) != len(self.ports):
            raise ValueError('wrong arity')
        if any(type(x) is not int or not 0 <= x < self.alphabet.size for x in symbols):
            return ring.zero
        leaves = {b.path:x for b,x in zip(self.bindings,symbols)}
        p, m = self.data['par'], self.data['m']
        def word(path):
            value = sum(leaves[path+(j,)]*self.codec.Q**j for j in range(self.codec.width))
            if value >= self.codec.modulus:
                raise ValueError('noncanonical Fourier word')
            return value
        try:
            frequency = {name:arith.dft([word(('hs',name,j,self.index)) for j in range(p['L2'])],
                                       pow(p['omega'],p['L1'],p['modulus']),p['modulus'])
                         for name in self.dependencies}
            public = {key:arith.dft(arith.limbs(value,p),p['omega'],p['modulus'])
                      for key,value in self.data['constants'].items()}
            vals = []
            for t2 in range(p['L2']):
                # Absent source names are provably outside this equation's
                # syntactic dependency set; they receive literal zero.
                source = {name:(frequency[name][t2] if name in frequency else 0)
                          for name in arith.source_names(m,'general')}
                constants = {key:row[self.index+p['L1']*t2] for key,row in public.items()}
                vals.append(arith.operations(source,constants,m,mode='general')[self.eq] % p['modulus'])
            wanted = arith.idft(vals,pow(p['omega'],p['L1'],p['modulus']),p['modulus'])
            actual = [word(('outs',self.eq,j,self.index)) for j in range(p['L2'])]
            return ring.val(int(actual == wanted))
        except ValueError:
            return ring.zero


class EquationClosure(Predicate):
    def __init__(self, base, eq, bindings):
        super().__init__(base.uid+':eq'+str(eq), base.data, base.codec, 'closure', base.index, bindings)
        self.eq = eq

    def entry(self, symbols, ring):
        if len(symbols) != len(self.ports):
            raise ValueError('wrong arity')
        if any(type(x) is not int or not 0 <= x < self.alphabet.size for x in symbols):
            return ring.zero
        leaves = {b.path:x for b,x in zip(self.bindings,symbols)}
        p, m = self.data['par'], self.data['m']
        try:
            raw = {}
            for j2 in range(self.index*p['R'], (self.index+1)*p['R']):
                row = []
                for k in range(p['L1']):
                    value = sum(leaves['outs',self.eq,j2,k,d]*self.codec.Q**d for d in range(self.codec.width))
                    if value >= self.codec.modulus:
                        raise ValueError('noncanonical Fourier word')
                    row.append(value)
                for j1,x in enumerate(arith.H_inverse_row(row,j2,p)):
                    raw[j1*p['L2']+j2] = x if x <= p['modulus']//2 else x-p['modulus']
            arithmetic = self.eq < arith.arithmetic_count(m)
            stop = p['L'] if arithmetic else (p['b']+p['k']-1)//p['k']
            def carry_at(t):
                bits = [leaves['carry_ports',self.eq,t,k] for k in range(p['carry_bits'])]
                if any(x not in (0,1) for x in bits):
                    raise ValueError('nonbinary carry')
                return arith.signed_read(bits,p['carry_bound'],p['carry_bits'])
            for start in range(self.index*p['R'],stop,p['L2']):
                carry = 0 if start == 0 else carry_at(start)
                end = min(start+p['R'],stop)
                for t in range(start,end):
                    value = raw[t]+carry
                    if arithmetic:
                        assert value % p['B'] == 0
                        carry = value//p['B']
                    else:
                        carry,digit = divmod(value,p['B'])
                        for bit in range(p['k']):
                            if self.data['hash_ell'] <= t*p['k']+bit < p['b']:
                                assert not ((digit >> bit) & 1)
                    assert abs(carry) < p['carry_bound']
                if end < stop:
                    assert carry == carry_at(end)
                elif arithmetic:
                    assert carry == 0
            return ring.one
        except (ValueError,AssertionError):
            return ring.zero


class Factored:
    def __init__(self, compiled):
        self.compiled = compiled
        d, p = compiled.data, compiled.data['par']
        m = d['m']
        deps = equation_dependencies(m)
        self.new_edges = {}
        tensors = []
        next_edge = len(compiled.net.edges)
        for base in compiled.net.tensors:
            if not hasattr(base,'kind'):
                tensors.append(base)
                continue
            if base.kind == 'first':
                _, dummy, _ = source_rows(compiled,row=base.index)
                pairs = [(('x1','x2'),2)] + [((f'X{j-1}',f'x{j}'),j) for j in range(3,m+1)]
                used = {name for names,j in pairs for name in names}
                groups = pairs+[((name,),None) for name in arith.source_names(m,'general') if name not in used]
                parts = []
                for number,(names,state) in enumerate(groups):
                    selected = {b.edge for b in base.bindings
                                if ((b.path[0] in ('hs','range') and b.path[1] in names)
                                    or (b.path[0] == 'nonzero' and b.path[1] == state))}
                    part = Pinned(base,{e:v for e,v in dummy.items() if e not in selected})
                    part.uid = base.uid+':part'+str(number)
                    assert set(part.ports) == selected
                    parts.append(part)
                assert Counter(e for t in parts for e in t.ports) == Counter(base.ports)
                tensors += parts
            elif base.kind == 'closure':
                for eq in range(arith.arithmetic_count(m)+m):
                    bindings = [b for b in base.bindings if b.path[1] == eq]
                    tensors.append(EquationClosure(base,eq,bindings))
            else:
                input_bindings = defaultdict(list)
                for binding in base.bindings:
                    if binding.path[0] != 'hs':
                        continue
                    consumers = [eq for eq,names in enumerate(deps) if binding.path[1] in names]
                    assert consumers
                    if len(consumers) == 1:
                        input_bindings[consumers[0]].append(binding)
                    else:
                        ports = [binding.edge]
                        for eq in consumers:
                            edge = next_edge
                            next_edge += 1
                            ports.append(edge)
                            self.new_edges[edge] = binding.edge
                            input_bindings[eq].append(Binding(edge,binding.path))
                        tensors.append(Equality('copy:'+str(binding.edge),ports,compiled.codec.alphabet))
                for eq,names in enumerate(deps):
                    bindings = input_bindings[eq]+[b for b in base.bindings if b.path[0] == 'outs' and b.path[1] == eq]
                    tensors.append(EquationMiddle(base,eq,names,bindings))
        endpoints = defaultdict(list)
        for tensor in tensors:
            for edge in tensor.ports:
                endpoints[edge].append(tensor.uid)
        assert all(len(ends) == 2 and ends[0] != ends[1] for ends in endpoints.values())
        # Explicit port incidence validates the graph in O(S); avoid the legacy
        # mixed() constructor's scan of all tensors for every edge.
        net = sc.Network.__new__(sc.Network)
        net.name='C38-factored'; net.rows=(); net.cols=()
        net.tensors=tuple(tensors); net.edges=tuple(sorted(endpoints)); net.alphabet=compiled.codec.alphabet
        assert len({t.uid for t in tensors}) == len(tensors)
        net.tail={e:v[0] for e,v in endpoints.items()}; net.head={e:v[1] for e,v in endpoints.items()}
        self.net = net
        self.dependencies = deps

    def lift(self, old):
        result = dict(old)
        result.update({edge:old[source] for edge,source in self.new_edges.items()})
        return result

    def accepts(self, assignment):
        if set(assignment) != set(self.net.edges):
            return False
        ring = sc.Ring(1)
        return all(t.entry(tuple(assignment[e] for e in t.ports),ring) == ring.one for t in self.net.tensors)

    def ledger(self):
        multiplicity = Counter(tuple(sorted((self.net.tail[e],self.net.head[e]))) for e in self.net.edges)
        p = self.compiled.data['par']
        fan = {name:sum(name in names for names in self.dependencies)
               for name in arith.source_names(self.compiled.data['m'],'general')}
        expected_extra = p['L']*self.compiled.codec.width*sum(v for v in fan.values() if v > 1)
        assert len(self.new_edges) == expected_extra
        return dict(original_edges=len(self.compiled.net.edges), factored_edges=len(self.net.edges),
                    added_equality_edges=len(self.new_edges), original_max_degree=max(len(t.ports) for t in self.compiled.net.tensors),
                    factored_max_degree=max(len(t.ports) for t in self.net.tensors),
                    max_parallel_multiplicity=max(multiplicity.values()), vertices=len(self.net.tensors),
                    max_equation_sources=max(map(len,self.dependencies)), max_source_consumers=max(fan.values()))
