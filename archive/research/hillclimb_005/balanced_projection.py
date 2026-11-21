"""Actual C38 copy factors inside an exactly balanced, fully occupied packing.

The test is a finite node identity, not the full asymptotic schedule.
All padding tensors are explicit pins; no old entry oracle is assumed.
"""
from collections import Counter
from itertools import product
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent/'hillclimb_003'))
from c38_slices import copy_restriction, Ones, sc, common


class Pin:
    pc = 1
    axes = ()
    def __init__(self, uid, edge, alphabet):
        self.uid, self.ports, self.alphabet = uid, (edge,), alphabet
    def entry(self, symbols, ring):
        return ring.one if symbols == (0,) else ring.zero


def fixture(g=3, H=8, a=3, head_word=(1,0)):
    old, _, pairs, accepted = copy_restriction(2)
    actual = old.rows[0]
    live = [e for pair in pairs for e in pair]
    rows, cols = [actual], []
    labels = {actual.uid:(0, 0)}
    edges = []
    fresh = max(live)+1
    classes = Counter()
    words = list(product(range(g), repeat=2))
    for left in words:
        for right in words:
            for k in range(H):
                if left == (0, 0) and right == head_word and k < len(live):
                    e = live[k]
                    tail = actual
                    head = Ones('freeX'+str(e), [e], old.alphabet)
                else:
                    e, fresh = fresh, fresh+1
                    tail = Pin('pinY'+str(e), e, old.alphabet)
                    head = Pin('pinX'+str(e), e, old.alphabet)
                    rows.append(tail)
                cols.append(head)
                labels[tail.uid], labels[head.uid] = left, right
                edges.append(e)
                classes[left, right] += 1
    assert len(classes) == g**4 and set(classes.values()) == {H}
    net = sc.Network('balanced-C38-copy-packing', rows, cols, edges, old.alphabet)
    packets = sc.geom.packetize(net, labels, a, 2)
    if a == 3 and H == 8:
        assert len(packets) == len(edges)//8
        assert all(set(p['positions']) == set(range(8)) for p in packets)
    pi = next(i for i, p in enumerate(packets) if set(live) <= set(p['positions'].values()))
    positions = {e:x for x,e in packets[pi]['positions'].items()}
    assert set(packets[pi]['positions']) == set(range(1<<a))
    assert [positions[e] for e in live] == ([0, 6, 4, 2] if a==3 else [0, 12, 8, 4])
    return net, labels, packets, pi, pairs, accepted


def polynomial_table(pairs, positions, forms, pi, ids, ring, fixed_rows=(), degree_limit=4):
    n = len(ids)
    zero = (0,)*n
    def add(a, b):
        out = dict(a)
        for exponent, value in b.items():
            out[exponent] = ring.add(out.get(exponent, ring.zero), value)
        return {e:v for e,v in out.items() if v != ring.zero}
    def scale(a, scalar):
        return {e:ring.scale(v, scalar) for e,v in a.items()}
    def multiply(a, b):
        out = {}
        for ea, va in a.items():
            for eb, vb in b.items():
                e = tuple(x+y for x,y in zip(ea, eb))
                if sum(e) > degree_limit:
                    continue
                out[e] = ring.add(out.get(e, ring.zero), ring.mul(va, vb))
        return {e:v for e,v in out.items() if v != ring.zero}
    def linear(x):
        return {tuple(int(j == k) for j in range(n)):forms[pi,x]['modes'].get(mode, ring.zero)
                for k, mode in enumerate(ids)}
    out = {zero:ring.one}
    for a,b in pairs:
        za, zb = linear(positions[a]), linear(positions[b])
        factor = add({zero:ring.one}, scale(add(za, zb), -1))
        factor = add(factor, scale(multiply(za, zb), 2))
        out = multiply(out, factor)
    for x in fixed_rows:
        out = multiply(out,add({zero:ring.one},scale(linear(x),-1)))
    return out


def direct_sparse_node(packets, pi, live, pairs, accepted, W, ring):
    """Independent 64-term sum of original Y/X assignments.

    Every omitted Y and X is pinned to zero. Other packets therefore
    contribute exactly one, with no coefficient or mode assumptions.
    """
    out = ring.zero
    for bits in product((0, 1), repeat=len(live)):
        X = dict(zip(live, bits))
        coordinates = {}
        for x,e in packets[pi]['positions'].items():
            value = ring.zero
            for y,f in packets[pi]['positions'].items():
                if X.get(f, 0):
                    value = ring.add(value, W[y-x])
            coordinates[e] = value
        for Y in accepted:
            term = ring.one
            for e, value in coordinates.items():
                term = ring.mul(term, value if Y.get(e, 0) else ring.sub(ring.one, value))
            out = ring.add(out, term)
    return out


def run_case(a=3, H=8, head_word=(1,0), precision=96, jet_length=9):
    start = time.monotonic()
    net, labels, packets, pi, pairs, accepted = fixture(H=H,a=a,head_word=head_word)
    live = [e for pair in pairs for e in pair]
    positions = {e:x for x,e in packets[pi]['positions'].items()}
    ring = sc.Ring(precision, (('balanced', jet_length),))
    W, modes, forms = sc.modes_and_forms(packets, a, 1, ring, 'balanced', 1)
    # All nonzero source variables occur in these four modes. Every other
    # mode has all its source ports pinned to zero, so positive exponents
    # vanish. This justifies reducing the full global join to these modes.
    ids = [j for j,m in enumerate(modes)
           if any(packets[m['packet']]['positions'][y] in live for y in m['coeffs'])]
    assert all(modes[j]['packet']==pi and modes[j]['target']==0 for j in ids)
    assert len(ids) == (4 if a==3 else 7)
    fixed_rows = [x for x,e in packets[pi]['positions'].items() if x%2==0 and e not in live]
    degree_limit = min(4+len(fixed_rows),jet_length-1)
    target_table = polynomial_table(pairs, positions, forms, pi, ids, ring,fixed_rows,degree_limit)
    engine = sc.Engine(net.alphabet, max_depth=0, groups=3)
    full, binary = ring.zero, ring.zero
    records = []
    checked_entries = source_terms = all_packet_checks = 0
    special={tuple([n]+[0]*(len(ids)-1)) for n in (0,2,4)}
    for wanted in product(range(degree_limit+1), repeat=len(ids)):
        if sum(wanted) > degree_limit:
            continue
        beta = [0]*len(modes)
        for j,b in zip(ids,wanted):
            beta[j] = b
        assert sum(m['weight']*b for m,b in zip(modes,beta)) < ring.p
        values = []
        for group in (0,1):
            child = engine.child(net, labels, packets, 0, group, beta, modes, forms, ring)
            packet = child.cols[pi]
            assignments = accepted if group == 0 else [dict(zip(live,bits)) for bits in product((0,1),repeat=4)]
            total = ring.zero
            for assignment in assignments:
                symbols = tuple(assignment.get(e,0) for side,e in packet.ports)
                value = packet.entry(symbols,ring)
                total = ring.add(total,value)
                if group == 1:
                    source_terms += 1
                if wanted in special:
                    for tensor in child.rows:
                        assert tensor.entry(tuple(assignment.get(e,0) for side,e in tensor.ports),ring) == ring.one
                    for other, tensor in enumerate(child.cols):
                        if other != pi:
                            assert tensor.entry((0,)*len(tensor.ports),ring) == ring.one
                            all_packet_checks += 1
            values.append(total)
        expected = target_table.get(wanted,ring.zero)
        assert values[0] == expected
        checked_entries += 1
        term = ring.mul(*values)
        full = ring.add(full,term)
        if max(wanted) <= 1:
            binary = ring.add(binary,term)
        if wanted in special and sum(wanted) in (2,4):
            assert expected != ring.zero
            first = next(i for i,x in enumerate(expected) if x)
            assert first == sum(wanted)
            if a==3:
                c1=F(-144,7)
                leading=c1**first*(F(36,35) if first==2 else F(4,105))
                assert expected[first]==ring.val(leading)[0]
            records.append(dict(exponents=wanted, leading_degree=first,
                                leading_valuation=sc.geom.v2(expected[first]),
                                coefficient=list(expected)))
    direct = direct_sparse_node(packets,pi,live,pairs,accepted,W,ring)
    assert checked_entries==comb(len(ids)+degree_limit,degree_limit)
    assert direct == full
    assert full != binary
    difference = ring.sub(full,binary)
    first = next(i for i,x in enumerate(difference) if x)
    result = dict(status='passed', edges=len(net.edges), vertices=len(net.tensors),
                  ordered_label_classes=81, load_per_class=H, packet_size=1<<a,
                  packets=len(packets), full_packets=sum(len(p['positions'])==1<<a for p in packets),
                  selected_packet_key=packets[pi]['key'], selected_packet_fully_occupied=True,
                  actual_C38_copy_pairs=2, live_positions=[positions[e] for e in live],
                  selected_mode_weights=[modes[j]['weight'] for j in ids],
                  original_global_count=4, target_entries_checked=checked_entries,
                  source_assignments_evaluated=source_terms, unit_padding_packet_checks=all_packet_checks,
                  fixed_target_rows=fixed_rows, precision=precision, jet_length=jet_length,
                  full_join_matches_independent_64_term_deformation=True,
                  repeated_mode_coefficients=records,
                  binary_mode_truncation_changes_full_join=dict(first_degree=first,valuation=sc.geom.v2(difference[first])),
                  compact_formula='product_pairs (1-u_a-u_b+2*u_a*u_b) * product_fixed_rows (1-u_x), with u_e=sum_j A[e,j]*z_j',
                  scope='Exactly balanced finite three-group, two-coordinate packing with a full selected packet and a genuine restricted C38 F0 factor. Not the full C38 network or asymptotic schedule; no lower bound against all symbolic algorithms.',
                  seconds=round(time.monotonic()-start,3))
    print(json.dumps({k:v for k,v in result.items() if k!='repeated_mode_coefficients'},indent=2))
    return result


def main():
    common.verify_sources()
    cases=[run_case(),run_case(a=4,H=4,head_word=(1,1),precision=128,jet_length=5)]
    result=dict(status='passed',cases=cases)
    (ROOT/'balanced_results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
