"""Actual C38 descriptors, relay faults, partial contraction and extraction."""
from collections import Counter
from copy import deepcopy
from itertools import product
import json
from pathlib import Path
import random
import time
from c38_network import Compiled, NonincidentRead, LocalMap, recover_by_counts, arith, sc, common
import adaptive_pair_join as adaptive


def main():
    common.verify_sources()
    rng = random.Random(2026110102)
    counts = Counter()
    cases = []
    fixtures = []
    for b, L, c in ((4, 8, 3), (4, 16, 3), (7, 16, 27)):
        for m in (2, 3, 4):
            data, values, factors = arith.make_instance(b, L, c, m, rng, 'general')
            witness = arith.witness(data, values)
            for t in (1, 2):
                compiled = Compiled(data, t)
                assignment = compiled.encode_witness(witness)
                assert compiled.accepts(assignment)
                assert compiled.decode_witness(assignment) == witness
                assert arith.recover_tuple(data, compiled.decode_witness(assignment)) == factors
                counts['public_data_compilations'] += 1
                counts['incident_predicate_programs_checked'] += sum(hasattr(x, 'kind') for x in compiled.net.tensors)
                counts['explicit_edges_checked'] += len(compiled.net.edges)
                counts['equality_relays_checked'] += sum(x.uid.startswith('relay:') for x in compiled.net.tensors)
                # Fault one endpoint of a relay without changing its other edge.
                for field in ('carry_ports', 'range', 'nonzero'):
                    e = next(e for e, path in compiled.paths.items() if path[0] == field)
                    bad = dict(assignment)
                    bad[e] = 1-bad[e]
                    assert not compiled.accepts(bad)
                    relay = next(x for x in compiled.net.tensors if x.uid.startswith('relay:') and e in x.ports)
                    assert relay.entry(tuple(bad[k] for k in relay.ports), sc.Ring(1)) == (0,)
                    counts['rejected_split_relay_states'] += 1
                # The local reader must reject a hidden read, not return a default.
                try:
                    LocalMap({})['not-an-incident-port']
                except NonincidentRead:
                    counts['nonincident_read_controls'] += 1
                else:
                    raise AssertionError('missing read guard')
                # Compare the explicit graph with independently reconstructed
                # global assignments under coherent (two-endpoint) mutations.
                for _ in range(5):
                    bad = deepcopy(witness)
                    eq = rng.randrange(len(bad['outs']))
                    j = rng.randrange(data['par']['L2'])
                    k = rng.randrange(data['par']['L1'])
                    old = arith.read_word(bad['outs'][eq][j][k], data['par'])
                    bad['outs'][eq][j][k] = arith.word((old+1) % data['par']['modulus'], data['par'])
                    assert compiled.accepts(compiled.encode_witness(bad)) == arith.accept(data, bad) == False
                    counts['coherent_mutation_agreements'] += 1
                if b == 4 and L == 8 and m == 2:
                    fixtures.append((compiled, assignment, factors))
                cases.append(dict(b=b, L=L, m=m, t=t, edges=len(compiled.net.edges),
                                  tensors=len(compiled.net.tensors), categories=dict(compiled.categories)))
    recursion = []
    original = sc.matrix_join.contract_group_tables
    try:
        sc.matrix_join.contract_group_tables = adaptive.contract_group_tables
        for compiled, assignment, factors in fixtures:
            # Different real parts of the arithmetic graph, plus both sides of
            # a relay. Everything else is explicitly pinned in this diagnostic.
            selected = []
            for field in ('hs', 'outs', 'carry_ports', 'range', 'nonzero'):
                selected.append((field, [next(e for e, path in compiled.paths.items() if path[0] == field)]))
            e0 = next(e for e, path in compiled.paths.items() if path[0] == 'range')
            selected.append(('relay_pair', [e for e, path in compiled.paths.items() if path == compiled.paths[e0]]))
            for label, free in selected:
                start = time.monotonic()
                pins = {e: x for e, x in assignment.items() if e not in free}
                net = compiled.pin(pins)
                p = (compiled.codec.Q**len(free)).bit_length()+1
                ring = sc.Ring(p)
                # Integer Boolean-product sum, independently of modular ring
                # accumulation and of the contraction recursion.
                expected = 0
                for symbols in product(range(compiled.codec.Q), repeat=len(free)):
                    trial = dict(pins)
                    trial.update(zip(free, symbols))
                    expected += compiled.accepts(trial)
                    counts['full_graph_completions_enumerated'] += 1
                checked = Counter()
                def observer(child, child_ring, answer, depth):
                    assert answer == child.direct(sc.Ring(child_ring.p, child_ring.axes))
                    checked['active_calls'] += 1
                engine = sc.Engine(compiled.codec.alphabet, max_depth=1, groups=9, observer=observer)
                got = engine.contract(net, p)[0]
                assert got == expected == 1
                count_stats = Counter()
                def count(residual):
                    precision = (compiled.codec.Q**len(residual.edges)).bit_length()+1
                    solver = sc.Engine(compiled.codec.alphabet, max_depth=1, groups=9)
                    value = solver.contract(residual, precision)[0]
                    count_stats['count_calls'] += 1
                    count_stats['active_calls'] += solver.stats['complete_simplex_contractions']
                    return value
                recovered = recover_by_counts(compiled, pins, count)
                assert recovered == assignment
                decoded = compiled.decode_witness(recovered)
                assert arith.recover_tuple(compiled.data, decoded) == factors
                counts['pinned_contractions'] += 1
                counts['pinned_witness_recoveries'] += 1
                counts['independently_checked_active_calls'] += checked['active_calls']
                counts['conditional_count_calls'] += count_stats['count_calls']
                counts['recovery_active_calls'] += count_stats['active_calls']
                recursion.append(dict(t=compiled.codec.alphabet.dimension, part=label, total_edges=len(compiled.net.edges),
                                      free_edges=len(free), pinned_edges=len(pins), expected=expected, got=got,
                                      recursion=dict(engine.stats), recovery=dict(count_stats),
                                      seconds=round(time.monotonic()-start, 3)))
                print(json.dumps({'part': label, 't': compiled.codec.alphabet.dimension, 'status': 'passed',
                                  'seconds': recursion[-1]['seconds']}), flush=True)
            # An inconsistent pin must return no witness through the same API.
            bad = dict(assignment)
            bad[e0] = 1-bad[e0]
            assert recover_by_counts(compiled, bad, count) is None
            counts['inconsistent_pin_rejections'] += 1
    finally:
        sc.matrix_join.contract_group_tables = original
    result = dict(status='passed', counts=dict(counts), cases=cases, recursion=recursion,
                  inherited_source_hashes=common.verify_sources(),
                  scope='Public-data C38 graph compilation, explicit relay equality, incident-only entry access, exact contraction and witness recovery with all but one or two edges pinned. No unpinned relation search, asymptotic packing implementation, or asymptotic runtime validation.')
    (Path(__file__).parent/'network_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': dict(counts)}, indent=2))


if __name__ == '__main__':
    main()
