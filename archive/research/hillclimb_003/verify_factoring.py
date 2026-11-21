"""Check complete C38 graph refactoring and its charged wire ledger."""
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import random
from c38_slices import Compiled, arith, sc, common
from factor_network import Factored, equation_dependencies


def main():
    common.verify_sources()
    rng = random.Random(2026110104)
    counts = Counter()
    records = []
    for b,L,c,m,t in ((4,8,3,2,1),(4,8,3,3,2),(4,16,3,5,1),(7,16,27,8,2)):
        data,values,points = arith.make_instance(b,L,c,m,rng,'general')
        witness = arith.witness(data,values)
        original = Compiled(data,t)
        factored = Factored(original)
        assignment = original.encode_witness(witness)
        lifted = factored.lift(assignment)
        assert original.accepts(assignment) and factored.accepts(lifted)
        assert arith.recover_tuple(data,original.decode_witness(assignment)) == points
        counts['complete_graph_refactorings'] += 1
        # Compare each original predicate with the product of its replacements.
        # Cover all source, middle and closure vertices and many mutations of
        # actual valid incidences, not just arbitrary zero-valued inputs.
        ring = sc.Ring(1)
        replacements = {base.uid:[new for new in factored.net.tensors if new.uid.startswith(base.uid+':')
                                   and not new.uid.startswith('copy:')]
                        for base in original.net.tensors if hasattr(base,'kind')}
        for base in original.net.tensors:
            if not hasattr(base,'kind'):
                continue
            for case in range(9):
                trial = dict(assignment)
                if case:
                    for edge in rng.sample(list(base.ports),min(case,len(base.ports))):
                        trial[edge] = (trial[edge]+1) % original.codec.Q
                lifted_trial = factored.lift(trial)
                expected = base.entry(tuple(trial[e] for e in base.ports),ring)
                actual = ring.one
                for new in replacements[base.uid]:
                    actual = ring.mul(actual,new.entry(tuple(lifted_trial[e] for e in new.ports),ring))
                assert actual == expected, (base.uid,case)
                counts['local_product_equivalences'] += 1
                if not case:
                    assert actual == ring.one
            counts['original_predicates_checked'] += 1
        for _ in range(12):
            trial = dict(assignment)
            edge = rng.choice(original.net.edges)
            trial[edge] = (trial[edge]+1) % original.codec.Q
            assert factored.accepts(factored.lift(trial)) == original.accepts(trial)
            counts['whole_graph_mutation_agreements'] += 1
        # Inconsistent new equality branches have no valid extension.
        for edge in list(factored.new_edges)[:8]:
            bad = dict(lifted)
            bad[edge] = (bad[edge]+1) % original.codec.Q
            assert not factored.accepts(bad)
            counts['rejected_new_wire_inconsistencies'] += 1
        ledger = factored.ledger()
        records.append(dict(b=b,L=L,m=m,t=t,**ledger))
    # The syntactic operation program admits constant dependency bounds.
    # Check a range including both parities and all endpoint special cases.
    dependency_scan = []
    for m in range(2,101):
        deps = equation_dependencies(m)
        fan = {name:sum(name in d for d in deps) for name in arith.source_names(m,'general')}
        assert max(map(len,deps)) <= 8 and max(fan.values()) <= 6
        dependency_scan.append(dict(m=m,equations=len(deps),sources=len(fan),
                                    max_sources=max(map(len,deps)), max_consumers=max(fan.values()),
                                    added_arrays=sum(x for x in fan.values() if x > 1)))
    counts['dependency_arities_checked'] = len(dependency_scan)
    result = dict(status='passed',counts=dict(counts),records=records,dependency_scan=dependency_scan,
                  scope='Full public-data graph refactoring with exact local-product checks, explicit copy-wire bijection, and charged wire ledger. No full unpinned recursive contraction, improved runtime bound, or closure under all recursive projected joins is claimed.')
    (Path(__file__).parent/'factoring_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':dict(counts),'records':records},indent=2))


if __name__ == '__main__':
    main()
