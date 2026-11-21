"""Complete small m=2 source searches versus independent curve-pair searches.

This is exhaustive arithmetic enumeration, not full-network contraction.
It includes independently chosen hashes, uncovered targets, and every
candidate source triple in F_13^3 for each fixture.
"""
from collections import Counter
from copy import deepcopy
from itertools import product
import json
from pathlib import Path
import random
from c38_network import Compiled, arith, common


def add_points(a, b, p):
    # Independent implementation, including identity and tangent branches.
    if a is None:
        return b
    if b is None:
        return a
    x, y = a
    u, v = b
    if x == u and (y+v) % p == 0:
        return None
    if x == u:
        slope = ((3*x*x+1)*pow(2*y, p-2, p)) % p
    else:
        slope = ((v-y)*pow(u-x, p-2, p)) % p
    xx = (slope*slope-x-u) % p
    return xx, (slope*(x-xx)-y) % p


def main():
    common.verify_sources()
    rng = random.Random(2026110103)
    p = 13
    points = [(x, y) for x in range(p) for y in range(p) if (y*y-x*x*x-x-1) % p == 0]
    template, _, _ = arith.make_instance(4, 8, 3, 2, rng, 'general')
    counts = Counter()
    records = []
    for fixture in range(12):
        data = deepcopy(template)
        target = points[(fixture*5) % len(points)]
        data['constants']['Xm'], data['constants']['Ym'] = target
        # These hashes are chosen without a planted decomposition.
        for slot in (1, 2):
            data['constants']['alpha'+str(slot)] = rng.randrange(1, 16, 2)
            data['constants']['beta'+str(slot)] = rng.randrange(16)
        data['hash_ell'] = 2 if fixture < 8 else 4
        def member(slot, x):
            c = data['constants']
            return (c['alpha'+str(slot)]*x+c['beta'+str(slot)]) % 16 < (1 << data['hash_ell'])
        expected = set()
        for a, b in product(points, repeat=2):
            counts['independent_curve_pairs_examined'] += 1
            if a[0] != b[0] and member(1, a[0]) and member(2, b[0]) and add_points(a, b, p) == target:
                expected.add((a, b))
        compiled = Compiled(data, 2)
        found = set()
        witnesses = set()
        for x, y, slope in product(range(p), repeat=3):
            counts['source_triples_examined'] += 1
            if x == y or not member(1, x) or not member(2, y):
                continue
            values = dict(x1=x, x2=y, l2=slope)
            residuals = arith.operations(values, data['constants'], 2)[:2]
            if any(v % p for v in residuals):
                continue
            for i, v in enumerate(residuals):
                values['q'+str(i)] = v//p+8*p
                assert 0 <= values['q'+str(i)] < 16*p
            original = arith.witness(data, values)
            encoded = compiled.encode_witness(original)
            assert compiled.accepts(encoded)
            decoded = compiled.decode_witness(encoded)
            recovered = tuple(arith.recover_tuple(data, decoded))
            assert recovered not in found
            assert tuple(encoded.values()) not in witnesses
            found.add(recovered)
            witnesses.add(tuple(encoded.values()))
            counts['complete_graph_witnesses_checked'] += 1
        assert found == expected, (fixture, found, expected)
        counts['complete_source_searches'] += 1
        counts['uncovered_targets'] += not found
        counts['covered_targets'] += bool(found)
        records.append(dict(fixture=fixture, target=target, hash_bits=data['hash_ell'],
                            hashes={k:v for k,v in data['constants'].items() if k.startswith(('alpha','beta'))},
                            witnesses=len(found)))
    assert counts['covered_targets'] and counts['uncovered_targets']
    result = dict(status='passed', counts=dict(counts), records=records,
                  scope='Every source triple for 12 independent-hash m=2 F_13 instances compared with every curve-point pair. All found witnesses pass the full explicit graph. Enumeration is used, not unpinned recursive contraction.')
    (Path(__file__).parent/'source_enumeration_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':dict(counts)}, indent=2))


if __name__ == '__main__':
    main()
