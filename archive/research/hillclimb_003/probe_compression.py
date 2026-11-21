"""Exact root Fourier minors and actual projected-child compression probes."""
from collections import Counter
from itertools import product, combinations
import json
from math import comb
from pathlib import Path
import time
from c38_slices import public_instance, fourier_minor, copy_restriction, sc, common


def binary_rank(rows, columns):
    rows = [sum((x & 1) << j for j,x in enumerate(row)) for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(rows)) if (rows[i] >> column) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank+1, len(rows)):
            if (rows[i] >> column) & 1:
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def sparse_child_value(child, ring, accepted, stats):
    """Exact sum over original-tensor support, independently evaluated first.

    This is a diagnostic enumeration, not a compact-table algorithm.
    """
    out = ring.zero
    for trial in accepted:
        value = ring.one
        for tensor in child.cols:
            symbols = tuple(trial[port[1]] for port in tensor.ports)
            value = ring.mul(value, tensor.entry(symbols, ring))
        out = ring.add(out, value)
        stats['sparse_child_assignments'] += 1
    return out


def projected_case(r):
    start = time.monotonic()
    net, labels, pairs, accepted = copy_restriction(r)
    stats = Counter()
    precision, jet_length = 6*r+5, 2*r+1
    ring = sc.Ring(precision, (('probe', jet_length),))
    packets = sc.geom.packetize(net, labels, 2, 2)
    assert all(set(p['positions']) == {0} for p in packets)
    W, modes, forms = sc.modes_and_forms(packets, 2, 1, ring, 'probe', 1)
    assert len(modes) == 2*r
    mode_for_edge = {}
    for index, mode in enumerate(modes):
        edge = packets[mode['packet']]['positions'][0]
        assert mode['target'] == 0 and mode['weight'] == 3
        assert mode['coeffs'][0] == ring.one
        mode_for_edge[edge] = index
    c_values = [forms[mode['packet'],0]['modes'][index] for index,mode in enumerate(modes)]
    assert all(c == c_values[0] for c in c_values)
    c = c_values[0]
    assert c[0] == 0 and sc.geom.v2(c[1]) == 3
    block = ((ring.one, ring.neg(c)), (ring.neg(c), ring.scale(ring.mul(c,c), 2)))
    engine = sc.Engine(net.alphabet, max_depth=0, groups=3)
    table = {}
    joined = ring.zero
    for beta_bits in product((0,1), repeat=2*r):
        left, right = beta_bits[:r], beta_bits[r:]
        beta = [0]*len(modes)
        for j,(a,b) in enumerate(pairs):
            beta[mode_for_edge[a]], beta[mode_for_edge[b]] = left[j], right[j]
        child = engine.child(net, labels, packets, 0, 0, beta, modes, forms, ring)
        got = sparse_child_value(child, ring, accepted, stats)
        expected = ring.one
        for a,b in zip(left,right):
            expected = ring.mul(expected, block[a][b])
        assert got == expected
        table[left,right] = got
        stats['actual_projected_entries'] += 1
        exterior = ring.one
        for group, bits in ((1,left), (2,right)):
            source_child = engine.child(net, labels, packets, 0, group, beta, modes, forms, ring)
            moment = source_child.direct(ring)
            assert moment == ring.val(1 << (r-sum(bits)))
            exterior = ring.mul(exterior, moment)
            stats['actual_source_children_checked'] += 1
        joined = ring.add(joined, ring.mul(got, exterior))
        if sum(beta_bits) in (0, 2*r) or (r <= 2 and sum(beta_bits) == 1):
            assert got == child.direct(ring)
            stats['full_child_cross_checks'] += 1
    protected = None
    if r % 2 == 0:
        subsets = list(combinations(range(r), r//2))
        matrix = []
        for a in subsets:
            left = tuple(int(j in a) for j in range(r))
            row = []
            for b in subsets:
                right = tuple(int(j in b) for j in range(r))
                coefficient = table[left,right][r]
                assert coefficient % (1 << (3*r)) == 0
                reduced = (coefficient >> (3*r)) & 1
                assert reduced == int(not set(a).intersection(b))
                row.append(reduced)
            matrix.append(row)
        rank = binary_rank(matrix, len(subsets))
        assert rank == comb(r,r//2)
        protected = dict(normalized_binary_rank=rank, valuation=3*r, coefficient_degree=r,
                         precision=precision, degree_cap=2*r, weight=3*r,
                         jet_rank_lower_bound=f'{rank}/{r+1}')
    # The target table alone sums to (1-2*c+2*c^2)^r. In the full
    # three-group join a source moment is 2 for exponent zero and 1
    # for exponent one. Those two exterior tables cannot be dropped.
    enumerated = ring.zero
    for value in table.values():
        enumerated = ring.add(enumerated, value)
    single = ring.add(ring.sub(ring.one, ring.scale(c,2)), ring.scale(ring.mul(c,c),2))
    factored = ring.power(single, r)
    assert factored == enumerated
    full_single = ring.add(ring.sub(ring.val(4), ring.scale(c,4)), ring.scale(ring.mul(c,c),2))
    full_factored = ring.power(full_single,r)
    assert joined == full_factored and joined != enumerated
    if r <= 2:
        assert joined == sc.direct_deformation(net, packets, W, ring)
        stats['full_deformation_cross_checks'] += 1
    return dict(r=r, original_tensor_support=len(accepted), projected_entries=len(table),
                counts=dict(stats), protected_slice=protected, exact_factored_join=True,
                target_table_sum='(1 - 2*c(Z) + 2*c(Z)^2)^r',
                factored_formula='(4 - 4*c(Z) + 2*c(Z)^2)^r',
                rejected_shortcut='Ignoring exterior source moments changes the constant term from 4^r to 1',
                seconds=round(time.monotonic()-start,3))


def main():
    common.verify_sources()
    root_minors = []
    for b,L,p in ((4,8,13), (17,16,131071), (31,64,2147483647)):
        for sources in (1,2):
            record = fourier_minor(public_instance(b,L,p), sources)
            root_minors.append(record)
            print(json.dumps({'root_minor':'passed', **record}), flush=True)
    projections = []
    for r in (1,2,3,4):
        record = projected_case(r)
        projections.append(record)
        print(json.dumps({'projected_case':'passed','r':r,'seconds':record['seconds']}), flush=True)
    abstract = []
    for r in (2,4,6,8,10,16,32,64):
        n = comb(r,r//2)
        # Complementation is an exact permutation of the r/2-subsets.
        # No large matrix needs to be materialized to certify this family.
        abstract.append(dict(r=r, protected_rank=n, jet_rank_lower_bound=f'{n}/{r+1}',
                             minimum_precision=3*r+1, minimum_jet_length=r+1,
                             scope='Algebraic consequence of the proved copy-factor formula'))
    result = dict(status='passed', root_minors=root_minors, actual_projected_cases=projections,
                  symbolic_family=abstract, inherited_source_hashes=common.verify_sources(),
                  scope='Natural Fourier cuts of actual C38 first tensors and actual projected children of explicit local restrictions. High exact rank coexists with a compact factored formula. These finite packet geometries are not the asymptotic balanced schedule; no full C38 projected-table closure or ECDLP improvement is claimed.')
    (Path(__file__).parent/'compression_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','root_minor_cases':len(root_minors),
                      'actual_projected_entries':sum(x['projected_entries'] for x in projections)}))


if __name__ == '__main__':
    main()
