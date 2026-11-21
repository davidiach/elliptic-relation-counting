"""Independent finite sums and real simplex recursion through adaptive joins."""
from collections import Counter
from itertools import product
import random
import common
import adaptive_pair_join as adaptive
import simplex_contraction as sc


def direct(groups, modes, P, limit, table, ring):
    """Plain Cartesian enumeration; no pair-sector or join implementation."""
    incident = {g: tuple(i for i, m in enumerate(modes) if g in (m['source'], m['target'])) for g in groups}
    out = ring.zero
    for assignment in product(range(limit+1), repeat=len(modes)):
        if sum(assignment) > limit or sum(x*m['weight'] for x, m in zip(assignment, modes)) >= P:
            continue
        term = ring.one
        for g in groups:
            term = ring.mul(term, table(g, tuple(assignment[i] for i in incident[g])))
        out = ring.add(out, term)
    return out


def main():
    common.verify_sources()
    counts = Counter()
    rng = random.Random(2026103102)
    def accumulate(stats):
        for key, value in stats.items():
            if key.startswith('peak_'):
                counts[key] = max(counts[key], value)
            else:
                counts[key] += value
    # Perfect-square pair sectors make the tree selector's objective exact.
    # Check its realized partition on dense weights as well as sparse modes.
    for _ in range(50):
        weights = tuple(rng.randrange(20) for _ in range(36))
        _, selected = adaptive.choose_tree(tuple(w*w for w in weights))
        objective = sum(sum(weights[e] for e in bundle)**2 for bundle in adaptive.BUNDLES[selected])
        assert 4*objective <= 27*sum(w*w for w in weights)
        counts['realized_partition_bound_checks'] += 1
    for case in range(48):
        n = 2+case % 8
        groups = tuple(range(n))
        modes = []
        for _ in range(case % 5):
            a, b = rng.sample(groups, 2)
            modes.append(dict(source=a, target=b, weight=rng.randrange(1, 4)))
        P, limit = 3+case % 4, case % 4
        axes = (('old', 2),) if case % 2 else ()
        ring = sc.Ring(5, axes)
        def table(g, assignment):
            v = g+1+sum((j+2)*b+b*b for j, b in enumerate(assignment))
            return ring.val(v) if not axes else (v % ring.N, (g+sum(assignment)) % ring.N)
        expected = direct(groups, modes, P, limit, table, ring)
        stats = Counter()
        got = adaptive.contract_group_tables(groups, modes, P, table, ring, stats, limit)
        assert got == expected, (case, got, expected)
        accumulate(stats)
        counts['independent_mode_sum_cases'] += 1

    # Two equally weighted coordinate copies on each triangle side force
    # nontrivial 2x2 Strassen products in an exact (weight,degree) sector.
    modes = [dict(source=a, target=b, weight=1, coordinate=c)
             for a, b in ((0, 1), (0, 2), (1, 2)) for c in (0, 1)]
    ring = sc.Ring(7, (('old', 2),))
    def table(g, assignment):
        return ((g+1+sum((i+1)*v for i, v in enumerate(assignment))) % ring.N,
                (1+sum(v*v for v in assignment)) % ring.N)
    stats = Counter()
    expected = direct((0, 1, 2), modes, 4, 3, table, ring)
    got = adaptive.contract_group_tables((0, 1, 2), modes, 4, table, ring, stats, 3)
    assert got == expected and stats['strassen_calls'] > 0
    accumulate(stats)
    counts['duplicated_coordinate_nonnumeric_ring_cases'] += 1

    # Run the unchanged simplex engine with the new join adapter in this process.
    # All full recursive calls are checked, including old coefficient axes.
    original = sc.matrix_join.contract_group_tables
    recursion_records = []
    def observer(net, ring, answer, depth):
        assert answer == net.direct(sc.Ring(ring.p, ring.axes))
        counts['active_recursions_checked_directly'] += 1
        counts['checked_depth_'+str(depth)] += 1
    try:
        sc.matrix_join.contract_group_tables = adaptive.contract_group_tables
        for t, depth, axes in ((2, 2, ()), (3, 1, (('jet1', 2),))):
            alphabet = sc.Alphabet(t, t)
            def left(symbols):
                a = symbols[0]+1
                return (a, 1) if axes else a
            def right(symbols):
                a = 2-symbols[0]
                return (a, symbols[0] % 2) if axes else a
            net = sc.network('adaptive_simplex_'+str(t), [([0], left)], [([0], right)],
                             alphabet=alphabet, p=3, axes=axes)
            # Independent root oracle: coefficient arithmetic on the stored
            # canonical integers, without the engine's ring multiplication.
            expected = [0]*(2 if axes else 1)
            for j in range(alphabet.size):
                a, b = net.rows[0].table[(j,)], net.cols[0].table[(j,)]
                expected[0] += a[0]*b[0]
                if axes:
                    expected[1] += a[0]*b[1]+a[1]*b[0]
            expected = tuple(x % 32 for x in expected)
            engine = sc.Engine(alphabet, max_depth=depth, groups=9, observer=observer)
            got = engine.contract(net, 5, axes)
            assert got == expected
            assert engine.stats['pair_sectors'] > 0
            recursion_records.append(dict(t=t, depth=depth, axes=axes, value=got,
                                          independent_expected=expected, stats=dict(engine.stats)))
    finally:
        sc.matrix_join.contract_group_tables = original
    assert counts['checked_depth_1'] > 0
    common.save('adaptive_simplex.json', dict(status='passed', counts=dict(counts), recursion=recursion_records,
        scope='Exact pair-sector join implementation compared with Cartesian sums, duplicated coordinate families, nilpotent coefficients, actual nontrivial Strassen products, and full small simplex recursion. The asymptotic matrix formulas and packing search are not executed.'))


if __name__ == '__main__':
    main()
