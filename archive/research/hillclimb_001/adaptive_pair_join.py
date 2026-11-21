"""Exact pair-sector contraction with a selected nine-group join tree.

Tree selection uses rational square-root approximations. Algebraic evaluation
is exact. Strassen is the finite backend; published rectangular asymptotic
formulas are NOT implemented by this module.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import isqrt
import common
import weighted_matrix_join as matrix

EDGES = tuple(combinations(range(9), 2))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
BETA = F(37134, 67867)
D = 1-F(32133, 100000)*BETA


def partitions(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for pair in combinations(vertices[1:], 2):
        block = (first,)+pair
        rest = tuple(v for v in vertices if v not in block)
        for tail in partitions(rest):
            yield (block,)+tail


PARTITIONS = tuple(partitions(tuple(range(9))))
assert len(PARTITIONS) == 280


def cross(left, right):
    return tuple(EDGE_INDEX[tuple(sorted((a, b)))] for a in left for b in right)


BUNDLES = tuple(tuple(cross(p[i], p[j]) for i, j in ((0, 1), (0, 2), (1, 2)))
                for p in PARTITIONS)


def rect(dimensions):
    a, b, c = sorted(dimensions)
    return max(b+c, c+D*b+BETA*a)


@lru_cache(maxsize=4096)
def choose_tree(pair_weights, root_bits=32):
    if len(pair_weights) != 36 or any(type(x) is not int or x < 0 for x in pair_weights):
        raise ValueError('36 nonnegative integer pair weights required')
    if root_bits < 1:
        raise ValueError('root precision must be positive')
    # Scaling by 2**root_bits is common, so integer roots suffice for ordering.
    w = tuple(isqrt(x << (2*root_bits)) for x in pair_weights)
    scores = [sum(sum(w[e] for e in bundle)**2 for bundle in bundles) for bundles in BUNDLES]
    selected = min(range(280), key=lambda i: (scores[i], i))
    partition = PARTITIONS[selected]

    def dimensions(left, right):
        rest = set(range(9))-set(left)-set(right)
        return tuple(sum(w[e] for e in cross(a, b))
                     for a, b in ((left, right), (left, rest), (right, rest)))

    trees = []
    for triple in partition:
        options = []
        for isolated in triple:
            pair = tuple(v for v in triple if v != isolated)
            cost = max(rect(dimensions((pair[0],), (pair[1],))),
                       rect(dimensions(pair, (isolated,))))
            options.append((cost, isolated, (pair, isolated)))
        trees.append(min(options)[2])
    return ((trees[0], trees[1]), trees[2]), selected


def contract_group_tables(groups, modes, P, value, ring, stats=None,
                          degree_limit=None, multiply=matrix.rectangular):
    """Sum exactly over W(beta)<P and |beta|<=degree_limit.

Pair sectors include exact weight and degree. The latter adds only a fixed
polynomial number of sectors and avoids mixing incompatible degree budgets.
"""
    groups = tuple(groups)
    if len(groups) > 9 or len(set(groups)) != len(groups):
        raise ValueError('at most nine distinct groups required')
    if degree_limit is None:
        degree_limit = P-1
    if P <= 0 or degree_limit < 0:
        return ring.zero
    if not groups:
        if modes:
            raise ValueError('modes without groups')
        return ring.one
    stats = stats if stats is not None else Counter()
    lookup = {g: i for i, g in enumerate(groups)}
    indexed = []
    pair_ids = defaultdict(list)
    for i, mode in enumerate(modes):
        if mode['weight'] <= 0 or mode['source'] == mode['target']:
            raise ValueError('invalid mode')
        ends = tuple(sorted((lookup[mode['source']], lookup[mode['target']])))
        pair = EDGE_INDEX[ends]
        indexed.append(ends)
        pair_ids[pair].append(i)
    incident = {g: tuple(i for i, ends in enumerate(indexed) if g in ends) for g in range(9)}
    options = {}
    for pair, ids in sorted(pair_ids.items()):
        sectors = defaultdict(list)
        for assignment, weight, degree in matrix.assignments(tuple(ids), modes, P, degree_limit):
            sectors[(weight, degree)].append(assignment)
        options[pair] = tuple(sorted(sectors.items()))
    pairs = tuple(options)
    cached_values = {}

    def projected(g, assignment):
        if g >= len(groups):
            return ring.one
        key = (g, assignment)
        if key not in cached_values:
            cached_values[key] = value(groups[g], assignment)
            stats['adaptive_projected_values'] += 1
        else:
            stats['adaptive_projected_cache_hits'] += 1
        return cached_values[key]

    chosen = {}
    weights = [0]*36

    def sectors(j, weight, degree):
        if j == len(pairs):
            yield dict(chosen), tuple(weights)
            return
        pair = pairs[j]
        for (wp, dp), assignments in options[pair]:
            if weight+wp >= P or degree+dp > degree_limit:
                continue
            chosen[pair] = assignments
            weights[pair] = wp
            yield from sectors(j+1, weight+wp, degree+dp)

    output = ring.zero
    for selected, sector_weights in sectors(0, 0, 0):
        tree, partition_index = choose_tree(sector_weights, max(32, P.bit_length()*2))
        stats['pair_sectors'] += 1
        stats['selected_partition_'+str(partition_index)] += 1

        @lru_cache(maxsize=None)
        def assignments_on(ids):
            if not ids:
                return ((),)
            used = tuple(pair for pair in pairs if pair_ids[pair][0] in ids)
            assert set(ids) == {i for pair in used for i in pair_ids[pair]}
            answers = []
            for parts in product(*(selected[pair] for pair in used)):
                mapping = {i: b for pair, part in zip(used, parts) for i, b in zip(pair_ids[pair], part)}
                answers.append(tuple(mapping[i] for i in ids))
            return tuple(answers)

        def combine(ids, xids, yids, x, y):
            mapping = dict(zip(xids, x))
            mapping.update(zip(yids, y))
            return tuple(mapping[i] for i in ids)

        def evaluate(node):
            if isinstance(node, int):
                ids = incident[node]
                return {node}, ids, {a: projected(node, a) for a in assignments_on(ids)}
            left, right = node
            lg, li, lv = evaluate(left)
            rg, ri, rv = evaluate(right)
            assert not lg & rg
            X = tuple(i for i in li if i in ri)
            Y = tuple(i for i in li if i not in X)
            Z = tuple(i for i in ri if i not in X)
            boundary = tuple(sorted(Y+Z))
            xs, ys, zs = assignments_on(X), assignments_on(Y), assignments_on(Z)
            A = [[lv[combine(li, X, Y, x, y)] for x in xs] for y in ys]
            B = [[rv[combine(ri, X, Z, x, z)] for z in zs] for x in xs]
            C = multiply(A, B, ring, stats)
            table = {combine(boundary, Y, Z, y, z): C[i][j]
                     for i, y in enumerate(ys) for j, z in enumerate(zs)}
            stats['adaptive_tree_joins'] += 1
            stats['peak_adaptive_partial_entries'] = max(stats['peak_adaptive_partial_entries'], len(table))
            return lg | rg, boundary, table

        _, ids, table = evaluate(tree)
        assert not ids and set(table) == {()}
        output = ring.add(output, table[()])
    return output
