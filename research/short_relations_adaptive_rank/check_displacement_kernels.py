"""Exact finite cross-mean ranks on the order-19 curve. Run with assertions on."""
from collections import Counter
from itertools import combinations, combinations_with_replacement, product

from check_short_relations import curve13, rank, short_fixture


R = 19
CLASSES = tuple(range(1, 10))  # representatives up to sign


def affine_block(cycle, parity_maps):
    block = []
    for index, (x, y) in enumerate(cycle[1:], start=1):
        a, b = parity_maps[y % 2]
        if (a * x + b) % 13 < 3:
            block.append(index)
    return tuple(block)


def transition_rows(blocks, shift, prime=R):
    """Signed incidences before pin averaging; use the known fixture indices."""
    starts = []
    n = 0
    for block in blocks:
        starts.append(n)
        n += len(block)
    tuples = []
    for digits in product(*(range(len(block)) for block in blocks)):
        value = sum(block[j] for block, j in zip(blocks, digits)) % prime
        columns = [start + j for start, j in zip(starts, digits)]
        tuples.append((value, columns))
    rows = []
    for x, source in tuples:
        for y, target in tuples:
            if (y - x) % prime == shift:
                row = [0] * n
                for i, j in zip(source, target):
                    row[i] -= 1
                    row[j] += 1
                rows.append(row)
    return rows


def contrasts(rows):
    return [[x - y for x, y in zip(row, rows[0])] for row in rows[1:]]


def score(row, witness):
    return sum(x * y for x, y in zip(row, witness)) % R


def check_witness(cycle, add, blocks, shifts, witness, slopes, expected_ranks):
    f = short_fixture(R, blocks, cycle, add)
    ranks = [rank(rows, R, f['n']) for rows in f['by_t']]
    rows = [row for s in shifts for row in f['by_t'][s]]
    stacked = rank(rows, R, f['n'])
    complete = rank(f['rows'], R, f['n'])
    assert ([ranks[s] for s in shifts], stacked, complete) == expected_ranks
    assert any(score(row, witness) for row in f['rows'])
    blind = {s for s, matrix in enumerate(f['by_t'])
             if all(score(row, witness) == 0 for row in matrix)}
    assert blind == {0, *shifts, *(R-s for s in shifts)}
    for s, slope in zip(shifts, slopes):
        for t in (s, R-s):
            assert all(score(row, witness) == 0 for row in f['by_t'][t])
        assert {score(row, witness) for row in transition_rows(blocks, s)} == {slope}
    print(f'Shifts {shifts}, ranks {[ranks[s] for s in shifts]}, '
          f'stack {stacked}/{complete}, raw score increments {slopes}.')
    return f


def named_examples(cycle, add):
    blocks = [(2, 7, 9, 14, 15), (4, 6, 9, 11, 18)]
    check_witness(cycle, add, blocks, (7,), [0,0,0,0,1,0,0,0,0,0], (0,), ([6],6,7))
    blocks = [(3, 6, 9, 14, 17), (5, 8, 11, 12, 16)]
    check_witness(cycle, add, blocks, (3,8), [0,0,0,0,0,0,0,0,1,0], (0,0), ([6,6],6,7))

    examples = [
        ([((1,0), (1,0)), ((2,6), (2,6))],
         (1,2), [0,0,0,0,1,1,0,0], (0,0), ([4,3],4,5)),
        ([((7,3), (5,3)), ((9,2), (6,11))],
         (2,3,9), [1,0,18,0,1,0,0], (18,18,0), ([2,2,3],3,4)),
    ]
    for masks, shifts, witness, slopes, expected in examples:
        blocks = [affine_block(cycle, pair) for pair in masks]
        f = check_witness(cycle, add, blocks, shifts, witness, slopes, expected)
        masses = [count for table in f['H'].values() for count, _ in table.values()]
        assert (min(masses), max(masses)) == (1, 6)
        assert len(f['d']) == (14 if len(shifts) == 2 else 10)
        if len(shifts) == 3:
            assert all(sum(a*s for a, s in zip(signs, shifts)) % R
                       for signs in product((-1, 1), repeat=3))


def modular_cancellation():
    prime = 7
    blocks = [(1, 3), (2, 4, 5), (2, 4, 5)]
    witness = [6, 0, 3, 6, 0, 1, 0, 0]
    f = short_fixture(prime, blocks)
    t = 1
    raw = transition_rows(blocks, t, prime)
    assert len(f['d']) == prime
    assert rank(contrasts(raw), prime, len(witness)) == 4
    assert rank(f['by_t'][t], prime, len(witness)) == 3
    assert {sum(x*y for x, y in zip(row, witness)) % prime for row in raw} == set(range(prime))
    for count, incidence in f['H'][t].values():
        assert (sum(x*y for x, y in zip(incidence, witness)) - 4*count) % prime == 0
    assert any(count % prime for count, _ in f['H'][t].values())
    print('C7: full tuple support, raw contrast rank 4, pin-mean rank 3; '
          'all seven raw scores occur, every centered pin sum is zero.')


def punctured_support():
    prime = 23
    blocks = [(1, 5, 6, 7, 8, 20), (8, 11, 14, 15, 18, 22), (2, 5, 11, 17, 19, 20)]
    m = len(blocks)
    for i in range(m):
        retained = blocks[:i] + blocks[i+1:]
        support = {sum(t) % prime for t in product(*retained)}
        assert len(support) == 22 and 6*m*(prime - len(support)) < prime
    n = sum(map(len, blocks))
    for s in range(prime):
        rows = transition_rows(blocks, s, prime)
        assert rank(contrasts(rows), prime, n) == n-m-1
    print('C23: all three punctured supports have size 22; raw rank 14 at all 23 shifts.')


def sweep(cycle, add):
    blocks = {affine_block(cycle, [(a,b)] * 2)
              for a, b in product(range(1,13), range(13))}
    blocks.discard(())
    accepted = deficient = raw_losses = 0
    failed_lists = Counter()
    for pair in combinations_with_replacement(sorted(blocks), 2):
        f = short_fixture(R, pair, cycle, add)
        if 12 * f['missing'] >= R or 2 * len(f['d']) <= R:
            continue
        accepted += 1
        n, matrices = f['n'], f['by_t']
        target = n - 3
        assert rank(f['rows'], R, n) == target
        # A failing list has only failing subsets. Full rank persists under stacking.
        failures = set()
        for s in CLASSES:
            rs = rank(matrices[s], R, n)
            assert rank(matrices[s] + matrices[R-s], R, n) == rs
            assert rank(matrices[R-s], R, n) == rs
            if rs < target:
                failures.add((s,))
                deficient += 1
                raw_rank = rank(contrasts(transition_rows(pair, s)), R, n)
                assert rs <= raw_rank <= target
                raw_losses += raw_rank > rs
        for length in range(2, 5):
            next_failures = set()
            for subset in combinations(CLASSES, length):
                if not all(sub in failures for sub in combinations(subset, length-1)):
                    continue
                rows = [row for s in subset for row in matrices[s]]
                if rank(rows, R, n) < target:
                    next_failures.add(subset)
            failed_lists.update(next_failures)
            failures = next_failures
    assert (len(blocks), accepted, deficient, raw_losses) == (60, 996, 636, 0)
    bad_pairs = sum(len(s) == 2 for s in failed_lists)
    bad_triples = [s for s in failed_lists if len(s) == 3]
    assert (bad_pairs, len(bad_triples)) == (36, 15)
    assert not any(len(s) == 4 for s in failed_lists)
    assert all(a + b == c or a + b + c == R for a, b, c in bad_triples)
    print(f'Symmetric masks: {accepted} accepted block pairs; {deficient} deficient '
          f'nonzero sign classes; {raw_losses} with further rank loss under pin averaging.')
    print('Fixed lists failing on at least one accepted pair: 36/36 pairs, '
          '15/84 triples, 0/126 quadruples.')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions required')
    curve, addition = curve13()
    named_examples(curve, addition)
    modular_cancellation()
    punctured_support()
    sweep(curve, addition)
