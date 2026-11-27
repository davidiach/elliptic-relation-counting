"""Exact rank regressions on small, explicitly indexed cyclic sets."""
from collections import Counter
from itertools import combinations, product

from check_displacement_kernels import contrasts, transition_rows
from check_short_relations import rank, short_fixture


def positive_pin_tables(r, m):
    # Each free slot contributes a nonzero group element. Reflection permutes
    # the allowed entries, so it has the same count array as a positive slot.
    powers = [[1] + [0] * (r - 1)]
    for _ in range(2*m - 1):
        previous = powers[-1]
        powers.append([sum(previous[(s-x) % r] for x in range(1, r)) for s in range(r)])
    one_fixed, two_fixed = powers[2*m-1], powers[2*m-2]
    labels = [(i, a) for i in range(m) for a in range(1, r)]
    tables = []
    for t in range(r):
        pins = []
        for i, a in labels:
            mass = one_fixed[(t-a) % r]
            row = []
            for j, b in labels:
                if j == i:
                    target = mass if b == a else 0
                else:
                    target = two_fixed[(t-a-b) % r]
                source = two_fixed[(t-a+b) % r]
                row.append(target-source)
            pins.append((mass, row))
        tables.append(pins)
    return tables


def cross_rows(pins, r):
    ref = next(((n, v) for n, v in pins if n % r), None)
    if ref is None:
        return []
    mass, vector = ref
    return [[(mass*x-n*y) % r for x, y in zip(v, vector)] for n, v in pins]


def invisible(rows, u, r):
    return all(sum(x*y for x, y in zip(row, u)) % r == 0 for row in rows)


def classification(r, m, direct=False):
    positive = positive_pin_tables(r, m)
    n = m*(r-1)
    blocks = [tuple(range(1, r))] * m
    fixture = short_fixture(r, blocks) if direct else None
    matrices = []
    for t in range(r):
        pins = positive[t] + [(mass, [-x for x in row]) for mass, row in positive[-t % r]]
        if fixture is not None:
            for j in range(n):
                for sign, idx in ((1, j), (-1, j+n)):
                    assert tuple(fixture['H'][t][sign, j]) == pins[idx]
        matrices.append(cross_rows(pins, r))

    k = (2*m-1) % r
    raw = n-m-1
    h = (r-1)//2
    loss = h-1 if k == 0 else m-1 if k == 1 else 0
    for t in range(1, r):
        assert rank(matrices[t], r, n) == raw-loss
        if k == 0:
            odd = [pow(a, 3, r) for _ in blocks for a in range(1, r)]
            assert invisible(matrices[t], odd, r)
        if k == 1 and m >= 2:
            parity = [1 if (a * pow(t, -1, r) % r) % 2 else -1 for a in range(1, r)]
            witness = parity + [-x for x in parity] + [0]*((m-2)*(r-1))
            assert invisible(matrices[t], witness, r)
    if k == 0:
        assert not matrices[0]
    pair_loss = h-1 if k == 0 else 0
    for s, t in combinations(range(1, h+1), 2):
        assert rank(matrices[s]+matrices[t], r, n) == raw-pair_loss
    union = rank((row for mat in matrices for row in mat), r, n)
    assert union == raw-pair_loss
    print(f'C{r}, m={m}, k={k}: raw target {raw}, rank at every nonzero shift '
          f'{raw-loss}, complete stacked rank {union}.')


def unequal_progressions():
    prime = 11
    blocks = [(2, 4, 6), (2, 4, 6, 8, 10, 1, 3)]
    f = short_fixture(prime, blocks)
    start = 0
    for block in blocks:
        for j, a in enumerate(block[:-1]):
            assert (block[j+1]-a) % prime == 2
            minus = f['H'][prime-1][1, start+j][0]
            plus = f['H'][1][1, start+j+1][0]
            assert minus == plus and minus % prime
        start += len(block)
    assert rank(f['by_t'][1], prime, start) == 7
    print('C11 unequal progressions: every consecutive edge has nonzero mass, rank 7/7.')


def pointwise_no_wrap():
    r = 101
    blocks = [(1, 1+d) for d in (1, 2, 4, 8, 16, 22, 50)]
    f = short_fixture(r, blocks)
    assert len(f['d']) == r and Counter(f['d'].values()) == {1: 74, 2: 27}
    masses = [mass for table in f['H'].values() for mass, _ in table.values()]
    assert (min(masses), max(masses)) == (64, 92)
    u = [0,0,78,0,9,0,1,0,2,0,1,0,0,0]
    ranks = [rank(rows, r, len(u)) for rows in f['by_t']]
    blind = {25, 26, 75, 76}
    assert all(value == (5 if t in blind else 6) for t, value in enumerate(ranks))
    assert rank(f['by_t'][25]+f['by_t'][26], r, len(u)) == 5
    assert {t for t, rows in enumerate(f['by_t']) if invisible(rows, u, r)} == blind
    for t in (25, 26):
        raw = transition_rows(blocks, t, r)
        scores = {sum(x*y for x, y in zip(row, u)) % r for row in raw}
        assert len(raw) == 169 and len(scores) > 1 and 0 not in scores
        assert rank(contrasts(raw), r, len(u)) == 6
        assert all(sum(x*y for x, y in zip(row, u)) % r == 0
                   for _, row in f['H'][t].values())
    print('C101: tuple counts 1..2, every pin mass 64..92; at shifts 25 and 26 '
          'all 169 scores per shift are nonzero and every pinned sum is zero. '
          'Raw ranks 6, individual and joint cross-mean ranks 5.')


def dense_binary_examples():
    # All coordinates are known cyclic fixture labels. Give the two slots of
    # each binary choice as 1 and 1+step, keeping every label nonidentity.
    examples = [
        (7, (1, 1, 2, 2, 3), (1, 2), [6,0,6,0,1,0,1,0,0,0], [3,2], 3, 4),
        (17, (1,1,1,3,3,7,7), (2, 3), [11,0,11,0,11,0,16,0,16,0,1,0,0,0], [2,5], 5, 6),
    ]
    for r, steps, shifts, witness, expected_ranks, joint, target in examples:
        blocks = [(1, 1+step) for step in steps]
        f = short_fixture(r, blocks)
        n = f['n']
        counts = list(f['d'].values())
        assert len(counts) == r and 0 < min(counts) <= max(counts) < r
        punctures = []
        for i in range(len(blocks)):
            d = Counter(sum(t) % r for t in product(*(blocks[:i]+blocks[i+1:])))
            assert len(d) == r
            punctures.extend(d.values())
        assert max(punctures) < r
        assert [rank(f['by_t'][s], r, n) for s in shifts] == expected_ranks
        assert rank([row for s in shifts for row in f['by_t'][s]], r, n) == joint
        assert rank(f['rows'], r, n) == target
        blind = {s for s, mat in enumerate(f['by_t']) if invisible(mat, witness, r)}
        assert blind == {*shifts, *(r-s for s in shifts)}
        for s in shifts:
            raw = contrasts(transition_rows(blocks, s, r))
            assert rank(raw, r, n) == target and not invisible(raw, witness, r)
        if r == 17:
            for s in (*shifts, *(r-s for s in shifts)):
                assert all(mass % r for mass, _ in f['H'][s].values())
            assert rank(f['by_t'][1], r, n) == target
        else:
            def sign_class(x):
                return min(x % r, -x % r)
            dilates = {tuple(sorted(sign_class(a*s) for s in shifts)) for a in range(1, r)}
            assert dilates == set(combinations((1, 2, 3), 2))
        print(f'C{r} binary blocks: tuple counts {min(counts)}..{max(counts)}, '
              f'punctured counts {min(punctures)}..{max(punctures)}, raw rank {target}, '
              f'shifts {shifts} have ranks {expected_ranks} and joint rank {joint}.')


def no_wrap_example():
    r = 5
    blocks = [(1, 4), (1, 4), (1, 2, 3, 4)]
    f = short_fixture(r, blocks)
    assert [f['d'][g] for g in range(r)] == [2, 4, 3, 3, 4]
    for i in range(3):
        a, b = blocks[:i]+blocks[i+1:]
        sums = {(x+y) % r for x in a for y in b}
        assert len(sums) == (3 if i == 2 else 5)
    u = [4, 0, 4, 0, 1, 1, 0, 0]
    for t in range(r):
        assert any(mass % r for mass, _ in f['H'][t].values())
        assert rank(f['by_t'][t], r, 8) == 3
        assert invisible(f['by_t'][t], u, r) == (t != 0)
        if t:
            rows = contrasts(transition_rows(blocks, t, r))
            assert rank(rows, r, 8) == 4
            assert not invisible(rows, u, r)
    nonzero = [row for t in range(1, r) for row in f['by_t'][t]]
    assert rank(nonzero, r, 8) == 3 and rank(f['rows'], r, 8) == 4
    print('C5 unequal blocks: tuple counts 2..4; all nonzero shifts together '
          'give rank 3, adding the identity restores rank 4.')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions required')
    for prime, slots, direct in [(5, 1, True), (5, 2, True), (5, 3, True), (5, 5, False), (7, 3, True),
                                 (7, 4, False), (5, 6, False), (7, 8, False),
                                 (11, 2, True), (11, 6, False)]:
        classification(prime, slots, direct)
    unequal_progressions()
    pointwise_no_wrap()
    dense_binary_examples()
    no_wrap_example()
