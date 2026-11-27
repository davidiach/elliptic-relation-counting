"""Verify affine-mask witnesses for consecutive targets with excess pinned nullity."""

import json
from fractions import Fraction
from itertools import combinations, product
from math import isqrt, prod
from pathlib import Path

import numpy as np
from flint import nmod_mat

from research.deterministic_endgame_barriers.check_pinning_barriers import curve_cycle


def affine_block(cycle, p, threshold, maps):
    logs = []
    for k, (x, y) in enumerate(cycle[1:], 1):
        alpha, beta = maps[y % 2]
        assert 0 < alpha < p and 0 <= beta < p
        if (alpha*x + beta) % p < threshold:
            logs.append(k)
    return np.array(logs, dtype=np.int64)


def sum_counts(blocks, r):
    counts = np.zeros(r, dtype=np.int64)
    counts[0] = 1
    for block in blocks:
        counts = sum(np.roll(counts, int(k)) for k in block)
    return counts


def pinned_matrices(blocks, r):
    # After pinning one slot, count completions with another slot fixed.
    indices = [i for i in range(len(blocks)) for _ in range(2)]
    singles, pairs = [], {}
    for i, xs in enumerate(blocks):
        rest = indices.copy()
        rest.remove(i)
        singles.append(sum_counts([blocks[k] for k in rest], r))
        for j in range(i, len(blocks)):
            others = rest.copy()
            others.remove(j)
            pairs[i,j] = (sum_counts([blocks[k] for k in others], r),
                          np.add.outer(xs, blocks[j]))

    def at(t):
        rows = [[None] * len(blocks) for _ in blocks]
        for (i,j), (ways, sums) in pairs.items():
            part = ways[(t-sums) % r]
            if i == j:
                part[np.diag_indices(len(blocks[i]))] += singles[i][(t-blocks[i]) % r]
            else:
                part *= 2
            rows[i][j], rows[j][i] = part, part.T
        return np.block(rows)
    return at


def rank(matrix, r):
    return nmod_mat(matrix.tolist(), r).rank()


def verify_two_point_contrast(source, blocks, shift, matrices, step, expected_counts):
    """Check the antisymmetric mode on two tagged translates of the same pair."""
    r = len(matrices)
    pair = source[0]
    assert len(pair) == 2 and np.array_equal(pair, source[1])
    v = np.zeros(sum(map(len, blocks)), dtype=np.int64)
    for x, sign in zip(pair, (1, -1)):
        v[np.flatnonzero(blocks[0] == (x+shift) % r).item()] = sign
        v[2 + np.flatnonzero(blocks[1] == x).item()] = -sign
    # Removing two pins from the four slots on the duplicate blocks leaves
    # two copies of the half-list source[1:]. Use integer cyclic convolution.
    remaining = sum_counts([*source[1:], *source[1:]], r)
    weights = np.roll(remaining, (2*shift + int(pair.sum())) % r)
    for weight, matrix in zip(weights, matrices):
        assert np.array_equal(matrix @ v, 2*weight*v)
    at_run = [int(weights[k*step % r]) for k in range(len(expected_counts))]
    assert at_run == expected_counts and all(w % r == 0 for w in at_run)
    print(f'shared contrast completion counts={at_run}')


def enumerate_target(t, buckets, logs, block_sizes, r):
    """Join all ordered half-tuples and independently count the pinned rows."""
    m, n = len(block_sizes), len(logs)
    total = sum(len(left)*len(buckets[(t-s) % r]) for s, left in enumerate(buckets))
    sample_indices = np.linspace(0, total-1, min(16*n, total), dtype=int)
    sample = np.empty((len(sample_indices), 2*m), dtype=np.int64)
    matrix = np.zeros(n*n, dtype=np.int64)
    pins = np.zeros(n, dtype=np.int64)
    offset = 0
    for s, left in enumerate(buckets):
        right = buckets[(t-s) % r]
        if not len(left) or not len(right):
            continue
        slots = np.column_stack((np.repeat(left, len(right), axis=0),
                                 np.tile(right, (len(left), 1))))
        assert np.all(logs[slots].sum(axis=1) % r == t)
        columns = np.ascontiguousarray(slots.T)
        for pin in columns[:m]:
            pins += np.bincount(pin, minlength=n)
            row_index = pin*n
            for slot in columns:
                matrix += np.bincount(row_index+slot, minlength=n*n)
        lo, hi = np.searchsorted(sample_indices, [offset, offset+len(slots)])
        sample[lo:hi] = slots[sample_indices[lo:hi]-offset]
        offset += len(slots)
    assert offset == total

    # A certified lower bound on the full incidence rank from genuine rows.
    # The m-1 block contrasts and the shifted logarithm give the upper bound.
    offsets = np.cumsum([0, *block_sizes])
    forced = np.zeros((n, m), dtype=np.int64)
    for j in range(m-1):
        forced[offsets[j]:offsets[j+1], j] = 1
        forced[offsets[-2]:, j] = -1
    forced[:, -1] = (2*m*logs-t) % r
    assert rank(forced, r) == m
    rows = np.zeros((len(sample), n), dtype=np.int64)
    for slot in sample.T:
        rows[np.arange(len(sample)), slot] += 1
    assert not np.any(rows @ forced % r)
    assert rank(rows, r) == n-m
    return matrix.reshape(n,n), total, int(np.count_nonzero(pins % r == 0))


def verify(case):
    p, a, b, r = (case[k] for k in ('p','a','b','r'))
    assert all(z > 1 and all(z % j for j in range(2,isqrt(z)+1)) for z in (p,r))
    cycle = curve_cycle(p, a, b)
    assert len(cycle) == r
    m = len(case['seed'])
    q = next(k for k in range(1,r) if k**m >= 16*r)
    threshold = (p*q + r-2)//(r-1)
    assert threshold == case['threshold']
    blocks = [affine_block(cycle,p,threshold,maps) for maps in case['seed']]
    sizes = list(map(len,blocks))
    n = sum(sizes)
    assert sizes == case['sizes']
    length = len(case['nullities'])
    step = cycle.index(tuple(case['P']))
    assert step and 1 <= length < r
    # Translating one block by H shifts every full target by 2H.
    # The new run has midpoint (length-1)P/2, so 4H=(length-1)P.
    shift = (length-1)*step*pow(4,-1,r) % r
    old = affine_block(cycle,p,threshold,[case['source_first_map']]*2)
    source = [old, *blocks[1:]]
    assert all(set(v) == set(-v % r) for v in source)
    assert set((old+shift) % r) == set(blocks[0])

    offsets = np.cumsum([0, *sizes])
    half = np.array(list(product(*(range(i,j) for i,j in zip(offsets,offsets[1:])))))
    logs = np.concatenate(blocks)
    sums = logs[half].sum(axis=1) % r
    hist = np.bincount(sums,minlength=r)
    assert np.array_equal(hist, np.roll(sum_counts(source,r),shift))
    buckets = np.split(half[np.argsort(sums)], np.cumsum(hist)[:-1])
    w = prod(sizes)
    assert w*w * (2*m)**2 < 2**63
    energy = int(hist @ hist)
    mu = Fraction((r-1)*threshold,p)**m
    phi = (r*energy + 64*(w-mu)**2)/mu**2
    assert phi == Fraction(case['phi']) and phi < Fraction(64,57)
    assert np.all(hist > 0)
    print(f'r={r}, m={m}, T={threshold}, sizes={sizes}, P={cycle[step]}')
    print(f'half-tuples={w}, energy={energy}, half-sum support={len(hist)}, Phi={phi}')

    matrix = pinned_matrices(blocks,r)
    matrices = [matrix(t) for t in range(r)]
    nullities = [n-rank(mat,r) for mat in matrices]
    # The shared-mode witnesses also have documented exceptions outside the run.
    expected_bad = set(case['exceptional_logs'])
    for t, value in enumerate(nullities):
        if t not in expected_bad:
            assert value == m, ("ordinary target nullity", t, value, m)
    exceptional = [t for t,v in enumerate(nullities) if v>m]
    assert exceptional == case['exceptional_logs']
    sampled = []
    for k in range(length):
        t = k*step % r
        direct, total, missing = enumerate_target(t,buckets,logs,sizes,r)
        assert np.array_equal(matrices[t],direct)
        assert nullities[t] == case['nullities'][k] > m
        assert missing == case['zero_pins_mod_r'][k]
        sampled.append(direct)
        print(f'k={k}, log={t}, tuples={total}, incidence rank={n-m}, '
              f'pinned rank={n-nullities[t]}, nullity={nullities[t]}, zero pins mod r={missing}')
    joint = n-rank(np.vstack(sampled),r)
    assert joint == case['joint_nullity']
    next_nullity = nullities[length*step % r]
    assert next_nullity == case['next_nullity']
    print(f'exceptional logs={exceptional}, joint nullity={joint}, nullity at {length}P={next_nullity}')
    if 'contrast_counts' in case:
        verify_two_point_contrast(source, blocks, shift, matrices, step, case['contrast_counts'])
    if 'every_pair_nullity' in case:
        pair_nullity = case['every_pair_nullity']
        for a, b in combinations(matrices, 2):
            assert n-rank(np.vstack((a,b)),r) == pair_nullity
        print(f'all {r*(r-1)//2} distinct target pairs have joint nullity {pair_nullity}')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run with assertions enabled.')
    for fixture in json.loads(Path(__file__).with_name('witnesses.json').read_text()):
        verify(fixture)
