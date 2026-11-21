"""Verify projected half-sum joins, intermediate ranks and signed zero counts."""
from collections import Counter
from itertools import product
from random import Random

from flint import nmod_mat

from fixtures import determinant, evaluation
from probe_signed_moment import Probe


def split_sums(slots, radix, order):
    halves = []
    for choices in product(*slots):
        address = 0
        for _, digit in choices:
            address = address*radix+digit
        halves.append((sum(point for point, _ in choices) % order, address))
    return halves


def main():
    probe = Probe(101, 17, (0,))
    p, rng, m = probe.p, Random(7861), 4
    # One admissible affine mask (multiplier 1) and three disjoint fixed-bit masks.
    slots = []
    for offset in range(0, p-1, 2):
        slot = [(i, 2*j+a) for j in (0, 1) for i, a in probe.by_x.get(offset+j, ())]
        if len(slot) == 2:
            slots.append(slot)
        if len(slots) == m:
            break
    used = {probe.curve[i][0] for slot in slots for i, _ in slot}
    actual = Counter(sum(i for i, _ in c) % probe.n for c in product(*slots))
    left = split_sums(slots[:2], probe.block, probe.n)
    tail = split_sums(slots[2:], probe.block, probe.n)
    off_curve = next((0, y) for y in range(3) if y*y % p != 1)
    coordinates = [off_curve, *probe.curve[1:]]
    assert len(set(coordinates)) == probe.n

    for count in (0, 1):
        target = next(t for t in range(1, probe.n) if actual[t] == count and probe.curve[t][0] not in used)
        right = [((target-i) % probe.n, address) for i, address in tail]
        forbidden = set()
        for (i, _), (j, _) in product(left, right):
            if i == j:
                continue
            x, y = coordinates[i]
            u, v = coordinates[j]
            forbidden.add(p if y == v else (u-x)*pow(y-v, -1, p) % p)
        best = None
        for c in range(p+1):
            labels = [(y if c == p else x+c*y) % p for x, y in coordinates]
            exact = all((i == j) == (labels[i] == labels[j]) for (i, _), (j, _) in product(left, right))
            assert exact == (c not in forbidden)
            size = min(len({labels[i] for i, _ in half}) for half in (left, right))
            if exact and (best is None or best[0] < size):
                best = size, labels
        assert len(forbidden) <= len(left)*len(right) and best
        s, labels = best
        a, b = (sorted({labels[i] for i, _ in half}) for half in (left, right))
        assert len(set(a) & set(b)) == count and s == 4
        def rank(entry, first=a):
            return nmod_mat([[entry(x-y) % p for y in b] for x in first], p).rank()
        assert rank(lambda d: int(d == 0)) == count
        for j in range(s):
            coefficients = [rng.randrange(p) for _ in range(j)]+[rng.randrange(1, p)]
            assert rank(lambda d: sum(v*pow(d, k, p) for k, v in enumerate(coefficients))) == j+1
        disjoint = [x for x in a if x not in b]
        assert rank(lambda d: pow(d % p, -1, p), disjoint) == min(len(disjoint), len(b))

        endpoint = evaluation(probe.curve[-target], m+1, p)
        hits = []
        for choices in product(*slots):
            rows = [evaluation(probe.curve[i], m+1, p) for i, _ in choices]+[endpoint]
            if determinant(rows, p) == 0:
                digits = [a for _, a in choices]
                hits.append(sum(digit*probe.block**(m-1-k) for k, digit in enumerate(digits)))
        assert len(hits) == count
        for r in range(probe.block**m):
            signed = sum((-1)**((r & address).bit_count()) for address in hits)
            weights = []
            for digits, half in zip(divmod(r, probe.block**2), (left, right)):
                table = Counter()
                for i, address in half:
                    table[labels[i]] += (-1)**((digits & address).bit_count())
                weights.append(table)
            result = sum(w*weights[1][x] for x, w in weights[0].items())
            assert result == signed
        print(f"F_{p}, {count} witness: {len(forbidden)}/{p+1} directions produce a false match; "
              f"all {probe.block**m} signed joins agree with the determinant test.")

    # Larger label lists separately exercise the full rank profile and the inverse shortcut.
    p, s = 1031, 28
    sample = rng.sample(range(p), 2*s-1)
    a, b = sample[:s], sample[s-1:]
    for j in range(s):
        matrix = nmod_mat([[pow((x-y) % p, j, p) for y in b] for x in a], p)
        assert matrix.rank() == j+1
    matrix = nmod_mat([[pow((x-y) % p, -1, p) for y in b] for x in a[:-1]], p)
    assert matrix.rank() == s-1
    print(f"F_{p}, {s} labels on each side, one match: polynomial ranks 1 through {s}, disjoint reciprocal rank {s-1}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    main()
