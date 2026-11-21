"""Exhaust the coarse nonlinear join score on three actual width-one slots."""
from itertools import product
from math import comb, prod, sqrt

from check_address_parity import character
from probe_signed_moment import Probe


def majority(values):
    return 1 if sum(values) >= 0 else -1


def majority_coefficient(size):
    return comb(size-1, (size-1)//2)/2**(size-1)


def signed_word(value, mask):
    return -1 if (value & mask).bit_count() % 2 else 1


def fixture():
    probe = Probe(101, 17, ())
    valid_x = [x for x in range(probe.p) if len(probe.by_x[x]) == 2]
    slots = [tuple(probe.by_x[x]) for x in valid_x[:3]]
    target, marked = next((t, w) for t in range(1, probe.n)
                          if (w := probe.singleton(slots, t)[0]) is not None)
    left = [(point, address << 2) for point, address in slots[0]]
    right = [((target-p-q) % probe.n, (a << 1) | b)
             for (p, a), (q, b) in product(slots[1], slots[2])]
    return probe, target, marked, left, right


def group_buckets(left, right, mask):
    labels = [[[], []] for _ in range(2)]
    for side, half in enumerate((left, right)):
        for index, (point, _) in enumerate(half):
            labels[(point & mask).bit_count() % 2][side].append(index)
    return [(u, v) for u, v in labels if u and v]


def check_score():
    probe, target, witness, left, right = fixture()
    assert len(left) == 2 and len(right) == 4
    counts = []
    tests = []
    for mask in range(1 << probe.n.bit_length()):
        buckets = group_buckets(left, right, mask)
        count = sum(len(u)*len(v) for u, v in buckets)
        counts.append(count)
        if len(buckets) == 2 and count > 2:
            tests.append((mask, buckets, count))
    mean = sum(counts)/len(counts)
    assert mean == 1+(len(left)*len(right)-1)/2
    print(f"p={probe.p}, N={probe.n}, guarded singleton T=[{target}]G, "
          f"witness address {witness}, {len(left)} left and {len(right)} right tuples.")
    print(f"  All {len(counts)} binary hashes: mean coarse matches {mean:.3f}.")
    mask, buckets, count = tests[0]
    masses = [len(u)*len(v)/(2*count) for u, v in buckets]
    edges = [(left[i][1] | right[j][1], b, i, j) for b, (u, v) in enumerate(buckets)
             for i in u for j in v]
    assert len({a for a, _, _, _ in edges}) == count
    assert witness in {a for a, _, _, _ in edges}
    aligned = [0.0]*count
    abs_values = [0.0]*count
    bits = len(left)+len(right)+len(buckets)
    words = 1 << bits
    challenges = probe.size
    for signs in range(words):
        u_sign = [signed_word(signs, 1 << i) for i in range(len(left))]
        v_sign = [signed_word(signs, 1 << (len(left)+j)) for j in range(len(right))]
        z_sign = [signed_word(signs, 1 << (len(left)+len(right)+b)) for b in range(len(buckets))]
        scores = []
        for r in range(challenges):
            factors = []
            for b, (u, v) in enumerate(buckets):
                left_vote = majority(u_sign[i]*character(left[i][1], r) for i in u)
                right_vote = majority(v_sign[j]*character(right[j][1], r) for j in v)
                factors.append(complex(sqrt(1-masses[b]),
                                       sqrt(masses[b])*z_sign[b]*left_vote*right_vote))
            scores.append(prod(factors).imag)
        assert max(map(abs, scores)) <= 1+1e-12
        coeffs = [sum(value*character(address, r) for r, value in enumerate(scores))/challenges
                  for address, _, _, _ in edges]
        assert sum(c*c for c in coeffs) <= 1+1e-12
        for e, ((_, b, i, j), coefficient) in enumerate(zip(edges, coeffs)):
            aligned[e] += z_sign[b]*u_sign[i]*v_sign[j]*coefficient/words
            abs_values[e] += abs(coefficient)/words
    for e, (_, b, _, _) in enumerate(edges):
        u, v = buckets[b]
        expected = (sqrt(masses[b])*prod(sqrt(1-x) for other, x in enumerate(masses) if other != b)
                    * majority_coefficient(len(u))*majority_coefficient(len(v)))
        assert abs(aligned[e]-expected) < 1e-10
        assert abs_values[e] >= aligned[e]-1e-10 >= 1/(4*sqrt(count))-1e-10
    assert sum(abs_values)/count <= 1/sqrt(count)+1e-10
    print(f"  Mask {mask}, {count} candidate addresses, {len(buckets)} occupied buckets: "
          f"all {words} coin strings, {challenges} challenges and all candidate alignments checked.")
    print(f"  Minimum mean absolute coefficient {min(abs_values):.6f}, "
          f"guaranteed floor {1/(4*sqrt(count)):.6f}, minimax ceiling {1/sqrt(count):.6f}.")


def check_prefix_occupancy():
    # The two exclusive label sets have distinct keys. Enumerating their
    # random bucket values suffices for this occupancy identity.
    buckets = 4
    assignments = buckets**(2*buckets)
    shared = 0
    for labels in product(range(buckets), repeat=2*buckets):
        shared += len(set(labels[:buckets]) & set(labels[buckets:]))
    expected = buckets*(1-(1-1/buckets)**buckets)**2
    assert shared/assignments == expected
    assert expected >= buckets/4
    print(f"  Four exclusive labels per side: all {assignments} hash assignments "
          f"gave mean {expected:.6f} jointly occupied prefixes out of four.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_score()
    check_prefix_occupancy()
