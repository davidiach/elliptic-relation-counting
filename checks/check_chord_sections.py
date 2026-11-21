"""Check streamed section summaries against full determinant products."""
from collections import Counter
from itertools import product
from math import prod
from random import Random

from check_chord_difference import literal, normalized_difference
from probe_signed_moment import Probe


def combine(left, right, p):
    h, d = left
    j, e = right
    if h + j >= 2:
        return 2, 0
    if h + j == 0:
        return 0, (d + e) % p
    return left if h else right


def streamed(slots, order, pinned, p, linear):
    """Fix the pin coordinates one section at a time. Keep all determinant rows."""
    if any(not slot for slot in slots):
        return (0, 0), (0, 0, 0, 0)
    result = (0, 0)
    max_side = 0
    sections = 0
    nonzero_empty = 0
    kinds = set()
    for assignment in product(*(slots[j] for j in pinned)):
        section = slots.copy()
        for j, choice in zip(pinned, assignment):
            section[j] = [choice]
        section = [section[j] for j in order] + [section[-1]]
        middle = len(order) // 2
        max_side = max(max_side, prod(map(len, section[:middle])),
                       prod(map(len, section[middle:])))
        h, d = normalized_difference(section, p, linear)
        kinds.add(min(h, 2))
        nonzero_empty += int(h == 0 and d != 0)
        result = combine(result, (min(h, 2), d), p)
        sections += 1
    return result, (max_side, sections, len(kinds), nonzero_empty)


def fixture(probe, m):
    slots = []
    for offset in range(0, probe.p - 1, 2):
        entries = [(j, 2 * digit + sign) for digit in (0, 1)
                   for j, sign in probe.by_x.get(offset + digit, ())]
        if len(entries) == 2:
            slots.append(entries)
        if len(slots) == m:
            return slots
    raise ValueError("insufficient disjoint width-two slots")


def check_cancellation(base, order, pinned, p, linear, reference):
    """Force a global zero by linearity and check that its individual sections disagree."""
    rng = Random(p)
    for _ in range(10):
        first, second = [
            [[(pt, rng.randrange(p)) for pt, _ in slot] for slot in base] for _ in range(2)
        ]
        _, d = normalized_difference(first, p, linear)
        _, e = normalized_difference(second, p, linear)
        if not (d and e):
            continue
        moving = [[(pt, (e*u - d*v) % p) for (pt, u), (_, v) in zip(a, b)]
                  for a, b in zip(first, second)]
        answer, (_, _, _, nonzero) = streamed(moving, order, pinned, p, linear)
        assert answer == (0, 0)
        if nonzero >= 2:
            assert literal(moving, p, linear, 5) == reference
            return nonzero
    raise AssertionError("no example with cancellation between sections")


def check(p, linear, m):
    probe = Probe(p, linear, (0,))
    seed = fixture(probe, m)
    used = {probe.curve[j][0] for slot in seed for j, _ in slot}
    counts = Counter(sum(j for j, _ in a) % probe.n for a in product(*seed))
    targets = {}
    for t in range(1, probe.n):
        if probe.curve[t][0] not in used:
            targets.setdefault(min(counts[t], 2), t)
    order = [0, 3, 4, 1, 2, 5] if m == 6 else [0, 2, 1, 3]
    pinned = [order[0], order[-1]]
    rng = Random(p + m)
    comparisons = Counter()
    mixed_sections = 0
    canceled = 0
    for kind, t in sorted(targets.items()):
        original = [[(probe.curve[j], label) for j, label in slot] for slot in seed]
        original.append([(probe.curve[-t], 0)])
        base = [[(point, 0) for point, _ in slot] for slot in original]
        reference = literal(base, p, linear, 5)
        if kind == 0:
            canceled = check_cancellation(base, order, pinned, p, linear, reference)
        for s in range(12):
            if s < 2 * m:
                moving = [[(point, (label >> (s % 2)) & 1 if i == s // 2 else 0)
                           for point, label in slot] for i, slot in enumerate(original)]
            else:
                moving = [[(point, rng.randrange(p)) for point, _ in slot] for slot in original]
            for variant in (moving, [[(pt, (speed+rng.randrange(p)) % p) for pt, speed in slot]
                                     for slot in moving]):
                (h, digit), (space, batches, ntypes, _) = streamed(variant, order, pinned, p, linear)
                expected_h, expected_digit = normalized_difference(variant, p, linear)
                assert (h, digit) == (min(expected_h, 2), expected_digit)
                assert h == kind and (digit == 0) == (literal(variant, p, linear, 5) == reference)
                assert space == 2 ** (m // 2 - 1) and batches == 4
                mixed_sections += int(ntypes > 1)
                comparisons[kind] += 1
    empty = [[]] + [seed[j] for j in range(1, m)] + [[]]
    assert streamed(empty, order, pinned, p, linear)[0] == (0, 0)
    print(f"F_{p}, {m} width-two slots: {sum(comparisons.values())} literal product comparisons, "
          f"hit classes {sorted(comparisons)}, {mixed_sections} comparisons mixed section classes. "
          f"Peak entries per section side {2 ** (m//2-1)} vs {2 ** (m//2)} in the balanced join.")
    assert canceled >= 2
    print(f"  Additional literal equality from cancellation of {canceled} nonzero empty-section digits.")
    return set(comparisons)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    kinds = check(101, 17, 4) | check(1031, 9, 6) | check(101, 17, 6)
    assert kinds == {0, 1, 2}
