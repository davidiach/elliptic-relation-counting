"""Exact GCDs of full-rank minors of the unweighted two-pin stacks.

Each fixture specifies the shift modulus (None means equality over the integers),
the number of occurrences per block per tuple, the blocks, and the expected GCD.
All Gram entries and row operations use integers.
"""
from collections import Counter, defaultdict
from itertools import product
from math import prod


FIXTURES = (
    (None, 1, [(0, 1)] * 3, 1),
    (None, 1, [(0, 1, 2)] * 3, 1),
    (None, 1, [(0, 1, 3), (0, 1, 2), (0, 2, 3)], 1),
    (None, 1, [(0, 1, 2, 3)] * 3, 1),
    (None, 1, [(0, 1)] * 4, 1),
    (None, 1, [(0, 1, 2)] * 4, 27),
    (None, 2, [(0, 1)] * 2, 1),
    (None, 2, [(0, 1, 2)] * 2, 3),
    (None, 2, [(0, 1, 3), (0, 1, 2)], 1),
    (None, 2, [(0, 1, 2, 3)] * 2, 1),
    (None, 2, [(0, 1, 2, 4), (0, 1, 3, 4)], 1),
    (7, 1, [(1, 2, 3)] * 3, 1),
    (7, 1, [(1, 2, 4)] * 3, 1),
    (7, 1, [(1, 2, 3), (1, 2, 4), (1, 3, 5)], 1),
    (11, 2, [(1, 2, 3)] * 2, 3),
    (11, 2, [(1, 2, 3, 4)] * 2, 1),
    (11, 2, [(1, 3, 4, 5, 9)] * 2, 1),
    (11, 1, [(1, 2, 3), (1, 2, 4), (1, 3, 6)], 1),
)


def two_pin_rows(modulus, copies, blocks):
    # In each incidence difference, a block's last coordinate is minus the sum
    # of the others. Its row and column in each Gram are thus redundant over Z.
    width = sum(len(block) - 1 for block in blocks)
    slot_blocks = [block for block in blocks for _ in range(copies)]
    fibers = defaultdict(list)
    for entries in product(*slot_blocks):
        value = sum(entries)
        if modulus is not None:
            value %= modulus
        counts = []
        pins = []
        for i, block in enumerate(blocks):
            selected = entries[i * copies:(i + 1) * copies]
            counts.extend(selected.count(label) for label in block[:-1])
            # Identically distributed repeated slots give the same pinned Gram.
            pins.append((i, selected[0]))
        fibers[value].append((counts, pins))

    grams = defaultdict(lambda: [[0] * width for _ in range(width)])
    for value, sources in fibers.items():
        successor = value + 1 if modulus is None else (value + 1) % modulus
        for (x, source_pins), (y, target_pins) in product(sources, fibers.get(successor, ())):
            delta = [(i, b - a) for i, (a, b) in enumerate(zip(x, y)) if a != b]
            for pin_pair in product(source_pins, target_pins):
                gram = grams[pin_pair]
                for i, di in delta:
                    for j, dj in delta:
                        gram[i][j] += di * dj
    return [row for gram in grams.values() for row in gram]


def maximal_minor_gcd(rows, width):
    """GCD of all width-by-width minors, including zero for insufficient rank."""
    rows = [row.copy() for row in rows if any(row)]
    pivots = []
    for col in range(width):
        pivot = next((i for i in range(col, len(rows)) if rows[i][col]), None)
        if pivot is None:
            return 0
        rows[col], rows[pivot] = rows[pivot], rows[col]
        for i in range(col + 1, len(rows)):
            # Euclidean row subtraction and swap preserve the ideal of minors.
            while rows[i][col]:
                top, lower = rows[col], rows[i]
                quotient = top[col] // lower[col]
                rows[col], rows[i] = lower, [a - quotient * b for a, b in zip(top, lower)]
        pivots.append(rows[col][col])
    # Only the first width rows survive, forming an upper triangular matrix.
    return abs(prod(pivots))


def main():
    counts = Counter()
    for modulus, copies, blocks, expected in FIXTURES:
        width = sum(len(block) - 1 for block in blocks)
        rows = two_pin_rows(modulus, copies, blocks)
        actual = maximal_minor_gcd(rows, width)
        assert actual == expected, (modulus, copies, blocks, expected, actual)
        counts[actual] += 1
    assert counts == {1: 15, 3: 2, 27: 1}
    print(f"Two-pin integer stacks: {len(FIXTURES)} fixtures; maximal-minor GCD "
          f"1 in {counts[1]}, 3 in {counts[3]}, 27 in {counts[27]}.")
    print("CHECK_INTEGER_MINOR_GCD: ALL CHECKS PASSED.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without Python optimization so assertions execute.")
    main()
