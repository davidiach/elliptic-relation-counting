"""Count false addresses and internal completions of hashed chord trees."""
from collections import Counter
from fractions import Fraction
from itertools import product

from check_block_support import group_tables


def binary_product(a, b, modulus):
    """Multiplication in the fixture's binary extension field."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a.bit_length() == modulus.bit_length():
            a ^= modulus
    return result


def hashes(n, bits, modulus):
    mask = (1 << bits) - 1
    size = 1 << (modulus.bit_length() - 1)
    return [tuple(binary_product(scale, x, modulus) & mask for x in range(n))
            for scale in range(size)]


def check():
    points, add, neg = group_tables(13, 4, 1)
    n = len(points)
    slot_x = (0, 2, 3, 4)
    slots = [tuple(i for i, p in enumerate(points) if p and p[0] == x) for x in slot_x]
    addresses = list(product(*slots))
    sums = [add[add[a][b]][add[c][d]] for a, b, c, d in addresses]
    counts = Counter(sums)
    target = next(t for t in sorted(counts) if counts[t] == 1 and t and points[t][0] not in slot_x)
    assert len(addresses) == 16 and n == 19
    for bits1, bits2, modulus in ((2, 3, 0b100101), (5, 5, 0b100101),
                                  (7, 7, 0b10000011)):
        left_hashes = hashes(n, bits1, modulus)
        right_hashes = hashes(n, bits2, modulus)
        q1, q2 = 1 << bits1, 1 << bits2
        for table, q in ((left_hashes, q1), (right_hashes, q2)):
            assert all(sum(h[x] == h[y] for h in table) * q == len(table)
                       for x in range(n) for y in range(x))
        alpha = (1 - Fraction(1, q1)) * (1 - Fraction(1, q2))
        mass = (1 + Fraction(n-1, q1)) * (1 + Fraction(n-1, q2))
        wrong_paths = (mass - alpha) / n
        for shape in ("balanced", "chain"):
            total_wrong = 0
            for (a, b, c, d), correct_sum in zip(addresses, sums):
                first_sum = add[a][b]
                path_total = 0
                address_total = 0
                single_jump_total = 0
                for h1 in left_hashes:
                    for h2 in right_hashes:
                        paths = 0
                        one_jump = False
                        for candidate in range(n):
                            if h1[first_sum] != h1[candidate]:
                                continue
                            if shape == "balanced":
                                second_sum = add[c][d]
                                consumed = add[target][neg[candidate]]
                            else:
                                second_sum = add[candidate][c]
                                consumed = add[target][neg[d]]
                            if h2[second_sum] == h2[consumed]:
                                paths += 1
                                one_jump |= (first_sum == candidate or second_sum == consumed)
                        if correct_sum == target:
                            assert paths >= 1
                        path_total += paths
                        address_total += paths > 0
                        single_jump_total += one_jump
                trials = len(left_hashes) * len(right_hashes)
                expected = wrong_paths + (alpha if correct_sum == target else 0)
                assert path_total == expected * trials
                if correct_sum != target:
                    assert single_jump_total == (1-alpha)*trials
                    assert single_jump_total <= address_total <= path_total
                    total_wrong += address_total
            if bits1 == 7:
                assert 3*total_wrong < trials
            print(f"F_13, {shape}, tags {bits1}+{bits2} bits: 1 witness kept for all {trials} "
                  f"hash pairs; mean wrong addresses {Fraction(total_wrong,trials)}, "
                  f"one-jump lower bound {15*(1-alpha)}, mean wrong completions {15*wrong_paths}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check()
