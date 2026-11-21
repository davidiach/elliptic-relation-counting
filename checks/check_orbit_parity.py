"""Check the binary address derivative, actual mask spans and fixed-target ranks."""
from collections import Counter
from itertools import product

from check_strided_isolation import curve


def rank(rows):
    basis = {}
    for word in rows:
        while word:
            pivot = word.bit_length()
            if pivot not in basis:
                basis[pivot] = word
                break
            word ^= basis[pivot]
    return len(basis)


def members(word):
    while word:
        low = word & -word
        yield low.bit_length()-1
        word ^= low


def check(p, a, c):
    points, sums = curve(p, a, c)
    n, modulus, width = len(points), 1 << p.bit_length(), 4
    symbols = [0]*modulus
    for i, (x, _) in enumerate(points[1:], 1):
        symbols[x] |= 1 << i
    inverse = [row.index(0) for row in sums]
    fixed = {i for i in range(n) if inverse[i] == i}
    pairs = [(1 << i) | (1 << j) for i, j in enumerate(inverse) if i < j]

    def convolve(left, right):
        word = 0
        for i, j in product(members(left), members(right)):
            word ^= 1 << sums[i][j]
        return word

    anchors = []
    for beta in range(modulus):
        full, high = 0, 0
        for x, symbol in enumerate(symbols):
            address = (x+beta) % modulus
            if address < width:
                full |= symbol
                if address >= width//2:
                    high |= symbol
        anchors.append((full, high))
    assert all(any(not symbols[x] for x in range(i, modulus, width//2))
               for i in range(width//2))
    assert rank(high for _, high in anchors) == sum(bool(word) for word in symbols)

    def boxes(free):
        return [sum(symbol for x, symbol in enumerate(symbols) if x & ~free == pattern)
                for pattern in range(modulus) if not pattern & free]

    products = set()
    singleton_checks = tangent_checks = 0
    for (full, high), second, third in product(anchors, boxes(5), boxes(3)):
        tail = convolve(second, third)
        high_counts = convolve(high, tail)
        products.add(high_counts)
        second_sign = sum(1 << i for i in members(second) if points[i][1] > (-points[i][1]) % p)
        joint_tangent = high_counts ^ convolve(convolve(full, second_sign), third)
        multiplicities, weighted = Counter(), 0
        for i, j, k in product(members(full), members(second), members(third)):
            target = sums[sums[i][j]][k]
            multiplicities[target] += 1
            address_bit = ((high >> i) ^ (second_sign >> j)) & 1
            weighted ^= address_bit << target
        assert weighted == joint_tangent
        singleton_checks += sum(count == 1 for count in multiplicities.values())
        tangent_checks += 1
    product_dimension = rank(products)
    assert rank(list(products)+pairs) == product_dimension
    assert len(pairs) == (n-len(fixed))//2

    ranks = {}
    split_ranks = []
    good_boxes = [box for box in boxes(5) if box and not fixed.intersection(members(box))]
    tails = sorted({convolve(second, third) for second, third in product(boxes(5), boxes(3))})
    actual = [[convolve(high, tail) for tail in tails] for _, high in anchors]
    for target in range(1, n):
        twice = sums[target][target]
        if not twice:
            continue
        subgroup = {0}
        point = twice
        while point:
            subgroup.add(point)
            point = sums[point][twice]
        images = [frozenset(sums[i][j] for j in subgroup) for i in range(n)]
        cosets = set(images)
        fixed_cosets = {images[i] for i in range(n) if sums[i][i] in subgroup}
        fixed_images = {images[i] for i in fixed}
        expected = len(pairs)-(len(cosets)+len(fixed_cosets))//2+len(fixed_images)
        gram = [sum(((convolve(x, y) >> target) & 1) << j for j, y in enumerate(pairs))
                for x in pairs]
        assert rank(gram) == expected
        assert expected >= n/4-2
        ranks.setdefault(len(subgroup), set()).add(expected)
        restricted = rank(sum(((convolve(pair, box) >> target) & 1) << j
                              for j, box in enumerate(good_boxes)) for pair in pairs)
        stable = sum(all(box & (1 << sums[i][twice]) for i in members(box)) for box in good_boxes)
        assert stable <= 4*width**2
        assert 2*restricted >= len(good_boxes)-stable
        split_rank = rank(sum(((word >> target) & 1) << j for j, word in enumerate(row))
                          for row in actual)
        assert split_rank >= restricted
        split_ranks.append(split_rank)
    print(f"p={p}, N={n}, rational two-torsion size={len(fixed)}: "
          f"{tangent_checks} mixed seeds, {singleton_checks} singleton target fibers, "
          f"product span {product_dimension} contains pair dimension {len(pairs)}.")
    print("  Fixed-target Gram ranks by order(2T): "
          + ", ".join(f"{order}->{sorted(values)}" for order, values in sorted(ranks.items())))
    print(f"  Exact anchor/tail seed matrix ranks across eligible targets: "
          f"{min(split_ranks)} to {max(split_ranks)}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    for parameters in ((13, 4, 1), (17, -1, 0), (29, 1, 1)):
        check(*parameters)
