"""Check simultaneous shifted abscissa tests against exact signed incidence."""
from itertools import product

from flint import nmod_mat

from check_block_support import group_tables
from fixtures import curve_cycle


def check_curve():
    points, add, neg = group_tables(13, 4, 1)
    n = len(points)
    assert n == 19
    shifts = [0]
    for _ in range(3):
        shifts.append(add[shifts[-1]][1])
    assert len(set(shifts)) == 4
    orbit = [min(i, neg[i]) for i in range(n)]
    false_counts = [0] * len(shifts)
    for a, b, t in product(range(n), repeat=3):
        true = add[a][b] == t
        survives = True
        for index, offset in enumerate(shifts):
            u, v = add[a][offset], add[b][offset]
            goal = orbit[add[t][add[offset][offset]]]
            survives &= goal in (orbit[add[u][v]], orbit[add[u][neg[v]]])
            assert survives or not true
            false_counts[index] += int(survives and not true)
        assert survives == true
    print(f"F_13, order 19: all {n**3:,} ordered point/point/target triples. "
          f"Wrong triples after each successive offset: {false_counts}.")


def check_determinant_obstruction(p, label, certify_order_bound):
    m = 3
    dim = m+1
    height = 16  # dim^(dim/2)
    radix = 2*m*height + 1
    limit = 2*m*height*sum(radix**a for a in range(dim+1))
    if certify_order_bound:
        assert limit < p
    offsets = [radix**a for a in range(dim+1)]
    signs = [s for s in product((-1, 1), repeat=m) if -1 in s]
    checked = 0
    for patterns in product(signs, repeat=dim+1):
        matrix = [[*s, -1] for s in patterns]
        augmented = [[*row, 2*s.count(-1)*offset]
                     for row, s, offset in zip(matrix, patterns, offsets)]
        assert nmod_mat(matrix, p).rank() < nmod_mat(augmented, p).rank()
        checked += 1
    print(f"{label}: all {checked:,} nontrivial "
          f"sign-pattern systems inconsistent across {len(offsets)} offsets at arity {m}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_curve()
    actual_order = len(curve_cycle(1031, 9))
    assert actual_order == 1019
    check_determinant_obstruction(actual_order, "F_1031, cyclic curve order 1019", False)
    check_determinant_obstruction((1 << 61) - 1, "Cyclic additive group of order 2^61-1", True)
