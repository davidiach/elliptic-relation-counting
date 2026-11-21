"""Verify balanced lifted-product equality on actual sampled abscissas."""
from itertools import combinations, product
from math import prod

from check_section_deformation import lift, tangent_row
from check_strided_isolation import curve
from fixtures import determinant, evaluation


def norm(values, p):
    total = 1
    for value in values:
        total = total*value % p
    return total


def check_grid(p, points, sums, slots, target, thorough=True):
    size = len(slots)+1
    endpoint = (points[target][0], -points[target][1] % p)
    sizes = list(map(len, slots))
    volume = prod(sizes)
    group_slots = [[points[i] for i in slot] for slot in slots]
    grid = list(product(*group_slots))

    def chord(left, right):
        return (left[1]+right[1])*pow(left[0]-right[0], -1, p) % p

    def reciprocal(slot, x):
        return sum(pow(x-Q[0], -1, p) for Q in slot)*pow(len(slot), -1, p) % p

    a = [-reciprocal(slot, endpoint[0]) % p for slot in group_slots]
    coeffs = [(volume*endpoint[1]*(v+sum(a))) % p for v in a]
    dets, endpoint_partials = [], []
    hits = sum(sums[sums[i][j]][k] == target for i, j, k in product(*slots))
    directions = [[int(i == j) for i in range(size-1)] for j in range(size-1)]
    labels = [lambda P, i=i: (i+1)*(P[0]+2*P[1])+7 for i in range(size-1)]
    local_direct = 0
    for chosen in grid:
        row_points = [*chosen, endpoint]
        rows = [evaluation(P, size, p) for P in row_points]
        delta = determinant(rows, p)
        partials = [determinant(rows[:i]+[tangent_row(P, size, 1, p)]+rows[i+1:], p)
                    for i, P in enumerate(row_points)]
        kval = [sum(chord(P, Q) for j, Q in enumerate(row_points) if j != i)
                for i, P in enumerate(row_points)]
        for i, j in combinations(range(size), 2):
            assert (partials[i]-partials[j]-(kval[i]-kval[j])*delta) % p == 0
        local_direct += sum(label(P)*(K-kval[-1]) for label, P, K in zip(labels, chosen, kval))
        dets.append(delta)
        endpoint_partials.append(partials[-1])
    assert dets.count(0) == hits
    full = norm(dets, p)
    local = 0
    for i, slot in enumerate(group_slots):
        for P in slot:
            inner = 2*chord(P, endpoint)
            for j, other in enumerate(group_slots):
                if i != j:
                    inner += P[1]*reciprocal(other, P[0])-endpoint[1]*reciprocal(other, endpoint[0])
            local += (volume//len(slot))*labels[i](P)*inner
    assert (local-local_direct) % p == 0

    def lifted_product(speeds, e, offset):
        lifted_slots = [[evaluation(lift(P, speeds[i](P), 1, p, offset), size, p*p)
                         for P in slot] for i, slot in enumerate(group_slots)]
        last = evaluation(lift(endpoint, e, 1, p, offset), size, p*p)
        return norm((determinant([*rows, last], p*p) for rows in product(*lifted_slots)), p*p)

    for offset in ((0, 3) if thorough else (0,)):
        base = lifted_product([lambda _: 0]*len(slots), 0, offset)
        assert base % p == full
        for j, speeds in enumerate(directions):
            shifted = lifted_product([lambda _, c=c: c for c in speeds], -1, offset)
            assert (shifted-base) % (p*p) == p*((full*coeffs[j]) % p)
            if coeffs[j]:
                assert (base == shifted) == (hits > 0)
        if thorough:
            weighted = lifted_product(labels, 5, offset)
            common = 0
            for k, chosen in enumerate(grid):
                speed = 5+sum(f(P) for f, P in zip(labels, chosen))
                common += speed*endpoint_partials[k]*norm(dets[:k]+dets[k+1:], p)
            assert (weighted-base) % (p*p) == p*((full*local+common) % p)
    return coeffs, hits


def run():
    p = 29
    points, sums = curve(p, 1, 1)
    target = next(i for i, P in enumerate(points) if P and P[0] == 10)
    exceptional = [
        (((24, 25), (6, 14), (8, 12)), 0),
        (((22, 27), (8, 12), (24, 25)), 1),
        (((22, 27), (6, 14), (8, 12)), 2),
    ]
    for raw, expected in exceptional:
        modulus = 1 << p.bit_length()
        left, right = raw[0]
        alpha = pow(right-left, -1, modulus)
        beta = -alpha*left % modulus
        assert all(((alpha*x+beta) % modulus < 2) == (x in raw[0]) for x in range(modulus))
        for low, high in raw[1:]:
            assert (low ^ high).bit_count() == 1
        slots = [[i for i, P in enumerate(points) if P and P[0] in pair] for pair in raw]
        coefficients, hits = check_grid(p, points, sums, slots, target)
        assert not any(coefficients) and hits == expected
        free_bit = raw[1][0] ^ raw[1][1]
        answers = []
        for bit in (0, 1):
            part = [i for i in slots[1] if bool(points[i][0] & free_bit) == bit]
            cs, count = check_grid(p, points, sums, [slots[0], part, slots[2]], target, False)
            assert any(cs)
            answers.append(count)
        assert sum(answers) == hits
        print(f"p={p}, target x=10, {hits} witnesses: all three balanced coefficients vanish. "
              f"Pinning one fixed-slot bit gives calibrated subgrid witness counts {answers}.")
    print("Every tuple passed the alternant identity for all row pairs. "
          "Two independent lift choices passed balanced and address-dependent product checks.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run()
