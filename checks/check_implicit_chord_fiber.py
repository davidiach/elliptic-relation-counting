"""Check target-fixed partial fractions, sign pooling and the occupied kernel rank."""
from collections import Counter
from itertools import product
from random import Random

from check_strided_isolation import curve
from check_direct_fiber import energy
from fixtures import curve_cycle


def matrix_rank(rows, p):
    basis = {}
    for row in rows:
        work = list(row)
        for i in range(len(work)):
            value = work[i]
            if value:
                if i not in basis:
                    scale = pow(value, -1, p)
                    basis[i] = [v*scale % p for v in work]
                    break
                work = [(v-value*b) % p for v, b in zip(work, basis[i])]
    return len(basis)


def check(p, linear, seed):
    pts, add = curve(p, linear, 1)
    neg = [row.index(0) for row in add]
    reps = [i for i in range(len(pts)) if i <= neg[i]]
    ordinary = [i for i in reps if i != neg[i]]
    fixed = set(reps)-set(ordinary)
    target = next(i for i in ordinary if pts[i][0] != 0)
    t, y = pts[target]

    def poles_and_constant(B, T):
        tx, ty = pts[T]
        residues = {}
        for r in ordinary:
            v = pts[r][1]*(B[add[r][T]]-B[add[r][neg[T]]])*pow(ty, -1, p) % p
            if v:
                residues[pts[r][0]] = v
        const = sum(m*pow(tx-pts[q][0], -1, p) for q, m in B.items()
                    if q not in (0, T, neg[T]))
        const += B[T]*(3*tx*tx+linear)*pow(2*ty*ty, -1, p)
        return residues, const % p

    def quadratic(s, u):
        if u is None:
            return 2, (t-s) % p
        mid = (s+u)*(s*u+linear)+2
        num = 2*(s-u)**2*t-2*mid
        den = (s-u)**2*t*t-2*mid*t+(s*u-linear)**2-4*(s+u)
        return num % p, den % p

    def jet(z, plus, minus):
        return [z, (plus-minus)*pow(2, -1, p) % p, (plus+minus-2*z)*pow(2, -1, p) % p]

    columns = 0
    for q in reps:
        B = Counter([q, neg[q]])
        poles, constant = poles_and_constant(B, target)
        u = pts[q][0] if q else None
        for s in range(p):
            raw = [quadratic(s+d, u) for d in (0, 1, -1)]
            num, den = (jet(*col) for col in zip(*raw))
            i = next((i for i, v in enumerate(num) if v), 3)
            j = next(j for j, v in enumerate(den) if v)
            if i < j:
                assert j-i == 1 and num[i]*pow(den[j], -1, p) % p == poles[s]
            else:
                assert s not in poles
                val = num[j]*pow(den[j], -1, p) % p
                assert val == (constant+sum(v*pow(s-z, -1, p) for z, v in poles.items())) % p
            columns += 1

    rng = Random(seed)
    modulus = 1 << p.bit_length()
    empty = singleton = bit_one = branch_only = rank_cases = tries = 0
    best_rank = largest_boundary = 0
    while tries < 800 and (rank_cases < 15 or empty < 50):
        tries += 1
        alpha, beta = rng.randrange(1, modulus, 2), rng.randrange(modulus)
        anchor = [(i, (alpha*P[0]+beta) % modulus) for i, P in enumerate(pts) if P]
        anchor = [(i, u) for i, u in anchor if u < 2]
        slots = [[i for i, _ in anchor]]
        for _ in range(3):
            bit = 1 << rng.randrange(p.bit_length())
            pattern = rng.randrange(modulus) & ~bit
            slots.append([i for i, P in enumerate(pts) if P and P[0] & ~bit == pattern])
        xx = [{pts[i][0] for i in slot} for slot in slots]
        if any(not slot for slot in slots) or sum(map(len, xx)) != len(set.union(*xx)):
            continue
        selected = [i for i, u in anchor if u == 1]
        if not selected:
            continue
        A = Counter(add[i][j] for i, j in product(selected, slots[1]))
        full = Counter(add[i][j] for i, j in product(slots[0], slots[1]))
        B = Counter(add[i][j] for i, j in product(*slots[2:]))
        for T in ordinary:
            tx, _ = pts[T]
            if any(tx in slot for slot in xx):
                continue
            hits = sum(v*B[add[T][neg[i]]] for i, v in full.items())
            if hits > 1:
                continue
            bit = sum(v*B[add[T][neg[i]]] for i, v in A.items())
            residues, const = poles_and_constant(B, T)
            largest_boundary = max(largest_boundary, len(residues))
            colliding = {pts[i][0] for i in ordinary if A[i] and
                         (B[add[i][T]] or B[add[neg[i]][T]])}
            branches = sum(A[i]*B[add[T][neg[i]]] for i in fixed)
            assert bool(bit) == (any(pts[i][0] in residues for i in ordinary if A[i]) or bool(branches))
            if hits:
                singleton += 1
                bit_one += bit
                branch_only += bool(bit and branches)
            else:
                empty += 1
                direct = sum(a*b*pow(tx-pts[add[i][j]][0], -1, p)
                             for i, a in A.items() for j, b in B.items() if add[i][j])
                pooled = sum(a*(const+sum(v*pow(pts[i][0]-z, -1, p) for z, v in residues.items()))
                             if i else a*const for i, a in A.items())
                assert (pooled-direct) % p == 0
            row_x = [pts[i][0] for i in ordinary if A[i] and pts[i][0] not in colliding]
            candidates = [q for q in ordinary if B[q] and
                          all(r not in fixed for r in (add[q][T], add[q][neg[T]]))]
            used = set()
            columns_away = []
            for q in candidates:
                pole, c0 = poles_and_constant(Counter([q, neg[q]]), T)
                assert len(pole) == 2
                if used.isdisjoint(pole):
                    used.update(pole)
                    columns_away.append([ (c0+sum(v*pow(s-z, -1, p) for z, v in pole.items())) % p
                                          for s in row_x])
            assert 3*len(columns_away) >= len(candidates)
            k = min(len(columns_away), max(0, (len(row_x)-1)//2))
            assert matrix_rank(columns_away[:k], p) == k
            rank_cases += int(k >= 2)
            best_rank = max(best_rank, matrix_rank(columns_away, p))
    assert empty and bit_one
    print(f"p={p}, N={len(pts)}: {columns} full-field kernel column checks, "
          f"{empty} empty regular targets, {singleton} singleton targets ({bit_one} bit one).")
    print(f"  {branch_only} selected witnesses required the branch-point check. "
          f"Largest reduced right denominator degree {largest_boundary}, "
          f"largest checked independent-column matrix rank {best_rank}.")


def check_fixed_half():
    points, add = curve(29, 5, 1)
    neg = [row.index(0) for row in add]
    slots = [[i for i, pt in enumerate(points) if pt and pt[0] == x] for x in (15, 21, 0, 1)]
    a, b = slots[0][0], slots[1][0]
    left = add[a][b]
    target = points.index((11, 16))
    right = Counter(add[i][j] for i,j in product(*slots[2:]))
    assert neg[left] == left and right[add[target][left]] == 1
    assert right[add[left][target]] == right[add[left][neg[target]]]
    print("p=29, full rational two-torsion: a regular singleton on the four "
          "pinned abscissas 15,21 | 0,1 needs the separate fixed-left-sum check.")


def check_slot_boundaries():
    p = 1031
    pts = curve_cycle(p, 9)
    n = len(pts)
    modulus = 1 << p.bit_length()
    rng = Random(34987)
    good = largest_size = largest_overlap = 0
    for _ in range(120):
        alpha, beta = rng.randrange(1, modulus, 2), rng.randrange(modulus)
        slot = [i for i, P in enumerate(pts) if P and (alpha*P[0]+beta) % modulus < 8]
        if not slot or energy(pts, slot)[1]:
            continue
        good += 1
        counts = Counter((i-j) % n for i,j in product(slot, repeat=2) if i != j)
        overlap = max(counts.values(), default=0)
        assert overlap <= 4
        largest_overlap = max(largest_overlap, overlap)
        largest_size = max(largest_size, len(slot))
    assert good and largest_size > 4
    print(f"p={p}: {good} nonempty affine masks with paired four-term relations, "
          f"up to {largest_size} points; maximum overlap across all nonzero shifts {largest_overlap}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check(29, 1, 982)
    check(101, 17, 3081)
    check_fixed_half()
    check_slot_boundaries()
