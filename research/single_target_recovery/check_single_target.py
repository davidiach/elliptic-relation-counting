"""Verify the single-target recovery claims on small exact instances.

Every rank is computed exactly modulo the prime group order. The pinned
matrices come from the slot-by-slot entry formula of RESEARCH.md, which is
checked against direct enumeration on small curves and against the
independent construction in research/sparse_target_runs.
"""

if not __debug__:
    raise RuntimeError("This verifier requires Python assertions.")

import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np
from flint import nmod_mat

from research.deterministic_endgame_barriers.check_pinning_barriers import curve_cycle
from research.sparse_target_runs.check_four_targets import affine_block, pinned_matrices

HERE = Path(__file__).resolve().parent


def rank(rows, r):
    return nmod_mat(np.asarray(rows, dtype=np.int64).tolist(), r).rank() if len(rows) else 0


def counts(blocks, r):
    """Integer number of ways to write each residue as one element from each listed block."""
    out = np.zeros(r, dtype=np.int64)
    out[0] = 1
    for block in blocks:
        out = sum(np.roll(out, int(x)) for x in block)
    return out


class Pinned:
    """Pinned rows R_{t,(i,A)} restricted to block columns, from the entry formula."""

    def __init__(self, blocks, r):
        self.blocks = [np.array(sorted(int(x) % r for x in b), dtype=np.int64) for b in blocks]
        self.r, self.m = r, len(blocks)
        self.sizes = [len(b) for b in self.blocks]
        self.offsets = np.cumsum([0, *self.sizes])
        self.n0 = int(self.offsets[-1])
        self.labels = np.concatenate(self.blocks)
        slots = [i for i in range(self.m) for _ in range(2)]
        self.E, self.phi = [], {}
        for i in range(self.m):
            rest = slots.copy()
            rest.remove(i)
            self.E.append(counts([self.blocks[k] for k in rest], r))
            for j in range(self.m):
                others = rest.copy()
                others.remove(j)
                self.phi[i, j] = counts([self.blocks[k] for k in others], r)

    def integer(self, t):
        """Return M_t over the integers and the pin counts n_h(t)."""
        r, M = self.r, np.zeros((self.n0, self.n0), dtype=np.int64)
        pins = np.zeros(self.n0, dtype=np.int64)
        for i in range(self.m):
            for a, x in enumerate(self.blocks[i]):
                h = self.offsets[i] + a
                pins[h] = self.E[i][(t - x) % r]
                M[h, h] += pins[h]
                for j in range(self.m):
                    c = 1 if j == i else 2
                    M[h, self.offsets[j]:self.offsets[j+1]] += c*self.phi[i, j][(t - x - self.blocks[j]) % r]
        return M, pins

    def nullity(self, t):
        return self.n0 - rank(self.integer(t)[0] % self.r, self.r)

    def gauge(self):
        """Block-constant vectors with zero sum: the space B_0."""
        B = np.zeros((self.m - 1, self.n0), dtype=np.int64)
        for k in range(self.m - 1):
            B[k, self.offsets[0]:self.offsets[1]] = 1
            B[k, self.offsets[k+1]:self.offsets[k+2]] = -1
        return B % self.r

    def decompose(self, g):
        """Lexicographically first representation of g with two slots per block."""
        r, slots = self.r, [i for i in range(self.m) for _ in range(2)]
        reach = [None]*(2*self.m + 1)
        reach[-1] = np.zeros(r, dtype=bool)
        reach[-1][0] = True
        for k in range(2*self.m - 1, -1, -1):
            reach[k] = np.zeros(r, dtype=bool)
            for x in self.blocks[slots[k]]:
                reach[k] |= np.roll(reach[k+1], int(x))
        columns, rest = [], g % r
        for k, i in enumerate(slots):
            a = next(a for a, x in enumerate(self.blocks[i]) if reach[k+1][(rest - x) % r])
            columns.append(self.offsets[i] + a)
            rest = (rest - self.blocks[i][a]) % r
        assert rest == 0
        return columns

    def explicit(self, x):
        """Normalisation and consistency rows for columns blocks, *P (log 1), *Q (log x)."""
        n = self.n0 + 2
        rows = [np.eye(1, n, self.n0, dtype=np.int64)[0]]
        for j, g in enumerate([*self.labels, 1, x]):
            row = np.zeros(n, dtype=np.int64)
            row[j] += 1
            for c in self.decompose(int(g)):
                row[c] -= 1
            rows.append(row % self.r)
        rhs = np.zeros(len(rows), dtype=np.int64)
        rhs[0] = 1
        return np.array(rows), rhs

    def target_rows(self, t):
        M, pins = self.integer(t)
        rows = np.hstack([M % self.r, np.zeros((self.n0, 2), dtype=np.int64)])
        return rows, (t*pins) % self.r


def solve(rows, rhs, r):
    """Return the unique solution of a consistent full-column-rank system, else None."""
    n = rows.shape[1]
    reduced, rk = nmod_mat(np.hstack([rows, rhs[:, None]]).tolist(), r).rref()
    if rk != n or any(int(reduced[i, i]) != 1 for i in range(n)):
        return None
    return np.array([int(reduced[i, n]) for i in range(n)], dtype=np.int64)


def enumerate_pinned(blocks, r, t):
    """Count pinned rows directly from all ordered tuples with two slots per block."""
    P = Pinned(blocks, r)
    M = np.zeros((P.n0, P.n0), dtype=np.int64)
    choices = [range(P.offsets[i], P.offsets[i+1]) for i in range(P.m) for _ in range(2)]
    for cols in product(*choices):
        if P.labels[list(cols)].sum() % r != t:
            continue
        for pin in cols[0::2]:
            for c in cols:
                M[pin, c] += 1
    return M


def check_small_curves():
    for p, a, b, picks in ((37, 2, 8, ((1, 4, 9), (2, 6, 11))), (41, 1, 5, ((1, 3), (5, 12, 20), (7, 9, 30)))):
        cycle = curve_cycle(p, a, b)
        r = len(cycle)
        assert all(r % k for k in range(2, int(r**0.5) + 1))
        P = Pinned(picks, r)
        for t in range(r):
            M, pins = P.integer(t)
            assert np.array_equal(M, enumerate_pinned(picks, r, t))
            assert np.array_equal(M @ np.ones(P.n0, dtype=np.int64), 2*P.m*pins)
            u = (P.labels - t*pow(2*P.m, -1, r)) % r
            assert not np.any(M @ u % r) and not np.any(M @ P.gauge().T % r)
        stack = np.vstack([P.integer(t)[0] % r for t in range(r)])
        assert P.n0 - rank(stack, r) == P.m - 1
        kernel = nmod_mat(stack.tolist(), r).nullspace()[0]
        basis = np.array([[int(kernel[i, j]) for j in range(P.m - 1)] for i in range(P.n0)]).T
        assert rank(np.vstack([basis, P.gauge()]), r) == P.m - 1
        print(f"p={p}, r={r}, m={P.m}, sizes={P.sizes}: entry formula equals direct enumeration at all "
              f"{r} targets; forced kernel checked; all-target stack nullity {P.m-1} = dim B_0.")


def recovery_report(P, x, targets, label):
    E, rhs0 = P.explicit(x)
    truth = np.array([*P.labels, 1, x], dtype=np.int64) % P.r
    assert not np.any((E @ truth - rhs0) % P.r)
    explicit_rank = rank(E, P.r)
    full = []
    for t in targets:
        R, rhs = P.target_rows(t)
        assert not np.any((R @ truth - rhs) % P.r)
        rows = np.vstack([E, R])
        if rank(rows, P.r) == P.n0 + 2:
            full.append(t)
    t0 = targets[0]
    R, rhs = P.target_rows(t0)
    solution = solve(np.vstack([E, R]), np.concatenate([rhs0, rhs]), P.r)
    assert solution is not None and np.array_equal(solution, truth)
    print(f"  {label}: explicit rank {explicit_rank}/{P.n0+2}; explicit rows plus M_t have rank "
          f"{P.n0+2} for {len(full)}/{len(targets)} targets; solving at t={t0} returns log Q={solution[-1]}.")
    return full


def check_curve_1009():
    case = json.loads((HERE / "fixture_1009.json").read_text())
    p, a, b, r, m = case["p"], case["a"], case["b"], case["r"], len(case["seed"])
    cycle = curve_cycle(p, a, b)
    assert len(cycle) == r and all(r % k for k in range(2, int(r**0.5) + 1))
    q = next(k for k in range(1, r) if k**m >= 16*r)
    threshold = (p*q + r - 2)//(r - 1)
    assert threshold == case["threshold"] and q < 4096*m
    blocks = [affine_block(cycle, p, threshold, maps) for maps in case["seed"]]
    assert [len(x) for x in blocks] == case["sizes"]
    half = counts(blocks, r)
    mu = Fraction((r - 1)*threshold, p)**m
    phi = (r*int(half @ half) + 64*(int(half.sum()) - mu)**2)/mu**2
    assert phi == Fraction(case["phi"]) and phi < Fraction(64, 57) and np.all(half > 0)
    P = Pinned(blocks, r)
    assert np.array_equal(P.integer(0)[0] % r, pinned_matrices(blocks, r)(0) % r)
    hist, bad = {}, []
    for t in range(r):
        k = P.nullity(t)
        hist[k] = hist.get(k, 0) + 1
        if k != m:
            bad.append(t)
    print(f"p={p}, r={r}, m={m}, q={q}, T={threshold}, sizes={P.sizes}, Phi={float(phi):.6f} < 64/57, "
          f"half-sum support {int(np.count_nonzero(half))}/{r}")
    print(f"  all-target scan: nullity histogram {dict(sorted(hist.items()))}; H_1(t) fails at "
          f"{len(bad)} targets {bad}.")
    full = recovery_report(P, r - 3, list(range(r)), "single-target recovery, Q=[r-3]P")
    assert len(full) == r


def check_sparse_witnesses(witness_file=HERE.parent / "sparse_target_runs" / "witnesses.json"):
    cases = json.loads(witness_file.read_text(encoding="utf-8"))
    for case in cases:
        p, a, b, r = (case[k] for k in ("p", "a", "b", "r"))
        cycle = curve_cycle(p, a, b)
        blocks = [affine_block(cycle, p, case["threshold"], maps) for maps in case["seed"]]
        P = Pinned(blocks, r)
        reference = pinned_matrices(blocks, r)
        for t in (0, 1):
            assert np.array_equal(P.integer(t)[0] % r, reference(t) % r)
        exceptional = [t for t in range(r) if P.nullity(t) > P.m]
        recorded = case["exceptional_logs"]
        assert all(type(t) is int for t in recorded) and sorted(recorded) == exceptional, (
            f"r={r}, m={P.m}: exceptional targets differ: recorded {recorded}, computed {exceptional}"
        )
        status = "fails" if 0 in exceptional else "holds"
        print(f"r={r}, m={P.m}, sizes={P.sizes}: H_1(0) {status} (nullity {P.nullity(0)}); "
              f"exceptional targets X={len(exceptional)} {exceptional}, N_0={P.n0}")
        for x in (2, r - 3):
            recovery_report(P, x, [0], f"Q=[{x}]P")


def digit_blocks(r, m, n, delta):
    D = 2*n - 1
    for s1 in range(1, r):
        s = [((1 - delta)*pow(2, -1, r) - s1*(m - 1)) % r] + [s1]*(m - 1)
        blocks = [[(s[i] + D**i*k) % r for k in range(n)] for i in range(m)]
        if all(0 not in block for block in blocks):
            return blocks


def check_digit():
    r, m, n, delta = 1087, 2, 17, 2
    D = 2*n - 1
    assert D**m - delta == r
    blocks = digit_blocks(r, m, n, delta)
    P = Pinned(blocks, r)
    top = P.offsets[m - 1]
    first = {}
    for t in range(1, r):
        M = P.integer(t)[0]
        for j in range(1, n):
            if j not in first and np.any(M[:, top + j]):
                first[j] = t
    assert all(first[j] == D**(m - 1)*j + 1 - delta for j in range(1, n))
    X = sum(P.nullity(t) > m for t in range(r))
    single = counts(blocks, r)
    covered = int(np.count_nonzero(single))
    assert X == r and covered <= n**m < Fraction(7, 8)*r
    print(f"digit family r={r}, m={m}, n={n}: column (m-1,j) first used at t=D^(m-1) j-1 for every j; "
          f"H_1(t) fails at all {X} targets; single-slot sumset covers {covered}/{r} < 7r/8.")


def check_subgroup():
    r, k = 541, 90
    g = next(g for g in range(2, r) if all(pow(g, (r - 1)//q, r) != 1 for q in (2, 3, 5)))
    assert r - 1 == 2**2 * 3**3 * 5
    h = pow(g, (r - 1)//k, r)
    H = sorted(pow(h, j, r) for j in range(k))
    P = Pinned([H, H], r)
    bad = [t for t in range(r) if P.nullity(t) > 2]
    single = counts([H, H], r)
    W, C = int(single.sum()), int(single @ single)
    assert Fraction(r*C, W*W) < Fraction(64, 57) and np.count_nonzero(single) > Fraction(7, 8)*r
    print(f"subgroup family r={r}, F_1=F_2=H with |H|={k}: N_0={P.n0}, coverage "
          f"{int(np.count_nonzero(single))}/{r}, rC/W^2={r*C/W**2:.6f}; H_1(t) fails at X={len(bad)} targets.")


def check_chevalley_warning():
    r, n, h = 5, 2, 2
    M = n + 2*h
    rng = np.random.default_rng(7)
    weights = rng.integers(0, r, size=(h, M))
    support = M - n + 1
    tau = next(v for v in product(range(r), repeat=support) if any(v)
               and all(sum(w[e]*v[e]**2 for e in range(support)) % r == 0 for w in weights))
    A = np.zeros((M, n), dtype=np.int64)
    A[:support, 0] = tau
    for c in range(1, n):
        A[support + c - 1, c] = 1
    assert rank(A, r) == n
    for w in weights:
        B = (A.T * w) @ A % r
        assert rank(B, r) < n
    print(f"Chevalley-Warning construction r={r}, N={n}, h={h}, M={M}: rank-N array singular "
          f"under both fixed weightings.")


if __name__ == "__main__":
    check_small_curves()
    check_curve_1009()
    check_sparse_witnesses()
    check_digit()
    check_subgroup()
    check_chevalley_warning()
