"""Independent comparisons and recurrence checks for the implicit resultant note."""

from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import prod
from pathlib import Path
from random import Random

from flint import nmod_mat

from hash_probe import raw_residue, slot_gcd
from implicit import implicit_curve_jet, implicit_jet, three_blocks
from scalar import form, half_forms, polynomial, split_norm


def load_reference():
    path = Path(__file__).resolve().parents[1] / "scalar_norm" / "verify.py"
    spec = spec_from_file_location("explicit_scalar_reference", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reference = load_reference()


def moment_check(slots, target, p, a, c, chart):
    """Explicit reference halves ONLY in this audit of the local-series proposal."""
    cut = len(slots)//2
    A, _, f = half_forms(slots[:cut], None, p, a, c, False)
    B, _, h = half_forms(slots[cut:], target, p, a, c, False)
    z = polynomial([chart, -1], p)
    denominator, numerator = A.compose(z), B.compose(z)
    assert denominator[0] and numerator[0]
    length = f+h
    ratio = (numerator * denominator.inverse_series_trunc(length)).truncate(length)
    moments = [int(ratio[j]) if j >= 0 else 0 for j in range(1+h-f, f+h)]
    hankel = nmod_mat([moments[i:i+f] for i in range(f)], p)
    rank = hankel.rank()
    infinity_common = min(f - A.degree(), h - B.degree())
    common = A.gcd(B).degree() + infinity_common
    value = (-1)**(f*(f-1)//2) * pow(int(denominator[0]), f+h, p) * int(hankel.det()) % p
    assert rank == f-common
    assert value == three_blocks(slots, target, p, a, c, False)[0]
    return f, h, rank, common


def small():
    p, a, c = 13, 4, 1
    points = reference.old.points(p, a, c)
    targets = [x for x, y in points if 0 < y < p-y]
    seed = ((3, 4, 2), (7, 9, 2), (1, 8, 4))
    status = Counter()
    for window in product(*(reference.old.children(s) for s in seed)):
        slots = [sorted({x for x, y in slot}) for slot in reference.old.hash_slots(points, 16, window)]
        for target in targets:
            expected = reference.dense_jet(slots, target, p, a, c)
            assert implicit_jet(slots, target, p, a, c) == expected
            assert three_blocks(slots, target, p, a, c) == split_norm(slots, target, p, a, c)
            value, (dz, ds, *_) = expected
            status["nonzero" if value else "simple" if dz else "one repeated" if ds else "several"] += 1
    chart = next(x for x in range(p) if x not in targets)
    moments = 0
    for slots in ([[0, 2]]*m for m in (2, 3, 4, 5, 6)):
        for target in targets:
            assert implicit_jet(slots, target, p, a, c) == reference.raw_product(slots, target, p, a, c)
            moment_check(slots, target, p, a, c, chart)
            moments += 1
    print(f"All 63 dyadic boxes, {sum(status.values())} value and full-gradient checks: {dict(status)}.")
    print(f"{moments} Hankel determinant/rank checks including repeated slots and infinite roots.")


def torsion():
    p, a, c = 7, -1, 0
    points = reference.old.points(p, a, c)
    slots = reference.old.hash_slots(points, 8, ((1, 0, 4), (3, 1, 4), (5, 3, 4)))
    eligible = [pt for pt in points if pt[1]]
    for t in eligible:
        assert implicit_curve_jet(slots, t, p, a, c, reference.old.point_add) == reference.curve_jet(
            slots, t, p, a, c, reference.old.point_add)
    print(f"Full rational two-torsion: {len(eligible)} signed targets.")


def local_hits():
    p, a, M, width = 13, 4, 16, 2
    points = reference.old.points(p, a, 1)
    target = (0, 1)
    hits = Counter()
    for alpha in range(1, M, 2):
        for beta in range(M):
            slot, = reference.old.hash_slots(points, M, [(alpha, beta, width)])
            roots = set()
            for pt in slot:
                shifted = reference.old.point_add(target, pt, p, a)
                roots.add(shifted[0] if shifted is not None else None)
            hits.update(roots)
    seeds = M*M//2
    assert max(hits.values())*M <= seeds*2*width
    print(f"All {seeds} right-slot seeds: a single candidate joined abscissa hit at most "
          f"{max(hits.values())}/{seeds} (bound {2*width}/{M}).")


def hash_probes():
    checks = 0
    fixtures = [(13, 16, 4, 1, range(1, 16, 2), range(16), (1, 2, 4, 8)),
                (65521, 65536, 44, 1, (509,), (8371,), (16384,))]
    for p, M, a, c, multipliers, shifts, widths in fixtures:
        x = polynomial([0, 1], p)
        probes = (x, (x-2)**2*(x-3), x**4+3*x+1, (x-2)**13*(x-5)**3)
        for alpha, beta, width in product(multipliers, shifts, widths):
            inverse = pow(alpha, -1, M)
            preimages = [inverse*(j-beta) % M for j in range(width)]
            canonical = [u for u in preimages if u < p]
            accepted = [u for u in canonical if pow((u**3+a*u+c) % p, (p-1)//2, p) == 1]
            for probe in probes:
                raw = polynomial([1], p)
                for u in canonical:
                    raw = raw*(x-u) % probe
                genuine = polynomial([1], p)
                for u in accepted:
                    genuine = genuine*(x-u) % probe
                assert raw_residue(p, M, alpha, beta, width, probe) == raw
                assert slot_gcd(p, M, alpha, beta, width, a, c, probe) == genuine.gcd(probe)
                checks += 1
    print(f"{checks} raw residues and genuine-slot gcds checked, up to width 16384.")


def twist_check():
    p, a, c = 13, 4, 1
    points = reference.old.points(p, a, c)
    genuine = {x for x, _ in points}
    nonsquare = next(x for x in range(p) if x not in genuine)
    targets = sorted(genuine)
    checked = different = 0
    for collection in ([[0, nonsquare], [0, nonsquare], [0, 2]],
                       [[nonsquare], [nonsquare], [0]], [[0, 1, 2, 3]]*3):
        real = [[x for x in slot if x in genuine] for slot in collection]
        for t in targets:
            counts = []
            for slots in (real, collection):
                total = 0
                linear = polynomial([-t, 1], p)
                for row in product(*slots):
                    f = form(row, p, a, c)
                    while not f(t):
                        total += 1
                        f = f // linear
                counts.append(total)
            assert counts[0] <= counts[1] and (counts[1]-counts[0]) % 2 == 0
            different += counts[0] != counts[1]
            checked += 1
    assert different
    print(f"Raw/curve multiplicities: {checked} parity agreements, {different} differing integers.")


def large():
    p, a, c, M = 65521, 44, 1, 65536
    rng = Random(18112026)
    squares = {y*y % p for y in range(1, p)}
    for _ in range(2):
        slots = []
        for width in (8, 16, 16, 16):
            alpha, beta = rng.randrange(1, M, 2), rng.randrange(M)
            roots = [pow(alpha, -1, M)*(i-beta) % M for i in range(width)]
            slots.append([x for x in roots if x < p and (x**3+a*x+c) % p in squares])
    assert tuple(map(len, slots)) == (4, 7, 10, 10)
    expected = (0, [43446, 16953, 13292, 47523, 12119, 9163])
    assert implicit_jet(slots, 0, p, a, c) == expected == reference.raw_product(slots, 0, p, a, c)
    chart = next(x for x in range(p) if (x**3+a*x+c) % p not in squares)
    f, h, rank, common = moment_check(slots, 0, p, a, c, chart)
    sizes = [*map(len, slots), 1]
    groups = [sizes[j::3] for j in range(3)]
    degrees = [prod(g)*2**(len(g)-1) for g in groups]
    section_cap = 2*sorted(degrees)[1]
    assert degrees == [80, 14, 10] and section_cap == 28 and (f, h, rank, common) == (56, 400, 55, 1)
    print(f"Order 65419: blocks {degrees}, maximum section degree {section_cap}; "
          f"old halves {f},{h}, local Hankel rank {rank}, common degree {common}.")
    print(f"Norm and first derivatives: {expected}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    small()
    torsion()
    local_hits()
    hash_probes()
    twist_check()
    large()
