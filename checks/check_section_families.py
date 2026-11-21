"""Batch section tests on a pencil, a conic and the complete tiny anchor law."""
from collections import Counter, defaultdict
from itertools import product
from math import prod

from check_nonlinear_branches import mixed_fixture
from check_projective_coefficient import directions
from probe_signed_moment import Probe


def normal(row, p):
    row = tuple(x % p for x in row)
    nonzero = next((x for x in row if x), 0)
    if not nonzero:
        return None
    scale = pow(nonzero, -1, p)
    return tuple(x*scale % p for x in row)


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def rows(probe, slot, target):
    endpoint = probe.vectors[-target][1:]
    return [(tuple((x-y) % probe.p for x, y in zip(probe.vectors[i][1:], endpoint)), addr)
            for i, addr in slot]


def check_bulk_sums():
    probe = Probe(13, 4, (0,))
    p = probe.p
    slots = [[(i, 2*j+a) for j, x in enumerate(xs) for i, a in probe.by_x.get(x, ())]
             for xs in ((0, 1), (0, 4), (0, 8))]
    slot_rows = [rows(probe, s, 1) for s in slots]
    tests = 0
    seen_base = False
    for axis in directions(p, 3):
        pencil = [v for v in directions(p, 3) if dot(v, axis) % p == 0]
        base_rows, exceptional = [], []
        for slot in slot_rows:
            base, by_param = [], defaultdict(list)
            for row, addr in slot:
                v = normal(cross(axis, row), p)
                if v is None:
                    base.append(addr)
                else:
                    by_param[v].append(addr)
            base_rows.append(base)
            exceptional.append(by_param)
        seen_base |= any(base_rows)
        exceptional_set = set().union(*exceptional)
        for r in range(probe.size):
            masks = [r//probe.block**k % probe.block for k in (2, 1, 0)]
            def mass(addrs, mask):
                return sum((-1)**((a & mask).bit_count()) for a in addrs)
            base = [mass(a, k) for a, k in zip(base_rows, masks)]
            bulk = (p+1-len(exceptional_set))*prod(base)
            bulk += sum(prod(u+mass(exc.get(v, ()), k) for u, exc, k in zip(base, exceptional, masks))
                        for v in exceptional_set)
            literal = sum(prod(mass([a for row, a in slot if dot(v, row) % p == 0], k)
                               for slot, k in zip(slot_rows, masks)) for v in pencil)
            assert bulk == literal
            tests += 1
    assert seen_base
    print(f"F_{p}: {tests} exact signed pencil sums including common base points and overlapping masks.")


def check_regular():
    probe = Probe(101, 17, (0,))
    p, target = probe.p, 12
    tail = mixed_fixture(probe)[:2]
    assert len({probe.curve[i][0] for i, _ in tail[0]} & {probe.curve[i][0] for i, _ in tail[1]}) == 0
    visible = [rows(probe, s, target) for s in tail]
    partial = Counter(normal(cross(v, u), p) for (v, _), (u, _) in product(*visible))
    assert len(partial) == 15 and max(partial.values()) == 2 and None not in partial
    witness_mass = Counter()
    seen_seeds = []
    excluded = {probe.curve[target][0]} | {probe.curve[i][0] for slot in tail for i, _ in slot}
    for alpha in range(1, probe.modulus, 2):
        inverse = pow(alpha, -1, probe.modulus)
        for beta in range(probe.modulus):
            slot = [(i, 2*u+sign) for u in range(probe.width)
                    for i, sign in probe.by_x.get((u-beta)*inverse % probe.modulus, ())]
            if any(probe.curve[i][0] in excluded for i, _ in slot):
                continue
            address, _ = probe.singleton([slot, *tail], target)
            if address is None:
                continue
            s = [next(v for v, a in row if a == address//probe.block**k % probe.block)
                 for row, k in zip(visible, (1, 0))]
            key = normal(cross(*s), p)
            witness_mass[key] += 1
            seen_seeds.append((alpha, beta, key))
    assert seen_seeds
    scale = probe.modulus*probe.width//2
    assert max(witness_mass.values()) <= scale
    axes = directions(p, 3)
    max_cover = max_success = 0
    for axis in axes:
        chosen = next(slot for slot in visible if all(normal(row, p) != axis for row, _ in slot))
        roots = {normal(cross(axis, row), p) for row, _ in chosen}
        candidates = {v for v in roots if all(any(dot(row, v) % p == 0 for row, _ in slot) for slot in visible)}
        assert candidates == {v for v in partial if dot(axis, v) % p == 0}
        assert len(candidates) <= len(chosen)
        mass = sum(witness_mass[v] for v in candidates)
        assert mass <= len(chosen)*scale
        max_cover = max(max_cover, len(candidates))
        max_success = max(max_success, mass)
    assert max_cover == 4
    print(f"F_{p}: all {len(axes)} pencils, {p+1} sections each; at most four candidates per pencil, "
          f"{len(partial)} distinct partial sections from 16 tuples.")
    print(f"Among {probe.modulus**2//2} anchor seeds, {len(seen_seeds)} regular guarded singletons; "
          f"largest fixed-pencil hit count {max_success} (bound {4*scale}).")

    # A nonsingular conic through three independent visible candidates. The
    # selected section may depend on both words of the affine anchor seed.
    a, b, c = next((a, b, c) for a, b, c in product(witness_mass, repeat=3) if dot(a, cross(b, c)) % p)
    conic = [normal(tuple(x+(y-x-z)*u+z*u*u for x, y, z in zip(a, b, c)), p) for u in range(p)]+[c]
    assert len(set(conic)) == p+1
    candidates = set(conic) & partial.keys()
    assert 3 <= len(candidates) <= 2*len(tail[0])
    hits = sum(conic[(alpha+beta) % p] == key for alpha, beta, key in seen_seeds)
    assert hits <= sum(witness_mass[v] for v in candidates) <= len(candidates)*scale
    print(f"Seed-selected point on a smooth conic: {len(candidates)} visible candidates (bound 8), "
          f"{hits} actual hits.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_bulk_sums()
    check_regular()
