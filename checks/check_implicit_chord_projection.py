"""Check a capped predictor of an unshifted product comparison over whole target populations."""
from collections import Counter
from itertools import product
from math import prod

from check_chord_difference import add_with_slope, literal, negative
from probe_signed_moment import Probe


def extend(table, slot, p, linear, cap=None):
    merged = {}
    for point, (n, u, c, uc) in table.items():
        for item, s in slot:
            new_point, slope = add_with_slope(point, item, p, linear)
            old = merged.get(new_point, (0, 0, 0, 0))
            merged[new_point] = (old[0]+n, (old[1]+u+n*s)%p,
                                 (old[2]+c+n*slope)%p,
                                 (old[3]+uc+s*c+slope*u+n*s*slope)%p)
            if cap is not None and len(merged) > cap:
                return None
    return merged


def prepare(slots, p, linear):
    active = {P for slot in slots for P, speed in slot if speed}
    if not active or not all(slots):
        return "equal", {}
    index = next(i for i, slot in enumerate(slots) if any(speed for _, speed in slot))
    cap = 1+len(active)
    table = {None: (1, 0, 0, 0)}
    for i, slot in enumerate(slots):
        if i != index:
            table = extend(table, slot, p, linear, cap)
            if table is None:
                return "pole", {}
    return "exact", extend(table, slots[index], p, linear)


def decide(mode, table, slots, endpoint, p, linear):
    if mode == "equal":
        return True
    if mode == "pole":
        return False
    final = extend(table, [(endpoint, 0)], p, linear)
    n, u, _, _ = final.get(None, (0, 0, 0, 0))
    if n:
        return n > 1 or u == 0
    weighted_chords = sum(values[3] for values in final.values())
    with_endpoint = [*slots, [(endpoint, 0)]]
    size = prod(map(len, with_endpoint))
    for i, a in enumerate(with_endpoint):
        for b in with_endpoint[:i]:
            weight = size//(len(a)*len(b))
            for P, u in a:
                for Q, v in b:
                    K = add_with_slope(P, negative(Q, p), p, linear)[1]
                    weighted_chords += weight*(u-v)*K
    return weighted_chords % p == 0


def check_population(p, linear):
    probe = Probe(p, linear, ())
    xs = sorted(x for x in range(p) if len(probe.by_x[x]) == 2)[:3]
    raw = [probe.by_x[x] for x in xs]
    eligible = [t for t in range(1, probe.n) if probe.curve[t][0] not in xs]
    seed_volume = prod(map(len, raw))
    arity = len(raw)
    assert arity*(arity+1)*seed_volume < p
    target_counts = Counter(sum(i for i, _ in word) % probe.n for word in product(*raw))
    totals = Counter()
    for name, mask in (("zero", 0), ("one ordinate", 1), ("three ordinates", 7)):
        slots = [[(probe.curve[i], sign if mask & (1 << j) else 0) for i, sign in slot]
                 for j, slot in enumerate(raw)]
        mode, table = prepare(slots, p, linear)
        assert mode == {"zero": "equal", "one ordinate": "pole", "three ordinates": "exact"}[name]
        active = {P for slot in slots for P, s in slot if s}
        residues = Counter()
        for chosen in product(*raw):
            speed = sum(s for j, (_, s) in enumerate(chosen) if mask & (1 << j))
            group_sum = sum(i for i, _ in chosen) % probe.n
            if group_sum:
                residues[(-group_sum) % probe.n] += speed
            for j, (i, s) in enumerate(chosen):
                if mask & (1 << j):
                    residues[i] -= s
        if mode == "pole":
            assert any(value % p and probe.curve[i] not in active for i, value in residues.items())
        errors = zero_errors = singleton_errors = singletons = 0
        truth = {}
        shifted_errors = [0]*arity
        shifted_empty_errors = [0]*arity
        simultaneous_singleton_errors = 0
        zero_sum_speeds = [sum(s for j, (_, s) in enumerate(word) if mask & (1 << j))
                           for word in product(*raw) if sum(i for i, _ in word)%probe.n == 0]
        for j in range(1, arity+1):
            finite = Counter(residues)
            for word in product(*raw):
                total = sum(i for i, _ in word)%probe.n
                if total:
                    finite[-total%probe.n] -= j
                for i, _ in word:
                    finite[i] -= j
            residue_at_infinity = -sum(finite.values())
            assert residue_at_infinity == j*(arity+1)*seed_volume+sum(k-j for k in zero_sum_speeds)
            assert j*arity*seed_volume <= residue_at_infinity < p
        for t in eligible:
            endpoint = probe.curve[-t]
            actual = literal([*slots, [(endpoint, 0)]], p, linear, 4)
            base_slots = [[(P, 0) for P, _ in slot] for slot in slots]
            baseline = literal([*base_slots, [(endpoint, 0)]], p, linear, 4)
            truth[t] = actual == baseline
            predicted = decide(mode, table, slots, endpoint, p, linear)
            wrong = predicted != truth[t]
            errors += wrong
            zero_errors += wrong and target_counts[t] == 0
            singleton_errors += wrong and target_counts[t] == 1
            singletons += target_counts[t] == 1
            observed_shifts = [0] if truth[t] else []
            for j in range(1, arity+1):
                matched = actual == literal([*base_slots, [(endpoint, j)]], p, linear, 4)
                if matched:
                    observed_shifts.append(j)
                    shifted_errors[j-1] += 1
                    shifted_empty_errors[j-1] += target_counts[t] == 0
            if target_counts[t] == 1:
                assert len(observed_shifts) == 1
                simultaneous_singleton_errors += (observed_shifts != ([0] if predicted else []))
        assert max(shifted_empty_errors) <= seed_volume+sum(map(len, raw))+1
        if mode == "pole":
            assert zero_errors <= seed_volume+len(active)+1
            assert errors <= 2*seed_volume+len(active)+1
            assert 2*singleton_errors == singletons > 0
            assert all(truth[t] != truth[-t % probe.n] for t in eligible if target_counts[t] == 1)
            assert simultaneous_singleton_errors == singletons
        else:
            assert not errors
        totals[mode] += 1
        print(f"  {name}, {mode}: errors {errors}/{len(eligible)} overall, "
              f"{zero_errors} on empty targets, {singleton_errors}/{singletons} on singletons.")
        print(f"    Positive shifts j=1,2,3 had {shifted_errors} errors; "
              f"all-shift failure on {simultaneous_singleton_errors}/{singletons} singleton targets.")
    assert len(totals) == 3


def check_stabilizer():
    n = 114
    subgroup = {0, 38, 76}
    selected = {(a+x)%n for a in (-2, 2) for x in subgroup}
    stabilizer = {z for z in range(n) if {(a+z)%n for a in selected} == selected}
    assert stabilizer == subgroup
    cases = 0
    for left, right in product(({0}, {38, 76}, {1, 113}, {2, 112}), repeat=2):
        predicted = all((x+y)%n in stabilizer for x, y in product(left, right))
        projected = Counter(min(z, -z%n) for a, x, y in product(selected, left, right)
                            if (z := (a+x+y)%n))
        slot = Counter({x: v*len(left)*len(right) for x, v in
                        Counter(min(a, -a%n) for a in selected).items()})
        assert predicted == (projected == slot)
        cases += 1
    print(f"  Group order {n}: checked {cases} abscissa-multiset identities with stabilizer order three.")


def check_walsh_transcript():
    p, linear = 8191, 1
    probe = Probe(p, linear, (0,))
    raw = []
    for start in range(0, p-1, 2):
        slot = [(point, 2*j+s) for j in (0, 1) for point, s in probe.by_x[start+j]]
        if len(slot) == 4:
            raw.append(slot)
        if len(raw) == 4:
            break
    arity = len(raw)
    assert arity == 4
    volume = prod(map(len, raw))
    assert arity*(arity+1)*volume < p
    for i in range(arity):
        sums = {sum(point for point, _ in word)%probe.n for word in product(*(raw[:i]+raw[i+1:]))}
        assert len(sums) > 1+sum(map(len, raw))
    used_x = {probe.curve[i][0] for slot in raw for i, _ in slot}
    counts = Counter(sum(i for i, _ in word)%probe.n for word in product(*raw))
    witness = None
    for word in product(*raw):
        index = sum(i for i, _ in word)%probe.n
        addr = sum(a << (2*j) for j, (_, a) in enumerate(word))
        if not index or counts[index] != 1 or probe.curve[index][0] in used_x or not addr:
            continue
        prefix = word[0][0]
        for point, _ in word[1:]:
            if not prefix or probe.curve[prefix][0] == probe.curve[point][0]:
                break
            prefix = (prefix+point)%probe.n
        else:
            witness = addr
            break
    assert witness is not None
    scores = []
    for r in range(256):
        slots = [[(probe.curve[i], ((a << (2*j)) & r).bit_count()%2) for i, a in slot]
                 for j, slot in enumerate(raw)]
        mode, _ = prepare(slots, p, linear)
        assert mode == ("pole" if r else "equal")
        scores.append(1 if mode == "equal" else -1)
    walsh_numerator = sum(val*(-1)**((r & witness).bit_count()%2) for r, val in enumerate(scores))
    assert walsh_numerator == 2
    print(f"  F_{p}: all 256 challenge transcripts on a regular guarded singleton; "
          f"measured nonzero witness Walsh coefficient 2/{volume}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    print("F_101, cyclic order 103: comparing with both literal precision-two products on each target.")
    check_population(101, 17)
    check_stabilizer()
    check_walsh_transcript()
