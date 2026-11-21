"""Check the short hash interfaces, selected target counts, and asymptotic word ledger."""

from copy import deepcopy
from collections import Counter
from fractions import Fraction
from itertools import product
from math import isqrt, prod
from pathlib import Path
from random import Random
import sys

from hash_blocks import InverseHash, Layout, arith, ceildiv, common


def main_x_rows(main, x):
    original = arith.H(arith.limbs(x, main), main)
    return {j*main["L2"]+row: value
            for row, words in enumerate(original)
            for j, value in enumerate(arith.H_inverse_row(words, row, main))}


def rejects(fn):
    try:
        fn()
    except AssertionError:
        return
    raise AssertionError("Corrupt input was accepted")


def retained_c38_checks():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hillclimb_003"))
    from c38_slices import Compiled, sc
    from factor_network import EquationClosure, EquationMiddle, equation_dependencies

    rng = Random(2026111702)
    for b, grid, difference, m, dimension in ((4, 8, 3, 2, 1), (7, 16, 27, 5, 2)):
        data, values, _ = arith.make_instance(b, grid, difference, m, rng, "general")
        compiled = Compiled(data, dimension)
        encoded = compiled.encode_witness(arith.witness(data, values))
        assert compiled.accepts(encoded)
        par, eq_count = data["par"], arith.arithmetic_count(m)
        # Remove every old hash-only wire before evaluating the retained predicates.
        retained = {e: symbol for e, symbol in encoded.items()
                    if not (compiled.paths[e][0] in ("outs", "carry_ports")
                            and compiled.paths[e][1] >= eq_count)}
        assert sum(compiled.paths[e][0] == "outs" for e in encoded if e not in retained) == (
            m * grid * compiled.codec.width)
        ring = sc.Ring(1)
        deps = equation_dependencies(m)
        x_rows = [{} for _ in range(m)]
        for tensor in compiled.net.tensors:
            if not hasattr(tensor, "kind"):
                if tensor.ports[0] in retained:
                    assert tensor.entry(tuple(retained[e] for e in tensor.ports), ring) == ring.one
            elif tensor.kind == "first":
                symbols = tuple(retained[e] for e in tensor.ports)
                assert tensor.entry(symbols, ring) == ring.one
                ass = tensor.local_assignment(symbols)
                row = tensor.index
                for j in range(m):
                    words = [arith.read_word(w, par) for w in ass["hs"][f"x{j+1}"][row]]
                    for i, digit in enumerate(arith.H_inverse_row(words, row, par)):
                        x_rows[j][i*par["L2"]+row] = digit
            else:
                # Individual checks inside one retained vertex can share their input wires.
                for eq in range(eq_count):
                    if tensor.kind == "middle":
                        binds = [d for d in tensor.bindings if (
                            d.path[0] == "hs" and d.path[1] in deps[eq]) or (
                            d.path[0] == "outs" and d.path[1] == eq)]
                        check = EquationMiddle(tensor, eq, deps[eq], binds)
                    else:
                        binds = [d for d in tensor.bindings if d.path[1] == eq]
                        check = EquationClosure(tensor, eq, binds)
                    assert check.entry(tuple(retained[e] for e in check.ports), ring) == ring.one
        layout = Layout(b, par["k"], data["hash_ell"], grid, par["L2"], 2)
        for j in range(m):
            x = values[f"x{j+1}"]
            alpha, beta = (data["constants"][f"{name}{j+1}"] for name in ("alpha", "beta"))
            g = InverseHash(par, layout, alpha, beta)
            t = (alpha*x+beta) % (1 << b)
            boundary = g.witness(t)
            for output in (boundary[0], *boundary[1]):
                for row in output:
                    for word in row:
                        assert compiled.codec.decode(compiled.codec.encode(word)) == word
            g.check(boundary, x_rows[j])
            # Keep the same valid core assignment but choose a hash excluding x.
            bad = InverseHash(par, layout, alpha, (-1-alpha*x) % (1 << b))
            for t in range(1 << layout.ell):
                rejects(lambda: bad.check(bad.witness(t), x_rows[j]))
        print(f"Retained C38 predicates: b={b}, m={m}, alphabet={compiled.codec.Q}, "
              f"old hash-only symbols removed={len(encoded)-len(retained)}, "
              f"replacement hash slots checked={m}")


def hash_checks():
    main = arith.parameters(4, 8, 3, "general")
    layout = Layout(4, main["k"], 2, 8, main["L2"], 2)
    xs = [main_x_rows(main, x) for x in range(16)]
    count = 0
    for alpha, beta in product(range(1, 16, 2), range(16)):
        gadget = InverseHash(main, layout, alpha, beta)
        for t in range(4):
            ports = gadget.witness(t)
            for x in range(16):
                if (alpha*x+beta) % 16 == t:
                    gadget.check(ports, xs[x])
                else:
                    rejects(lambda: gadget.check(ports, xs[x]))
                count += 1
    print(f"Exhaustive 4-bit forward/inverse tests: {count}")

    rng = Random(2026111701)
    checks = corruptions = overflows = 0
    for b, ell, grid, nfft, run in ((73, 17, 128, 32, 4), (136, 20, 128, 32, 2),
                                      (111, 9, 128, 32, 8), (151, 19, 128, 64, 4)):
        main = arith.parameters(b, grid, 1, "general")
        layout = Layout(b, main["k"], ell, nfft, main["L2"], run)
        for trial in range(12):
            modulus = 1 << b
            alpha, beta = rng.randrange(1, modulus, 2), rng.randrange(modulus)
            gadget = InverseHash(main, layout, alpha, beta)
            t = rng.randrange(1 << ell)
            inverse = pow(alpha, -1, modulus)
            untrimmed = inverse*t + (-inverse*beta) % modulus
            x = untrimmed % modulus
            assert (alpha*x+beta) % modulus == t
            overflows += untrimmed >= modulus
            ports = gadget.witness(t)
            rows = main_x_rows(main, x)
            gadget.check(ports, rows)
            source_limbs = arith.encode(t, layout.short, gadget.radix)
            padded = source_limbs + [0]*(nfft-layout.short)
            assert arith.H_inverse(ports[0], gadget.par) == padded
            for block, boundary in zip(gadget.public_blocks, ports[1]):
                assert arith.H_inverse(boundary, gadget.par) == arith.conv(block, padded)
            if trial == 0:
                changed = deepcopy(ports)
                changed[1][0][0][0] = (changed[1][0][0][0]+1) % gadget.modulus
                rejects(lambda: gadget.check(changed, rows))
                changed = deepcopy(ports)
                changed[2][run] ^= 1
                rejects(lambda: gadget.check(changed, rows))
                changed = deepcopy(ports)
                changed[0][0][0] = gadget.modulus
                rejects(lambda: gadget.check(changed, rows))
                rejects(lambda: gadget.check(ports, main_x_rows(main, x ^ 1)))
                corruptions += 4
                if ell % main["k"]:
                    bad_t = 1 << ell
                    bad_ports = gadget.witness(bad_t)
                    bad_x = (inverse*(bad_t-beta)) % modulus
                    rejects(lambda: gadget.check(bad_ports, main_x_rows(main, bad_x)))
                    corruptions += 1
            checks += 1
    assert overflows
    print(f"Block-interface tests: {checks}, corruptions rejected: {corruptions}, "
          f"discarded high overflow: {overflows}")


def check_extended_equations(factors, extra, p, a):
    m = len(factors)
    partial = factors[0]
    sources = {f"x{i}": pt[0] for i, pt in enumerate(factors, 1)}
    for j, point in enumerate(factors[1:], 2):
        assert point[0] != partial[0]
        sources[f"l{j}"] = (point[1]-partial[1])*pow(point[0]-partial[0], -1, p) % p
        partial = arith.point_add(partial, point, p, a)
        if j < m:
            sources[f"X{j}"] = partial[0]
            if j % 2 == 0:
                sources[f"Y{j}"] = partial[1]
    constants = dict(Xm=partial[0], Ym=partial[1], a=a)
    for j in range(1, m+1):
        constants[f"alpha{j}"], constants[f"beta{j}"] = 1, 0
    assert all(z % p == 0 for z in arith.operations(sources, constants, m)[:arith.arithmetic_count(m)])
    y_first = sources["l2"]*(sources["x1"]-sources["X2"])-sources["Y2"]
    y_last = sources[f"l{m}"]*(sources[f"x{m}"]-partial[0])-partial[1]
    assert sources["x1"] == sources[f"x{m}"] and (y_first+y_last) % p == 0
    if extra is not None:
        j = m-1
        if j % 2 == 0:
            y = sources[f"l{j}"]*(sources[f"x{j}"]-sources[f"X{j}"])-sources[f"Y{j}"]
        else:
            y = sources[f"l{j}"]*(sources[f"x{j}"]-sources[f"X{j-1}"])+sources[f"Y{j-1}"]
        assert (sources[f"x{j}"], y % p) == extra


def target_checks():
    p, a = 211, 1
    points = arith.points(p, a, 1)
    size = len(points)+1
    assert size == 223
    main = arith.parameters(8, 16, 45, "general")
    layout = Layout(8, main["k"], 1, 8, main["L2"], 2)
    def add(u, v):
        return arith.point_add(u, v, p, a)
    def neg(u):
        return (u[0], -u[1] % p)
    def guarded(factors):
        out = factors[0]
        for point in factors[1:]:
            if out is None or out[0] == point[0]:
                return False, None
            out = add(out, point)
        return True, out
    def slot(alpha, beta):
        g = InverseHash(main, layout, alpha, beta)
        xs = []
        for t in range(2):
            x = pow(alpha, -1, 256)*(t-beta) % 256
            if x < p:
                g.check(g.witness(t), main_x_rows(main, x))
                xs.append(x)
        result = [pt for pt in points if pt[0] in xs]
        assert result == [pt for pt in points if (alpha*pt[0]+beta) % 256 < 2]
        return result
    P, Q = points[0], points[8]
    assert P[0] != Q[0]
    S = add(P, Q)
    cyclic = [P, Q, neg(S)]
    specs = [(h, (-h*point[0]) % 256) for h, point in zip((5, 15, 23), cyclic)]
    cases = [(specs, None), (specs, add(S, S)), ([(1, 3), *specs[1:]], Q),
             (specs[:2], S),
             ([(5, -5*P[0] % 256)]*2, add(P, P)),
             ([(5, -5*P[0] % 256)]*2, None)]
    for parameters, target in cases:
        sets = [slot(*s) for s in parameters]
        m = len(sets)
        extra = Q if target is None else None
        endpoint = Q if target is None else target
        literal = original_guarded = extended = checked_equations = 0
        for factors in product(*sets):
            total = None
            for point in factors:
                total = add(total, point)
            literal += total == target
            valid, out = guarded(factors)
            original_guarded += valid and out == target
            for offset in points:
                path = (offset, *factors, *((extra,) if extra is not None else ()), neg(offset))
                valid, out = guarded(path)
                if valid and out == endpoint:
                    extended += 1
                    if checked_equations < 8:
                        check_extended_equations(path, extra, p, a)
                        checked_equations += 1
        upper = min(prod(len(s) for j, s in enumerate(sets) if j != i) for i in range(m))
        assert (2*m+7)*upper < size
        assert 0 <= literal*size-extended <= (2*m+7)*literal
        assert ceildiv(extended, size) == literal
        print(f"Target={target}, slots={list(map(len,sets))}, ordinary={literal}, "
              f"original guarded={original_guarded}, extended={extended}, recovered={ceildiv(extended,size)}")


def ledger_checks():
    # Exact size arithmetic only. These are distinct from finite transform executions.
    for m in (16, 64, 256, 1024, 4096):
        log_b = m
        b = 1 << log_b
        k = log_b**2
        n = ceildiv(b, k)
        main_n = 1 << ((2*n+4)-1).bit_length()
        v = 1 << ((main_n.bit_length())//2)
        ell = ceildiv(b, m) + 2
        nfft = 1 << (ceildiv(n, isqrt(m))-1).bit_length()
        nfft = max(nfft, 2*v)
        layout = Layout(b, k, ell, nfft, v, 1)
        assert layout.payload+layout.short-1 <= nfft
        assert 0 <= layout.payload*layout.blocks-n < layout.payload
        ratio = Fraction((layout.blocks+1)*nfft, n)
        assert ratio >= 1
        print(f"m={m}, (forward+inverse hash words)/number of long limbs={float(ratio):.6f}")


def all_targets_with_two_torsion():
    p, a = 71, 2
    points = arith.points(p, a, 1)
    group = (None, *points)
    torsion = [point for point in points if point[1] == 0]
    assert len(torsion) == 3
    # Two identical dyadic windows include a rational two-torsion point.
    fixed = next(pt for pt in points if pt[1] and (pt[0]-torsion[0][0]) % 2)
    alpha = pow(fixed[0]-torsion[0][0], -1, 128)
    slot = [pt for pt in points if (alpha*(pt[0]-torsion[0][0])) % 128 < 2]
    assert len(slot) == 3 and 11*len(slot) < len(group)
    expected, extended = Counter(), Counter()
    for first, second in product(slot, repeat=2):
        total = arith.point_add(first, second, p, a)
        expected[total] += 1
        for use_fixed in (False, True):
            body = (first, second, fixed) if use_fixed else (first, second)
            for offset in points:
                partial = offset
                for factor in (*body, (offset[0], -offset[1] % p)):
                    if partial[0] == factor[0]:
                        break
                    partial = arith.point_add(partial, factor, p, a)
                else:
                    if use_fixed:
                        if partial == fixed:
                            extended[None] += 1
                    else:
                        extended[partial] += 1
    for target in group:
        assert ceildiv(extended[target], len(group)) == expected[target]
    print(f"All {len(group)} targets checked on curve with full rational two-torsion, "
          f"slot size={len(slot)}, identity count={expected[None]}")


if __name__ == "__main__":
    common.verify_sources()
    hash_checks()
    retained_c38_checks()
    target_checks()
    all_targets_with_two_torsion()
    ledger_checks()
