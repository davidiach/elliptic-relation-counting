"""Finite Walsh prefix-energy and two-block bounds on mixed signed slot addresses."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import prod, sqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_strided_isolation import finite_laws
from fixtures import curve_cycle


def transform(scores):
    """Integer Walsh coefficients. Divide by the input length for the uniform transform."""
    words = list(scores)
    n = len(words)
    assert n and n.bit_count() == 1
    width = 1
    while width < n:
        for start in range(0, n, 2*width):
            for j in range(start, start+width):
                a, b = words[j], words[j+width]
                words[j], words[j+width] = a+b, a-b
        width *= 2
    assert sum(v*v for v in words) == n*sum(v*v for v in scores)
    return words


def character(index, mask):
    return -1 if (index & mask).bit_count() % 2 else 1


def check_prefixes():
    bits, secret = 6, 43
    n = 1 << bits
    for noisy in (False, True):
        score = [character(secret, r)*(-1 if noisy and r % 8 == 0 else 1) for r in range(n)]
        coefficients = transform(score)
        assert coefficients[secret] == (3*n//4 if noisy else n)
        live = [0]
        for length in range(1, bits+1):
            tail = 1 << (bits-length)
            energies = []
            for prefix in range(1 << length):
                lo = prefix*tail
                energy = sum(a*a for a in coefficients[lo:lo+tail])
                conditional = [sum(score[v*tail+y]*character(v, prefix) for v in range(1 << length))
                               for y in range(tail)]
                assert energy == tail*sum(a*a for a in conditional)
                energies.append(energy)
            live = [2*v+bit for v in live for bit in (0, 1) if energies[2*v+bit] >= (n//2)**2]
            assert len(live) <= 4
        assert live == [secret]
    print("Exact prefix energies and Parseval verified for pure and noisy signed-address parities.")


@dataclass
class Block:
    mass: int
    score: list
    coefficients: list
    # Valid signed tuples with (point labels, packed address, block sum).
    marks: list


def make_blocks(laws, curve, address_bits):
    n = len(curve)
    size = 1 << (address_bits*len(laws))
    blocks = []
    for masks in product(*laws):
        score = [-1]*size
        marks = []
        used = set()
        for choices in product(*masks):
            labels = tuple(i for i, _ in choices)
            code = 0
            for _, address in choices:
                code = (code << address_bits) | address
            assert code not in used
            used.add(code)
            subtotal = sum(labels) % n
            # A nonconstant toy score using only these masks and the local mask input.
            # Lookup of the group-law result is charged only to this finite check.
            if subtotal == 0 or curve[subtotal][0] % 2 == 0:
                score[code] = 1
            marks.append((labels, code, subtotal))
        coefficients = transform(score)
        mass = prod(law[mask] for law, mask in zip(laws, masks))
        blocks.append(Block(mass, score, coefficients, marks))
    return blocks


def block_energy(blocks, n, rho, arity):
    scale = len(blocks[0].score)
    seeds = sum(block.mass for block in blocks)
    counts, absolute, squared = (Counter() for _ in range(3))
    for block in blocks:
        for _, code, subtotal in block.marks:
            value = block.coefficients[code]
            counts[subtotal] += block.mass
            absolute[subtotal] += block.mass*abs(value)
            squared[subtotal] += block.mass*value**2
    atom_bound = n**(arity-1)*rho**arity
    means = [Fraction(absolute[g], seeds*scale) for g in range(n)]
    for g in range(n):
        assert Fraction(counts[g], seeds) <= atom_bound
        assert absolute[g]**2 <= counts[g]*squared[g]
    assert sum(a*a for a in means) <= atom_bound
    return means


def check_mixed():
    p, positions, m = 13, (0,), 3
    curve = curve_cycle(p, 4)
    n = len(curve)
    affine, boxes = finite_laws(p, curve, positions)
    bits = len(positions)+1
    left = make_blocks([affine], curve, bits)
    right = make_blocks([boxes, boxes], curve, bits)
    rho = Fraction(1 << len(positions), 1 << p.bit_length())
    lam = n**(m-1)*rho**m
    a = block_energy(left, n, rho, 1)
    b = block_energy(right, n, rho, 2)
    total_seeds = actual_seeds = 0
    good, all_values, selected = Counter(), Counter(), Counter()
    scale = 1 << (m*bits)
    for ib, jb in product(left, right):
        mass = ib.mass*jb.mass
        actual_seeds += 1
        total_seeds += mass
        full = transform([v*w for v in ib.score for w in jb.score])
        tensor = [v*w for v in ib.coefficients for w in jb.coefficients]
        assert full == tensor
        hits = Counter()
        allowed = {}
        for (i,), li, g in ib.marks:
            for (j, k), rj, h in jb.marks:
                code = (li << (2*bits)) | rj
                target = (g+h) % n
                hits[target] += 1
                all_values[target] += mass*abs(full[code])
                subtotal = (i+j) % n
                if subtotal and curve[i][0] != curve[j][0] and curve[subtotal][0] != curve[k][0]:
                    allowed[target] = code
        for target, code in allowed.items():
            if hits[target] == 1:
                selected[target] += mass
                good[target] += mass*abs(full[code])
    exact = []
    for t in range(1, n):
        singleton = Fraction(good[t], total_seeds*scale)
        ordinary = Fraction(all_values[t], total_seeds*scale)
        assert ordinary == sum(a[g]*b[(t-g) % n] for g in range(n))
        assert selected[t] and singleton <= ordinary and ordinary**2 <= lam/n
        exact.append(singleton)
    assert total_seeds == 8192 and any(exact)
    print(f"Field {p}: {total_seeds:,} weighted seeds ({actual_seeds} distinct mask combinations), "
          f"all {scale} joint Walsh coefficients checked on each. Across {n-1} fixed targets, "
          f"singleton-weighted magnitudes {float(min(exact)):.5f}..{float(max(exact)):.5f}, "
          f"upper bound {sqrt(lam/n):.5f}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    check_prefixes()
    check_mixed()
