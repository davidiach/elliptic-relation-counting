"""Remainders of canonical dyadic hash candidates using short interval products."""

from bisect import bisect_left
from collections import defaultdict
from math import isqrt

from flint import nmod_poly


def translate(f, s, p):
    """Return f(X+s) by divide and conquer, in arbitrary characteristic."""
    if f.degree() < 1 or s % p == 0:
        return f
    linear = nmod_poly([s % p, 1], p)
    powers = {1: linear}
    k = 1
    while 2*k <= f.degree():
        powers[2*k] = powers[k]*powers[k]
        k *= 2

    def apply(g):
        if g.degree() < 8:
            return g.compose(linear)
        split = 1 << (g.degree().bit_length()-1)
        return apply(g.truncate(split)) + powers[split]*apply(g.right_shift(split))

    return apply(f)


def batch_remainders(f, moduli):
    """Remainder tree. Its moduli may share factors."""
    levels = [moduli]
    while len(levels[-1]) > 1:
        layer = levels[-1]
        parents = [layer[i]*layer[i+1] for i in range(0, len(layer)-1, 2)]
        if len(layer) % 2:
            parents.append(layer[-1])
        levels.append(parents)
    residues = [f % levels[-1][0]]
    for layer in reversed(levels[:-1]):
        residues = [residues[i//2] % modulus for i, modulus in enumerate(layer)]
    return residues


def raw_residue(p, modulus, alpha, beta, width, probe):
    """Residue of the raw canonical preimage polynomial, before the curve-lift test."""
    k = probe.degree()
    if k < 1 or width < 1 or width & (width-1) or width > modulus:
        raise ValueError("Use a positive-degree probe and an admissible dyadic width.")
    if modulus < p or modulus & (modulus-1):
        raise ValueError("The hash modulus must be a power of two at least p.")
    step = pow(alpha, -1, modulus)
    start = -step*beta % modulus
    desired = max(1, isqrt(width*k))
    block = min(width, 1 << (desired-1).bit_length())
    baby = sorted((start+j*step) % modulus for j in range(block))
    shifts = defaultdict(list)
    for j in range(0, width, block):
        offset = j*step % modulus
        for s in (offset, offset-modulus):
            lo = block + bisect_left(baby, -s)
            hi = block + bisect_left(baby, p-s)
            while lo < hi:
                if lo % 2:
                    shifts[lo].append(s)
                    lo += 1
                if hi % 2:
                    hi -= 1
                    shifts[hi].append(s)
                lo //= 2
                hi //= 2
    tree = [None]*block + [nmod_poly([(-x) % p, 1], p) for x in baby]
    for j in range(block-1, 0, -1):
        tree[j] = tree[2*j]*tree[2*j+1]
    result = nmod_poly([1], p)
    for j, offsets in shifts.items():
        transformed = [translate(probe, s, p) for s in offsets]
        for s, remainder in zip(offsets, batch_remainders(tree[j], transformed)):
            result = result*translate(remainder, -s, p) % probe
    return result


def slot_gcd(p, modulus, alpha, beta, width, a, c, probe):
    """The exact monic gcd for abscissas with two rational lifts."""
    common = probe.gcd(raw_residue(p, modulus, alpha, beta, width, probe))
    if common.degree() == 0:
        return common
    rhs = nmod_poly([c % p, a % p, 0, 1], p)
    return common.gcd(rhs.pow_mod((p-1)//2, common)-1)
