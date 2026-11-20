"""Direct zero-multiplicity certificates for two pooled three-slot sign orbits."""

from math import comb

from flint import nmod_poly


def pooled_certificate(p, moments, orders=(8, 9, 10, 12)):
    """Map each prefix length to the compatible zero counts, or None if inconclusive.

    Assumes 16 signed determinants in two full three-slot sign orbits from
    regular elliptic slots of sizes 2, 2, 4 in any order. Only
    the field-valued pooled moments enter the readout. A singleton set certifies
    the count; a set containing only positive integers also certifies existence.
    """
    top = max(orders)
    if p <= 7 or top < 8 or min(orders) < 8 or top > len(moments):
        raise ValueError("Need characteristic >7 and at least eight supplied moments")
    zero, one, h = (nmod_poly(c, p) for c in ([], [1], [0, 1]))
    t = h*h
    inv2 = pow(2, -1, p)
    raw = [1] + [int(v)*pow(16, -1, p) % p for v in moments[:top]]
    c = raw[1]
    centered = [sum(comb(j, i)*raw[i]*pow(-c, j-i, p) for i in range(j+1)) % p
                for j in range(top+1)]

    # The central moment of order 2q in either block is A_q +/- B_q/h.
    A, B = [one], [zero]
    for q in range(1, 4):
        a, b = one*centered[2*q], one*centered[2*q+1]
        for r in range(q):
            a -= comb(2*q, 2*r)*t**(q-r)*A[r]
            b -= comb(2*q+1, 2*r)*t**(q-r)*B[r]
        A.append(a)
        B.append(b*pow(2*q+1, -1, p))

    profiles, jets = [], []
    for sign in (1, -1):
        # Scaling every centered determinant by h clears the Laurent denominators.
        evens = [one] + [h**(2*q)*A[q] + sign*h**(2*q-1)*B[q] for q in range(1, 4)]
        _, s1, u4, u6 = evens
        s2 = (3*s1*s1-u4)*inv2
        s3 = (u6-15*u4*s1+30*s1**3)*pow(16, -1, p)
        e2 = (s1*s1-s2)*inv2
        e3 = (s1**3-3*s1*s2+2*s3)*pow(6, -1, p)
        const = s1*s1-4*e2
        # Coefficients of the scaled block octic at z^0,z^2,...,z^8.
        octic = [const*const, -4*s1*const-64*e3,
                 4*s1*s1+2*const, -4*s1, one]
        for q in range(4, top//2+1):
            evens.append(-sum((octic[r]*evens[q-4+r] for r in range(4)), zero))
        profiles.append(evens)
        minus_center = -h*c - sign*t
        jets.append([sum((comb(2*j, r)*octic[j]*minus_center**(2*j-r)
                          for j in range((r+1)//2, 5)), zero) for r in range(5)])
    # The two smaller slots give four ordered pairs. Each has at most one
    # completing point. Thus only coefficients 0 through 4 are ever needed.
    at_zero = [sum((jets[0][i]*jets[1][r-i] for i in range(r+1)), zero)
               for r in range(5)]
    answers = {}
    g = zero
    for j in range(8, top+1):
        predicted = sum((comb(j, 2*r)*(sign*t)**(j-2*r)*profile[r]
                         for sign, profile in zip((1, -1), profiles)
                         for r in range(j//2+1)), zero)*inv2
        residual = h**j*centered[j] - predicted
        assert residual.degree() <= 2*j
        if not residual.is_zero():
            g = g.gcd(residual)
        if j not in orders:
            continue
        if g.is_zero() or not any(centered[1:j+1:2]):
            # Equal centers have not been excluded, or all separations survive.
            answers[j] = None
            continue
        admissible = (g.gcd(pow(h, p-1, g)-one) if g.degree() else one)
        counts = set()
        for r, coeff in enumerate(at_zero):
            if admissible.degree() <= 0:
                break
            remaining = admissible.gcd(coeff)
            if remaining.degree() < admissible.degree():
                counts.add(r)
            admissible = remaining
        if not counts:
            raise ValueError("Moment vector is inconsistent with the promised two orbits")
        answers[j] = frozenset(counts)
    return answers
