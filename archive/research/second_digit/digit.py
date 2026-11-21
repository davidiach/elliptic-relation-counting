"""One coefficient of a group-algebra Frobenius correction from runs of equal step."""

from collections import Counter
from itertools import combinations_with_replacement, groupby
from math import comb, factorial


def multiple(point, n, add, identity):
    out = identity
    while n:
        if n & 1:
            out = add(out, point)
        point = add(point, point)
        n >>= 1
    return out


def scalar_digit(weight, prime):
    modulus = prime*prime
    return ((pow(weight, prime, modulus)-weight) % modulus)//prime


def target_digit(runs, step, target, prime, add, identity):
    """Return ((h**prime - Phi(h))/prime)[target] mod prime.

    Each (origin, length, weight) adds weight on origin + j*step for
    0 <= j < length, with multiplicity. Points must be hashable. The
    caller supplies a prime and addition in a finite abelian group.
    """
    runs = list(runs)
    if any(length < 0 for _, length, _ in runs):
        raise ValueError("Run lengths must be nonnegative")
    runs = [(a, length, w) for a, length, w in runs if length and w]
    if not runs:
        return 0
    max_degree = prime*(max(length for _, length, _ in runs)-1)
    stride = multiple(step, prime, add, identity)

    def hits(origin, degree, ramp_length=None):
        if degree > max_degree:
            return 0
        point = add(origin, multiple(step, degree, add, identity))
        total = 0
        for i in range((max_degree-degree)//prime+1):
            if point == target:
                total += 1 if ramp_length is None else min(i+1, ramp_length)
            point = add(point, stride)
        return total

    # The lifted numerator (1-z)H(z) has two formal endpoints per run.
    endpoints = Counter()
    for origin, length, weight in runs:
        endpoints[origin, 0] += weight
        endpoints[origin, length] -= weight
    atoms = [(point, degree, w) for (point, degree), w in endpoints.items() if w]

    # Division by 1-z**prime propagates each term along one progression.
    answer = 0
    for point, degree, weight in atoms:
        value = scalar_digit(weight, prime)
        if value:
            answer += value*hits(multiple(point, prime, add, identity), prime*degree)

    # Equal endpoints have already been combined. A coefficient divisible by
    # prime contributes only to the pure scalar case handled above.
    units = [(point, degree, w) for point, degree, w in atoms if degree <= max_degree and w % prime]
    inverse_factorials = [pow(factorial(j), -1, prime) for j in range(prime)]
    for indices in combinations_with_replacement(range(len(units)), prime):
        if indices[0] == indices[-1]:
            continue
        origin, degree, coefficient = identity, 0, -1  # Wilson's theorem.
        for index, repetitions in groupby(indices):
            count = sum(1 for _ in repetitions)
            point, shift, weight = units[index]
            degree += count*shift
            origin = add(origin, multiple(point, count, add, identity))
            coefficient *= pow(weight, count, prime)*inverse_factorials[count]
            coefficient %= prime
        answer += coefficient*hits(origin, degree)

    # Subtract Phi(H) * D(1-z)/(1-z**prime). The last term occurs
    # only at prime two. Computing it explicitly handles that case.
    denom_digit = [(j, (-1)**j*(comb(prime, j)//prime)) for j in range(1, prime)]
    if scalar_digit(-1, prime):
        denom_digit.append((prime, scalar_digit(-1, prime)))
    for point, length, weight in runs:
        shifted = multiple(point, prime, add, identity)
        for degree, coefficient in denom_digit:
            answer -= weight*coefficient*hits(shifted, degree, length)
    return answer % prime
