"""Independent exact Smith elimination and candidate-enumeration barrier."""
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path


def v2(x):
    x = F(x)
    if not x:
        return float('inf')
    def integer(n):
        n = abs(n)
        return (n & -n).bit_length()-1
    return integer(x.numerator)-integer(x.denominator)


def smith_exponents(matrix):
    """Minimum-valuation pivots and integral row/column operations over Z_(2).

    No inherited mode implementation or claimed pivot formula is used.
    """
    a = [[F(x) for x in row] for row in matrix]
    n = len(a)
    exponents = []
    for k in range(n):
        _, i, j = min((v2(a[i][j]), i, j) for i in range(k, n) for j in range(k, n))
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        pivot = a[k][k]
        assert pivot
        exponents.append(v2(pivot))
        for i in range(k+1, n):
            factor = a[i][k]/pivot
            assert v2(factor) >= 0
            for j in range(k, n):
                a[i][j] -= factor*a[k][j]
        for j in range(k+1, n):
            factor = a[k][j]/pivot
            assert v2(factor) >= 0
            for i in range(k, n):
                a[i][j] -= factor*a[i][k]
    assert exponents == sorted(exponents)
    return exponents


def bit_count_sum(k):
    """Sum popcount(q) for 0 <= q < k, without enumerating k integers."""
    total = 0
    bit = 1
    while bit < k:
        cycles, rest = divmod(k, 2*bit)
        total += cycles*bit+max(0, rest-bit)
        bit *= 2
    return total


def weight_sum(k, shift):
    return shift*k+2*k*(k-1)-2*bit_count_sum(k)


def main():
    records = []
    for a in range(1, 6):
        L = 1 << a
        for i in range(a):
            R = 1 << i
            for bit in (0, 1):
                slots = [x for x in range(L) if (x >> i) & 1 == bit]
                for perturbation in (1, 3):
                    T = R+(2*L)*perturbation
                    matrix = [[F(R, T+x-y) for y in slots] for x in slots]
                    actual = smith_exponents(matrix)
                    expected = [4*(k//R)-2*(k//R).bit_count() for k in range(L//2)]
                    assert actual == expected
                    records.append(dict(a=a, i=i, bit=bit, perturbation=perturbation, exponents=actual))
    for k in range(256):
        assert bit_count_sum(k) == sum(q.bit_count() for q in range(k))
        assert weight_sum(k, 7) == sum(7+4*q-2*q.bit_count() for q in range(k))
    barriers = []
    for u in (4, 8, 16, 32, 64, 128):
        D, M = 1 << (2*u), 1 << (10*u)
        Q = M//D
        root = isqrt(u)
        ceiling = root + (root*root < u)
        k = D//ceiling
        shift = 2*u+1
        P = 6*M*D//u
        used = Q*weight_sum(k, shift)
        bits = Q*k
        assert k <= D//2 and used < P and bits <= M
        # Every subset of these Q*k distinct modes is a feasible tuple.
        assert 16*bits*bits*u >= M*M
        barriers.append(dict(u=u, log2_M=10*u, log2_D=2*u, families=str(Q),
                             ranks_per_family=str(k), precision=str(P), used_weight=str(used),
                             max_degree=str(bits), log2_feasible_subset_count=str(bits),
                             normalized_squared=str(F(bits*bits*u, M*M))))
    root = Path(__file__).resolve().parents[2]
    sources = [root/'_review/sharper-complete-ecdl-v6-reproducibility/sharper-complete-ecdl-v6/reproducibility/paper'/name
               for name in ('global_modes.tex', 'entropy.tex', 'model.tex', 'uniform_implementation.tex',
                            'orientation_packet_rounding.tex', 'composition.tex')]
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    result = dict(status='passed', independent_smith_cases=len(records), spectra=records,
                  digit_sum_controls=256, candidate_subset_certificates=barriers, source_hashes=hashes,
                  scope='Independent rational Smith spectra for finite scalar crossing matrices and explicit large feasible subsets of the unrestricted mode universe. This is not a lower bound for actual nonzero C38 support or for algorithms that avoid enumeration.')
    (Path(__file__).parent/'spectrum_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status':'passed', 'independent_smith_cases':len(records), 'subset_certificates':len(barriers)}))


if __name__ == '__main__':
    main()
