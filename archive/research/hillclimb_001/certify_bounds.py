"""Exact rational certificates for conditional composition constants.

The inherited contraction/encoding theorems and cited matrix existence
results are premises. Passing this program does not prove those premises.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial, isqrt, comb
import common
import intervals as I
import adaptive_pair_join as joins


def support(m, radius2):
    assert m >= 0 and radius2 >= 0
    if m*m < radius2/6:
        # The maximum is inside the interval, not at the boundary. This
        # branch matters when changing the inherited matrix parameters.
        return I.sqrt(I.point(18*m*m+radius2))
    return I.add(I.mul(I.point(m), I.sqrt(I.point(F(27, 2)))),
                 I.mul(I.point(F(1, 2)), I.sqrt(I.point(radius2))))


def candidate(alpha, omega):
    beta = (omega-2)/(1-alpha)
    d = 1-alpha*beta
    assert 0 < beta < d < 1
    assert beta*beta+d*d < 1
    assert beta*beta+7*d*d+7 < 14
    sqrt14 = I.sqrt(I.point(14))
    exceptional = ((beta, F(1), d), (d, F(1), beta))
    forms = []
    for a, b, c in permutations((beta, d, F(1))):
        m = (2*b+c)/3
        r2 = 2*a*a+4*(b-c)**2
        bound = support(m, r2)
        if (a, b, c) not in exceptional:
            assert bound[1] < sqrt14[0]
        forms.append(dict(coefficients=list(map(str, (a,b,c))),
                          m=str(m), radius_squared=str(r2),
                          maximum='boundary' if m*m>=r2/6 else 'interior',
                          bound=I.decimal_interval(bound), exceptional=(a,b,c) in exceptional))
    # Derive all eight averaged forms from their six coefficients, rather
    # than importing the closed k-formula from the old certificate.
    averages = []
    for bits in product((0, 1), repeat=3):
        internal_coefficients = [beta if bit else d for bit in bits]
        external_coefficients = [d if bit else beta for bit in bits]
        aa = [(internal_coefficients[i]+internal_coefficients[j])/3 for i,j in combinations(range(3),2)]
        tt = [(2+x)/3 for x in external_coefficients]
        m = sum(tt)/3
        r2 = sum(a*a for a in aa)+6*sum((x-m)**2 for x in tt)
        interval = support(m, r2)
        averages.append(interval)
    ceiling = (max([sqrt14[0]]+[x[0] for x in averages]),
               max([sqrt14[1]]+[x[1] for x in averages]))
    nine = I.parameters(9)
    K = I.mul(I.div(ceiling, I.point(6)), nine['B'])
    J = I.mul(I.div(sqrt14, I.point(6)), nine['B'])
    assert nine['child'][1] < J[0] <= K[1]
    general = I.mul(I.point(2), I.sqrt(I.mul(I.point(19), K)))
    c44 = I.mul(I.point(2), I.sqrt(I.mul(I.point(22), K)))
    return dict(alpha=str(alpha), omega=str(omega), beta=str(beta), d=str(d),
                six_support_forms=forms,
                eight_averages=[I.decimal_interval(x) for x in averages],
                constants={name:I.decimal_interval(x) for name,x in
                           (('local_join',ceiling),('K',K),('J',J),('child',nine['child']),
                            ('C38_simplex_time',general),('C38_simplex_space',I.div(general,I.point(2))),
                            ('C44_simplex_time',c44))}), K


def partition_certificate():
    matrix = [[0]*36 for _ in range(36)]
    for bundles in joins.BUNDLES:
        for bundle in bundles:
            for e in bundle:
                for f in bundle:
                    matrix[e][f] += 1
    for e in range(36):
        assert sum(matrix[e]) == 1890
        for f in range(36):
            expected = 210 if e==f else 60 if set(joins.EDGES[e]) & set(joins.EDGES[f]) else 40
            assert matrix[e][f] == expected
    return dict(partitions=280, matrix_entries=1296, average_row_sum='27/4',
                proof='For nonnegative symmetric M, 2xy<=x^2+y^2 gives w^T M w <= max_row_sum*||w||^2.')


def finite_schedule(lam, t):
    S = 1 << lam
    gamma = (lam**12-1).bit_length()
    previous_a = 0
    mass = F(S)
    threshold = S//lam
    p = 2*t*S+1
    rows = []
    while mass > threshold:
        lo, hi = 0, lam
        while lo < hi:
            u = (lo+hi+1)//2
            if (1 << (previous_a+gamma))*90**u*lam**48 <= S:
                lo = u
            else:
                hi = u-1
        u = lo
        if u < 4:
            return dict(lam=lam,t=t,fallback=True,active_levels=len(rows),reason='coordinate guard')
        a = previous_a+u+gamma
        assert (1 << a)*45**u*lam**48 <= S
        assert (1 << (previous_a+gamma))*90**(u+1)*lam**48 > S
        ceil_mass = -(-mass.numerator//mass.denominator)
        # P_i upper bound uniformly over all retained valuations, using
        # q_i <= 2tM and exact integer sum j/2^j < 2.
        nu = (factorial(t) & -factorial(t)).bit_length()-1
        qmax = 2*t*ceil_mass
        child_p_bound = p+nu*ceil_mass+2*(1 << a)*qmax+qmax*(a+1)+1
        assert child_p_bound <= 32*t*S*lam*lam*(1 << a)
        rows.append(dict(u=u,a=a,packet_gap=lam-a))
        p = child_p_bound
        mass *= F(2*(lam+1),9*lam)
        previous_a = a
    return dict(lam=lam,t=t,fallback=False,active_levels=len(rows),
                min_u=min(x['u'] for x in rows), min_packet_gap=min(x['packet_gap'] for x in rows))


def main():
    common.verify_sources()
    baseline, K0 = candidate(F(8,25),F(239,100))
    published, Kpub = candidate(F(32133,100000),F(118567,50000))
    recent, Krecent = candidate(F(1,2),F(9,4))
    assert Krecent[1] < Kpub[0] < Kpub[1] < K0[0]
    assert F(published['beta']) == joins.BETA and F(published['d']) == joins.D
    # Explicit counterexample to reusing the old monotonicity shortcut with
    # new matrix parameters: an interior point exceeds the endpoint bound.
    m, radius2, position = F(7,12), F(9,4), F(18,5)
    assert m*m < radius2/6
    endpoint = I.add(I.mul(I.point(m),I.sqrt(I.point(F(27,2)))),I.point(F(3,4)))
    interior = I.add(I.point(m*position),I.sqrt(I.point(radius2*(1-position*position/18))))
    assert interior[0] > endpoint[1]
    scan = []
    for alpha in (F(3,10),F(31,100),F(8,25),F(321,1000),F(3213,10000),F(32133,100000)):
        for omega in (F(239,100),F(238,100),F(2375,1000),F(2372,1000),F(118567,50000)):
            try:
                data, K = candidate(alpha,omega)
            except AssertionError:
                scan.append(dict(alpha=str(alpha),omega=str(omega),status='certificate_rejected'))
            else:
                scan.append(dict(alpha=str(alpha),omega=str(omega),status='certified_support',K=I.decimal_interval(K)))
                assert Kpub[1] <= K[1]
    schedules = [finite_schedule(lam,t) for lam in (64,256,1023,1024,1025,2047,2048,2049,4095,4096,4097,16384) for t in (1,2,3,8,16)]
    assert all(not row['fallback'] for row in schedules if row['lam']==16384)
    alphabet = []
    for t in (1,2,3,4,8,16,32,64,128):
        Q=comb(2*t,t)
        assert Q <= 1 << (2*t)
        normalized=I.div(I.point(t),I.log2_integer(Q))
        assert normalized[0] >= F(1,2)
        alphabet.append(dict(t=t,Q=str(Q),per_bit_factor=I.decimal_interval(normalized)))
    # Compare with a deliberately outward enclosure of the supplied simplex
    # coefficient, not with rounded point data used as an exact theorem.
    old = (F('32.15695165969816'), F('32.15695165969817'))
    new = I.mul(I.point(2), I.sqrt(I.mul(I.point(19),Kpub)))
    gain = I.mul(I.point(100), I.sub(I.point(1),I.div(new,old)))
    common.save('bounds.json',dict(status='passed',counts=dict(matrix_parameter_candidates=len(scan),
        integer_schedule_cases=len(schedules),selector_sizes=len(alphabet)),
        partition=partition_certificate(),original_matrix_input=baseline,
        published_matrix_input=published,recent_preprint_input=recent,
        matrix_parameter_scan=scan,finite_schedules=schedules,finite_alphabet_penalties=alphabet,
        gain_percent_vs_supplied_simplex=I.decimal_interval(gain,digits=13),
        rejected_monotonicity_shortcut=dict(m=str(m),radius_squared=str(radius2),
            tested_interior_point=str(position),interior_value=I.decimal_interval(interior),
            endpoint_value=I.decimal_interval(endpoint)),
        premises=['Inherited finite contraction, padding, entropy and encoding proofs',
                  'Uniform simplex lifting with per-call scalar precision, as extended in proof.tex',
                  'Published rational matrix-formula existence; recent-preprint row has an additional dependency'],
        scope='Exact support inequalities, rational numerical enclosures, finite schedule guards and conditional transfer arithmetic. This is not a formal proof of the inherited foundations or an implementation of their asymptotic runtime.'))


if __name__ == '__main__':
    main()
