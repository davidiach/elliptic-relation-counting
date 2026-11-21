"""Decode small regular pools using the target quartic and its odd Gram block.

Requires python-flint and NumPy. Run from the project root with assertions on.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product

import numpy as np
from flint import nmod_mat, nmod_poly

from probe_target_jets import curve_cycle, gram_data, hash_families, moment_polynomial, trial


def curve_value(x, a, c, p):
    return (x**3 + a*x + c) % p


def roots(poly):
    return sorted(int(value) for value, multiplicity in poly.roots())


def pair_sum(u, v, a, c, p):
    """Polynomial for the two abscissas x(P+R), x(P-R), assuming u != v."""
    return nmod_poly([(u*v-a)**2-4*c*(u+v), -2*((u+v)*(u*v+a)+2*c), (u-v)**2], p)


@lru_cache(maxsize=8192)
def kernel(u, v, a, c, p):
    """Coefficients of K(u,v,z,X) by interpolating in z at 0,...,4."""
    X = nmod_poly([0, 1], p)
    f = X**3 + a*X + c
    d = (u-v)**2
    rows = []
    for z in range(5):
        g = d*(u-z)**2*(v-z)**2
        poly = (curve_value(u,a,c,p)*(v-z)**2*(X-v)**2*(X-z)**2
                + curve_value(v,a,c,p)*(u-z)**2*(X-u)**2*(X-z)**2
                + curve_value(z,a,c,p)*d*(X-u)**2*(X-v)**2 + f*g)
        rows.append([int(poly[j]) for j in range(5)] + [g])
    vandermonde = nmod_mat([[pow(z,j,p) for j in range(5)] for z in range(5)], p)
    matrix = (vandermonde.inv() * nmod_mat(rows,p)).transpose()
    return matrix.tolist() + [[1,0,0,0,0]]


def support_from_powers(k, sums, p):
    """Recover a monic support polynomial via Newton's identities."""
    coeff = [1]
    for j in range(1,k+1):
        accum = sum(coeff[j-i]*sums[i] for i in range(1,j+1))
        coeff.append(-accum*pow(j,-1,p) % p)
    return roots(nmod_poly(coeff[::-1],p))


def consistent_moments(reduced, rank, k, p):
    """Resolve the one free parameter with the first unused Newton identity."""
    if rank == 5:
        return [[int(reduced[j,5]) for j in range(5)]]
    if rank != 4 or k == 4:
        return None
    pivot = [next(j for j in range(5) if reduced[i,j] != 0) for i in range(rank)]
    free, = set(range(5)) - set(pivot)
    sums = [None]*5
    sums[free] = nmod_poly([0,1],p)
    for i,j in enumerate(pivot):
        sums[j] = nmod_poly([int(reduced[i,5]),-int(reduced[i,free])],p)
    coeff = [1]
    for j in range(1,k+1):
        coeff.append(-sum(coeff[j-i]*sums[i] for i in range(1,j+1))*pow(j,-1,p))
    residual = sum(coeff[k+1-i]*sums[i] for i in range(1,k+2))
    if residual.is_zero():
        return None
    return [[int(poly(t)) for poly in sums] for t in roots(residual)]


def decode(shape, coefficients, odd, target_x, a, c, p):
    """Return True, False, or None. Inputs contain no masks, hashes or labels."""
    if p <= 5 or shape[:2] != (1,1) or not 1 <= shape[2] <= 4:
        return None
    k = shape[2]
    scale = pow(8,-1,p)
    q = [int(val)*scale % p for val in coefficients]
    g = int(odd)*scale % p
    W = nmod_poly(q,p) - 2*g*nmod_poly([c,a,0,1],p)
    if W.is_zero():
        return None

    def valid(x):
        return x != target_x and pow(curve_value(x,a,c,p),(p-1)//2,p) == 1

    answers = set()
    for u,v in combinations([x for x in roots(W) if valid(x)],2):
        equations = kernel(u,v,a,c,p)
        rhs = q + [g,k]
        augmented = nmod_mat([list(row)+[entry] for row,entry in zip(equations,rhs)],p)
        reduced, rank = augmented.rref()
        if any(all(reduced[i,j] == 0 for j in range(5)) and reduced[i,5] != 0
               for i in range(rank)):
            continue
        candidates = consistent_moments(reduced,rank,k,p)
        if candidates is None:
            return None
        for sums in candidates:
            support = support_from_powers(k,sums,p)
            if (len(support) != k or any(x in (u,v) or not valid(x) for x in support)
                    or sums != [sum(pow(x,j,p) for x in support) % p for j in range(5)]):
                continue
            left = pair_sum(u,v,a,c,p)
            hit = any(left.gcd(pair_sum(w,target_x,a,c,p)).degree() > 0 for w in support)
            answers.add(hit)
    return answers.pop() if len(answers) == 1 else None


def direct_feature(u, v, support, a, c, p):
    """Independent reference using the four Vandermonde cofactors."""
    X = nmod_poly([0,1],p)
    q = nmod_poly([],p)
    g = 0
    for w in support:
        xs = [u,v,w,X]
        for i,t in enumerate(xs):
            others = xs[:i]+xs[i+1:]
            det = (others[1]-others[0])*(others[2]-others[0])*(others[2]-others[1])
            q += (t**3+a*t+c)*det**2
            if i == 3:
                g += det**2
    return tuple(int(8*q[j]) % p for j in range(5)), 8*g % p


def exhaustive(p, a, width, target_index=1):
    curve = curve_cycle(p,a)
    target_x = curve[target_index][0]
    masks, weights = hash_families(curve,1<<p.bit_length(),width)
    support = [sorted({curve[i][0] for i in mask}) for mask in masks]
    gr = np.array([gram_data([curve[i] for i in mask],p) for mask in masks])
    eligible = []
    for ids in combinations(range(len(masks)),3):
        i,j,k = sorted(ids,key=lambda i:len(masks[i]))
        U,V,S = (support[t] for t in (i,j,k))
        if (len(U),len(V)) != (1,1) or not 1 <= len(S) <= 4:
            continue
        if len(set(U+V+S+[target_x])) != 3+len(S):
            continue
        eligible.append((i,j,k))
    idx = np.array(eligible)
    coeff,odd = moment_polynomial([gr[col] for col in idx.T],p,a,include_odd_block=True)
    stats = Counter()
    for (i,j,k),q,g in zip(eligible,coeff,odd):
        U,V,S = (support[t] for t in (i,j,k))
        assert (tuple(q),g) == direct_feature(U[0],V[0],S,a,1,p)
        truth = any((target_index-r-s) % len(curve) in masks[k] for r in masks[i] for s in masks[j])
        decision = decode((1,1,len(S)),q,g,target_x,a,1,p)
        assert decision is None or decision == truth
        mass = 6*weights[i]*weights[j]*weights[k]
        stats[len(S),"all"] += mass
        stats[len(S),"certified"] += mass*(decision is not None)
        stats[len(S),"positive"] += mass*truth
        stats[len(S),"found"] += mass*(decision is True)
    print(f"Exact F_{p} width {width}, target {curve[target_index]} (ordered seed multiplicities):")
    for k in sorted({k for k,label in stats}):
        print(f"  k={k}: "+", ".join(f"{s} {stats[k,s]:,}" for s in ("all","certified","positive","found")))


def fresh(p, a, width, seed):
    def audit(sizes, coefficients, odd, rows):
        tally = Counter()
        cached = lru_cache(maxsize=None)(decode)
        for sizes3, q, g, (eligible,unseen,truth,memorized) in zip(sizes,coefficients,odd,rows):
            shape = tuple(sorted(map(int,sizes3)))
            decision = cached(shape,tuple(map(int,q)),int(g),0,a,1,p) if eligible else None
            assert decision is None or decision == truth
            hybrid = memorized if decision is None else decision
            tally["standalone_errors"] += (decision is True) != truth
            tally["hybrid_errors"] += hybrid != truth
            tally["certified"] += decision is not None
            if unseen:
                tally["novel_certified"] += decision is not None
                tally["novel_true_positives"] += decision is True
                tally["novel_false_positives"] += decision is True and not truth
                if shape[:2] == (1,1) and 1 <= shape[2] <= 4:
                    tally["novel_in_shape"] += 1
                    tally["novel_in_shape_hits"] += truth
        print("  Structural readout: "+", ".join(f"{s} {v:,}" for s,v in tally.items()))
    assert curve_cycle(p,a)[1][0] == 0
    trial(p,a,width,seed,audit=audit)


def collision():
    """Two seeds with identical entire unweighted slot Gram data and opposite labels."""
    from fixtures import determinant
    p,a,c = 37,27,2
    target = (13,16)
    points = [(x,y) for x,y in product(range(p),repeat=2) if y*y % p == curve_value(x,a,c,p)]
    assert len(points)+1 == 47
    def slot(alpha,beta,width):
        return [P for P in points if (alpha*P[0]+beta) % 64 < width]
    def add(P,R):
        if P is None or R is None:
            return R if P is None else P
        x,y = P
        z,w = R
        if x == z and (y+w) % p == 0:
            return None
        num,den = (3*x*x+a,2*y) if P == R else (w-y,z-x)
        slope = num*pow(den,-1,p) % p
        new_x = (slope**2-x-z) % p
        return new_x,(slope*(x-new_x)-y) % p
    first,second = slot(49,37,8),slot(45,21,8)
    thirds = [slot(23,9,32),slot(11,24,32)]
    assert [sorted({x for x,y in S}) for S in [first,second]+thirds] == [
        [2,32],[1,21],[3,6,8,9,14,17,22,23,36],[4,5,6,16,17,18,22,23,27]]
    axis = [{x for x,y in S} for S in [first,second]+thirds]
    assert all(target[0] not in s for s in axis)
    assert all(axis[i].isdisjoint(axis[j]) for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3)])
    assert gram_data(thirds[0],p) == gram_data(thirds[1],p)
    counts = []
    weighted = []
    for S in thirds:
        triples = list(product(first,second,S))
        count = sum(add(add(P,R),U) == target for P,R,U in triples)
        zero_count = sum(determinant([[1,x,x*x,y] for x,y in (P,R,U,(target[0],-target[1]))],p)==0
                         for P,R,U in triples)
        assert count == zero_count
        counts.append(count)
        moments = [np.array([gram_data(F,p)]) for F in (first,second)]
        powers = [sum(pow(x,i+1,p) for x,y in S) % p for i in range(5)]
        data = powers+[sum(x*y*y for x,y in S) % p]
        q,odd = moment_polynomial(moments+[np.array([data])],p,a,c,include_odd_block=True)
        weighted.append((tuple(map(int,q[0])),int(odd[0])))
    assert counts == [0,5] and weighted[0] != weighted[1]
    print(f"Collision F_37: shared third-slot data {gram_data(thirds[0],p)}, counts {counts}, "
          f"x-weighted readouts {weighted}.")
    singles = [slot(1,6,8),slot(7,54,8)]
    assert [sorted({x for x,y in F}) for F in singles] == [[1],[2]]
    other_target = (30,5)
    assert all(other_target[0] != x for F in singles+thirds for x,y in F)
    single_counts = [sum(add(add(P,R),U) == other_target for P,R,U in product(*singles,S))
                     for S in thirds]
    assert single_counts == [1,0]
    print(f"Same third-slot collision at {other_target}, single-abscissa anchors: {single_counts}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    collision()
    for params in ((13,4,2),(23,5,2),(23,5,4,2)):
        exhaustive(*params)
    # Freeze the additional curve by selecting the first non-anomalous
    # prime-order curve y^2=x^3+a*x+1 with a>=1 at the next prime after 101.
    for params in ((23,5,2,14372),(101,17,4,94275),(103,33,4,103041),(103,33,8,103081)):
        fresh(*params)
