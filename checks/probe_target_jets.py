"""A slotwise evaluator of the full second-moment polynomial and a planted lookup trial."""
from collections import Counter
from itertools import combinations, product
from math import comb
from random import Random

import numpy as np

from fixtures import curve_cycle, hash_families


def gram_data(points, p):
    return [sum(pow(x,j,p) for x,y in points) % p for j in range(5)] + [
            sum(y*y for x,y in points) % p]


def moment_polynomial(slots, p, a, c=1, *, include_odd_block=False):
    """Increasing x-coefficients, optionally with the coefficient of y(T)^2."""
    coeff = np.zeros((len(slots[0]),5),dtype=np.int64)
    if include_odd_block:
        odd = np.zeros(len(slots[0]), dtype=np.int64)
    for n in range(1,4):
        for subset in combinations(slots,n):
            s0,s1,s2,s3,s4,yy = (sum(subset) % p).T
            a00 = (s2*s4-s3*s3) % p
            a01 = (s2*s3-s1*s4) % p
            a02 = (s1*s3-s2*s2) % p
            a11 = (s0*s4-s2*s2) % p
            a12 = (s1*s2-s0*s3) % p
            a22 = (s0*s2-s1*s1) % p
            det = (s0*a00+s1*a01+s2*a02) % p
            if include_odd_block:
                odd += (-1)**(3-n)*det
            even = np.column_stack([a00,2*a01,2*a02+a11,2*a12,a22])
            coeff += (-1)**(3-n)*(yy[:,None]*even+det[:,None]*np.array([c,a,0,1,0]))
    if include_odd_block:
        return coeff % p, odd % p
    return coeff % p


def taylor(coeff, x, p):
    return np.column_stack([sum(comb(j,s)*pow(x,j-s,p)*coeff[:,j] for j in range(s,5)) % p
                            for s in range(5)])


def trial(p, a, width, seed, *, audit=None):
    curve = curve_cycle(p,a)
    order, modulus = len(curve), 1 << p.bit_length()
    masks, _ = hash_families(curve,modulus,width)
    to_index = {mask:i for i,mask in enumerate(masks)}
    mask_of = {(u,v):to_index[tuple(j for j,P in enumerate(curve[1:],1)
                                     if (u*P[0]+v) % modulus < width)]
               for u,v in product(range(1,modulus,2),range(modulus))}
    gr = np.array([gram_data([curve[i] for i in mask],p) for mask in masks])
    valid = [bool(mask) and 1 not in mask for mask in masks]
    mask_sets = [set(mask) for mask in masks]
    mu_num, mu_den = ((order-1)**3+1)*width**3, order*modulus**3

    def draw(num, plant, rng):
        samples = []
        for _ in range(num):
            anchors = []
            if plant:
                while not anchors:
                    i,j = rng.randrange(1,order),rng.randrange(1,order)
                    k = (1-i-j) % order
                    if k:
                        anchors = [i,j,k]
            idx = []
            for t in range(3):
                u = rng.randrange(1,modulus,2)
                v = ((rng.randrange(width)-u*curve[anchors[t]][0]) % modulus
                     if plant else rng.randrange(modulus))
                idx.append(mask_of[u,v])
            samples.append(tuple(sorted(idx)))
        return samples

    def observations(samples):
        idx = np.array(samples)
        jets, odd = moment_polynomial([gr[col] for col in idx.T],p,a,include_odd_block=True)
        for key, jet, g in zip(samples,jets,odd):
            i,j,k = key
            if (not all(valid[t] for t in key) or not mask_sets[i].isdisjoint(mask_sets[j])
                    or not mask_sets[i].isdisjoint(mask_sets[k]) or not mask_sets[j].isdisjoint(mask_sets[k])):
                yield None
            else:
                shape = tuple(sorted(len(masks[t]) for t in key))
                yield (shape+tuple(jet), key, shape+tuple(jet)+(int(g),))

    num = 30000
    fits = []
    for plant, offset in [(False,0),(True,1)]:
        bins = [Counter(),Counter(),Counter()]
        for obs in observations(draw(num,plant,Random(seed+offset))):
            if obs:
                for bucket,key in zip(bins,obs):
                    bucket[key] += 1
        fits.append(bins)
    fresh, planted = fits
    positive = [{key for key,num_plants in pos.items()
                 if 2*mu_num*num_plants > mu_den*unlabeled[key]}
                for pos,unlabeled in zip(planted,fresh)]
    samples = draw(20000,False,Random(seed+2))
    error = [0,0,0]
    total_hits = novel_hits = novel = 0
    novel_correct_pos = novel_false_pos = 0
    augmented_novel_tp = augmented_novel_fp = 0
    if audit is not None:
        audit_rows = []
    for key,obs in zip(samples,observations(samples)):
        i,j,k = key
        hit = any((1-u-v) % order in mask_sets[k] for u,v in product(masks[i],masks[j]))
        total_hits += hit
        preds = [bool(obs and obs[c] in positive[c]) for c in range(3)]
        for c in range(3):
            error[c] += preds[c] != hit
        unseen = bool(obs and key not in planted[1] and key not in fresh[1])
        if audit is not None:
            audit_rows.append((obs is not None, unseen, hit, preds[1]))
        if unseen:
            novel += 1
            novel_hits += hit
            novel_correct_pos += preds[0] and hit
            novel_false_pos += preds[0] and not hit
            augmented_novel_tp += preds[2] and hit
            augmented_novel_fp += preds[2] and not hit
    print(f"F_{p}, {num:,} plants + {num:,} unlabeled seeds; 20,000 independent fresh tests:")
    print(f"  Errors: absent {total_hits:,}, jet {error[0]:,}, memorized mask triple {error[1]:,}.")
    print(f"  Eligible new mask triples: {novel:,} draws, {novel_hits:,} positives; "
          f"jet called {novel_correct_pos:,} true and {novel_false_pos:,} false positives on them.")
    if audit is not None:
        print(f"  Jet plus odd-block lookup: {error[2]:,} errors; unseen true/false positives "
              f"{augmented_novel_tp:,}/{augmented_novel_fp:,}.")
        indices = np.array(samples)
        coeff, odd = moment_polynomial([gr[col] for col in indices.T],p,a,include_odd_block=True)
        sizes = np.array([len(mask)//2 for mask in masks])[indices]
        audit(sizes, coeff, odd, audit_rows)


if __name__ == "__main__":
    trial(23,5,2,14372)
    trial(101,17,4,94275)
