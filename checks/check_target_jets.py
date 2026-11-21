"""Exact fresh-seed feature-collision counts for target derivatives on toy curves."""

import numpy as np

from fixtures import curve_cycle, hash_families
from probe_target_jets import gram_data, moment_polynomial, taylor


def compress(array, indicators):
    for _ in range(3):
        array = np.tensordot(array, indicators.T, axes=(0, 0))
    return array


def study(p, a, width):
    curve = curve_cycle(p, a)
    points = curve[1:]
    xt, yt = points[0]
    order = len(curve)
    modulus = 1 << p.bit_length()
    masks, weights = hash_families(curve, modulus, width)
    indicators = np.array([[j in mask for j in range(1,order)] for mask in masks], dtype=np.int64)
    x, y = np.array(points).T
    rows = [np.ones_like(x), x, x*x % p, y]
    cf = []
    for omitted in range(4):
        u,v,w = [row for j,row in enumerate(rows) if j != omitted]
        cross = (np.outer(v,w)-np.outer(w,v)) % p
        minor = ((u[:,None,None]*cross[None,:,:]) - u[None,:,None]*cross[:,None,:]
                 + u[None,None,:]*cross[:,:,None])
        cf.append(((-1)**(3+omitted)*minor) % p)
    const, linear, quad, ycoef = cf
    # Expand in s = x - x(T). The target row at -T contains -y(T).
    coeff = [(const+linear*xt+quad*xt**2) % p, (linear+2*xt*quad) % p, quad]
    f = [yt**2 % p, (3*xt**2+a) % p, 3*xt % p, 1, 0]
    jets = []
    for n in range(5):
        power = sum(coeff[j]*coeff[n-j] for j in range(3) if 0 <= n-j < 3)
        jets.append(compress((power + ycoef*ycoef*f[n]) % p,indicators) % p)
    delta = (coeff[0]-ycoef*yt) % p
    raw = []
    power = np.ones_like(delta)
    for _ in range(5):
        power = power*delta % p
        raw.append(compress(power,indicators) % p)
    assert np.array_equal(raw[1],jets[0])
    # Independent seed-level computation via the Gram determinant, with no
    # determinant tuples serving as input. Compare every mask triple in slabs.
    gr = np.array([gram_data([curve[j] for j in mask],p) for mask in masks])
    jmask, kmask = np.indices((len(masks),len(masks))).reshape(2,-1)
    for i in range(len(masks)):
        slab = moment_polynomial([np.broadcast_to(gr[i],(len(jmask),6)),gr[jmask],gr[kmask]],p,a)
        assert np.array_equal(taylor(slab,xt,p),np.column_stack([jet[i].ravel() for jet in jets]))
    # Validate all five Taylor coefficients by direct evaluation on many targets.
    for xc,yc in points:
        predicted = sum(t*pow(xc-xt,i,p) for i,t in enumerate(jets)) % p
        at = (const+linear*xc+quad*xc**2-ycoef*yc) % p
        assert np.array_equal(predicted,compress(at*at % p,indicators) % p)
    labels = np.arange(1,order)
    actual = (labels[:,None,None]+labels[None,:,None]+labels[None,None,:]) % order == 1
    counts = compress(actual,indicators)
    sizes = indicators.sum(axis=1)
    disjoint = (indicators@indicators.T)==0
    allowed = (sizes>0) & (indicators[:,0]==0)
    eligible = (allowed[:,None,None] & allowed[None,:,None] & allowed[None,None,:]
                & disjoint[:,:,None] & disjoint[:,None,:] & disjoint[None,:,:])
    assert np.array_equal(counts[eligible],compress(delta==0,indicators)[eligible])
    idx = np.nonzero(eligible)
    sizes3 = np.column_stack([sizes[i] for i in idx])
    weights3 = np.prod(np.column_stack([np.array(weights)[i] for i in idx]),axis=1)
    truth = counts[idx]>0
    pool = np.prod(sizes3,axis=1)//8
    print(f"F_{p}, order {order}, widths {width}, target {points[0]}")
    for name, subset in (("regular nonempty",pool>=1),("at least three sign orbits",pool>=3)):
        w,h = weights3[subset], truth[subset]
        ids=tuple(i[subset] for i in idx)
        print(f"  {name}: mass {w.sum():,}, positive {w[h].sum():,}")
        features = {"Z2":[raw[1]], "raw2":raw[:2], "jet2":jets[:2],
                    "raw3":raw[:3], "raw2+dZ2":[*raw[:2],jets[1]],
                    "raw5":raw, "jet5":jets, "raw5+jet":raw+jets[1:]}
        for label, feature in features.items():
            data = np.column_stack([sizes3[subset],*(v[ids] for v in feature)])
            _, inv = np.unique(data,axis=0,return_inverse=True)
            masses = np.zeros((inv.max()+1,2),dtype=np.int64)
            np.add.at(masses,(inv,h.astype(int)),w)
            print(f"    {label:10s} {masses.min(axis=1).sum():,}")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    for params in ((13,4,2),(13,4,4),(23,5,2)):
        study(*params)
