"""Plant-only field-equation learning with independent fresh tests on toy curves."""

from collections import Counter
from itertools import product
from math import comb, log2, prod
from random import Random

import numpy as np
from flint import nmod_mat
from scipy.stats import beta

from pooled_zero import pooled_certificate

from fixtures import curve_cycle


def monomials(moments, size, p):
    if size % p == 0 or moments[0] == 0:
        return None
    inv_size = pow(size, -1, p)
    center = int(moments[0])*inv_size % p
    raw = [1] + [int(v)*inv_size * pow(center, -j, p) % p for j, v in enumerate(moments[:6], 1)]
    central = [sum(comb(j, i)*(-1)**(j-i)*raw[i] for i in range(j+1)) % p for j in (2,4,6)]
    # All monomials in the scale-normalized even central moments of weight at most eight.
    return [prod(pow(c, e, p) for c, e in zip(central, exps)) % p
            for exps in MONOMIALS]


def raw_quadrics(sums, p):
    return [1, *sums, *(x*y % p for i, x in enumerate(sums) for y in sums[i:])]


MONOMIALS = [a for a in product(range(5), range(3), range(2))
             if 2*a[0]+4*a[1]+6*a[2] <= 8]


class Fixture:
    def __init__(self, p, a, widths):
        self.curve = curve_cycle(p, a)
        self.p, self.widths = p, widths
        self.order = len(self.curve)
        self.modulus = 1 << p.bit_length()
        self.target = 1
        xt, yt = self.curve[self.target]
        self.x = [None] + [v[0] for v in self.curve[1:]]
        self.fibers = [[] for _ in range(self.modulus)]
        for j in range(1, self.order):
            self.fibers[self.x[j]].append(j)
        self.rows = np.array([[0,0,0]] + [[(x-xt) % p, (x*x-xt*xt) % p, (y+yt) % p]
                                            for x, y in self.curve[1:]], dtype=np.int64)
        self.inverse = {a: pow(a, -1, self.modulus) for a in range(1, self.modulus, 2)}
        self.lam = self.order**2 * prod(widths)/self.modulus**3
        assert 1/(2*log2(p.bit_length())) < self.lam <= 1/log2(p.bit_length())

    def observe(self, rng, mode, truth=False):
        if mode == "plant":
            while True:
                i, j = rng.randrange(1, self.order), rng.randrange(1, self.order)
                k = (self.target-i-j) % self.order
                if k:
                    anchors = (i,j,k)
                    break
        elif mode == "matched":
            anchors = tuple(rng.randrange(1, self.order) for _ in self.widths)
        else:
            anchors = None
        occupied = {self.x[self.target]}
        slots, regular = [], True
        for i, width in enumerate(self.widths):
            alpha = rng.randrange(1, self.modulus, 2)
            beta = (rng.randrange(self.modulus) if anchors is None
                    else (rng.randrange(width)-alpha*self.x[anchors[i]]) % self.modulus)
            ids = []
            for u in range(width):
                x = (self.inverse[alpha]*(u-beta)) % self.modulus
                if self.fibers[x]:
                    regular &= x not in occupied
                    occupied.add(x)
                    ids.extend(self.fibers[x])
            slots.append(np.array(ids, dtype=int))
        size = prod(map(len, slots))
        hits = None
        if truth:
            left, mid, right = slots
            hits = int(np.sum((left[:,None,None] + mid[None,:,None] + right[None,None,:])
                              % self.order == self.target))
        if not regular or not size:
            return None, hits
        left, mid, right = (self.rows[s] for s in slots)
        delta = (-np.cross(left[:,None,:], mid[None,:,:]).reshape(-1,3) @ right.T).ravel() % self.p
        if truth:
            assert int(np.sum(delta == 0)) == hits
        pool = size//8
        sums, power = [], delta.copy()
        for _ in range(9):
            sums.append(int(power.sum() % self.p))
            power = power*delta % self.p
        return (pool, size, sums), hits


def learn(fixture, count, rng, mode):
    collections = {s: [] for s in ("one", "two", "many", "many raw")}
    for _ in range(count):
        obs, _ = fixture.observe(rng, mode)
        if obs is None:
            continue
        pool, size, sums = obs
        row = monomials(sums, size, fixture.p)
        if row is not None:
            collections["one" if pool == 1 else "two" if pool == 2 else "many"].append(row)
        if pool > 2:
            collections["many raw"].append(raw_quadrics(sums, fixture.p))
    models = {}
    for s, data in collections.items():
        if not data:
            continue
        mat = nmod_mat(data, fixture.p)
        kernel, dimension = mat.nullspace()
        d = len(data[0])
        models[s] = [[int(kernel[i,j]) for i in range(d)] for j in range(dimension)]
        print(f"  {mode}, {s}: {len(data):,} samples, feature rank {mat.rank()}/{d}.")
        if mode == "plant" and s in ("many", "many raw"):
            # If any hyperplane contains a fraction h of the planted law, a
            # block of d independent rows is singular with probability >= h**d.
            blocks = len(data)//d
            if blocks:
                singular = sum(nmod_mat(data[i*d:(i+1)*d], fixture.p).rank() < d
                               for i in range(blocks))
                upper = (beta.ppf(.995, singular+1, blocks-singular)
                         if singular < blocks else 1.)
                print(f"  Singular plant blocks {singular}/{blocks}, simultaneous "
                      f"upper bound on one-equation recall {upper**(1/d):.4f} (99.5%).")
    return models


def run(p, a, widths, seed, n_train=12000, n_test=16000):
    fixture = Fixture(p, a, widths)
    print(f"F_{p}, order {fixture.order}, widths {widths}, lambda={fixture.lam:.4f}")
    sequences = np.random.SeedSequence(seed).generate_state(3)
    plant_rng, matched_rng, test_rng = [Random(int(i)) for i in sequences]
    model = learn(fixture, n_train, plant_rng, "plant")
    learn(fixture, n_train, matched_rng, "matched")
    tot = Counter()
    cache = {}
    for _ in range(n_test):
        obs, hits = fixture.observe(test_rng, "fresh", truth=True)
        present = hits > 0
        tot["hits"] += present
        answer = False
        if obs:
            pool, size, sums = obs
            tot["regular nonempty"] += 1
            name = "one" if pool == 1 else "two" if pool == 2 else "many"
            tot[name] += 1
            tot[name, "hits"] += present
            row = monomials(sums, size, p)
            kernel = model.get(name)
            if row is not None and kernel:
                answer = all(sum(a*b for a, b in zip(row, basis)) % p == 0 for basis in kernel)
            if name == "many" and model.get("many raw"):
                raw = raw_quadrics(sums, p)
                answer |= all(sum(a*b for a, b in zip(raw, basis)) % p == 0
                              for basis in model["many raw"])
            if pool == 2:
                key = tuple(sums)
                if key not in cache:
                    cache[key] = pooled_certificate(p, sums, orders=(9,))[9]
                possible = cache[key]
                if possible is not None:
                    assert hits in possible
                    if len(possible) == 1:
                        tot["exact pooled counts"] += 1
                    if possible == frozenset({0}) or 0 not in possible:
                        answer = 0 not in possible
            tot[name,"errors"] += answer != present
        elif present:
            # A hit requires a nonempty slot product, so this mask is irregular.
            tot["irregular hits"] += 1
        tot["errors"] += answer != present
        tot["false positives"] += answer and not present
    print(f"  {n_test:,} independent fresh seeds, {tot['hits']:,} positives, "
          f"{tot['regular nonempty']:,} regular nonempty.")
    for name in ("one","two","many"):
        print(f"  {name}: {tot[name]:,} seeds, {tot[name,'hits']:,} hits, {tot[name,'errors']:,} errors.")
    assert tot["hits"] == tot["irregular hits"] + sum(tot[name, "hits"] for name in ("one", "two", "many"))
    print(f"  Irregular masks: {tot['irregular hits']:,} fresh hits.")
    print(f"  Nine-moment exact counts on {tot['exact pooled counts']:,}/{tot['two']:,} two-orbit seeds.")
    print(f"  Total readout errors {tot['errors']:,}/{n_test:,}, including {tot['false positives']:,} false positives. "
          f"Always absent makes {tot['hits']:,} errors.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run(101, 17, (4,4,4), 94201)
    run(1031, 9, (8,16,16), 94202)
