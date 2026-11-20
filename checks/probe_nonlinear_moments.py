"""Nonlinear signed-moment probes against wrong legal addresses on the same seed."""
from itertools import combinations, product
from math import prod, sqrt
from random import Random

import numpy as np

from check_strided_isolation import finite_laws
from fixtures import curve_cycle, evaluation
from probe_signed_moment import elimination


def walsh(values):
    """Unnormalized Walsh transforms along the last axis, in integer arithmetic."""
    out = np.array(values, dtype=np.int64, copy=True)
    dimension = out.shape[-1]
    half = 1
    while half < dimension:
        sections = out.reshape(-1, dimension//(2*half), 2, half)
        left = sections[:, :, 0].copy()
        right = sections[:, :, 1].copy()
        sections[:, :, 0], sections[:, :, 1] = left+right, left-right
        half *= 2
    return out


class FieldProbe:
    def __init__(self, p, linear, count, positions, target=1):
        self.p, self.count, self.positions, self.target = p, count, positions, target
        self.curve = curve_cycle(p, linear)
        self.order = len(self.curve)
        self.modulus = 1 << p.bit_length()
        self.width = 1 << len(positions)
        self.block = 2*self.width
        self.size = self.block**count
        sign_mask = sum(self.block**j for j in range(count))
        self.abscissa_only = np.arange(self.size) & sign_mask == 0
        # Subtract the endpoint row, then expand the constant column.
        endpoint = np.array(evaluation(self.curve[-target], count+1, p)[1:])
        self.eta = [None]+[(np.array(evaluation(v, count+1, p)[1:])-endpoint) % p
                           for v in self.curve[1:]]
        self.xs = {}
        for v, (x, y) in enumerate(self.curve[1:], 1):
            self.xs.setdefault(x, []).append((v, int(y > (-y) % p)))
        self.chi = np.full(p, -1, dtype=np.int64)
        self.chi[0] = 0
        self.chi[[a*a % p for a in range(1, p)]] = 1
        self.checks = 0

    def plant(self, rng):
        while True:
            points = [rng.randrange(1, self.order) for _ in range(self.count-1)]
            last = (self.target-sum(points)) % self.order
            if last:
                points.append(last)
                break
        alpha = 1+2*rng.randrange(self.modulus//2)
        beta = (rng.randrange(self.width)-alpha*self.curve[points[0]][0]) % self.modulus
        inv = pow(alpha, -1, self.modulus)
        images = [[inv*(i-beta) % self.modulus for i in range(self.width)]]
        free = sum(1 << j for j in self.positions)
        for v in points[1:]:
            base = self.curve[v][0] & ~free
            images.append([base+sum(((i >> j) & 1) << bit for j, bit in enumerate(self.positions))
                           for i in range(self.width)])
        return [tuple((v, 2*i+sign) for i, x in enumerate(xs) for v, sign in self.xs.get(x, ()))
                for xs in images]

    def witness_and_domain(self, slots):
        domain, witness, hits = [], None, 0
        for choices in product(*slots):
            address, partial, guarded = 0, 0, True
            for j, (v, digit) in enumerate(choices):
                address = address*self.block+digit
                if j and (not partial or self.curve[partial][0] == self.curve[v][0]):
                    guarded = False
                partial = (partial+v) % self.order
            domain.append(address)
            if partial == self.target:
                hits += 1
                if guarded:
                    witness = address
        if hits != 1 or witness is None:
            return None
        return witness, np.array(domain)

    def grids(self, slots):
        # This enumerative reference constructs exact values on the entire binary
        # grid. Only the independent per-challenge formula has the advertised cost.
        minors = {(): np.ones(1, dtype=np.int64)}
        lane = np.zeros(1, dtype=np.int64)
        for j, slot in enumerate(slots, 1):
            rows = np.zeros((self.block, self.count), dtype=np.int64)
            occupied = np.zeros(self.block, dtype=np.int64)
            for point, addr in slot:
                rows[addr] = self.eta[point]
                occupied[addr] = 1
            lane = np.add.outer(lane, occupied).ravel()
            next_minors = {}
            for cols in combinations(range(self.count), j):
                terms = sum((-1)**(j-1-k)*np.multiply.outer(minors[cols[:k]+cols[k+1:]], rows[:, col])
                            for k, col in enumerate(cols))
                next_minors[cols] = terms.ravel() % self.p
            minors = next_minors
        delta = (-1)**self.count*next(iter(minors.values()))
        return np.stack([lane == self.count, delta, (delta**2) % self.p])

    def check_factorization(self, slots, r, moments):
        means, grams = [], []
        mass = 1
        for i, slot in enumerate(slots):
            digit = r//self.block**(self.count-i-1) % self.block
            weights = np.array([(-1)**((digit & addr).bit_count()) for _, addr in slot])
            rows = np.array([self.eta[v] for v, _ in slot])
            means.append(weights @ rows)
            grams.append(rows.T @ (weights[:, None]*rows))
            mass *= sum(int(w) for w in weights)
        total = 0
        for k in range(1, self.count+1):
            for block in combinations(grams, k):
                g = np.sum(block, axis=0) % self.p
                total += (-1)**(self.count-k)*elimination(g.tolist(), self.p)
        first = (-1)**self.count*elimination([row.tolist() for row in means], self.p)
        assert tuple(moments[:, r]) == (mass % self.p, first % self.p, total % self.p)
        self.checks += 1

    def labelled_moments(self, slots):
        hit = self.witness_and_domain(slots)
        if hit is None:
            return None
        witness, legal = hit
        grid = self.grids(slots)
        assert grid[1, witness] % self.p == 0
        assert np.count_nonzero(grid[0]) == len(legal) and np.all(grid[0, legal])
        moments = walsh(grid) % self.p
        if self.checks < 16:
            self.check_factorization(slots, self.checks*(self.size-1)//15, moments)
        clean = np.count_nonzero((grid[1] % self.p == 0) & (grid[0] == 1)) == 1
        return moments, witness, legal[legal != witness], clean

    def contrast(self, witness, wrong):
        # Every wrong legal address is equally weighted on this seed.
        marks = np.zeros(self.size, dtype=np.int64)
        marks[witness], marks[wrong] = len(wrong), -1
        return walsh(marks)/len(wrong)

    def root(self, slots):
        labelled = self.labelled_moments(slots)
        if labelled is None:
            return None
        moments, witness, wrong, clean = labelled
        mass, first, second = moments
        disc = (mass*second-first**2) % self.p
        mix = (second+first**2) % self.p
        fixed = np.stack((self.chi[first], self.chi[second], self.chi[mix], self.chi[disc]))
        if self.count >= 4 and self.order % 2:
            assert np.all(first[self.abscissa_only] == 0) and np.all(mass[~self.abscissa_only] == 0)
            reduced = self.chi[-1]*(first != 0)+self.abscissa_only*self.chi[mass]*self.chi[second]
            assert np.array_equal(fixed[3], reduced)
        cells = sum((fixed[i]+1)*3**i for i in range(4))
        return fixed, cells, self.contrast(witness, wrong), witness, wrong, clean

    def sample(self, rng, count):
        for _ in range(count):
            while True:
                root = self.root(self.plant(rng))
                if root is not None:
                    yield root
                    break


def learn(probe, seed, n):
    total, covariance = np.zeros(81), np.zeros((81, 81))
    choices = Random(seed+314159)
    parity = np.array([(-1)**a.bit_count() for a in range(probe.size)])
    challenges = np.arange(probe.size)
    for _, cells, contrast, witness, wrong, _ in probe.sample(Random(seed), n):
        total += np.bincount(cells, weights=contrast, minlength=81)
        # A mean of eight uniformly selected wrongs is unbiased for the
        # same-seed second-moment control. Evaluation later uses every wrong.
        controls = choices.sample(wrong.tolist(), min(8, len(wrong)))
        hist = np.array([np.bincount(cells, weights=parity[challenges & v], minlength=81)
                         for v in [witness, *controls]])
        covariance += np.outer(hist[0], hist[0])-hist[1:].T @ hist[1:]/len(controls)
    simple = np.sign(np.round(total, 10)).astype(np.int64)
    # With binary cell weights, the difference of squared coefficients is a
    # quadratic form. Use deterministic starting points and bit flips.
    best_value = float("-inf")
    for initial in range(8):
        weights = np.array([1 if initial == 0 or choices.randrange(2) else -1 for _ in range(81)])
        for _ in range(100):
            changed = False
            for cell in range(81):
                gain = -4*weights[cell]*(covariance[cell] @ weights-covariance[cell, cell]*weights[cell])
                if gain > 1e-8*n:
                    weights[cell] *= -1
                    changed = True
            if not changed:
                break
        value = weights @ covariance @ weights
        if value > best_value:
            best_value, best = value, weights.copy()
    if best @ total < 0:
        best = -best
    return simple, best


class Results:
    def __init__(self):
        self.mass = 0
        self.sums = np.zeros((6, 5))
        self.gap2 = np.zeros(6)

    def add(self, scores, witness, wrong, weight):
        self.mass += weight
        signed = walsh(scores)/scores.shape[1]
        spectrum = np.abs(signed)
        true = spectrum[:, witness]
        others = spectrum[:, wrong]
        contrast = true-others.mean(axis=1)
        rank = ((true[:, None] > others)+0.5*(true[:, None] == others)).mean(axis=1)
        oriented = signed[:, witness]-signed[:, wrong].mean(axis=1)
        self.sums += weight*np.stack((true, others.mean(axis=1), contrast, rank, oriented), axis=1)
        self.gap2 += weight*contrast**2

    def display(self, cohort, exact):
        print(f"  {cohort}: {self.mass} singleton seeds")
        for i, name in enumerate(("chi(M1)", "chi(M2)", "chi(M2+M1^2)", "chi(M0*M2-M1^2)",
                                  "fit81-contrast", "fit81-power")):
            true, other, gap, rank, oriented = self.sums[i]/self.mass
            se = sqrt(max(0., self.gap2[i]/self.mass-gap**2)/(self.mass-1))
            error = "" if exact else f", paired SE {se:.5f}"
            print(f"    {name}: witness {true:.5f}, wrong {other:.5f}, gap {gap:+.5f}{error}, "
                  f"rank {rank:.4f}, signed gap {oriented:+.5f}")


def run(p, linear, count, positions, train, heldout):
    probe = FieldProbe(p, linear, count, positions)
    fit_signed, fit_power = learn(probe, 7471+p+count, train)
    tiny = p == 13
    if tiny:
        affine, boxes = finite_laws(p, probe.curve, positions)
        laws = [affine]+[boxes]*(count-1)
        roots = ((probe.root(slots), prod(law[slot] for law, slot in zip(laws, slots))) for slots in product(*laws))
    else:
        roots = ((root, 1) for root in probe.sample(Random(97417+p+count), heldout))
    all_seeds, clean_seeds = Results(), Results()
    for root, weight in roots:
        if root is None:
            continue
        fixed, cells, _, witness, wrong, clean = root
        scores = np.vstack((fixed, fit_signed[cells], fit_power[cells]))
        all_seeds.add(scores, witness, wrong, weight)
        if clean:
            clean_seeds.add(scores, witness, wrong, weight)
    assert probe.checks == 16
    mode = "exact conditional law" if tiny else "independent conditional draws"
    print(f"p={p}, m={count}, w={probe.width}, {mode}, training {train}; {probe.checks} joint moment checks")
    all_seeds.display("all", tiny)
    clean_seeds.display("one determinant zero", tiny)


def confirm(target, seed, count):
    """Independently repeat the largest-field mixture following the initial comparison."""
    probe = FieldProbe(65521, 44, 6, (0,), target)
    binary = np.arange(probe.size)
    abscissa = sum((((binary//probe.block**j) % probe.block)//2)*probe.width**j
                   for j in range(probe.count))
    values, projected = [], []
    for fixed, _, _, witness, wrong, _ in probe.sample(Random(seed), count):
        spectrum = walsh(fixed[2])/probe.size
        true = abs(spectrum[witness])
        control = np.abs(spectrum[wrong])
        rank = np.mean((true > control)+0.5*(true == control))
        values.append((true-control.mean(), rank, spectrum[witness]-spectrum[wrong].mean()))
        # On the all-even slice of the variance score, M1 vanishes and the
        # remaining function is the abscissa-only score chi(M0*M2).
        candidate = abscissa[witness]
        rivals = np.unique(abscissa[wrong])
        rivals = rivals[rivals != candidate]
        if rivals.size:
            coefficients = walsh(fixed[3, probe.abscissa_only])
            assert np.all(coefficients[candidate] == coefficients[rivals])
            orbit = np.abs(coefficients)/probe.width**probe.count
            value, others = orbit[candidate], orbit[rivals]
            rank = np.mean((value > others)+0.5*(value == others))
            projected.append((value, others.mean(), value-others.mean(), rank))
    mean = np.mean(values, axis=0)
    se = np.std(values, axis=0, ddof=1)/sqrt(count)
    print(f"p=65521, m=6, w=2, independent repeat T=[{target}]G, n={count}, chi(M2+M1^2): "
          f"magnitude gap {mean[0]:+.5f} (SE {se[0]:.5f}), rank {mean[1]:.4f} (SE {se[1]:.4f}), "
          f"signed gap {mean[2]:+.5f} (SE {se[2]:.5f})")
    mean = np.mean(projected, axis=0)
    se = np.std(projected, axis=0, ddof=1)/sqrt(len(projected))
    print(f"  chi(M0*M2), abscissa-only T=[{target}]G, n={len(projected)} with alternatives: "
          f"witness {mean[0]:.5f}, other {mean[1]:.5f}, gap {mean[2]:+.5f} (SE {se[2]:.5f}), "
          f"rank {mean[3]:.4f} (SE {se[3]:.4f})")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run(13, 4, 3, (0,), 1000, 0)
    run(101, 17, 3, (0,), 1000, 1200)
    run(1031, 9, 4, (0, 3), 800, 1200)
    run(1031, 9, 5, (0,), 800, 1200)
    run(65521, 44, 6, (0,), 800, 1200)
    confirm(1, 816905, 5000)
    confirm(2, 315480, 5000)
