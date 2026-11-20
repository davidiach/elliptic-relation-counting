"""Dense raw-residue nonlinear readouts of two signed determinant moments."""
from math import sqrt
from random import Random

import numpy as np

from probe_nonlinear_moments import FieldProbe, learn, walsh


def residue_bin(x, p):
    """A separate zero bin and eight contiguous bins partitioning the nonzero residues."""
    return np.where(x == 0, 0, 1+8*(x-1)//(p-1))


def centered_bit(x, p):
    return np.where(x == 0, 0, np.where(x <= p//2, 1, -1))


class DenseProbe(FieldProbe):
    def __init__(self, p, linear, count, positions, target):
        super().__init__(p, linear, count, positions, target)
        self.nonsquare = next(z for z in range(2, p) if self.chi[z] == -1)
        assert pow(self.nonsquare, (p-1)//2, p) == p-1
        # Fixture acceleration. A single actual query only needs one field inverse.
        inverses = [0, 1]
        for i in range(2, p):
            inverses.append(-(p//i)*inverses[p % i] % p)
        self.inverses = np.array(inverses)
        assert np.all((np.arange(1, p)*self.inverses[1:]) % p == 1)

    def root(self, slots):
        labelled = self.labelled_moments(slots)
        if labelled is None:
            return None
        (_, first, second), witness, wrong, clean = labelled
        cells = 9*residue_bin(first, self.p)+residue_bin(second, self.p)
        mixture = (second+first**2) % self.p
        denominator = (first**2-self.nonsquare) % self.p
        ratio = (second*self.inverses[denominator]) % self.p
        fixed = np.stack((centered_bit(mixture, self.p), centered_bit((mixture+self.p//4) % self.p, self.p),
                          centered_bit(ratio, self.p)))
        return fixed, cells, self.contrast(witness, wrong), witness, wrong, clean


def run(p, linear, m, positions, target, train, n):
    probe = DenseProbe(p, linear, m, positions, target)
    base = p+977*m+3841*target
    oriented, power = learn(probe, 93531+base, train)
    table = power.reshape(9, 9)
    # The alternating difference vanishes on either coordinate axis. It isolates
    # the genuinely two-variable part relative to the public zero baselines.
    interaction = (table-table[:, :1]-table[:1, :]+table[0, 0]).ravel()
    assert np.all(interaction.reshape(9, 9)[0] == 0)
    assert np.all(interaction.reshape(9, 9)[:, 0] == 0)
    scales = np.array([1, 1, 1, 1, 1, 4])[:, None]
    names = ("center(M2+M1^2)", "center(shift+M2+M1^2)", "center(M2/(M1^2-d))",
             "joint-contrast", "joint-power", "power interaction")
    results = []
    occupancy = np.zeros(81, dtype=np.int64)
    for fixed, cells, _, witness, wrong, _ in probe.sample(Random(892313+base), n):
        scores = np.vstack((fixed, oriented[cells], power[cells], interaction[cells]))
        signed = walsh(scores)/probe.size/scales
        real = np.abs(signed[:, witness])
        competitors = np.abs(signed[:, wrong])
        control = competitors.mean(axis=1)
        rank = ((real[:, None] > competitors)+0.5*(real[:, None] == competitors)).mean(axis=1)
        signed_gap = signed[:, witness]-signed[:, wrong].mean(axis=1)
        results.append(np.stack((real, control, real-control, rank, signed_gap), axis=1))
        occupancy += np.bincount(cells, minlength=81)
    assert probe.checks == 16 and occupancy.sum() == probe.size*n
    means = np.mean(results, axis=0)
    errors = np.std(results, axis=0, ddof=1)/sqrt(n)
    joint = occupancy.reshape(9, 9)[1:, 1:].sum()/occupancy.sum()
    print(f"p={p}, m={m}, w={probe.width}, T=[{target}]G, d={probe.nonsquare}, train={train}, fresh={n}, "
          f"occupied cells {np.count_nonzero(occupancy)}/81, both moments nonzero {joint:.4f}")
    for name, mean, error in zip(names, means, errors):
        real, control, gap, rank, signed_gap = mean
        print(f"  {name}: witness {real:.5f}, wrong {control:.5f}, gap {gap:+.5f} (SE {error[2]:.5f}), "
              f"rank {rank:.4f} (SE {error[3]:.4f}), signed gap {signed_gap:+.5f}")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run(1031, 9, 4, (0, 3), 1, 1200, 1600)
    run(1031, 9, 4, (0, 3), 2, 1200, 1600)
    run(65521, 44, 6, (0,), 1, 1200, 1600)
    run(65521, 44, 6, (0,), 2, 1200, 1600)
