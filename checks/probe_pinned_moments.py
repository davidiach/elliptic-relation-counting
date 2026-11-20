"""Compare nonlinear moment scores after projecting onto a single sign vector."""
from collections import defaultdict
from copy import copy
from itertools import product
from math import prod, sqrt
from random import Random

import numpy as np

from probe_nonlinear_moments import FieldProbe, walsh


def unique_relation(probe, slots):
    """Count all original signed solutions before choosing a sign vector."""
    split = len(slots)//2
    left = defaultdict(list)
    for choices in product(*slots[:split]):
        left[sum(v for v, _ in choices) % probe.order].append(choices)
    hit = None
    for choices in product(*slots[split:]):
        key = (probe.target-sum(v for v, _ in choices)) % probe.order
        matches = left.get(key)
        if matches:
            if len(matches) != 1 or hit is not None:
                return None
            hit = matches[0]+choices
    if hit is None:
        return None
    partial = hit[0][0]
    for v, _ in hit[1:]:
        if not partial or probe.curve[partial][0] == probe.curve[v][0]:
            return None
        partial = (partial+v) % probe.order
    return hit


class PinnedProbe(FieldProbe):
    def __init__(self, p, linear, count, positions, target):
        super().__init__(p, linear, count, positions, target)
        # Each pinned point uses its ordinary abscissa address, with no sign bit.
        self.block = self.width
        self.size = self.width**count
        self.projection_checks = 0

    def check_projection(self, slots, hit, moments):
        original = copy(self)
        original.block = 2*self.width
        full = walsh(original.grids(slots)) % self.p
        signs = [digit % 2 for _, digit in hit]
        inverse = pow(2**self.count, -1, self.p)
        for xi in (0, self.size//3, self.size-1):
            total = np.zeros(3, dtype=np.int64)
            for bits in product((0, 1), repeat=self.count):
                r = 0
                for i, bit in enumerate(bits):
                    address = xi//self.width**(self.count-i-1) % self.width
                    r = r*original.block+2*address+bit
                total += (-1)**sum(s*b for s, b in zip(signs, bits))*full[:, r]
            assert np.array_equal(total*inverse % self.p, moments[:, xi])
            self.projection_checks += 1

    def draw(self, rng, anchor_rng):
        while True:
            slots = self.plant(rng)
            hit = unique_relation(self, slots)
            if hit is None:
                continue
            # A regular original grid has distinct abscissas between slots
            # and avoids the endpoint, so the unique witness is its only zero.
            xs = [{self.curve[v][0] for v, _ in slot} for slot in slots]
            if (sum(map(len, xs)) != len(set.union(*xs))
                    or any(self.curve[self.target][0] in part for part in xs)):
                continue
            break
        pinned, witness = [], 0
        for group, (_, digit) in zip(slots, hit):
            pin = digit % 2
            pinned.append([(v, address//2) for v, address in group if address % 2 == pin])
            witness = witness*self.width+digit//2
        grid = self.grids(pinned)
        legal = grid[0].astype(bool)
        assert legal[witness] and grid[1, witness] % self.p == 0
        assert np.count_nonzero(legal & (grid[1] % self.p == 0)) == 1
        mass, first, second = moments = walsh(grid) % self.p
        if not self.projection_checks:
            self.check_projection(slots, hit, moments)
        if self.checks < 8:
            self.check_factorization(pinned, self.checks*(self.size-1)//7, moments)
        scores = np.stack((self.chi[first], self.chi[second],
                           self.chi[(second+first**2) % self.p],
                           self.chi[(mass*second-first**2) % self.p]))
        raw_coefficients = walsh(scores)
        # Recenter at an independent uniform legal address of this root.
        # For any old address v the new address is v XOR anchor. All signed
        # input moments acquire the same Walsh character of the anchor.
        anchor = 0
        for slot in pinned:
            anchor = anchor*self.width+anchor_rng.choice(slot)[1]
        shifts = np.arange(self.size) ^ anchor
        plus = np.array([((r & anchor).bit_count() % 2 == 0) for r in range(self.size)])
        reflected = self.chi[(first**2-second) % self.p]
        centered_mix = np.where(plus, scores[2], reflected)
        centered = np.stack((walsh(centered_mix)[shifts], raw_coefficients[3, shifts]))
        offset_hull = False
        hull_size = 1
        for slot in pinned:
            first_address = slot[0][1]
            span = {0}
            for _, address in slot:
                diff = address ^ first_address
                span.update([v ^ diff for v in span])
            offset_hull |= first_address not in span
            hull_size *= len(span)
        if offset_hull:
            even_rows = [0, 1, 3] if self.p % 4 == 1 else [3]
            assert not np.any(raw_coefficients[even_rows][:, legal])
            assert not np.any(walsh(scores[2]+reflected)[legal])
        coefficients = np.abs(np.vstack((raw_coefficients, centered)))/self.size
        true = coefficients[:, witness]
        legal[witness] = False
        rich = min(map(len, pinned)) >= 2
        assert legal.sum() == prod(map(len, pinned))-1
        if not legal.any():
            return None, rich, 1, offset_hull
        wrong = coefficients[:, legal]
        ties = (true[:, None] == wrong)
        rank = np.mean((true[:, None] > wrong)+0.5*ties, axis=1)
        average = wrong.mean(axis=1)
        data = np.stack((true, average, true-average, rank), axis=1)
        # Integrate out the choice of origin, using v XOR c in the even
        # part and v in the odd part. Put the witness in the first row.
        occupied = np.r_[witness, np.flatnonzero(legal)]
        relative = np.bitwise_xor.outer(occupied, occupied)
        reflected_coeffs = walsh(reflected)
        even = raw_coefficients[2]+reflected_coeffs
        odd = raw_coefficients[2]-reflected_coeffs
        all_mix = np.abs(even[relative]+odd[occupied, None])/(2*self.size)
        all_disc = np.abs(raw_coefficients[3, relative])/self.size
        chosen_col = np.flatnonzero(occupied == anchor)[0]
        for row, matrix in ((4, all_mix), (5, all_disc)):
            assert np.array_equal(matrix[:, chosen_col], coefficients[row, occupied])
            marked, rivals = matrix[0], matrix[1:]
            a, b = marked.mean(), rivals.mean()
            rank = np.mean((marked > rivals)+0.5*(marked == rivals))
            if row == 5 and hull_size == len(occupied):
                assert a == b and rank == 0.5
            data = np.vstack((data, [a, b, a-b, rank]))
        return data, rich, wrong.shape[1]+1, offset_hull


def run(p, linear, count, positions, target, seed, draws):
    probe = PinnedProbe(p, linear, count, positions, target)
    rng = Random(seed)
    anchor_rng = Random(seed+1879)
    cases = [probe.draw(rng, anchor_rng) for _ in range(draws)]
    assert probe.checks == 8 and probe.projection_checks == 3
    print(f"p={p}, m={count}, w={probe.width}, T=[{target}]G, {draws} regular singleton draws; "
          f"{sum(v is None for v, _, _, _ in cases)} with no wrong pinned tuple, "
          f"{sum(flag for _, _, _, flag in cases)} in an affine hull excluding zero.")
    for name, rich_only in (("with alternatives", False), ("at least two per slot", True)):
        kept = [(v, size) for v, rich, size, _ in cases if v is not None and (not rich_only or rich)]
        assert len(kept) > 10
        values = np.stack([v for v, _ in kept])
        means = values.mean(axis=0)
        errors = values.std(axis=0, ddof=1)/sqrt(len(values))
        median = np.median([size for _, size in kept])
        print(f"  {name}: n={len(values)}, median pinned tuples {median:.0f}.")
        for i, formula in enumerate(("chi(M1)", "chi(M2)", "chi(M2+M1^2)", "chi(M0*M2-M1^2)",
                                     "recentered mix", "recentered discriminant",
                                     "mix all legal origins", "discriminant all legal origins")):
            a, b, gap, rank = means[i]
            print(f"    {formula}: witness {a:.5f}, wrong {b:.5f}, gap {gap:+.5f} "
                  f"(SE {errors[i,2]:.5f}), rank {rank:.4f} (SE {errors[i,3]:.4f}).")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run with assertions enabled.")
    run(1031, 9, 4, (0, 3), 1, 75403, 600)
    run(65521, 44, 6, (0, 3), 1, 451029, 600)
    run(65521, 44, 6, (0, 3), 2, 942791, 600)
