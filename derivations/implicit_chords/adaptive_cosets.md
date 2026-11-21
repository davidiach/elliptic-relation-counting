# How much of a sampled half can adaptive cosets cover?

An isogeny can replace a complete kernel coset by one quotient evaluation in the chord calculation. Could one choose different kernels after seeing the seeds, covering the half by whichever cosets are most useful? The following estimate allows every rational subgroup at once, including on composite-order curves. It also allows cosets with holes. [Exact trace formula](algebra.md)

Let $G=E(\mathbb F_p)$ have order $N$, $b=\lceil\log_2 p\rceil$, and take one balanced block $I$ of $k=\lfloor m/2\rfloor$ or $\lceil m/2\rceil$ independent mixed slots. Use the thinned widths, growing $m=O(\sqrt b)$ and the free-bit run bound $s=o(b)$ in the isolation regime. Let $S_I$ be the set of **distinct** block sums. Choose $0<\eta\le1$ with $\log_2(1/\eta)=o(b)$. Count a point in $S_I$ as covered if it belongs to *any* coset of a nontrivial subgroup $H$ containing at least $\max(2,\eta|H|)$ **distinct** points of $S_I$.

With probability $1-o(1)$ conditional on the fixed-target guarded singleton event, the total number of covered points is at most $2^{b/4+o(b)}$. The block itself has at least $2^{b/2-o(b)}$ distinct sums with the same conditional probability. In particular the conclusion holds if a method can mix different isogenies, search for cosets adaptively, or use a density as small as $2^{-o(b)}$. The count concerns coverage by dense rational subgroup cosets. It charges neither arbitrary implicit evaluations nor identities between overlapping translates. [Widths and isolation probability](../sampled_slots/isolation.md) · [Distinct sums on a successful half](../direct_parity/scalar_join.md)

## A displacement estimate for any block

For a block $J$ let $\nu_J=(Nw/2^b)^{|J|}$. For each displacement $W\in G$, the expected number of ordered pairs of **different tuples** in $J$ with difference of sums $W$ is at most

$$
\beta_J=C^{|J|}|J|\,2^s\,\nu_J^2/N .
$$

For $W=0$ the identical tuple is removed. For every remaining pair some slot displacement is nonzero. Mark one such slot, fix the other displacements, and bound the forced displacement at the marked slot by its uniform nonzero-displacement bound. Sum the masses in the other slots. This is the same estimate for a prescribed nonzero $W$ as for a collision. For any fixed subgroup of size $e$, let $J_H$ denote the number of ordered different tuples of this block agreeing modulo $H$. Summing over its $e$ displacements gives $\mathbb E J_H\le e\beta_J$. On the balanced block, $\beta_I\le2^{o(b)}$. [Slot displacement proof](../sampled_slots/prefix_dictionary.md)

## Large subgroups require collisions in both quarters

Split $I$ into two independent blocks $A,B$ whose arities differ by at most one. Let $n_A,n_B$ be their actual tuple counts. Each is at most $L=(2w)^{\lceil k/2\rceil}=2^{b/4+o(b)}$. In $G/H$ let $a_z,b_z$ be their tuple frequencies, and write $J_A=\sum_z a_z(a_z-1)$ and $J_B=\sum_z b_z(b_z-1)$. The number of full tuples in **any** coset $c$ obeys

$$
\sum_z a_z b_{c-z}\ \le\ n_A+n_B+\sqrt{J_A J_B}.
$$

Indeed, write each frequency as its occupancy indicator plus its excess over one. Terms with at least one indicator have total at most $n_A+n_B$. The remaining inner product is bounded by Cauchy–Schwarz, since $(a-1)_+^2\le a(a-1)$ for every nonnegative integer $a$.

Take $e>4L/\eta$. A coset containing at least $\eta e$ distinct half sums also contains at least that many tuples, so $J_A J_B\ge(\eta e/2)^2$. The seeds of the two quarters are independent **before conditioning**. Consequently the probability that even one such coset exists for this fixed $H$ is at most

$$
\frac{4\,\mathbb E J_A\,\mathbb E J_B}{\eta^2 e^2}
\ \le\ \frac{4\beta_A\beta_B}{\eta^2}
\ \le\ 2^{-b+o(b)}.
$$

Here $\nu_A\nu_B=\nu_I=2^{b/2+o(b)}$. No union over coset representatives is needed since the same inequality holds for every $c$. Only $e\le (2w)^k/\eta=2^{b/2+o(b)}$ can have a sufficiently dense coset.

A rank-at-most-two finite abelian group has at most $\sigma_1(e)\le e\tau(e)$ subgroups of order $e$. For example, realize the group as a quotient of $\mathbb Z^2$ and pass to the dual group: the annihilator of an order-$e$ subgroup has index $e$. Its preimage is an index-$e$ sublattice, and the two-dimensional Hermite normal forms count $\sum_{d\mid e}d$ such lattices. Here $\tau$ is the divisor count and $\tau(N)=2^{o(b)}$. Summing over divisors $e\mid N$ below this size gives at most $2^{b/2+o(b)}$ candidate subgroups. The probability that any large subgroup has such a coset is therefore at most $2^{-b/2+o(b)}$.

## Small subgroups can cover few points in total

For a fixed subgroup of order $2\le e\le4L/\eta$, let $t_c$ count the distinct half sums in its coset $c$. In a qualifying coset, $t_c(t_c-1)\ge t_c\,\eta e/2$. Each ordered pair of distinct sums supplies at least one ordered tuple pair with a nonzero displacement in $H$. The expected number of covered points for that subgroup is therefore at most $2(e-1)\beta_I/(\eta e)\le2\beta_I/\eta$. Counting **all** subgroups in this size range costs at most $(4L/\eta)\tau(N)^2=2^{b/4+o(b)}$. The expected total coverage, even with double counting, is at most $2^{b/4+o(b)}$.

The fixed-target event used above has inverse probability $2^{o(b)}$. Dividing the unconditional small-subgroup expectation and large-subgroup failure bound by this probability preserves their displayed exponential orders. Applying Markov to the small-subgroup contribution with an additional growing subexponential factor proves the claim.

One can also choose the balanced split **after seeing the seeds**. Sum the small-subgroup expectations and large-subgroup failure bounds over the at most $2^m=2^{o(b)}$ choices before conditioning. The total still has the same exponents. For the support lower bound, the event that every slot has at least $\nu/m^3$ points is common to all splits, with $\nu=Nw/2^b$. The total expected number of pairs of different tuples colliding at zero across all balanced blocks is at most $\sum_I\beta_I=2^{o(b)}$. Conditioning and Markov bound this total by $2^{o(b)}$ with probability $1-o(1)$. Removing at most that many repeated tuples leaves $2^{b/2-o(b)}$ distinct sums in every balanced block on the same event.

The finite check enumerates every subgroup on curves of orders 19 and 20. It compares coset occupancies with the product of the two quarter collision counts, checks three density thresholds, and computes the independent expectations with their seed multiplicities. [Checker](../../checks/check_quarter_cosets.py) · [Transcript](../../evidence/check_quarter_cosets.txt)
