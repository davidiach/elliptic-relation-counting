# False addresses from fingerprinting intermediate group sums

Consider a binary tree doing exact complete group additions. At each of $t$ specified internal links let the producer output $U_i$ and let the consumer guess any $V_i$ in the full rational group $G$ of size $N$. They need only agree under a fingerprint $h_i$ of range size $q_i$. The root checks its exact target $T$. The distinct fingerprints are independently sampled, with $\Pr[h_i(X)=h_i(Y)]=1/q_i$ for each fixed unequal pair. Random binary linear maps on an injective point encoding have this property for power-of-two ranges and polynomial-size descriptions. All seeds, the target and the proposed signed address may be public. The fingerprint coins are independent of them.

For one fixed tuple of true sum $S$, define $d_i=V_i-U_i$. Given the $d_i$, all produced and consumed group values are uniquely determined from the leaves. The required root value gives $\sum_i d_i=T-S$, even when two modified links lie on the same path. Each assignment of jumps has acceptance probability $\prod_{i:d_i\ne O}q_i^{-1}$, since all of its group values are fixed before the fingerprint coins. Let
$$A=\prod_i(1-1/q_i),\qquad B=\prod_i(1+(N-1)/q_i).$$
The expected number of internal completions of the tuple is exactly
$$ A\,1_{S=T}+(B-A)/N.$$
Indeed each link contributes the group function $(1-1/q_i)\delta_O+q_i^{-1}{\bf1}_G$. Its convolution with the next such function retains a point mass and a constant. After all links the point-mass coefficient is $A$ and total mass is $B$. This count includes different auxiliary completions of the same address.

There is a complementary bound for **distinct** wrong addresses. If $S\ne T$, assign the entire error $T-S$ to one chosen link and zero to all others. Only one fingerprint collision is required. The $t$ such events are independent because each uses a different hash and its own deterministic pair of points. At least one occurs with probability $1-A$. Therefore for any fixed singleton grid of $D$ legal signed addresses, the number $F$ of accepted wrong addresses satisfies
$$
(D-1)(1-A)\ \le\ \mathbb E_h F\ \le\ (D-1)(B-A)/N.
$$
All witnesses survive by taking all jumps zero. The bounds apply to this unrestricted group-state verifier. Intersections with other independent soundness checks, states restricted to reachable sums by extra local tests, and fingerprints designed from the particular seed and target have different constraints.

The upper bound can meet the requested statistical budget. For instance take $q_i=2^{\lceil\log_2 N\rceil}\ge N$ and $t=m-2$. Then $B\le2^{m-2}$. At the project's thinned widths $D/N\le(2w)^m/N\le(1+o(1))\,2^{-m}/m^2$, so the expected number of wrong addresses is $O(m^{-2})$ on every fixed singleton. The full widths remain essential in this model. For each $i$, $1-A\ge1/q_i$, so the desired conditional mean below $1/3$ requires $q_i>3\,\mathbb E[D-1\mid U]$ for fixed range sizes. On typical guarded singletons $D\ge2^{b-o(b)}$, by the same good-slot event as in the [half-size estimate](../direct_parity/scalar_join.md). Thus each unrestricted link needs $b-o(b)$ fingerprint bits. Even an arbitrarily complicated implementation of the local additions would still have these label widths if it retains this model. The bound does not apply to general joint seed arithmetic.

The fixture enumerates complete hash families and all legal signed addresses on a regular order-19 curve example. It checks two tree shapes, the expected number of auxiliary completions, and the probability of acceptance by at least one single jump. [Checker](../../checks/check_candidate_splices.py) · [Output](../../evidence/check_candidate_splices.txt)
