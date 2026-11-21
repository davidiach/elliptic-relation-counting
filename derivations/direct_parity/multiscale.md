# A bounded nonlinear parity score from coarse half-sum matches

Suppose a coarse test of the group equality leaves $K$ legal candidate addresses, including the true address. When those candidates come from joining two bucketings, a bounded nonlinear score achieves expected absolute Walsh correlation at least $1/(4\sqrt K)$ at **each** candidate. Evaluating the score needs the two half lists and their bucket memberships. The score itself can avoid expanding the $K$ pairs. The constant-factor rate is optimal if all $K$ candidates need the same guarantee. Uniform hashing quantifies how fine the buckets must be for this guarantee to meet the outer decoder's correlation requirement.

## Combining the half-address votes

Split the slots into two blocks. A left half tuple $u$ has group label $P(u)$, and a right half tuple $v$ has group label $Q(v)=T-\sum_{i\text{ right}}P_i$. Their signed binary addresses occupy disjoint blocks in the complete address. Assign the two sets of labels to buckets by the same function. Drop buckets lacking a half tuple on either side. For a surviving bucket $j$ write $l_j,r_j$ for its two numbers of tuples, $k_j=l_jr_j$ and $K=\sum_j k_j$. Repeated labels within a half still count as separate tuples.

Give every **half tuple** an independent uniform sign $e_u$ or $d_v$, and give every surviving bucket an independent uniform sign $z_j$. Sample the signs once for this score and reuse them for its Walsh challenges. At challenge $r$, let $s_j(r)$ be the product of the two votes in bucket $j$. A left vote is the sign of $\sum_{u\in j}e_u\chi_{a(u)}(r)$, with ties assigned $+1$. Use $d_v$ for the right vote. Put $t_j=k_j/(2K)$ and set

$$
f(r)=\operatorname{Im}\prod_j\left(\sqrt{1-t_j}+i z_j s_j(r)\sqrt{t_j}\right).
$$

The factors have absolute value one, so $|f(r)|\le1$. A grid with $K=0$ can return zero. The products and vote counts use time proportional to the total input tuples and occupied buckets, up to polynomial factors in the problem bit length. Computing the factor square roots and complex product to $O(\log(2K))$ bits of absolute precision per operation suffices for uniform error at most $1/(8\sqrt K)$: there are at most $K$ unit-modulus factors, and the elementary product-error bound sums their absolute errors up to a factor of two if their total is small. Clip the final approximation to $[-1,1]$. This preserves at least half the correlation guarantee below.

Fix **any** candidate $(u,v)$ in bucket $j$. Its address $a$ has $\chi_a=\chi_{a(u)}\chi_{a(v)}$. For independent signs $X_1,\ldots,X_n$ the correlation of $X_1$ with the sign of their sum is

$$
c_n=2^{1-n}\binom{n-1}{\lfloor(n-1)/2\rfloor}\ge 1/\sqrt{2n}.
$$

Condition on the other signs to obtain the equality. The lower bound follows from $\binom{2k}{k}/4^k\ge1/(2\sqrt k)$ for $k\ge1$, by induction on $k$. For every fixed challenge the Walsh factors just flip independent signs. Multiplying $f$ by $z_j$ and averaging over all the bucket signs retains just the linear term in that factor. Averaging its two votes over the half signs gives the exact identity

$$
\mathbb E_{e,d,z}[z_j e_u d_v\,\widehat f(a)]
=\sqrt{t_j}\prod_{j'\ne j}\sqrt{1-t_{j'}}\;c_{l_j}c_{r_j}.
$$

Here $\widehat f$ is the normalized address-Walsh coefficient. Since $\sum_j t_j=1/2$, the product of square roots is at least $1/\sqrt2$. Substitution yields $\mathbb E|\widehat f(a)|\ge1/(4\sqrt K)$. The aligning sign is used only in this proof. The absolute-correlation decoder never receives the identity of the true tuple. The expectation includes the score coins fixed for each run. On each root a randomly chosen score can therefore act as the fixed oracle for a run of the [address decoder](../sampled_slots/address_parity.md).

For any bounded score on the same challenge space, Parseval gives $\sum_{a\text{ candidate}}|\widehat f(a)|^2\le1$, since the $K$ addresses are distinct. Thus the average absolute coefficient over these candidates is at most $1/\sqrt K$, even with score randomness. This upper bound concerns a guarantee for every candidate or for a uniformly hidden candidate. It permits higher correlation at the actual witness when other information in the seeds distinguishes it.

## Hash width, setup cost and a one-slot-at-a-time baseline

On a fixed ordinary singleton exactly one left-right pair has equal group labels. Let $D$ be the total number of legal tuples. Hash the labels with a fresh pairwise independent uniform $B$-bucket map. For example, when $B$ is a power of two, use a random binary affine map on an injective binary encoding of the group points. Each of the $D-1$ false pairs becomes a coarse match with probability $1/B$. Thus $\mathbb E K=1+(D-1)/B$. Jensen's inequality proves for **each actual singleton**

$$
\mathbb E_{\text{hash},e,d,z}|\widehat f(a_*)|
\ge\frac{1}{4\sqrt{1+(D-1)/B}}.
$$

On regular guarded singletons this is correlation with the signed lift of $M_0-M_{p-1}$. At the thinned widths $D\le(2w)^m=2^{b+O(m)}$. To obtain a guaranteed $2^{-O(b/m)}$ correlation through this inequality on a root of size $D$, it suffices to use $B\ge D/2^{O(b/m)}$. The guarantee describes a weak predictor used with list decoding, not a high-accuracy parity classifier. Even with $B$ this large, the occupied buckets can be stored sparsely. Building the two half lists and assigning their independent signs still takes $2^{b/2-o(b)}$ time on typical isolated roots at balanced splits. Each half has that many distinct sums under the [sampled-slot displacement estimate](scalar_join.md). The calculation provides no lower bound on implicit arithmetic that recovers the useful buckets or exploits the true witness's distribution.

There is also an inexpensive explicit score obtained by taking the **product of one independently signed majority per individual slot**. For every legal tuple, the same sign-alignment calculation yields expected absolute correlation at least $2^{-m/2}/\sqrt D$. It uses $O(\sum_i |F_i|)$ signed terms per challenge once the slots are available. With $D\le(2w)^m$ this gives a fixed-target singleton lower bound of $2^{-b/2-O(m)}$. The [separated-score inequality](../sampled_slots/address_parity.md) supplies the matching upper scale $2^{-b/2+O(m+\log m)}$ for this product on the actual sampler. Thus the inexpensive baseline attains the optimal leading correlation exponent for these separated scores, which is too small for the decoding budget. The address parity remains computable exactly with a balanced half-list join. Neither construction improves the end-to-end exponent.

## Visiting every jointly occupied prefix

Refining hashes a few bits at a time can also require a large intermediate list. On a singleton, the sets of distinct left and right group labels share exactly one label. Delete it from both. Let $B$ be the largest power of two bounded by the smaller of the two remaining cardinalities. From each side retain any $B$ distinct labels, selected before the new hash. These $2B$ keys are disjoint. At depth $\log_2 B$, let $X,Y$ count hits from these two sets in a given prefix. For a four-wise independent prefix hash,

$$
\mathbb E[XY]=1,\qquad
\mathbb E[X^2Y^2]=(2-1/B)^2\le4 .
$$

To obtain the second equality, expand the four indicators. The only repeated indices can be two draws on the same side. This gives exactly the product of the two binomial second moments. Cauchy–Schwarz now gives $\Pr(XY>0)\ge1/4$. Summing over prefixes, at least $B/4$ are jointly occupied in expectation. Counts can only increase when omitted labels are restored. A four-wise independent finite hash family with polynomial seed length is obtained by evaluating a uniformly random cubic over a sufficiently large binary field at distinct encoded keys and taking output bits.

At a balanced split on typical successful roots $B=2^{b/2-o(b)}$. An exhaustive traversal that visits every jointly occupied prefix up to an eventual fine hash incurs this expected width even if each occupancy test is free. A search that can identify which shared prefixes contain the true equality can skip other branches. The result does not bound that selection task. [Finite checks of the score and occupancy identities](../../checks/check_multiscale.py) · [Output](../../evidence/check_multiscale.txt)
