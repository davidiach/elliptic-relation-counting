# Fresh uniform point weights conceal the whole moment vector

One way to combine the signed moments is to first randomize the measurements. Give each retained **point** $P$ in slot $i$ an independent uniform coefficient $u_{i,P}\in\mathbb F_p$. Multiply the ordinary Walsh sign $(-1)^{r_i\cdot v(P)}$ by that coefficient. In this variant, when the coefficients are private to the evaluator and refreshed on each call, even an arbitrary nonlinear readout of all the raw field values gains exponentially little information about the isolated address. The deterministic signed moments and coefficients revealed to the readout have different observation laws.

There is an exact sign invariance on **every** fixed root and target. For any fixed challenge, multiplying the array of private weights by its array of Walsh signs is a bijection preserving their product-uniform law. Hence the joint distribution of all moment orders is identical for all challenges. A rule seeing only the vector (or also the fixed root and target) has a challenge-independent mean and its coefficient is zero at every nonzero address. The near-uniformity below gives a stronger comparison for readouts that may also see the challenge: on typical singleton roots the vector can be replaced by independent uniform field values even while conditioning on the root, target and challenge.

## Uniformity on a fixed witnessed grid

Let $d=m+1$. Fix a smooth curve, a nonidentity target $T$, and slots whose abscissas are disjoint between slots and avoid $x(T)$. These are the regular grids. Assume the slots contain at least one relation of sum $T$, and some slot $i$ contains more than $kd$ points, for $k\ge1$. Fix **any** sign challenge $r$. For $0\le j\le k$, let

$$
Z_j=\sum_{\mathbf P}\Delta_{\mathbf P}(T)^j
       \prod_h u_{h,P_h}(-1)^{r_h\cdot v(P_h)}
$$

with $\Delta^0=1$, including at zero. On these fixed data, the vector $(Z_0,\ldots,Z_k)$ is within $ (k+1)(m-1)/p$ in total variation of uniform on $\mathbb F_p^{k+1}$. The same statement for $(Z_1,\ldots,Z_k)$ has bound $k(m-1)/p$.

Freeze the coefficients in all slots except $i$. As a function $B_j(P)$ of the point in slot $i$, take the sum over the other slots with their weights, leaving out the sign and weight in slot $i$. Then $B_j$ belongs to $L(jdO)$. In the elliptic evaluation basis there is exactly one element with pole order $d$ at $O$. The leading coefficient of $\Delta$ there, as a function of $P$, is up to a fixed nonzero scalar the evaluation minor after deleting that basis element and the row at $P$.

Choose the other $m-1$ points from a witnessed tuple and include the endpoint $-T$. These $m=d-1$ points are distinct by regularity. Their sum is $-P_i\ne O$. Therefore their evaluation matrix on $L((d-1)O)$ is nonsingular. A vanishing determinant would give a nonzero function with exactly those $d-1$ zeros and pole divisor $(d-1)O$, implying their group sum is $O$. The retained tuple consists of affine points, so its omitted point cannot be the identity. [Elliptic determinant criterion](../joint_parity/target_incidence.tex)

For every $j\ge1$ that tuple contributes a nonzero coefficient to the pole-$jd$ term of $B_j$, viewed as a polynomial in the other slots' independent weights. Each tuple contributes a different monomial, including the nonzero Walsh signs. The coefficient polynomial has total degree $m-1$. The corresponding constant coefficient of $B_0$ also contains a nonzero monomial. By the polynomial root bound and a union bound, every one of these coefficient polynomials is nonzero at the sampled weights except with probability at most $(k+1)(m-1)/p$. On the good event the functions $B_j$ have distinct pole orders $0,d,\ldots,kd$ and are linearly independent. A nontrivial linear combination has at most $kd$ zeros, so the $k+1$ evaluation rows remain independent on $F_i$. Multiplying columns by the fixed signs in slot $i$ preserves rank. The still-independent uniform vector $(u_{i,P})_{P\in F_i}$ therefore maps onto $\mathbb F_p^{k+1}$ with equal-size fibers. The claimed total variation bound follows.

## Exceptional roots and adaptive queries

Write $n$ for the number of uniformly masked eligible nonidentity subgroup targets. The target is independent of the root. Use the common marginal retention probability $\rho$ of the mixed sampler. The unconditional chance of an irregular root is at most $\binom m2 p\rho^2+m\rho$ by independence across slots and the single-abscissa marginal. If every slot has at most $kd$ points, then for fixed seeds the full grid reaches at most $(kd)^m$ targets. Write $U$ for a unique guarded relation with no other ordinary relation and $q=\Pr(U)>0$. For each challenge this gives the following bound on the average conditional distance from an independent uniform vector:

$$
\frac{(k+1)(m-1)}p+
 \frac{\binom m2p\rho^2+m\rho+(kd)^m/n}{q}.
$$

With $k\le 2^{b/m^2}$, $\log_2 n=b-O(1)$, $\ell=b/m+O(1)$, and the isolating law $q^{-1}=2^{O(m+\log m)}$, this distance is $2^{-b+o(b)}$ for growing $m=O(\sqrt b)$. Any bounded real readout loses at most twice this distance when its observation is replaced by an independent uniform field vector. This comparison permits other side information determined by the root, target and challenge in both experiments; it bounds only the additional information supplied by the observation. If its measurement coefficients are hidden, the conclusion applies even when reading arbitrary field bits or testing target-specific buckets. If the rule sees only the measured vector, its uniform-input baseline is independent of the challenge, so its normalized Walsh coefficient vanishes at each nonzero address. The all-zero singleton address contributes at most $1/(nq)$ conditional on $U$ by uniqueness of its target for each seed.

For a sequence of $Q$ adaptive calls on the same root and target, use fresh independent coefficients on every call. On good roots couple each call to independent uniform output until the first discrepancy. This multiplies the first term by $Q$ and leaves the exceptional-root term unchanged. Reusing a single set of hidden coefficients across challenges or making the coefficients visible changes the joint law and lies outside this argument. [Finite rank and uniform-projection check](../../checks/check_hidden_weights.py)
