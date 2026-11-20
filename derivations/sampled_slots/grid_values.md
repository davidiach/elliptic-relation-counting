# What the sampled grid preserves about the elliptic relation polynomial

Fix a smooth short Weierstrass curve in characteristic $p>3$, an affine target $T$ with $2T\ne O$, and $m\ge2$. Let $F_T(X_1,\ldots,X_m)$ be an irreducible equation for the abscissas of ordered points with sum $T$, over the algebraic closure. Fix any one nonzero scalar normalization. If $T$ is rational the polynomial can be normalized over $\mathbb F_p$. We examine its **actual field values** on the product of the sampled abscissa sets. [Semaev's summation polynomials](https://eprint.iacr.org/2004/031)

## A top coefficient survives restriction to legal points

Put $d=2^{m-1}$. The degree of $F_T$ in every variable is $d$, and its coefficient at $\prod_i X_i^d$ is nonzero. Here is a proof covering the fixed-target specialization. Fix generic abscissas in the first $m-1$ slots. Their independent sign flips form a group of order $d$. Its orbit on the last abscissa consists of the values $x(T-\sum_{i<m}\epsilon_iP_i)$. Equality of two such values requires the two points to be equal or opposite. In the former case a nonempty signed subset has double zero. In the latter the doubled signed sum of the unchanged points equals $2T$. A nonempty signed subset varies nontrivially with one point, and the empty unchanged set would imply $2T=O$. Thus the orbit has size $d$. Symmetry gives this degree in all variables.

The complete sum fiber is isomorphic to $E^{m-1}$. Its image under the finite abscissa map to $(\mathbb P^1)^m$ is closed and irreducible. Its equation is the separate homogenization of $F_T$, since the exact individual degrees ensure no boundary divisor divides that equation. The point $(\infty,\ldots,\infty)$ is absent: its only possible preimage is the all-identity tuple. Evaluating the homogenization there gives precisely the claimed top coefficient. In particular $\deg F_T=md$.

Let $A_i\subset\mathbb F_p$ be the *distinct, liftable* abscissas in slot $i$. The choices may depend arbitrarily on the public seeds. If every $|A_i|>d$, the unique representative of $F_T$ reduced modulo the grid polynomials $\prod_{a\in A_i}(X_i-a)$ is $F_T$ itself. Successive univariate division gives a representative of degree below $|A_i|$ in each variable and never raises total degree. A reduced polynomial vanishing on the full grid is zero, by repeated univariate interpolation. Thus any polynomial agreeing with $F_T$ at every grid point has total degree at least $md$.

There is also an agreement bound. Let $h>0$ and suppose all $|A_i|\ge h$. For any polynomial $G$ of total degree less than $md$, even one chosen after seeing the whole seed, independent uniform $a_i\in A_i$ satisfy

$$
\Pr\!\left[G(a_1,\ldots,a_m)=F_T(a_1,\ldots,a_m)\right]\le\min(1,md/h).
$$

Indeed $F_T-G$ has a nonzero top coefficient. A nonzero polynomial of degree $D$ on independent uniform finite sets of sizes at least $h$ vanishes with probability at most $D/h$. Induct on the variables, first excluding zeros of the leading coefficient in the last variable and then applying its univariate root bound. This argument concerns abscissas chosen uniformly from a *realized* grid, without a condition that the tested tuple is a witness.

## Probability on seeds that actually contain a guarded relation

Let $N=\#E(\mathbb F_p)$ and let the slot seeds be independent. Each affine point is included with the same marginal probability $\rho$ in each slot. Put $\nu=N\rho$ and $\lambda=N^{m-1}\rho^m$. Let $H$ be the ordinary ordered signed relation count, $X\le H$ the count after the addition guards, and $\kappa=\Pr(X>0)>0$. For every integer $h\ge1$,

$$
\Pr\!\left[\min_i |A_i|\le h\mid X>0\right]\le
\frac{2mh\lambda}{\nu\kappa}.
$$

Condition just on slot $j$. Each of its at most $2|A_j|$ affine points has at most $N^{m-2}$ choices of the other points summing to the prescribed target. Their joint inclusion probability is $\rho^{m-1}$ by independence. On $|A_j|\le h$ this gives $\mathbb E[H\,1_{\{|A_j|\le h\}}]\le2hN^{m-2}\rho^{m-1}$. Sum in $j$, use $1_{\{X>0\}}\le H$, and divide by $\kappa$. This proof requires no independence of two points in the same slot.

At balanced dense widths with $\log_2 \nu=b/m+O(1)$, $\lambda\le 2^{O(m)}$ and $\kappa^{-1}\le 2^{O(m)}$, take $h=\lfloor\sqrt{\nu}\rfloor$. When $m^2=o(b)$, outside conditional seed probability $2^{-b/(2m)+O(m+\log m)}$, every degree-below-$md$ predictor has field-value agreement at most $2^{-b/(2m)+O(m+\log m)}$. These width and coverage conditions include the paid regime of the fixed-bit sampler with one independent affine anchor. They have been stated explicitly so the probability conclusion can be used with other seed laws.

## Zero-only checks can be much shorter

There is a finite-field restriction even on a zero-only filter if it must choose from a small predeclared family. Fix the curve and target and let $\mathcal G$ contain at most $2^a$ polynomials over $\overline{\mathbb F}_p$ of total degree at most $e$. Suppose none is divisible by $F_T$ in $\overline{\mathbb F}_p[X_1,\ldots,X_m]$, so the identically zero polynomial is excluded. A seed-dependent procedure can select any family member for free. For independent slots with affine point marginals $\rho_i$, the probability that a genuine target tuple in the sampled sets passes even one of the filters $G=0$ is at most

$$
2^a\,2^m\, (md)e\, p^{m-2}\prod_i\rho_i .
$$

Indeed, for every nonconstant $G$ in the family, the intersection $V(F_T,G)$ has dimension at most $m-2$ because $F_T$ is irreducible over the algebraic closure. Its cumulative degree is at most $(md)e$ by affine Bézout. An affine variety of cumulative degree $\Delta$ and dimension $r$ has at most $\Delta p^r$ rational points, including for varieties defined over the algebraic closure. A nonzero constant has empty zero set. Each rational abscissa tuple lifts to at most $2^m$ signed tuples and each genuine one survives with probability $\prod_i\rho_i$. Apply the first-moment bound and sum over the family. Independent randomization of the family before the draws of the slot seeds is allowed by conditioning. [Lachaud and Rolland, Theorems 1.1 and 2.1](https://arxiv.org/abs/1405.3027)

For $\lambda=N^{m-1}\prod_i\rho_i=2^{o(b)}$, $m=o(b)$, and $a+\log_2(1+e)=o(b)$, the bound is $2^{-b+o(b)}$ by Hasse's estimate. In particular every nonzero filter of degree below $md$ meets the nondivisibility hypothesis. This result limits finding even one genuine hit with such a family. It assigns no complexity penalty to individual polynomial evaluations. Allowing a full field word of seed-dependent advice, an unrestricted operation on the seeds, or a polynomial divisible by the true relation can circumvent its hypotheses.

Suppose exactly one abscissa tuple on a sampled grid hits. For each $a_i\in A_i$, use its public slot address $0\le t_i<w_i$; every legal abscissa in a slot has a unique address. In the finite model a single affine equality on those addresses can isolate the marked tuple. If $\prod_i w_i\le p$, pack $J=\sum_i t_i\prod_{j<i}w_j$ into $\mathbb F_p$ and test $J=J_*$, where $J_*$ is the packed address of the hit. More generally use a bounded number of base-$p$ limbs when the address range exceeds $p$. Obtaining $J_*$ already determines all the marked addresses by ordinary mixed-radix division. This example explains the limitation of the value theorem: exact numerical evaluation, zero-only tests tailored to the hits, and finding the hits are different tasks. No lower bound for arbitrary arithmetic circuits, gates with auxiliary witnesses, or charged running time follows from the grid degree.

The executable fixture extracts the full mixed leading coefficient from actual grids of liftable abscissas over $\mathbb F_{101}$, cross-checks every zero by signed point addition, and exhausts a small mixed seed law to check the conditional estimate. [Finite checks](../../checks/check_grid_values.py)

A distinct argument applies directly to zero-only verification if a sequential chain decodes and transmits the exact partial sum at each prefix. [Fixed-target dictionary bound](prefix_dictionary.md)
