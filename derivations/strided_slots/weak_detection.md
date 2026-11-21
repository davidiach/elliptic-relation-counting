# Weak planted detection on the original affine hashes

A bounded statistic whose mean shifts under planting can suffice for witness recovery, even if thresholding it is a poor predictor of fresh target existence. The comparison must be made against independently *marked* slots, since planting also biases slot sizes. This note gives the advantage and cost conditions. A statistic meeting them at growing arity is still missing.

## Statement

Use the original product law $P=\bigotimes_i P_i$ of independent odd affine hashes with dyadic widths $w_i$ and $\rho_i=w_i/2^b$. Fix a target with $2T\ne O$. Write $N=\#E(\mathbb F_p)$, $H$ for its ordinary ordered signed count, $X$ for the guarded count, $\lambda=N^{m-1}\prod_i\rho_i$, and $\mu=\mathbb E_P H$. Use balanced sparse widths with $\lambda\asymp1/\log b$, $\mu\sim\lambda$, $\mathbb E_P H^2=O(\lambda)$ and $\mathbb E_P(H-X)/\mu=O(m/N)$. These estimates follow from the existing factorial and guard bounds when $m/\min_i N\rho_i=o(1)$ and $2^m/N=o(1)$. For square-root-in-bitlength arity they hold as well. [Factorial moment estimates](../joint_parity/fresh_windows.tex)

Let $K_{i,x}$ be the seed law conditional on accepting abscissa $x$. Under $Q$, draw an ordered tuple of uniform affine points conditional on sum $T$, then independently draw the seeds from their $K_{i,x_i}$. Under $R$, draw the marked points independently uniformly affine and use the same conditional seeds. If $D=\prod_i|F_i|$, then $dQ/dP=H/\mu$ and $dR/dP=D/\mathbb ED$. Let $R_i$ be the one-slot matched law, so $R=\bigotimes_i R_i$.

Suppose a bounded score $0\le f\le1$ costs $D_m$ per call and
$\left|\mathbb E_Q f-\mathbb E_R f\right|\ge\delta$ for each required target. Assume $m2^m/N=o(\delta)$. Fresh private coins are permitted. Use a known lower bound on $\delta$. There is a verified witness attempt with inverse success and cost bounded respectively by

$$
\mathrm{poly}(\lambda^{-1},2^m,\delta^{-1}),\qquad
(1+\textstyle\sum_i w_i)(1+D_m)\,
\mathrm{poly}(b,m,\lambda^{-1},2^m,\delta^{-1}).
$$

All simulated calls, candidate enumeration and verification are included. In particular, for $m\asymp\sqrt b$, $D_m,\delta^{-1}\le2^{O(\sqrt b)}$ would suffice for the complete $2^{O(\sqrt b)}$ target under the inherited outer reduction. This implication concerns a single scalar score with a uniform targetwise advantage. It requires an evaluator at those parameters. The new coefficient in the [fixed-bit theorem](results.md) does not use it.

## Isolating the joint signal

Write $\Pi_T$ for the planted law of the marked tuple before drawing seeds. Regard the probability measures as averaging operators acting on their own axes, and put

$$
B_{i,x}=(K_{i,x}-R_i)/2,\qquad
g_k(x_1,\ldots,x_k;s_{>k})=
\left(\bigotimes_{i\le k} B_{i,x_i}\right) f(s_{>k}).
$$

The fully mixed coefficient obeys

$$
\left|2^m\mathbb E_{\Pi_T}g_m-(\mathbb E_Q f-\mathbb E_R f)\right|
\le 2\left[(1+1/(N-1))^m-1\right].
\tag{1}
$$

Expand the product of differences. The full $K$ term has law $Q$. Every proper term averaging over $s$ planted coordinates and $t=m-s$ fresh matched ones differs from $R$ by at most $2/(N-1)^t$ in total variation. To see this, the sum of $j$ independent nonidentity points has mass $[1+(-1)^j(N1_{u=O}-1)/(N-1)^j]/N$ at $u$. Conditioning the sum of all $m$ marks to $T$ changes the $s$-mark density according only to whether their sum is $T$. If that event originally has probability $a_s$, its total variation change is $Na_s(1-a_s)/[(N-1)^t(1-(-1)^m/(N-1)^m)]$. This gives the displayed per-term bound for sufficiently large $N$. Passing marks through a sampling kernel cannot increase it. The proper terms' coefficients sum to $-1$ and the empty term has zero error. Summing the remaining errors proves (1). [Proper-marginal calculation](../joint_parity/pooled_zero.tex)

## Keeping one whole planted path

Set $h=\delta/(16\,2^m)$. By (1), $|\mathbb E_{\Pi_T}g_m|\ge8h$ eventually. Draw the root seed vector under $Q$ with its marked tuple. For each slot independently flip a fair coin. With sign $+1$, copy the root seed to a working vector. With sign $-1$, independently resample from $R_i$. Let $U=f(s^{\rm work})\prod_i\sigma_i$. For each fixed marked tuple, expose these coin signs and working seeds from the last slot backwards. At the stage where $k$ prefix slots remain hidden, the conditional expectation of $U$ is $(\prod_{i>k}\sigma_i) g_k$.

Stop at the first $k=m,m-1,\ldots,1$ with $|g_k|<4h$, if one exists. Conditioning on that first visit gives $|\mathbb E[U1_{\rm stop}]|\le4h$. Since $|U|\le1$ and $|\mathbb E U|\ge8h$, the entire path clears the barrier with probability at least $4h$. The chance the planted mark fails its addition guards is $O(m/N)=o(h)$.

Evaluate the same working-vector construction on a fresh root under $P$. Let $G$ count guarded signed root tuples that clear every barrier. The mark weighting identity gives $\mathbb E_P G\ge c\mu h$, for an absolute $c>0$. As $G\le H$, the second moment bound gives $\Pr_P(G>0)\ge c'\lambda h^2$. Denote a guaranteed lower bound of this order by $\zeta$. This argument allows multiple relations in a root and charges full-path survival across all $m$ slots.

## Bounding the number of surviving prefixes

For one slot write $I_x=1_{x\ {\rm accepted}}$ on the set of liftable abscissas. The covariance matrix under $P_i$ is $\rho_i I-\rho_i^2 LL^\top$, where $L$ is the zero-one membership matrix of residue classes modulo $w_i$. Distinct abscissas in the same class cannot both be retained. At different classes their joint acceptance probability is exactly $\rho_i^2$. Thus the covariance operator is bounded above by $\rho_i I$. Taking the squared norm of the adjoint of the centered-incidence operator gives, for any real one-slot function $v$,

$$
\sum_x |K_{i,x}v-P_i v|^2\le\rho_i^{-1}\mathbb E_{P_i}|v|^2 .
$$

The point-mark weights $\pi(x)$ satisfy $1/(N-1)\le\pi(x)\le2/(N-1)$. Centering at $R_iv=\sum_x\pi(x)K_{i,x}v$ minimizes weighted squared deviation. Hence $\sum_x|B_{i,x}v|^2\le(2\rho_i)^{-1}\mathbb E_{P_i}|v|^2$. With no rational nonidentity two-torsion the weights are uniform and the constant improves to $(4\rho_i)^{-1}$. Tensorizing, for every fixed working suffix,

$$
\sum_{x_1,\ldots,x_k}|g_k(x_{\le k};s_{>k})|^2
\le\prod_{i\le k}(2\rho_i)^{-1}.
$$

Under a fresh root the first $k$ masks are independent of the working suffix. The expected number of root abscissa prefixes with $|g_k|\ge h$ is at most $2^{-k}h^{-2}$. Summing over $k$ bounds the expected total by $h^{-2}$. This controls the cost of the tree even when the score is highly irregular.

## Recovery and costs

Generate the working vector once. Extend the live root-abscissa prefixes one slot at a time. Estimate each $g_k$ using independent fair signs and conditional $K_{i,x}$ or matched $R_i$ samples on the prefix, keeping the working suffix fixed. Average the score times the product of these fresh signs. To sample $K_{i,x}$, choose a uniform odd multiplier and a uniform legal address and solve for the offset. To sample $R_i$, first draw a uniform affine point, then use its abscissa. An exact uniform affine sampler proposes a uniform field abscissa and a random index in the two possible sorted roots, rejecting an unavailable index. Hasse's bound gives constant expected proposals. These routines take polynomial expected bit cost independently of the root.

Keep a prefix when its estimate has absolute value at least $2h$. Abort if the total number of retained prefixes exceeds $J=\lceil 10/(\zeta h^2)\rceil$. With every estimate accurate to $h$, this abort event has probability at most $\zeta/10$ by the preceding expectation bound. The cap gives at most $1+(J+1)\sum_i w_i$ tests. Hoeffding concentration with $O(h^{-2}\log((2+(J+1)\sum_i w_i)/\zeta))$ independent draws per test controls all errors by $\zeta/10$, using conditional error bounds on adaptively queried prefixes. Every path above $4h$ survives an accurate test. At full length verify the signed choices, at most $2^m$ per abscissa tuple. Thus with probability $\Omega(\zeta)$ a guarded relation is returned. The simulation setup can be given a global time cap by Markov's inequality at an extra factor $O(1/\zeta)$. The cost bounds in the statement absorb this factor.

If the evaluator has only a fresh mean-time bound, deterministically cap each call at $O(D_m/(\lambda\delta^2))$. Markov bounds fresh timeout mass, Cauchy–Schwarz with $H/\mu$ bounds its planted mass, and $\mathbb E_P(D/\mathbb ED)^2=O(1)$ bounds its matched mass. Choosing the implied constant large preserves at least half the original gap. This gives a worst-case evaluator for the simulation, whose conditional seed calls may be far from fresh. The one-pilot theorem permits the recovered relation to be chosen randomly with any one fixed rule shared by training and online use. Its factor generation, verified linear solve and repetitions cost $2^{O(b/m)}$ in total whenever $m=O(b/m)$, $D_m,\delta^{-1}\le2^{O(b/m)}$, and $\lambda^{-1}=\mathrm{poly}(b)$.

## A tradeoff for independent thinning

Consider a readout that independently retains each root abscissa in slot $i$ with probability $\theta_i$, then tests whether the thinned grid contains any target tuple. Write $\Theta=\prod_i\theta_i$. Under $Q$ the planted tuple survives with probability $\Theta$. The union bound over all target tuples gives $Q(\mathrm{hit})\le\Theta\mathbb E_QH=O(\Theta)$, since $\mathbb E_QH=\mathbb E_PH^2/\mu$. Likewise $R(\mathrm{hit})\le\Theta\mathbb E_RH=O(\Theta\sqrt{\lambda})$ by Cauchy–Schwarz and the bounded second moment of $D/\mathbb ED$. Thus in the balanced sparse regime the gap lies between $(1-O(\sqrt{\lambda}))\Theta$ and $C\Theta$.

To preserve advantage at least $2^{-O(b/m)}$, this scheme may discard only $O(b/m)$ bits of grid size in expectation. Indeed the fresh expected number of all surviving signed tuples is exactly $\Theta(N-1)^m\prod_i\rho_i\sim\Theta N\lambda$. At square-root arity this is at least $2^{b-O(\sqrt b)}$. Simply enumerating that whole thinned grid misses the required budget. The argument leaves compressed existence tests open.

## Finite tests

On the order-19 curve over $\mathbb F_{13}$, complete small seed laws were summed for a synthetic noisy witness control at three and four slots and for the quadratic character of the actual Gram moment $Z_2$ at three slots. The last statistic had planted-minus-matched gap 0.002734 and fully mixed coefficient $2^m\mathbb E g_m=0.002409$. Each case checked the energy inequality and the whole-prefix survival bound. These finite sums use floating arithmetic for the probabilities. [Small-law check](../../checks/check_weak_recovery.py)

The independent larger-field test fixes its readout to whether $Z_2(T)$ is a nonzero square. The Gram formula uses the even and odd blocks of the sign-symmetric moment matrix, and is checked against direct sums of squared determinants on small cases of each arity. This permits sampling the two full seed laws at larger field size without a label oracle. The arithmetic comparison is one scalar readout and does not measure the information in the full vector of moments. [Larger-field experiment](../../checks/probe_weak_gate.py)

For $p=65521,\ a=44,\ T=(0,1)$ on $y^2=x^3+ax+1$, a larger independent run used one million draws per law in each row. Rates and gaps below are rounded separately.

| Slot widths | Matched acceptance | Planted acceptance | Planted minus matched |
|---|---:|---:|---:|
| 8, 8, 16, 16 | 0.49941 | 0.50031 | +0.00090 |
| 4, 4, 4, 4, 8, 8 | 0.49977 | 0.49972 | −0.00005 |

With $n$ independent draws under each law, Hoeffding gives an individual 95% half-width $\sqrt{\ln(40)/n}=0.001921$ for each difference. Both intervals contain zero. This experiment supplies no positive advantage certificate at either arity. Smaller gaps and other readouts are unresolved.
