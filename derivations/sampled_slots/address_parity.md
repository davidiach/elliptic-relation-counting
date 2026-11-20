# Predicting random parities of an isolated signed address

Fix a target $T$ and independently sample the slots. Give each admissible affine point in a slot a canonical binary address consisting of its abscissa address within the dyadic window and a sign bit distinguishing the two ordinates. Reserve only one sign code at a point of order two. A concatenated legal address then denotes exactly one ordered signed tuple. Write $D$ for the full address length. On the isolation event $U=\{H=X=1\}$ the ordinary count and the guarded count are both one. Let $v_*(\sigma)\in\{0,1\}^D$ be that tuple's address at root seed $\sigma$ and let $q_U=\Pr(U)>0$. The mixed affine/fixed-bit sampling law can achieve $q_U^{-1}=2^{O(m)}\operatorname{poly}(m)$ at its thinned widths. [Isolation proof](isolation.md)

## A weak parity signal suffices

A score queried on $\sigma,T$ and a binary challenge $r$ returns a random real number in $[-1,1]$. Its mean over fresh private coins is $f_\sigma(r)$. Define the normalized Walsh coefficient

$$
\widehat f_\sigma(v)=\mathbb E_r f_\sigma(r)(-1)^{r\cdot v}
$$

with uniform $r$ and binary dot product. Suppose a known $\delta>0$ satisfies $\mathbb E[|\widehat f_\sigma(v_*)|\mid U]\ge\delta$. Then an attempt using the score returns a verified guarded tuple with probability at least $3q_U\delta/8$. It uses $\operatorname{poly}(D,\delta^{-1})$ score calls with a deterministic cap on every seed. The sign of the correlation may vary with the seed. This is a sufficient condition on a joint score, and supplies no such score by itself.

Indeed, conditional on $U$, at least a $\delta/2$ fraction of seeds have $|\widehat f_\sigma(v_*)|\ge\tau=\delta/2$, since every coefficient has magnitude at most one. At each seed list the Walsh coefficients of magnitude at least $\tau$ with probability at least $3/4$. The argument is a version of Goldreich–Levin list decoding. [Goldreich and Levin, STOC 1989](https://doi.org/10.1145/73007.73010)

Here are explicit query and failure controls. For a prefix $\alpha$ of length $j$ define its weight $W(\alpha)=\sum_z\widehat f_\sigma(\alpha,z)^2$. Orthogonality over the omitted challenge bits gives

$$
W(\alpha)=
\mathbb E_{u,u',t}\left[f_\sigma(u,t)f_\sigma(u',t)(-1)^{\alpha\cdot(u+u')}\right].
$$

The two independent length-$j$ prefixes share the same uniform challenge suffix $t$. Use independent private evaluator coins on the two calls, including when both challenges coincide. Each sample in this estimator has magnitude at most one. Extend each live coefficient prefix by each bit and keep it if the estimated weight is at least $3\tau^2/4$. Abort the attempt if more than $\lceil 2/\tau^2\rceil$ prefixes remain at any depth. At most $2D\lceil2/\tau^2\rceil$ estimates are attempted, even on a bad seed. Hoeffding's inequality with $O(\tau^{-4}\log(2+D/\tau^2))$ independent samples per estimate makes the chance of any error exceeding $\tau^2/4$ at most $1/4$, conditional on each fixed seed. This uses fresh samples at every adaptive query. On correct estimation all prefixes of a heavy coefficient survive, and every live prefix has weight at least $\tau^2/2$. Parseval gives total weight at most one at every depth, so the cap never aborts in this case.

Decode each listed full address, apply the field-range, ordinate and addition guards, and verify the prescribed sum. A fixed selection policy on verified outputs can be used in both training and online attempts. The success estimate follows without knowing the event $U$ or testing it. The score requires a stable conditional mean on repeated calls. If the model has per-attempt shared randomness, one may include its independently sampled setup in the seed and apply the hypotheses to that larger seed law.

All costs used to construct and query the score must be charged. For example, suppose the independent fresh-seed, uniform-challenge mean call time including per-call setup is at most $C$ for each required target. Truncate a call after $4C/(q_0\delta)$ steps, returning zero on timeout, where $q_0$ is a known uniform lower bound for $q_U$. Changing the score changes its conditional expected absolute coefficient by at most the timeout probability given $U$, at most $\delta/4$ by Markov. The truncated oracle therefore retains at least $3\delta/4$ of absolute correlation and gives a per-call cap on every seed. Use that reduced threshold in the decoder. A common seed-preprocessing phase must likewise be bounded and charged for each fresh setup at which it is needed. With $D=\operatorname{poly}(b)$, $q_0^{-1},\delta^{-1},C\le2^{O(b/m)}$ and all verification and table costs charged, the attempt cost and inverse success are each $2^{O(b/m)}$. To invoke the outer reduction these guarantees must hold in its target and seed laws and with the same selection policy on every attempt.

## What separated scores can see

For a fixed partition into two nonempty blocks $I,J$, suppose $f_\sigma(r)=g_{\sigma_I}(r_I)h_{\sigma_J}(r_J)$, where both real functions have absolute value at most one. They can use $T$ and arbitrary computation on their own seeds. Let $\nu_i=N\rho_i$, $N=|E(\mathbb F_p)|$, $\nu_I=\prod_{i\in I}\nu_i$ and $\lambda=\prod_i\nu_i/N$. Then

$$
\mathbb E\left[|\widehat f_\sigma(v_*)|\mid U\right]
\le \frac{\sqrt{\lambda/N}}{q_U}.
$$

For each group point $y$, let $a_I(y)$ be the expectation over the block seeds of the sum of $|\widehat g_{\sigma_I}(v)|$ over all legal signed block addresses $v$ with group sum $y$. The expected number of such addresses is at most $\nu_I/N$: the last group point is forced by the others and each specified affine tuple is included with probability $\prod_{i\in I}\rho_i$. Cauchy–Schwarz over the block seed and the legal address gives

$$
a_I(y)^2\le (\nu_I/N)
\mathbb E\!\sum_{\substack{v\ {\rm legal}\\\operatorname{sum}(v)=y}}|\widehat g_{\sigma_I}(v)|^2.
$$

The address sets for different $y$ are disjoint. Parseval on the full binary challenge cube yields $\sum_y a_I(y)^2\le\nu_I/N$. The other block satisfies the analogous bound. The joint Walsh coefficient factors as the product of the two block coefficients. Sum their absolute products over **all** legal full tuples of sum $T$ before conditioning. This nonnegative quantity bounds $1_U|\widehat f_\sigma(v_*)|$, and independent block seeds make its expectation $\sum_y a_I(y)a_J(T-y)$. Cauchy–Schwarz bounds this by $\sqrt{\nu_I\nu_J}/N=\sqrt{\lambda/N}$. Dividing by $q_U$ proves the claim.

For a linear combination of such products with predeclared coefficients of total absolute weight $B$, the numerator bound scales by $B$. Different fixed partitions per term are permitted. Independent public coins can be conditioned on before applying the bound. At the isolating widths the result is $B\,2^{-b/2+O(m+\log m)}$, too small for $\delta\ge2^{-O(b/m)}$ when $m\to\infty$, $m=o(b)$ and $B=2^{o(b)}$. The calculation controls the score linearly. A subsequent nonlinear readout, seed-dependent choice using both blocks, or joint arithmetic need not admit the same decomposition.

## Rank and communication across a split

For one fixed partition there is a sharper bound. The block argument also works for unbounded real factors with finite second moment. Parseval now gives $\sum_y a_I(y)^2\le(\nu_I/N)\|g\|_2^2$, where the norm averages over both the block seeds and the uniform challenge bits. Thus a product contributes at most $\sqrt{\lambda/N}\|g\|_2\|h\|_2$ to the unconditional numerator.

View the bounded joint score as a matrix on the two block inputs (seeds and challenge bits together), using their product probability measure. If this matrix has rank $R$, its singular value decomposition in the weighted inner products writes it as a sum of $R$ normalized products with nonnegative weights $\gamma_j$, with $\sum_j\gamma_j^2=\|f\|_2^2\le1$. The triangle inequality and Cauchy–Schwarz give

$$
\mathbb E[1_U|\widehat f_\sigma(v_*)|]
   \le \sqrt{\lambda/N}\sum_j\gamma_j
   \le \sqrt{R\lambda/N}.
$$

In particular a score attaining conditional mean correlation $\delta$ must have $R\ge (N/\lambda)(q_U\delta)^2$. The matrix includes the seeds. A rank bound just on the challenge inputs after fixing all the seeds would allow the unknown witness to be hard-coded.

Give each of two parties its block of seeds and challenge bits and let them compute the score by exchanging at most $c$ bits in a binary protocol tree. Fix their coins. Each transcript leaf is a rectangle in the two inputs, since it imposes a separate consistency condition on each party's messages. If one party returns a bounded score there, absorb that local output into its rectangle factor. There are at most $2^c$ leaves, hence rank at most $2^c$ for these fixed coins. The bound on the numerator survives averaging over independent coins by convexity even when the mean score has larger rank. Consequently

$$
c\ \ge\ \log_2\!\left((N/\lambda)(q_U\delta)^2\right)
    \ \ge\ b-O(m+\log m+b/m)
$$

at the isolating widths when $\delta\ge2^{-O(b/m)}$. This statement charges all communication about the initially separate seeds and challenge bits, including any joint preprocessing. Giving the parties the complete challenge or both seed blocks at no cost changes the hypothesis. [Finite Walsh and isolation checks](../../checks/check_address_parity.py)
