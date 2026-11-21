# Decoding unseen masks from the target quartic

The second determinant moment, supplemented by one scalar already computed by the slot Gram formula, gives a training-free certificate on three-slot pools of shape $(1,1,k)$, $1\le k\le4$. Sizes here count distinct abscissas with both signs retained. The certificate may abstain. It gave correct answers on previously unseen mask triples on three fields. A pair of genuine hash seeds also shows that these data can lose the existence bit, even when the positive input has exactly one witness.

## How to read the six field values

Work on $y^2=f(x)=x^3+ax+c$ over a prime field with $p>5$. Assume disjoint slot abscissas, none equal to the target abscissa $t$, and full pairs of nonzero ordinates. Write $Q(X)=Z_2(X)$ and $G$ for the coefficient of $y(T)^2$ in the Gram formula **before** substituting $y(T)^2=f(X)$. This scalar is the odd diagonal entry of the mixed adjugate in the ordered basis $1,x,x^2,y$. It is additional information beyond $Q$. The slot sizes and regularity are known from setup. The decoder receives only the size shape, $Q,G,t$, and the curve.

For four abscissas $u,v,w,X$, write $K(u,v,w,X)=\sum_i f(x_i)V_i^2$, where $V_i$ is the Vandermonde on the other three. Let $d(u,v,w)=(u-v)^2(u-w)^2(v-w)^2$. The [sign expansion](target_jets.md) gives $q=Q/8=\sum_w K$ and $g=G/8=\sum_w d$ when the first two slots are singletons $u,v$. Set

$$
W(X)=q(X)-2g f(X).
$$

At $X=u$ the first and last summands of $K$ are each $f(u)d$ and the two others vanish. Hence $W(u)=0$, and likewise $W(v)=0$, for any third-slot size. If $W$ is nonzero, its rational roots give at most six unordered candidate anchor pairs. A zero polynomial causes abstention.

For each candidate pair expand $K(u,v,w,X)$ in $w$. Its degree is at most four in each variable. Matching its five $X$-coefficients to $q$, matching $d$ to $g$, and supplying $s_0=k$ give seven linear equations in the five sums $s_j=\sum_w w^j$, $0\le j\le4$. These equations require only interpolation at $w=0,\ldots,4$. In a full-rank system the sums are determined. For a consistent rank-four system parameterize the sums by a single field element. Newton's identities give the monic degree-$k$ support polynomial. When $k<4$, its next power-sum recurrence gives a univariate residual of degree at most $k+1$ in that parameter. The decoder checks all its field roots. It abstains on an identically zero residual or a consistent lower-rank case it cannot resolve.

Reject a candidate if its support fails to split into $k$ distinct admissible abscissas or fails any of the five recovered power sums. For each remaining support, decide the group-sum condition. For different $u,v$ the quadratic

$$
 (u-v)^2 z^2-2((u+v)(uv+a)+2c)z+(uv-a)^2-4c(u+v)
$$

has as its roots the abscissas of the two signed pair sums (up to overall sign). For a candidate $w$, a positive-degree gcd of the quadratics for $(u,v)$ and $(w,t)$ detects the relation. All input abscissas lift to rational points, and the two members of each pair have different abscissas. These conditions exclude an intermediate point at infinity and ensure that matching abscissas allows the signs to be chosen. This is the classical third summation polynomial. [Semaev, 2004](https://eprint.iacr.org/2004/031.pdf)

Return a bit only if all remaining candidates agree. The true anchors and support satisfy every equation, which proves soundness even when the augmented quartic does not uniquely identify the masks. Bounded-degree field factoring and linear algebra take expected $\operatorname{poly}(\log p)$ bit operations after the six inputs are available, at fixed $m=3$. Slot construction and the Gram readout still have to be charged. The actual singleton abscissas could of course be read directly from the original slots. The result here concerns the information retained by the pooled values.

## Fixed-seed holdout against memorization

The draws, labels, and novelty definition are those of the [previous holdout](target_jets.md). Each field used 30,000 plants and 30,000 unlabeled seeds to fit comparison tables, and 20,000 independent fresh draws at $T=(0,1)$ to score them. A mask triple is unseen if its unordered triple of exact masks occurred in neither training sample. The structural decoder receives no training examples. Off its supported shape and on abstentions its stand-alone prediction is absence. The hybrid uses the mask memorizer as fallback. All errors below are on the entire fresh sample, including irregular inputs.

| Field, width | Fresh positives | Mask memorizer errors | Lookup on the same six values, errors | Structural errors | Hybrid errors | Structural true positives among unseen eligible masks | Unseen false positives of structural rule |
|---|---:|---:|---:|---:|---:|---:|---:|
| 23, 2 | 3,165 | 1,975 | 1,976 | 1,955 | 1,752 | 7 / 254 | 0 |
| 101, 4 | 4,763 | 4,744 | 4,744 | 4,143 | 4,138 | 599 / 3,284 | 0 |
| 103, 4 | 4,230 | 4,216 | 4,216 | 3,622 | 3,620 | 587 / 2,914 | 0 |
| 103, 8 | 14,797 | 14,797 | 14,797 | 14,730 | 14,730 | 67 / 6,262 | 0 |

The six-value lookup called no unseen positives in any run. Within the supported shape the unseen positive counts were 7, 599, 587 and 68 respectively. At $p=101$ there were 3,632 unseen eligible draws in shape and the readout certified 3,631 of them. The first two rows reuse previously inspected fixtures. The last two rows test the unchanged decoder on an additional curve chosen after the decoder was fixed, by taking the next prime after 101 and the first positive $a$ for which $y^2=x^3+ax+1$ has non-anomalous prime order. This gives $a=33$, order 97. At doubled window width only 310 of 8,337 unseen eligible draws were in shape. Masks outside the supported shape account for most of the missed hits.

Exhaustive tests on actual hash masks, weighted by ordered seed multiplicity, also used a second target with nonzero abscissa. All certificates agreed with independently computed elliptic group sums. For the third row block, $T=(12,15)$. The other two blocks use $(0,1)$.

| Field, width | $k$ | Eligible seed mass | Certified | Positive mass | Certified positives |
|---|---:|---:|---:|---:|---:|
| 13, 2 | 1 | 146,784 | 146,784 | 75,840 | 75,840 |
| 13, 2 | 2 | 168,480 | 168,480 | 133,008 | 133,008 |
| 23, 2 | 1 | 10,721,424 | 10,658,784 | 3,060,720 | 3,025,728 |
| 23, 2 | 2 | 11,211,984 | 11,202,384 | 5,644,272 | 5,638,128 |
| 23, 4 | 1 | 1,654,080 | 1,648,080 | 466,416 | 466,416 |
| 23, 4 | 2 | 5,308,224 | 5,298,816 | 2,654,592 | 2,650,944 |
| 23, 4 | 3 | 3,328,608 | 3,304,800 | 2,205,552 | 2,193,264 |
| 23, 4 | 4 | 95,040 | 88,512 | 81,456 | 77,232 |

## Exact collision on independently hashed slots

On $y^2=x^3+27x+2$ over $\mathbb F_{37}$, the group has prime order 47. Use the actual hash law $(\alpha x+\beta\bmod64)<\tau$. Two third slots with $(\alpha,\beta,\tau)=(23,9,32),(11,24,32)$ retain respectively the abscissas

$$
 \{3,6,8,9,14,17,22,23,36\},\qquad
 \{4,5,6,16,17,18,22,23,27\}.
$$

Over the signed points both have power sums $(18,17,11,11,28)$ in degrees zero through four, and both have $\sum y^2=25$. Consequently they have identical entire unweighted slot Gram matrices. With the other slots fixed they give the same $Q,G$ and **all** target derivatives at **every** target.

* With first slots $(1,6,8),(7,54,8)$, the singleton abscissas are 1 and 2. At $T=(30,5)$ the respective numbers of witnesses are **1 and 0**.
* With first slots $(49,37,8),(45,21,8)$, their abscissas are $\{2,32\},\{1,21\}$. At $T=(13,16)$ the respective counts are **0 and 5**, verified by group addition and by counting literal zero determinants.

All four triples are regular at the quoted targets. This proves failure of exact universal decoding from the six values and shape. It gives no asymptotic average-case error lower bound by itself. The widths in these examples are $(8,8,32)$. Weighting just the third slot by $x$ distinguishes the second pair, with $G_{\rm weighted}=30,27$ respectively. That observation identifies higher slot information erased by the original Gram statistic. Its ability to decode typical fresh inputs remains open.

## What the root test measures at general arity

Write $C=\sum_{S\subseteq[m]}(-1)^{m-|S|}\operatorname{adj}(\sum_{i\in S}A_i)$ for the [mixed Gram matrix](target_incidence.tex). Give each signed point in each slot its own formal weight $h_{i,P}$. For $T=(X,Y)$, polarization of the mixed Cauchy–Binet identity gives

$$
 B_h(X)=v(T)^\mathsf{t} C v(-T)
   =\sum_{\mathbf P}\left(\prod_i h_{i,P_i}\right)
       \Delta_{\mathbf P}(T)\Delta_{\mathbf P}(-T).
$$

Symmetry of $C$ cancels the mixed even/odd terms. Substituting $Y^2=f(X)$ gives a polynomial of degree at most $m+1$, even with different weights on opposite signs. For $m=3$ and unit weights this polynomial is $Q-2Gf=8W$.

**Claim.** Suppose $m\ge2$, the rational group has odd prime order, and the slots are nonempty, sign-paired and pairwise disjoint in abscissa. At an affine rational target $T$, $B_h(x(T))$ is identically zero as a polynomial in the separate point weights exactly when some slot has the single abscissa $x(T)$.

If a slot forces that abscissa, every tuple repeats either $T$ or $-T$ in one of the two determinants. Conversely, when no slot forces it, choose one other abscissa in each slot. Fix all signs except those of two points $P,R$ from distinct slots. The four possible sums are different: doubling is injective in odd order and $P\ne\pm R$. At least one sum avoids both $T$ and $-T$. All the selected abscissas avoid $x(T)$ and are pairwise distinct, so the elliptic determinant incidence criterion makes both determinants nonzero for that tuple. Its formal weight monomial is distinct from every other tuple's, giving a nonzero coefficient.

With independent uniform field weights, the probability of accidentally vanishing at any prescribed non-forced rational abscissa is at most $m/p$ by the polynomial zero bound. With two independent sets of weights the probability of any spurious common rational curve abscissa is at most $m^2/p$ by a union bound over at most $p$ abscissas. In the growing-arity regime the common roots can be read via a bounded-by-$m+1$ degree gcd and field factoring if at least one polynomial is nonzero. A simultaneous identically zero output can be flagged. The computation uses two existing Gram accumulations at the same exponential scale. The direct check compared the actual polynomials with signed determinant products at all rational abscissas on seven fixed fixtures spanning $m=2,\ldots,6$.

## End-to-end consequence

At the balanced tuning the expected signed size of each slot is $\mu_i=(N-1)\rho_i=2^{b/m+O(1+(\log L)/m)}$ with $L=\log_2 b$. The project's variance estimate gives $\operatorname{Var}|F_i|\le2\mu_i$. For $\min_i\mu_i>16$, Chebyshev and a union bound give

$$
 \Pr(\exists i:\ |F_i|\le8)\le \sum_i {2\mu_i\over(\mu_i-8)^2}
       =O(m/\min_i\mu_i)=o(L^{-2}).
$$

A predictor that calls presence only on such sparse shapes can recover at most that much fresh probability mass. Its false-negative mass is at least $\lambda-O(\lambda^2)$, using the project's targetwise hit estimate. The required total error is $O(\lambda^2)$. The general cross-target root test detects forced abscissas, so it cannot provide the missing existential test on typical regular inputs. With the same packet and reduction assumptions the best complete log-time upper bound remains $(26.2933775156+o(1))b/L^{1/4}$.

From the project root, run the following with Python, NumPy and python-flint.

~~~sh
python -B checks/check_quartic_decoder.py
python -B checks/check_opposite_targets.py
~~~
