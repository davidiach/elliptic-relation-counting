# Adaptive nonlinear equality tests on signed moments

A bounded score may use all public seeds and branch on joint polynomial identities among any of the signed moments. If the total pole degree of the equalities it actually visits on the generic target is small, its correlation with the correct address is exponentially small on the masked-target distribution. The same argument bounds a procedure that directly returns a short list of addresses. This strengthens the earlier target-averaged existence bound to address recovery, adaptive branching, and limited arbitrary target advice. [Earlier target-incidence argument](../joint_parity/target_incidence.tex)

## Target dependence

Let $\mathcal T$ be a set of $n$ distinct nonidentity targets, sampled uniformly and independently of the slot seeds $\sigma$. The outer reduction uses the nonidentity elements of its search subgroup. All addresses and their curve points are fixed before drawing $T$. Write $U=U(\sigma,T)$ for the event that exactly one signed tuple has sum $T$ and passes the guards, with no further ordinary tuple of that sum. Let $q=\Pr(U)>0$ and $v_*=v_*(\sigma,T)$ be the unique address on $U$. The challenge $r$ is uniform on the binary address cube. For every fixed seed, two different targets having unique witnesses have different witness addresses because each tuple has a single sum. [Address convention and sparse success](../sampled_slots/isolation.md)

At a fixed seed and challenge, the $j$-th signed moment $Z_j(T)$ is a rational function of the endpoint with its only possible pole at $O$, of order at most $jd$, where $d=m+1$. Expand each evaluation determinant in its endpoint row, whose entries belong to $L(dO)$, then raise to the $j$-th power and sum with the target-independent challenge weights. The same statement holds for any target-independent point weights chosen from the seeds and challenge, and for different weights at different nodes of a decision tree.

Consider a rule whose nodes test whether a polynomial in signed moments is zero. Coefficients, weights, and the choice of the next test may depend on all the seeds, the challenge, independent randomness and the previous comparison bits. A polynomial with monomials $\prod_j Z_j^{e_j}$ has pole order at most $d\max\sum_j je_j$ if its coefficients are independent of the target. Direct target-coordinate polynomials with accounted pole degree are allowed. A rational comparison can first branch on each required denominator and then test its cleared numerator with the respective degree bounds. A leaf returns any score in $[-1,1]$ determined by the data and transcript used so far. Extra arbitrary information about $T$ is allowed as a label $a_\sigma(T)$ with at most $J$ values, supplied before the challenge is drawn. An unrestricted field lookup indexed by a moment value or a seed-and-challenge-dependent target bit is outside this access model.

For each fixed seed, challenge, coins and possible label, follow the tree formally in the function field of the curve. Take the equality branch precisely when the tested function vanishes identically. This defines a **generic transcript** for that label without selecting any particular target. Assume these transcripts terminate. Let $B$ bound, or upper-bound in expectation over seeds, challenge and coins, the sum of pole degrees for the nonidentically-zero tests on these paths **over all $J$ labels**. The tree can be arbitrarily large away from those paths.

## Advantage and list-size bounds

Let $f_{\sigma,T}(r)$ be the conditional mean score. With the normalized Walsh transform,

$$
\mathbb E_{\sigma,T}\!\left[1_U
       |\widehat f_{\sigma,T}(v_*)|\right]
 \ \le\ \sqrt{Jq/n}+2B/n .
$$

To prove it, fix the random tape. For a fixed seed and challenge, every target with a given label follows its generic path unless it is a zero of a nonzero function tested on that path. A nonzero rational function on a smooth complete curve has as many zeros as poles, counting multiplicity. Thus across all labels the number of exceptional targets is at most the charged sum of degrees. This counts the **first** departure; later tests on the exceptional path incur no charge in this argument. The two scores differ by at most two. Averaging over uniform $T$ and over $r$ bounds their contribution to the expected absolute Walsh coefficient by $2B/n$. [Principal divisors have degree zero, Sutherland, Lecture 23](https://math.mit.edu/classes/18.783/2022/LectureNotes23.pdf)

For fixed $\sigma$, let $h$ be its number of singleton targets and $g_a(r)$ its generic score for label $a$. The sum of the squares of $|\widehat g_{a_\sigma(T)}(v_*)|$ over singleton targets is at most $J$. Indeed, within each label the witness frequencies are distinct, and each bounded $g_a$ has Parseval energy at most one. Cauchy–Schwarz therefore bounds the sum of magnitudes by $\sqrt{Jh}$. Take expectations and use $\mathbb E h=nq$. Averaging over the independent random tape preserves the result by convexity. Freezing the tape means using one coupled realization of the score procedure across inputs; it imposes no constraint on how calls may use fresh coins.

For a procedure that returns a list of at most $L$ full addresses at a leaf, the same generic-path argument gives

$$
\Pr(\text{the returned list contains a target tuple}) \le (JL+B)/n .
$$

For a fixed seed, challenge and tape the $J$ generic lists contain addresses of at most $JL$ distinct targets. All additional successes require an exceptional target. This version needs no uniqueness or guard assumption, and allows exact external verification of the returned candidates. The bound also applies to an adaptive multi-query procedure if its whole transcript fits the tree model and its whole generic path is charged.

On a large subgroup with $\log_2 n=b-O(1)$ and the isolating law with $q=2^{-O(m+\log m)}$, $J,B=2^{o(b)}$ imply conditional mean magnitude at most $2^{-b/2+o(b)}$ for growing $m=o(b)$. On this target class, attaining conditional magnitude $\delta$ requires $B\ge (nq/2)(\delta-\sqrt{J/(nq)})$ whenever the right side is positive. This is an average with respect to the actual, success-weighted fresh target. A uniform targetwise signal would also have to satisfy it. Repeated squaring can create enormous pole degree with few operations, so this degree bound does not imply an arithmetic circuit time lower bound.

For example, a Hankel rank test or a decision tree testing products, sums, and ratios of low orders has an immediate degree charge even with challenge-specific weights. For $Z_0Z_2-Z_1^2$ it is $2d$, and for $Z_1Z_3-Z_2^2$ it is $4d$. Tests against a short table of prescribed values can charge the number of entries times the comparison degree, or test the product of the corresponding differences. The field bits, multiplicative characters of high-degree expressions, and arbitrary algorithms for processing the full moment vector remain open. [Finite check](../../checks/check_nonlinear_branches.py)
