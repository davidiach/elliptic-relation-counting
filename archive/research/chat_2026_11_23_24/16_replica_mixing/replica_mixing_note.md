# Uniform replica mixing and Poisson relation counts for translated elliptic slots

Research checkpoint, 24 November 2026. Continues `affine_seed_spectrum_note.md` (especially its higher-order incidence section), the earlier affine-digit sampler, and v3 of *Counting relations in independently hashed elliptic-curve slots*. **No strict general-prime ECDLP exponent improvement is established.** The previous candidate coefficient 11.31486814 is not re-audited. This note proves a stronger sampling/statistics result, rather than a faster per-seed relation evaluator. No claim of priority or independent peer review is made.

## 1. Model and theorem

Let `p>3` be prime, let `E/F_p: y^2=f(x)=x^3+a x+c` be smooth, and suppose `G=E(F_p)` has odd prime order `r`. Set `b=ceil(log_2 p)`. Take any fixed sets `S_i ⊂ F_p`, with sizes `w_i`, and any fixed nonzero multipliers `t_i`. Choose translations `a_i` independently uniformly in `F_p`. The slot `F_i` retains all rational affine points with `x(P) ∈ a_i+t_i S_i`, including both signs. Thus this result applies conditional on *every* choice of the multipliers of the previous field-affine sampler. No random multipliers are necessary. This is not the original paper's mod-2^b odd-affine hash model.

Write `W=∏_i w_i`, `lambda=W/r`, and `H_T= #{(P_i)∈∏F_i : ∑ P_i=T}` for the ordinary ordered signed relation count. The convention `W/r` differs slightly from the paper's `r^(m-1)∏(w_i/p)`; the ratio is `(r/p)^m=1+O(m/sqrt(p))` in the range considered.

**Theorem.** Assume `m→∞`, `m≤sqrt(b)`, `w_i=2^((1+o(1))b/m)` uniformly in i, and `b^(-A)≤lambda≤C` for fixed positive constants A,C. For any positive integer `k` up to

    K=floor(log(m)/(4 log log(m)))

(for large m; logs in this expression have any fixed common base), let `T_1,...,T_k` be nonidentity targets, possibly equal, but with `T_j != -T_l` for every distinct pair of indices. If the multiplicity of target U in this list is `k_U`, then, uniformly in targets, curve, supports and multipliers,

    E ∏_U (H_U)_(k_U) = lambda^k + O(2^(-b/(2m))).        (1)

Here `(h)_j=h(h-1)...(h-j+1)` is the falling factorial. The implied bound holds eventually along the stated width and arity regime. The constants in the explicit error bound below can be used instead of this asymptotic form. For fixed k the same theorem is understood once K≥k. In particular at polynomially sparse density the relative error for any fixed factorial moment is o(1).

If `lambda→lambda_*>0`, then `H_T` converges to Poisson(lambda_*) uniformly over targets and all allowed samplers. Counts at any fixed number of distinct, pairwise nonopposite targets converge jointly to independent Poisson(lambda_*) variables. Opposite targets cannot be included in the independence statement: the sign symmetry gives `H_T=H_{-T}` exactly.

In particular, for a fixed draw of the slots, the fraction of nonidentity targets with a relation converges in probability to `1-exp(-lambda_*)`. This is a statement about target coverage using the *same factor bases*, not just a mean over fresh slots for one target. The theorem also gives the sparse expansion `Pr(H_T=1)=lambda-lambda^2+O(lambda^3+2^(-b/(2m)))` when lambda→0. See §5.

## 2. A uniform one-slot replica-character estimate

The higher-order incidence formula in the preceding checkpoint counts, for a prescribed collection of abscissas, how many affine masks contain them. To sum the terms arising in factorial moments collectively, interchange that sum with the choice of ordered distinct offsets in each mask. The following estimate is the load-bearing result.

For any fixed ordered d distinct field offsets `s_1,...,s_d`, a fixed `t!=0`, and a uniform translation `a`, form the finite measure on `G^d` that counts every tuple of rational affine lifts of the abscissas `a+t s_j`, with weight 1/p. Let its mass be `u_d`. For a product of group characters `(chi_1,...,chi_d)`, its Fourier coefficient is the corresponding weighted sum of `∏ chi_j(P_j)`. We show, for `d≤k`,

    |u_d-1| ≤ C_k/sqrt(p),
    |Fourier coefficient| ≤ C_k/sqrt(p)  if some chi_j is nontrivial,
    C_k = 16 k 2^k.                                      (2)

A more precise bound for the nontrivial coefficient is `((3d-3)2^(d-1) sqrt(p)+2^d)/p`.

### 2.1 All branch loci are disjoint, for every nonzero slope

Because the group order is odd, `f` has no F_p root (a root gives a rational point of order two). It is therefore an irreducible cubic. Its three roots in F_(p^3) are cyclically permuted by Frobenius. No two distinct roots can differ by a nonzero element `h∈F_p`: if `beta=alpha+h=Frob^j(alpha)` with j=1 or 2, then iterating Frob^j three times gives `alpha=alpha+3h`, impossible in characteristic >3. Hence the roots of `f(A+t s_j)` are disjoint across j, since an overlap would give such a root difference. This is stronger than the generic bound excluding at most `9 binom(d,2)` exceptional multipliers. It explains why *all* multipliers can be fixed in the prime-order case.

### 2.2 The fiber-product curve

Over the algebraic closure consider the function field

    k(A)(sqrt(f(A+t s_1)),...,sqrt(f(A+t s_d))).

The d square classes are independent: a root of the j-th cubic has odd valuation in that factor and zero valuation in all the others. Thus the smooth projective normalization C of this fiber product is geometrically connected and has a Galois map to the A-line of degree `2^d` with group `(Z/2)^d`. At the 3d finite branch points the inertia has order two; at infinity it also has order two because the cubics have the same odd pole order and their ratios are squares in the completed algebraically closed local field. No other point ramifies. The tame Riemann–Hurwitz formula gives

    2g(C)-2 = 2^(d-1)(3d-3).                             (3)

The projections `pi_j:C→E` given by `(A+t s_j,y_j)` have degree `2^(d-1)` and are morphisms on the smooth projective models. Their coordinate sign changes are independent deck transformations on C. The affine fiber product is nonsingular: at each branch point at most one square-root coordinate vanishes, and the corresponding cubic has a simple root. Hence the finite rational tuples above are exactly the rational points of C away from the fiber at infinity. There are at most `2^d` omitted rational points. Hasse–Weil on C proves the mass part of (2), including d=1.

### 2.3 Geometric nontriviality of the character

For each group character `chi_j`, use its unramified rank-one Artin character/local system on E supplied by the finite étale Lang covering `Frob_p-1:E→E`, with kernel E(F_p). A nontrivial character has order r and is geometrically nontrivial, as the cover has geometrically connected domain and degree r. Its trace on a rational point is the group character (choosing the inverse convention if necessary). This construction also works when r=p, since the Lang map is étale.

The product character on C is the tensor product of the `pi_j^*chi_j`. Suppose it were geometrically trivial and some `chi_j` nontrivial. Pulling back by the automorphism which flips only y_j leaves every other factor invariant and inverts the j-th factor. Comparing the two characters forces `(pi_j^*chi_j)^2` to be geometrically trivial. But the pullback of a connected degree-r cyclic covering through a map of degree `2^(d-1)` remains connected: its field extension cannot intersect the degree-power-of-two extension nontrivially, since r is an odd prime. Thus the pullback character has order r, and its square is nontrivial. Contradiction.

Apply the classical Weil bound for a geometrically nontrivial unramified rank-one character on a smooth projective genus-g curve: the absolute character sum over rational points is at most `(2g-2)sqrt(p)`. Removing the points above infinity costs at most `2^d`. Equations (2) and (3) follow. The result uses the standard character-sum theorem, not an unproved estimate on growing-order elliptic characters. For the general Lang character-sheaf construction see Sawin and Templier, *On the Ramanujan conjecture for automorphic forms over function fields I*, Lemma 2.11, https://preprints.ihes.fr/storage/OnRamanujan.pdf . For the curve character bound see the character-sum section and Lemma 1 at https://www.normalesup.org/~tibouchi/papers/welldistributed.pdf . Only the unramified, geometrically nontrivial case of that bound is used here.

### 2.4 Average over offsets without paying for their enumeration

For slot i let `mu_{i,d}` be the finite measure of ordered d retained signed points having pairwise distinct abscissas, averaged over its uniform translation. Each such tuple determines its d distinct offsets uniquely for any realized seed. Therefore `mu_{i,d}` is the sum of the preceding measures over the `(w_i)_d` ordered offset choices. We do NOT need to construct these choices to establish the estimates, because (2) holds uniformly for all of them.

Put `a0=C_k/sqrt(p)` and `eta=2a0`; suppose a0≤1/2. We have

    mass(mu_{i,d}) = (w_i)_d(1+theta), |theta|≤a0,
    |normalized Fourier coefficient| ≤ eta  at all nontrivial characters.  (4)

The normalized measure is a probability measure. These are estimates on the true signed abscissa masks, not a replacement with independent points. In particular the previous higher-order interval incidence formula is not needed at this step: the resulting collective estimates are support-shape independent.

## 3. Signed collision patterns and an exact factorial-moment decomposition

For k replica tuples, each slot contributes a partition of replica indices according to which points have equal abscissas. Within a block the points are the representative point or its negative. Choose the lowest-index replica as representative with sign +; all other signs are then uniquely determined because the group has odd order. A signed partition P has d(P) blocks. Encode it by a k×d matrix A_P whose columns have disjoint supports and entries ±1 in their own blocks, with exactly one nonzero per row. Its image V_P in `G^k ≅ F_r^k` is a d-dimensional linear subspace. The isomorphism is only used in the proof and does not require computing any elliptic discrete logarithms.

There are at most `B=(2k)^k` signed patterns (the number is `∑_d {k brace d} 2^(k-d)`, with Stirling numbers of the second kind). Each type's expected count measure on `G^k` is precisely the pushforward of `mu_{i,d}` through A_P. In particular its mass is as in (4), its support is contained in V_P, and the normalized Fourier coefficient of a character is 1 if it annihilates V_P, of modulus ≤eta otherwise.

For a sequence of types `P_1,...,P_m`, define it admissible if no two replica indices have *identical* rows (same block and sign) in every slot. This is precisely the condition that the resulting full signed tuples are pairwise distinct. A repeated full tuple can only have a repeated target. The desired mixed factorial moment is therefore the sum, over the admissible type sequences, of the convolutions of the expected slot measures at the target vector `T=(T_1,...,T_k)`. This equality is exact: slot measures are independent and every tuple of replicas has one unique type sequence. It applies at arbitrary finite widths and arities, including opposite targets (though the asymptotic theorem excludes those targets).

A useful elementary structural fact: `V_P` determines P uniquely. Indeed the coordinate functionals restricted to V_P are proportional up to a sign precisely for indices in the same block, and their relative signs determine the normalized signed pattern. In particular, if a subspace V has dimension R then at most one type P with `V_P⊆V` has d(P)=R; all others have d(P)≤R-1.

If `R<k` and the maximal type (when it exists) occurs in every slot, at least two replica indices in some block have a fixed relative sign. If the target vector has nonzero probability, the two target points must have that relative sign. Opposite targets were excluded. Thus those replicas are identical full tuples and this sequence cannot contribute to the mixed factorial moment. This is the key step suppressing the potentially large diagonal terms.

## 4. Collective estimate with rare patterns separated

Set `h=2k+4`, `J=(h-1)B`, and suppose m>J. Call a type frequent if it occurs at least h times among the m slots. Let V be the span of the frequent-type spaces, dimension R≥1. Every frequent slot lies inside V, and at most J slots are outside V (meaning `V_{P_i}` is not a subspace of V). Include all inside slots in their convolution. For any nontrivial character of V, at least one frequent type has nontrivial restriction, so the product of inside-slot Fourier coefficients has modulus ≤eta^h. Fourier inversion on V then yields pointwise probabilities

    r^(-R)(1 ± delta),   delta = r^k eta^h.               (5)

Here and below an upper bound with delta remains valid if its lower counterpart is negative. For the all-singleton-type sequence a sharper error is `r^k eta^m`. The mass factors of the unnormalized measures are bounded above by `(1+a0)^m∏ w_i^(d_i)`.

Write `v=min_i w_i`, `U=max_i w_i` and `A_lambda=max_{1≤R≤k} lambda^R`.

### Case A: R=k

All slots are inside. The unique all-singleton type sequence (d_i=k for every i) contributes

    lambda^k [∏_i ((w_i)_k/w_i^k)] (1±a0)^m (1±r^k eta^m).

The notation means a value in the product of the indicated intervals. The error from lambda^k is `O(lambda^k (m k^2/v + ma0 + r^k eta^m))` when these three quantities tend to zero. Every other type sequence has at least one deficit `k-d_i≥1`. Summing over ALL possible type sequences with a deficit (an upper bound), the total contribution of this case apart from the main sequence is at most

    (1+delta)(1+a0)^m lambda^k ((1+B/v)^m-1).             (6)

### Case B: R<k and no outside slots

There are at most `2^B` possible V, each a span of a subset of the type spaces. For a fixed V of dimension R, by §3 at most one inside type has d_i=R, and an all-maximal sequence is inadmissible or impossible for the target vector. Hence a contributing sequence has at least one deficit `R-d_i≥1`. The sum over this case is at most

    2^B (1+delta)(1+a0)^m A_lambda ((1+B/v)^m-1).          (7)

This bound still holds if no maximal type exists; adding a hypothetical maximal type only increases the estimate on sequences with a deficit. Importantly, the factor `B^m` is NOT charged to sequences with just one deficit. That would lose the useful asymptotic estimate at large m.

### Case C: R<k and some outside slot

Pick one outside slot j. There is a linear functional of `F_r^k` vanishing on V and nonzero on `V_{P_j}`. Every nonzero scalar multiple of it induces a nontrivial character on that slot. Equation (4) and one-dimensional Fourier inversion give maximum point probability `1/r+eta` for the projected variable. Condition on the other outside slots. For the sum of the outside values to land in the desired coset of V, at least this one projected equality must hold. Given the outside values, the inside sum has point mass at most `(1+delta)/r^R` by (5). Thus each sequence's normalized target probability is at most `(1+delta)r^(-R)(1/r+eta)`.

There are at most J outside slots, each with d_i≤k; inside slots have d_i≤R. Consequently `∏_i w_i^d_i ≤ W^R U^(kJ)`. Summing even over all `B^m` type sequences gives the upper bound

    (1+delta)(1+a0)^m A_lambda B^m U^(kJ) (1/r+eta).      (8)

These three cases exhaust all admissible sequences. Equations (6)-(8), together with the stated all-singleton error, are explicit quantitative bounds for the factorial moment. They allow arbitrary fixed multipliers and different support shapes. No factorial-moment calculation enumerating the elliptic relation tuples was assumed.

### Asymptotic accounting

For `k≤log(m)/(4 log log(m))`, eventually `B≤m^(1/3)`, `J≤m^(1/2)` and `kJ=o(m)`. In the given regime, `log_2 v` and `log_2 U` are `(1+o(1))b/m`, Hasse gives `log_2 r=b+O(1)`, and `ma0=o(1)`. Also `delta≤p^(-2+o(1))` because h=2k+4 and `log C_k=O(k+log k)`. The all-singleton term's `r^k eta^m` is negligible as m≫k.

For (6) and (7), `(1+B/v)^m-1=O(mB/v)`; `2^B=2^o(b/m)` since `m≤sqrt(b)`. Moreover `A_lambda≤max(1,C^k)=2^o(b/m)`. Thus these terms and the main-term error are `2^(-b/m+o(b/m))` (the main-term `ma0` is even smaller). For (8), `m log_2 B=o(b)` and `kJ log_2 U=o(b)`, while `1/r+eta=p^(-1/2+o(1))`. Hence (8) is `2^(-b/2+o(b))`. Taking a safe exponent of `-b/(2m)` gives (1) uniformly for all k in the stated growing range.

## 5. Consequences and quantitative probability formulas

For a fixed limit `lambda_*`, convergence of all fixed mixed factorial moments implies joint Poisson convergence. One can alternatively use the finite Bonferroni inequalities and tightness from the first moments. Uniformity follows by applying the same bounds to every possible sequence of targets, supports and fixed slopes. For fixed k and lambda≥b^(-A), the error in (1) is negligible relative to lambda^k.

For example at density lambda=1+o(1), the limiting probabilities are

    Pr(H_T=0) → e^(-1) = 0.3678794412,
    Pr(H_T=1) → e^(-1) = 0.3678794412,
    Pr(H_T≥2) → 1-2/e = 0.2642411177,
    Pr(H_T>0) → 1-1/e = 0.6321205588.

These improve the earlier lower bound 1/2-o(1) for coverage at density one, and give the whole limiting distribution. They are asymptotic probabilities, not empirical cryptographic-scale measurements.

For sparse lambda→0, the first three factorial moments and ordinary inclusion-exclusion give

    Pr(H_T>0) = lambda-lambda^2/2+O(lambda^3+e_b),
    Pr(H_T=1) = lambda-lambda^2+O(lambda^3+e_b),
    Pr(H_T=2) = lambda^2/2+O(lambda^3+e_b),
    Pr(H_T≥3) = O(lambda^3+e_b),
    e_b=2^(-b/(2m)).

For example use the inequalities H-(H)_2 ≤ 1_(H=1) ≤ H-(H)_2+(H)_3/2 for integer H≥0. The other formulas follow from the analogous alternating binomial inequalities. No unproved assumption that the relation count is Poisson was made in these expansions.

For the same-base coverage assertion, let `Z=(r-1)^(-1)∑_{T≠O}1_(H_T>0)`. The expectation tends to `1-e^(-lambda_*)`. Joint convergence for pairs T,U distinct and nonopposite makes their covariances tend uniformly to zero. The O(r) equal/opposite pairs contribute O(1/r) to the variance. Therefore Var Z=o(1). More quantitatively, one can use the first K factorial moments of H_T, and of H_T+H_U (the binomial expansion of its factorial moments costs at most `2^K e_b` error), to bound the zero-count probabilities by truncation of the exponential series with error `O((C'/K)^K+2^K e_b)`, for constant C' depending on the density bound. This also yields uniform coverage concentration. The statement concerns random translations with fixed support and slope choices; it does not concern arbitrary adversarial fixed seeds.

Regularity (distinct abscissas across slots, avoiding x(T)) can be enforced at cost of an o(1) probability loss: each abscissa has uniform marginal w_i/p under translation alone and different slots are independent, so the failure bound is at most `sum_(i<j) w_i w_j/p + sum_i w_i/p`. The ordinary and guarded counts also coincide with probability 1-o(1) at the stated densities. For a fixed target, each bad equal-abscissa or identity partial-sum event is a nontrivial linear constraint on group coordinates in addition to the total-sum equation (or impossible for T!=O), so there are O(m r^(m-2)) bad tuples before sampling. With affine marginal `∏(w_i/p)` the expected bad count is O(m lambda/r) in this regime. Markov's inequality suffices. These arguments are specifically for the prime-order group and nonidentity target; composite-order extension and the original modulo-power-of-two sampler are not proved here.

## 6. What changes computationally, and what does not

- The full low-order factorial moments can be approximated collectively to the error in (1) by the explicit numbers `lambda^k`, with no tuple enumeration. This is an **ensemble statistic**, averaged over the random translations, not the count of a given input seed.
- Multipliers may be fixed to 1, even in every slot. The sampler therefore has only one random field word per slot. The fixed translation sampler has stronger coverage than the earlier average-target-only guarantee, on the prime-order subclass. At the prior compiler widths one can in principle choose whatever public support fits the arithmetic representation most conveniently; deterministic slopes do not harm the above coverage. Because those slopes were already public constants before contraction, eliminating their randomness by itself does not lower the leading wire count or the asymptotic relation-evaluation cost.
- Higher-order incidence probabilities can depend strongly on support shape (as shown in the preceding note). Nevertheless, the limiting relation-count distribution does not, under the hypotheses here. Searching for a support with much better constant-density relation coverage cannot yield an asymptotic improvement through the coverage statistic alone. This does NOT prove equality of nonlinear computational complexity, or of per-input likelihoods.
- For density one, inverse success of an exact verified finder improves at most by a constant using this coverage proof. Even removing retry losses entirely does not change the inherited balance of `b/m` against `m b/sqrt(log b)`.

The remaining task is evaluation conditional on the actual given translations. It may use the low entropy and deterministic common multipliers, but the collective Fourier argument here integrates the translations out, removing exactly the input-specific information a finder needs. The desired `o(b/(log b)^(1/4))` general-prime exponent, including polynomial-cofactor curves, remains unproved. Neither a new compiler ledger nor a new stronger contraction theorem is claimed. The previous 11.31486814 coefficient remains a conditional candidate, not a newly validated baseline. No practical discrete logarithm computation was attempted.

## 7. Reproduction and verification

Run `OPENBLAS_NUM_THREADS=1 python -B verify.py` with Python 3 and NumPy. The verifier:

1. Generates all signed partitions of sizes 2 through 6 and checks distinct image subspaces using exact modular RREF, plus the duplicate-pair condition on compatible nonopposite target vectors.
2. Constructs expected finite slot arrays independently by direct replica-tuple classification and by distinct-representative expansion. On four small actual curve fixtures, compares the complete signed-pattern convolution with literal seed-class enumeration and integer falling factorial counts. The convolution uses exact NTTs in auxiliary prime fields, and checks that the product of moduli exceeds an a priori upper bound for the integer counts. Thus the resulting equality checks are exact, not floating point.
3. For fixed offsets of sizes 1 through 4, checks shifted cubic coprimality by polynomial gcd, the mass bound, and the nontrivial character sums for every slope by numerical FFT on the small group arrays. These FFT comparisons are numerical regressions, not exact certificates; the geometric proof above is needed for the asymptotics. As predicted, every tested slope is good. These toy bounds can be loose and the tests are not large-scale Poisson experiments.

The successful transcript is in `verification.txt`. The tests cover the finite signed-pattern and group-character identities, not the complete contractor or ECDLP runtime. The growing-arity results rely on the classical Weil bound, the elementary geometric nontriviality argument, and the quantitative pattern summation proof above.
