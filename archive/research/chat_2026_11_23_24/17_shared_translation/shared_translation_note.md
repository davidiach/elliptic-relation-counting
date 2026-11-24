# Correlated translation samplers, cofactor transfer, and inverse relation queries

Research checkpoint, 24 November 2026. Continues `replica_mixing_note.md`, `affine_seed_spectrum_note.md`, and version 3 of *Counting relations in independently hashed elliptic-curve slots*.

**Status.** This checkpoint proves a sampler theorem that permits many slots to share the same translation, covers polynomial-cofactor curves, and gives a constructive inverse-translation query in expected `2^{O(m)} poly(b)` bit time. It does **not** establish a faster forward relation evaluator or beat `b/(log b)^(1/4)`. The previous candidate coefficient 11.31486814 is not independently re-audited. The proof invokes the classical Weil character-sum theorem and standard elliptic-curve degree identities. No claim of research priority or external peer review is made.

## 1. Model and main results

Let `p>3` be prime and `E/F_p: y²=f(x)=x³+a x+c` smooth. Suppose `N=#E(F_p)=h r`, where `r` is an odd prime, `gcd(h,r)=1`, and `h=poly(b)` for `b=ceil(log2 p)`. Thus `r=2^{b-O(log b)}`. Let `G_r=[h]E(F_p)`, of order r, and use `π(P)=[h]P` (not an idempotent projection, but a surjection). This includes h=1 and does not assume r differs from p.

Partition m slots into q nonempty balanced blocks `B_1,...,B_q`, their sizes differing by at most one. Use fixed finite offset sets `S_i⊂F_p`, `w_i=|S_i|`. Within each block the sets must be pairwise disjoint, and the shifted cubics `f(A+s)` for every s in their union must be pairwise coprime over the algebraic closure. Across distinct blocks there is no disjointness condition. A construction of these public supports is given in §2. Choose just q independent uniform translations `A_j∈F_p`. In slot i belonging to B_j retain both signs of every point with `x(P)∈A_j+S_i`, excluding zero ordinates. Removing zero ordinates makes signed-pattern encoding unambiguous; it costs a negligible probability at the prescribed widths.

Let `W=∏ w_i`, `λ=W/r`, and

    H_T = #{(P_1,...,P_m)∈∏ F_i : π(∑ P_i)=T}, T∈G_r.

This is a count of raw ordered signed tuples, not of distinct projected tuples. Assume `m→∞`, `m≤sqrt(b)`, `log2 w_i=(1+o(1))b/m` uniformly, and `b^(-A)≤λ≤C` for fixed positive A,C.

**Theorem A (one translation).** For q=1 and a uniform independent nonidentity target T∈G_r, the probability `H_T>0` is at least `2^{-O(m)} / poly(b,λ^{-1})`. The proof is independent of the prior independent-slot sampler. If the supplied query algorithm can find a witness whenever one exists (after the negligible regularity checks below), the one-pilot outer reduction of the preceding large-arity note applies. The single shared random seed has b bits, rather than mb bits. One may use q=1 even on general primes in the polynomial-cofactor class.

**Theorem B (five translations).** For q=5, uniformly for nonidentity T,U∈G_r with U≠-T,

    E H_T = λ + O(e_b),
    E (H_T)_2 = λ² + O(e_b),
    E H_T H_U = λ² + O(e_b) if U≠T,
    e_b = 2^{-b/(2m)}.

In particular at λ→1 the targetwise positive probability is at least `1/2-o(1)`, and at polynomially sparse λ→0 the singleton probability is `λ+O(λ²+e_b)`. This is a second-moment statement; five translations alone are not claimed to give the full Poisson law.

**Theorem C (all fixed moments with very few translations).** For a moment of order k, set `s=2k+4`, `J=k(s-1)` and take `16kJ≤q≤m`. For any k distinct replica positions carrying nonidentity targets `T_1,...,T_k`, with no opposite-target pair, denote by k_U the multiplicity of target U. Then

    E ∏_U (H_U)_(k_U) = λ^k + O(e_b).

The bound holds uniformly for `1≤k≤floor(log(m)/(4 log log(m)))` with the stated q condition. Thus `q=16K²(2K+3)=O(K³)` shared translations suffice for all moments through K, whenever q≤m (eventually true in that range). If q→∞ at any rate, q≤m, the fixed-moment conclusions hold for each fixed k eventually. In particular q may grow as slowly as ceil(log log m). For λ→λ_*>0 this gives uniform Poisson(λ_*) convergence, joint independence at any fixed number of distinct nonopposite targets, and same-factor-base coverage concentration at `1-exp(-λ_*)`. In the sparse regime, the first three moments give the expansion `Pr(H_T=1)=λ-λ²+O(λ³+e_b)`.

**Theorem D (inverse query).** Given one block of d≤m fixed distinct offsets obeying the branch condition and a specified target, all translations that make the chosen address tuple support a signed relation can be found in expected `2^{O(d)} poly(b,h)` bit time. The algorithm uses a univariate norm polynomial of degree `2^{O(d)}`. It does not iterate over p translations or over all choices of offsets. It is an inverse query: the offsets are given, the translation is unknown. The desired forward problem fixes a translation and searches over W offset tuples. No collective forward evaluator is established here.

## 2. Deterministic public branch-avoidance construction

Let D consist of the nonzero differences of distinct roots of f that lie in F_p. There are at most six such differences. They can be found by factoring the cubic in its splitting field (degree at most three for a separable cubic over a finite field), with expected polynomial bit cost. Equivalently, `v∈D` iff `f(A)` and `f(A+v)` have a common root. For any public pool S with distinct elements, the forbidden nonzero slopes are

    Bad = { v/(u-u') : v∈D, u,u'∈S, u≠u' }.

Thus `|Bad|≤6 |S-S\{0}|≤6|S|(|S|-1)`. If this is less than p-1, one can choose one common nonzero slope t outside Bad and replace every S_i by t S_i. The branch conditions hold simultaneously on every block. Construct Bad by enumerating all differences, or enumerate a containing difference set, and choose the least missing nonzero residue. This is deterministic conditional on the computed D, and has expected bit cost `O(|S|² poly(b))` including the small-degree finite-field factorization. At |S|=mw=2^{(1+o(1))b/m} this fits `2^{O(b/m)}`. An interval pool S={0,...,mw-1} has at most 2mw-2 nonzero differences, so its avoidance setup only costs `O(mw poly(b))`. For an irreducible cubic in characteristic >3, D is empty by the earlier Frobenius argument; no multiplier search is needed. For the general cofactor case, D need not be empty (see tests).

This construction is also compatible with the previous polynomial-ring thin digit representation: choose one injective public thin coefficient box of size at least mw, split it by consecutive address intervals of lengths w_i, and apply one public nonzero scalar t to all its field images. Use the prior resultant injection argument on the box, not an assumption that every arbitrary subset has a cheap circuit. Enlarging the box by the factor m costs O(log m) address bits. The original owner-ordered comparisons can impose the two public interval boundaries. This is a qualitative compatibility observation, not a new audited numerical wire ledger.

## 3. One-block geometric estimates

Fix d distinct offsets from one block. The affine tuples in a uniform translation are the points on the fiber-product curve

    y_j²=f(A+s_j), j=1,...,d.

Remove the points above infinity and the finite branch points (zero ordinates). Over the algebraic closure all d square classes are independent since their finite branch loci are disjoint. The smooth projective normalization C is geometrically connected, the map to the A-line has degree 2^d, and its sign involutions act independently. At each of the 3d finite branch points inertia has order 2. At infinity inertia also has order 2, because all cubics have the same odd pole order and their ratios have square roots in the algebraically closed completed local field. Hence Riemann–Hurwitz gives exactly

    2g(C)-2 = 2^{d-1}(3d-3).

The morphisms `ρ_j:C→E`, `(A+s_j,y_j)`, have degree 2^{d-1}. At most `(3d+1)2^d` points are removed. Give each retained rational point of C weight 1/p, mapped to the tuple `(π ρ_1,...,π ρ_d)`. Call this finite measure μ_d. For all characters of G_r^d,

    |mass(μ_d)-1|≤a_d,
    |μ̂_d(χ)|≤a_d when χ is nontrivial,
    a_d = 16 (d+1) 2^d / sqrt(p).

The mass estimate follows from Hasse–Weil. For the character estimate use the unramified rank-one Artin character on E from the connected finite étale Lang map `Frob_p-1:E→E` and the nontrivial order-r character χ_j∘π when χ_j is nontrivial. Its connected degree-r quotient cover remains connected after pullback through ρ_j, since deg ρ_j is a power of 2. Thus its pullback has order r. If the tensor product of the d pulled-back characters were geometrically trivial, pull back by the involution flipping only one nontrivial coordinate. It inverts that factor and leaves all others unchanged, forcing its square to be geometrically trivial. This contradicts r odd. The classical Weil bound for a geometrically nontrivial unramified rank-one Artin character is `(2g-2)sqrt(p)` on the complete curve. Adding the omitted rational points in the error gives the estimate above. The argument remains valid for polynomial cofactor, including rational 2-torsion, and for r=p, because the Lang map is étale.

For a block of slots, a collection of signed collision patterns (defined below) and ordered distinct offsets in each slot defines a pushforward of μ_d under a fixed linear map to G_r^k. Sum over all `(w_i)_(d_i)` ordered distinct offset choices in each slot. Every chosen offset set obeys the same estimate. Consequently the unnormalized block mass is

    (1±a0) ∏_{i in B_j}(w_i)_(d_i),

where `a0=16 (k m+1) 2^{k ceil(m/q)}/sqrt(p)` bounds all a_d. For a0≤1/2, the normalized measure is supported in the block's image subspace V_j; the Fourier coefficient is 1 on its annihilator and has modulus at most `η=2a0` elsewhere. This is a collective bound; its proof does not compute any offset tuple enumeration.

### 3.1 Exact degree of projected sum morphisms

For integers c_j not all zero, put `F=[h]∑ [c_j]ρ_j:C→E`. Then

    deg F = h² 2^{d-1} ∑ c_j².                 (D1)

Use the parallelogram identity for the degree of the pullback of the symmetric divisor (O) on E (the theorem of the square): `deg(f+g)+deg(f-g)=2 deg f+2 deg g`. A constant morphism is assigned pullback degree zero. The involution changing only coordinate j identifies `deg(f+[c_j]ρ_j)` and `deg(f-[c_j]ρ_j)` when f involves only the other coordinates. The degree is therefore additive for these particular maps. Multiplication by an integer t has degree t², including inseparable multiplication in characteristic p. Since deg ρ_j=2^{d-1}, induction proves (D1). In particular F is nonconstant. Every geometric fiber of a nonconstant curve morphism has at most its degree many geometric points. Hence for a specified group value and uniform translation the unnormalized signed-lift count for `π∑ c_j P_j` is at most `h² 2^{d-1} ∑ c_j²/p`. Normalizing a measure of mass ≥1-a0 costs a factor ≤2. For d≤2m and |c_j|≤2 this is `2^{O(m)} poly(b)/r` (use N≈p, h poly b). We use this stronger point-probability estimate only for the small integer coefficient choices in Theorems A and B.

## 4. Signed partitions and exact block convolution

For k raw replicas in a slot partition their indices by equal abscissas. In each part choose the lowest-index replica's point as representative and attach ± signs relative to it, first sign +. Zero ordinates were excluded, so this is unique even when several raw points have the same projected group value. The type P is a signed partition, with d(P) parts. The associated `k×d(P)` matrix has one signed unit entry per row and disjoint column supports; its rank over F_r is d(P). The image V_P determines P uniquely: coordinate functionals restricted to V_P are proportional up to sign exactly in the same part, and determine the relative signs. There are at most `B=(2k)^k` types. A block type sequence has image subspace `V_j=∑_{i∈B_j} V_{P_i}`; its dimension may exceed any individual slot's. The normalized block measure from §3 is supported in V_j and has the stated Fourier estimate there. Different blocks are independent.

For the mixed factorial moment, a type sequence over all slots is admissible if no two replica indices have identical rows in every slot. This is precisely the condition of distinct raw signed tuples. If a pair of tuples were identical then their projected target would be identical. The mixed falling-factorial moment is exactly the sum, over admissible sequences, of the convolutions of the unnormalized block measures at the target vector. This uses block independence, not slot independence within a block. It includes cofactor multiplicities and possible projected identities of individual raw leaves.

If all slots have a maximal type contained in the same proper subspace V of dimension R<k with d(P)=R, they have the same type (by uniqueness of V_P). Some pair of replica indices then has a fixed relative sign in every slot. If it is minus, the prescribed targets are opposite and excluded; if it is plus, the replicas are identical and inadmissible. Thus such an all-maximal sequence has zero contribution. This elementary observation is needed below to suppress the full diagonal.

## 5. Proof of Theorem A: one shared shift

With q=1, the expected total number of raw tuples is `W(1±a0)` by the mass estimate. Let T be uniform in G_r. Then `E H_T=λ(1±a0)`. For a pair of distinct raw tuples, shared abscissas within the same slot may have equal or opposite signs; use k=2 signed partitions. The difference of their sums is a linear combination of at most 2m representative points with coefficients 0,±1,±2, not all zero (a raw nonidentity pair has some differing point). The degree estimate (D1) bounds the expected number of rational lift pairs with equal projected sums, for a fixed offset pattern, by `2^{O(m)} poly(b)/r`. There are at most 3^m signed type sequences, and at most W² choices of offsets for any type. Summing gives

    E_A ∑_{T∈G_r} H_T(H_T-1) ≤ W² 2^{O(m)} poly(b)/r.

Divide by r for the uniform target average. To exclude T=O, note separately that for every fixed offset tuple the sum morphism is nonconstant and (D1) gives `E_A H_O ≤ W 2^{O(m)} poly(b)/r`, negligible compared with W. Therefore, for T uniform nonidentity,

    E H_T=λ(1+o(1)),
    E H_T² ≤ O(λ)+λ² 2^{O(m)} poly(b).

By Cauchy–Schwarz `Pr(H_T>0) ≥ (E H_T)² / E H_T²`, which proves the asserted average coverage. This does not prove a targetwise lower bound for q=1.

## 6. Proof of Theorem B: five shared shifts

For k=1 all five block spaces have dimension one and the same full target span. Their Fourier convolution is pointwise `r^{-1}(1±rη^5)`. The first moment follows from the mass factors.

For k=2 there are just three slot types: distinct offsets (rank 2), same offset with same sign (line L_+), same offset with opposite signs (line L_-). A rank-one block must consist entirely of the same rank-one type, since the two lines are distinct (r odd). A rank-two block may contain any combination of slot types.

If all five block spaces have rank 2, the normalized convolution is pointwise `r^{-2}(1±r²η^5)`. The all-distinct slot sequence contributes

    λ² ∏_i(1-1/w_i) (1±a0)^5 (1±r²η^5).

All other sequences of full-rank blocks are bounded in total by

    λ² (1+a0)^5 (1+r²η^5) [ (1+2/v)^m-1 ],

where v=min w_i; use the unique rank-2 slot type and two rank-1 types. This upper bound sums more sequences than necessary, which is safe. Since `r²η^5=2^{-b/2+o(b)}`, these terms give the claimed main term and error.

If some block space is a line L and total span has dimension 2, choose a second block not contained in L. The linear functional ℓ annihilating L has coordinates (1,±1) (up to a nonzero scalar) and is nonzero on the second block. Condition on every other block. The second block must first meet a prescribed ℓ-value, which by (D1) has probability at most `2^{O(m)}poly(b)/r`: its linear sum uses coefficients of magnitude at most two and not all zero. Given that block's entire value, the line block must meet one specified point in L. Its probability has the same bound, using either coordinate on L and (D1). Independence gives a joint point probability `2^{O(m)}poly(b)/r²`. A rank-one block loses one offset choice in *every slot of that block*, at least floor(m/5) factors of v relative to W². Even summing over all 3^m sequences, the total contribution is bounded by

    λ² v^{-floor(m/5)} 2^{O(m)} poly(b) = 2^{-b/5+o(b)}.

If total span has rank one, all blocks are the same line. For nonopposite targets the minus line cannot contribute, and the plus line means identical raw tuples, excluded from the factorial count. This exhausts the cases. The small first-moment and second-moment errors are `2^{-b/m+o(b/m)}+2^{-Ω(b)}`, bounded by e_b eventually. At λ→1, Paley–Zygmund gives 1/2-o(1) coverage; for sparse λ the inequalities `E H-E(H)_2≤Pr(H=1)≤E H` give the stated singleton estimate. No use of independence between slots of the same block was made.

## 7. Proof of Theorem C: robust-span pruning replaces frequent-type counting

For a fixed full type sequence start with all q block spaces and their span V. If a nonzero linear functional on V is nonzero on fewer than s of the current block spaces, delete the blocks on which it is nonzero and recompute the span. Each step deletes at most s-1 blocks and drops the dimension by at least one: all remaining spaces lie in its kernel. After at most k steps stop. Let V be the remaining span of dimension R; the deleted set contains at most J=k(s-1) blocks. Since q>J, some block remains and R≥1. Reinsert any deleted block lying in V. The number outside V is still at most J; every nontrivial character on V is nontrivial on at least s of the inside blocks. No enumeration of hyperplanes is needed by a relation algorithm: this is a proof device to bound the moment. It improves on counting all possible block subspace types, whose number may be much larger.

For the inside convolution Fourier inversion gives normalized point probabilities

    r^{-R}(1±δ), δ=r^k η^s.

Since `s=2k+4`, `δ≤p^{-2+o(1)}` in the stated regime (the exponent error is O(k²m)=o(b)). The product of all block mass upper factors is `(1+a0)^q`. Use `v=min w_i`, `U=max w_i`, `Aλ=max(1,C^k)`.

**Case I: R=k.** The all-singleton slot-type sequence (all k offsets distinct in each slot) has full space in every block. Its contribution is

    λ^k ∏_i ((w_i)_k/w_i^k) (1±a0)^q (1±δ).

Its error from λ^k is `O(λ^k(m k²/v + q a0 + δ))`. Every other sequence has a deficit `k-d_i≥1` in some slot. Sum their contribution using the unique full type and at most B deficient types per slot to bound by

    (1+δ)(1+a0)^q λ^k [ (1+B/v)^m -1 ].             (C1)

**Case II: R<k, no outside blocks.** For fixed V, each slot type space lies in V. At most one has d_i=R, and the all-maximal sequence cannot contribute by §4. All remaining sequences lose at least one factor of v. Every possible V is spanned by at most k signed type columns from `{0,±1}^k`, so there are at most `(k+1)3^{k²}` possible spans, regardless of r. The case is bounded by

    (k+1)3^{k²}(1+δ)(1+a0)^q Aλ [ (1+B/v)^m -1 ].    (C2)

The precise number of candidate spans is immaterial; its logarithm is O(k²), rather than a factor exponential in the number of types.

**Case III: outside blocks exist (R<k).** Pick one outside block and a functional on F_r^k annihilating V but nonzero on that block's space. Every nonzero scalar multiple is a nontrivial character on that block. Its projected point probability is at most `1/r+η` by one-dimensional Fourier inversion. Condition on the other outside values; membership in the desired coset of V requires that projected equality. The inside sum has point probability at most `(1+δ)/r^R` given the outside values. Thus the joint point probability is at most `(1+δ)(1/r+η)/r^R`. There are at most J outside blocks containing at most `J ceil(m/q)` slots. In the inside slots d_i≤R, and in the outside slots d_i≤k. Consequently `∏ w_i^{d_i} ≤ W^R U^{kJ ceil(m/q)}`. Sum over all B^m sequences, even inadmissible ones, to get the bound

    (1+δ)(1+a0)^q Aλ B^m U^{kJ ceil(m/q)}(1/r+η).   (C3)

The exponents are as follows. For `k≤log m/(4 log log m)`, B≤m^{1/3} eventually, `m log B=o(b)`, `k²=o(b/m)`, and `kJ=o(m)`. Because q≥16kJ,

    log2(U^{kJ ceil(m/q)}) ≤ (1/16+o(1))b.

Also `η=2^{-b/2+o(b)}`, `q a0=2^{-b/2+o(b)}`, `δ≤2^{-2b+o(b)}` and `log2 Aλ=o(b/m)`. Thus (C3) is `2^{-7b/16+o(b)}`. Equations (C1),(C2) and the main-term error are `2^{-b/m+o(b/m)}` since `(1+B/v)^m-1=O(mB/v)`. Together these prove the uniform error e_b. The proof is uniform for balanced partitions and arbitrary public supports satisfying §2, not a statement for adversarially chosen realized translations.

For λ→λ_*, convergence of each fixed mixed factorial moment gives the joint Poisson law by the usual method of moments; the Poisson distribution is moment-determinate. Opposite targets are excluded because `H_{-T}=H_T` exactly. To get same-base coverage concentration let Z be the fraction of nonidentity projected targets with H_T>0. Its mean tends uniformly to `1-exp(-λ_*)`; its variance tends to zero using pairwise convergence for distinct nonopposite T,U and bounding the O(r) exceptional ordered pairs by one. Sparse expansion follows from the first three factorial moments and the Bonferroni inequalities for integer counts, as in the prior checkpoint.

## 8. Constructive inverse translation query

Take d fixed offsets s_i in one shared block. Work in the polynomial algebra

    R=F_p[A,y_1,...,y_d]/(y_i²-f(A+s_i), i=1,...,d).

Its localization at nonzero polynomials in A is a field, by the disjoint branch loci and independent square classes. Each element of R is represented in the basis of 2^d squarefree y monomials with univariate coefficients in A. Polynomial arithmetic on this basis is explicit, using the reduction rules y_i²=f(A+s_i).

Let P_i=(A+s_i,y_i) be the generic points. Construct their generic sum in mixed Jacobian coordinates `x=X/Z², y=Y/Z³`. Start `(X,Y,Z)=(A+s_1,y_1,1)`. On adding the next affine point `(u,v)=(A+s_j,y_j)`, set

    H=u Z²-X, D=v Z³-Y, V=X H²,
    X'=D²-H³-2V,
    Y'=D(V-X')-Y H³,
    Z'=Z H.

The identity follows from the ordinary chord formula on H Z≠0. Generically H≠0: otherwise the preceding partial sum equals ±P_j, contradicting the nonconstancy of the corresponding signed sum morphism (D1), or directly its independent sign involutions. Z is nonzero in the generic field. All actual specialization paths that hit an exceptional chord have Z=0 at that step and thereafter, so the final Z detects those paths even though the other coordinates may no longer represent points there.

For a finite target U, form the nonzero algebra element

    R_U = Z (X-x(U) Z²).

It is nonzero since the sum morphism is nonconstant. For U=O use `R_U=Z`. Compute its field norm to F_p(A) by multiplying with its image under each successive independent sign involution; after each multiplication that variable's sign disappears. The result lies in F_p[A] and is a nonzero polynomial. Its specialization is the product over all 2^d sign choices (even when some roots coincide); consequently every translation supporting a genuine signed relation at U is among its field roots. The same holds for exceptional paths because of the factor Z. Spurious roots (opposite sign, exceptional but wrong sum, nonrational lifts) are allowed and verified individually.

Use weighted degree with wt(A)=2, wt(y_i)=3; reduction of y_i² never increases weighted degree, and the mixed-addition formulas multiply the maximum degree by at most an absolute constant and add a constant each step. Therefore the polynomial degree of every coefficient after d additions is 2^{O(d)}. Norm computation increases it by at most another factor 2^d, and intermediate algebra sizes are 2^{O(d)}. Direct arithmetic, polynomial factoring over F_p (first gcd with A^p-A), and checking every real root using all 2^d signed lifts cost `2^{O(d)}poly(b)` expected bit operations. The exact nonzero-norm argument ensures root enumeration is finite; no assumption of generic specialized behavior is used. For projected nonidentity target T, try all h full-group preimages U of T, with polynomial extra overhead. This proves Theorem D. At d≈sqrt b, a fixed-address inverse query meets a `2^{O(sqrt b)}` scale.

**Scope of the primitive.** The inverse procedure fixes the address and finds translations. Forward relation search fixes the shared translation and needs an address; there are W≈r λ such choices. Taking a product of all the inverse polynomials or calling the inverse evaluator once per address has exponential-in-b cost. We have not established a fast structured aggregation of these polynomials. Nor can one simply generate a new factor base tailored to each target using the inverse algorithm and then reuse the previous linear-system reduction: that reduction requires fixed labeled columns across its training rows.

## 9. Cofactor handling, regularity, and cost accounting

The target equations are for `π=[h]` as defined in §1. The kernel of π on E(F_p) has h elements because gcd(h,r)=1. All preimages of T∈G_r can be listed from one preimage `[h^{-1} mod r]T` and the kernel. Obtain a uniform kernel sample by multiplying a uniform rational point by r; coupon collection lists the entire polynomial-sized kernel in expected polynomial time. A uniform rational point can be sampled by uniform x and a random root index with rejection, adding the identity with the appropriate rational probability (or uniform sampling of the affine points followed by the known identity mass 1/N). To find a projected relation, try each of these h raw targets using an exact original-group evaluator and verify the sum. This adds only polynomial overhead. In the one-pilot linear-algebra reduction choose T=`[h]([u]G+[v]Q)` with independent uniform u,v modulo r, excluding identity. Since [h]G generates G_r, v remains uniform conditional on T. Relation rows refer to the projections of the same labeled factor-base points. This is the prior size-biased pilot argument with a new sampling law; the sampler uses one random translation per block during base construction and keeps it fixed through all queries on that base.

Within a block the x sets are disjoint exactly. Between blocks a collision has probability at most `∑_{i<j}w_iw_j/p`; avoiding the abscissas of all raw target preimages costs at most `h∑_i w_i/p`. Exclusion of rational 2-torsion is already part of the sampler, and its failure probability if tested before enumeration is ≤`3∑ w_i/p`. These bounds are `2^{-b+o(b)}` and negligible even relative to the average coverage in Theorem A. With the section norm's regular-support predicate, one may reject whole seeds/targets failing these checks. For the older guarded chord chain, common-shift correlations change the independent-slot guard estimate; (D1) bounds a specified nontrivial zero-sum guard failure in one shared block, but we do not claim the old guard coverage proof transfers unchanged for every block partition. The regular section predicate or an exact ordinary relation evaluator is the appropriate assumed interface here.

The one-pilot ratio only needs average success. Theorem A shows that using one translation adds at most `O(m+log b+log λ^{-1})` to the logarithm of retry cost. Theorem B retains a targetwise constant (or sparse polynomial) success guarantee using five independent translations, and Theorem C gives the full limiting law with q=o(m). Branch-avoidance preprocessing has log cost at most `O(b/m)` at balanced widths. With the inherited O(mb)-wire arithmetic graph, however, the log query cost remains `O(mb/sqrt(log b))` under the generic contractor; optimizing b/m against it still yields `O(b/(log b)^{1/4})`. Seed entropy is **not** the number of transmitted certificate bits. Do not infer a smaller wire count from q alone. The older 11.31486814 leading coefficient was conditional and is neither improved nor revalidated here.

The independent-slot seed-local slice-rank and communication obstructions of earlier checkpoints cannot simply be applied to this shared-translation model: there are q independent seeds rather than m. In particular the lower bounds proved by varying one independent slot mask while holding the others fixed do not have their original hypotheses. This is an opening for a block-correlated arithmetic compiler, not a proof of small multilinear degree or a universal computational lower bound. The inverse norm primitive identifies a different target for collective computation. None of the results here solves that computation or achieves the stronger general-prime exponent.

## 10. Verification and dependencies

Run from the package directory with Python 3 and NumPy:

    OPENBLAS_NUM_THREADS=1 python -B verify.py
    OPENBLAS_NUM_THREADS=1 python -B inverse_translation_verify.py

The first verifier compares literal shared-translation enumeration, exact signed pattern measures, and complete falling-factorial arrays using auxiliary-prime NTTs; a product-of-moduli bound ensures equality as integers. It includes polynomial-cofactor fixtures, 2-torsion exclusion, and the five-block configuration. It also checks the robust-span deletion identities, the exact forbidden-slope set, rational projected fibers against (D1), and numerical FFT character bounds. The second verifier implements the squarefree multiquadratic polynomial arithmetic, mixed Jacobian sum, and norm calculation, checking every rational lift and literal root set in four small cases. Its toy root search exhausts the small field; the asymptotic algorithm's root search uses standard finite-field factoring. The recorded transcripts distinguish numerical regressions from exact comparisons. No large-scale forward relation evaluator or discrete logarithm computation was run.

Sources: `replica_mixing_note.md` for the preceding independent-translation proof and its Lang-cover argument; v3 of the uploaded paper for the inherited reduction and target exponent; Michel–Tibouchi character-sum exposition, https://www.normalesup.org/~tibouchi/papers/welldistributed.pdf (classical unramified Weil estimate); Kohel–Shparlinski, *On Exponential Sums and Group Generators for Elliptic Curves over Finite Fields*, https://www.i2m.univ-amu.fr/perso/david.kohel/pub/character.pdf (elliptic group characters); standard theorem of the square for the symmetric line bundle (O) on an elliptic curve. The robust-span pruning, shared-block factorial moment analysis, branch-avoidance transfer, one-shift average coverage proof, and inverse norm query are the new derivations in this checkpoint. Novelty relative to all literature has not been established.
