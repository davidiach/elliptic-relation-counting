# A slice-rank obstruction to low-degree relation compilers

24 November 2026. Research checkpoint following `contractor_degree_note.md` and `affine_digit_sampler_note.md`. This is a negative result for an explicit model of exact relation compilation, not a lower bound for ECDLP or arbitrary tensor contraction. No improved general-prime ECDLP exponent is claimed. The result is self-contained; priority relative to unpublished work has not been checked.

## 1. Main statement

Let `p` be an odd field prime, `b=ceil(log2 p)`, and `G=E(F_p)` have a supplied odd prime-order subgroup `C` of order `r`, with polynomially bounded cofactor. Fix `T in C\{O}`. For each of `m>=2` slots, choose a fixed nonempty proper set `S_i subset F_p`. Its seed `(a_i,t_i)` ranges over `F_p x F_p^*`; the sampled point set contains every affine point with abscissa in `a_i+t_i S_i`, retaining both signs. Let `H_T` be the ordinary ordered relation count, including all intermediate group-law cases, as a tensor indexed by the independent seeds.

A **seed-local network** is a fixed bipartite Boolean tensor network whose vertex functions may depend on public curve/target data and on at most one slot seed each. It may have seed-free auxiliary tensors, arbitrary local computations, and dependencies on the seeds on either side. Its edge count is `M`. On the chosen row side, the product of the row tensors has a multilinear polynomial of degree at most `D` in the edge bits, uniformly for every seed tuple. Degree may be measured modulo two, or over a dyadic coefficient ring at the precision used by the degree-sensitive contractor. No restriction on the rank or topology of any single vertex is needed.

**Theorem.** If such a network computes `H_T mod 2` exactly on every seed tuple, and `s=(r-1)/2>2m`, then

    M D > (m^2/16) [log2(s/(2m))]^2.

The same conclusion with `log2(s/(2Km))` holds for at least one of `K` networks if a fixed, seed-independent linear combination of their outputs gives exact parity, with `s>2Km`. It also holds for a fixed linear combination equal to `2^(P-1) H_T mod 2^P` over a dyadic unramified ring, using the degree at precision `P` (see Section 5). Thus a fixed power-of-two normalization does not remove the obstruction. More general output normalizations must be justified separately.

For `m=polylog(b)`, `K=2^o(b)`, and polynomial cofactor, the right side is `Omega(m^2 b^2)`. In particular an inflation of the wire count to reduce the degree cannot reduce the product `M D` by a growing factor in this model. The theorem applies to ordinary count parity, not automatically to the parity of arbitrary nonunique certificates, guarded-only counts, exact existence, or a fallible verified-witness oracle.

## 2. Full paired Fourier diagonal from the actual affine seed family

All rank arguments in this section are over an algebraic closure `k` of `F_2`, so primitive roots of unity of any odd order exist. This is an analysis field; no large extension arithmetic is required in the proposed algorithm.

### 2.1 Affine mask span

Write `1_S` for a nonempty proper support's indicator vector on `F_p`. Its translation Fourier transform has at least one nonzero nontrivial coefficient: otherwise inversion would make the vector constant. Multiplicative dilations permute the `p-1` nontrivial additive characters transitively. The linear span of all affine images is translation invariant. Applying the character projection, a linear combination of translations using Fourier orthogonality (p is a unit), shows that it contains every nontrivial character line. Thus it contains the entire zero-sum subspace of `k^{F_p}`. When `|S|` is odd it contains the constant character line as well. This also proves the corresponding rank over `F_2`, since rank does not change upon extension of the coefficient field.

Projection of the zero-sum space onto any proper subset of coordinates is surjective: adjust one omitted coordinate to make the sum zero. The set of rational abscissas is proper for all sufficiently large p: Hasse's bound gives at most `(p+3+2 sqrt(p))/2` up to the constant contribution of two-torsion, which is less than p eventually. More precisely, it is `(N-1+t)/2`, where `N=|G|` and `t` is the number of nonidentity rational two-torsion points, at most 3.

It follows that linear combinations of the actual seed vectors realize **every** function on the rational abscissas. In particular they realize functions supported on the pairs `{P,-P}` of the given odd subgroup `C`, zero on all other rational points. This is a linear restriction of the count tensor, not an assertion that a single mask can be chosen to isolate every such pair.

### 2.2 Fourier diagonal on the odd subgroup

Identify `C` abstractly with `Z/r`, for analysis only. Let `V` be its space of symmetric functions vanishing at O. It has dimension `s=(r-1)/2`. For `f in V`, the group Fourier transform has `fhat(1)=0` in characteristic two and `fhat(chi)=fhat(chi^{-1})`. The invertible group Fourier transform therefore identifies `V` with all arbitrary coordinates indexed by the s nontrivial inverse-character pairs.

Ordinary fixed-target count is convolution, hence

    H_T(f_1,...,f_m)
      = (1/r) sum_chi chi(-T) product_i fhat_i(chi)
      = (1/r) sum_{chi modulo inverse, chi!=1}
                    (chi(T)+chi(-T)) product_i fhat_i(chi).

Every displayed coefficient is nonzero. Indeed, equality `chi(T)=chi(-T)` would imply `chi(2T)=1`, impossible for a nontrivial character on a prime-order subgroup and `T!=O`. This establishes a diagonal tensor of length s with nonzero entries. For a prime-order full group this is an exact diagonalization of the seed tensor up to injective local maps. In the cofactor case it is a linear restriction, sufficient for the lower bound.

### 2.3 Slice rank

A slice is `u(x_i) v(x_1,...,x_{i-1},x_{i+1},...,x_m)` for some choice of mode i. Slice rank is the smallest number of slices whose sum gives the tensor. Applying separate linear maps on the modes cannot increase slice rank. A diagonal tensor of length s with nonzero diagonal entries has slice rank s. Here is a short proof. For m=2 this is ordinary matrix rank. For m>=3, suppose a decomposition has k_1 slices in mode one. Choose a vector orthogonal to their first-mode factors, with at least s-k_1 nonzero coordinates. Such a vector exists: their annihilator has dimension >=s-k_1; choose a set of independent pivot coordinates and prescribe 1 on them. Contract along this vector. The first-mode slices vanish, and the remaining tensor has a diagonal restriction of length at least s-k_1, with at most k-k_1 slices. Induction yields k>=s. The reverse inequality is the trivial sum of s diagonal slices.

Consequently the actual ordinary relation count as a tensor of affine seeds satisfies

    slice_rank(H_T mod 2) >= s.

No generic-position assumption on the support S_i is needed, and the width can be even or odd. The same lower bound holds when any disjoint groups of slots are treated as single modes, provided at least two groups remain: after the preceding local linear restrictions and Fourier maps, the coordinate vector for character chi in each grouped mode is `e_chi^{tensor group_size}`. These vectors are independent and can be mapped back to the s standard coordinate vectors. The diagonal is preserved.

## 3. Degree-to-slice-rank lemma for a general seed-local network

Let the row polynomial product be expanded as

    A(X; seeds) = sum_{J subset edges, |J|<=D} a_J(seeds) product_{e in J} X_e.

The row port sets are disjoint, so `a_J` is the product of the selected local multilinear coefficients. The columns have disjoint port sets as well. Summing each column tensor times its specified monomial gives the local moment associated to the subset J of that column's incident edges. The entire contraction therefore has an expansion

    Z(seeds) = sum_{|J|<=D} c_J product_i g_{i,J intersect E_i}(seed_i).

Here `E_i` is the union of edges incident to a tensor that reads seed i; `c_J` includes all seed-free factors. The factor for seed i depends only on the indicated subset, because both the coefficient selections and the local moment selections are determined by that subset. This works even when the two ends of an edge depend on two different seeds. In that case it can be counted twice, but never more. Write `c_i=|E_i|`. Then

    sum_i c_i <= 2M,   sum_i |J intersect E_i| <= 2|J| <= 2D.

For arbitrary nonnegative integer thresholds `k_i<=c_i`, if `sum_i(k_i+1)>2D`, every term has at least one i for which `|J intersect E_i|<=k_i`. Assign the term to the least such i. Group terms with the same i and the same local subset. Each resulting group is a single slice. This yields the **exact bound**

    slice_rank(Z) <= sum_i sum_{j=0}^{k_i} binom(c_i,j).

Without any condition on D, one also has `slice_rank(Z)<=2^c_i` by grouping every monomial using a fixed i. No enumeration of these monomials is required for the lower-bound argument. Isolated seed-local factors cause no problem. In the one-sided seed model, each edge occurs in at most one E_i; the factors 2 in the two edge-count bounds can be removed.

### 3.1 An elementary entropy bound

For integers `1<=k<=c`, take `x=k/c` in `(0,1]`. Then

    sum_{j<=k} binom(c,j) <= x^(-k)(1+x)^c <= (e c/k)^k.

The maximum of `ln(e t)/sqrt(t)` for t>=1 is `2/sqrt(e)`, at t=e. Therefore

    log2(sum_{j<=k} binom(c,j)) <= 2 sqrt(c k).

It also holds for k=0.

### 3.2 The product lower bound, including multiple contractions

Let `B=log2(s/(2Km))>0`. For a particular network, if any `c_i<=B`, its slice rank is at most `2^B`. Otherwise set

    k_i=floor(B^2/(4 c_i)),  so 0<=k_i<c_i.

All local subset counts are at most `2^B`. Moreover Cauchy--Schwarz gives

    sum_i (k_i+1) > (B^2/4) sum_i 1/c_i
                   >= B^2 m^2/(8M).

If `M D <= B^2 m^2/16`, this strictly exceeds `2D`. Hence the contraction has slice rank at most `m 2^B=s/(2K)` in either case. If every one of K networks obeyed that degree/edge product bound, their fixed linear combination would have slice rank at most s/2 by subadditivity. This contradicts Section 2. Taking K=1 proves the main statement. In the one-sided seed model the constant 1/16 can be strengthened to 1/4 by the two sharper edge inequalities.

The same reasoning works with slots replaced by disjoint blocks. The relevant arity in the lower bound is then the number of blocks, not the number of original slots.

## 4. Joint seed dependence: how much is necessary?

Here consider the odd prime-order full-group case and the common support `S` of size w (as in the affine digit-block sampler). Allow an arbitrary subset of the slot seeds at each tensor. Form the **seed dependency graph** on the m slots: connect two if some tensor reads both. For several networks use the union graph. If I is an independent set of size q>=2, every tensor reads at most one active seed after all outside seeds are fixed. Choose the *same* nonempty actual affine mask F for every outside slot. Such a mask exists by translating a nonempty support onto any rational abscissa. This restriction need not have high probability under the fresh law; the exact compiler is assumed correct on every tuple.

Let `f=1_F` on the group. Its Fourier transform has at least `r/|F|` nonzero coordinates: the circulant matrix of translations of f has rank equal to Fourier support size. The supports of a basis of its rows must cover all r coordinates (every column is nonzero), and each such row has support size |F|. This proves the bound over the characteristic-two splitting field without invoking a characteristic-zero uncertainty theorem. Since the mask is symmetric, excludes O, and `|F|<=2w`, there are at least `r/(4w)` live nontrivial inverse-character pairs. The coefficient of each such pair in the restricted q-mode diagonal is the earlier nonzero coefficient times `fhat(chi)^(m-q)`, still nonzero. If no slots are fixed, use s itself.

It follows from Section 3 that for at least one network

    M D > q^2/16 * [log2((r/(4w))/(2Kq))]^2

whenever the logarithm is positive (use s instead if all slots are active). At balanced width `log2 w=b/m+o(b)`, `K=2^o(b)`, and m tending to infinity, the logarithm is `(1-o(1))b`. In particular, for dependency maximum degree Delta with a greedy independent set of size at least `m/(Delta+1)` and at least two, this yields `MD=Omega(m^2 b^2/(Delta+1)^2)`. Adjacent-only joint seed processing (bounded Delta) does not avoid the asymptotic obstruction. A dense dependency graph is not sufficient for a fast algorithm; it is only a way to evade this proof.

For M=O(mb), the degree-sensitive promise `D <= M/(L^(1/3) log L)` with `L=log2 b` therefore requires `Delta+1=Omega(L^(1/6) sqrt(log L))` in this model, or requires leaving the exact-count assumptions. This is a quantitative necessary degree of joint dependence, not a sufficient construction.

## 5. Dyadic normalization and multiple precisions

The degree-sensitive contractor uses dyadic rings. Suppose at some common precision P, after clearing a fixed scalar normalization, a seed-independent linear combination of the networks equals `2^(P-1) H_T mod 2^P`. Allow row polynomial degrees measured over the ring at P. The slice decomposition of Section 3 is purely ring algebra, so still applies.

Take a sufficiently large finite unramified extension of `Z/2^P` for the rank proof, with residue field containing the finite splitting field used in Section 2. This extension is a proof device. Lift all local linear restriction and diagonalization matrices from the residue field to the extension ring. Since the target tensor is multiplied by `2^(P-1)`, the resulting tensor is diagonal with nonzero entries `2^(P-1)` times unit lifts of the diagonal coefficients. Off-diagonal entries vanish exactly at precision P. Separate linear maps do not increase slice rank over the ring.

A diagonal tensor with nonzero diagonal entries over this finite chain ring has slice rank equal to its diagonal length. The same induction as before works with a unit-support vector in the kernel of the first-axis slice factors. For completeness, a matrix with k_1 rows over the ring has a kernel containing a free direct summand of rank at least s-k_1: use Smith normal form over the principal ideal ring and the unconstrained coordinates. The reductions of the free generators are linearly independent over the residue field. A linear combination can be chosen with at least s-k_1 unit coordinates (use an invertible minor of their reductions); lift its coefficients and retain exact kernel membership. Contract and restrict to those unit coordinates. Nonzero diagonal entries remain nonzero after multiplication by units. For m=2 the same kernel argument shows that fewer than s rank-one terms are impossible. Induction for higher m then proves the claim.

This proves the same product lower bound for the scaled output. If network outputs are used at different precisions, the statement applies only when the final readout can be expressed as a fixed linear combination over a common ring after clearing denominators and using sufficient precision. Arbitrary nonlinear decoding, seed-dependent combination coefficients, and insufficient-precision lifts are not included.

## 6. Consequence for the proposed low-degree route

The earlier degree-sensitive root bound has an exponential term of order

    sqrt(M D log r_groups / log M),

along with child-contraction, precision-floor and entry-access terms. In the old relation setting `M=Theta(mb)` and `log M=Theta(log b)`. Our new theorem forces the quantity `sqrt(MD/log M)` to be at least order `mb/sqrt(log b)` for slot-local exact parity compilers, even if their wire counts grow polynomially. This is a statement about the *parameter expression of that upper-bound method*, not a lower bound on the true running time of every contraction algorithm. The contractor's analysis could be sharpened by a different mechanism.

The disjoint-block variant gives a specific further limitation. Group the m original slots into q blocks, with largest block containing at least m/q slots. If a proposed compiler explicitly enumerates all tuples in each block, at balanced original width it pays log time at least `b/q` up to lower-order terms. The grouped diagonal bound forces `MD=Omega(q^2 b^2)` (for exact local parity over the q blocks) so the degree-sensitive expression is at least `q b/sqrt(log b)` when M is polynomial in b. Balancing these two quantities gives scale `b/(log b)^(1/4)`. Thus merely enumerating blocks to get a smaller effective degree does not certify a stronger logarithmic saving by this method. The conclusion does not cover joint block computations that avoid explicit block enumeration.

For comparison, in a genuinely binary additive group `(F_2)^b`, the exact sum test on b-bit labels is a product of b parity predicates, each of Boolean degree one on its m incident bits. The natural representation has M=mb and D=b. This illustrates why the desired low-degree promise would be possible for an additive binary group. It is not an ECDLP compiler: odd prime-order group convolution is Fourier-diagonal over characteristic two, whereas the binary group algebra in its own characteristic is not semisimple. The binary example also does not by itself satisfy every local degree-cap hypothesis of the contractor for arbitrary mask tensors.

## 7. What is and is not ruled out

The result covers arbitrary seed-local graph redesigns and edge inflation, not just per-word gauge transformations of the existing Fourier compiler. It also covers several subexponentially many exact local parity calls combined linearly, and sparsely overlapping joint seed computations via the independent-set restriction. It requires exact ordinary parity for every seed tuple. It does not rule out the following:

* A verified relation finder that succeeds on an appropriate fresh distribution and does not give exact ordinary count parity.
* A compiler with substantive joint seed preprocessing, or a dense dependency graph, that cannot be expanded into a subexponential-size sum of slot-local contractions.
* A relation count with variable certificate multiplicities whose exact parity is not the ordinary parity, or a nonlinear readout.
* Different deformations, coefficient-pairing algorithms, or contraction estimates that do not pay the present `sqrt(MD/log M)` term.
* All algorithms for ECDLP or all low-degree Boolean networks.

A next viable research gate is to construct a **joint** representation or an average-case detector outside the theorem's model, and prove charged runtime and verified recovery. An exact local low-degree reformulation of the present count is no longer a credible route to a stronger exponent without breaking a specific hypothesis. The best previously proposed coefficient remains 12.39481703 at the old 1/4 saving power, with its inherited assumptions and compiler audit qualifications; no new coefficient or exponent is claimed here.

## 8. Finite checks and reproduction

`python -B seed_local_slice_rank_audit.py` (standard library only) regenerates the retained transcript. The tests use four smooth small curves of odd prime group order, including even-width and odd-width affine supports. They verify mask spans, exact paired Fourier diagonals and fixed-mask restrictions in explicitly constructed binary extension fields, and directly compare small seed-local Boolean networks with their monomial-to-slice decompositions. The large arity-4 tuple sets use 1,000 deterministic pseudorandom examples where exhaustive enumeration would exceed 1,000 tuples; smaller sets are exhaustive. These tests establish finite identities, not the asymptotic theorem, which is proved above. The exact diagonal slice-rank argument and the degree-product inequality do not depend on experiments.

Sources: uploaded v3 `Counting relations in independently hashed elliptic-curve slots`, especially the signed slot model; prior checkpoint notes `affine_digit_sampler_note.md` and `contractor_degree_note.md`; the project's `derivations/strided_slots/results.md` and earlier seed-cut and separated-rank bounds, revision `8171080efb0d9024f49eb47e1e311914f0f7a4c8`. The diagonal slice-rank lemma, Fourier orthogonality, and Smith normal form are standard; proofs/applications needed here are included. This is AI-assisted mathematical research, not independently peer reviewed.
