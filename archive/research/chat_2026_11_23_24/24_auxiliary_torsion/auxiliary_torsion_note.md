# Auxiliary elliptic torsion supports and compact symbolic fiber elimination

Continuation of `coordinate_deformed_resultant_note.md`, `fiber_gcd_note.md` and `shared_translation_note.md`.

**Status.** This checkpoint removes the exponential slot-width factor from symbolic support reduction and from the construction of a fiber-indicator circuit. It does not give a fast evaluation of the collective trace of that circuit, nor a smaller general-prime ECDLP exponent. In particular, a gate over a remaining product algebra is not a scalar field operation. The previously proposed coefficient 11.31486814 and the full sampler-to-ECDLP reduction are not re-audited here. No novelty or priority claim is made.

## 1. A large support with a logarithmic-size membership circuit for arbitrary field primes

Let E/F_p be a smooth short Weierstrass curve y²=f(x)=x³+ax+c, p>3, of order hr, where r is an odd prime, h=poly(b), b=ceil(log_2 p). Take m→∞, m≤sqrt(b). In particular r=p/poly(b) up to a constant factor by Hasse. All asymptotic assertions are for sufficiently large b.

Choose a prime ell in [r^(1/m)/8, r^(1/m)/4] with ell not dividing p(p²−1). The prime number theorem in dyadic intervals supplies such a prime: the excluded prime divisors in that interval number O(b/log ell)=O(m), whereas the interval contains asymptotically many more primes. Enumerating candidates with polynomial-time primality tests costs 2^{O(b/m)} poly(b). Put w=(ell−1)/2 and W=w^m. Then lambda=W/r=2^{−Theta(m)}. This choice does not require factoring p−1 or p+1.

**Auxiliary-curve construction.** Independently of E, sample nonsingular short Weierstrass curves E'/F_p with uniformly random coefficients (a',c'). The resulting isomorphism class distribution is proportional to 1/|Aut_Fp(E')| because each class has (p−1)/|Aut_Fp(E')| coefficient models. Howe's weighted curve-order estimate for ell not dividing p−1 gives probability 1/(ell−1)+O(ell/sqrt(p)) of ell dividing #E'(F_p). Here ell=p^{o(1)} and ell²=o(sqrt(p)), so this is Omega(1/ell). Count points using Schoof's polynomial-in-b algorithm until finding such a curve. This costs expected O(ell poly(b)) bit operations. Sampling nonsingular curves and rational points by rejection has constant expected overhead. This invokes the classical point-counting and curve-order results, not a new estimate of the distribution of elliptic curve orders.

Frobenius on E'[ell] has characteristic polynomial with roots 1,p modulo ell, which are distinct by our choice of ell. Neither is −1. Therefore E'[ell](F_p) is cyclic of order ell. Its ell-primary rational group is cyclic (its subgroup killed by ell has dimension one). Given N'=#E'(F_p), multiplying a uniformly random rational point by N'/ell produces a nonzero point H of order ell with probability 1−1/ell.

Define S={x([j]H):1≤j≤(ell−1)/2}. Its size is exactly w. The standard odd division polynomial psi_ell of E' has the crucial property

    S = {u∈F_p : psi_ell(u)=0}.

Indeed a geometric nonzero ell-torsion point with rational abscissa is carried by Frobenius to itself or its negative. The latter is impossible for a nonzero point because −1 is not an eigenvalue. Conversely every rational ell-torsion abscissa occurs in the list. Notice that the *full* division polynomial has degree (ell²−1)/2, not w. It has other roots over the algebraic closure. We must retain the rational-root equation U^p−U=0 when reducing it modulo another polynomial. In particular its resultant alone is **not** the scalar norm over S.

Enumerating S costs O(w poly(b)) bit operations, needed for factor-base access and public setup. The symbolic membership circuit below does not need the expanded polynomial with roots S.

References: E. W. Howe, *On the group orders of elliptic curves over finite fields*, Compositio Math. 85 (1993), 229–247 (weighted distribution theorem); R. Schoof, *Elliptic curves over finite fields and the computation of square roots mod p*, Math. Comp. 44 (1985), 483–494. Standard division polynomial identities, e.g. Silverman, *The Arithmetic of Elliptic Curves*, Chapter III, and the corresponding Kummer-line multiplication formulas.

## 2. Coverage with a shared translation, including the thinner widths

We give a fresh average-coverage argument because the previous shared-translation estimate had a different nominal-density range. Use the public set S above in each slot, with disjoint translates c_i+S. They can be chosen greedily: at stage i fewer than (i−1)w² translates can intersect a previous one. Choose a public common multiplier t∈F_p^* such that the shifted cubics f(A+t(c_i+u)), for all distinct public offsets c_i+u, have disjoint finite branch sets. For each pair of distinct offsets, a bad t is determined by a nonzero difference of two geometric roots of f; at most six choices lie in F_p. Thus O((mw)²) bad choices suffice to bound all failures, and this is o(p). Public enumeration/checking fits 2^{O(b/m)} poly(b), and the finite branch differences can be computed over the degree-at-most-six splitting field of the cubic. Sample a single uniform A∈F_p. Slot i has abscissas A+t(c_i+u) for u∈S. Keep rational lifts with nonzero ordinates. The actual factor base can be enumerated and checked within the same setup bound.

For q distinct offsets d_1,...,d_q, the multiquadratic fiber product y_j²=f(A+d_j) is geometrically integral: disjoint finite branch loci show independence of all q square classes. Its smooth projective model C has genus satisfying 2g−2=2^(q−1)(3q−3), by the tame degree-2 ramification at each of 3q finite branch points and at infinity. Its coordinate maps rho_j:C→E have degree 2^(q−1). For integer coefficients k_j not all zero, the theorem of the square and the sign involutions on C give

    deg([h]∑[k_j]rho_j) = h² 2^(q−1) ∑ k_j².

For clarity, use the symmetric degree-one line bundle O_E(O); its pullback degree is a quadratic form in maps C→E, with constant maps of degree zero. Flipping rho_j to −rho_j while fixing other coordinates preserves degree by the corresponding automorphism of C, so all mixed coefficients vanish. Multiplication by h scales degree by h², proving the identity. A nonconstant map of smooth projective curves has at most its degree many geometric points in any fiber, also in the inseparable case.

For each of W address tuples, Hasse–Weil therefore gives an average of 1+O(2^{O(m)}poly(b)/sqrt(p)) rational signed lifts over A; removing zero ordinates or A=∞ costs at most 2^{O(m)}poly(b)/p. The map to the projected group sum [h]∑P_i is nonconstant and its zero fiber is bounded by the degree formula, so excluding target O has negligible additional effect. Let T be independently uniform among the r−1 nonidentity subgroup points. Write H_T for the number of rational signed tuples summing after projection to T. Then E H_T=(1+o(1))W/(r−1).

For two distinct signed tuples, combine their union of at most 2m public offsets. Their projected sum difference is a map of the displayed type, with coefficients in {−2,−1,0,1,2}, not all zero. To bound pairs, for any fixed pair of address tuples, split possible lifts at repeated addresses according to their relative signs (at most 2^m possibilities); identical signed tuples are the diagonal and are not counted again. Each nontrivial pattern has at most 2^{O(m)}poly(b) geometric lifts over all A for which the projected sums agree. Averaging additionally over T divides by p(r−1). Summing over W² pairs proves

    E H_T² ≤ O(lambda)+2^{O(m)}poly(b) lambda².

Paley–Zygmund yields Pr(H_T>0)≥2^{−O(m)}/poly(b), since lambda=2^{−Theta(m)}. This is an *average* coverage result, not a uniform targetwise Poisson theorem. Its inverse success cost fits the intended exponent scale. Filtering nonzero rational ordinates is charged explicitly; no heuristic dismissal of mixed-twist branches is necessary.

## 3. Reduce the auxiliary support predicate in a small quotient

For the remaining m−1 coordinates let B be the product algebra of F_p-valued functions on their retained address grid. Equivalently, for their distinct-address monic polynomials Q_j, B=⊗_{j≠i} F_p[U_j]/(Q_j). The symbol U is the eliminated address. Other coordinates are restricted to abscissas with rational nonzero ordinates, using the enumerated factor base. Let F∈B[U] be the raw projected elliptic relation polynomial for target x(T). Its declared U degree d=h²2^(m−1) is an upper bound, and every lane is a nonzero polynomial. This follows from the homogeneous resultant factorization of the previous checkpoint: the unprojected sum form has primitive linear factors for its signed orbits, including at infinity; for fixed other points, the possible final points reaching a finite projected target are a finite set because [h] and translation are finite morphisms. No specialized lane can be identically zero. Individual degrees may drop.

For a monic lane of degree d'≤d, work in B[U]/(F) (or use a degree-specific artificial monic lane below). Reduce psi_ell(U) by standard normalized division polynomial recurrences using only polynomially many operations in d and b. Here are explicit polynomial formulas over arbitrary rings. Write phi_n=psi_n for odd n, phi_n=psi_n/(2y) for even n; the latter notation denotes a polynomial, not a division performed in the algorithm. Define D=4(U³+a'U+c'). The bases are

    phi_0=0, phi_1=phi_2=1,
    phi_3=3U⁴+6a'U²+12c'U−a'²,
    phi_4=2(U⁶+5a'U⁴+20c'U³−5a'²U²−4a'c'U−8c'²−a'³).

For j≥3 use

    phi_(2j)=phi_j (phi_(j+2) phi_(j−1)² − phi_(j−2) phi_(j+1)²).

For odd index 2j+1≥5, take

    phi_(j+2) phi_j³ − phi_(j−1) phi_(j+1)³,

multiplying the first term by D² if j is even, and the second term by D² if j is odd. These are polynomial identities. At each recursion depth the required indices are in an interval of constant length: the next indices are floors of halves shifted by at most two. Memoization gives O(log ell)=O(b) distinct states, with a constant number of ring operations each. Thus zero divisors and vanishing 2y do not impede the computation. Reducing in a monic quotient costs poly(d) underlying B-gates per operation. Compute U^p−U and f(A+t(c_i+U))^((p−1)/2)−1 by O(b) further quotient multiplications.

The simultaneous gcd of F and those three remainders is monic and squarefree, because U^p−U is squarefree. Its roots are precisely valid eliminated addresses with rational nonzero ordinate. The other leaves are rational already. Since T has nonzero ordinate, a signed relation with projected x-coordinate x(T) can be oriented to T by simultaneous negation. Sign enumeration for a recovered address costs 2^m poly(b).

Degree drops can be handled by the fixed masked Euclidean circuit from the fiber-gcd checkpoint. If f_k denotes the kth coefficient of F, use masks e_k=f_k^(p−1) ∏_{j>k}(1−f_j^(p−1)) and artificial monic U^k+∑_{j<k} f_j f_k^(p−2) U^j, run the degree-specific circuits, and recombine using e_k. The zero-polynomial lane, if allowed in a general input, has the full filtered support, with its moments precomputed in O(w poly(b)); actual elliptic lanes need no such fallback. All degree choices and Euclidean steps have poly(d,b) gate cost. Degree and root power sums then give local counts and moment tags using the previous Newton-identity construction. This improves the prior interval support reduction bound ~O(d+sqrt(dw)) to poly(d,b) B-operations, with no slot-width term in the symbolic circuit. Setup still depends on w.

## 4. Build the summation polynomial as a circuit over B, not a dense multivariate expansion

The earlier dense coefficient-construction bound 2^{O(m²+m log h)} is unnecessary for this interface. For any block of k leaves construct the homogeneous binary sum form of declared degree 2^(k−1) in the output coordinate, sequentially joining another leaf with S3, using the previous checkpoint's declared-degree resultant identity. Its coefficients lie in B. Compute the resultant at 2^k+1 distinct output field coordinates and interpolate the polynomial of degree at most 2^k. Each scalar evaluation is a Sylvester determinant of dimension 2^{O(k)} over B, computable by a polynomial-size division-free determinant algorithm, such as Berkowitz. Interpolation denominators are units (public distinct field elements), and the homogeneous declared degree is retained even when some specialized leading coefficients vanish. The operation count for a join is 2^{O(k)}, and summation over the joins costs 2^{O(m)}. Projection uses the degree-h² Kummer multiplication form U_h−x(T)V_h, computed in poly(h) field operations; its resultant is degree at most h²2^(m−1) in the eliminated leaf coordinate. Evaluate the projected polynomial at d+1 distinct eliminated addresses (not required to belong to the support), and interpolate all its U coefficients over B. All degree bounds are p^{o(1)} in the specified range, so enough public interpolation points exist.

Consequently construction of F, support reduction, masked gcds and a bounded number of tagged moments costs 2^{O(m)}poly(b,h,K) symbolic gates over B for a moment cap K. No dense expansion in the remaining variables occurs. This statement charges arithmetic gates, not collective trace evaluation. At m≈sqrt(b) it removes the previous m² construction term. The public setup has expected bit cost 2^{O(b/m)}poly(b,h).

## 5. A branch-free alternative: one small determinant per fiber

Here is a simpler randomized projected support circuit that needs neither Euclidean degree selectors nor a variable matrix dimension. It applies to every nonzero actual F lane, of declared degree d, on a remaining grid of N lanes. Extend coefficients to K=F_{p²}. Pick alpha uniformly in K. For each lane form

    F'(V)=V^d F(alpha+1/V).

The coefficient of V^d is F(alpha). With probability at least 1−dN/p² it is nonzero on *all* lanes, by a union bound. Normalize it by raising it to p²−2; this is an arithmetic circuit even on bad choices, which are counted as failures rather than assumed detected. In a good chart F' is monic of degree d in the product algebra. Degree drops in F are represented as roots V=0, not discarded.

In the quotient by F', write X=1+alpha V, Z=V, corresponding to U=X/Z. Compute the homogeneous normalized division polynomial Psi_ell(X,Z) via the same recurrences, replacing D by 4(X³+a'XZ²+c'Z³) and homogenizing the bases. The degrees (n²−1)/2 for odd n and (n²−4)/2 for even n give exactly homogeneous recurrences. Define

    P1=Psi_ell(X,Z),
    P2=X^p−XZ^(p−1),
    L=(A+t c_i)Z+tX,
    P3=(L³+a LZ²+cZ³)^((p−1)/2)−Z^(3(p−1)/2).

All operations take place modulo the monic F'. At V=0, P1=Psi_ell(1,0)=ell, a unit. At finite points the three equations P1=P2=P3=0 select precisely the rational auxiliary support and rational nonzero original ordinates. Pick eta uniformly in K and form P=P1+eta P2+eta²P3. Let M_P denote multiplication by P in the d-dimensional monic quotient, and define Delta=det M_P over B⊗K.

For every positive lane, Delta=0 for all eta. For every negative lane, each of the at most d geometric roots of F' gives a nonzero degree-at-most-two polynomial in eta (its coefficients are not all zero). The determinant vanishes iff P and F' have a common root, including nonreduced quotients, because a multiplication map is invertible iff the multiplier is a unit. Thus at most 2d values of eta can create a false positive in each negative lane. Including the chart failure, probability of any error on the full remaining grid is at most 3dN/p². Since N≤w^(m−1)≤r=O(p), d=2^{O(m)}poly(b), this is 2^{−b+O(m+log b)}, negligible compared with the average coverage in §2. The indicator is 1−Delta^(p²−1). Its circuit construction has 2^{O(m)}poly(b,h) B⊗K-gate cost, using a polynomial-size division-free determinant (the matrix entries are quotient multiplications of P by 1,V,...,V^(d−1)). The test implementation uses lane-wise Gaussian elimination only to check the identity; it does not implement the asymptotic division-free circuit.

If a projected positive lane is located, factor its original polynomial F of degree ≤d over F_p, apply all support/liftability tests, and check the signed group sum. This is polynomial in d,b times the sign enumeration. If the projected grid has exactly one positive lane, its count trace and the m−1 coordinate-tag traces recover it. Multiple projected lanes require the previously derived bounded-count moment recovery; no automatic isolation guarantee beyond the average coverage above is asserted here.

## 6. Representation barrier that survives the new sampler

For an unfiltered set S of w distinct elements of F_p with w<p, consider the cardinal polynomial L_s(U−A) for the translated grid A+S. Its coefficient of U^j is a polynomial in A of degree w−1−j with nonzero leading coefficient binom(w−1,j)(−1)^(w−1−j)/Q'(s), as p>w. Union over w choices of s and all j shows that, with failure probability at most w²(w−1)/(2p), every cardinal polynomial has every coefficient nonzero. A union over axes works also with a common shift. A nonzero singleton indicator in the ordinary reduced polynomial basis on this *unfiltered shifted grid* is the product of cardinal polynomials, hence has full support with high probability. This is a statement in translated abscissa coordinates, not raw auxiliary-address coordinates, and not a theorem about coefficients after liftability filtering. Maximal individual degree for a singleton holds in any distinct-address grid. It is a barrier to dense or sparse expanded monomial representations, not to arithmetic circuits or fast trace evaluation.

## 7. Actual cost and precise missing theorem

The new setup is 2^{O(b/m)}poly(b,h) expected bit operations. The complete single-fiber projected indicator circuit has 2^{O(m)}poly(b,h) gates *over the remaining product algebra*. These are two different complexity measures. If the B values are explicitly stored, their lane count can be w^(m−1), and computing the trace by scanning lanes restores the exponential cost. Neither the sparse membership circuit nor the small determinant is a fast algorithm for the product over all lanes or for their collective trace. Even a cheap arithmetic-circuit description of a function does not imply that summing it over a Cartesian grid is cheap.

The needed result is an algorithm for the collective tagged traces of these specific determinant-zero circuits, without enumerating the remaining lanes, with total bit cost 2^{O(b/m+m)}poly(b,h) (or any bound strictly beating the current asymptotic scale after optimization). At m≈sqrt(b), that would fit the square-root bitlength log-time scale, conditional on a full relation-reuse and outer-reduction audit. No such algorithm is established. No exponent reduction is claimed.

## 8. Executed exact checks and reproducibility

Run from this directory: `python verify.py`, `python verify_chart.py`, `python verify_nonrational.py`. The file `previous_fiber_verifier.py` is the prior checkpoint's helper copied for reproducibility. All scripts use the standard library and exact finite-field arithmetic; the finite fixture point counts use exhaustive enumeration, not Schoof, and the determinant tests use small direct algorithms. No collective trace evaluator is benchmarked.

- Seven auxiliary curves over p=101,...,4093: 8,099 abscissa checks, 504 group-multiplication comparisons, 130 curve proposals; normalized division-polynomial recurrence also tested with 512-bit indices. These are not 512-bit-field benchmarks.
- 399 masked gcd lanes, including 270 degree drops, 50 zero-polynomial lanes, repeated roots, and 2,394 root-moment comparisons; the same 399 lanes tested with original-curve square filtering.
- 32 actual elliptic relation grids, 422 address checks, 24 positive controls, 13 projected cases. The initial random fixtures happened to give no positive cases; the updated run includes 16 planted controls and all final assertions passed. Initial and intermediate transcripts are retained.
- 420 homogeneous recurrence checks; 313 checks of the odd-index leading term at infinity; 316 projected chart fibers and 3,792 randomized determinant checks, including 1,200 true-positive comparisons and zero false positives in that run. Another 4,225 eta values exhaustively checked the error bound on small fields, with 22 negative-fiber root-bound checks.
- 40 product-ring summation-polynomial circuit fixtures: 720 scalar comparisons and 1,062 coefficient lanes; toy division-free determinant implementation.
- Deliberate nonrational-root negative control over F_13 with ell=5: the unfiltered division polynomial falsely accepts an irreducible quadratic. The filtered combined determinant has 3 bad eta values among 169, within its rigorous 2d=4 bound. The randomized test is not perfectly sound at small p, as the theorem predicts.
- 106 coefficient-root-bound checks, 31,660 translation evaluations, and four union-bound fixtures for the expanded singleton indicator.

All final exact component assertions passed. Transcripts and checksums accompany the scripts. These checks validate identities and edge cases on finite fixtures, not the asymptotic setup rate or ECDLP runtime.
