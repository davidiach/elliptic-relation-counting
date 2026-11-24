# Block norm elimination, one-call witness extraction, and the grid lifting ideal

Research checkpoint, 24 November 2026. Continuation of `shared_norm_aggregate_note.md` (§§4–5) and `shared_translation_note.md`. This checkpoint gives (i) a division-free block norm algorithm, (ii) exact recovery of a bounded set of zero addresses from one perturbed ring-valued norm, (iii) an obstruction to standard small-root lifting for the full summation polynomial, (iv) the exact index and low-degree information of the grid quotient ideal, and (v) a parameterized hardness comparison for general succinct Kronecker determinants. It does not solve the scalar structured norm evaluation problem, beat the general-prime exponent, or independently audit the earlier conditional coefficient 11.31486814. No novelty claim for the classical tools or exhaustive priority search is made.

## 1. Exact problem and assumptions

Let E/F_p be the nonsingular curve in short Weierstrass form, p>3, N=hr with r odd prime, gcd(h,r)=1, and h=poly(b), b=ceil(log2 p). For the cofactor statements below take p>h, which holds eventually in this regime. Fix a nonidentity subgroup target T, so 2T≠O. Let F be the universal projected fixed-target abscissa relation polynomial R_{m,h}(A+U_1,...,A+U_m,x(T)) from the previous checkpoint, with scalar normalization immaterial for its zero addresses. The supports S_i consist of distinct field elements; w_i=|S_i|, W=prod w_i. Define B=⊗_i F_p[U_i]/(Q_i(U_i)), where Q_i(U)=prod_{s in S_i}(U-s). CRT identifies B with F_p^W and identifies multiplication by F with the diagonal operator with eigenvalues F(s) for grid points s. Its norm N(F)=prod_s F(s) is its determinant. A zero is an algebraic abscissa relation; enumerating signs and verifying [h] times the rational sum recovers the genuine relations. The existing reduction's coverage and identical-selection-law assumptions remain necessary.

## 2. Exact block norm algorithm (one slot)

**Theorem 1.** Let R be any commutative F_p-algebra, F(U)∈R[U] have degree at most d, 1≤d<p, and 1≤w≤p. For a,t∈F_p, t≠0, the norm over {a,a+t,...,a+(w-1)t} can be computed in soft-O(d+sqrt(dw)) R-operations. For d=0 use ordinary powering. This counts operations in R, not bit time for an arbitrary succinct representation of R.

Proof. Put f(X)=F(a+tX). Let H=max(1,ceil sqrt(w/d)), clipped at w. Form B_H(X)=prod_{j=0}^{H-1} f(X+j). Its degree is at most dH. The norm is the product of B_H(qH) over q=0,...,floor(w/H)-1 and the remaining fewer than H scalar f-evaluations. The affine substitution and each Taylor shift are fast polynomial convolutions: since d<p, the factorials through d are invertible field scalars. Explicitly, coefficients after a shift c have formula [X^k]f(X+c)=1/k! sum_{i>=k} f_i i! c^{i-k}/(i-k)!. Reversing one array makes this convolution. Product trees build B_H, and a subproduct/remainder tree evaluates it at the required qH. Monic polynomial divisions use only inversion of a reversed series with constant coefficient 1, which is valid in R without any assumption that nonzero entries of R are units. Tail evaluations and products fit the same cost. Fast polynomial multiplication over arbitrary algebras (Cantor–Kaltofen) and the standard monic remainder tree give soft-O(dH+w/H+d) ring operations, or soft-O(d+w) when d>w. Optimize H as above. No inverses of F(s) are used. This is an application of the classical baby-step/giant-step factorial technique, not a claim of a new general factorial algorithm. ∎

References: A. Bostan, P. Gaudry, E. Schost, *Linear Recurrences with Polynomial Coefficients and Application to Integer Factorization and Cartier–Manin Operator*, SIAM J. Comput. 36 (2007), https://doi.org/10.1137/S0097539704443793. D. G. Cantor, E. Kaltofen, *On fast multiplication of polynomials over arbitrary algebras*, Acta Informatica 28 (1991), https://doi.org/10.1007/BF01178683. FLINT documentation, https://flintlib.org/doc/algorithms/polynomials/subproduct-tree.html (monic remainder evaluation over arbitrary commutative rings).

For d=h²2^(m-1), log w=b/m+O(log b) at balanced widths, the first slot norm needs at most 2^(b/(2m)+O(m+log b)) coefficient-ring operations, provided d≤w and d<p. This is an exact local speedup, not an end-to-end complexity bound. When d>w the bound is instead soft-O(d+w), still covered by the displayed exponential upper bound with appropriate O(m+log b).

### Degree growth in successive eliminations

For the unprojected fixed-target relation, let P_{m,T} be its irreducible projective image equation. Its separate affine degrees are d=2^(m-1), by fixing generic abscissas of m-1 points: the sign choices force 2^(m-1) distinct final abscissas (2T≠O). The image of the sum fiber E^(m-1) under the finite abscissa map is irreducible and closed in (P^1)^m. For m=1, use x_1-x(T). The image omits the point with all abscissas infinity, since each such point is O. Hence the corner monomial prod x_i^d has nonzero coefficient. Semaev's symbolic resultant recurrence, with declared homogeneous degrees, specializes to a nonzero scalar multiple of P_{m,T}: the specialized polynomial vanishes on the image, is nonzero at a generic tuple outside it, and has the same separate degree upper bounds. This uses the generic partial-sum argument and degree bounds of the previous checkpoint. In particular there are no hidden extra factors at the chosen fixed T.

Taking the leading coefficient in x_j is specialization of the multihomogeneous equation at x_j=∞. The restricted zero set is exactly the lower-arity image for target T, since the j-th point is O. The restriction is nonzero at the remaining all-infinity point. It has multidegree d in every other coordinate, whereas the lower-arity image equation has multidegree d/2. Thus it is a nonzero scalar multiple of P_{m-1,T}² for m≥2. This also follows by factoring the homogeneous polynomial: any nonconstant factor defines a divisor on the product of projective lines, and the only allowable divisor is the lower-arity image. No coordinate-at-infinity factors can be inserted. The corner coefficient is nonzero for every fixed eligible T, not just a generic target.

For the cofactor case h<p, the preimages [h]^{-1}(T) are h² distinct points over the algebraic closure, none of order dividing two. Their abscissa images are distinct: for a generic tuple on a fixed sum fiber, an independent nonconstant point prevents mixed sign flips from changing its sum to a different fixed point; global negation changes T to -T. Equivalently the two fixed-sum images agree only for identical or opposite sums. Products over these h² preimages (with the equivalent ±T fibers identified by global negation) give the projected relation polynomial. The symbolic projected resultant has the same separate degree bounds h²2^(m-1), is nonzero generically, and contains these factors; thus they exhaust it. The same corner and leading-coefficient-square identities hold with projected degrees. Affine substitutions x_i=A+U_i, and nonzero scalar slopes on coordinates, preserve the corner and the degrees.

Suppose k slots have been norm-eliminated, with at least two slots still free. Each individual specialized polynomial in the product has degree exactly d in each remaining coordinate, as a polynomial over the other free coordinates. Its leading coefficient is the lower-arity projected relation polynomial squared, which cannot vanish identically after the eliminated slots are fixed: one remaining free variable other than the variable being eliminated from the leading term can only take finitely many values satisfying a prescribed sum, after fixing all other variables and enumerating algebraic signs. Multiplying the factors gives exact remaining separate degrees

    d * prod_{eliminated i} w_i.

In particular a pipeline that eliminates slots sequentially and explicitly expands a dense coefficient vector in each next variable must pay 2^(b/2-o(b)) positions already at k=floor(m/2) (for growing m and balanced widths, ignoring lower-order terms). This bound applies to that representation, not to sparse or implicit circuits, other orderings of fundamentally different eliminations, or the time of an unrestricted norm algorithm. The square-root-in-width first step does not prevent this growth.

## 3. Recover every low-multiplicity witness from one norm (new transfer identity)

**Theorem 2.** For any polynomial F on a Cartesian grid S, let H count its distinct zero addresses, and choose K≥0. There is one perturbed norm call over a ring of dimension (m+1)(K+1)^2 over F_p which recovers all zero addresses if H≤K and a random linear projection separates them. It distinguishes H=0 from H>K. Conditional on H≤K, the projection fails with probability at most H(H-1)/(2p).

Choose a uniform linear form ell(U)=sum_i c_i U_i. Let

    R = F_p[Z,Y,eta_1,...,eta_m]/(Z^(K+1), Y^(K+1), eta_i eta_j for every i,j).

Evaluate the collective norm of the multiplication polynomial

    F(U) + Z (Y - ell(U) - sum_i eta_i U_i)

over this ring, with the same slot grid. This requires an evaluator for the perturbed polynomial over R, not just a scalar zero-test oracle specialized to the unperturbed elliptic polynomial. It adds only m+1 linear or constant terms, but the evaluator's ring and polynomial input requirements must explicitly allow this.

By CRT the norm is the product over all addresses. Each zero address contributes a factor Z; each nonzero address has nonzero constant term F(s). If H>K the truncated norm vanishes identically. If H≤K, its first nonzero Z coefficient is at order H, and is

    P(Y,eta) = C prod_{F(s)=0} (Y-ell(s)-sum_i eta_i s_i),
    C = prod_{F(s)≠0} F(s) ≠0.

Y-degree at order j in Z is at most j, so the Y truncation cannot affect this formula. H=0 yields a nonzero constant and no addresses. For H>0 write P=P_0(Y)+sum_i eta_i P_i(Y). The coefficient of Y^H in P_0 is C, hence the first nonzero Z order is detectable even if some projection value is zero. If ell is injective on the zero set, P_0 splits into H distinct linear factors. Factor it over F_p. For every root y=ell(s), recover the individual coordinates by

    s_i = - P_i(y)/P_0'(y).

Indeed the epsilon derivative on that root differentiates exactly its own linear factor; the other terms vanish. This requires only ordinary univariate polynomial factorization and field operations, expected poly(b,m,K) time and space. Verify that reconstructed residues lie in S_i and satisfy F, and verify signed rational group sums separately. For distinct grid addresses s,s', the uniform ell(s-s') is uniform in F_p. Union bound over pairs proves the probability assertion. Detect a repeated root by gcd with the derivative, and retry or abstain. The formula is exact for arbitrary eigenvalue multiplicities of the nonzero values: no characteristic-polynomial or full trace recurrence is needed. ∎

For H=1, a smaller square-zero ring suffices: use F(U)+eta_0+sum_i eta_i U_i. The norm constant distinguishes no zero, and its linear coefficients are C and Cs_i on a singleton. If there are ≥2 zero addresses all linear coefficients vanish. This obtains the complete address from one norm evaluation, with m+2 base-field dimensions, rather than O(b) separate zero queries. It does not itself give an efficient norm evaluator.

### Transfer to shared-translation average coverage

For each address and any translation, at most 2^m nonidentity subgroup targets can be algebraically reachable: the ordinate signs determine at most 2^m projected sums (global negation pairs targets). Thus for independent uniform A,T,

    E H_alg ≤ W 2^m/(r-1).

Write κ for the previous checkpoint's lower bound on probability of a genuine hit, κ≥2^(-O(m))/poly(b,λ^-1). At balanced λ=W/r=poly(b)^(-1)..O(1), choose K≥4 W 2^m/((r-1)κ), rounded up. Markov gives Pr(H_alg>K)≤κ/4. For m≤sqrt b, log K=O(m+log b) and K²/p=o(κ). On a nonoverflowing true positive, Theorem 2 recovers all algebraic addresses with probability 1-o(1). Enumerating 2^m signs and verifying each recovered address preserves success probability Ω(κ). Mixed-twist or exceptional algebraic addresses can be rejected without needing a separate false-positive bound for the evaluator. The ring has 2^O(m+log b) coefficients; even schoolbook arithmetic is within this exponential scale. If a uniform exact norm algorithm for these perturbations had time 2^O(b/m+m²) times poly(ring dimension,b), this would preserve the b^(2/3) sufficient objective. If its time were 2^O(b/m+m) times the same factor, the square-root objective would follow, subject to the inherited outer reduction. These are still sufficient hypotheses, not obtained algorithms. A scalar-only algorithm that depends on special properties of the unperturbed relation polynomial need not extend to this perturbation.

## 4. Full small-root lifting obstruction for the actual summation polynomial

**Theorem 3.** Let F∈Z[X_1,...,X_m] be monic at its componentwise maximal monomial X^d=prod_i X_i^{d_i}, all d_i>0, with all other exponents componentwise ≤d. Put I=(p,F). For any integer t≥1 and bounds B_i≥1 satisfying prod_i B_i^{d_i}≥p, every nonzero g∈I^t satisfies

    ||g(B_1 X_1,...,B_m X_m)||_∞ ≥ p^t.

The same holds for its coefficient l_1 and l_2 norms. In particular no polynomial in this ideal meets the ordinary strict coefficient-norm small-root lifting certificate at these bounds. This does not require a restriction to triangular shift sets or any estimate of the lattice determinant.

Proof. In the integer polynomial ring, the family of polynomials

    X^beta F^j,    j≥0, beta not componentwise ≥d

is a Z-basis. Indeed every exponent alpha has the unique decomposition alpha=beta+j d with j=min_i floor(alpha_i/d_i). The corresponding polynomial has leading monomial X^alpha with coefficient 1 in any monomial order, and all other terms are smaller. Triangular reduction, terminating by the monomial well-order, proves spanning and independence. Multiplication by F increases j by one. It follows that the coefficient of a basis member with index j in any g∈I^t is divisible by p^max(t-j,0): I^t is generated by p^(t-k)F^k for 0≤k≤t, and arbitrary polynomial multipliers can be expressed in this same basis. Conversely such divisible coefficient arrays give elements of I^t. For the largest nonzero monomial exponent alpha of g, triangularity shows its coefficient equals that basis coefficient. Thus its scaled absolute value is at least

    p^max(t-j,0) prod_i B_i^(beta_i+j d_i) ≥ p^t.

This proves the result, including arbitrary cancellations among shift generators and arbitrary integer linear combinations. ∎

Application: Section 2 proves that the actual fixed-target summation polynomial (also projected, under the cofactor assumptions) has a nonzero full corner monomial with separate degrees d_i=d. Normalize it over F_p, then lift with corner coefficient exactly one over Z. The ideal (p,F) is independent of the integer lift and of multiplication of its residue by a nonzero field scalar. At balanced interval widths, even using centered bounds B_i roughly w_i/2, log prod B_i=b-O(m+log b), hence prod_i B_i^d≥p for d≥2 and sufficiently large b. A modular root of F makes all polynomials in I^t vanish modulo p^t, but no nonzero polynomial from this ideal can meet the strict small-coefficient criterion. In particular simply generating more monomial shifts of p^(t-j)F^j and using a better lattice basis or exact shortest-vector computation does not bypass this norm obstruction. This concerns the stated coefficient norm, not all methods for small modular roots: adding grid equations, auxiliary algebraic relations, changing variable representation, or using a different way to certify small evaluations lies outside it. It is not a runtime lower bound or a claim that multivariate Coppersmith methods cannot work in a lifted system. For background on the usual shift-polynomial and determinant analysis, see K. Ryan, *Solving Multivariate Coppersmith Problems with Known Moduli*, EUROCRYPT 2025, https://doi.org/10.1007/978-3-031-91095-1_13.

## 5. Add exact grid equations: a genuinely different ideal

Let the integer supports be the canonical distinct representatives of S_i and Q_i be their monic products over Z. Form the free rank-W algebra

    A_Z = Z[U_1,...,U_m]/(Q_1,...,Q_m),   J=(p,F) in A_Z.

**Theorem 4.** If H=|{s∈grid : F(s)=0 mod p}|, then for every t≥1,

    [A_Z:J^t] = p^(t H).

Proof. Modulo p^t the evaluation matrix on the standard monomial basis is a tensor product of Vandermonde matrices. Its determinant is a p-unit, since all nodes are distinct modulo p in each slot, so evaluation is an algebra isomorphism onto (Z/p^t)^W. At any address where F(s) is a unit, F(s)^t generates the full coordinate ring. At any other address, every generator p^(t-j)F(s)^j of J^t vanishes modulo p^t. As p^t A_Z⊂J^t, the quotient has exactly H independent coordinates of size p^t. ∎

In particular the quotient ideal has index p on a singleton, despite the full-grid rank W≈p. Reducing modulo p, its image ideal is exactly all functions vanishing on the zero set. With w_i≥2, the affine-linear monomial subspace V=span(1,U_1,...,U_m) has dimension m+1. On a singleton s*, J mod p intersect V is the hyperplane

    c_0 + sum_i c_i s_i*=0.

A basis of this intersection suffices to reconstruct s* by linear algebra: its annihilating row, normalized at the first coordinate, is (1,s_1*,...,s_m*). Equivalently, extracting the ideal's low-degree part would solve relation search without reconstructing the full determinant. This is a precise computational target, not a proven efficient ideal-membership algorithm. The compact description of J does not supply its intersection with V for free.

There is a complementary degree obstruction for algorithms that try to build the right kernel of multiplication by F in the standard monomial basis. On a singleton its kernel is one-dimensional and spanned by the tensor-product Lagrange polynomial

    L_s(U)=prod_i [Q_i(U_i)/((U_i-s_i) Q_i'(s_i))].

CRT proves this directly. Its leading corner coefficient is nonzero and its total degree is sum_i(w_i-1), the maximum total degree in the standard monomial basis. Hence every nonzero representative in the kernel has this full degree. An algorithm searching for a kernel vector within total-degree truncated subspaces of that basis sees no kernel until the truncation includes all W basis monomials. This is not an obstruction to finding the factorized Lagrange form directly or to algorithms using a different basis. The degree-one **image ideal** and the high-degree **kernel** are different computational objects. The norm-perturbation theorem computes the same zero addresses by a different route.

## 6. Generic Kronecker determinant evaluator would solve parameterized clique

**Theorem 5 (scope warning, not an elliptic lower bound).** A general algorithm deciding whether multiplication by an explicitly supplied polynomial on a product of arbitrary split monic slot algebras is singular in time f(m)(w b)^C, where C is absolute, would give an FPT algorithm for k-clique. This holds for polynomials of total degree two with only O(m) separable monomials, equal slot widths w, and pairwise disjoint supports. With randomized prime generation it also holds, as a randomized FPT consequence, under the balanced promise b=m log2 w+O(log(k+2)). A randomized consequence continues to hold for an evaluator guaranteed only on instances with zero or one singular grid address.

Proof. Use k≥3 color classes with n vertices each; ordinary clique reduces by making k copies of its vertices and allowing cross-color edges corresponding to distinct adjacent originals. Set m=k+binom(k,2), w=n². For each color assign distinct positive integer vertex codes in [L], and introduce a variable V_i ranging over those n codes, padded to w with distinct dummy codes in [L+1,L+w]. Set B>L+w. For each i<j introduce an edge variable E_ij whose allowed set consists of u+Bv for each allowed cross-color edge with codes u,v, padded to width w by distinct dummies in [B²,B²+w-1]. Since all vertex codes are <B, base-B uniqueness makes E_ij=V_i+B V_j possible exactly for allowed edges (including the requirement that the vertices not be dummy). Every vertex participates in at least one edge. Define the quadratic

    F = sum_{i<j} (E_ij - V_i - B V_j)^2.

Each residual is an integer of absolute value O(B²). Choose p>8mB⁴, so the integer sum is less than p and F≡0 mod p iff each consistency condition holds. Clique choices are in bijection with grid zeros. Disjoint supports are obtained by adding separate public offsets κ_i to the variable supports and substituting the shifted variable U_i-κ_i in F. Offsets of size O(mB²) fit below p; they do not change the polynomial's degree or its O(m) monomial bound. For the ordinary deterministic FPT reduction choose L=kn, B=L+w+1, and search deterministically for a prime just above the polynomial-size lower bound (Bertrand and polynomial-time primality suffice, even a full scan here is polynomial in n with an absolute exponent). Multiplication by F is a sum of O(m) tensor products of companion-matrix powers of degrees 0,1,2. Singular iff the original graph has a clique. Input construction is f(k)n^O(1). The standard W[1]-completeness of k-clique then gives the FPT implication (not an unconditional contradiction).

For the balanced version choose p randomly from primes between N0 and 2N0, where N0=max(w^m,8mB⁴+1). Random prime search with exact primality testing has polynomial expected cost in log N0, by the prime number theorem. With k≥3, B=poly(k)n² at most in the versions here, N0≤poly(k)w^m, so b=m log2 w+O(log(k+2)). As n grows for fixed k, N0=w^m. This prime-generation step is randomized; do not claim a deterministic FPT reduction under the balanced promise from that step.

For the zero/one promise assign independently random codes in [L], L=8(kn)², and reject a collision among any of the kn actual vertex codes. Collision probability ≤ binom(kn,2)/L. The isolation lemma gives a unique minimum weight clique under weight sum of its vertex codes, if any clique exists, with probability at least 1-kn/L before collision rejection. To see this, fix all weights except one vertex. The minimum weights among sets with and without that vertex tie for at most one value of its weight. A nonunique minimum implies such a tie for some vertex. Union bound proves the assertion. With constant probability both distinctness and isolation hold. Add to F the term (sum_i V_i-τ)^2 and try τ=1,...,kL. The same p bound, enlarged by an absolute constant if needed, ensures no modular wrap, because the vertex codes and τ have magnitude O(kB). This adds only O(k²)=O(m) monomials, preserves total degree two, and yields exactly one zero for the isolated minimum weight. A no-instance always yields zero roots for every τ, so the promised evaluator never produces a false positive there. On a yes-instance with nonisolated queries the evaluator may answer arbitrarily; the isolated query suffices on the successful random trial. If its running time is guaranteed only on promised inputs, cap each call at its promised bound. Repeat independently for constant success. This is a randomized FPT algorithm if the hypothesized promise evaluator exists. ∎

The standard complexity-theoretic statement about CLIQUE is W[1]-completeness with parameter k (see e.g. https://iuuk.mff.cuni.cz/~koutecky/paralg2627-graph/lecture-2-proving-the-hard-side.html). The balanced and isolated variants here give randomized FPT implications, not deterministic class equalities or unconditional lower bounds. In particular this argument does not rely on ETH, and it says nothing against a specialized algorithm for the elliptic summation polynomial, whose slot supports and coefficients are not arbitrary clique inputs. It shows why low Kronecker rank, low polynomial degree, balanced widths, and even singleton promises alone are insufficient structural explanations for a generic FPT evaluator.

## 7. Exact component regression

From this package directory run (assertions enabled, Python with SymPy):

    python -B verify.py
    python -B verify_perturbation.py

The tests use schoolbook polynomial arithmetic and exact direct enumeration on small grids. They verify identities, not the asymptotic speed of fast multiplication or a cryptographic-scale forward solver. The actual Semaev tests use independent declared-degree resultant determinants and modular interpolation, not a symbolic assumption about the leading coefficient. The integer grid ideal index is checked by exact Smith normal form of its generator matrix. The perturbed norm verifier explicitly multiplies all factors in the truncated ring, then recovers coordinates; its tiny-field factoring uses enumeration only for regression. The clique tests exhaust small grids, with optional specified total weight, against the direct multicolored clique predicate. They do not measure a parameterized hardness result experimentally.

Recorded results: 600 block-norm comparisons, including 265 nonzero zero-divisor outputs and 219 tail cases; 400 nonzero lifting-ideal elements, 4,710 F-adic coefficient checks and 400 cancellation controls; 180 exact Smith-normal-form index checks, including 25 singleton annihilators, 9 empty and 26 multiple-zero base fixtures; 216 leading-coefficient comparisons on Semaev polynomials, 27 collective degree/leading comparisons, 9 curve/arity checks; 44 complete clique-grid comparisons over 180,224 addresses (23 singleton, 4 empty); 340 full truncated norm checks including 100 actual summation-polynomial cases, 95 successful recovery fixtures with 179 recovered zero addresses, 113 empty, 124 overflow and 8 detected projection collisions. All assertions passed. The finite tests do not separately benchmark the cofactor extension of the degree proof.

## 8. Frontier and limitations

The strongest positive new route is the perturbation identity: it turns a ring-valued collective norm into a complete bounded-multiplicity recovery oracle with only 2^O(m+log b) ring overhead at the intended shared-translation parameters. It does not require enumerating subgrids or perfectly isolating true hits. The missing computation is still the actual collective norm over that ring, or the degree-one intersection of the grid quotient image ideal. We have a fast individual slot elimination but cannot yet compose it at the necessary cost; expanding its intermediate coefficients is too expensive. The full lifting obstruction is rigorous within its coefficient-norm model and the grid ideal shows exactly where that model ceases to apply. The clique comparison prevents treating a generic small Kronecker description as a sufficient algorithmic explanation. Further work should target special addition-polynomial identities that act on the grid quotient or the perturbed norm without explicit expansion. No improved end-to-end exponent or conditional constant is claimed.
