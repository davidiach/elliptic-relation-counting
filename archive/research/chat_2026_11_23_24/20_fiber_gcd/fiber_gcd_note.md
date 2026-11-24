# Fiber-gcd compression, trace-only witness recovery, and sharp image-certificate degrees

Research checkpoint, 24 November 2026. Continuation of `kronecker_determinant_note.md` §§2, 3, 5 and `shared_norm_aggregate_note.md` §§2–5. The aim is a strict asymptotic improvement over the conditional general-prime log-time scale `b/(log b)^(1/4)`, not a constant change. This checkpoint gives a deterministic, small-fiber construction and an exact transfer to trace evaluation, plus lower bounds on explicit image-certificate methods. **It does not establish the missing fast trace evaluator or a faster ECDLP algorithm.** The earlier 11.31486814 coefficient remains an unrevalidated conditional candidate. No literature-priority or independent review claim is made.

## 1. Setting and scope

Let `p>3`, `b=ceil(log2 p)`, `E/F_p` a nonsingular short Weierstrass curve with order `hr`, `r` odd prime, `gcd(h,r)=1`, `h=poly(b)` and `p>h`. Let `T` be a nonidentity point of the subgroup of order `r`. Use the projected fixed-target relation polynomial `F` of the preceding checkpoints (separate degrees `d=h² 2^(m-1)`), whose algebraic zeros represent the abscissas of signed tuples with projected sum `T` or `-T`. This checkpoint uses that earlier relation-polynomial result rather than independently reauditing its cofactor extension. A recovered algebraic address must still be tested for rational lifts and verified using group addition. All results on finite product algebras below are independent of the elliptic interpretation.

The grid is `S=∏ S_j`, with `w_j=|S_j|`, `W=∏w_j`. Let `Q_j(U)=∏_{s∈S_j}(U-s)` and `A=F_p[U_1,...,U_m]/(Q_1,...,Q_m)`. CRT identifies it with the product of `W` base fields. Fix a slot `i` whose support is an affine interval `a+t{0,...,w-1}`, `t≠0`, `w≤p`; apply the public affine change of variable so its index is `U∈{0,...,w-1}`. The other coordinates form `B=⊗_{j≠i} F_p[U_j]/(Q_j)` with `W_-=W/w` lanes. The polynomial in the index coordinate is `F(U)∈B[U]`, of degree at most `d<p`. We never assume `B` is a domain. A ring operation in a succinct description of `B` is NOT asserted to cost polynomial bit time in `b`.

## 2. Monic remainder compression

**Theorem 1.** For a commutative `F_p`-algebra `B` and a monic `f∈B[U]` of degree `1≤k<p`, the remainder of the falling factorial `Q(U)=∏_{j=0}^{w-1}(U-j)` modulo `f` can be computed in soft-O(`k+sqrt(kw)`) operations in `B`, using fast polynomial multiplication over algebras and monic division.

Proof. Take `H=min(w,max(1,ceil sqrt(kw)))`, `J=floor(w/H)`, and `P_H(X)=∏_{j=0}^{H-1}(X-j)`. For each `q=0,...,J-1`, define the monic polynomial `f_q(X)=f(X+qH)`. Build the subproduct tree of the `f_q` and compute all remainders `R_q(X)=P_H(X) mod f_q(X)` by a remainder tree. Substituting `X=U-qH` gives

    R_q(U-qH) ≡ P_H(U-qH) (mod f(U)).

Multiply these remainders modulo `f` and include the last fewer than `H` factors. The two shift operations cost soft-O(k) per block by factorial convolution (factorials through k are units). Building `P_H`, the monic tree and the remainders costs soft-O(H+kJ+k). No relative primeness between the `f_q` is needed, and their coefficients may be zero divisors. Monic remainder division is valid by inversion of a reversed series with constant coefficient 1. Multiplying all the reduced factors costs soft-O(kJ), and forming/reducing the final tail costs soft-O(H+k). Balance `H` with `kw/H`; when k>w use direct construction of Q and division. ∎

This differs from the preceding block norm method: there we multiplied evaluations of a fixed degree polynomial over an interval; here we compute the interval's entire annihilating polynomial *modulo a short monic polynomial*. In particular for monic F,

    B[U]/(F,Q) ≅ (B[U]/F)/(Q mod F).

The auxiliary free module dimension over B is k, instead of w. The underlying small matrix is multiplication by the remainder modulo F. At a fixed remaining-grid lane its cokernel dimension equals `deg gcd(F,Q)`, which is exactly the number of distinct roots in the interval because Q is squarefree. This is a local fiber compression; it does not compress the number W_- of lanes. The complexity uses the classical fast multiplication/monic remainder infrastructure referenced in the preceding note (Cantor–Kaltofen; product/remainder trees and fast factorial methods of Bostan, Gaudry and Schost). The finite verifier uses schoolbook multiplication and an actual monic remainder tree; it verifies the identities, not fast-arithmetic runtime.

## 3. Deterministic masked Euclidean algorithm on all fibers

**Theorem 2.** Given F of degree at most d, one can construct a fixed straight-line circuit of size

    soft-O(poly(d,b) * (1+sqrt(w)))

in B-operations that represents the monic gcd of F with Q on every nonzero-polynomial lane. The coefficient arrays have length at most d+1. A separate idempotent marks lanes on which F is the zero polynomial. In particular no uniform leading-coefficient or degree assumption across the grid is necessary.

For a coefficient vector f=(f_0,...,f_d), define

    e_k = f_k^(p-1) ∏_{j>k}(1-f_j^(p-1)),
    e_zero = ∏_{j=0}^d(1-f_j^(p-1)).

These are pairwise orthogonal idempotents summing to 1 in B: by CRT, they identify the true degree of each lane, or the zero polynomial. For every k≥1 form the artificial monic polynomial

    f^(k)=U^k+∑_{j<k} f_j f_k^(p-2) U^j.

It equals the monic normalization of F wherever e_k=1. Its behavior elsewhere is harmless because the final answer is multiplied by e_k. Compute `Q mod f^(k)` via Theorem 1. The k=0 nonzero constant lane has gcd 1.

To compute gcds in parallel, start with monic `A=f^(k)`, `C=Q mod f^(k)`, pad to degree bound k. In each of k+1 iterations, compute degree idempotents of C. If C is zero on a lane, retain A and set next C=0. For every possible nonzero degree l form the artificial monic normalization C^(l) and compute `A mod C^(l)` by monic long division. Combine the next pair as `(C^(l),A mod C^(l))` on its mask and `(A,0)` on the zero mask. Every lane either has finished or its second degree decreases strictly. The result is the usual monic gcd in every lane after at most k+1 iterations. Each fixed iteration costs poly(k,b) B-operations, using binary powering for all Fermat masks and inverses, and no nonunit divisions. Sum the resulting coefficients across the original e_k masks. On the e_zero lanes the gcd is Q itself, represented by that flag, not by a d-degree array. This proves the circuit bound, including the cost of all d possible original degrees and all degree branches per iteration. ∎

On actual elliptic fibers the zero-polynomial case does not occur: after all other finite abscissas are fixed over the algebraic closure, each sign choice gives a fixed partial sum, and the possible last point is in a finite fiber of [h]. The abscissa image condition bounds the finite solutions by d. The specialization cannot vanish for every index. This argument uses the previous exact projected polynomial, including its boundary cases. The masked method also works for general polynomials with genuinely whole-slot fibers, so this assumption is not needed for the finite algebraic algorithm. The method is not a new faster elliptic scalar last-slot completion algorithm: if all the other abscissas are explicitly known, direct enumeration of their signs and group addition is already cheap. The advantage here is a uniform, branch-free ring-valued construction without enumerating the other-grid lanes or constructing a full determinant.

## 4. Recover all low-multiplicity witnesses using ordinary traces, without a norm or nilpotent perturbation

**Theorem 3.** With the input of Theorem 2, fix a bound K with `K+2<p`. If the total number H of distinct zero addresses is at most K, then a circuit of size

    soft-O(poly(d,b,m,K) * (1+sqrt(w)))

in B-operations plus O(m(K+1)) applications of the ordinary B-trace yields all zero addresses, with failure probability at most `H(H-1)/(2p)` from a random linear projection. The trace means `Tr_B(v)=∑_z v(z)` over the remaining grid. If W<p, the trace also detects H>K exactly. Otherwise the H≤K promise is required for the count readout; a modular count by itself does not distinguish H from H+p. Reconstructed addresses are checked directly.

Proof. Let `G_z(U)=gcd(F_z(U),Q(U))` on each nonzero-polynomial lane z, monic of variable degree `n_z≤d`. Use the degree idempotents of its coefficient array to reverse it into `C_z(t)=t^{n_z}G_z(1/t)` with constant 1. On zero-polynomial lanes set C=1 and handle separately. The logarithmic derivative yields local power sums

    P_0(z)=n_z,
    ∑_{j≥1}P_j(z)t^j = -t C_z'(t)/C_z(t).

Truncated inverse series with constant 1, or the coefficient recurrence, computes all P_j through K+1 using poly(d,K) B-operations. Squarefreeness of Q makes these power sums count distinct roots, including repeated roots of F only once. On a zero-polynomial lane add the known sums `∑_{u=0}^{w-1}u^j`, and count w at j=0. They can be constructed in poly(K,log p) operations: telescope `(u+1)^(j+1)-u^(j+1)` to solve for the j-th power sum using the previous lower ones, since j+1<p. All these procedures use only the base coefficient field and its product-ring operations.

Choose uniform independent c_i∈F_p and let `ell(U)=∑ c_i U_i`, with the affine index coordinate used in the eliminated slot. Write `ell=c_i U+ell_-(z)`. By the binomial expansion construct in B the local sums of ell^j for 0≤j≤K and the local sums of each original coordinate times ell^j for 0≤j<K. For the eliminated coordinate use P_{k+1}; for the others multiply by the corresponding coordinate in z. Take ordinary traces of these arrays, obtaining

    M_j = ∑_{s:F(s)=0} ell(s)^j,
    M_{l,j} = ∑_{s:F(s)=0} s_l ell(s)^j.

M_0 counts H modulo p. On the stated promise its canonical value is H. Newton identities, dividing only by integers at most H<p, construct the monic polynomial with roots ell(s). Factor it over F_p and detect projection collisions by a repeated-root test. For distinct projected values y_1,...,y_H, each tagged vector solves a Vandermonde system

    M_{l,j}=∑_{a=1}^H y_a^j s_{a,l},  0≤j<H.

Solving these systems recovers every original address. A pair of distinct addresses has collision probability exactly 1/p over the chosen c_i. Union bound proves the failure estimate. All extraction work takes expected poly(b,m,K) bit time (standard finite-field polynomial factorization), separately from the missing trace evaluation. No large characteristic polynomial, perturbed norm, or determinant evaluation is required in this transfer. ∎

For a singleton one only needs the count and m first-coordinate moments, without random projection: the trace of P_0 is 1 and the coordinate traces give the witness. Thus the degree-one image ideal in the preceding checkpoint is recovered immediately as the annihilator of `(1,s_1,...,s_m)`.

### Relation to the intended exponent and exact limitation

Using the preceding separate-degree bound d=h²2^(m-1), balanced `log w=b/m+O(log b)` and K with `log K=O(m+log b)` as in the previous shared-translation coverage/overflow argument, the **circuit construction** for the B-moment arrays has log size

    b/(2m)+O(m+log b)

apart from producing the original F-coefficients. A fully explicit coefficient tensor of F can be obtained by interpolation with `(d+1)^m` evaluations of the addition polynomial; its log cost is `O(m²+m log h+log b)` using the earlier resultant representation. This gives a total construction bound of `2^O(b/m+m²+m log h+log b)` with circuit nodes representing field operations on the remaining coordinates. It is not an estimate of the cost of tracing that circuit. In the previous coverage regime W<p may be arranged by thinning widths at the stated lower-order retry expense; otherwise use the promised H bound and explicitly check the returned roots. Recovered algebraic addresses still require sign enumeration and verified projected group addition. No new outer-reduction claim is made.

Directly evaluating all W/w lanes costs, even ignoring coefficient setup, `soft-O((W/w) poly(d,b,K,m) sqrt(w))`, of log scale `b-b/(2m)+O(m+log b)` at balanced widths. This implementation is worse than the old balanced half-list approach at growing arity. A trace evaluator for these succinct circuits, not a faster implementation of individual gcds, is the missing operation. Moreover if there is exactly one grid zero and no whole-slot fiber, the fiber count function is exactly the Lagrange indicator of its remaining coordinates. In the standard remaining-grid monomial basis it has degree `w_j-1` in EVERY remaining variable and maximum total degree. Thus explicitly expanding this intermediate does not preserve a small-degree invariant for repeated use of Theorem 1. Succinct circuits and alternative trace algorithms are not excluded by this fact.

## 5. Image-certificate lower bound for the actual elliptic polynomial

The preceding checkpoint established maximal degree of the kernel generator but left open whether the *image* ideal could admit short-degree certificates. The following result directly addresses that question.

**Theorem 4.** Suppose F is a polynomial of exact separate degree d_i≥2 in each variable, with standard grid widths `w_j>1+max_i d_i`. If a nonzero affine-linear polynomial L is in the image ideal F A, every standard-basis multiplier B_L satisfying

    F B_L ≡ L mod (Q_1,...,Q_m)

obeys

    deg_{U_i} B_L ≥ w_i-d_i     for every i.

Proof. The leading coefficient of F as a polynomial in U_i is a nonzero polynomial in the other variables with separate degrees at most d_j. If L has a nonzero coefficient of U_i, choose a remaining-grid lane where this leading coefficient is nonzero; such a lane exists by Cartesian interpolation. Otherwise multiply that leading coefficient by the nonzero L, which is independent of U_i. Its separate degrees are at most d_j+1<w_j. Hence again there is a lane where both are nonzero. On the chosen lane F has degree d_i and L is nonzero of degree at most 1. If the separate degree of the multiplier were less than w_i-d_i, the specialized product minus L would have degree less than w_i and vanish at all w_i distinct slot nodes. It would be the zero polynomial. But a nonzero product with the degree-d_i polynomial cannot have degree at most 1, and the zero product cannot equal the nonzero L. Contradiction. ∎

The preceding elliptic note supplies exact separate degree d_i=h²2^(m-1) for its polynomial, including under affine substitutions. At arities with `w_min≫d` (e.g. m=b^(1/3) and polynomial h), Theorem 4 holds for that actual polynomial without a generic-coefficient assumption. A rectangular monomial cutoff must therefore include degrees through w_i-d_i along every axis; a dense total-degree cutoff must include degrees through at least w_min-d. For equal widths w the latter contains at least

    binom(w-d+m,m) ≥ (w-d+1)^m/m!

standard monomials, of log size `b-O(m log m + m*d/w + log b)` when log W=b+O(log b). More generally with log w_min=(1+o(1))b/m, m log m=o(b), w_min/d→∞, this is `2^(b-o(b))`. Thus dense Macaulay-style searches for a degree-one image certificate cannot be made affordable merely by increasing the multiplier cutoff until it succeeds. This is a representation lower bound, not a lower bound on sparse multipliers, alternative bases, or implicit ideal algorithms.

## 6. Sharp degree-profile theorem via reciprocal grid moments

Let `D=∑(w_j-1)` and `ω_x=1/∏Q_j'(x_j)` for grid point x. For standard representatives, the linear functional

    Λ(B)=∑_x ω_x B(x)

is the coefficient of `∏ U_j^(w_j-1)`, by tensor-product Lagrange interpolation. The pairing Λ(AB) is nondegenerate by CRT (all weights are nonzero). If `A_{≤k}` denotes the total-degree-truncated standard monomial space, its orthogonal complement is `A_{≤D-k-1}`. Inclusion follows because reducing by the Q_j cannot increase total degree, and equality follows by dimension: exponent complementation in the standard box exchanges degrees j and D-j. This also proves the statement for k=1.

Suppose a singleton s is the only zero of F. For nonzero affine L with L(s)=0, write `L(x)=c·(x-s)`. Every multiplier is uniquely fixed off s as `L(x)/F(x)` and is arbitrary at s. Choose its value at s to kill Λ(B_L); this is always possible and yields degree at most D-1. Every other choice differs by a nonzero multiple of the Lagrange indicator at s and has degree D.

**Theorem 5.** Define the m×m symmetric reciprocal moment matrix

    R_{ij}=∑_{x≠s} ω_x (x_i-s_i)(x_j-s_j)/F(x).

Then the normalized multiplier for L has degree at most D-2 exactly when `Rc=0`. In particular, if R is nonsingular, *every* nonzero affine-linear image certificate has minimum multiplier degree exactly D-1.

Proof. A degree-D-2 multiplier is orthogonal to span(1,U_1,...,U_m). Its Λ(B_L)=0 condition is already imposed. The remaining conditions can be written `Λ((U_j-s_j)B_L)=0`. The term at s vanishes and all others are the stated reciprocal sums. This is exactly Rc=0. ∎

More generally the analogous threshold D-t-1 is characterized by the rectangular matrix with rows indexed by nonconstant standard centered monomials `(U-s)^α` of degree at most t and columns i, entries `∑_{x≠s}ω_x (x-s)^α (x_i-s_i)/F(x)`. Centering leaves the low-degree monomial span unchanged; duality gives the result. This identifies the additional sums needed for lower degree, rather than assuming that a codimension-one image ideal has a cheap low-degree certificate.

**Theorem 6 (representation-dependent sharpness after unit preconditioning).** Over a field of size q containing F_p, choose a uniformly random polynomial G of total degree at most two. With probability at least `1-(2W-1-m)/q`, G is nonzero on the entire grid and the reciprocal moment matrix for the equivalent generator FG is nonsingular. Thus, over F_{p²} with W=O(p), this statement holds with probability `1-O(1/p)` for *any* singleton instance, including the actual elliptic polynomial. Multiplication by G does not change the image ideal when G is a unit, but it can change certificate degrees. This is not a lower bound for the unpreconditioned F or for all unit preconditioners.

Proof. By Cauchy–Binet the determinant of R(FG) is the sum over m-element subsets J of grid points other than s of

    det([(x-s)]_{x∈J})² ∏_{x∈J} ω_x/[F(x)G(x)].

Multiplying by `∏_{x≠s}G(x)` gives a polynomial in the coefficients of G of degree W-1-m. Show this polynomial is not identically zero. For each i choose t_i∈S_i distinct from s_i. The simultaneous zero set of the quadratics `(U_i-s_i)(U_i-t_i)` and `(U_i-s_i)(U_j-s_j)` for i<j is exactly s and its m single-coordinate neighbors n_i (n_i differs from s only at i, taking t_i). Take a formal generic linear combination G_0 of these quadratics; with coefficients in a rational function field it vanishes at exactly these points on the finite grid. Substitute `G=G_0+z` with z another indeterminate. The coefficient of z^-1 in R(FG) is the diagonal matrix whose i-th entry is `ω_{n_i}(t_i-s_i)²/F(n_i)`, all nonzero. Its determinant has a nonzero z^-m leading Laurent coefficient. This proves the polynomial is nonzero over the base field, even in finite characteristic. Schwartz–Zippel bounds its vanishing probability by (W-1-m)/q. Each G(x) is a nonzero linear functional in the random coefficients (the constant coefficient is free); union bound charges W/q for failure of unitness. ∎

The sharpness result is intentionally scoped. It does not justify genericity of the actual addition polynomial's reciprocal matrix under shared translations. It establishes that even an equivalent generator can force maximal-degree certificates, while the direct separate-degree obstruction of Theorem 4 applies to the actual unpreconditioned polynomial.

## 7. Exact component regression

Run `python -B verify.py` (Python standard library only). The block algorithm uses a true monic subproduct/remainder tree with schoolbook arithmetic, not the fast algorithms used for the theoretical operation count. Every B lane is explicit in the verifier, not an implementation of the open trace evaluator. The elliptic polynomial is computed by independent declared-degree Sylvester determinant and interpolation helpers adapted from the preceding verifier. Actual elliptic regressions are for h=1; no cofactor-specific experiment is claimed in this checkpoint. Polynomial factoring on the tiny fields is performed by scanning roots only in the tests. Whole-slot moments use the telescoping power-sum recurrence; multiple-root and degree-drop cases are included. The random quadratic unit trial reports nonunits and the singular unit case, rather than retrying until success. The exact diagonal pole certificates check the nonvanishing proof independently of random trials.

Recorded output in `verification.txt`:

- 200 block remainder comparisons over product rings, including 117 tail cases.
- 821 masked gcd/moment lanes; 188 mixed-degree fixtures, 144 whole-slot lanes, 56 repeated-root input lanes.
- 135 actual Semaev grids, 1,097 slices including 10 degree drops, 3,752 literal address comparisons. Across two projection choices per grid: 246 empty responses, 21 successful full recoveries, 3 correctly detected projection collisions, 26 recovered addresses counting repeats across tests.
- 754 exact degree/reciprocal-matrix equivalence comparisons; 2 singular unpreconditioned matrices among 100 fixtures.
- 100 exact diagonal pole certificates; 87 nonsingular random unit trials, 12 nonunit trials and one singular unit trial (none retried).
- 50 actual Semaev affine-linear image generators on 20 singleton grids satisfy the separate-degree lower bound (217 sampled grids tested to find those singleton instances).

All assertions passed. These are component identities and finite comparisons. They do not benchmark a faster forward relation search.

## 8. Frontier

The new constructive object is a branch-free circuit for fiber gcd moments followed by ordinary trace, with one coordinate compressed from w to at most d and no determinant or square-zero coefficient ring required for witness extraction. The direct image-certificate theorem rules out affordable *dense* multiplier cutoff searches, even before the sharp generic total-degree obstruction is applied. The remaining trace is still over W/w lanes. Expanding the projected indicator restores full remaining-grid degrees. The research gate is a collective trace algorithm exploiting the particular gcd circuit, or a different simultaneous elimination that does not create high-degree remaining-grid functions. No strict asymptotic improvement, new numerical coefficient, or unconditional ECDLP algorithm is claimed.
