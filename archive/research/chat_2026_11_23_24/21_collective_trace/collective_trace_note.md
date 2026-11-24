# Projected norms, batched rational traces, and a resolvent readout obstruction

Research checkpoint, 24 November 2026. Continues `fiber_gcd_note.md` §§2–4 and `shared_norm_aggregate_note.md` §§2–5. Objective: beat the conditional general-prime log-time scale `b/(log b)^(1/4)`. **This checkpoint does not establish that improvement.** It simplifies the missing collective evaluator, removes unnecessary per-fiber gcd operations, gives a rational alternative and quantifies two limitations. The coefficient 11.31486814 from earlier checkpoints is still an unrevalidated conditional candidate. No priority or independent peer-review claim is made.

## 1. Finite-grid model and hypotheses

Let p be an odd prime, b=ceil(log2 p), S=∏_{i=1}^m S_i a finite Cartesian grid of distinct field elements within each slot, w_i=|S_i|, and W=∏w_i. Write `A=⊗_i F_p[U_i]/(Q_i(U_i))`, where `Q_i(U)=∏_{u∈S_i}(U-u)`. This algebra is a product of W copies of F_p. The ordinary trace is the sum of the values at the lanes. Let F be a polynomial on this grid, and write H for its number of distinct zero addresses.

Eliminate one interval slot, normalized by its public affine substitution to U=0,...,w-1, leaving z in the remaining grid and its algebra B. There are N=W/w lanes of B. The degree of F in U is at most d (degrees may differ by lane). A B-operation is a node in a succinct arithmetic circuit, NOT a bit operation on a field element. Explicit B arithmetic still requires processing N lanes.

For the elliptic applications we use the preceding checkpoints' projected summation polynomial and their common-translation sampler, not a new proof of every upstream claim. The actual elliptic finite regressions below have h=1. Cofactor statements in §5 depend on the earlier bidegree and branch-locus results. Algebraic addresses may have nonrational lifts, so one must still perform the rational-lift and projected-group-sum verification, under the previously stated sampling/coverage assumptions. We do not claim a new outer-reduction audit.

## 2. Projected support before factoring: no masked Euclidean algorithm

Define

    D(z) = ∏_{u=0}^{w-1} F(u,z),
    J(z) = 1 - D(z)^(p-1) ∈ B.

By the CRT and Fermat, J(z)=1 iff at least one address above z is a zero. This holds even if F_z=0 as a polynomial, if F_z has degree drops or repeated roots, or if several witnesses share the same remaining address. J is NOT the number of zeros above z; it counts positive fibers once.

**Theorem 1 (simple projected-support recovery).** Suppose the number C of positive fibers is promised at most K<p. Constructing D and J from the coefficient circuit of F takes soft-O(d+sqrt(dw)+b) B-operations for d≥1; for d=0 the constant polynomial product is handled by ordinary binary powering. If one can compute the traces of J and of its low-degree address-tagged powers as described below, every projected positive fiber, and then every actual zero address under a bound H≤K, is recoverable in expected time polynomial in b,m,d,K, with projection failure probability at most C(C-1)/(2p). An identically zero polynomial lane contributes all w original addresses, and these can be reported whenever the output-size promise permits it.

Proof of the construction cost. This is the earlier block norm algorithm rather than a new fast polynomial-arithmetic theorem. Put `B_H(X)=∏_{j=0}^{H-1}F(X+j,z)`. Taylor shifts, a balanced polynomial product tree, and monic multipoint evaluation at qH for 0≤q<floor(w/H) give the full block products. Include the final fewer than H factors. Cost is soft-O(dH+w/H+d) operations, using fast arithmetic over commutative algebras; no denominator in the coefficient ring must be inverted by the product and remainder trees. Choose H of order sqrt(w/d), clipped to [1,w]; if d>w use evaluation on the whole interval. Taylor shifts by factorial convolution use d<p, as in the intended application; ordinary shift algorithms also suffice in bounded-degree regimes. A final binary power costs O(b) B multiplications. The exact finite verifier uses schoolbook arithmetic and a genuine monic multipoint tree. It verifies correctness, not fast-arithmetic timings.

Recovery. Let `ell(z)=∑_{j≠i} c_j z_j` for uniform independent c_j. Obtain

    M_k = τ(J ell^k),         0≤k≤C,
    M_{j,k} = τ(z_j J ell^k), 0≤k<C.

Here τ denotes ordinary trace in B. M_0 gives C, by the promise C<p. It suffices to first obtain M_0, or to compute moments up to K. Newton identities with denominators 1,...,C construct the polynomial whose roots are ell(z) for the C positive fibers. Each distinct pair collides with probability 1/p. Factor, reject projections with repeated roots, and invert the C×C Vandermonde matrix to recover every z_j. Evaluate F_z in the eliminated variable ONLY for the recovered positive fibers; use finite-field factorization, discard non-base-field roots and roots outside the interval, and handle a zero polynomial separately. This produces the actual roots, even when different original addresses have the same projected address. Under the output-size promise the polynomial setup, factoring, and recovery work is expected poly(b,m,d,K), in addition to any cost of evaluating the original F coefficient circuit. The cost of the collective traces is not included in that recovery overhead and is still open. On a singleton, no random projection is needed: τ(J)=1 and τ(z_j J) is the coordinate. ∎

The key simplification relative to the previous masked Euclidean algorithm is that no per-lane leading coefficient masks, degree tests, gcd coefficients, or modular division are needed. The final circuit uses one high-power indicator D^(p-1). All tag baselines τ(ell^k) and τ(z_j ell^k) can themselves be calculated from low-order slot power sums via elementary polynomial arithmetic and multinomial identities (or enumerating each individual slot), without traversing the Cartesian product. Thus the new sufficient evaluator is a weighted trace of one high power, not general trace access to the entire previous gcd circuit. This simplification is NOT a proof that this power trace is efficient. The original unprojected zero count is already a power trace of F; our result reduces one coordinate before the power test, and trades a more expensive block norm circuit for that lower lane count.

If one explicitly enumerates N remaining fibers, the log of the cost, even with the fast block construction, is about `b-b/(2m)+O(m²+m log h+log b)` for the previous balanced widths. This is worse than the older half-list method asymptotically. At m≈b^(1/3), setup alone may fit 2^O(b^(2/3)), but the trace is not known to fit. No gain to the full exponent is claimed.

## 3. Audit of the initially proposed projective normalization

A projective degree normalization does work algebraically, but its determinant is just a predictable scalar times the simpler block norm. It is not a second, faster evaluator. We record the identity to avoid paying for an unnecessarily complicated circuit.

For any scalar fiber f(U) of degree ≤d, d≥1, choose α in an extension K of F_p, and set `a=f(α)` and

    t(V)=V^d f(α+1/V)=∑_{j=0}^d t_j V^j,  t_d=a,
    g(Z)=Z^d+∑_{j<d}t_j a^(d-1-j) Z^j.

This is a monic polynomial with no inverses, even when a vanishes. If a≠0, it equals `a^(d-1)t(Z/a)`. Nonzero roots correspond to finite roots of f via `Z=a/(U-α)` (when α is not a grid point). If deg f<d, the extra roots at Z=0 correspond to the infinite point of the chart. Let

    R(Z)=∏_{u∈S_i} (a+(α-u)Z).

Let M_R be multiplication by R modulo g in the d-dimensional monic quotient. Then

    det M_R = a^((d-1)w) ∏_{u∈S_i} f(u).                 (1)

Proof. Determinants multiply. For each factor with α≠u, the determinant of multiplication by a+(α-u)Z is `(-1)^d(α-u)^d g(-a/(α-u))`, equal to a^(d-1)f(u) after substitution. This identity extends to a=0 and α=u as a polynomial identity, or directly by the determinant of the linear factor. Multiplication over u proves (1). The case d=1 uses a^0=1, including a=0. All arithmetic can be realized in B by the same polynomial identity in the coefficients. ∎

If α has extension degree d+1, then no nonzero base-field polynomial of degree ≤d can vanish at α, so this normalization can be made uniform across every nonzero fiber; an irreducible extension polynomial can be found in expected polynomial time in d,b by standard irreducibility testing and random sampling. But (1) shows why this additional extension and determinant are superfluous for the projected-support test. We use the simpler product from §2. The tests explicitly cover a=0 and degree drops; the projective recovery interpretation assumes the valid chart, while the determinant identity does not.

## 4. Rational trace alternative: one dual-number block product for all tags

Here is a different exact representation, avoiding Fermat zero-test powers at the circuit-construction level. It replaces them with one unit inverse in the remaining algebra and a final global residue problem.

For arbitrary polynomials P_j(U,z), 1≤j≤k (tags can include 1, coordinates, or powers of a random linear form), form over B with a formal parameter t the block norm

    N(t,η) = ∏_{u∈S_i} (F(u,z)+t+∑_{j=1}^k η_j P_j(u,z)),
    η_iη_j=0 for every i,j.

Write `N=N_0+∑η_j N_j`. Then, in the localization where N_0 is a unit,

    N_j/N_0 = ∑_{u∈S_i} P_j(u,z)/(t+F(u,z)).             (2)

Proof: expand the product and retain the single η_j term. Divide by the common product. ∎

Every denominator is nonzero in the generic field F_p(t); at any fixed t in F_{p²}\F_p it is nonzero on every lane, even zero-polynomial fibers. No leading coefficient or squarefreeness assumption is needed. Compute the block product with all tags in the square-zero algebra. If D=max(1,deg_U F,max deg_U P_j), the operation count is soft-O((k+1)(D+sqrt(Dw))) in B (or in its scalar extension), since square-zero multiplication has linear overhead in k. Use one common unit inverse of N_0 in B followed by k multiplications. Over a finite extension this unit inverse may itself be implemented by powering with O(b times extension degree) B multiplications; it is not a free collective inverse. Formal rational operations do not imply small numerator/denominator descriptions after expansion.

Taking traces gives the rational function

    R_j(t) = τ(N_j/N_0) = ∑_{x∈S} P_j(x)/(t+F(x)).

Its residue at t=0 is EXACTLY `∑_{x:F(x)=0} P_j(x)` over F_p. The poles are simple even if many addresses share an evaluation. Repeated roots of f_z cause no difficulty because the sum counts distinct grid addresses rather than algebraic multiplicities. Using all desired low-degree tags recovers the same witness moments as the previous gcd method. On a singleton, the residue with P=1 is 1 and that with P=U_i is its address coordinate. More generally apply the Newton/Vandermonde recovery from the previous checkpoint on the full address grid. This alternative does not count projected fibers: its residues count all zero addresses. It can avoid final Fermat tests *if* the residue of the collective rational function can be obtained directly. We do not know such a fast aggregate algorithm.

An independent single-fiber reciprocal identity may also be useful. For monic f of degree d≥1, Q monic grid polynomial, gcd(f,Q)=1, and deg P<d,

    ∑_{u:Q(u)=0} P(u)/f(u)
       = -[U^(d-1)] (P Q' Q^(-1) mod f).                 (3)

Here Q is squarefree, inverses are in the small quotient, and the equality holds without requiring f to be squarefree. To see this, the rational function P Q'/(f Q) has vanishing residue at infinity. Its residues at the grid roots equal P(u)/f(u), and the sum at the roots of f is the coefficient of U^-1 of its proper part with denominator f, equivalently the indicated top coefficient modulo f. The proof works over the algebraic closure and by CRT over product algebras with the required invertibility. Q' mod f can be obtained by differentiating Q mod f²; the earlier block remainder theorem therefore applies with degree 2d when its hypotheses hold. Polynomial parts of a numerator of degree ≥d can be treated separately by elementary power sums. In this checkpoint the implemented rational evaluator is (2), which is simpler, also handles degree drops, and does not require monic charts. Equation (3) is an alternative mathematical identity, not separately benchmarked.

## 5. Large projected spectra: an obstruction on the actual shared-translation inputs

The preceding shared-norm note established the large spectrum of the original multiplication operator. We extend that argument to *any partial block elimination*. This is the key limitation for scalar filters applied after the new projected norm.

Assume the shared-translation model and prior branch-locus hypotheses: common random A∈F_p, disjoint public S_i with pairwise coprime shifted cubics; subgroup order r odd prime, cofactor h=poly(b). For each full offset tuple s, use the bivariate specialized relation polynomial F_s(A,X) of the preceding note. We use its asserted irreducibility, nonassociation across different s, and bidegrees

    a_0 = h² m 2^(m-1),   x_0 = 2^(m-1).

These follow there from the independent-sign cover, the degree of the sum morphism, and the finite branch-locus fingerprint. Multiplying any individual F_s by a nonzero field scalar does not affect the argument. In particular use the normalization of the universal relation polynomial used in the algorithm, conditional on the earlier result that each such specialization is a scalar multiple of the primitive polynomial. No independent cofactor re-audit is claimed here.

Take any nonempty proper subset I of slots to eliminate. Write `L=∏_{i∈I}w_i`, `N=∏_{i∉I}w_i`, W=LN, and for a remaining offset tuple z define

    D_z(A,X)=∏_{u∈∏_{i∈I} S_i} F_(u,z)(A,X).

The products for distinct z are pairwise nonassociate. Their irreducible factors correspond to disjoint sets of full address tuples, so they cannot be equal up to a scalar. They have bidegrees L a_0, L x_0. Fix T uniformly from the nonidentity subgroup and set X=x(T). Since r is odd, X is uniform on (r-1)/2 distinct abscissas.

**Theorem 2 (projected spectral energy).** Define E(A,T)=∑_v n_v², where n_v counts how many of the N projected norm values D_z(A,x(T)) equal v. Let

    e = min(1, L a_0/p + 2 L x_0/(r-1)).

Then E[E]≤N+N(N-1)e. On any event U of probability κ>0, with conditional probability ≥1-ε, the number s of distinct projected norm values satisfies

    s ≥ κ ε N / (1+(N-1)e).                            (4)

Proof. For z≠z', the difference D_z-D_z' is a nonzero polynomial of the stated bidegrees. The exceptional translations at which the whole polynomial in X vanishes are bounded by its degree in A (use one nonzero coefficient); for other translations its number of X roots is bounded by its degree in X. The resulting collision probability is ≤e. Sum over ordered pairs. Conditional expectation is at most E[E]/κ, and Markov's inequality combined with `N²≤s E` proves (4). ∎

In the prior balanced regime W=O(r), the denominator in (4) is at most 2^O(m)poly(b) since `N L = W` and r≤#E≈p. If one conditions on an actual successful rational relation event with κ≥2^-O(m)/poly(b), using the earlier common-translation coverage theorem, ε=1/b yields

    s ≥ N / (2^O(m) poly(b))

with conditional probability ≥1-1/b. For one eliminated balanced slot, log2 s ≥ b-b/m-O(m+log b), with the earlier width precision log2 w=b/m+O(log b); with weaker width precision write log2 N-O(m+log b) instead. The proof remains valid for deterministic slopes and for block choices fixed before drawing the random seed. It is not asserted for a seed-adaptive selection of I unless the union-bound cost is charged.

**Corollary.** On a successful instance with s distinct projected norm values, ANY rational scalar filter A(X)/B(X) that is defined at all these values and equals one at zero and zero at every nonzero norm value must have `deg A ≥ s-1`. Its coefficients may be chosen after seeing the instance. This follows because A(0)≠0 and A has s-1 distinct nonzero roots. Thus the new single zero-test power cannot simply be replaced by a low-degree polynomial or globally defined rational function of the projected norm alone. The bound concerns degrees of pointwise filters, not arithmetic circuit size (binary powering already gives a high degree in O(b) gates), and not a lower bound for evaluating the sum of those filters by another representation.

It also follows that the rational Cauchy transform of the projected spectrum, `∑_v n_v/(t+v)`, has almost s distinct poles when represented as a reduced rational function: any cancellations correspond to n_v=0 mod p, and at most floor(N/p) distinct values can have such positive integer multiplicities. In particular if N<p there is no such cancellation. For fixed actual input, explicitly reconstructing this reduced rational denominator still costs exponential degree. This does not exclude implicit extraction of its residue.

## 6. Extension-valued resolvents and jets: a precise generic readout lower bound

The actual spectral theorem motivates using rational evaluations rather than whole-spectrum reconstruction. There is a rigorous limitation even on adaptive resolvent-only algorithms if they must work for arbitrary field-valued histograms on a known set of distinct poles. The model in this section is intentionally stronger in its input family than the actual elliptic histograms.

Let S⊂F_p be s distinct values containing 0, and `n=(n_v)_{v∈S}∈F_p^S` arbitrary. The goal is n_0. Give the decoder the total ∑n_v for free. It may make adaptive queries at extension-field points t_j in F_{p^{e_j}}, excluding all -v, and obtain the entire first k_j+1 Hasse-Taylor coefficients of `R_n(t)=∑_v n_v/(t+v)` at t_j. The coefficient of order k is `(-1)^k∑ n_v/(t_j+v)^(k+1)`. Each such query returns e_j(k_j+1) base-field words.

**Theorem 3 (exact readout).** Any deterministic algorithm correct for every n must use, on its all-zero histogram transcript, total query budget at least s-1 base-field words, counting extension coordinates and jets as above. This holds for arbitrary nonlinear postprocessing and adaptive query choices. If the input total is prescribed, the same statement holds on any fixed reference histogram transcript within the affine space of that total. It does not assert a randomized bounded-error lower bound or an arbitrary computational lower bound on actual elliptic input histograms.

Proof. Freeze the reference histogram (zero if the total is zero), and follow its transcript. The chosen queries are fixed along this path and define F_p-linear measurement rows. If n_0 were in their row span together with the total row, its selector v↦1_{v=0} would equal a function on S of the form

    c + ∑_{j,k≤k_j} Tr_{F_{p^{e_j}}/F_p}(a_jk/(t_j+v)^(k+1)).

This uses the nondegenerate trace pairing to express any linear functional on each returned extension value. Expanding each field trace gives a rational function in v with poles among the conjugates of -t_j and order at most k_j+1 at those conjugates. None is in S, because every v∈S lies in the base field. A common denominator has degree at most E=∑ e_j(k_j+1); after multiplication the numerator has degree at most E. To be the zero selector off v=0 it needs at least s-1 distinct roots while remaining nonzero at 0. Therefore E≥s-1.

If E<s-1, the row-span membership fails. Linear duality then supplies a vector Δn with total zero and all measurement rows zero but Δn_0≠0. Add it to the reference histogram. Every successive returned value agrees with the reference path, so adaptation cannot distinguish the two histograms, and a common output cannot be correct. ∎

Matching control: if p-s≥s-1, the total and s-1 distinct admissible base-field resolvents have full rank. A linear combination vanishing on all S has numerator of degree ≤s-1 with s roots, hence is identically zero. Distinct simple poles then force each coefficient to be zero. Thus the order lower bound is exact for this generic measurement model. Counting an extension query as one operation without its e coordinates would hide the cost.

Theorem 3 is for arbitrary F_p-valued histograms, not necessarily nonnegative integer multiplicities of a fixed-size elliptic grid. In particular we do not infer that an algorithm exploiting correlations among the actual histograms needs s-1 measurements. The actual-input conclusions established independently in §5 are the large spectrum, large degree of pointwise scalar filters, and large reduced denominator. A circuit-aware implicit evaluator or a different nonlinear statistical model remains open.

## 7. Verification and reproducibility

Run `python -B verify.py` in this directory. Standard library only. The file `previous_fiber_verifier.py` is an unmodified copy of the previous checkpoint's verifier, included because the new tests reuse its small declared-degree Semaev resultant, interpolation, determinant and linear-system helpers. No previous test output is counted as a fresh regression. The new script uses exact modular arithmetic, explicit finite product-ring lanes, a quadratic extension, and square-zero tags. Polynomial arithmetic is schoolbook; no claim of a practical crossover or fast collective trace implementation is made.

Recorded results in `verification.txt`:

- 220 block products, 128 nonempty tails, 759 exact projected masks; 107 zero-polynomial fibers and 20 repeated-root inputs. 220 simultaneous dual-number norms, 2,277 independent tagged reciprocal comparisons and 660 residue comparisons.
- 300 exact projective determinant identity checks, including 231 degree drops and 38 a=0 cases; 130 valid-chart zero tests.
- 108 actual Semaev grids, 863 fibers, 3,200 address comparisons. 22 positive fibers in the random grids; 211 successful readouts (including empty ones), five detected projection collisions, 32 recovered addresses counted with repetitions across readouts. An explicit mixed whole-slot/multiple-hit control recovered seven addresses. In addition, 288 constructed genuine elliptic two-hit fibers (576 valid addresses) were verified against independent direct group addition on the order-19 curve over F_13.
- 500 exact dual collisions for the generic resolvent-readout obstruction; 1,792 base-field jet rows and 1,240 independent zero-transcript comparisons; 460 full-rank positive controls.
- 14 distinct projected-polynomial pair comparisons using 5,132 signature evaluations at sufficient bivariate interpolation grids, and 31,068 exact projected spectrum observations on two small prime-order curves, of which 2,030 were algebraically positive. Observed mean collision energies are logged; these finite cases do not establish the asymptotic spectrum theorem, which follows from the bidegree argument.

All assertions passed. No failure was suppressed by repeating until success, except where a deliberate construction is identified (e.g. finding a prime-order toy curve). No independent elliptic cofactor regression is claimed. The proof of (3) is not a separate implemented test; the implemented rational trace is the directional derivative (2).

## 8. Frontier and candid accounting

The useful new objects are a single projected-support power trace and a one-inverse, all-tag rational trace with a final residue. The projective matrix idea is algebraically correct but redundant, as the exact determinant identity demonstrates. The earlier gcd circuit can be retained for other local outputs, but it is not necessary for forward witness location.

The projected spectral theorem and generic resolvent readout theorem rule out several naive scalar-decoding approaches without asserting a general computational lower bound. The large spectral degree does not preclude short arithmetic circuits. We have not established an efficient method for collectively evaluating a tagged power trace of the block norm, or extracting the rational residue without reconstructing the spectrum. In particular, the previously proposed 2^O(b^(2/3)) construction-cost bound still excludes this load-bearing evaluation. The previous conditional general-prime coefficient has not been re-audited and the asymptotic exponent has not been beaten.

The next mathematically relevant gate would be a circuit-aware algorithm for the trace of a **single** block-norm power, or one which extracts a specific residue directly from the factorized norm circuit instead of learning its Cauchy transform from point queries. The generic resolvent route alone no longer offers a credible argument for a saving.
