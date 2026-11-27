# Enumerable prescribed fibers at the quarter-power scale

27 November 2026. Baseline `main` at `04a8b3d8bb904af8a30aa996f1af4027b52f6dd4`; compare open PRs #9–#12, in particular #11 at `e833763d99b9e17e3c3544e0e01ad409bc111ed6` and #12 at `8a2ee0b78d1ea37e74a8d2ff700ed3bb819213f0`. We change the seed and measurement architecture. The principal new result is a deterministic construction of explicitly enumerable, real-nondegenerate prescribed fibers at the original quarter-power scale. A second result improves the cost of the earlier *globally* pointwise-uniform seed via batched polynomial evaluation. A third result supplies a positive response to #12. None establishes a uniform improved end-to-end ECDLP theorem. Resource claims assume LZR ePrint 2026/3020, Theorem 2.2, the same closed exact-counting interface and canonical arithmetic gadgets as in the repository's `MAIN_RESULT.md` and `COMPLEXITY_AUDIT.md`.

## 1. A pointwise second moment around the empirical mass

Let G be a cyclic group of odd prime order r, and let F_i be independent random subsets of G minus its identity. Within a block the membership indicators are pairwise independent with marginal theta. Put q=(r-1)theta. This law is available on a whole prime-order elliptic curve over F_p, p>3, using independent parity-specific affine masks ax+c, with a,c uniform over F_p **including zero slopes**. Retain a nonidentity point exactly when its canonical evaluation is in [0,T), theta=T/p. Distinct points of the same y parity have different x coordinates. Evaluation at two such x values is a bijective map of the two affine coefficients, and the parities are independent. Thus the stated pairwise independence is exact, not the inequality for a nonzero-slope affine law. A size guard below handles rare constant masks.

For any k>=2 retained blocks let d(g) be their tuple representation count, W their actual product of sizes, and mu_k=q^k/r. Then, for every fixed g,

    E[(d(g)-W/r)^2] <= 2 mu_k.                         (1)

Proof: use the mathematical complex characters of G, with unnormalized Fourier transforms. There is no character evaluation in the algorithm. For nontrivial characters chi,psi and one block indicator f, pairwise independence and the omission of identity give

    E[hat f(chi) conj(hat f(psi))]
      = theta(1-theta)(r 1{chi=psi}-1)+theta^2.

Indeed the mean of either nontrivial transform is -theta. The centered covariance is the sum over the r-1 nonidentity points, since the off-diagonal membership covariances vanish. The displayed quantity is at most q in absolute value on the diagonal and at most theta off the diagonal. The zero-frequency term of the product of transforms is exactly W. Fourier inversion, squaring, independence across blocks, and the triangle inequality therefore give

    E[(d(g)-W/r)^2] <= (r-1)q^k/r^2 + (r-1)(r-2)theta^k/r^2
                      <= q^k/r + theta^k <= 2q^k/r.

The last step uses (r-1)^k>=r for k>=2. In contrast with a moment around the nominal expectation q^k/r, this removes the common random block-size fluctuations. This distinction is important for the selected-pin potential.

## 2. Prescribed fibers: deterministic seed with only affine masks

Let b=ceil(log2 p), L=log2 b, m>=4 growing as Theta(L^(1/4)). Let S be an arbitrary nonempty explicit list of ell group values (known as points, not necessarily with known scalars), ell<=2^{C b/m} for any fixed C. Include O,P,Q if desired. Duplicate values can be removed. Put

    K=2^20 m^2(ell+1), q*=ceil((Kr)^(1/(m-2))),
    T=ceil(p q*/(r-1)), q=(r-1)T/p, mu=q^m/r.

For sufficiently large b, T<=p, q>=q*, q^(m-2)>=Kr, and mu=(1+o(1))K q^2. Hasse gives log2 r=b+O(1). In this range q,mu,ell all have logarithm O(b/m). The bounded initial range of parameters is not an asymptotic claim; direct finite counting is possible.

For a set I of pinned blocks of size zero, one, or two, and a tuple of selected pins A in product_{i in I} F_i, let

    e_I(g,A)=d_{-I}(g-sum A)-W_{-I}/r,
    E_I(g)=sum_{A in product F_i} e_I(g,A)^2.

For I empty this is just the full fiber squared deviation. The pinned labels and their selection indicators are independent of the remaining blocks. Applying (1) at each possible group value and taking the expected number q^|I| of selected pin tuples gives the uniform inequality E[E_I(g)]<=2mu, including I empty.

Use the nonnegative rational potential

    Phi = m^4 sum_i ((|F_i|-q)/q)^2
          + sum_{g in S} [ E_empty(g)/mu^2
                        +(q^2/mu^2) sum_{1<=|I|<=2} E_I(g) ].       (2)

The exact variance of a block size is q(1-theta), and mu>=Kq^2. Consequently

    E Phi <= m^5/q + 2 ell/mu + 2 ell m^2 q^2/mu
          <= m^5/q+4 ell m^2/K < (1/64)^2

for sufficiently large b. All constants are deliberately loose.

### Conditional expectation is affordable

For each g,I, the exact integer-scaled quantity is

    r^2 E_I(g) = sum_{X,Y: X_I=Y_I}
                 (r 1{val X=g}-1)(r 1{val Y=g}-1),                (3)

where X,Y are full retained short tuples. Expand the product into four ordinary counts with the specified shared slots, one or both total constraints, or neither. In particular no selected-pin list has to be enumerated inside the expectation. Size terms use single-point or two-point membership counts. Compute each expectation conditional on the current seed prefix by including the remaining affine coefficients as variables, restricting them to their canonical intervals and dividing by the number of seed completions. The coefficient ranges and the prefix completion numbers have polynomial bit length. All signs and rational arithmetic occur after exact nonnegative counts. Each network has O(mb) Boolean edges, O(b^(3/4)) maximum degree using the inherited arithmetic gadgets, and O(mb) exact-count precision; the two mask coefficients per parity are shared via equality trees. In particular the zero slope is not inverted. There are O(ell m^2) count types per prefix and O(mb) prefix decisions. Enumerating ell explicitly and all polynomial work on these rational quantities is charged. Choose the smaller nonempty prefix expectation; the final seed has Phi<epsilon^2 with epsilon=1/64.

Write n_i=|F_i|, W=product n_i, d=d(g). The size term in (2) gives |n_i/q-1|<=epsilon/m^2 and W/q^m in [1-epsilon/m,1+2epsilon/m]. Every selected g has mu/2<=d<=2mu by the empty-pin term. Thus the final block sizes are O(q), even though slopes zero were permitted. Exact prefix counting enumerates each block and **the entire short fiber at every supplied g** with at most O(mb) times the number of nonempty leaves plus their siblings. No high-degree external score lookup occurs. Given the LZR closed premise the seed, all enumeration, and any polynomial-time exact linear algebra on these explicit arrays fit

    T,S <= 2^{O(b/m + mb/sqrt(log b))} = 2^{O(b/L^(1/4))}.       (4)

The bound also charges the ell input points and at most O(ell mu) explicit tuples. This is a new seed selected for the prescribed group values, not an assertion that the old seed satisfies (2).

## 3. All real directions detected in each prescribed fiber

For an arbitrary real labeling u_i(A), subtract its uniform mean within each block. Let V_i=(1/n_i)sum_A u_i(A)^2 and V=sum_i V_i. On the uniform distribution of short tuples of sum g, write Z=sum_i u_i(X_i). We claim

    V/2 <= Var(Z | val X=g) <= 2V.                              (5)

This is a statement simultaneously for every labeling, not a concentration test for one residual direction. The conditional one- and two-slot masses are exactly (W_{-I}/r+e_I(g,A))/d. Write c=W/(rd); then |c-1|<=2epsilon. For one-slot diagonal contributions, n_i|e_{i}(g,A)|/d<=4(q/mu)sqrt(E_i(g))<=4epsilon. Hence the sum of the conditional second moments of the individual slots differs from V by at most 6epsilon V. The independent contribution to each off-diagonal term is zero by centering. By Cauchy--Schwarz

    |E[u_i(X_i)u_j(X_j)|g]|
       <= a_ij sqrt(V_i V_j), a_ij=sqrt(n_i n_j E_ij(g))/d.

The sum of a_ij^2 over i<j is at most 16epsilon^2 by (2). Also sum_{i<j} V_i V_j<=V^2/2. Thus the doubled sum of cross terms in E Z^2 is bounded in absolute value by 4sqrt(2)epsilon V<6epsilon V. Finally

    |E[Z|g]|^2 <= V sum_i n_i E_i(g)/d^2 <= 8epsilon^2 V/q.

This proves (5), since 1-12epsilon-8epsilon^2>1/2 and 1+12epsilon<2. In particular the proof needs an aggregate squared error for one- and two-pin counts, not their individual pointwise uniformity on all group values.

Let D=sum_i(n_i-1). Delete one reference column per block and let c_X be the reduced incidence vector. For any g in S take one reference tuple X_g. The explicitly listed contrasts c_X-c_Xg over the complete g fiber have rank exactly D over Q: a rational null vector has conditional score variance zero, so by (5) every block label is constant, and the references are zero. Every contrast has squared row norm <=2m. A rational basis and Hadamard therefore yield the modular rank lower bound

    rank_r(complete g contrasts) >= D-floor(D log_r(2m)/2).      (6)

The true nonzero log-difference vector lies in the modular kernel of the contrasts, so their modular rank is at most D-1. We do not infer modular nondegeneracy from (5). For actual rows the integer lattice rank, modular rank, and determinantal divisors can all be computed on the explicitly enumerated matrix within (4), not inferred from a compressed Gram. Stacking contrasts at any other g in S only improves rank. With O,P in S append c_XP-c_XO of known RHS one, which is independent of every zero-shift contrast. This gives at least the bound (6) plus one (capped at D) for the calibrated system.

If a single selected fiber attains D-1 contrast rank, or the stacked fibers do, that calibrated system has full modular rank. With Q in S the target is then obtained from any enumerated Q-tuple relative to X_O. More generally the full enumerated list gives an **exact** test of multi-representation identifiability: compute W from the calibrated rows and all available zero-shift contrasts, including those from the Q fiber, and test c_XQ-c_XO in W. This is equivalent to testing whether the affine hull of all the projected target difference vectors contains zero, and needs no further implicit discrepancy selection inside those fibers. A positive test certifies a target using the known right sides; a negative test is honestly reported. A remaining obstruction is whether the actual selected seeds necessarily give sufficient modular rank or a sufficient target certificate. The covariance estimate does not decide it.

A further useful property of these fibers is an explicitly constructible matching of at least n_min/(2m) disjoint-label tuples (rounding down is harmless). Each marginal is <=2/n_i by the diagonal bound. Greedily select a tuple and discard tuples sharing any label in the same block: at most 2d sum_i 1/n_i of the d original tuples can be discarded per step. Thus the matching claim follows. The matching is much larger than the relative-rank-defect scale D log m/b at the chosen m, but does not by itself give new directions in an arbitrary supplied quotient.

## 4. Positive response to PR #12 by explicit sparse fibers

On the globally pointwise-flat blocks of #12 the short count is at most 2K=2^{O(b/m)} at any supplied group value. Under an affordable ordinary count-prefix interface its entire identity fiber is therefore enumerable in a small number of calls. If the quotient of a supplied calibrated subsystem has dimension one and the scores in that fiber are nonconstant, one genuine zero-shift contrast raises rank to full. Equivalently two distinct short-fiber scores may be transported to a prescribed target by adjoining a common short tuple of the required complementary group value. For higher corank, all enumerated contrasts give the exact quotient span and a checkable sufficient condition; no arbitrary table lookup is required.

For the random-augmentation existence construction of #12, the condition at the **single prescribed identity fiber** holds with probability 1-o(1). Consider only tuples consisting entirely of newly retained points. Their expected count at this value is (1+o(1))K. To see this, each allowed ambient set is the complement of the small dictionary plus identity; expand their convolution around the constant functions as in #12 section 3. The nonconstant contribution has magnitude at most 2^{O(m)} before multiplying by theta^m, and the main product of masses gives (n-L)^m/r=(1+o(1))K. For two distinct such tuples sharing j>=1 slots and having that same prescribed sum, the two linear group equations on 2m-j distinct ambient variables have rank two, since each tuple has an unshared slot. Thus the expected number of overlapping ordered distinct pairs is at most

    (q')^(2m)/r^2 * sum_{j=1}^{m-1} binom(m,j)/(q')^j
       = O(m K^2/n) = o(1),

where q'=r theta=(1+o(1/m))n. The disjoint tuple indicators are independent. The variance of the all-new count is at most its mean plus this overlapping-pair expectation, hence O(K). Chebyshev and Markov give at least K/2 all-new tuples with pairwise disjoint indexed labels at the fixed identity value with probability 1-o(1). Conditional on the chosen point sets, #12's random c-labels on these disjoint tuples are independent uniform residues. Their scores fail to vary with probability at most r^(1-floor(K/2)). Therefore the identity-fiber rank-increasing contrast exists and is selectable by full short-fiber enumeration on a probability 1-o(1) portion of the very augmentation family that defeats both prior criteria. This proof is consistent with intersecting all of #12's probability events. It is not a universal statement for every existential block family or every actual mask output. The same elementary argument with independent uniform vector labels in F_r^d and M disjoint tuples bounds failure to affinely span the quotient by r^(d+1-M); it does not supply independent random quotient labels for the actual extractor.

PR #12's asymptotic arbitrary augmented sets still lack a demonstrated affordable membership network, so this is a count-query resolution of that construction, not an asserted bit-time algorithm for those sets. The new targeted affine seed of sections 1–3 does have the ordinary arithmetic membership circuit. Its modular-rank guarantee remains (6), not the independent-random-label conclusion.

## 5. Faster globally pointwise seed via batched mask evaluation

Independently of the targeted second-moment construction, revisit POINTWISE_SEEDS.md's centered high-moment seed. There the same degree-h polynomial was evaluated separately at up to h private x coordinates in each parity and block. This cost O(h^2) field operations per block, giving the 1/8 seed exponent. For any fixed epsilon'>0 we can compile **all** evaluations using O(h^(1+epsilon')) field operations on the same field, with no variable inversions.

For completeness choose fixed integer k>=2 with log_k(2k-1)<1+epsilon'/2. For characteristic p>2k-2 (eventually), split coefficient arrays into k chunks, evaluate the chunk polynomial at 2k-1 distinct public constants, multiply recursively, interpolate with the fixed inverse Vandermonde matrix, and combine shifted chunks. This is ordinary Toom-Cook with recurrence (2k-1)M(ceil(n/k))+O_k(n). Padding to powers of k gives O(n^(log_k(2k-1))) field operations. The interpolation uses only constants precomputed from p in polynomial bit time. Form the balanced product tree of (X-x_j). Descend with monic remainders; for division by monic polynomials reverse the divisor, invert its constant-one series by Newton iteration g <- g(2-fg) mod X^(2s), and reverse the truncated quotient. Subtract quotient times divisor. Degree bounds can be fixed in advance, padding with zeros; no conditional degree scan or private inversion is needed. The product and remainder trees cost O(M(h) log h) field operations, thus O(h^(1+epsilon')). Repeated evaluation points are allowed. For the degree-h mask and h leaves, the root remainder is just ordinary monic reduction, including a possible leading term. Evaluate both parity polynomials on the same list of x values and select the prescribed output using the parity bit. For j<h in a moment query pad with public dummy x coordinates and ignore their results. All inputs and polynomial coefficients have canonical b-bit field representations; wires in the arithmetic straight-line program have bounded-fanout copies implemented by equality trees. This does not replace the contractor by a fast dense field-linear operation inside a high-degree tensor.

Consequently the h-copy count on m blocks has O(m h^(1+epsilon') b) Boolean edges under the inherited linear-edge arithmetic compilation, O(b^(3/4)) maximum local degree, and count precision O(m h^(1+epsilon') b). Network generation and public constants have polynomial cost. This is a circuit-size claim for simultaneous evaluations at variable points, not the usual preprocessing assumption that the points are fixed. The underlying product/remainder-tree and Toom-Cook methods are classical; the application to the shared private mask is the contribution. Background reference: https://flintlib.org/doc/algorithms/polynomials/subproduct-tree.html (monic remainder tree, repeated points permitted). The finite test uses Toom-Cook and Newton monic division, not a tensor contractor.

For a fixed nonnegative integer s, choose q* = ceil((Kr)^(1/(m-s))) with K=2^floor(b/(64m)), h=1024m even, and delta=2^-floor(b/(1024m)), as in POINTWISE_SEEDS.md except for the redundancy change. Take the sum of its nonnegative normalized h-th centered potentials over every retained block subset of size at least m-s. There are O(m^s) such subsets. Each has nominal mean at least K and q'>=K, so the centered-partition proof in that note bounds the total expectation by 2^(-3b+o(b)); the same binomial expansion and two ambient mean values apply to each subset. Conditional expectation with exact ordinary counts fixes a seed with relative pointwise error O(delta) simultaneously for all group values and all these retained subsets. The full fiber size is O(K q^s), with log size O_s(b/m); enumeration is therefore affordable. For s>=2 the same one- and two-slot marginal argument as section 3 gives full real contrast rank on every prescribed fiber (in fact the relative pin errors themselves are uniformly small). We do not require universal ordered-pair pin no-wrap for s>=1, which conflicts with the larger redundancy as noted in POINTWISE_SEEDS.md.

With h=Theta(m), the improved seed and enumeration time envelope is

    log2 T = O(b/m + m^(2+epsilon') b/sqrt(log b)).

Balancing at m=Theta(L^(1/(6+2epsilon'))) gives T,S=2^{O(b/L^delta')} for delta'=1/(6+2epsilon'). Thus any fixed delta'<1/6 is available for the **globally** pointwise seed with any fixed number of deletions and enumerable individual fibers, improving the previous 1/8 seed exponent. It is not a complete deterministic ECDLP exponent. The new targeted second-moment construction is faster if only a prescribed subexponential collection of anchors is needed.

## 6. Exact finite verification and success audit

Run `python -B check_fibers.py`, `python -B check_seed_prefix.py`, and `python -B check_prior_obstruction.py /path/to/PR12/research/dual_lift_obstruction`. Requires NumPy; the imported #12 fixture also needs SymPy. `verification.txt`, `prefix_verification.txt` and `prior_verification.txt` record the outputs. The tests use exact integer/rational arithmetic, and a fraction-free positive-definiteness check, not floating-point eigenvalue judgments. The actual small order-seven curve prefix experiment tests exact conditional expectation via convolutions of integer joint-membership tables, checks the parent-child weighted averages, and verifies every enumerated tuple using the curve group law. It cannot satisfy the enormous asymptotic K guard, and is not presented as an implementation of the contractor. The six dense controls satisfy the explicit numerical covariance guard and pass both exact matrix inequalities; all 66 moderate-density synthetic fibers attain maximal modular contrast rank, whereas none of the 22 sparse controls do. This is encouraging but not evidence of a universal modular-rank theorem. All three #12 finite fixtures are repaired by a single contrast from the identity fiber (16,115,6 tuples, respectively). The polynomial checks use arbitrary field inputs and repeated evaluation points; the Python routine trims polynomials for convenience, whereas the proof describes a fixed-shape padded straight-line circuit.

**Strongest new result:** a deterministic seed and explicit complete prescribed-fiber measurement architecture at quarter-power cost, with uniform real covariance and exact checkable modular/target-rank tests; a second component improves globally pointwise seed cost to any fixed saving below 1/6. The short-fiber enumeration also bypasses both sufficient-criterion failures on the random #12 examples in the stated count-query model.

**Missing load-bearing lemma:** prove that an affordable deterministic collection of explicit short fibers for the actual selected affine or higher-degree masks has sufficient modular rank or contains an affine target certificate. We have not established this from conditional covariance, prescribed fiber cardinality, disjoint-label matching, or the seed potential. The upstream contractor is assumed, not independently machine-verified here. Preliminary literature review covered Freiman homomorphisms and the classical multipoint evaluation algorithm, not priority of the combined prescribed-fiber theorem. No full independent referee review. No claimed new end-to-end delta, and the mandatory success checklist does not pass. The new work changes the explicit representation and cost of the rank-selection problem; it is not another assertion that near-complete rank implies target identifiability.
