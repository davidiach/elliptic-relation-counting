# Anchored translations and deterministic support selection

Follow-up to RESEARCH.md, PR #4, parent c40ab298a621eafd532752ef0e5d75dcea09e648. 27 November 2026. These are conditional mathematical results on the common prime-order elliptic-curve class, not a completed deterministic exponent improvement. Use the closed Boolean contraction premise in LZR Theorem 2.2. The finite tests below do not implement that contractor.

## 1. Exponentially dense support and pointwise no-wrap pin counts

Use the seed potential in RESEARCH.md section 1, but increase its density parameter to K=2^floor(b/(16m)), where m=Theta((log b)^(1/4)). In the asymptotic range K>=4096m, r-1>=K and the threshold q* in that section is ceil((Kr)^(1/m)); q*=max(Km,ceil((Kr)^(1/m))) is the uniform definition. The bounded initial range can be handled by exhaustive group arithmetic. The mean block size q=(r-1)T/p is q*+O(1); mu=q^m=Theta(Kr). The same expectation calculation and bitwise minimization give Psi<=5/K. Thus W=prod n_i=(1+O(K^(-1/2)))mu and, for the one-per-block sum count d(g) and probability p0=d/W,

    epsilon := 20/K,
    r sum_g (p0(g)-1/r)^2 <= epsilon,
    C := sum_g d(g)^2 = O(K^2 r).

All parameters and rational comparisons have polynomial bit length. Since each parity-split block has size n_i<=2T=O(q), the product lower bound gives min n_i >= W/(2T)^(m-1) >= 2^(b/m-O(m+1)). Consequently K^2/sqrt(min n_i)=o(1). The enumeration and seed-construction time remains 2^{O(b/(log b)^(1/4))}; the number N of indexed points, including P and Q, is 2^{O(b/m+log b)}. The stronger density does not change the asymptotic scale. The earlier short-row count changes from poly(m)r to Theta(K^2 r); this is still less than r^2 for large b. The original smaller K is sufficient for the rank theorem in section 3 below.

Consider a pair of short tuples X,Y with displacement t=val(Y)-val(X). Pin a point A in one slot i of Y. Let n_h(t) be the integer count, and let v_h(t) be its summed signed incidence column. Its total count M_h=W^2/n_i. Writing d_-i for the unnormalized convolution of the other blocks, n_h(t) is a translate of d_-i convolved with reflected d. Thus Cauchy-Schwarz and nonnegativity give

    n_i ||d_-i||_2^2 <= ||d||_2^2=C,
    0 <= n_h(t) <= C/sqrt(n_i) < r     for every t.

The same holds for a pin on the source side. In particular nonzero integer mass is exactly nonzero field mass; no modular wrap anywhere. The normalized displacement distribution contains the full unpinned opposite half as a convolution factor, hence its chi-square distance from uniform is <=epsilon (convolution contracts centered l2). At most epsilon r displacements have mass zero.

More generally, pin any consistent nonempty subset of slots on a single side, possibly the entire tuple. Its mass is at most a one-pin mass; the complete opposite half is still present. These two conclusions therefore hold with total mass W^2 divided by the product of the pinned block sizes. This is a uniform bound for each such pin pattern, not a claim that every pattern has support at every displacement.

## 2. Deterministic simultaneous support selection

Here is a constructive improvement over the union-bound support arguments in RESEARCH.md. A constraint is a specified one-side pin pattern h and a known group translation S. It asks that n_h(t+S) be a nonzero field element. Let mu_h=M_h/r>0 and

    phi_h,S(t) = ((n_h(t+S)-mu_h)/mu_h)^2.

For uniform group t its expectation is at most epsilon. If the constraint fails, its penalty is exactly one (by the pointwise no-wrap result). For L prescribed constraints with L epsilon<1 there is an algorithm finding a common good start t, within the quarter-power bit budget. Compute the conditional expectation of the sum of penalties while fixing the canonical point encoding of t bit by bit. For every prefix, the exact sum of n_h(t+S) is a closed arithmetic count of one pinned short pair; the sum of its square counts two independent pairs sharing the same displacement, translated by S. Include the point-validity conditions and the prefix restriction on t. The number of valid t extending the prefix is another closed count. All networks have O(mb) edges and maximum degree O(b^(3/4)) after the LZR equality/subdivision construction; integer output precision is O(mb). Use exact rational averages, discard empty branches, and take the smaller conditional expectation. At the end the total penalty is <1, so all constraints hold. No candidate labels or field inversions are needed by the counting oracle.

A list version is useful for many more constraints. Let there be J explicitly specified bundles of at most q0 constraints each, with q0 epsilon<1. A bundle is satisfied if all its constraints are satisfied. Apply the same conditional expectation selection to the sum of penalties of every remaining bundle. A failed bundle contributes >=1, so after the selection at most J q0 epsilon bundles remain. Remove covered bundles and repeat. The process takes at most 1+floor(log J/log(1/(q0 epsilon))) selections, interpreting this bound as one when J=1. If J=2^{O(b/m+log b)} and q0=poly(b), the stronger K above yields O(1) starts. Evaluating every prescribed mass at each selected start, and generating all prefix counts, stays in the quarter-power time/space budget. The bound allows any explicitly enumerable polynomial-size list in N,b; it does not enumerate the exponentially many possible tuples or arbitrary paths.

For example it covers all prescribed test shifts for a single anchor, requiring that both its mass at t and at t+S be nonzero, using a constant-size deterministic start list. It also applies to prescribed full source tuples. This is a hitting list for **supports of counts**, not a hitting list for modular discrepancies or determinants. Those are different predicates.

## 3. A one-reference translation certificate with constant gap

Fix just one valid positive-slot pin h. Write n(t),v(t) for its count and incidence vector; pad v with zero dedicated P,Q coordinates. For the true log table lambda,

    v(t)^T lambda = n(t) log_P(t)  (mod r).

This identity is mathematical; neither the count oracle nor the proposed machine knows log_P(t). For a shift S define the row and scalar weight

    B_S(t) = n(t)v(t+S)-n(t+S)v(t),
    c_S(t) = n(t)n(t+S).

Choose the explicit scalar small-bias multiset from RESEARCH.md section 1, with nontrivial character bias at most 1/4. For each scalar a in it, add the family of equations

    B_[a]P(t)^T y = a c_[a]P(t).

For every indexed base point or dedicated point A_j (including P,Q), add the family

    (B_Aj(t)-c_Aj(t)e_j)^T y=0.

These are all genuine equations: multiplying out the true mean identities cancels the unknown log_P(t). There are only O(N+poly(b)) families, each indexed by r displacements. A zero-mass endpoint yields the zero equation because n=0 implies v=0. At a nonzero-mass endpoint, v/n is a field-valued linear form. No arbitrary label lookup, no explicit consistency rows, and no cross comparison between two different pins is required for completeness.

**Theorem.** If the anchor has nonzero mass on at least (1-delta)r displacements with delta<=1/256, then for every false candidate y at least one of these fixed families has a violated equation at >=r/4 displacement values. Thus the complete equation family has full column rank. The family can depend on y; it is not identified by this existential proof. The density is measured over uniform t, not by weighting with the number of tuple pairs. The proof also allows counts with modular zeros if those are excluded from the good set, but section 1 gives the clean no-wrap version.

Proof: on the good set define f_y(t)=v(t)^T y/n(t) in F_r, and let e(t)=f_y(t)-log_P(t). Extend e arbitrarily to the bad set. Let p_c be the proportion of all t with e(t)=c. For a scalar shift a the uniform-t equality probability e(t+[a]P)=e(t) averaged over the small-bias list is at most

    1/4 + (3/4) sum_c p_c^2.

Indeed apply Fourier orthogonality to each indicator of an error value; the averaging operator fixes constants and has all other eigenvalues of absolute value <=1/4. This estimate needs no representation of the unknown error function for the algorithm.

Suppose all families violate at <r/4 positions. Scalar-shift mismatch on the entire group is at most the genuine violation fraction plus 2 delta (endpoints outside support). Averaging implies sum p_c^2 > 1-(1/4+2 delta)/(3/4) >=21/32. Thus a modal error c0 occurs on >21/32 of all group values. Pick a wrong column j, which exists for y!=lambda. For shift A_j, the two error values both equal c0 on >2(21/32)-1=5/16 of all t by the union bound, and both endpoints are good except for a fraction at most 2 delta. At those good equal-error endpoints the conditional means differ by the true log_P(A_j), not by y_j. Hence the point-shift equation violates on >5/16-2 delta>1/4 of all t, contradiction.

The true table is a solution, so applying the theorem to lambda+u for any nonzero kernel vector u proves full column rank. In particular normalization and consistency equations from LZR can be included as useful extra rows but are not needed in this proof. With K=2^floor(b/(16m)) the density assumption holds eventually by section 1; already the original polynomial K suffices. This is stronger than the prior 1/(12m) density over the union of cross-mean families in RESEARCH.md: one fixed shift family has constant-density violations, and only one reference pin is used. The numerical constants are not optimized.

## 4. Evaluation, fast paths and the remaining gap

Each specified B_S(t),c_S(t) requires two anchor count/first-moment queries. First moments expand into signed slot pins; for a given incidence entry only one additional pin is needed. Combine counts outside the graph. The explicit shift points [a]P and A_j cost polynomial bit time per shift. Prefix sums of the rows use two independent pinned pairs whose displacements differ by S; Gram entries use four such pairs, with a shared canonical t, and first-moment incidences expanded as additional pins. All counts are of nonnegative integer terms with precision O(mb); apply signs and reduce modulo r afterward. Thus one evaluation, all family sums, the support selector of section 2, and the adaptive first/quadratic tests of RESEARCH.md fit the target time 2^{O(b/(log b)^(1/4))} and space bound, conditional on the closed contractor. The same holds for any fixed-order higher moment. A moment of growing order k requires O(k) pair copies and N^{O(k)} expanded summands and is not free.

An explicit certified fast path is to start with the normalization and consistency rows (optional), append the translation equations at the support-selected starts, append all family sums over t, and then run adaptive prime-adic rank extraction on each family. At every step solve or update only exact genuine equations, check rank over F_r, and verify any returned log table by scalar multiplication. If full rank is reached this is a deterministic algorithm at the target scale for that input. If the test stalls return the residual kernel; no faster worst-case guarantee follows. The constant-gap theorem gives many individual violated equations but not a deterministic way of locating one, and a nonzero mass at a selected start does not imply a nonzero discrepancy there.

There is a simple obstruction to replacing the complete stack by the unweighted sum in every family. The test on y^2=x^3+4x+9 over F_29, with r=37 and generator (0,3), uses indexed log coordinates for the illustrative blocks

    F_1: 1,3,5,10,18,19,27,32,34,36
    F_2: 18,19,21,31
    target: 35; anchor: first positive slot, point P.

Every anchor mass is between 2 and 7. The complete stack has rank 16, but even using every scalar shift (not merely the small-bias list) the summed families have rank 15. Adding normalization and a specific valid doubled decomposition for every column (first in the script's indexed-log order) still gives rank 15. This is an actual curve-group arithmetic fixture. It does not meet the strengthened seed potential, nor is its decomposition asserted to be LZR's lexicographic coordinate choice; thus it is not a counterexample to a possible theorem that uses those further promises. It does refute universal aggregation on the sole full-support/no-wrap promise. At start t=0 the signed powers of two and the point-shift equations give full rank in this same fixture. That is an experimental rescue, not a general theorem.

For the all-displacement sums there is also a useful algebraic expression. In F_r[G], denote reflection by star, write F_i=sum_{A in F_i} [A], F_i^u=sum u_i,A[A], L_i=F_i^u/F_i (each F_i is a unit since n_i<r), and N=sum_t n(t)[t]. The generating function of the averaged row evaluated on u is

    N N^* sum_i k_i (L_i-L_i^*),

where k_i=1 in the anchored block and 2 elsewhere. This follows by expanding the signed first moment as V_u/N=u_anchor-L_anchor^*+sum_{i != anchor}(L_i-L_i^*) and subtracting its reflection; the constant u_anchor cancels. This identity explains why changing only the pinned point within a block does not change its all-displacement sums. It offers structure for future rank selection, but it does not imply that evaluations at the selected scalar/point shifts preserve rank. The group-algebra units are formal, not cheaply computable dense arrays.

## 5. Verification and current frontier

Run check_anchored_translation.py from this directory (standard library, Python 3.10+). See anchored_translation_verification.txt. It enumerates small cyclic groups and actual order-19 and order-37 elliptic-curve group laws; it does not reproduce the asymptotic seed construction or contractor. Five full-support fixtures test exact pin bounds, normalized energy, complete rank, summed-family rank and a short list at a single start. The 244 false tables include single-coordinate errors and wrong dedicated columns. The minimum observed best violation rate is 0.945946, against the proven 1/4; this finite experiment is not evidence of that rate at large p. The prefix-based support selector is separately exercised on 1,028 translated two-count bundles with actual zeros in abstract count arrays, and 171 bundles on the order-19 curve. The order-37 negative fixture is the exact regression above. All checks are deterministic, including the fixed-seed candidate generation.

The next target is a uniform noncancellation/rank-selection theorem for this one-anchor family, exploiting its special displacement-convolution structure or a stronger seed promise. Support selection alone cannot establish it. A positive theorem for all eligible curves would finish the deterministic quarter-power reduction; an intermediate positive logarithmic saving would also be new relative to LZR if all costs and assumptions are charged. Neither is proved here. The upstream contraction assumptions and literature priority require independent scrutiny.

Sources: RESEARCH.md at the parent commit, LZR https://eprint.iacr.org/2026/3020 (Theorem 2.2, Theorem 3.1 and Lemma 3.2), and the small-bias construction and background reference in RESEARCH.md. Classical Fourier identities and conditional expectations are not claimed as new.
