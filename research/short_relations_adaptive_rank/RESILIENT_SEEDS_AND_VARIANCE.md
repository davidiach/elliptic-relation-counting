# Deletion-resilient seeds and variance-sensitive separation

Follow-up to PR #4 (`ANCHORED_TRANSLATIONS.md`, parent `14a87f197a24bc01d2c9eb13c4e64ef0e50c2352`). 27 November 2026. This is AI-assisted conditional mathematical research, not an independent referee report. The uniform faster deterministic end-to-end exponent remains **open**. Use the closed Boolean contraction premise of Liu–Zhang–Ren (LZR), Theorem 2.2, ePrint 2026/3020. The finite tests do not implement that contractor.

## 1. Joint deterministic seed selection for all punctures

Let G=<P> be the whole prime-order elliptic-curve group of odd prime order r, over F_p, p>3, with b=ceil(log2 p), L=log2 b. Let m grow as a positive fixed power of L (m=Theta(L^(1/4)) for the original target). Fix any deletion fraction 0<=alpha<1/2. Choose a fixed integer h>=2 satisfying

    1/(1-alpha) < 2-1/h.

This is possible for every fixed alpha<1/2. Let a=m-floor(alpha*m), K=2^floor(b/(32m)), n=r-1, q*=ceil((Kr)^(1/a)), T=ceil(p q*/n), theta=T/p, q=n theta, q'=r theta. The following proof applies for sufficiently large b, so that 1<=T<=p and K>=2. A bounded initial range can be handled by exhaustive group arithmetic, with constants depending on alpha,h and the fixed arity schedule.

Replace each LZR parity-specific affine map with a polynomial over F_p of degree h, with nonzero leading coefficient, independent uniform lower h coefficients, and independent seed for each block and parity. Retain an affine point when the polynomial on its abscissa has canonical value in [0,T). Exclude O. In each parity class x is injective. At any j<=h distinct abscissas the j evaluations are independent uniform field elements, even after fixing the leading coefficient: the lower-coefficient evaluation map is surjective by Vandermonde rank. Opposite parities are independent. Thus any j<=h distinct points of the ambient nonidentity set are jointly retained in a block with probability theta^j. Every mask contains at most 2hT points, since the polynomial is nonconstant and every value has at most h preimages. All seed words, membership checks and modular arithmetic have a closed BTN encoding using O_h(mb) edges, degree O_h(b^(3/4)), with constant h. This is a different seed family from the affine masks of LZR and the earlier PR.

For I a retained block set with s=|I|>=a, write d_I(g) for the count of one-point-per-block tuples summing to g, W_I=sum_g d_I(g), mu_I=q^s, and C_{h,I}=sum_g d_I(g)^h. The claim is that one deterministic seed makes **all** these retained subsets satisfy, uniformly,

    W_I = (1+o(1)) mu_I,
    r^(h-1) C_{h,I} / W_I^h = 1+2^(-Omega(b/m)),
    r sum_g (d_I(g)/W_I - 1/r)^2 <= 2^(-Omega(b/m)),
    0 <= d_I(g) < r  for every g.

In the second line the error term denotes a nonnegative quantity bounded above by 2^(-Omega(b/m)), not a matching lower bound. All implicit constants may depend on fixed alpha,h.

### Higher collision-moment bound

An ordered h-tuple of tuples contributing to C_h has an equality partition pi_i of [h] in each block i: indices in the same class selected the same point. Let ell_i be its number of classes, ell=sum ell_i, and R the rank over F_r of the equality-of-sums equations on the ell class-point variables. In abstract group coordinates these are h-1 homogeneous linear equations (compare the first h-1 sums with the last). No logarithm is computed by the algorithm. The indicator columns of each single partition are independent, and their span contains the all-ones vector. The difference map has kernel exactly that vector's span. Hence R>=max_i(ell_i-1), and R<=h-1. Ignoring the nonidentity and distinctness conditions for the upper bound, the number of solutions for this partition pattern is at most r^(ell-R). For each actual distinct-point assignment the joint retention probability is theta^ell. Its contribution to E C_h is therefore at most (q')^ell/r^R.

The pattern of all singleton classes in every block has R=h-1 and contributes at most M=(q')^(hs)/r^(h-1). If R=h-1 but some class merges, the ratio to M is <=1/q'. If R<h-1, then ell_i<=R+1 for every i, and the ratio is <=(r/(q')^s)^(h-1-R)<=1/K, since q^s>=Kr. There are at most B_h^s patterns, where B_h is the h-th Bell number. Consequently

    E C_{h,I} <= (q')^(hs)/r^(h-1) * [1+B_h^s(1/q'+1/K)].

This estimates all coincidence patterns, including duplicate tuples; no independence across repeated points has been assumed. Since q'/q=r/(r-1), the normalized leading bound differs from one by O_h(m/r). The loss from B_h^m is harmless because h is fixed and b/m grows much faster than m.

### One potential for all deletion patterns

Let

    Phi_I = (W_I/mu_I-1)^2 + (r^(h-1) C_{h,I}-W_I^h)/mu_I^h.

Both terms are nonnegative (power mean). By pairwise independence inside a block and full independence across blocks, E W_I=mu_I and Var(W_I)/mu_I^2 <= (1+1/q)^s-1. Jensen gives E W_I^h>=mu_I^h. The previous bound then gives

    E sum_{|I|>=a} Phi_I <= O_h(2^m * (m/q + m/r + B_h^m(1/q'+1/K)))
                         = 2^(-Omega(b/m)).

The equality denotes an asymptotic upper bound. There are at most 2^m subsets, not 2^N. Apply exact conditional expectations to the sum, fixing the canonical seed words bit by bit. To compute each conditional mean use exact counts of W_I, W_I^2, W_I^h and C_{h,I} over valid seed completions; the latter uses h tuples with their group sums equated. Unused seed words can be counted separately, or included with range and prefix checks. The valid completions count is the product of the allowed integer-prefix interval sizes. Polynomial masks require O_h(1) modular operations per membership. The full group law, canonical point validity, parity branch and equality gadgets are the LZR arithmetic encodings. Each count has O_h(mb) edges, degree O_h(b^(3/4)), and O_h(mb)-bit exact output; constants absorb h copies. All rational averages can be compared exactly. The different denominators depend on s=|I| and rational q, so a common denominator has polynomial bit length in m,b. Discard empty prefix branches. The conditional expectation decreases, producing a seed whose nonnegative sum is at most its initial mean eta=2^(-Omega(b/m)).

For each I, W_I/mu_I lies within sqrt(eta) of one, so W_I>0 eventually, and the normalized h-th moment ratio is <=1+eta/(1-sqrt(eta))^h. Power means applied to the nonnegative density r d_I/W_I bound its normalized second moment by the h-th moment raised to 2/h. This proves the chi-square bound. The maximum count is at most C_{h,I}^{1/h} <= (1+o(1)) W_I / r^(1-1/h). Since q=q*+O(1), Hasse gives

    q^m = r^(1/(1-alpha)+o(1)),
    max_{I,g} d_I(g) <= r^(1/(1-alpha)-1+1/h+o(1)) < r.

The strict inequality is exactly the choice of h. No modular wrap can occur in these single-tuple fiber counts. A zero fiber contributes 1/r to chi-square, hence each retained support has density >=1-2^(-Omega(b/m)). For alpha=1/4, h=2 suffices and the count exponent is at most 5/6+o(1). For alpha=0.45, h=6 suffices (higher constants); for every fixed alpha<1/2 a fixed h suffices. No claim uniform in growing h is made.

### Bit resources

Seed length O_h(mb); conditional-expectation calls O_h(2^m mb), with polynomial overhead. The factor bases can be enumerated by the same exact prefix counting as LZR (membership of a specified block); no arbitrary lookup of a seed-dependent point list inside the graph is needed. Each block size <=2hT=2^{O_alpha,h(b/m)}. All point-count setup and enumeration are charged. With the assumed closed contractor, total log time is

    O_alpha,h(b/m + mb/sqrt(log(mb)) + m + log b).

All storage is at most 2^{O_alpha,h(b/m+log b)}. The original tuning m=Theta(L^(1/4)) yields the same quarter-power bit-time bound. This is a structural improvement of the seed theorem, not a numerical leading-constant improvement. Degree h fixed can substantially increase constants, and the threshold is asymptotic.

### Near-optimality of the half-deletion threshold

For any arbitrary blocks, if all retained subsets of cardinality a have support at least (1-epsilon)r, then W_I>= (1-epsilon)r. Multiplying the inequalities over all size-a subsets (each index occurs in fraction a/m of them) gives W_full >= ((1-epsilon)r)^(m/a). If all full fiber counts are integers <r, then W_full <=r(r-1)<r^2. For any fixed a/m<1/2 and epsilon=o(1) these requirements contradict each other as r grows. Even at a/m=1/2, demanding W_I/r grow without bound (as in our redundant dense construction) is impossible. This is a cardinality obstruction to deletion-resilience plus **no-wrap**, not a lower bound on discrete-log algorithms or on constructions allowing larger counts.

### Deterministic support selection for the punctured fibers

Any predeclared consistent partial assignment of at most floor(alpha m) full-tuple slots leaves a retained subset I; completion counts are precisely translates of d_I. Thus all such patterns, including ones with a growing number of pinned slots, have the same density/no-wrap guarantees. For a known translation t0 and pattern, use penalty ((d_I(t+t0)-W_I/r)/(W_I/r))^2. Its uniform-t expectation is the chi-square bound and it equals one on failed support. The prefix counts for one or two copies of I, plus the validity count for t, are ordinary closed contractions. The exact conditional-expectation bundle selector in ANCHORED_TRANSLATIONS.md now applies to *single tuples*, not just pinned displacement pairs. For J=poly(N,b) explicitly enumerable bundles of poly(b) patterns each, it produces O(1) starts covering every bundle, since the failure penalty expectation decays as 2^(-Omega(b/m)). This charges all generation and evaluation; it does not enumerate all exponentially many possible long assignments. The selector only guarantees nonzero support, not nonzero candidate discrepancy.

## 2. A one-sided punctured-mean rank certificate

As an application of the joint seed, write d(g),v(g) for full one-per-block counts and integer first incidence moments, and d_-i(g),v_-i(g) for the subset omitting block i. This section requires alpha>0 so at least one deletion is covered eventually. Pad vectors to the full indexed factor columns. For A in F_i define a homogeneous row

    R_i,A(g) = d(g) [v_-i(g-A)+d_-i(g-A)e_i,A] - d_-i(g-A)v(g).

Every row annihilates the true log table lambda: both normalized means would have value log_P(g) whenever the denominators are nonzero. This follows by direct multiplication without division at zero counts. Each supplied row is computable by ordinary single-tuple counts with at most one extra pin for a first moment. It uses two single-tuple copies for its prefix sums, rather than two pairs. The Gram uses four single tuples. Both signs and scalar products with candidate coefficients are computed outside the graph.

**Rank theorem.** The kernel of the full stack over all i,A,g consists of u_i,A=kappa log_P A+c_i. If some nonaffine u survives the previously selected equations, the set of g where at least one of these rows fails on u has density >=(1/(3m)-3epsilon)/2, where epsilon bounds the chi-square of every retained subset. This is >=1/(12m) eventually. Restoring LZR's explicit doubled consistency and normalization rows (possible because the dense full sumset plus itself covers G) gives full column rank, including P,Q.

Proof: at a good full fiber define f(g)=v(g)^T u/d(g) in F_r; similarly f_-i(x) on the support of the punctured count. Suppose the set of g with any row violation has proportion eta. At a starting x where d_-i(x), d(x+A), d(x+B) are nonzero and both full endpoints avoid that bad set, the rows yield f(x+A)-f(x+B)=u_i,A-u_i,B. Each difference translation therefore has its predicted increment on a set of starts of density >=1-3epsilon-2eta. The support of the full sum count is >r/2, so each s in G is a sum of m within-block differences. Telescoping gives each such predicted increment on density >=1-m(3epsilon+2eta). If 3m(3epsilon+2eta)<1, any two representations of s have a common valid start; their predicted increments agree. Direct and composed translations have intersecting domains by the same bound, so the common increments define a homomorphism G->F_r. This forces all within-block label differences to be kappa times the true group log difference. Conversely any such affine labels annihilate the rows, by the fiber identity. For u in the homogeneous kernel of the explicit consistency equations, block constants satisfy c_i=2 sum_j c_j, and r>2m-1 forces all c_i=0. The normalization at P fixes kappa=0; consistency fixes Q. This proves rank and the density lower bound.

This rank proof is simpler than the cross-multiplied displacement means in RESEARCH.md and uses fewer tuple slots per query. It does **not** select an affordable spanning subfamily from the exponentially many group values. The simultaneous support theorem does not remove this residual modular cancellation problem.

**Aggregation obstruction.** On the abstract cyclic group C_101, take blocks (in known scalar coordinates for the finite test only)

    [2,14,39,61,72,82,85,89]
    [3,4,17,24,33,42,44,63]
    [1,3,7,18,32,41,42,44]
    [11,37,45,54,58,72,81,87].

All full fiber counts lie in [28,53], and all four punctured distributions have full support with positive counts below 101. The complete stack has the predicted rank 32-4-1=27, but replacing every family R_i,A(g) by its sum over g has rank 26. This demonstrates that full support and no-wrap alone do not guarantee noncancellation under aggregation, even for the simpler one-sided system. These are illustrative blocks, not an asserted counterexample to the much stronger asymptotic joint seed potential. The verifier checks the complete ranks on 32 qualifying full-support no-wrap fixtures.


## 3. Variance-sensitive exact polynomial separation

This is an alternative sufficient condition for deterministic rank progress that does not require arbitrary label lookups. It works on the doubled-tuple constant-gap family (known scalar shifts and point shifts of the indexed columns), using the present or earlier seed. It is a restricted fast path, not a universal derandomization.

Let u!=0 be a supplied homogeneous direction in the kernel of already collected genuine rows. For doubled tuples X,Y containing two points from each block, with displacement S, the homogeneous discrepancy is Lambda_u(Y)-Lambda_u(X)-c, where c=0 for a known scalar shift and c=u_j for point shift A_j. In the true inhomogeneous equations the right side is the known scalar for a scalar shift, or the dedicated variable moved to the left for a point shift. We will extract one genuine row with nonzero dot product with u. The doubled distribution nu has pointwise bounds (1-epsilon)/r <=nu(g)<=(1+epsilon)/r by Cauchy-Schwarz on the convolution of the centered single-tuple law. Take epsilon<=1/16. The conditional distribution of the initial group value at any fixed displacement has density at least 1/2 relative to uniform; each individual tuple marginal has density at most 2 relative to its independent unconditioned tuple law.

For completeness, the fixed family has a test with violation rate >=1/8. Apply the small-bias scalar-shift argument in RESEARCH.md to the false table lambda+u, or equivalently to the conditional distribution of Lambda_u(X) at each group value. If every actual violation rate were <1/8, every uniform-group-value mismatch would be <1/4. The scalar small-bias bound implies the squared l2 norm of the average conditional error distribution exceeds 2/3. A modal error has mass >2/3. For a wrong column, both endpoint errors equal that mode on a uniform-group probability >1/3, hence actual violation >1/6. Contradiction. All shifts and multipliers below use explicit polynomial-size small-bias lists as in RESEARCH.md; no true logarithms are needed by the separator.

For a multiplier z in a nontrivial-character-bias <=1/4 multiset, choose rational representatives x_i,A of z u_i,A/r modulo 1 in any interval of length one per block (for example the canonical interval [0,1)). Different cuts can be chosen deterministically by scanning the explicit tables; it suffices to use any specified cuts. Let V_z be the sum over blocks of their within-block variances under uniform points. Each mean cancels from the signed sum since the multiplicities in X,Y agree. For the scaled right side choose c_z in [-1/2,1/2] congruent to z c/r. Write

    D_z = sum_Y x_i,A - sum_X x_i,A - c_z.

The phase score 1-cos(2 pi D_z) vanishes on legal homogeneous rows and is in [0,2]. For every nonzero original residue its average over the multipliers is >=3/4. Fix a variance budget v. If at least 7/8 of the multiset (counted with multiplicity) has V_z<=v, then for each illegal row the total phase score averaged over the whole list but restricted to such multipliers is >=3/4-2/8=1/2. The test with violation rate >=1/8 thus has some low-variance multiplier with conditioned expected phase score >=1/16. This condition is directly checkable from u, without needing to estimate any transcendental expectation. More generally the fast path works whenever a tried low-variance multiplier has that gap.

### Global rational minorant

Let delta=1/4096. Choose a certified rational omega with |omega-2 pi|<=delta/(4m+1), 0<omega<6.5, using for example the alternating series in Machin's arctangent formula. All actual |D_z|<=4m+1. For even k set

    P_k(x) = sum_{j=1}^k (-1)^(j+1) (omega x)^(2j)/(2j)! - delta.

The alternating Taylor polynomials are global bounds, not merely bounds on a small interval: the error C_k(x)-cos(x) for C_k=sum_{j=0}^k (-1)^j x^(2j)/(2j)! has second derivative -(C_{k-1}(x)-cos(x)), and vanishing value and derivative at zero for k>=1; start with C_0-cos>=0 and use evenness. This proves for even k that 1-cos(omega x) minus the unshifted P_k lies between zero and the magnitude of its next term. The Lipschitz estimate for cosine then gives P_k(D_z)<=1-cos(2 pi D_z) for every actual row. In particular P_k<=0 on every legal row. The global bound, including tails, is essential for exact prefix extraction.

### Degree depends on variance, not arity

For a centered slot label Z, |Z|<=1, E Z=0. Taylor expansion gives E exp(sZ)<=exp((exp(|s|)-1-|s|)Var Z). The centered sum U_X of a doubled tuple has E exp(s U_X)<=exp(2(exp(|s|)-1-|s|)V_z) under the unconditioned independent law. Conditional marginal domination by 2 and Cauchy-Schwarz on exp(tU_Y)exp(-tU_X) (and its inverse) imply

    E_cond exp(t |D_z|) <= 4 exp(t/2+2(exp(2t)-1-2t)V_z)   (t>0).

At t=7, the expected Taylor remainder is bounded by

    E_cond |omega D_z|^(2k+2)/(2k+2)!
      <= (13/14)^(2k+2) * 4 exp(7/2+2(exp(14)-15)V_z).

For example take the even k=2 ceil(2^25 (v+1)); the crude bounds exp(14)<2^22 and (13/14)^(2k+2)<=exp(-(k+1)/7) make the remainder <delta. The constants here are deliberately very loose. Therefore a tested multiplier with expected phase score >=1/16 has E_cond P_k(D_z)>=1/16-3delta>0. A polynomial of degree O(v+1) suffices. No small-variance assumption has been made about the true log labels.

### Algorithm and complexity

For each test shift and eligible multiplier compute the **exact signed rational sum** of P_k over valid doubled tuple pairs. Expand the powers of D_z into signed terms with at most 2k slot occurrences; for every specified occurrence enumerate its indexed point, and count the unweighted arithmetic relation with these pins. Repeated pins on one slot must agree, and otherwise the term is zero. This is one ordinary pair network per summand, not k copies of it. Apply candidate coefficients outside the graph; clear denominators exactly (including r, omega, factorials). For a positive total, prefix restrict the coordinate encodings and retain a child of positive sum until a single pair remains. Its positive score implies a nonzero homogeneous discrepancy; append its genuine row, increasing modular rank. No approximate oracle, random choice, or unspecified high-arity lookup is used in this fast path.

The arithmetic graph has O(mb) edges, O(b^(3/4)) degree; even with many pins the same fixed slots are pinned only once. Number of expanded summands <=(O(m)(N+1))^{O(k)}; rational coefficient and count bit lengths are polynomial in k,m,b. Thus at a supplied direction with v<=m^theta, 0<=theta<1, log time is

    O((v+1)b/m + mb/sqrt(L) + polylog(b)).

At the original arity this gives a saving exponent (1-theta)/4. Retune m=Theta(L^(1/(2(2-theta)))) and use the same seed construction at that arity; balancing the two terms gives the stronger restricted exponent

    delta_theta = (1-theta)/(2(2-theta)).

For example theta=0 gives 1/4, theta=1/2 gives 1/6, theta=3/4 gives 1/10. These are sufficient bounds for **rank progress at a supplied low-variance direction**, not a universal deterministic DLP bound. Up to N such increments give the same scale if at each iteration a suitable direction is obtained by charged deterministic means (e.g. try a computed kernel basis and accept a qualifying vector). No proof of that condition on every residual kernel is supplied. In particular uniform-looking residues can have variance proportional to m for every tested multiplier. The certificate is different from the modal-error certificate of RESEARCH.md: clustered scalar labels in a direction need not imply clustered group error points. It does not supersede that other fast path.

## 4. Obstruction to a concentration-only universal moment argument

The following exact probability-measure counterexample explains why merely knowing that the normalized discrepancy has subgaussian width O(sqrt(m)) cannot give a universally positive low-degree polynomial minorant. This strengthens the earlier unstructured moment-matching example with a concentration promise, but is **not** an example realized by the prescribed curve relation family.

Let d be a positive multiple of four, r an odd prime >=16d, and use the d+1 strictly increasing rational nodes

    x_j = j/2 - (j mod 2)/(2r) - d/4,  0<=j<=d.

The even nodes are integers, the odd nodes are noninteger multiples of 1/r. In fact r x_j is always an integer. Put w_j=1/prod_{l!=j}|x_j-x_l|. Lagrange interpolation of any polynomial of degree <d, and comparison of its leading coefficient of degree d, gives sum_j (-1)^j w_j f(x_j)=0 (d even). In particular the positive even and odd weight sums are equal. Normalize each parity separately to probability measures mu_even and mu_odd. Their moments match through degree d-1; all even nodes are legal normalized discrepancies, all odd nodes illegal.

For the equally spaced nodes the weights would be w_j^0=2^d/(j!(d-j)!). Every spacing ratio differs from one by at most 1/(r|j-l|). Using |log(1+t)|<=2|t| for |t|<=1/2 gives |log(w_j/w_j^0)|<=4 H_d/r<=1/4. After parity normalization the probability ratios relative to the binomial(d,1/2) law conditioned on the corresponding parity lie in [exp(-1/2),exp(1/2)] subset [1/2,2]. Each measure is thus dominated by four times the unconditioned binomial law; its node differs from (j-d/2)/2 by <=1/(2r). Hoeffding gives, for either measure and t>=0,

    Pr(|X| >= t+1/(2r)) <= 8 exp(-8 t^2/d).

This yields subgaussian tails with variance scale O(d) and support inside [-d/4,d/4]. Nevertheless any polynomial of degree <d nonpositive on the legal nodes has nonpositive expectation under the all-illegal measure too, by exact moment equality. The proof uses rational weights; the required integer multiplicities need not fit the actual relation-row count. It rules out a degree-o(m) universal separation proof based only on this concentration envelope and scalar moments when d=Theta(m), not on more detailed arithmetic/conditional independence properties or other algorithms.

## 5. Verification, scope and next steps

Run `python check_resilient_seed_phase.py` (standard library). `resilient_seed_phase_verification.txt` records the finite checks: evaluation bijections of polynomial masks with fixed nonzero leading coefficient, exact collision-partition ranks and assignment counts, joint potential minimization in three small restricted seed libraries (one uses polynomial parity masks on an actual elliptic curve), rational low-variance polynomial/prefix extraction, exact concentrated moment equality, and one-sided punctured-mean rank and aggregation checks. The seed-library checks are not the asymptotic seed search. The phase fixtures use a smaller even Taylor truncation that succeeds on these toy tables, not the enormous sufficient constant in the proof. They do not demonstrate practical complexity or universal rank success. The probability obstruction is verified at five odd prime orders through degree 63.

The next bottleneck remains uniform modular rank selection. Jointly controlling every punctured subset's high moments removes many support and wraparound failures, but it does not imply that a small selected stack of the mean equations is full rank. The variance-sensitive polynomial provides a different exact progress condition, and the concentration obstruction limits a naive universal extension. Future work should use the exact arithmetic structure of the punctured counts to derive noncancellation, or prove a cheap low-variance annihilator selection lemma for the actual residual kernels. No such lemma is established here. No universal intermediate logarithmic saving is claimed. Independent proof review and novelty assessment remain outstanding.

Sources: LZR https://eprint.iacr.org/2026/3020 (Theorems 2.2, 3.1 and affine mask method); PR #4 https://github.com/davidiach/elliptic-relation-counting/pull/4, RESEARCH.md and ANCHORED_TRANSLATIONS.md at the parent commit. The polynomial limited-independence construction, Vandermonde argument, Lagrange identity, power means, Taylor bounds and conditional-expectation machinery are classical. No novelty claim for those ingredients individually.
