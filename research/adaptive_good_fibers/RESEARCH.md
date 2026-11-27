# Adaptive robust fibers from actual affine masks

**Research note — 27 November 2026.** This note returns from the abstract additive-avoidance problem to the arithmetic seed and relation architecture. Its main contribution is a single deterministically constructed affine seed that supports *adaptive* queries for small, approximately high-wise-uniform translated fibers, using only the original indexed columns. A uniform rank-or-distance theorem follows. Full modular rank and a faster complete discrete-logarithm algorithm remain unproved.

All resource estimates below are conditional on the closed Boolean contraction theorem (LZR, Theorem 2.2) and the ordinary bounded-degree, linear-edge arithmetic compilation used in the repository's fixed-shift manuscript. No stronger contractor or independent verification of that premise is claimed. All finite tests are small mathematical controls, not cryptographic performance experiments.

## 1. Arithmetic setting and parameters

Let E/F_p be a nonsingular elliptic curve whose full rational group G has odd prime order r. Write b=ceil(log_2 p), L=log_2 b. By Hasse, log_2 r=b+O(1). Take m growing as a fixed positive power of L, integers h>=16 and k>=0 with h+k=o(m), and put

    d0 = m-k-2h >= 2,
    a = floor((h+k+1)b/(128m)), K=2^(2a),
    q_star = ceil((Kr)^(1/d0)), T0 = ceil(p q_star/(r-1)),
    theta = T0/p, q = (r-1)theta, q' = r theta.

Only the sufficiently large parameter regime is claimed; in particular T0<=p. All thresholds can be computed with integer powers and comparisons. We have q=q_star+O(1), log q=O(b/m), q^d0=(1+o(1))Kr, and log K=o(b). Denote by Ical the family of retained sets I of block indices with |I|>=m-k. For each I put mu_I=q^|I|/r, so mu_I>=K q^(2h).

The indexed block F_i is obtained by independently selecting two pairs of affine coefficients (one per canonical ordinate parity), each uniformly in F_p^2, *including zero slopes*. Keep each nonidentity rational point A when its corresponding affine evaluation ax(A)+c, reduced canonically modulo p, is below T0. Within a parity class, distinct rational points have distinct abscissas: the other point with the same abscissa is its negative and has the opposite canonical parity. (There is no nonidentity point of order two.) The two evaluations are therefore uniform independent when their abscissas differ. Across parities the seeds are independent. Thus the inclusion indicators are exactly pairwise independent within each block, each with marginal theta, and independent across blocks. This is the actual affine-mask law already used in PR #13, rather than an abstract random-subset hypothesis.

Write n_i=|F_i|, W_I=product_{i in I} n_i, and d_I(g) for the integer number of one-per-block tuples from I with group sum g. The selected seed will have no empty block. For J subset I with |J|<=h let

    E_(I,J)(g) = sum_{A in product_(j in J) F_j}
                (d_(I\J)(g-sum A)-W_(I\J)/r)^2.
    Phi_I(g) = sum_{J subset I, |J|<=h} q^(2|J|) E_(I,J)(g)/mu_I^2.
    barPhi_I = (1/r) sum_g Phi_I(g).

The empty J is included, with one empty pin assignment. No independence of h different points within one mask is needed.

## 2. A low-order seed potential with high-order marginal consequences

### 2.1 The exact mean inequality

For any retained set R of at least two blocks and any group value g,

    E[(d_R(g)-W_R/r)^2] <= 2 q^|R|/r.                    (1)

This is the second-moment lemma of PR #13, reproduced for clarity. Use unnormalized Fourier coefficients of the mask indicator f_i on G. For nontrivial complex characters chi, psi, pairwise inclusion independence gives

    E[fhat_i(chi) conjugate(fhat_i(psi))]
      = theta(1-theta)(r 1_(chi=psi)-1)+theta^2.

The diagonal absolute value is at most q, and the off-diagonal absolute value is at most theta. Expand the centered convolution at g, use independence between the blocks and the triangle inequality on the off-diagonal terms. The resulting bound is q^|R|/r+theta^|R|<=2q^|R|/r. The random W_R is the exact trivial-character contribution; it is not replaced by its expectation.

Each pin assignment on J depends only on the seed of J, independently of the remaining masks. Apply (1) to I\J, then sum over ambient assignments with the probability theta^|J| of inclusion. Since |I\J|>=m-k-h>=2,

    E[E_(I,J)(g)] <= 2 mu_I,
    E[Phi_I(g)] <= 2 * 2^m / K.                         (2)

The same bound holds for barPhi_I. Thus controlling arbitrarily many selected pin slots (one point per block) uses only two tuple copies in the counting oracle, not h copies.

### 2.2 A short signed-collision constraint

Let t0=floor(3d0/4). Let B count all nonempty signed patterns of total length t<=t0, with at most two distinct indexed points in any block and coefficients in {-2,-1,1,2}, whose group sum is the identity. Distinctness is imposed within each block; points in different indexed blocks need not initially be different. A pattern means choices of occupied slots and their coefficients. For t>=1, fixing t-1 ambient point values leaves at most one choice for the last because its coefficient is a unit modulo odd r. Its expected count is at most (q')^t/r, by the same pairwise independence. There are at most 21^m patterns (1+4+16 possibilities per block). Consequently

    E[B] <= 21^m (q')^t0 / r = r^(-1/4+o(1)).           (3)

This is a statement about bounded coefficients and at most two points per indexed block, not arbitrary linear dependencies in the later quotient.

### 2.3 Deterministic selection

Use the nonnegative potential

    Psi = m^4 sum_i ((n_i-q)/q)^2 + 2^a sum_(I in Ical) barPhi_I + B.

The expected size term is at most m^5/q, by pairwise independence. Equations (2) and (3) give

    E[Psi] <= m^5/q + 2^(2m+1-a) + r^(-1/4+o(1)) < 1/256

for sufficiently large b in the stated regime. The same proof permits predeclared anchors: add 2^a sum_(I,g in S) Phi_I(g) for any explicit set S of size at most 2^(a/4). The expectation of that extra term is at most 2^(2m+1-3a/4). This optional term is not needed in the adaptive theorem.

For the exact conditional-expectation implementation use

    r^2 E_(I,J)(g)
      = sum_(X,Y: X_J=Y_J) (r [val X=g]-1)(r [val Y=g]-1).   (4)

Expand into four ordinary nonnegative integer counts; signs and rational coefficients are applied outside the network. In a sum over g the g variable must remain present even in terms with no equality constraint on it, to preserve its multiplicity. The size penalty uses one- and two-point counts, and B uses ordinary group-law constraints and within-block inequalities with at most two occurrences per block. Include the variable mask coefficients, their canonical range constraints, and any fixed prefix bits in each count. Normalize by the number of seed completions, computable from those range constraints. The inherited arithmetic compilation uses O(mb) edges per count, local degree O(b^(3/4)), precision O(mb); the number of terms is 2^O(m) times the number of predeclared anchors plus one. There are O(mb) seed bits. The rational weights, sums and divisions are exact and have polynomial bit length. At each nonempty prefix choose a child with no greater conditional expectation. Empty prefixes are skipped.

We obtain an actual affine seed satisfying

    |n_i/q-1| < 1/(16m^2),
    3/4 <= W_I/q^|I| <= 5/4  (all I),
    sum_I barPhi_I < 2^(-a)/256,
    B=0.                                                   (5)

The product bounds follow from the size guard. The zero integer B follows from Psi<1. In particular, two distinct tuples that differ in at most floor(t0/2) slots cannot have the same group sum: cancellation of their identical slots would give a forbidden signed pattern. Every partial one-per-block sum map on at most floor(t0/2) selected blocks is injective. This extra arithmetic constraint is built into the same seed, not assumed of the output of an arbitrary additive avoidance instance.

## 3. Adaptive robust-fiber access on the fixed seed

Put tau=2^(-h-14)/m^2 and delta=8 tau. Call t good for I if Phi_I(t)<tau^2. For every J of size at most h, and every specified pin assignment A, goodness gives

    |d_(I\J)(t-sum A)-W_(I\J)/r| < tau mu_I/q^|J|.

Using (5), this error is at most 2 tau times the nominal pinned count. The same estimate with J empty bounds d_I(t), in particular

    mu_I/2 <= d_I(t) <= 2 mu_I.

Dividing the pinned count by the total gives, for the uniform distribution P_t on the actual fiber,

    (1-delta)/product_(j in J)n_j
      <= P_t[X_J=A] <=
    (1+delta)/product_(j in J)n_j.                        (6)

This is pointwise multiplicative approximate h-wise uniformity, simultaneously for all selected pins. It implies every pin assignment of strength <=h has a completion. Good fibers contain at most 2 mu_I tuples; log mu_I=O((h+k+1)b/m). No bound of this size is asserted for every untranslated fiber.

For a retained set I define p_I=d_I/W_I. The empty-pin term gives

    chi_I^2 = r sum_g(p_I(g)-1/r)^2
      <= (q^|I|/W_I)^2 barPhi_I <= 4 barPhi_I.

Thus p_I*p_I has every point probability <=(1+chi_I^2)/r <=2/r. More generally, take a uniform doubled tuple T from I (two independent one-per-block tuples), possibly with some slots preset so that the remaining free index set of EACH copy is still in Ical. Its sum law is a translate of p_U*p_V for those free sets. Cauchy-Schwarz gives the bound (1+sqrt(chi_U^2 chi_V^2))/r <=2/r for each point of that law.

**Adaptive translation theorem.** After the seed has been fixed, choose any retained I and any explicit batch g_1,...,g_ell of group points, ell<=2^(a/4). They may depend on any earlier computations on the seed. There exists a doubled tuple T from I, respecting the stated optional pins, such that every g_s-val(T) is good. Moreover T can be selected deterministically using the ordinary counting oracle, with the same O(mb)-edge size per query as the seed selection.

Proof: the expected nonnegative quantity sum_s Phi_I(g_s-val T) is <=2 ell barPhi_I < tau^2 by (5). Apply conditional expectation to the slot encodings of T. Each term expands using (4), with the group test now val X+val T=g_s or val Y+val T=g_s. This adds only two one-copy tuples and their mask tests; prescribed slots, canonical range conditions and prefix constraints are exact. In terms where a constraint is absent keep the random variables contributing its completion multiplicity. Use exact counts for the valid completions at a prefix. This selects a concrete T satisfying the sum bound, hence every summand is <tau^2. Its incidence vector is recorded. Each good main fiber is then enumerated completely by ordinary prefix counting, costing O(mb) times its number of leaves in counting calls (plus the usual polynomial enumeration overhead). Adding the recorded T to every main tuple gives explicit relations at g_s with three occurrences per retained block and no new unknown columns. The within-fiber contrasts cancel the translation exactly. This is an adaptive construction on ONE arithmetic seed, not reseeding after each new point. QED.

The logarithm of the total conditional bit-time/space bound, including explicit blocks, batch size, enumeration, and polynomial work in the returned tuples, is

    O((h+k+1)b/m + mb/sqrt(L)).                            (7)

Here and below the polynomial overhead in m,b and the 2^O(m) term are absorbed. The same exponent holds for any 2^{O((h+k+1)b/m)} separate adaptive queries (with constants allowed to depend on that bound). For fixed h,k and m=Theta(L^(1/4)), this is the exact saving 1/4 for adaptive robust-fiber access. If h+k+1=Theta(m^theta0), 0<theta0<1, balancing gives component saving (1-theta0)/(2(2-theta0)). These are resource estimates for the stated operation, not the complete relation-selection theorem. In particular allowing *all* possible doubled translations would create a much larger tuple product and does not give the small-product promise of planar avoidance.

## 4. A fixed finite translation family preserving the near-linear tuple product

When an explicit small-product family is needed, it is possible to construct a common list of translations *before* knowing the target. For each retained I separately put f_I=Phi_I/tau^2, and take s=ceil(4b/a)=O(m/(h+k+1)) independent doubled tuples from I. The potential

    Theta = sum_(g in G) product_(j=1..s) f_I(g-val T_j)

has expectation <= r (2 barPhi_I/tau^2)^s <= r 2^(-a s/2)<1 eventually. Indeed a dominates h+log m, and (5) bounds the averages. Fix every translation tuple by exact conditional expectation. Expand every f_I as the finite sum (4); products create at most (4*2^m)^s terms, with O(smb) edges and precision O(smb) per ordinary count, and a shared group-value variable. All normalization factors have polynomial bit length in smb; no new indexed logarithm column is introduced. Empty branches are excluded. One can use the same procedure independently for all I, costing a factor <=2^m.

For each g at least one branch has f_I(g-val T_j)<1. On query evaluate the s potentials exactly, choose a good branch, and enumerate its main fiber. The translated family has at most s W_I tuples (duplicates may be retained as tagged translation choices); its sum covers G and

    s W_I <= r^(1+o(1)) < r^2.

Every recorded translate for I uses only blocks in I; it is not necessary to reintroduce deleted blocks. Treat the shift as a known-incidence choice rather than a new free column. This family also has affordable membership/decomposition by trying each recorded translation. Its exact enumeration guarantee is for a selected good branch, NOT all tuples in a fiber across all branches.

Including this fixed-list preprocessing, the logarithm of the cost is

    O((h+k+1)b/m + m^2 b/((h+k+1)sqrt(L))).                (8)

For fixed h,k choose m=Theta(L^(1/6)), giving component saving 1/6 without a growing-degree mask. For h+k+1=Theta(m^theta0) the balance gives (1-theta0)/(2(3-2theta0)). PR #13 proves globally small, pointwise-flat fibers for any fixed saving below 1/6 with batched higher-degree masks. Our 1/6 estimate has a different guarantee: there is a selected good branch for each target, not pointwise flatness for the entire untranslated tuple family. PR #20 gives full coverage with constant-degree masks at 1/4 but does not furnish this small robust-fiber access on a fixed near-linear-product family. The two distinctions are essential.

## 5. What the robust fibers guarantee for modular rank

For each retained I use reduced incidence coordinates (one reference label per block removed), dimension D_I=sum_(i in I)(n_i-1). A label direction u is in the kernel of all contrasts in a selected fiber exactly when its additive score sum_i u_i(X_i) is constant on that fiber. Translation T does not affect the kernel. Constants per block are quotiented out. Below, score labels may belong to any abelian group where indicated; modular rank uses F_r.

### 5.1 Real covariance and the inherited partial-rank bound

For real labels center within each block under the independent uniform law, and write V_i=E u_i^2, V=sum V_i. One- and two-pin estimates (6) show that the conditional second moment differs from V by at most delta |I| V: off-diagonal contributions are bounded by delta sqrt(V_i V_j). The squared conditional mean is at most delta^2 |I| V. Thus the conditional variance lies between V/2 and 2V. The full set of fiber contrasts therefore has rational rank D_I. As in PR #13, each contrast has squared norm <=2m; Hadamard on a full-rank minor and the prime divisibility of its determinant bound its modular deficiency by floor(D_I log_r(2m)/2). This part is the earlier prescribed-fiber argument now available at an adaptive target; it is not the new modular theorem.

### 5.2 A robust cylinder lemma

Fix a set P of at most h/4 pins. For each remaining block let C_i be an arbitrary allowed subset with excluded fraction e_i<=1/2 and E=sum_i e_i<=h/16. The good fiber contains a tuple satisfying the pins and all the allowed-set restrictions. In fact its conditional probability of satisfying the restrictions is at least (1/2) product_i(1-e_i).

Proof: conditioning on the pins divides two quantities from (6), so all remaining cylinder probabilities of order <=h-|P| are within multiplicative 1 +/- 3delta of the independent law. Let p0=product_i(1-e_i)>=exp(-2E)>=exp(-h/8). Write S_j=sum_(|J|=j) product_(i in J)e_i<=E^j/j!, with S_0=1. Let l be the largest odd integer <=h-|P|. Bonferroni gives a lower bound on the desired probability by the first l terms of inclusion-exclusion. Its independent value is >=p0-S_(l+1). Using l+1>=3h/4, we have S_(l+1)<=(e/12)^(3h/4)<=p0/4 (or zero beyond the number of events). The cylinder perturbation changes the partial sum by at most 3delta sum_j S_j<=3delta exp(E)<=p0/4 with the stated tau. This proves the claim. It is a support result for arbitrary allowed sets; their membership predicates need not be cheaply represented by the oracle because the *whole good fiber* is already explicitly listed. QED.

### 5.3 Uniform modular rank-or-distance theorem

For any label direction u put

    E(u)=sum_i (1-max_c |{A in F_i:u_i(A)=c}|/n_i).

If u is nonconstant within at least one block and E(u)<=h/16, choose modal classes in all blocks. Pin every diffuse block (excluded fraction >1/2) to a modal point. There are at most 2E(u) such blocks. Also pin one nonconstant block j, successively to two labels of different scores. The number of pins is <=2E(u)+1<=h/4 for h>=16. On all unpinned blocks impose their modal classes. Apply the cylinder lemma twice. The resulting two tuples belong to the same fiber, their total scores differ by the nonzero difference of the chosen labels in block j, and every other block has the same score in the two tuples. Thus the direction cannot be in the contrast kernel.

Consequently **every nonconstant direction in the residual modular kernel of a single good fiber has E(u)>h/16**, simultaneously, not merely for a direction selected before the construction. In particular nonconstant directions supported on at most floor(h/16) blocks are excluded. The true logarithm label table, with block constants allowed, belongs to the homogeneous kernel of zero-shift contrasts. Subtracting any affine slope from a surviving nonaffine direction gives another nonconstant kernel direction; the same strict distance lower bound holds against all such slopes. This result works over any abelian score group. Adding more valid constraints preserves the exclusion. It does not imply that no nonaffine directions survive.

Earlier work gave deterministic separators for individually supplied low-modal or low-variance directions (PR #4 residual-kernel notes) and conditional rank-completion branches (PR #11). Here the whole exclusion follows from enumerating one good fiber and holds without finding, classifying or even testing such directions. The effect is a new unconditional structural guarantee on the constructed subsystem, subject to the arithmetic seed/decomposition cost premise, not a full-rank theorem.

### 5.4 Trimmed short lifts: a stronger exclusion with growing h

Assume additionally h>=16+8 ceil(log_2(m+1)). Suppose u is constant in score modulo r on a good fiber. Permit p pinned blocks, and subsets C_i in the unpinned blocks with e_i<=1/2, such that p+E<=h/16. Choose arbitrary integer lifts of u_i on C_i (and of the fixed pinned labels). If the sum of the widths max u_i - min u_i of the chosen lifts on the retained subsets is <r, then u must be blockwise constant on the original blocks.

To prove this, first strengthen the cylinder lemma to a relative estimate for a rectangle under the pins, including after at most two additional fixed slots. With p'<=h/16+2, the largest even integer l<=h-p' satisfies l>=3h/4. The odd/even Bonferroni brackets and (6) bound the absolute error from the independent probability p0' by S_l+3delta exp(E). Here p0'>=exp(-h/8), E<=h/16, S_l<=(e/12)^(3h/4); hence the relative error is <=2^(-h/2). This estimate also applies if some events are absent or all are exhausted in the inclusion-exclusion sum. It guarantees positive probability. Dividing the rectangle probabilities with and without the additional pins, and using (6) conditional on P, shows that the one- and two-slot marginals under the trimmed conditional law are within relative error 8*2^(-h/2) of the independent uniform law on the allowed sets. The extra bound on h makes |I| times this error <1/16. The covariance proof in 5.1 applies to the trimmed law, for any real labels on its allowed sets.

All integer tuple scores for the trimmed law lie in an interval of width <r. Since their residues are identical, they are the same integer. Zero variance and the preceding real covariance bound imply that each score is constant on its allowed subset. The original sum of modal deficits is now at most p+E<=h/16. The rank-or-distance theorem forces original blockwise constancy. QED.

This exclusion is not limited to repeated or modal labels. A direction can have mostly distinct labels and large total circular width but admit narrow intervals after a small trim; such a direction cannot survive these good-fiber contrasts. All possible scalar multiples, affine slopes and lift choices are covered by the statement when their respective trim hypotheses hold. No efficient algorithm for finding arbitrary optimal trims is required for the exclusion, and no universal existence of a suitable trim in the residual kernel is claimed. With h=Theta(log m) and fixed k, (7) gives log-cost O(b sqrt(log L)/L^(1/4)) by choosing m on the scale L^(1/4) sqrt(log L). This yields any fixed saving strictly below 1/4 for this stronger component, not the exact fixed-h quarter-power bound.

## 6. Verification and limitations

Run `python -B check_adaptive_fibers.py` with standard Python and assertions enabled. The verifier rejects optimized mode at startup. The recorded output is in `verification.txt`. The tests cover 31,080 shared-pin energy identities, 720 signed-pattern comparisons, 14 girth controls, 497 mean inequalities and 30 affine bit-prefix identities on the order-seven curve, and 2,574 higher-pin mean inequalities and 96 pattern bounds on the order-thirteen curve. A separate two-slot fixture has two equal sums, 1+29 and 5+25. Its first signed collisions have length four. Changing 29 to 125 gives a collision-free control.

The fixed cover tests include 399 prefix identities, 463 cover checks and 3,250 translated tuples checked by the curve group law. In the five-block fixture each of the five four-block subsets gets its own translation family, covering all 19 targets using only the retained slots. Both selectors' callers check the length and block membership of each returned witness. The fixed-family checks also validate each reported shift, and the adaptive checks verify the prescribed slots. The adaptive tests have 43 successful selections, 91 target requests, 667 exact prefix identities, 209 expanded mean terms and 19 real-curve branches.

Small mask fixtures use relaxed finite goodness thresholds and do not meet the huge sufficient asymptotic guard. In particular the initial toy potential guard fails honestly. An attempted more demanding small h=2 deletion fixture did not pass the relaxed threshold in the tested search, so the real-curve deletion and higher-pin tests are separate controls. The numeric large-b guard checks use conservative bit bounds rather than asserting that a curve of a specified prime order was constructed; one guard fails and four pass. No tensor contractor is implemented or benchmarked. The fixed-list preprocessing estimate and adaptive cost estimates are mathematical conditional bounds, not measured practical performance.

The robust-cylinder fixtures perturb the h-wise-uniform parity distribution at the scale `delta = 8 tau = 1/(2^(h+11) m^2)`. On all 120 rectangles the verifier compares the measured conditional probability with half the independent probability using exact arithmetic. For 120 witness pairs it checks the allowed sets, parity and distinct scores modulo q. Of these pairs, 83 involve diffuse blocks. Negative controls cover a positive mass below the bound, inadmissible words and unequal integer scores with equal residues.

The trimmed controls check 960 exact one/two-pin marginal ratios and 24 variance inequalities. In the distinct-label example the original circular-width sum is 5,388, but a trim of total loss 2 reduces it to 378<1,009 and two attained scores are different. These are synthetic support and distribution checks of the abstract lemmas, not claims about actual elliptic-curve kernel statistics. The retained two-constraint abstract orthogonal array has full (m-1)-wise marginal uniformity yet two independent modular score relations: marginal uniformity in itself cannot prove uniqueness of the modular relation. That example is not claimed to arise from a prime-order elliptic fiber or to obstruct all rank selection.

Two runs of the final verifier reproduce the same transcript byte for byte. Proofs and code remain open for independent review. The universal modular-rank problem remains: eliminate surviving high-distance, wide-lift directions on the actual affine seed, or prove that a charged adaptive choice of good fibers suffices. The new result is neither a solution of general additive avoidance nor a new complete deterministic DLP exponent. The stronger statements about the original arithmetic seed narrow the possible obstructions without assuming a universal affine/tensor representation of arbitrary additive lists.

## Sources and relationship to earlier work

- Liu, Zhang, Ren, ePrint 2026/3020, especially Theorem 2.2; used only as the inherited closed-contraction premise: https://eprint.iacr.org/2026/3020 . The resource arguments rely on the existing arithmetic compilation in the repository rather than a new proof of this premise.
- [PR #13](https://github.com/davidiach/elliptic-relation-counting/pull/13), `research/enumerable_prescribed_fibers/RESEARCH.md` (head bef28551c5fdd352cc3e641cd25e7b3f6ec08eb8): equation (1), prescribed-anchor potential and real-rank bound. We extend to adaptive translated access, high-strength pinned constraints, and modular distance exclusions.
- [PR #18](https://github.com/davidiach/elliptic-relation-counting/pull/18), `research/universal_quotient_compression/RESEARCH.md` (merged head cf5780e7cebd140392aa625d91abaecc6e3dfd7d): quotient reduction and torsion presentation. No new quotient avoidance result is asserted here.
- [PR #20](https://github.com/davidiach/elliptic-relation-counting/pull/20), `research/translation_cover/RESEARCH.md` (merged head 528851643f507e1a92670ae1bcccfcb1f0a560d2): common translation covers; our recorded translations and small good-fiber property are additional.
- [PR #4](https://github.com/davidiach/elliptic-relation-counting/pull/4), residual-kernel, anchored-translation and pointwise-seed notes: existing low-modal separators and higher-moment methods.
- [PR #29](https://github.com/davidiach/elliptic-relation-counting/pull/29), `research/sparse_target_runs/SIX_TARGETS.md`: reduced-floor finite common-kernel controls. They do not establish failure under the new asymptotic seed conditions.
- [PR #39](https://github.com/davidiach/elliptic-relation-counting/pull/39): obstruction for generic affine and tensor representation classes, motivating the return to arithmetic seeds.

Fourier second moments, conditional expectation, and Bonferroni bounds for product/rectangle events are classical; the contribution is their combination with the actual mask family and adaptive enumeration and the simultaneous kernel exclusions. No independent literature-priority or machine-checked proof claim is made.
