# Carry-aware additive probes and the high-rate remote-point problem

27 November 2026. Follow-up to PR #22 (`06b2f51d7c7ba6d2cb67d1f399cf8514ca5b539b`) and the reductions in PRs #18 and #20. We work on the **abstract finite-group additive-cover problem**, not on deployed cryptographic systems. The main new result is a deterministic reduction from the prime-field plane to small-characteristic, high-dimensional additive probes, including an uncolored sparse-span/remote-point formulation with a proved volume gap. A separate exact polynomial method and a degree obstruction delimit that route. No sublinear-in-arity universal probe generator or improved end-to-end uniform deterministic ECDLP exponent is established.

## 1. The target and why a digit lift needs correction

Fix an odd prime r and m nonempty explicit lists A_i in F_r^2, with bounds n_i, N=sum_i n_i, M=prod_i n_i. Assume fixed product slack M<=r^(1+sigma), where 0<=sigma<1. PR #22 shows that a list E which is not contained in any such sumset suffices for its probe formulation; the generator need not identify the absent member. We consider both that input-independent version and an input-dependent list for one collection of A_i.

It is tempting to replace each coordinate by its base-q digits, working over F_q rather than F_r. This is not a group homomorphism, either because of ordinary digit carries or because of reduction modulo r. Ignoring those effects would invalidate the reduction. We give a one-sided, exact overapproximation that preserves the product slack and original arity scale.

## 2. Carry-aware digit-lift theorem

Let q>=2 be prime with q<=r, and k=floor(log_q r)>=1. Define phi(x) to be the k least significant base-q digits of the canonical integer representative 0<=x<r. Use the same low-digit extraction on the integer r itself, and put d=2k. Define Phi(x,y)=(phi(x),phi(y)) in F_q^d. The output decoding map iota maps every vector in F_q^d to the pair of nonnegative integers represented by its two k-digit words. Both decoded coordinates are below q^k<=r. Phi need not be injective on the whole original plane, but Phi(iota(z))=z for every z in the digit space.

Let C be the following correction sumset in F_q^d:

    R = {(-t phi(r),-u phi(r)): 0<=t,u<m},
    C0 = {(c_0,...,c_(k-1),e_0,...,e_(k-1)):
           c_0=e_0=0, all other digits in [-m,m] modulo q},
    C = R+C0.

Then the exact inclusion

    Phi(A_1+...+A_m) subset Phi(A_1)+...+Phi(A_m)+C          (1)

holds. In particular, if E' is not contained in the expanded sumset on the right, then iota(E') is not contained in the original sumset. An absent digit vector always decodes to an absent plane point. No inverse assertion is made: the enlarged sumset can contain spurious digit vectors.

**Proof.** For one coordinate write x_i for the canonical input integers and z=(sum_i x_i) mod r, so sum_i x_i=z+t r, with 0<=t<m. Let s_j=sum_i phi(x_i)_j-t phi(r)_j-phi(z)_j for 0<=j<k. We have sum_j q^j s_j=0 mod q^k and |s_j|<=m(q-1). Define c_0=0 and c_(j+1)=(s_j+c_j)/q. Integrality follows inductively by reducing the remaining low-digit congruence modulo q and dividing by q. The bound |c_j|<=m follows by induction, since |s_j+c_j|<=mq. Reducing the recurrence modulo q gives

    phi(z)_j = sum_i phi(x_i)_j - t phi(r)_j + c_j  (mod q).

The second coordinate has its own independent wrap count and carry sequence. Both belong to the advertised correction set. This proves (1), including original sums with one or both canonical coordinates >=q^k. For a decoded candidate, Phi(iota(z))=z, so the stated implication follows by contradiction. QED.

The correction need not be materialized as one exponentially large list. Partition the 2(k-1) carry coordinates into chunks of at most ell coordinates. For each chunk form one list supported there, allowing all values in [-m,m] for its coordinates and zero elsewhere. The sum of these lists is precisely C0. Together with R this is an exact factorization of the correction box, with at most m^2 choices in R and at most (2m+1)^ell choices in each chunk. It makes no assumption about actual correlations between carry states.

### Preserving the subexponential parameter budget

Let b=ceil(log_2 r), let m grow polylogarithmically in b with m log m=o(b), and take a fixed positive integer c>8/(1-sigma). By Bertrand's postulate choose a prime q in [(2m+1)^c,2(2m+1)^c]. Since c is fixed, even trial division through that interval costs polynomial time in m (not in r). Eventually q<r and q>2m+1. Choose the largest ell satisfying

    (2m+1)^ell <= 2^ceil(b/m).

This is an exact integer powering/comparison, not an uncharged logarithmic approximation. For sufficiently large inputs ell>=1. Let j=ceil(2(k-1)/ell) be the number of carry chunks. The lifted list count m'=m+1+j satisfies

    k=(1+o(1)) b/[c log_2(2m+1)],
    ell=(1+o(1)) b/[m log_2(2m+1)],
    m'=(1+2/c+o(1))m.

There are only O(m) extra lists, each with at most 2^ceil(b/m) entries, and the explicit correction entries have d log_2 q=O(b) bits. If M>=r, then N>=r^(1/m), so the correction sizes are O(N), and total explicit bit length is polynomial in N,b. More generally the time bound is polynomial in N+2^ceil(b/m),b; in the regime used by the prior coverage reduction these are equivalent on the logarithmic scale. No need to enumerate all corrections jointly.

Let M' be the product of the lifted size bounds. It satisfies

    M' <= M m^2 (2m+1)^(2(k-1))
       <= r^(1+sigma+2/c+o(1)).

The ambient cardinality is q^d=q^(2k)>r^2/q^2=r^(2-o(1)). Since c>8/(1-sigma), the gap is bounded away from zero in the exponent. In particular for all sufficiently large parameters M'<q^d, with fixed slack M'<= (q^d)^(1-epsilon) for some epsilon>0 depending on sigma,c. The field size q is polynomial in m, while d=Theta(b/log m) grows. The larger field-size bound q=poly(b) would permit materializing C in a single list at the cost 2^{O(b log m/log b)}; the chunk factorization is stronger because it makes q=poly(m) possible while keeping O(m) lists and the same explicit-size budget.

A universal list generator for these high-dimensional small-field inputs therefore supplies a universal plane probe list by iota. This preserves its uniform running time/output-size exponent N^{O(m^theta)} for any fixed theta>=0, up to polynomial factors in b,N, and fixed constants in the O. For an input-dependent generator the reduction is equally valid. It introduces no random choice or nonuniform advice. It neither claims that a small-field generator is known nor requires proving any nonmembership within the original sumset.

Finite parameter caveat: all comparisons should be guarded in an implementation. A conservative uncolored sparse-span bound below fails for one of the tested small-arity large-bitlength combinations due to integer rounding. The asymptotic theorem assumes growing m; it is not a claim that the asymptotic constants automatically certify each finite parameter choice.

## 3. A high-rate remote-point problem with actual volume slack

Previously the planar formulation could not be strengthened to the usual sparse-linear-span avoidance: in dimension two, two independent columns already span the plane. After the digit lift there is a different, stronger total-search problem that really has solutions under the promised volume bound.

Use the balanced bound relevant to PR #20's coverage setup:

    max_i n_i <= r^((1+sigma+o(1))/m),   M>=r,

allowing a bounded number of extra translating lists (absorb them into m and the o(1)). Keep repeated columns if necessary, and form a dictionary matrix H of dimension d by K from every lifted list including the correction lists. Then K<=N+O(m)*2^ceil(b/m)+m^2 and

    log_2 K <= (1+sigma+o(1)) b/m.

Put w=m'. Every expanded tuple sum is a linear combination of at most w dictionary columns, with coefficients one (or repeated columns can be combined). The set of ALL linear combinations of at most w columns, now permitting arbitrary coefficients in F_q, has size at most (1+qK)^w by encoding w ordered choices of either zero or a column and coefficient. Since log q=O(log m),

    log_2 (1+qK)^w <= [(1+2/c)(1+sigma)+o(1)] b < (2-o(1)) b.

The strict inequality is ensured by c>8/(1-sigma). Hence this stronger sparse span also fails to cover F_q^d. If rank(H)<d, finding a vector outside its whole column span is elementary. If rank(H)=d, the code L=ker H has length K and dimension K-d. By the classical syndrome characterization, a syndrome outside all w-sparse combinations of columns is equivalent to a coset of L whose Hamming distance from L exceeds w. Keeping repeated columns is legal and makes K grow much faster than d in the stated regime. Thus the code rate tends to one. A list containing at least one such syndrome suffices; an absence certificate for the syndrome need not be produced in PR #22's probe architecture.

This establishes a new route to an ordinary high-rate remote-point / remote-set problem on explicit matrices over q=poly(m), with d=Theta(b/log m), target remoteness w=O(m), and a certified constant exponent gap in the sphere-covering count. In particular no arbitrary-coefficient span promise that is false in the original plane has been smuggled in. For arbitrary unbalanced lists the colored variant (at most one dictionary element per list with arbitrary coefficients) still has volume <=q^m' M'<q^d for large inputs; the uncolored assertion uses the balanced bound above.

This is a reduction, not an application of an already sufficient remote-point theorem. Ben Lee Volk, *Improved Algorithms for the Remote Point Problem*, ECCC TR26-203 (2026), Theorem 1.2 (and its extension for a fixed rate bounded away from one) does not cover these high-rate codes. Its rational-field construction does not transfer to these finite fields. The syndrome equivalence itself is classical (see that paper's Lemma 2.1). In the earlier plane formulation the sparse span could already be the entire plane, so even existence of a remote point was unavailable; (1) repairs that precise obstruction. Finding such a syndrome or remote set in the desired uniform N^{o(m)} time remains open here. The generic theorem has no automatic stronger rate-dependent implication merely by inserting K-d for k in a formula outside its hypotheses.

## 4. A complete exact test for bounded-degree algebraic separators

Consider arbitrary nonempty lists B_i in F_q^d (including the lifted lists if desired), choose references b_i, and write s0=sum_i b_i, D_i=B_i-b_i. For t>=0 let S_t be the sums obtainable by choosing nonreference deviations in at most t blocks, with all other blocks fixed at their references. S_t is explicitly enumerable using at most

    sum_{I subset [m'], |I|<=t} prod_{i in I} |B_i|

choices. For every polynomial P in d variables of total degree <=t, its restriction to the full sumset is determined by its values on S_t. In particular

    P vanishes on sum_i B_i iff P vanishes on S_t.           (2)

**Proof.** For an arbitrary tuple of deviations v_i define Delta_v P(x)=P(x+v)-P(x), which lowers total degree by at least one unless it is already zero. The translation operators commute, so P(s0+sum_i v_i)=sum_I (prod_{i in I} Delta_(v_i))P(s0). All terms of order >t vanish. Every surviving term is an alternating sum of evaluations at partial deviation sums with at most t active blocks. If P is zero on S_t all terms are zero. Necessity follows from S_t being a subset. QED. This is an exact finite-difference identity over any field and does not assume independent samples, positive measures, or noncancellation of a modular Gram matrix.

Use reduced monomials with each exponent <q, total degree <=t, at most binom(d+t,t) of them. Build the evaluation matrix on S_t and solve for a nonzero vector in its kernel. A reduced nonzero polynomial over F_q is a nonzero function on the full field grid; uniqueness follows by repeated one-variable interpolation. A point where it is nonzero is constructed by fixing variables one at a time, keeping a nonzero coefficient polynomial at each step. At most q trials per variable suffice, or test only degree+1 distinct values. By (2) that point lies outside the full sumset. If the evaluation matrix is full column rank, report that no separator exists in this degree class, not that the full sumset equals the ambient space. The procedure takes polynomial work in the number of monomials, the explicit S_t bound, d, log q (and q if scanning all residues), hence (N'+d+t)^{O(t+1)} poly(log q,m') with degree-bounded trial sets and explicit list sizes N'. The count of monomials is charged, not replaced by an ideal-membership oracle.

This is an actual finite deterministic avoidance algorithm for the stated class. For instance products of the bilinear graphs {(x,y,xy): x,y in F_q} on disjoint coordinate triples have full affine span and **every nonzero scalar projection is surjective**, but a quadratic separator exists. Each nonzero linear form on one graph is surjective: if its xy coefficient is nonzero, fix x so the coefficient of y is nonzero; otherwise use a nonzero linear coefficient. The graph has no affine-linear annihilator. The polynomial z-xy vanishes exactly on the graph in that triple. Arbitrarily many disjoint triples are allowed, so this is not an implication of small scalar support. This example concerns the abstract digit-space problem, not a claim that the carry-expanded lists of arbitrary planar inputs preserve that quadratic graph.

### Two obstructions to using (2) as a universal breakthrough

**Hamming-ball obstruction.** In F_q^d let B(t0) be all vectors of Hamming weight at most t0. The sum of m lists B(t0) is exactly B(h), h=min(d,m t0), since a vector of weight <=m t0 can partition its support among the lists. For h<d, the reduced polynomial functions with every monomial supported on at most h coordinates restrict injectively to B(h): for each coordinate set J of size <=h, restrict outside coordinates to zero; vanishing on the full F_q^J forces all coefficients with support contained in J to be zero. Thus the minimum total degree of a nonzero reduced polynomial vanishing on B(h) is exactly h+1 (the product of h+1 distinct coordinate variables attains it). This is stronger than a dimension-count heuristic and also excludes high exponents on only h coordinates.

In the lifted asymptotic regime q=poly(m), d=2 floor(log_q r), choose t0=floor(d/(8m)). For growing m and b/(m log m) tending to infinity, t0>=d/(16m) eventually. The ball size is <=(t0+1)(e d(q-1)/t0)^t0, so log_2 |B(t0)|<=log_2(d+1)+(d/(8m))log_2(16e m(q-1)) < b/m eventually for q=(2m+1)^c up to factor two and c>8. The original list-size caps of at least r^(1/m) can accommodate these balls; pad to the stated balanced caps if desired. Since the digit decode of every vector is a legal plane point, one can realize these as the digit lifts of actual plane lists, and the carry correction includes zero, so the expanded sumset contains B(h). Their product satisfies the earlier volume guard with padding chosen within the fixed slack. Thus a universal low-degree annihilator guarantee in the lifted class cannot hold below Omega(d)=Omega(b/log m), much larger than the desired o(m) degree. This example is not asserted to be a projected factor base arising from the prescribed elliptic extractor, nor a lower bound for all algorithms.

**Carry-box obstruction (independent of Hamming balls).** In the digit lift with q>2m+1 and k>=2, separate the two least significant coordinates from the d-2 others. If the pre-box sumset sum_i Phi(A_i)+R projects surjectively onto those two coordinates, then the final expanded sumset contains, above each low pair, a translate of a full (2m+1)^(d-2)-point Cartesian grid in the other coordinates. Any reduced polynomial P with degree at most 2m in each of the latter coordinates which vanishes there is identically zero: fix a low pair, interpolate on its grid, then interpolate the coefficient polynomials on all q^2 low pairs. In particular no nonzero reduced total-degree <=2m separator exists. Surjectivity of the low pair is already attained by a single list containing the q^2 points with both coordinates in [0,q); this fits the huge asymptotic caps and can preserve the fixed volume slack. The examples in the verifier have product strictly below ambient even after the correction. The low-degree barrier holds for this intentionally enlarged support, not necessarily for the original plane support. Reducing the number of original slots by first enumerating groups of s lists only changes the naive carry width to O(m/s); a degree-by-pin method of that order on grouped lists of size N^s still has an N^{O(m)} budget, not a proved sublinear exponent. Exploiting correlated, rather than independent, carries could avoid this particular obstruction, but no affordable general construction is given.

## 5. Verification, literature and exact status

Run `python -B check_digit_lift.py` (standard library, assertions enabled). `verification.txt` records the completed run. The verifier checks 64,645 exhaustive scalar carry identities and 5,000 seeded larger cases; four direct-vs-factored correction sets, twelve full synthetic plane fixtures, 147 original-support inclusions, and **all 581,359 absent digit vectors** in those fixtures for safe decoding. It checks 192 exact finite-difference evaluation-span identities, 186 with strictly smaller restricted sets; four Hamming-ball interpolation ranks; thirty toy missing-syndrome checks by exhaustive sparse search; three actual quadratic selector successes on graph products (with every scalar projection surjective); and two carry-box negative controls. The sparse-syndrome experiments do not implement a fast high-rate remote-point algorithm. Large-parameter checks use conservative integer bit bounds, not actual huge primes or block enumeration. One initial parameter guard failed on the conservative uncolored volume estimate at m=4; it is deliberately retained as a negative control rather than presented as a pass. Four growing-arity large-bitlength controls pass both volume inequalities. Finite checks illustrate the theorems, not their asymptotic proof or practical cryptographic performance.

Literature inspected: Ben Lee Volk, *Improved Algorithms for the Remote Point Problem*, ECCC TR26-203, https://eccc.weizmann.ac.il/report/2026/203/ (Theorem 1.2, Lemma 2.1); Alon–Panigrahy–Yekhanin, *Deterministic approximation algorithms for the nearest codeword problem*, https://doi.org/10.4230/DagSemProc.09421.4. The connection of sparse factorization to the stronger arbitrary-coefficient "sumset-evasive" problem occurs in the literature on linear data structures and rigidity; see https://eccc.weizmann.ac.il/report/2018/188/revision/2/ . Those results concern a stronger notion than one-hot additive containment; none was used as a solver. Sumset recognition hardness (Abboud–Fischer–Safier–Wallheimer, https://arxiv.org/abs/2410.18661) and sumset covering literature do not directly give a lower bound for the promised avoidance task. The elementary carry, sparse volume, and finite-difference ingredients are classical in isolation; priority of this combined reduction is not established by this preliminary search.

**Status:** the new digit-lift theorem, preserving arity, field-size and volume slack, is unconditional. So are (2) and the degree obstructions. An affordable high-rate remote-set construction or a direct small-field universal probe generator remains a missing lemma. Through PR #22 and PR #20, a successful generator at N^{O(m^theta)} time for fixed theta<1 would retain their conditional exponent implication delta=(1-theta)/(2(2-theta)); no such runtime is established here. The upstream contractor remains an assumption; no independent proof-chain or novelty review is claimed. The mandatory end-to-end breakthrough checklist fails. No real key recovery, curve benchmarks, or deployment experiments are included.
