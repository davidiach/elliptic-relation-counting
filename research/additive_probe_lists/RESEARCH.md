# Universal additive probe lists and barriers to algebraic probe generators

Research follow-up, 27 November 2026, to PR #18 (`9d3ae7ab3e9253391c392bf0a4a06e1e935a2631`), PR #20 (`ea22af035570b6577c02423bd8186617ed134144`), and the merged compact-cover results in #15. This note works on the abstract additive avoidance problem itself, not a deployed cryptographic system. The strongest end-to-end deterministic ECDLP objective remains **unproved**. The constructive reduction to a universal probe list below inherits the closed LZR Theorem 2.2 and the short-coverage construction in #20. The finite-plane results and the generic avoidance baseline do not use the contraction premise.

## 1. We need a probe list, not a certified absent point

Fix an odd prime r, nonempty indexed sets A_i in F_r^2 of sizes at most n_i, i=1,...,m, and write

    M = prod_i n_i < r^2,  N = sum_i n_i,  D_* = 1 + sum_i(n_i-1).

A universal probe list E for the size vector (n_i) means a finite set E subset F_r^2 such that for every such collection A_i,

    E is not a subset of A_1+...+A_m.

E depends on r and the size bounds, but not on the contents of the A_i. A weaker, instance-dependent version suffices too: construct E from the lists, with a proof that E is not contained in their sumset. The task of identifying which member is absent is **not required** in the application of PR #18. This is a distinct weaker deliverable than outputting one certified missing point of AdditiveAvoid(2).

To see this, use the notation of PR #18: V=F_r^D, pi:V->G the public group-valued map, W the row space of the known valid equations including a calibrated row of value P, U=W intersect ker(pi), H=V/U of dimension d+1, where d is the residual modular dimension. If d>0 choose a surjective linear T:H->F_r^2. The projected tuple lists A_i=T(C_i) are explicit. Their size bounds can be kept at the original n_i even when projected points coincide. For every z in the probe list, lift it to some x in V with T(x+U)=z, compute the point g_0+pi(x), and use the ordinary short-tuple decomposition interface at this point to obtain X. Compute

    v=x-c_X.

Check pi(v)=0 by group arithmetic, and T(v+U) by explicit elimination. If T(v+U) is nonzero, the verified row v is not in U, and hence not in W. Append it; the residual dimension falls by one. If z is absent from the projected tuple support, **every** valid decomposition at that point gives nonzero T(v+U). Some present z may also give progress because the supplied decomposition is not one realizing z. The procedure therefore needs neither a coNP nonmembership verifier nor exact additive-sumset membership queries. A universal E can be reused in every quotient round, up to D rank increases. If no probe gives progress, then every tried z is in the sumset, contrary to universality. There is no assumption about the canonical decomposition's behavior other than validity and totality.

Using #20's short-coverage/translation seed, the bit-cost consequence is as follows. With b the curve field bit length, L=log_2 b, and m the main block count, the explicit list sizes have log N=O(b/m+log b), the total product is <r^2, and ordinary decomposition has the inherited O(mb/sqrt(L)) log-time cost, up to polynomial overhead. If a uniform deterministic generator for the requisite probe list has time and output size N^{O(m^theta)} poly(b) for any fixed 0<=theta<1, the total log bit-time bound including testing every probe, at most D rounds, and all linear algebra is

    O(b/m + m b/sqrt(L) + m^theta b/m).

Balancing gives the same conditional implication delta=(1-theta)/(2(2-theta)) as in #20. In particular an explicit polynomial-time universal probe generator would imply the quarter-power objective. The point is to remove the need to certify which vector lies outside the sumset. The generator itself is still a load-bearing unresolved obligation. Random or nonuniform existence is not a uniform construction.

## 2. The optimal scale of a universal probe list

**Theorem.** Let H(r;n_1,...,n_m) denote the least size of a universal probe list for these size bounds. Then

    D_* + 1 <= H(r;n_1,...,n_m)
              <= min(r^2, floor(2 D_* log r / log(r^2/M)) + 1).

Here log may be in any fixed base. In particular, for any fixed sigma<1, if M<=r^(1+sigma), then

    D_*+1 <= H <= floor(2 D_*/(1-sigma))+1.

The bounds hold for arbitrary distinct input sets, not just random or low-doubling sets. In the balanced growing-arity regime D_*=Theta(N); thus the existential list size is Theta(N). The upper bound is not a constructive time bound.

Lower bound: let E be a nonempty set of at most D_* points and fix e in E. Distribute its remaining |E|-1 point differences z-e among the m blocks, with at most n_i-1 differences assigned to block i. Form each block from its assigned differences and zero, and translate the first block by e. Every z in E now has a representation using its difference in the appropriate block and the zero choices in the others; e uses all zero choices before translation. Enlarge sets if exact sizes n_i are desired, without exceeding the size bounds. The empty E is trivially not universal. Thus every E of size <=D_* can be covered. The product promise guarantees all n_i<=r^2, so the padding causes no obstruction.

Upper bound: normalize each A_i by subtracting a chosen reference a_i. A sumset is described by the single aggregate shift sum_i a_i and at most sum_i(n_i-1) other field-plane points, padding by repeated zeros when necessary. Hence there are at most r^(2D_*) distinct possible sumsets (a count of ordered descriptions, which may overcount). Each contains at most M points. For K independent uniform probes the probability that some legal sumset contains every probe is at most

    r^(2D_*) (M/r^2)^K.

If K>2 D_* log r/log(r^2/M), the probability is strictly less than one. A successful sequence has a universal underlying set of at most K distinct points. The whole plane is also universal. This proves the theorem. No sample-based absence certificate, uniform generator or novel general union-bound method is asserted.

This distinguishes the additive structure from the arbitrary direction-oblivious fiber model of #14. That model requires Omega(sqrt(r)) probes even with one value per line, whereas actual balanced additive sumsets have universal lists of size O(N)=r^{O(1/m)} under fixed product slack. There is no contradiction: arbitrary thin-fiber supports can contain many sets that are not expressible as these small additive sumsets.

### An exact incidence-matrix formulation of the construction problem

Normalize the lists as above and allow repeated entries. Let X in F_r^(D_* x 2) hold the total shift and all nonreference point differences. Let R be the M possible row vectors with leading coordinate one and a one-hot choice (or zero for the reference) in each of the m blocks. A sequence Y of K plane points is coverable if and only if there are a K x D_* matrix B with each row in R and some X with

    B X = Y.

The two columns of Y must lie in col(B) simultaneously. Therefore constructing a universal list is exactly the task of constructing Y that is not in col(B) x col(B) for every such one-hot incidence matrix B. This is a precise restricted subspace-avoidance problem; it is not generic avoidance of all D_*-dimensional subspaces (which would be impossible for a two-column matrix when D_*>=2).

A useful exact potential is the integer number of normalized list descriptions **with chosen witnesses** for every probe:

    Psi(Y) = sum_{B in R^K} r^(2(D_*-rank B)) * 1{both columns of Y in col(B)}.

A coverable Y has Psi>=1, and a noncoverable Y has Psi=0. This follows by solving B X=Y columnwise; repeated list entries give legal smaller sets. For independent uniform Y the expectation is M^K r^(2D_*) r^(-2K). It supplies the same size bound and an exact conditional-expectation formulation. Direct evaluation still enumerates M^K witness assignments, or their huge equality-partition patterns; we have no small contraction for this potential. Its nonnegativity does not by itself make the conditional expectations affordable. This formulation identifies the specific subspace family on which a useful explicit generator must be proved, rather than invoking generic circuit range avoidance as a black box.

## 3. An interpolation obstruction to low-degree algebraic probes

The following is a separate constructive obstruction. It shows that checking a long segment of a simple algebraic curve (or a small collection of such curves) can fail even though a much smaller universal list exists. The construction uses positive-size ordinary unweighted sets and does not require the sets to contain all points of a chosen curve as individual elements.

**Interpolation-cover theorem.** Let F_1,...,F_ell be polynomial maps F_r^s -> F_r^2 of total degree at most d, where d>=1 and s>=1. Fix integers q>=2 and k>=1. Put L_0=d(k-1)+1 and assume L_0<=m and L_0<=r. For j=0,...,L_0-1 set

    R_j=(q-1) sum_{a=0}^{k-1} j^a

(with j^0=1). Then the union of all probes

    F_h(t_1,...,t_s),  1<=h<=ell,  0<=t_i<q^k,

where integers are reduced modulo r for evaluation, is contained in an m-fold sumset with the jth nontrivial list of size at most ell(R_j+1)^s. In particular a uniform cap n is sufficient if

    ell [1+(q-1) k max(1,L_0-1)^(k-1)]^s <= n.

With n^m<r^2 this is a legal family of size bounds for AdditiveAvoid(2). Pad the remaining lists with singletons, or enlarge to exact sizes when requested.

Proof. Expand each t_i in base q, t_i=sum_{a=0}^{k-1} t_(i,a) q^a, 0<=t_(i,a)<q, and let T_i(Z)=sum_a t_(i,a) Z^a. The polynomial vector F_h(T_1(Z),...,T_s(Z)) has degree at most d(k-1). Lagrange interpolation at the distinct nodes 0,...,L_0-1 over F_r gives public weights

    omega_j = prod_{v!=j} (q-v)/(j-v) in F_r

such that its value at q is the sum over j of omega_j times its value at j. Each integer T_i(j) belongs to [0,R_j], independently of the other parameters' digit choices. Set

    A_j = {omega_j F_h(x_1,...,x_s) : h in [ell], x_i in [0,R_j]}.

Every probe is the sum of the indicated one point from each A_j, using the same map h and consistent digit polynomials. The sets may contain many additional mixed choices; that is allowed for containment. Their sizes obey the bound, even if some omega_j=0 or some evaluations repeat. This establishes the claim. The proof uses classical interpolation; its application to the restricted probe family is the new observation. QED.

**Quantitative consequence.** Take the balanced regime n^m<r^2, log n >> m^2 log m, m growing. For fixed degree d and fixed parameter count s, choose k=floor((m-1)/d)+1 and ell<=n^eta for any fixed 0<=eta<1. Choose q as the integer floor of (n/ell)^(1/s)/(2 k max(1,m)^(k-1)), which is >=2 eventually and satisfies the size guard. Then

    log(q^(ks)) >= k(log n-log ell) - O(s k^2 log m)
                 >= (1-eta-o(1)) k log n.

For ell=1 and d fixed, a single such map on a parameter box with up to n^((1-o(1))m/d) grid points can be completely covered by the legal sumset. For the parabola (t,t^2), this includes all consecutive-parameter probes of length n^{o(m)}, and extends to length n^{(1-o(1))m/2}, subject to the stated field and size conditions. In the near-linear-product regime n^m=r^{1+o(1)} this is r^{1/2-o(1)}. In the fixed-slack regime n^m=r^{1+sigma+o(1)}, it is r^{(1+sigma)/2-o(1)}. These are parameter ranges for the displayed construction, not a claim that every subset of that many arbitrary plane points is coverable.

More generally, if an algebraic universal probe generator is based on one polynomial map of degree d, evaluated on a consecutive interval of length n^{C m^theta}, theta<1, then the preceding construction excludes all such generators for which k(log q)/log n >= C m^theta (in particular fixed d, and sufficiently small d of order m^{1-theta}, with the relevant constants). This is a restriction on that particular generator class, not a lower bound on the degree of an arbitrary algorithm. It also applies to small unions of curves and to small integer parameter boxes in several variables. An arbitrary set of parameter residues not contained in such a box after a public change of variables need not fall under this theorem. High-degree maps, rational maps, irregular point sets, and input-dependent probe selection are not ruled out.

### Even whole affine lines may be covered

A simpler and stronger result for degree one covers entire lines. Let q be the least integer with q^m>=r, and let ell q<=n. For any ell affine lines a_h+F_r v_h in the plane, let

    A_i = union_{h=1}^ell {1{i=0}a_h + j q^i v_h : 0<=j<q},
         i=0,...,m-1.

Any t in F_r has a representative in [0,r) with m base-q digits. Thus all points a_h+t v_h have representations in the sumset, and each list has size <=ell q. If n^m<r^2, this size estimate guarantees coverage of any family of at most floor(n/q) fixed lines, with q=ceil(r^(1/m)). At n^m=r^(1+sigma) the threshold before rounding is r^(sigma/m).

For comparison, section 2 gives a universal set of O(D_*) points under fixed product slack. Taking one line through each point gives a universal union of O(D_*) lines, hence O(N) lines when N is the sum of the size caps. This is an existence result. The efficient uniform construction of a probe list remains open.

## 4. Why short dependency arguments for universal probes are inadequate

A natural approach to the incidence formulation of section 2 is to choose Y with no short additive dependencies, then show that every large collection of one-hot rows has a short linear dependence. The latter claim is false in a very strong sense.

**Theorem (high-girth witness incidence).** Suppose all m blocks have n possible labels, n>=2. For integers K>=2 and s_0>=2 satisfying

    sum_{s=2}^{s_0} K^s (s/sqrt(n))^(m s) < 1,

there exist K distinct tuple-index rows such that every nonempty subset of at most s_0 rows has a singleton label in at least one block. Consequently every subset of at most s_0 of their full one-hot incidence rows is linearly independent over every field. The same holds after deleting the reference columns and adjoining the leading constant one, as in section 2.

Proof. Draw K independent uniform index tuples. For a fixed subset of s rows, if a block has no singleton among their labels, its equality partition has at most s/2 parts. There are at most s^s partitions and at most n^(s/2) assignments of distinct labels to those parts. The probability of no singleton in that block is thus at most (s/sqrt(n))^s. Different blocks are independent. Union over at most K^s choices of rows of size s and over s=2,...,s_0. The displayed guard makes the total probability strictly less than one. A successful array has no repeated tuple because a repetition would fail the s=2 singleton test. In a linear dependence with nonzero coefficients on its support, a singleton coordinate has just one such coefficient and cannot vanish. Deleting references and retaining the constant column preserves this argument: vanishing of all nonreference sums and of the constant sum implies vanishing of every omitted reference coordinate as well. QED.

In particular K<=n^(m/8), s_0<=floor(n^(1/4)), and n^(m/8)>=2 give a sum bounded by sum_{s>=2} 2^-s=1/2. For growing m this can make K much larger than the total number of columns, while no dependency exists on even n^(1/4) rows. All rows cannot be independent when K>D_*, so their linear dependencies necessarily have large support. A universal probe family of length n^{o(m)} may therefore face witness-incidence patterns with no local dependency at any of these subexponential support scales. This is not a construction of covering list values for a specified Y: the eventual large dependencies can still prevent that. It rules out a proof that tries to force a short linear dependency solely from the tuple-index incidence combinatorics.

There are also explicit elementary fixtures of this phenomenon. Over a small prime field F_q with q>=m, take all evaluations at m distinct points of polynomials of degree <=d<m. The K=q^(d+1) tuple-index rows pairwise agree in at most d slots. In any subset of s rows, if a fixed row's label were never a singleton, it would need at least m agreements with the other s-1 rows, whereas their total number is at most d(s-1). Thus every subset with d(s-1)<m has a singleton, independently of the characteristic used for the final incidence rank. The probabilistic bound above is stronger asymptotically; the explicit polynomial code is useful for finite tests.

## 5. An exact direct avoidance algorithm and charged baseline

The main constructive result on AdditiveAvoid(2) without extra promises in this pass is a meet-in-the-middle algorithm that returns an actual absent point without enumerating the full sumset or scanning all r^2 output points. This is a baseline improvement over full product enumeration, **not** a sublinear-in-arity exponent. It uses classical split-and-list and orthogonal range counting.

Let P=prod_i |A_i| be the number of input tuples, so P<=M=prod_i n_i. Partition the lists into two groups with products P_1,P_2 of their actual lengths, assigning each next list to the group with smaller current product. Then P_1 P_2=P. The ratio of the larger to the smaller product stays at most n_max=max_i n_i. Multiplying the smaller side by a list length reduces the old imbalance if it remains smaller, and otherwise the new ratio is at most that list length. Hence max(P_1,P_2)<=sqrt(P n_max)<=sqrt(M n_max).

Enumerate the two half-sum arrays A,B, retaining multiplicities. For a rectangle R=[a,b)x[c,d) in canonical plane coordinates 0,...,r-1, define C(R) as the number of pairs (x,y) in A x B with x+y modulo r in R. For each x in A, shift R by -x modulo r. It is a disjoint union of at most four ordinary axis-aligned rectangles. Counting y in B for each shifted rectangle is a batch of O(|A|) standard 2D orthogonal range counting queries on B. Include multiplicities throughout. The counts are exact nonnegative integers. Reduce rectangle queries to four dominance queries using inclusion-exclusion and compute all of them in O((|A|+|B|) polylog(|A|+|B|) poly(log r)) bit time by sorting and an offline sweep (or the standard divide-and-conquer dominance algorithm). On sequential tapes, a recursion by x with merge-sort and y-ordered cross-dominance sweeps uses polylogarithmically many passes; the boundary tie order is strict for half-open intervals. All arithmetic is on O(log M+log r)-bit integers. The finite implementation uses a coordinate-compressed Fenwick sweep on the word-RAM, not an assertion that arbitrary random access is free on a multitape machine.

The root rectangle contains every tuple, so its exact count is C(R)=P_1 P_2=P<=M<r^2=area(R). Bisect one non-singleton coordinate interval. Evaluate C on one child and subtract that count from the exact parent count to get the other. At least one child has count strictly below its area because the two child areas add to the parent area. Continue on that child. After O(log r) splits both intervals have length one. The count is a nonnegative integer strictly below one, hence zero. Return that plane point. A point has a representation if and only if its count is positive.

The total time and space bound on the explicit input is

    (P_1+P_2+N) poly(m,log r,log M)
      <= (sqrt(M n_max)+N) poly(m,log r,log M).

For balanced n_i this is n^{ceil(m/2)} up to polynomial factors. Because global group coverage in the original application forces M>=r, this is still exponential in b at the relevant parameters, not the desired 2^{o(b)} result. The proof does not hide an NP oracle or a scan of r^2 outputs. It provides a reproducible exact benchmark against which faster avoidance candidates can be tested.

**Optional output-sensitive variant.** If the two half-sum supports themselves are small, replace each enumerated half by its distinct support; the pair count of these supports still certifies absence and its total is <=M. Published deterministic sparse nonnegative convolution algorithms compute two integer-set sumsets in time near-linear in output size on the specified word-RAM (Bringmann, Fischer, Nakos, SODA 2022). Encode (x,y) as x+2r*y; adding two canonical points produces no cross-coordinate carry, so decoding and modular reduction yield the desired cyclic plane sumset. Each reduced sum has at most four unreduced lifts. As each intermediate nonempty partial support has size no greater than the final half support, the iterative output-sensitive construction is near-linear in m times the half-support sizes, up to polynomial factors in log r. This is a useful word-RAM bound when there is substantial intermediate collapse, not a universal small-support theorem or a new sparse convolution algorithm. The finite verifier uses literal partial sums, not an implementation of the cited fast convolution algorithm.

## 6. Research assessment and limitations

Run `python -B check_avoidance.py` from this directory (standard library; assertions enabled). The recorded transcript includes exact modular rectangle counts, complete finite avoidance trials, exhaustive universal-list minima on four tiny planes, synthetic quotient completion with universal probe reuse, polynomial interpolation coverage, whole-line coverage, and incidence-code checks. The larger polynomial interpolation check uses arithmetic identities and explicit bounds without listing its enormous lists. All such instances are mathematical controls and do not implement the elliptic tensor contractor. The asymptotic bounds are established by the preceding proofs, not inferred from the finite experiments.

The strongest new structural result is the tight order of the universal probe number together with the explicit interpolation cover showing that all short consecutive-parameter probes on bounded-degree algebraic maps can be hidden in valid additive lists. The useful new reduction is that the rank increment can be certified by ordinary decomposition without certifying any probe's nonmembership. This focuses the remaining task on *explicit* universal additive hitting lists (or an input-dependent list generator) rather than the stronger total-search problem of outputting one certified missing vector. We do not prove a sublinear-in-arity construction time for either. The incidence-matrix formulation and short-dependence obstruction identify why a proof based only on short additive trades cannot establish such a generator.

Several other attempted routes were considered but did not establish a theorem: coarse modular residue bucket counts lose the joint constraints at finer scales; low-order Fourier marginals do not certify any specified missing point; sparse convolution outputs can be as large as the entire input product; a direct ROABP encoding retains enormous partial-sum state spaces; and the straightforward conditional-expectation potential in section 2 has M^K witness patterns. None is claimed as a computational lower bound on all additive avoidance algorithms. No runtime improvement is inferred from the low-degree probe obstruction, which concerns a restricted **input-independent** family of probe generators.

Literature checked in this round includes Guruswami–Lyu–Wang, *Range Avoidance for Low-Depth Circuits and Connections to Pseudorandomness* (APPROX/RANDOM 2022, https://doi.org/10.4230/LIPIcs.APPROX/RANDOM.2022.20); Gajulapalli et al., *Range Avoidance for Constant Depth Circuits: Hardness and Algorithms* (APPROX/RANDOM 2023, https://doi.org/10.4230/LIPIcs.APPROX/RANDOM.2023.65); *Range Avoidance and Remote Point: New Algorithms and Hardness* (ITCS 2026, https://doi.org/10.4230/LIPIcs.ITCS.2026.79); and Bringmann–Fischer–Nakos, *Deterministic and Las Vegas Algorithms for Sparse Nonnegative Convolution* (SODA 2022, https://doi.org/10.1137/1.9781611977073.119, https://arxiv.org/abs/2107.07625). These verify the relevant distinctions between explicit hitting lists, algorithms for bounded-locality Boolean circuits, and output-sensitive convolution; they do not provide the needed generator for arbitrary explicitly listed modular summands. Classical interpolation, the elementary union bound, meet-in-the-middle and range counting are not claimed as new. Independent priority review of the precise additive-probe statements and independent proof-chain assurance for the inherited elliptic contraction premise remain outstanding. The ePrint PDF could not be freshly fetched via web in this pass; the reduction was checked against the exact repository notes cited above rather than treating the upstream theorem as independently reverified.

**Success verdict:** no fixed positive end-to-end deterministic saving exponent has been established. The uniform generator is still missing, so the mandatory breakthrough checklist fails. No claim of practical cryptographic performance, of a new computational lower bound on ECDLP, or of a stronger theorem for the actual selected elliptic masks is made.
