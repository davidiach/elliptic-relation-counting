# High-precision lifting, degree-aware mode contraction, and the remaining packing gap

24 November 2026. Continuation of `section_norm_research_note.md`, Section 7. The goal is to improve the asymptotic general-prime relation-search exponent, not merely its leading coefficient.

**Result and scope.** The compressed mode network can now be contracted rigorously with its high inherited precision and high-degree vertices charged, by adapting the public Boolean contraction proof. Its resulting bound is still `2^O(M/sqrt(log M))` in the original edge count M. We give a degree-ratio version of that proof, an exact precision-removal lemma, and a sharper per-group entropy bound valid for uneven budgets. A sparse-interaction sufficient condition would improve the logarithmic ECDLP saving power to 1/3 (up to log-log factors), but the dense quotient graph of the current Fourier compiler fails that sufficient condition. We prove barriers for explicit dense mode elimination and for pair-only reciprocal-edge packet geometry. Neither is a lower bound on arbitrary implicit algorithms. No improved general-prime exponent or numerical coefficient is claimed.

Public theorem used throughout: *The Exponential-Time Hypothesis is False*, 25 October 2026, Boolean tensor contraction theorem and its proof, at
https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-reduction.tex
and
https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex .
The earlier numerical finite-packet assumptions are not invoked to derive a new constant. The claims adapting the public proof are subject to its correctness. Finite tests below do not run that asymptotic contractor.

## 1. Exact precision lifting

Let a Boolean network have S edges, and t tensors whose entries are not restricted to literal 0 or 1. Remove or multiply isolated scalar factors as appropriate. Work over the ring

    A = (Z/2^P Z)[U,V]/(f(U), V^ell),

where f is monic of degree e (the lemma also works without an unramified hypothesis on this original ring). The high precision P, e, ell, t, and S are polynomial in an ambient size M. A tensor entry is accessed through its canonical coefficient representatives. The other tensors may be equality predicates or other 0/1 tensors and are left unchanged.

**Lemma (digit lift).** The contraction modulo the original ring can be recovered from polynomially many contractions on the identical Boolean graph, at scalar dyadic precision

    B = 1 + S + ceil(t log2(P e ell)),

in an unramified extension of degree O(log(M+2)). The additional time and space, except for the contraction calls and original entry access, are polynomial in M and the description length. Each transformed entry calls its original entry procedure once and then performs polynomial-time arithmetic. When `t log(P e ell) = O(S)`, the input precision of all new contractions is O(S), no matter how large the original polynomially bounded P was.

**Proof.** Replace the canonical coefficient a_ij of an entry by its binary digit polynomial

    d_ij(Z) = sum_{k=0}^{P-1} bit_k(a_ij) Z^k.

For every entry use the *unreduced* polynomial `sum d_ij(Z) X^i Y^j` with integer coefficients. Do not reduce X by f or Y by Y^ell yet. Each polynomial has at most P e ell monomials with coefficients 0 or 1. Let G be the contraction over Z[Z,X,Y]. It has nonnegative integer coefficients, and each coefficient is at most

    2^S (P e ell)^t < 2^B.

The total degree in each variable separately is at most t(P-1), t(e-1), t(ell-1), respectively. Choose mixed-radix bases `A_Z=t(P-1)+1` and `A_X=t(e-1)+1` and substitute

    Z=W, X=W^A_Z, Y=W^(A_Z A_X).

No two monomials of G collide. The resulting degree bound is

    D_* = (t(P-1)+1)(t(e-1)+1)(t(ell-1)+1)-1.

Choose an unramified extension of Z/2^B with residue field of size at least D_*+1. If the next contraction routine needs more interpolation points, take the extension degree large enough for those as well. Since all parameters are polynomial in M, its degree is O(log M); a monic irreducible polynomial modulo two can be found by the same polynomial-time exhaustive procedure used in the public proof. Select D_*+1 distinct residue lifts. Their differences are units. For each point, evaluate the transformed entries locally and contract the unchanged graph. Lagrange or Newton interpolation over the extension ring recovers G's mixed-radix coefficient vector modulo 2^B, without any denominator loss. By the strict coefficient bound the recovered coefficients are their exact nonnegative integers, in the base-ring coordinate. Decode the exponents, substitute Z=2, X=U, Y=V, and reduce in A. Ring homomorphism shows this is precisely the original contraction. All interpolation objects have polynomial size; sequential calls require only polynomial additional space. QED.

This avoids a potentially huge initial precision penalty when the outer graph has very few nontrivial tensors. Equality subdivisions do not increase t. We use this only for the outer mode graph; it is not a claim that arbitrary dense tensors become easier just because their coefficient bit length shrinks. The scalar precision reduction also applies when the original output is over an ordinary dyadic ring (e=ell=1).

## 2. A degree-ratio version of the Boolean contractor

Here is a useful direct adaptation of the public proof. Let the graph have at most Cs edges and maximum degree Delta, with polynomial-size descriptions and uniform original entry-access time E and polynomial space. Define

    Delta_0 = max(Delta, ceil(s^(1/2))),
    R = max(2, s/Delta_0),
    U = log2 R.

Constant factors in the graph size can be included in s. Assume scalar input precision O(s) and polynomial-size coefficient arithmetic. The coefficient ring can already be an unramified extension with sufficiently many distinct residues, rather than just its scalar subring: the proof never needs scalar-only entries before the final coordinate extraction. In the applications here the extension degree and the size of its operations are polynomial in the ambient M.

**Proposition.** The public proof gives, with modified parameters,

    time <= poly(ambient description size) * (1+E)
                  * 2^O(s/sqrt(log(2+s/Delta))).

The space is polynomial if original entry access has polynomial space. The usual s/sqrt(log s) bound is recovered when Delta<=s^(1-epsilon) for fixed epsilon>0. The expression becomes an ordinary exponential bound when Delta is a constant fraction of s.

**Proof of the parameter changes.** If R is bounded, use direct enumeration. Otherwise fix the public number of groups r0=64 and a sufficiently small absolute epsilon, e.g. 1/10000. Stop a recursive branch once M'<=H=ceil(s/sqrt(U)). At depth d choose

    u_d=floor(epsilon U/2^d),
    a_d=2 sum_{j<=d} u_j,
    L_d=2^a_d, D_arr=r0^(2 a_d).

The unchanged packing and child construction make every child have at most M'/8 wires; the old degrees are unchanged and each new packet degree is at most 2 L_d. On an active path d<= (log2 U)/6+O(1), so u_d=Omega(U^(5/6)) and a_d<=4 epsilon U. Also `D_arr<=R^(48 epsilon)` and `L_d<=R^(4 epsilon)`. The public packing inequalities follow since

    D_arr^2 Delta_0 <= s R^(-1+96 epsilon) = o(H),
    L_d D_arr <= R^(52 epsilon) = o(H),
    2 L_d <= Delta_0

for sufficiently large R. The degree cap is thus maintained. A bounded exceptional R range can be handled by direct enumeration, with constants chosen uniformly.

Use the unchanged tapered Hermite readout and crossing-mode expansion. The inherited precision at a depth-d call is at most `s U^O(1) 2^a_(d-1)`. For a retained valuation i<u_d, put D=L_d/2^i. Then `D>=2^(a_(d-1)+u_d)`, and the extra 2^u_d dominates every fixed power of U. Since M'>H, the public precision calculation gives working precision O(M'D/u_d). The number of exponent lists is therefore `exp(O(M'/sqrt(u_d)))`, exactly as in the public partition calculation. Summing this logarithm along a branch is O(s/sqrt U), because M' decreases by at least 8 and sqrt(1/u_d) grows by at most a fixed multiple of sqrt 2 per level. Polynomial interpolation factors over O(log U) levels contribute O(log s log U) bits to log time, which is lower order (or poly in ambient M).

The new packet entry calculation costs `exp(O(L_d log s))`, which is `exp(o(s/sqrt U))` because L_d<=s^(4 epsilon). Original tensor evaluations have their charged time E and polynomial space. At leaves the direct sum costs at most 2^H times these entry costs and polynomial overhead. The immutable precision-lift convention, interpolation, and parent-to-child identity are exactly those in the public proof. This yields the stated bound with log R comparable to log(2+s/Delta), including the clamp Delta_0>=sqrt(s). QED.

**Variant for a high-precision child.** For the particular application below, the first-stage child has size s'=O(M/r) with r polylogarithmic in M and original maximum degree at most M^(1-epsilon0). Its initial precision is at most M^(1+o(1)). More precisely it is at most `M * 2^O(log M/log r)` for growing r. The same proof can absorb this instead of using the digit lift on all its many vertices. Choose its own initial u'_0=Theta(log M) with a sufficiently small fixed coefficient depending on epsilon0. Then the smallest first-level D' among retained nodes is at least `2^u'_0 = M^Omega(1)`, which dominates the ratio of the inherited precision to s' times all polylogarithmic factors. At every active call its precision is consequently O(M' D'/u'_0), just as above. Later levels satisfy the same cumulative bound. All precision remains polynomial in ambient M. If the child is already a leaf, evaluate directly. For the polynomial-time original entry predicates and root packet entry costs `2^{M^o(1)}`, its time is `poly(M) * 2^O(M/(r sqrt(log M)))`. This variant is specific to the scale separation between the two packet widths; it is not claimed for arbitrarily huge initial precision on a high-degree outer graph.

## 3. Exact-budget compression with per-vertex entropy bounds

At a retained node of the first stage write M for the original wire count, r for its group count, u for its retained valuation levels, L for packet width, i for node valuation, and D=L/2^i. The public expansion bounds the working precision P by O(M D/u) and its mode multiplicity per weight by O(M/D). Every mode has two group endpoints.

Group modes by unordered endpoint pair e, and fix exact budgets h_e with sum h_e<P, as in the preceding section-norm note. Let c_e bound the number of modes of each weight belonging to pair e. Write N_e(h) for the number of lists in that block with exact weighted sum h. The generating function gives, for h>0 and any theta>0,

    ln N_e(h) <= theta h + c_e sum_{j>=1} -ln(1-exp(-theta j))
               <= theta h + c_e pi^2/(6 theta).

Thus, also for h=0 with N_e(0)=1,

    log2 N_e(h) <= K sqrt(c_e h), K = pi sqrt(2/3)/ln 2.

Skip impossible budgets, and encode feasible lists by ceil(log2 N_e(h_e)) bits, padding invalid codes by zero at either endpoint. For any budget vector the total outer edge count and the degree of a group g obey

    S <= K sqrt((sum_e c_e) P) + #nonempty pairs,
    Delta_g <= K sqrt((sum_{e incident to g} c_e) P) + (r-1).

These follow by Cauchy--Schwarz. The incident-budget sum can be used in place of P for sharper instance-specific bounds. No assumption of evenly distributed h_e is needed.

In the public packing proof, all ordered label-word classes are padded to the same load. At a fixed coordinate, the fraction of padded slots whose unordered symbol pair contains g is at most 2/r; for a specified pair of distinct symbols it is 2/r^2. The same fractions apply to packet counts because every packet has the same L slots and the same unordered pair at that coordinate. Each packet contributes at most 2^{i+1} tile modes per weight. Removing padding cannot increase these counts. Consequently

    sum_e c_e = O(M/D),
    sum_{e incident to g} c_e = O(M/(rD)).

This proves the uniform bounds

    S = O(M/sqrt u + r^2),
    Delta_mode = O(M/sqrt(r u) + r).

The previously noted lower bound Delta>=2S/r is not sufficient to run a contractor; the new upper bound is what we need. The square-root loss in r reflects the possibility that a large fraction of the entire weight budget sits on edges incident to one group. Actual balanced budgets may do better.

Subdivision of mode edges for bipartite form doubles S, adds only binary equality tensors, and does not increase the nontrivial degree. The exact-list unranking and budget enumeration have the same polynomial access and `P^O(r^2)` branch-count costs as in the preceding note.

## 4. A fully charged two-stage composition, without a new saving power

Assume the original graph has degree at most M^(1-epsilon0) and polynomial entry access. Take r=(log M)^c for any fixed c>0, choose a=2u with u a sufficiently small fixed multiple of log M/log r, and use one level of the public packing and tapered readout. The public inequalities hold because the label-array size is M to a sufficiently small fixed power and L=M^o(1). The child graphs have O(M/r) edges. By the high-initial-precision variant in Section 2, the cost of a group entry is

    E = poly(M) * 2^O(M/(r sqrt(log M))).

It includes exact decoding of its incident lists, creation-time packet tensor entry access, and the child's own contraction. It depends on no uncharged table of child entries. The outer graph has only r nontrivial tensors. Use the precision lift with P polynomial in M and extension/truncation sizes polynomial in M: the new scalar precision is O(S+r log M), hence O(s_0) with ambient `s_0=Theta(M/sqrt u + r^2 + r log M)`. The latter two terms are lower order for the chosen r. Each transformed group entry calls the charged E routine. The degree bound of Section 3 yields `s_0/Delta_0=Omega(sqrt r)` after the degree clamp. Apply the ratio-aware contractor, giving

    log time = O(s_0/sqrt(log r) + M/(r sqrt(log M)) + r^2 log M)
             = O(M/sqrt(log M)).

The last term charges all exact pair budgets, the polynomial factors charge unranking and precision interpolation. Root Hermite nodes and evaluations also number polynomially many. Sequential calls use polynomial space. The mode network is now an actual alternative contraction route within the public proof, rather than an uncosted oracle for an enormous-precision high-degree network.

But the factors still cancel:

    sqrt(u) sqrt(log r) = Theta(sqrt(log M)).

This is an upper-bound accounting, not a time lower bound. The new proof does not show that no better mode-network algorithm exists.

## 5. When a sparse interaction graph would improve the ECDLP exponent

For a fixed budget vector consider the graph on the r groups, each edge weighted by its binary index length. Let its ordinary vertex treewidth be t_w and maximum weighted degree Delta_mode. A tree decomposition of width t_w gives a carving decomposition with maximum cut weight O((t_w+1) Delta_mode): attach each graph vertex as a leaf to one bag containing it, and replace high-degree tree nodes by binary trees. A graph edge crossing a cut of the original decomposition has an endpoint in the separator bag intersection, since the two endpoint subtrees intersect at a bag containing that edge. Edges crossing cuts created by refining a node similarly touch its bag. Leaf cuts cost at most Delta_mode. This bounds every cut. Given this decomposition, contract bottom up, storing tables on each cut. Combining two tables requires at most the sum of their cut weights (and at most a constant multiple of the maximum cut weight) in free binary indices, so time is `poly(M)(1+E) 2^O((t_w+1) Delta_mode)`, and space at most `2^O((t_w+1) Delta_mode)` times polynomial factors. We allow exponential storage for this sufficient condition. For bounded t_w a decomposition can be found polynomially in r by fixed-width treewidth methods, or supplied; even an exhaustive search over bag trees of size poly(r) has log cost poly(r,log r), lower order here.

**Sufficient structural condition (not achieved).** Suppose one can pack the original relation network with the first-stage entropy, degree, entry, interpolation and balance guarantees above, and make the mode interaction graph have bounded treewidth at every required node and budget. Then choose `r=Theta((log M)^(1/3))`. We obtain

    log time = O(M sqrt(log log M)/(log M)^(2/3)),

including exact budgets and charged child entries. This would transfer, for M=O(m b) and the inherited verified relation reduction, to

    m = Theta((log b)^(1/3)/(log log b)^(1/4)),
    log expected ECDLP time = O(b (log log b)^(1/4)/(log b)^(1/3)).

The transfer assumes the existing relation compiler and coverage proof at that arity. It is a strict asymptotic improvement over b/(log b)^(1/4). The section-norm polynomial expansion has 2^m=poly(b) cost in this regime and m log b=o((log b)^2), so the earlier precision accounting has room for that arity; the graph sparsity condition is the missing new ingredient. More generally t_w=o(sqrt(r/log r)) is sufficient for this separator method to improve over M/sqrt(log M) at polylogarithmically growing r. This is not a claim that low treewidth is necessary for a different implicit contraction algorithm.

**Why ordinary packing of the Fourier compiler does not provide this.** The interaction graph of the *mode catalogue* contains the crossing quotient graph of the original edge assignments: an original edge from group g to h with g!=h supplies real row and column ports at the same packet position, hence a tile mode with endpoints h and g. Even if some particular interpolation evaluation or tensor entry kills it, the generic catalogue includes it. A complete bipartite Fourier boundary between two populated sides with balanced loads gives a complete bipartite graph between their nonempty group-label sets, omitting equal-label pairs. When that boundary has a fixed positive fraction of the original edges, both sets have Omega(r) elements by the group port caps (or directly for the fully populated uniform boundary). Choose disjoint subsets of size Omega(r) on the two sides, allowing overlap of the original sets. This supplies a K_{Omega(r),Omega(r)} subgraph of the catalogue. Its treewidth is Omega(r) (e.g. the minimum degree of that induced subgraph is Omega(r), and a graph of treewidth k is k-degenerate). Thus the natural catalogue is not of bounded treewidth. Some budget branches may delete many edges, and actual nonzero terms might admit compression. This argument does not prove they cannot.

## 6. Dense budget branches and explicit intermediate table obstruction

In a fully occupied balanced idealization, assume that each pair e has c_e>=c modes at each weight 1,...,J, with c=Theta(M/(r^2 D)), J=floor(c_* D/sqrt u), and a small fixed c_*>0. These are explicit model hypotheses about the mode catalogue, not a claim about nonzero coefficients of every input. For polylogarithmic r and the retained packet scales, choose c_* small enough that the total weight of selecting every one of these modes once is below P=Theta(M D/u). The exact-budget representation has a branch with simultaneously large alphabets: for each pair choose an h_e attaining the largest cardinality among the 0/1 exponent subsets. There are 2^(cJ) subsets and at most `1+cJ(J+1)/2` possible weights. Thus

    log2 N_e(h_e) >= cJ - log2(1+cJ(J+1)/2)
                  = Omega(M/(r^2 sqrt u))

for the growing parameters in question. All budgets are feasible simultaneously. This gives a complete weighted graph with each bundle carrying Omega(M/(r^2 sqrt u)) bits and total size Theta(M/sqrt u). A binary contraction tree grouping the original group vertices must have an internal cut with r/3<=|A|<=2r/3 (take the first node of size >2r/3 and descend to a suitably sized child). Its boundary contains |A|(r-|A|)=Omega(r^2) pairs, so materializing the entire intermediate tensor requires `2^Omega(M/sqrt u)` entries. In particular, this exact-budget method cannot rely on all its branches having small explicit boundary tables.

There is an analogous rank statement for a *universal* multilinear sum-of-products join. Treat each local table entry as an independent formal input. The coefficient tensor of the join is the product of equality tensors on the pair variables. Flatten across a balanced group cut. For each fixed assignment of the crossing edge labels, the two sides have disjoint supports in their local table-entry coordinates. Every such support is nonempty (choose arbitrary internal labels), and the corresponding block is a rank-one all-ones outer product. Different crossing assignments give disjoint blocks. Hence flattening rank is exactly the product of the crossing alphabet sizes, exponential in the cut weight. Any representation as a sum of fully separated products of linear forms in the independent group tables needs at least that rank. This does not bound general arithmetic circuits, nonlinear entry operations, algorithms exploiting actual group factors, or the public dyadic contractor (which uses entrywise canonical lifts).

## 7. A pair-only packet geometry shortcut fails

One tempting alternative to the label-array blowup is to put the reciprocal edges of each unordered pair of groups into their own wide packet. With a cyclic integer-shift deformation, encode the two orientations by a coloring c: Z/2^a -> {0,1}. A row in position x is in the same group as a column at x+t precisely when the orientation has flipped. All matched ports are local across the two groups iff

    c(x+t)=1-c(x) for all x.

Call such t an anti-period. The periods form a subgroup H of Z/2^a. If anti-periods exist, they form one coset t0+H, and 2t0 is a period while t0 is not. Since a cyclic 2-group has one subgroup of each order, H consists of multiples of 2^k for some k>=1 and the anti-period coset is exactly `2^(k-1) + H`. Every anti-period has valuation k-1. Thus a single two-group cyclic packet can localize its crossing pair at only one valuation level, no matter how many positions it has.

For the ordinary noncyclic interval [0,L), a related statement holds. If two positive shifts s,t have different dyadic valuations, write g=gcd(s,t). The indices 0,g,...,s+t-g form a cycle under rotation by t modulo s+t, in steps of either +t or -s, of odd length (s+t)/g. If L>=s+t-g+1, the full cycle is in the packet. An anti-period condition for both shifts would require alternating two colors along this odd cycle, impossible. In particular two shifts both at most L/2 cannot be anti-periods at different levels. This linear statement only concerns shifts with the stated overlap; the cyclic statement is unconditional for full packets.

Finally any two nontrivial partitions of the vertex set of a complete graph have a pair separated by both. Otherwise every pair crossing the first partition would be inside a block of the second partition, and the connected complete multipartite crossing graph would put all vertices in one block. Therefore in the cyclic reciprocal-pair-only scheme for a complete mode graph, two different valuation levels cannot both yield nontrivial globally localized partitions: some group pair would need a packet to have anti-periods at both levels. This identifies why many large reciprocal-edge packets do not by themselves supply the missing independent grouping levels. More general packets involving many group pairs, nontranslation kernels, and analytic approximations remain outside the claim.

## 8. Verification and research status

Run `python -B verify_mode_contraction.py` (Python standard library). The tests use independent direct contraction, exact sparse integer polynomial convolution, evaluation and Newton interpolation over an explicitly constructed unramified ring, and final reduction into the original ring. They verify scalar, extension and nilpotent cases. The remaining tests use exact colored-partition counts, uneven pair budgets, dense-budget witnesses, exhaustive small packet colorings and intersecting partition cuts. The transcript records 8 precision networks, 561 exact evaluations and recovered coefficients, 2,541 partition inequalities, 1,770 incident-degree inequalities, 4 dense-budget witnesses, 65,812 cyclic colorings, 8,188 linear colorings and 10,548 complete-graph cut pairs. These are finite identity checks; the general bounds are proved above.

**Next gate.** High inherited precision is no longer an unresolved obstacle in the two-stage route, and the group-degree bound now makes the second contraction fully chargeable. The remaining target is a specialized join of the group contractions, or a genuinely different packet geometry, that beats `S/sqrt(log r)` with S=Theta(M sqrt(log r/log M)). Equivalently, for r=polylog M one needs an extra unbounded saving beyond the ratio-aware generic contractor; a factor `(log M)^(1/6)` (up to log-log factors) would suffice for the 2/3 contraction saving power. Low-treewidth mode graphs are a precise sufficient route but the present dense catalogue fails it. The last proposed numerical coefficient 11.31486814 remains a conditional candidate, unchanged. The asymptotic general-prime objective remains open.
