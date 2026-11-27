# Constructive additive probes and a stronger covering lower bound

27 November 2026. Abstract finite-plane research, continuing PRs [#18](https://github.com/davidiach/elliptic-relation-counting/pull/18), [#20](https://github.com/davidiach/elliptic-relation-counting/pull/20), [#22](https://github.com/davidiach/elliptic-relation-counting/pull/22) and [#23](https://github.com/davidiach/elliptic-relation-counting/pull/23). Reviewed main at `ea9bd4101039c4d83deedb52a20794e11a2d0120`. The earlier probe-size bounds and polynomial interpolation obstructions in #22–23 are not claimed new here. We give constructive probe sets in several additional regimes and strengthen the universal size lower bound. These are finite-group statements. The general small uniform probe construction, general deterministic avoidance, and an improved end-to-end exponent for elliptic-curve discrete logarithms remain open. Independent proof and literature-priority review remain outstanding.

## 1. Setup

Let G=F_r^2 for an odd prime r, let m>=1, and let n_i>=1 be integer list-size bounds. Write M=product n_i<r^2, N=sum n_i, D=1+sum(n_i-1). A universal probe set H is not contained in A_1+...+A_m for any nonempty lists A_i of sizes <=n_i. We do not require that the generator identify which probe is missing. In the abstract problem, a certified missing point is a separate, stronger output requirement.

The known bound D+1 <= H_min follows by translating one probe to zero and distributing every other probe among the lists with capacities n_i-1, putting zero in each list. The known union-bound upper bound is min(r^2,1+floor(2D log r/log(r^2/M))). We improve on the constructive side under explicit conditions, not on that nonuniform upper bound in general.

## 2. Exact-minimum radix probes when the field is sufficiently large

**Theorem 1.** Set C=floor(sqrt((m+1)^D)) and B=C+1. If

    r >= B^ceil(D/2),

then H_min=D+1, with an explicit set constructible in polynomial bit time in D,m,log r. For m=2 the same conclusion holds with C=1, B=2. The case m=1 is immediate by selecting n_1+1 distinct points.

**Proof.** If H has a representation in the sumset, normalize each list by one reference point, collect the total shift, and assign a witness tuple to each point of H. Pad the lists by repeating references if they are smaller than their bounds. As in #22, this yields an integer incidence matrix T with one leading constant column and n_i-1 columns for each block, at most one 1 in each block per row, such that the two columns of H belong to col(T) over F_r. Thus T has D columns, and every row has squared norm at most m+1.

If there are D+1 witness rows, their rational row rank d is at most D. Choose d independent rows, d columns on which they are independent, and one additional row. The signed d-by-d minors on these d+1 rows furnish a nonzero integer dependence c on the selected columns, hence on all columns since the selected columns span the rational column space. The coefficient at the extra row is nonzero. Extend c by zeros on other rows. Hadamard gives |c_j|<=floor((m+1)^(d/2))<=C. The leading constant column gives sum c_j=0. We use the integer relation c^T T=0, not a claim that some minor remains a unit modulo r.

For m=2 there is a simpler, sharper argument: each witness is an edge between the corresponding two list labels, in a bipartite graph with n_1+n_2=D+1 vertices. Equal witness edges would give equal output points, so distinct probes have distinct edges. The D+1 edges cannot form a forest. A simple even cycle supplies an alternating signed dependence with coefficients in {-1,0,1}, which annihilates both full incidence blocks and the constant column. Thus C=1 suffices.

Choose H to consist of (0,0), k=ceil(D/2) points (B^j,0), 0<=j<k, and D-k points (0,B^j), 0<=j<D-k. They are distinct in F_r^2. For any nonzero integer c of sum zero and |c_j|<=C, a nonzero coefficient must occur away from the origin. In an axis with some such coefficient, take the largest exponent j with nonzero coefficient. Its magnitude is at least B^j, while all smaller contributions have total absolute value at most C*(B^j-1)/(B-1)=B^j-1. The coordinate sum is nonzero over Z and has absolute value <=B^k-1<r. Thus it cannot be zero mod r. This contradicts c^T H=0 mod r. The radix powers are constructed using integer arithmetic within the field-size guard. The earlier lower bound proves exact optimality. QED.

This argument is useful in the small-dimension regime, not for D vastly exceeding log r. The general sufficient guard is roughly D^2 log(m+1)=O(log r) with the appropriate constant. The two-list guard is substantially better. The proof does not say that these guards are necessary.

## 3. Sidon probes for a bipartition of the lists

**Theorem 2.** Let a,b>=1 bound the two partial-sum supports of a partition of the input lists (one may simply use the two products of their size bounds). Define

    E(a,b)=min(ab,
          floor((b+sqrt(b^2+4b*a*(a-1)))/2),
          floor((a+sqrt(a^2+4a*b*(b-1)))/2)).

If K=E(a,b)+1<=r, then H={(t,t^2):0<=t<K} is an explicit universal probe set for this size bound. In particular for two equally sized lists of size at most n it uses at most floor((n+n*sqrt(4n-3))/2)+1=O(n^(3/2)) points.

**Proof.** The parabola is Sidon for sums of unordered pairs in odd characteristic: equality of sums gives equality of the sum and product of the two parameters, hence the same unordered pair. For each distinct point of H that belongs to A+B choose one pair (a_i,b_j). Make these pairs the edges of a simple bipartite graph on the two lists. A four-cycle would give four distinct probe points satisfying the parallelogram identity; the Sidon property prohibits it. With e edges and right degrees d_j, every pair of left vertices has at most one common neighbor, so sum_j binom(d_j,2)<=binom(a,2). By Cauchy–Schwarz, e^2/b-e<=a(a-1). Reversing sides gives the second inequality; the trivial e<=ab also holds. Solve the quadratic inequalities to obtain E. Padding the lists for the size-bound interpretation does not change the argument. QED.

The underlying parabola Sidon property and C4-free bipartite graph bound are classical. This gives a short explicit universal list in the two-factor setting even when the radix theorem's field-size guard fails. For a balanced split with products near r^(1/2), the bound is near r^(3/4), too large for the growing-arity goal. It is not a better general algorithm than the input-dependent meet-in-the-middle baseline in #22. It also does not contradict #22's interpolation covers of short parabola prefixes.

## 4. Explicit near-linear probes for a certified low-growth class

For a positive integer h let hA denote the h-fold sumset, not scalar multiplication. There is a simple digit-basis H with hH=G, constructed as follows. Set q=ceil(r^(2/h)), k=floor(h/2). Include the k x-axis digit sets {(j q^i,0):0<=j<q} and k y-axis digit sets {(0,j q^i):0<=j<q}, for i=0,...,k-1. For even h this suffices because q^k>=r. For odd h add one joint set {(a q^k,b q^k):0<=a,b<s}, s=ceil(r/q^k)<=ceil(sqrt(q)). Reduce all coordinates mod r. Every pair of canonical residues has a representation using one point per digit set (the odd case uses the joint high digit), and all sets contain zero. Their union H has size at most hq+2 sqrt(q)+1 (at most 1+h(q-1) for even h). Conversely hH=G requires binom(|H|+h-1,h)>=r^2, hence |H|>=(h/e) r^(2/h)-h+1. The construction is optimal in order as r^(2/h) tends to infinity. It is an ordinary additive-basis construction, not a universal probe list for all inputs with M<r^2.

**Theorem 3 (growth criterion).** For S=A_1+...+A_m, suppose a partition of the lists into groups with partial-sum supports B_j is supplied, and for some h>=1 there are certified bounds U_j>=|hB_j| with product U_j<r^2. Then the explicit digit-basis H above is not contained in S. The same principle works for a signed pattern of h summands: negate the corresponding digit sets before taking their union to obtain the prescribed signed basis, and use upper bounds on the signed growth of the B_j.

**Proof.** hS=sum_j hB_j, with cardinality at most product U_j. If H were contained in S, then G=hH would be contained in hS, contradiction. The signed version follows by distributing the specified signs across each factor. QED.

This gives a new sufficient branch distinct from the rectangle covers in #15 and #23. The following two checks are charged explicitly:

* Compute higher growth by a capped exact support dynamic program. Begin with {0}; update T to T+B_j for h steps, retaining distinct points. Abort as soon as its cardinality exceeds cap C0. This is sound since adding any nonempty set cannot decrease support size. If the run completes, it gives U_j=|hB_j|. The update can use a deterministic sorted merge or bounded balanced search trees with polynomial overhead in log r, log C0, log |B_j|. At most h*C0*|B_j| pairs are examined for each group. With groups of at most g input factors, their B_j can be generated by enumerating at most N^g sums per group. No enumeration of S is needed in this test. For a fixed cap N^c and g=O(1) the work is polynomial in the explicit list length, h and log r. Trying larger g and h requires charging those parameters; no claim that a successful grouping always exists is made. Grouping can materially improve the cardinality product bound, as the finite check demonstrates.
* Alternatively compute K_j=|B_j+B_j|/|B_j| by explicit pairs. The classical Plunnecke–Ruzsa inequality (e.g. Petridis, arXiv:1101.2532; see also Ben Green, Additive Combinatorics 2025, Theorem 3.7) implies |h B_j|<=|B_j| K_j^h, and the same bound for a sum with any fixed pattern of plus and minus signs. Thus product_j |B_j| K_j^h<r^2 is a checkable sufficient condition using integer comparisons and no large support enumeration.

**Near-linear size consequence.** If M>=r and h=2m, the radix is q=ceil(r^(1/m)). By AM-GM, N>=m r^(1/m). Thus the explicit H has at most 1+2m(q-1)=O(N) points, irrespective of equality of the list sizes. In the sparse-product setting M<=r^(1+sigma), 0<=sigma<1, the ungrouped doubling condition is certainly met if product_i K_i<=r^((1-sigma-epsilon)/(2m)) for any fixed 0<epsilon<1-sigma. For constant K_i and m^2=o(log r), this holds eventually. This does not assume or imply bounded doubling of the actual complete list family in the elliptic-curve project. It is a rigorous uniform positive result for a subclass with a cheap observable certificate, not a general probe generator.

**Failure boundary.** The digit-axis input from #23 has m lists A_i={(j q^i,0),(0,j q^i):0<=j<q}, q=ceil(r^(1/m)). Its sumset contains both complete coordinate axes. For any two nonzero scalars alpha,beta, alpha S+beta S=G by taking a point from each axis. If m grows slowly and q is large, product |A_i|<r^2 still holds. Hence no test relying on a proper h-fold sumset, including signed or arbitrary nonzero scalar coefficients with h>=2, can succeed on all permitted inputs. This is a limitation of this criterion, not a lower bound on all avoidance algorithms.

## 5. A stronger universal-size lower bound from balanced-relation packing

This complements the explicit upper bounds in the first sections. Let all size bounds equal n>=2, let m>=2 and n^m<r^2. Choose an integer k>=2, and let s be the least integer >=k with binom(s,k)>r^2. Set

    g=max(0, floor((2(n-1)-s+1)/(2k))).

**Theorem 4.** Every set of at most

    1+m(n-1)+floor(m/2) g

points is coverable by m lists of size at most n. Therefore H_min>=2+m(n-1)+floor(m/2)g. If the bound exceeds the ambient size, it should be interpreted together with the original product promise; in fact the construction forces the coverable size to be at most the list product. The result is elementary and independent of any conjectures on addition circuits.

**Proof.** Normalize a nonempty target set by one of its points; the origin is covered by all zero choices. First consider a chunk T of at most 2(n-1)+g nonzero points assigned to a pair of lists. When |T|>=s, the k-subset sums of T cannot all be distinct. Two equal sums, after deleting common elements, give disjoint equal-size lists of j>=2 points U=(u_1,...,u_j), V=(v_1,...,v_j) with sum U=sum V and j<=k. (Distinct points exclude j=1.) Realize all 2j points by a cycle between two sets: start a_1=0 and recursively set b_i=u_i-a_i, a_(i+1)=v_i-b_i. The final a_(j+1)=0 by the equality of the sums. Thus u_i=a_i+b_i and v_i=a_(i+1)+b_i, using at most j-1 new left labels and j new right labels in addition to the always-available zero on each side. Swap the two roles if desired. A cycle covers 2j distinct target points at a capacity cost at most 2j-1, a saving of one slot.

Greedily remove disjoint cycles until fewer than s points remain. If the initial chunk size t exceeds 2(n-1), the number q_c of cycles satisfies q_c >= (t-s+1)/(2k) >= g when g>0, since each cycle removes at most 2k points. Thus the total cost of all cycles plus the remaining singleton points is at most t-q_c<=2(n-1). Orient successive cycles alternately so the total upper-bound capacity costs of the two sides differ by at most one. Each side's cycle cost is then <=n-1; distribute the remaining points as singleton labels to the unused capacity on either side. Coincident labels can only save more space. If t<=2(n-1), use just the singleton construction. This covers every chunk up to the stated capacity while leaving zero on both lists.

Partition the nonzero normalized targets among floor(m/2) pairs with those capacities and, for odd m, one last singleton-only list with capacity n-1. Every chunk is contained in its pair's sumset and every other list contains zero. Translate the first list back by the initial target point. The prescribed target set is covered. QED.

**Asymptotic consequence.** Suppose m tends to infinity, r<=n^m<r^2, and log n/(m log m) tends to infinity. Take k=2m+1. The inequality binom(s,k)>r^2 is ensured by s=k*(ceil(r^(2/k))+1), using binom(s,k)>=(s/k)^k. Hence the minimal s=o(n) under these hypotheses. It follows that g=(1-o(1)) n/(2m+1) and

    H_min >= 2+m(n-1)+(1-o(1))*n/4.

This strengthens the earlier linear covering threshold by an unbounded additive term. It does not close the factor-of-two gap to the nonuniform upper bound when m grows. Searching for collisions here is a proof of coverability, not a proposed fast general covering algorithm. The relation packing does not contradict the high-girth witness-incidence theorem in #22: we are constructing a covering representation for small target sets, not asserting a short dependence in every possible witness array.

## 6. Verification and assessment

Run `python -B verify.py` with Python 3.10+ (standard library only, assertions enabled). The two completed runs produced byte-identical `verification.txt`. The script tests 68,944 witness arrays for the radix construction, including exhaustive distinct-row arrays for four small fixtures and sampled arrays for the fifth; 15,400 Sidon probe exclusions (exhaustive up to translation for the first fixture, sampled for the others); 80 constructed cycle covers with 902 disjoint cycles; 28 complete plane-basis constructions; three positive higher-growth checks, one rescue from grouping and one explicit fast-growth failure. The test of the first version's largest radix fixture correctly failed its numerical field-size guard; the final fixture uses a verified larger prime. The tests are finite controls, not evidence that a generic large-list construction has been solved. They use no elliptic-curve computation, tensor contractor, or cryptographic instances.

Sources for the external classical facts: G. Petridis, *Plunnecke's Inequality*, [arXiv:1101.2532](https://arxiv.org/abs/1101.2532); Ben Green, [Additive Combinatorics, 2025](https://people.maths.ox.ac.uk/greenbj/papers/C3.10-Oxford-2025.pdf), Theorem 3.7. Related conceptual background includes Guruswami, Lyu and Wang, [Range Avoidance for Low-depth Circuits and Connections to Pseudorandomness](https://eccc.weizmann.ac.il/report/2022/102/), and Chattopadhyay and Liao, [Extractors for Sum of Two Sources](https://arxiv.org/abs/2110.12652). The latter concerns XOR-source extraction, not complement hitting for odd-prime modular sumsets; its guarantees do not solve this problem. Hadamard, elementary Sidon and C4-free bounds, radix bases, and Plunnecke–Ruzsa are classical; novelty of the precise combined statements has not been independently established.

**Remaining frontier:** explicit small probes for general high-growth inputs when the field is large but the total list dimension is far larger than log r. No universal sublinear-in-arity construction or improved end-to-end exponent follows from this note.
