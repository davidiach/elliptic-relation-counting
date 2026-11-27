# Higher-order explicit probes and high-growth prescribed covers

27 November 2026. Abstract finite-plane follow-up to [PR #22](https://github.com/davidiach/elliptic-relation-counting/pull/22), [PR #23](https://github.com/davidiach/elliptic-relation-counting/pull/23), and [PR #24](https://github.com/davidiach/elliptic-relation-counting/pull/24). Reviewed main at `ea9bd4101039c4d83deedb52a20794e11a2d0120` (merged PR #20), and the three cited open research branches. This note concerns additive sets in the abstract plane over a prime field. It gives a constructive near-linear two-factor probe theorem, and a high-growth restricted-class lower bound using short-seed randomized covers of arbitrary prescribed points. It does not solve uniform additive avoidance for arbitrary growing arity, improve an end-to-end ECDLP exponent, or implement an elliptic-curve algorithm. Independent proof and priority review remain open. All finite tests are mathematical toy controls.

## 1. Explicit higher-order probes for two factors

Let G=F_r^2 with r an odd prime. A probe set H is universal for two lists of size at most n if no A+B with |A|,|B|<=n contains H. The Sidon/parabola result of PR #24 uses O(n^(3/2)) points. Higher-order additive uniqueness yields a stronger constructive bound when n is small relative to r.

**Theorem 1.** Let h>=2 be an integer, and set

    K = 1 + n(1+ceil(n^(1/h))).

Choose a prime q >= max(K,h+1,3), and put B=h(q-1)+1. If

    r >= B^ceil(h/2),

there is an explicit universal probe set H of size K. The construction runs in polynomial bit time in n,h,log r when q is chosen as the least prime above the indicated lower bound. No search over candidate sets or list descriptions is needed. The same theorem applies to a partition of any number of input factors into two groups whose actual partial supports both have size at most n; finding or certifying those bounds is a separate task.

**Proof (graph bound).** A B_h-set means that every sum of h elements (with repetitions allowed) uniquely determines the unordered multiset. It is also a B_j-set for every j<=h by padding both sums with the same h-j elements. If H is a B_h-set, choose one representation of each point in H intersect (A+B). Regard these as edges in a simple bipartite graph on the labels of A and B. Distinct covered points give distinct edges. Any cycle of length 2j<=2h gives an equality between the sums of the two alternating j-element multisets of distinct points of H, because each vertex label occurs once on each side. This contradicts the B_j property. Hence the graph has girth at least 2h+2.

The classical even-girth Moore bound for irregular graphs (Alon–Hoory–Linial, *The Moore bound for irregular graphs*, Graphs and Combinatorics 18, 2002, [DOI](https://doi.org/10.1007/S003730200002)) gives, for a graph with minimum degree >=2, average degree dbar, and girth >=2h+2,

    v >= 2 sum_{j=0}^h (dbar-1)^j.

If the original average degree is less than two, the edge count is less than 2n and the desired bound holds. Otherwise prune vertices of degree zero or one until the nonempty minimum-degree-two core remains. This cannot reduce average degree, and the core still has at most 2n vertices. The last term of the Moore inequality implies that the original average degree is at most 1+n^(1/h). The number of edges, including the original graph before pruning, is consequently at most n(1+n^(1/h)) < K. Thus any B_h-set of size K is a universal probe set. This use of high-girth extremal bounds is a generalization of the C4-free argument in PR #24, not a new Moore bound.

**Proof (explicit B_h-set without logarithm computation).** For t=0,...,K-1 form h canonical field digits

    d_j(t) = t^j mod q in [0,q-1],   j=1,...,h.

In F_q^h, the vector of power sums through degree h of h parameter values determines their elementary symmetric functions via Newton's identities, since all 1,...,h are invertible in F_q. The resulting monic polynomial determines the multiset of h parameter roots, including multiplicities. Thus the moment-curve vectors are a B_h-set in F_q^h.

Pack the first ceil(h/2) digits into one integer using base B, and the other digits into a second integer using the same base. Reduce the pair modulo r. In the sum of any h such pairs, each base-B digit sum is at most h(q-1)=B-1; there are no carries, and either full coordinate sum is at most B^ceil(h/2)-1<r. Hence equality of h-term sums modulo r implies equality of every integer digit sum, then equality of every power sum modulo q, and therefore equality of the multisets. The packed set has K distinct points since its first digit distinguishes t. This proves the theorem. Bertrand's postulate bounds the chosen q by twice max(K,h+1,3), and even trial division prime search is polynomial in n,h; the powers and packing require polynomial bit operations per output point. QED.

For concreteness, with b=ceil(log_2 r), L=ceil(log_2(n+b+2)), take the even order

    h = 2 floor((b-1)/(12 L))

whenever this is >=2. The chosen q is O(n^(3/2)+h+1), B<=h q, and log_2 B<=5 L (using L>=2 and h<=b). Thus (h/2) log_2 B <=5(b-1)/12<b-1, verifying the guard. For n=r^o(1), h tends to infinity and the probe size is

    O(n exp(O(log n * log(n+b+2) / b))) = n^(1+o(1))

when n also tends to infinity; if log n=O(sqrt(b)), the size is O(n). In particular for two factors each of size n=r^(1/m) with growing m and log n>>log b, the bound is O(n^(1+O(1/m))). This is a deterministic polynomial-time construction in this regime, not merely a nonuniform count. The two-factor lower bound from the translated-singleton argument is 2n for equal sizes, so the new bound approaches the optimal exponent for small n. The numerical constants in the guard are intentionally conservative.

**Arity limitation.** This does not make a balanced split of a general m-list high-growth input cheap: when each half-support is close to r^(1/2), the field guard does not admit growing h. Also the witness-incidence structure for m factors is a hypergraph, not a bipartite graph. As shown in #22, high-girth one-hot witness arrays may avoid short dependencies at much larger sizes. We do not claim a generic short-dependency theorem or a uniform sublinear-in-arity construction. Higher-order Sidon sets and Newton identities are classical; background on B_h sets: Kevin O'Bryant, *Constructing Thick B_h-Sets*, [arXiv:2308.12406](https://arxiv.org/abs/2308.12406), JIS 2024. We did not establish literature priority for the combination and the stated parameter regime.

## 2. High-growth, small-bias inputs can still cover almost N prescribed probes

The remaining case of #24 cannot be disposed of just by assuming that the sets look random under several strong additive diagnostics. We construct a restricted input class with simultaneously near-maximal support at every retained subset of blocks, near-maximal individual doubling, small Fourier bias, and complete scalar projections of the total sumset. Any prescribed probe set of size nearly the total list length can be covered by some member of this class. This strengthens the high-growth side of the universal-size lower bound; it does not increase the unrestricted lower bound of #24.

**Theorem 2.** Let r tend through odd primes and m=m(r) tend to infinity with m^2 log(m+1)=o(log r). Put

    q = ceil(r^(1/m)), n = m^3 q, d = 2q-1, ell=n-d,
    N=mn, M=n^m, Q=r^2, P=ell(ell+1)/2,
    Delta=2m, t=16m, delta=sqrt((M+mP)/Q).

For all sufficiently large r, for every prescribed H subset G with |H| <= 1+(m-2)ell, there are m distinct-point lists A_i of exactly n points each such that:

1. H is contained in S=sum_i A_i, and S contains a translate of both complete coordinate axes. In particular every nonzero scalar linear projection of S is all F_r.
2. M=r^(1+o(1))<Q. For every nonempty I subset [m], the support of S_I=sum_{i in I} A_i obeys

       |S_I| / n^|I| >= 1 - m d/n - 2(m-2)/ell - delta = 1-o(1).

3. Each list satisfies

       |A_i+A_i| / binom(n+1,2) >= 1-2d/n-delta = 1-o(1).

   This is nearly maximal doubling, not a low-doubling hypothesis.
4. For every nontrivial additive character chi of G and every i,

       |(1/n) sum_{a in A_i} chi(a)|
            <= d/n + 2t ell^(1/8) sqrt((2Delta-1)ell)/n = o(1).

There is a randomized algorithm taking polynomial time in N,m,log r that constructs these covers with probability 1-o(1), from O(m) independent uniform field elements (O(m log r) bits of entropy and expected fair bits using rejection sampling). H is supplied explicitly. This is a success-probability guarantee, not a deterministic certification of the realized seed's global support or Fourier properties. The statements concern abstract input lists, not the distribution of the actual projected relation blocks.

### 2.1 A sparse graph supplies many prescribed representations

Write g=max(4,2(m-2)). Construct a simple bipartite graph with ell vertices in each part, (m-2)ell edges, maximum degree <=Delta and girth >g. This is possible by a direct greedy argument under the sufficient guard ell>4 Delta^g. As long as fewer than (m-2)ell edges have been inserted, fewer than ell/2 vertices in either part have degree Delta. Pick an unsaturated left vertex. The number of vertices at graph distance <=g-1 is at most 1+Delta+...+Delta^(g-1)<=2 Delta^(g-1)<ell/2. Thus an unsaturated right vertex at distance >=g (or in another component) exists. Insert the edge. The distance condition prevents a new cycle of length <=g; it also forbids repeated edges. Repeat. This costs polynomial time in ell,m by ordinary breadth-first search and scanning vertices. Our asymptotic guard holds since log ell >> m log m. See also Chandran, *A high girth graph construction*, [DOI](https://doi.org/10.1137/S0895480101387893), for related distance-based constructions; no use of a novel graph existence theorem is claimed.

Partition the edges into m-2 indexed groups E_i, each of ell edges. There are two auxiliary lists, one for each part of the graph, and one pivot list for each edge group. Associate to vertices formal random labels X_u and Y_v in G. For each edge e=(u,v) in a pivot group set its pivot label equal to

    c_e - X_u - Y_v,

where the constants c_e are chosen to include all prescribed target points other than one reference. The two auxiliary random parts consist of all X_u and Y_v. Each pivot random part has ell labels, one per edge. If H is nonempty, pick h0 in H and assign each element of H-h0 other than zero to a different edge; pad unused constants arbitrarily. For empty H take h0=0. Every list also receives the d fixed digit-axis labels

    D_i = {(j q^i,0):0<=j<q} union {(0,j q^i):0<=j<q},  i=0,...,m-1,

reducing mod r. They are distinct since q<r and q^i is a unit, except for the one shared origin. Choosing X_u,Y_v and the pivot on their edge, and choosing zero in every other list, represents c_e. Zero itself is represented using all zeros. The sum of the fixed D_i contains both full coordinate axes by base-q expansion of integers 0,...,r-1. Finally translate the entire first list by h0. This preserves cardinalities, growth ratios and character magnitudes and translates the total sumset by h0. Hence the cover and projection claims hold as soon as the list cardinalities are correct. No randomness or membership test is used to choose the targets.

### 2.2 Injective formal signatures on most tuples

For a retained subset I choose one random-part label per retained list, excluding choices of fixed labels. Call the tuple good if neither selected auxiliary vertex (when its list is retained) is an endpoint of any selected pivot edge. Given the pivot choices, at most m-2 vertices are forbidden in each auxiliary part. Thus there are at least ell^|I|(1-2(m-2)/ell) good tuples, by the union bound over the two possible auxiliary choices.

Ignore the constants c_e and take the coefficient vector of a tuple sum in the 2ell formal variables X_u,Y_v. In a good tuple the positive coordinates exactly identify its selected auxiliary vertices. Its negative coordinates give the degrees of the selected edges in the graph. Suppose two good tuples for the same retained I have identical coefficient vectors. The auxiliary choices agree, and their edge sets have identical degrees. Removing the common edges gives two disjoint edge sets with equal degrees at every vertex. If nonempty, the underlying symmetric difference has minimum positive degree at least two and hence contains a cycle with length at most the number of its edges, which is at most 2(m-2)<=g. This contradicts the girth. So their edge sets are identical, and the partition into pivot groups identifies each tuple choice. Therefore all good tuples for each fixed I have distinct coefficient vectors. We do not need injectivity across different retained subsets.

There is a similar injectivity for unordered pairs of random-part labels from a *single* list, including repeated labels. For an auxiliary list this follows immediately from its formal variables. For a pivot list two unordered pairs of edges giving the same coefficient vector have the same endpoint multiplicities. Cancel common edges as multisets. If a nonempty signed difference remains, the underlying simple support has no degree-one vertex: equal signed degrees force at least two incident distinct edges at any vertex. Hence it contains a cycle of length <=4. The graph forbids this. Thus all P unordered pair signatures are distinct in each list. This explains the use of g>=4 in the graph construction.

### 2.3 Short-seed randomness, simultaneous collision estimates, distinctness

Choose two independent random polynomials over F_r, each of degree <kappa, with kappa=4t=64m; all coefficients are independent and uniform. Evaluate them at 2ell distinct indices in F_r, and interpret the pair of evaluations at each index as its formal plane label. This is kappa-wise independent over the plane: the Vandermonde matrix at any s<=kappa distinct indices has row rank s, so the s evaluations of either polynomial are uniform and independent. We have 2ell<r eventually. The maximum number of distinct variables in a difference of two good tuple signatures is at most 4m, and its coefficients have magnitude <=m. They cannot all vanish modulo r when r>m. Accordingly every such distinct pair of good tuples collides with probability exactly 1/Q, regardless of the constants. The same conclusion holds for the unordered double-sum signatures, with at most eight variables, and for the difference of any two different intended labels in a list when at least one label is random. The two fixed axes have no duplicate entries beyond the shared origin, already removed.

Let C_I be the number of unordered colliding good-tuple pairs for subset I, and let J_i be the number of unordered colliding pairs among its P double-sum keys. The expectations obey

    E C_I <= ell^(2|I|)/(2Q),   E J_i <= P^2/(2Q).

Crucially, one can combine every retained subset in a single normalized collision potential, not lose a factor 2^m in the final bound:

    Z = sum_{I nonempty} C_I/ell^|I| + sum_i J_i/P,
    E Z <= ((1+ell)^m-1 + mP)/(2Q) <= (M+mP)/(2Q).

The last inequality uses d>=1. Markov gives Z<=delta except with probability at most delta/2. For each retained subset the number of distinct sums from its good tuples is at least the number of good tuples minus C_I, because a collision class of size a loses a-1 distinct values, at most binom(a,2). The same argument for the double-sum keys gives at least (1-delta)P distinct values. Therefore on this event,

    |S_I| >= ell^|I| (1-2(m-2)/ell-delta),
    |A_i+A_i| >= (1-delta) P.

Use (ell/n)^|I| >= 1-m d/n and P/binom(n+1,2) >= (ell/n)^2 to obtain items 2 and 3 of the theorem. An additional union bound ensures that every list has exactly n distinct points, except with probability at most m binom(n,2)/Q, using the label-difference uniformity above. Translation of a list does not change any of these cardinalities or collision counts within a fixed retained subset.

### 2.4 Simultaneous Fourier bias

The uniform distribution on the actual list is the uniform distribution on its intended indexed labels once distinctness is established. For a nontrivial additive character chi:G->C^*, the fixed d labels contribute at most d in absolute value to the unnormalized character sum. Partition the full edge set into at most 2Delta-1 matchings by greedy edge coloring: an edge conflicts with at most 2Delta-2 previously colored edges. Restrict the colors to each pivot group. Within a matching of L edges, any 2t labels depend on at most 4t disjoint formal vertices; by the chosen limited independence, those labels are jointly uniform and independent. Thus the character values Z_j have unit magnitude, mean zero, and are 2t-wise independent. The same holds for the random part of an auxiliary list as a single group. The constants c_e do not alter this property.

Expansion of E|sum_{j=1}^L Z_j|^(2t) gives only terms in which no index occurs just once, because a singleton's independent mean is zero. There are at most (2t)^(2t) L^t such terms: partition the 2t positions into at most t nonempty classes (at most (2t)^(2t) possible partitions), then assign distinct indices. All surviving term magnitudes are at most one. Consequently

    Pr(|sum Z_j| > 2t sqrt(L) ell^(1/8)) <= ell^(-t/4).

There are at most m(2Delta-1) matching groups across all lists, and Q-1 nontrivial characters. The union bound gives simultaneous success except with probability at most m(2Delta-1)(Q-1) ell^(-t/4). On success sum the bounds across the matching groups of each pivot list, using sum sqrt(L_j) <= sqrt((2Delta-1)ell). Dividing by n and adding the fixed labels proves item 4. The Fourier bound is purely a probability argument; computing the maximum over all characters is not part of the sampler.

### 2.5 Asymptotics, resources and lower bound

The hypotheses give q=r^(1/m)(1+o(1)), n=m^3 q, ell=n-d, and

    log M = log r + O(m log m + m/q) = (1+o(1)) log r.

Thus M/Q and mP/Q tend to zero. We have ell>=q, so ell^(-t/4)=ell^(-4m)<=r^(-4), and all three exceptional probabilities above tend to zero. Also m d/n<=2/m^2, 2(m-2)/ell=o(1) and d/n+2t ell^(1/8) sqrt((2Delta-1)ell)/n=O(m^(-3)+m^(3/2)ell^(-3/8))=o(1). This proves simultaneous satisfaction of the properties with probability 1-o(1). Constructing the sparse graph, evaluating the two degree-(64m-1) polynomials at the 2ell indices, and creating all the lists have polynomial cost in N,m,b=ceil(log_2 r). The random coefficients have 128m log_2 r bits of entropy, and rejection sampling requires O(m b) expected fair bits. No enumeration of any large sumset, or of every Fourier character, is needed to generate the random instance. Conversely we do not claim that a successful coefficient seed can be *found deterministically* with these same resources, or that the global guarantees are cheaply checkable. The proof is independent of any contraction or other algorithmic premise. QED.

**Corollary (near-linear lower bound in a strongly generic class).** Consider only the inputs satisfying the quantitative promises of Theorem 2. Every set H of size at most 1+(m-2)ell is covered by at least one such input. Thus the minimum size of a fixed universal probe set for just this restricted class is at least

    2+(m-2)ell = (1-o(1)) N.

The known union-bound upper bound for all lists from #22–23 remains (2+o(1))N, and hence also bounds this restricted class. This shows that even simultaneous near-maximal growth of all retained subsets, almost maximal doubling, Fourier small bias in every input list, and full scalar coverage do not force a universal probe set dramatically smaller than linear in the explicit list length. It is not a lower bound for input-dependent candidate selection, nor does it rule out explicit near-linear universal probes. It does not improve the unrestricted size lower bound of #24. In particular none of the pseudorandomness promises just established can be assumed to resolve avoidance by itself. Moreover since the total support contains translates of both full axes, any sum of two nonzero scalar multiples of that support is the whole plane; hence the low-growth/signed-growth branch of #24 fails on these instances. They are a separate obstacle from the geometric box-resistant construction of #23.

## 3. Verification and research assessment

Run `python -B verify.py` (Python 3.10+, standard library, assertions enabled). See [verification.txt](verification.txt). The verifier checks:

- 113,235 distinct multiset sums in six packed B_h fixtures, including odd h and an intentionally failed field-size guard;
- all 66,048 simple bipartite graphs on 3+3 and 4+4 vertices for the rounded Moore edge caps at h=2,3; 124 finite two-factor cover exclusions, with an explicit large intersection of 2n-1 probes per fixture;
- 32,730 good-tuple signatures and 234,957 double-sum signatures, with exhaustive enumeration on two small graph fixtures and the indicated sampling on larger ones (the large m=5 double-sum check is exhaustive); exact graph girths and matching decomposition bounds;
- 25 polynomial-evaluation uniform marginals, 820 nonzero finite-field linear forms and 16 combinatorial character-moment checks;
- four complete finite covering examples, with exact sumsets on every retained subset and every within-list doubling, and all intended good-tuple collision counts. Two include full translated axes; the two larger-field controls omit axes and illustrate only the high-growth mechanism. All four cover their prescribed translated parabola points.

Some first finite greedy graph fixtures failed to reach the requested size at their chosen small parameters. The proved sufficient asymptotic size guard was not met by those fixtures. The final graph fixtures use exactly the required girth threshold max(4,2(m-2)) and independently check girth; the larger sixth-block case uses ell=1000. These adjustments are not evidence that the theorem's guard is necessary. Finite controls do not exhaust the randomness of the large construction or empirically certify its Fourier maximum. The character estimates are proved analytically. Pseudorandom Python seeds are used solely for reproducible finite tests, not claimed to implement the exact mathematical uniform distribution at scale.

Related literature checked: Chattopadhyay and Liao, [*Extractors for Sum of Two Sources*](https://arxiv.org/abs/2110.12652), 2021, studies extraction from XOR sumset sources rather than complement hitting for odd-prime modular sums; Guruswami, Lyu and Wang, [*Range Avoidance for Low-Depth Circuits and Connections to Pseudorandomness*](https://doi.org/10.4230/LIPIcs.APPROX/RANDOM.2022.20), gives range-avoidance constructions under different circuit restrictions; Huang, Li and Zhong, [*Range Avoidance and Remote Point: New Algorithms and Hardness*](https://doi.org/10.4230/LIPIcs.ITCS.2026.79), studies general complexity-theoretic connections. None of these results is used to claim a general small deterministic generator for the present additive model.

**Assessment.** The positive result substantially improves the explicit two-factor probe size where the field-size guard allows growing order; the other result gives a short-seed construction of high-growth counterexamples to the idea that generic local additive behavior alone drastically reduces the universal probe requirement. It makes the high-growth frontier more precise, but the desired uniform small probe construction for arbitrary growing-arity inputs remains missing. No new end-to-end complexity claim, no cryptographic-scale experiment, and no general lower bound on deterministic avoidance follows. The classical Moore bound, B_h construction by power sums, graph greediness and limited-independence sampling are not claimed novel individually; independent literature priority of the combined statements remains open.
