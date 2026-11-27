# Universal candidate sets and geometric barriers for planar additive range avoidance

27 November 2026. Follow-up to the abstract additive-avoidance problem isolated in [PR #18](https://github.com/davidiach/elliptic-relation-counting/pull/18), the box-cover tests in [PR #15](https://github.com/davidiach/elliptic-relation-counting/pull/15), and the improved coverage reduction under review in [PR #20](https://github.com/davidiach/elliptic-relation-counting/pull/20). Main reviewed at `1616184a390f3c93844fd327c152d392783195d0`. This note concerns explicit sets in the abstract additive plane, not any deployed cryptographic system. No general deterministic additive-avoidance algorithm, complexity improvement for discrete logarithms, or independent priority claim is made.

## 1. Near-optimal nonuniform candidate sets

Let G=F_r^2, Q=r^2, and fix nonempty list sizes n_1,...,n_m with M=product n_i<Q and N=sum n_i. Treat the lists as sets (deduplication can only help the cardinality promise). A universal candidate set H is a subset of G not contained in ANY sumset of m sets of these specified sizes. It is a hitting set for the complements, not necessarily an algorithm for identifying which member is missing. Let h_min be its minimum size.

**Theorem 1.** With any consistent logarithm base,

    N-m+2 <= h_min <= min(Q, 1+floor((N-m+1) log Q / log(Q/M))).

In particular, when M=r^(1+o(1)) and m=o(N), the minimum universal candidate size is between (1-o(1))N and (2+o(1))N. These are bounds for a set fixed across all possible inputs, not lower bounds on input-dependent candidate selection.

Proof of the lower bound. If 1<=|H|<=1+sum_i(n_i-1), choose h0 in H. Distribute the other elements of H-h0 among the m sets with capacities n_i-1, and include zero in each set. Their sumset contains H-h0 by taking zero in every other set. Translate one set by h0 and extend each set to its prescribed size. Thus H is covered. The empty set cannot be universal either. The lower bound does not require field structure.

Proof of the upper bound. Translate each input list to contain zero, storing the sum of the translations separately. This gives at most

    Q * product_i binom(Q-1,n_i-1) <= Q^(N-m+1)

possible sumsets, counting duplicates harmlessly. Each has size at most M. For t independent uniform candidate points, the probability that all lie in any fixed sumset is at most (M/Q)^t. The union bound is strictly below one as soon as t>(N-m+1) log Q/log(Q/M). Deduplicating the successful tuple does not destroy the property. Alternatively the entire plane always works. QED.

### 1.1 Graph candidates with prescribed first coordinates

A related nonuniform result fixes the first coordinates in advance. Let X be any t distinct residues (t<=r). There exists a universal candidate graph {(x,y_x):x in X} whenever

    t*r > M,     t log(t*r/M) > 2(N-m+1) log r.

Indeed for any given sumset S let s_x=|{y:(x,y) in S}|. Independent uniform choices y_x are all in S with probability product_x(s_x/r), at most (M/(tr))^t by AM-GM. Use the same count of sumsets as above. For instance, in the equal-size regime n_i=4 ceil(r^(1/m)), m growing with m^2=o(log r), one can take t=ceil(4mN)<r eventually: log(tr/M)=(1+o(1))log(r)/m. Thus even a graph on a short prescribed interval has a universal nonuniform witness. This proof neither supplies an explicit graph function nor selects a missing value from an explicit sumset.

## 2. A digit-expansion obstruction to algebraic candidate lists

**Theorem 2.** Let f:F_r -> F_r^2 be a polynomial map of total degree at most d>=1. For integers q>=2,s>=1, the prefix H={f(t):0<=t<q^s}, interpreting t modulo r, is contained in a sum of B(s,d)=sum_{j=1}^{min(d,s)} binom(s,j) sets, each of cardinality at most q^d. If m>=B(s,d), n>=q^d, and n^m<r^2, this yields a valid covering input with exactly m sets of size n by padding the sets and adding unused lists containing zero. It holds even when the coefficients of f are arbitrary field elements.

Proof. Express t=sum_{j=0}^{s-1}q^j z_j, z_j in {0,...,q-1}. For each nonempty subset I of digit positions, let

    F_I(z_I) = sum_{J subset I} (-1)^(|I|-|J|) f(sum_{j in J} q^j z_j).

Boolean inclusion-exclusion, or ordinary monomial expansion, gives f(t)=f(0)+sum_{I nonempty} F_I(z_I). A monomial in the expanded polynomial depends on at most d distinct digits; hence F_I is identically zero for |I|>d. The image of each surviving F_I has at most q^|I|<=q^d points. Translate any one component list by f(0). Every prefix point is then a sum of one element per component. QED. The same argument applies to multivariate polynomial maps when their input coordinates are linear combinations of a total of s radix digits. More generally it applies to any function with an additive decomposition into at most m components, each depending on at most d radix digits.

For fixed d and growing m take q=floor(n^(1/d)) and s of order m^(1/d) small enough that B(s,d)<=m. This covers prefixes of length n^(Omega_d(m^(1/d))) when n grows. For d=2, a parabola prefix of length n^(Omega(sqrt(m))) can be covered. In particular fixed-degree polynomial prefixes of size polynomial in n are not universal in the growing-arity sparse-product regime. This is a lower bound on the suitability of this PARTICULAR family of candidate sets, not a lower bound on polynomial-time construction of all candidate sets. It does not apply to arbitrary subsets of parameter values chosen across the entire field, to arbitrary rational parametrizations, or to all polynomial-based pseudorandom generators.

The graph upper bound in section 1.1 and the polynomial-cover theorem together isolate an interesting explicitness gap: nonuniform graph witnesses on a short interval exist, but bounded-degree polynomial graph functions on that interval can fail entirely.

## 3. A near-exponential combination barrier for factorwise rectangle covers

Call a cyclic affine rectangle any set B={x in G: ell_1(x) in I_1, ell_2(x) in I_2}, where ell_1,ell_2 are independent linear forms and I_j are cyclic intervals of residues (singletons and the whole field allowed). Each rectangle may use its own frame. For explicit lists L_i, consider a certificate which first covers each L_i by a union U_i of J_i such rectangles and then certifies that U_1+...+U_m is not the full plane. The uncompressed expansion has product_i J_i choices of rectangles. We bound this parameter, not the number of distinct sums after simplification or the time of any algorithm that might process the cover without expanding it.

**Theorem 3.** Let m=m(r) tend to infinity with m log log r=o(log r), b=ceil(log_2 r), A=ceil(r^(1/m)), n=4A, J=floor(A/(36b)). For all sufficiently large r there exist equal-sized lists |L_i|=n, whose product is r^(1+o(1))<r^2, such that:

1. Their sumset contains both full coordinate axes. In particular every nonzero scalar linear projection of the sumset is all F_r.
2. Every factorwise rectangle certificate of the stated kind must satisfy

       product_i J_i >= (J+1)^(m-1) = r^(1-o(1)).

   The conclusion holds even if different rectangles use different affine frames.
3. The lists can additionally be required to cover any prescribed candidate set H of size at most m(A-1)+1. Alternatively they can cover a polynomial prefix as in Theorem 2 whenever B(s,d)<=m and q^d<=A.

Proof. The number of nonempty cyclic affine rectangles is at most r^8: choose two linear forms (at most r^4 choices) and two cyclic intervals (at most r^4 choices). Empty rectangles are allowed for padding. The number of unions of at most J rectangles is at most (1+r^8)^J<=r^(9J). Choose A points uniformly without replacement from G. A fixed union of size at most Q/2 contains all A points with probability at most 2^(-A). Hence the probability that any such union covers the sample is at most r^(9J)2^(-A)<=2^(-3A/4)<1. Fix one good sample S of size A. Every union of at most J rectangles covering S has density strictly greater than 1/2. The same S may be included in every list.

For i=0,...,m-1 include the digit set

    D_i = {(j A^i,0):0<=j<A} union {(0,j A^i):0<=j<A},

with coordinates reduced modulo r. This costs at most 2A-1 points in each list. S union D_i costs at most 3A-1. The sum of all D_i contains each axis because all integers 0,...,r-1 have an m-digit base-A representation (A^m>=r). Thus every nonzero linear form maps at least one of those axes onto the field.

There is room for an additional set C_i of at most A points in each list. To cover H in item 3, use the lower-bound construction of Theorem 1 with capacities A, or the component lists of Theorem 2, with unused C_i={0}. Take the unions S union D_i union C_i and pad to exactly 4A distinct points. All promises are preserved. Since A grows faster than every fixed power of b, n^m=r^(1+o(1)) and is below r^2. The actual sumset is consequently proper regardless of the auxiliary H.

Any successful factorwise cover must have at most one factor with J_i<=J. Otherwise two U_i each have size >Q/2, so their sum is all G: for any g, U_i intersects g-U_j. Adding the other nonempty factors preserves the whole plane. Therefore at least m-1 factors must use J_i>=J+1. The resulting product bound follows from

    (m-1) log(J+1) >= (1-1/m)log r - O(m log b).

Both losses are o(log r). QED.

This also shows a linear universal-candidate lower bound even within a restricted class of full-projection, box-resistant inputs: every H with at most m(A-1)+1 elements can still be covered by a member of this class. The randomized existence argument does not identify a cheap box-resistant sample. The next section gives a deterministic construction for the narrower common-frame version.

## 4. Explicit common-frame obstruction via small-bias sets

A common frame within each factor is the original coordinate-box model; different factors may choose different frames for the following theorem. Use the same asymptotic parameters m,A as in Theorem 3 and fix 0<eta<1/4. Set epsilon=A^(-1/2+eta) and J'=floor(A^(1/4-eta)). Rounding epsilon downward by a constant factor, e.g. to a dyadic rational, is harmless.

By Theorem 1 of Jalan and Moshkovitz, [Near-Optimal Cayley Expanders for Abelian Groups](https://doi.org/10.4230/LIPIcs.FSTTCS.2021.24), FSTTCS 2021, there is an explicit epsilon-biased distribution on F_r^2 with support size at most poly(b)*epsilon^(-2-o(1))=A^(1-2eta+o(1))<=A eventually. Here use the cyclic group of order r and its square; its characters are the ordinary additive characters. The construction in that paper is the external algorithmic input, not implemented by our verifier. We only need its support, even if it is presented as a weighted sample.

For completeness the discrepancy transfer is elementary. Normalize the Fourier transform by uniform expectation on G. An interval indicator on F_r has Fourier l1 norm at most 2+ln r: the zero coefficient is at most one and the geometric-series bound for a nonzero frequency k is at most 1/(2 min(k,r-k)). A rectangle has Fourier l1 norm at most (2+ln r)^2, by tensor products and invariance of this norm under an invertible linear change of coordinates. A union of J' rectangles in one frame is a disjoint union of at most (4J'+1)^2 grid cells obtained by splitting all circular interval endpoints in each axis. Its Fourier l1 norm is thus at most (4J'+1)^2(2+ln r)^2. For any such union U, the epsilon-biased distribution mu obeys

    |mu(U)-|U|/r^2| <= epsilon (4J'+1)^2(2+ln r)^2 = o(1).

This follows by Fourier inversion and the bound on all nontrivial character expectations. Because log A/log b tends to infinity, the right side tends to zero. Thus if a union covers the support of mu, its density is greater than one half for sufficiently large r. Pad the support to A points if needed and repeat the construction in Theorem 3. Every successful factorwise certificate using a common frame within each factor needs

    product_i J_i >= (J'+1)^(m-1) = r^(1/4-eta-o(1)).

This yields deterministically specified obstructing input lists using the cited explicit small-bias construction, although the stronger near-r lower bound for arbitrary orientations in section 3 is only existential. The discrepancy proof does NOT handle mixed frames within one union of J' boxes. Related discrepancy constructions for Bohr sets appear in Adi Akavia, [Explicit small sets with epsilon-discrepancy on Bohr sets](https://doi.org/10.1016/j.ipl.2014.05.002), 2014. No claim of novelty for Fourier discrepancy or small-bias constructions is made.

## 5. Reproducibility and next mathematical target

Run `python -B verify.py` (Python 3.10+, NumPy for the accelerated exact integer scan in the three-block control; no optimizers). The recorded finite run checks four exact small universal-candidate minima, 3,264 lower-bound cover constructions, 780 exact polynomial prefix representations, 1,060 higher mixed differences, 3,166,776 linear-form arc checks in four full-projection examples, and an exhaustive arbitrary-orientation two-box test on F_7^2. The first patterned two-box sample failed the proposed property (a cover of area 24 was found); it is retained as a negative control. The replacement sample passes. These finite controls do not implement the asymptotic small-bias construction or reproduce its parameter bounds. The exact minimum-rectangle calculation checks all nonzero linear forms, including all nonzero scalings, and then takes the two smallest minima from distinct projective directions. The combined three-block example at r=997 contains both coordinate axes and a nine-point quadratic prefix. All four plane sumsets are proper by their list-product bounds. The standalone two-box sample is deliberately not claimed to be an input with sparse list product. See [verification.txt](verification.txt).

The broad total-search setting is related to Huang, Li and Zhong, [Range Avoidance and Remote Point: New Algorithms and Hardness](https://doi.org/10.4230/LIPIcs.ITCS.2026.79), ITCS 2026, and references in PR #18. The general results for Boolean circuit classes do not directly imply an efficient construction for these additive set families. The universal-set counting argument is classical, and wider priority of the particular boundaries here is unverified.

The new, focused open problem is to construct a universal complement-hitting set in explicit small size (or to obtain input-dependent certified avoidance) without using a candidate family susceptible to the digit expansion above. A set of near-linear size exists nonuniformly. The geometric obstruction shows that full scalar coverage and small list product alone cannot force an affordable factorwise rectangle certificate. This is not a lower bound on all algorithms: direct joint covers, different compressed representations, and richer input-dependent procedures remain possible. Neither the universal explicit hitting-set problem nor general planar additive avoidance is solved by this note. The previous conditional resource claims and their open upstream assumptions remain unchanged.
