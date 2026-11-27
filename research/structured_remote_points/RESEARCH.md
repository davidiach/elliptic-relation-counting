# Exact remote syndromes from hidden rank geometry and separator profiles

27 November 2026. Abstract finite-field/coding-theory follow-up to PR #25 (`911799c4f2c7d9fde822590db9391c4b32440d98`). We study explicitly supplied dictionaries, not real cryptographic instances. Two independently constructive results are proved: (i) recognition of a hidden bounded-rank matrix dictionary, including a noisy version and rank-metric probes for arbitrary outliers; (ii) exact min-plus profile elimination on a supplied tree of small linear interfaces. They are complementary **structured** algorithms, not a uniform algorithm for the general high-rate remote-point problem. No new end-to-end ECDLP time bound is claimed.

## 1. Abstract problem and definitions

Let q be prime, V=F_q^d, and H a finite explicitly listed dictionary of K columns in V. Repeated columns and zero may be removed when computing its sparse span. Define

    S_w(H) = { sum_(i in I) a_i h_i : |I| <= w, a_i in F_q }.

A remote syndrome is a vector outside S_w(H). If H spans V, and L=ker(H) is the length-K linear code with parity-check matrix H, the minimum Hamming weight in a syndrome coset equals the minimum number of dictionary columns in a representation of that syndrome. In particular the maximum of these minimum weights equals the covering radius. These equivalences are classical. The interesting regime from #25 has q polynomial in w, d very large compared with w^2, K exponential in d/w, and a certified volume slack. The crude bound |S_w(H)| <= (1+qK)^w supplies existence, not a deterministic efficient constructor. If H does not span V, ordinary elimination already finds a missing vector.

## 2. Hidden matrix rank: an algorithm not limited by polynomial degree

Write d=s^2 and let B_k be the full set of s-by-s matrices over F_q of rank at most k. Suppose first that H is B_k after an unknown invertible F_q-linear change of coordinates, where q>=7, 1<=k and s>=8k. The radius parameter w is independent of k. We prove a deterministic polynomial-time reconstruction of a rank coordinate system from the *list of vectors alone*. No coordinate arrangement or hidden basis is supplied. The numerical guard is sufficient, not asserted sharp.

### 2.1 Overlap counts reveal the rank-one set

For any finite set B write C_B(v)=|{x in B: x+v in B}|. For B=B_k, the count depends only on rank(v): multiplication on the left and right by invertible matrices acts transitively on each rank stratum and preserves B_k. Write C_1 for the count at any rank-one shift. We prove the quantitative gap

    C_B(v) <= C_1/16    for every nonzero v of rank >=2.       (1)

Let R_h be the number of rank-h matrices. Counting full-rank factors X (s-by-h) and Y (h-by-s), and dividing by the number of choices of a basis of the h-dimensional column space, gives

    R_h = prod_(i=0)^(h-1) (q^s-q^i)^2 / prod_(i=0)^(h-1) (q^h-q^i).

For q>=7, prod_(i>=1)(1-q^(-i)) >= 1-1/(q-1) >= 5/6, so

    (1/2) q^(h(2s-h)) <= R_h <= 2 q^(h(2s-h)).

Summing the geometric upper bounds for h<k yields

    sum_(h<k) R_h / R_k <= 8 q^(-2(s-k)-1).                    (2)

For a uniformly chosen rank-k matrix X, the column and row spaces are independent uniform k-dimensional subspaces. This follows from the independent transitive left/right action, or from the uniform number of invertible maps between fixed factor spaces. For rank-one shift uv^T, whenever u belongs to the column space of X, X+uv^T has rank at most k. Thus

    C_1/R_k >= (q^k-1)/(q^s-1) >= (1-1/q) q^(-(s-k)).       (3)

Now let rank(v)=j>=2. If j>2k the two balls cannot intersect. Otherwise write v=PQ with rank-j factors and X=UV with rank-k factors. Let a and b be the dimensions of the intersections of their respective column and row spaces. The ranks of [U P] and [V;Q] are k+j-a and k+j-b. Sylvester's rank inequality for their product gives rank(X+v)>=k+j-a-b. Hence a+b>=j is necessary for X+v to have rank at most k. In particular a+b>=2.

The probability that a uniform k-subspace contains a specified nonzero vector is at most 2q^(k-s); for a specified ordered independent pair it is at most 4q^(2(k-s)). Union over at most q^j nonzero vectors or q^(2j) ordered pairs in the fixed j-subspace. Independence of the two random spaces and a union bound for a>=2, b>=2, or both a,b>=1 give

    Pr(a+b>=2) <= 12 q^(-2(s-k-j)) <= 12 q^(-2(s-3k)).

We allow every lower-rank X to contribute. By (2),

    C_B(v)/R_k <= 12 q^(-2(s-3k)) + 8 q^(-2(s-k)-1).

Divide by (3). With s>=8k and q>=7 the ratio is at most

    (7/6) [12 q^(-(s-5k)) + 8 q^(-(s-k)-1)] < 1/16,

as claimed. This proof uses only elementary rank and subspace counting, not an imported assertion about the intersection numbers of the bilinear forms association scheme.

Consequently for a complete hidden dictionary all nonzero rank-one vectors are exactly the shifts with maximal nonzero overlap. Counts can be obtained by tabulating ordered differences of distinct input vectors in polynomial time in the explicit dictionary size and d,log q. Sorting gives a deterministic implementation. For the complete ball it suffices to compare the maxima, so even in smaller finite examples where the 1/16 gap fails the exact maximum method can succeed. It need not always succeed in small parameters: the F_2, 2-by-2, k=1 fixture has C(v)=6 at **all fifteen** nonzero shifts. The larger-parameter proof does not include that exception.

### 2.2 Reconstructing the two tensor factors

Suppose the full set R of nonzero rank-one vectors has been identified in hidden coordinates. Put A=R union {0}. For any h in R, define N(h)={x in A : x+h in A}. In the true matrix coordinates, write h=u v^T. The sum of two nonzero rank-one matrices has rank at most one exactly when their left factors are proportional or their right factors are proportional. In the nonproportional case independent left and right factors give rank two. Thus N(h) is the union of two s-dimensional linear spaces L_h={u z^T} and U_h={z v^T}, intersecting in the line F_q h.

Choose any a in N(h) outside that line (s>=2). The set of x in N(h) with x+a in A is exactly the component containing a. Call it L; the other component is (N(h) minus L) union F_q h, called U. This identification is valid in the unknown coordinates because the original change of coordinates is linear. Compute bases h=a_1,a_2,...,a_s of L and h=b_1,b_2,...,b_s of U. Their original matrices can be written a_j=u_1 v_j^T and b_i=u_i v_1^T after possibly interchanging the two factors; the u_i and v_j form bases of their respective factors.

For i,j>=2, identify the other factor space at a_j and at b_i by the same adjacency construction, removing L or U and restoring the corresponding common line. Their intersection is the one-dimensional line representing u_i v_j^T. Choose its unique nonzero element c_ij for which

    h + a_j + b_i + c_ij belongs to A.

In the two-by-two factor coordinates the matrix of this sum is [[1,1],[1,t]], which has rank one precisely at t=1. Thus the scalar normalization is forced. Set c_1j=a_j, c_i1=b_i. The resulting s^2 vectors are a basis, and the matrix coordinates in this basis preserve rank. This is a constructive linear-algebra proof, without assuming access to the classification of rank-distance isometries.

All set tests, basis calculations, and intersections are polynomial in the size of the explicit rank-one set, d and q (or in the number of line elements scanned). In the promised hidden ball this size is no greater than K. The reconstruction can be implemented with explicit failure checks. On an arbitrary input, even if overlap classification or reconstruction gives a spurious candidate layout, **validate** that its basis is invertible and calculate the actual matrix rank of every dictionary column in that layout. If every column has rank <=k and w k<s, the sum of its diagonal basis matrices is a remote syndrome, by rank subadditivity. No unverified inference about hidden coordinates is needed for the output. If the validation or resource guard fails, this branch reports failure.

### 2.3 A noise-tolerant recognition result

Let S be an independently p-thinned subset of B_k (retain each point with probability p), and then add any set O of at most L arbitrary outliers, even chosen after observing S. All statements can be transported by an unknown invertible linear map. Put mu=p^2 C_1. For any fixed nonzero v the directed pairs (x,x+v) with both endpoints in B_k have underlying undirected graph of maximum degree two, with no duplicate unordered edges for odd q. Partition the edges into at most three matchings. Within each matching the indicators that both endpoints are retained are independent Bernoulli variables of mean p^2. The total mean is p^2 C_B(v)<=mu. Applying the elementary additive Chernoff bound to each matching with tolerance mu/24, and a union bound, gives

    Pr(|C_S(v)-p^2 C_B(v)| > mu/8) <= 6 exp(-mu/2000).

For example the two-sided bound for a sum of independent Bernoulli variables of mean nu is 2 exp(-a^2/(2(nu+a/3))); here nu<=mu and a=mu/24, so the displayed constant is conservative. Union over all q^(s^2)-1 nonzero shifts gives simultaneous success probability at least 1-6q^(s^2)exp(-mu/2000). This expression may be negative in small cases and is interpreted as a useful lower bound only when positive.

Adding L outliers changes each overlap count by at most 2L, because for any outlier there is at most one ordered pair starting there and one ending there for a fixed shift. Suppose L<=mu/32. On the simultaneous event, every rank-one count lies in [7mu/8,19mu/16], and every nonzero higher-rank count is <=mu/4. Thus the set of observed differences with count greater than **half the largest nonzero count** is exactly the full rank-one set. In particular the rank-one shifts need not themselves have been retained as dictionary entries. Tabulation of ordered differences of the supplied set, including the outliers, recovers the needed set without knowing p, C_1, or the hidden factors. A deterministic alternative allows at most C_1/16 arbitrary deletions from the full ball and C_1/32 additions, since deleting one element affects a fixed count by at most two; take p=1 in the same inequalities.

For fixed p>0 the bound tends to one extremely rapidly in the growing-parameter regime in section 4, since (3) gives C_1>= (1/2)(1-1/q) q^(k(2s-k)+k-s). The random statement describes an input ensemble, **not randomness required by the recognition algorithm**, nor a promise that arbitrary dictionaries have this form. Even on a noisy input outside the theorem, correctness of a subsequent output follows from the independently checked ranks and the outlier certificate below.

## 3. Certified remote probes with arbitrary high-rank outliers

Suppose an explicit linear layout identifies all but L dictionary atoms as matrices of rank at most k. The outliers can be arbitrary matrices. Any w-sparse dictionary combination is a sum of a matrix of rank at most R=wk and a member of the w-sparse span of the outliers. Let C be an upper bound on the number of distinct outlier centers, for example

    C = min(q^t, sum_(j=0)^min(w,L) binom(L,j)(q-1)^j),

where t is the rank (linear span dimension) of the outlier vectors. Both bounds are computable from the explicit data. The first follows by allowing all elements of the outlier span, the second by counting choices of the actual sparse combination. The center enumeration itself is not needed to output a probe list.

Assume 2R<s and C<q^(s(s-2R)). There is a completely explicit list of C+1 matrices, pairwise rank distance >=2R+1. Indeed set ell=s-2R, use the degree-s extension F_(q^s), and for each coefficient tuple (a_0,...,a_(ell-1)) define the F_q-linear map

    f_a(x) = sum_(i=0)^(ell-1) a_i x^(q^i).

Represent it by its s-by-s matrix in a fixed extension basis. A nonzero difference polynomial has ordinary degree <=q^(ell-1), so at most q^(ell-1) roots. Its F_q-linear kernel has dimension <=ell-1, and the difference matrix has rank >=s-ell+1=2R+1. This is the classical linearized-polynomial rank-metric construction, with a self-contained distance proof. Its full size is q^(s ell); take the first C+1 coefficient tuples. By rank subadditivity, each rank ball of radius R about any center contains at most one of these probes. Thus at least one probe is outside S_w(H). This gives a certified *list* guarantee even without enumerating centers.

If a single certified vector is desired, enumerate the smaller of the two types of centers in the bound C (discard duplicate centers if helpful), test the rank distance of each probe to every listed center, and keep one whose rank distances are all >R. This uses polynomial work in C and the explicit data. If only the list is wanted, output it without any such verification of absence. The proof covers all coefficients in the original small field, not only one-hot sums.

All extension-field setup is charged: for the prime fields considered here, one elementary deterministic method enumerates monic degree-s polynomials and tests division by all monic polynomials of degree <=floor(s/2); it finds an irreducible polynomial in q^O(s) field operations. Thereafter representing and evaluating a fixed linearized polynomial is polynomial in s,log q. An already supplied irreducible polynomial can be checked and used instead. Thus the full cost is polynomial in K,d,q,q^s,C (with absolute exponents); it is not silently polynomial in log(q^s) for the exhaustive construction. The algorithm succeeds within a cap only when those quantities fit it. In particular, if the input is a constant fraction of the rank ball in section 4, q^s is polynomial in K. If L=poly(d) and w grows polylogarithmically in d, then log_q C=O(w log_q(d+q))=o(s) at the displayed scales, so the rank-metric list capacity and explicit enumeration costs are comfortably within polynomial-in-K time. More generally the exact inequalities above, not the qualitative claim of a few outliers, are the success criteria.

### 3.1 Unknown affine offsets are harmless under a stronger radius guard

The overlap statistic is unchanged by translating the entire hidden rank ball or its thinned subset. Thus the recognition theorem of section 2 identifies the linear tensor layout even if the inliers are an unknown affine translate a+B_k (after an invertible linear transformation). We do not have to identify a itself.

In the recovered coordinates, choose an actual supplied atom h_0 maximizing the number of supplied atoms h for which rank(h-h_0)<=2k. As long as there is at least one inlier, a choice of an inlier gives such a neighborhood containing every inlier, by rank subadditivity. Hence the maximizing choice leaves at most L exceptional difference atoms if there were at most L outliers in the input. These differences need not coincide with the original outliers. This step is an explicit pairwise rank calculation.

Write every original atom as h_0+(h-h_0). In any w-sparse combination the coefficient of h_0 is a single arbitrary scalar; the other part consists of a matrix of rank <=2wk plus a w-sparse combination of the exceptional difference atoms. All centers therefore lie in F_q h_0 plus that sparse span. Their number is bounded by

    C_aff <= min(q^(t+1), q * sum_(j=0)^min(w,L) binom(L,j)(q-1)^j),

where t is the span dimension of the exceptional differences. The first term can also use the actual rank of their span together with h_0. The rank-metric probe construction above applies with R=2wk, requiring 4wk<s and C_aff<q^(s(s-4wk)). Its correctness is again independently checkable from the explicit ranks and center bound. In the full affine-ball case there are no exceptional differences, and at most q+1 probes suffice under that radius guard. This is a genuine affine extension, not the erroneous assumption that sparse spans are invariant under translation.

## 4. A high-width, high-degree family with polynomial-time construction

These results are not subsumed by the earlier low-degree polynomial separator or low-width routes. To see this take q>=7, s tending to infinity with w=o(s), k=floor(s/(4w))>=1 (and impose s>=8k, e.g. w>=2). Use the full rank ball B_k as a dictionary of distinct matrices, possibly hidden by any invertible change of coordinates. Write d=s^2. Every matrix of rank <=k factors as an s-by-k times k-by-s product, so K<=q^(2sk); the rank-k count also shows log_q K=Theta(sk)=Theta(d/w). The uncolored volume bound is

    log_q (1+qK)^w <= w(2sk+2) <= d/2+2w <= 3d/4

for sufficiently large s. The code rate 1-d/K tends to one and the required volume slack is real, not an artefact of the actual support collapsing. All nonzero scalar linear projections of the dictionary are surjective because it contains every scalar multiple of every matrix unit.

**Exact sparse span and polynomial degree.** The sparse span is exactly B_(wk) when wk<s. One inclusion is rank subadditivity. For the converse, factor a matrix of rank <=wk into rank-one terms and group them into at most w blocks of <=k terms. Every grouped matrix is an available atom. The minimum total degree of a nonzero reduced polynomial function vanishing on S_w(H) is **wk+1**. Indeed that set contains every coordinate vector of Hamming support at most wk, since the corresponding matrix rank is no greater than its number of nonzero entries. Restrict a putative annihilator to every coordinate subspace of dimension at most wk, and apply finite-field interpolation to the reduced polynomial. Every coefficient whose monomial uses at most wk coordinates must vanish, and every monomial of total degree <=wk has that property. A minor of order wk+1 supplies a matching nonzero polynomial function. This degree bound is invariant under an invertible linear change of coordinates (compose and reduce exponents; the degree does not increase). The #25 degree-by-pin algorithm therefore cannot establish a separator below order wk+1 on this family, whereas the reconstruction and rank test find one remote syndrome in polynomial time. In the displayed regime wk=Theta(s)=Theta(sqrt(d)), much larger than w.

**Maximal individual-column tree separator width.** Every proper hyperplane contains strictly less than one third of the dictionary. To prove this first draw independent uniform s-by-k and k-by-s factors X,Y. For any nonzero matrix A defining a linear functional on their product, conditional on X the value <A,XY> is a linear form on Y. It is identically zero exactly when A^T X=0, an event of probability q^(-k rank(A))<=q^(-k). Otherwise its probability of being zero is 1/q. Hence

    Pr(<A,XY>=0) <= 1/q+(1-1/q)q^(-k).

The probability that both factors have full rank is at least 1-2q^(k-s)/(q-1), by a union bound on failures to extend independent vectors. Conditional on this event the product is uniform over the rank-k matrices: each rank-k matrix has exactly |GL_k(F_q)| full-rank factorizations. The lower-rank tail from (2) contributes at most 8 q^(-2(s-k)-1) relative to R_k. Consequently the fraction of the full dictionary in the hyperplane is bounded by

    [1/q+(1-1/q)q^(-k)]/[1-2q^(k-s)/(q-1)] + 8q^(-2(s-k)-1) < 1/3

under the stated guards. Any subset of at least K/3 columns therefore spans the entire matrix space. In any binary tree with the individual columns as leaves, descend along a heavier child until the leaf count is at most 2K/3. Its size is >K/3, giving a cut with both sides of size at least K/3. Both sides span V, so their linear interface intersection has dimension d. This holds even after hiding the matrix coordinates. The exact profile algorithm in the next section does not provide a small-width bound for such individual-column trees. This is not a computational lower bound; the rank recognition method is precisely an alternative efficient algorithm on the family.

**Limit of the matrix-rank method.** It cannot be claimed that all sparse dictionaries admit such a layout, even allowing arbitrary surjective linear maps to rectangular matrices. For a given d,w put t=floor((floor(sqrt(d))-1)/w)+1, and take as dictionary every vector of coordinate Hamming weight at most t. Suppose a surjective linear map T: F_q^d -> M_(a,b)(F_q) had rank(Tx)<=k for every dictionary element, with w k<min(a,b). Then k+1<=t, since ab<=d. A (k+1)-minor of T(x) is a polynomial function of total degree <=k+1; choose a minor that is not identically zero by surjectivity. The polynomial would vanish on the entire Hamming ball of radius t, contrary to the interpolation argument above (reducing exponents does not increase degree). Thus every such rank-only certificate fails. The original sparse span here is the Hamming ball of radius min(d,wt), so when wt<d an all-nonzero-coordinate vector is remote by a simpler method. The dictionary can also satisfy a strong volume slack as d grows. This exact obstruction is to the rank-only layout criterion, not to all rank-metric methods with outliers or all algorithms.

## 5. Exact weighted profile elimination on a small-interface tree

We give a different algorithm for arbitrary explicit *weighted* leaf supports. It handles ordinary additive avoidance as well as sparse syndrome distance. Let V_i be the span of the vectors in the i-th explicit local table A_i. A table associates a nonnegative integer cost to each listed vector; unlisted vectors have infinite cost. Empty tables are allowed mathematically with V_i supplied (or defined as zero if there are no columns); the proof also allows supplied larger local spans containing their supports. For example one dictionary atom h has table {0:0, a h:1 for a!=0} with minimum cost used at duplicates. A one-choice additive list has cost zero on all listed elements. A group of columns can be a leaf with its exact capped local cost table, whose construction cost must be included. If its local dimension l has q^l affordable, elementary dynamic programming by one column at a time computes that table in polynomial-in-q^l time.

Fix a rooted binary tree on the tables. At a node v let V_v be the span of its descendant local spans, let O_v be the span of all outside local spans, and put S_v=V_v intersect O_v. Its interface width is dim S_v. Let kappa be the maximum width of the supplied tree. The root has S_root=0. Bases for these spaces are computable bottom-up and top-down by ordinary field elimination. In a column-wise tree this parameter is the usual cut-rank connectivity of the represented column matroid, up to the conventional additive one in branchwidth definitions. We charge the supplied tree and its width verification; we do not assume that an optimal tree is cheaply found.

Let c_v(y) be the minimum total cost to express y as a sum of one allowed choice from each descendant table, or infinity if none. Fix a budget w and use the cap w+1. For every y in V_v consider its **boundary distance profile**

    F_(v,y)(s) = min(w+1,c_v(y+s)),   s in S_v.

Store one actual witness y for every distinct profile. There are at most (w+2)^(q^kappa) possible profiles. The crucial point is that the number of cosets V_v/S_v, potentially much larger, never has to be enumerated.

**Leaf.** Let T_i be the vectors of table i whose cost is <=w. Every non-all-capped profile has a witness among y=x-s, x in T_i and s in S_v. Enumerate these at most |T_i| q^kappa witnesses and compute the profile by table lookup. There is an all-capped profile exactly when some coset of S_v in V_v misses T_i. Compute quotient coordinates, sort the observed coset indices as base-q integers, and find a gap in [0,q^(dim V_v-dim S_v)). If one exists, lift its coordinate vector to obtain a witness. The integer bound q^dim is represented in O(d log q) bits; no scan of the quotient is required. This also handles tables whose finite costs all exceed the cap.

**Internal combination.** Let the two children be a,b. For any y_a in V_a,y_b in V_b set y=y_a+y_b. For every s in S_v the exact capped identity is

    F_(v,y)(s) = min_(s_a in S_a,s_b in S_b, s_a+s_b=s)
                  min(w+1, F_(a,y_a)(s_a)+F_(b,y_b)(s_b)).      (4)

The minimum of an empty set is w+1. To prove (4), each pair of child representations on the right supplies a parent representation. Conversely, suppose a parent representation splits into z_a in V_a and z_b in V_b with z_a+z_b=y+s. Write s_a=z_a-y_a and s_b=z_b-y_b. Then s_a+s_b=s. As s belongs to O_v, s_a lies in V_a intersect (V_b+O_v)=S_a, and similarly s_b lies in S_b. Thus every parent representation is accounted for. Nonnegative costs justify capping at every operation. Formula (4) depends only on the two profiles, not the witnesses used. Enumerate all pairs of stored child profiles, compute (4), and keep one witness y_a+y_b for each resulting profile. Every y in V_v has a decomposition into child vectors because V_v=V_a+V_b, so by induction **every** actual profile is represented, and each stored witness realizes its profile.

**Root.** S_root=0, so the profiles are precisely the attainable capped minimum costs. If V_root is a proper subspace of the ambient space, a vector outside it is immediately missing. Otherwise a stored root profile (w+1) is a syndrome with distance/cost greater than w. If there is no such profile, the support with budget w covers V_root. To compute the exact covering radius for a spanning ordinary dictionary, use w=d (or the rank of the whole dictionary); a basis represents every syndrome with at most d atoms. For additive avoidance take w=0 and tables of cost zero; the root decides and constructs a missing point even without a volume promise. Formula (4) also works for general nonnegative integer local costs, subject to the chosen cap.

The cost is polynomial in the total explicit leaf table size, d, q and tree size, times q^O(kappa)(w+2)^(2 q^kappa). Indeed each leaf uses at most |T_i|q^kappa candidate witnesses plus one gap search, each internal node enumerates at most (w+2)^(2q^kappa) profile pairs, and each combine can precompute at most q^(2kappa) valid interface pairs. Store explicit capped vectors and actual witness coordinates; deterministic sorting/deduplication has polynomial overhead. The bit cost of integer caps and vector coordinates is polynomial in their explicit lengths. This is not an efficient generic high-width algorithm: the term exponential in q^kappa must be paid. With fixed small interface width it can process huge sumsets or highly redundant dictionaries without enumerating the ambient space. At the high-rate scales of #25 an input-dependent branch with q^kappa log(w+2)=O(log K) is polynomial in the explicit dictionary size. The rank-ball family in section 4 demonstrates that the volume slack does not force this interface bound on individual-column trees.

## 6. Exact verification and scope

Run `python -B check_structured_remote.py` with assertions enabled (standard library). The finite verifier is not an implementation of any asymptotic curve algorithm. It checks:

* 84 exact rational/integer parameter triples for the analytic overlap and hyperplane guards;
* five complete balls in unknown linear coordinates, including a rank-two ball, with 1,774 individual rank comparisons and exact sparse-support exclusion of the constructed full-rank syndrome;
* a retained binary 2-by-2 negative overlap example (all nonzero shifts tie), plus a noisy binary 5-by-5 example with four deletions, two arbitrary high-rank additions and complete recovery of its 961 rank-one differences; the rank-metric outlier witness is independently checked against the actual two-atom span;
* five additional rank-metric outlier controls, 183 pairwise rank-distance comparisons, one further unknown-affine noisy fixture (961 supplied atoms, two recentered outliers, six centers and seven probes, exact one-atom exclusion), 400 projective hyperplanes in the F_7 2-by-2 rank-one example, and all 20,160 invertible binary four-dimensional reshapes for the Hamming-ball obstruction;
* three binary random-thinning toy runs, outside the asymptotic q,s guard; all three recover the full rank-one set;
* 48 random weighted/sparse tree instances, 6,668 stored node profiles independently checked against full small-span computations; a 65-dimensional colored chain with 3^64 tuple choices and width 2; and 7,764 nonzero atoms in overlapping local dictionaries of dimension 59 over F_3, for which the exact profile algorithm returns covering radius 12 without enumerating 3^59 syndromes.

The finite fixtures do not reproduce the huge asymptotic rank-ball parameters, nor do they demonstrate that arbitrary projected dictionaries satisfy either new structure. The first attempted maximum-overlap verification failed on the small binary exception, which lies outside the theorem guard; the failure is now retained in the transcript as a negative control. All completed assertions and the final transcript are reproducible. The large structural examples are supported by proofs rather than extrapolation of finite timings.

Related literature: Volk, *Improved Algorithms for the Remote Point Problem*, https://arxiv.org/abs/2609.25765, gives the classical syndrome formulation and algorithms for code rates bounded away from one; its hypotheses do not cover a general high-rate dictionary here. Dvir, Golovnev and Weinstein, *Static Data Structure Lower Bounds Imply Rigidity*, https://golovnev.org/papers/ds.pdf, identify sparse-span-evasive sets with sparse linear factorizations, explaining why an explicit universal remote set is substantially stronger than its counting existence proof. For the bilinear forms association scheme and its intersection numbers see *Strengthening the Balanced Set Condition for the Distance-Regular Graph of the Bilinear Forms*, https://doi.org/10.1007/s00026-026-00857-3; related rank-distance isometry results appear at https://publications.hse.ru/en/articles/1103759761. The linearized-polynomial rank-metric construction is classical (Delsarte/Gabidulin/Roth; see https://doi.org/10.3934/amc.2022034). Matroid pathwidth/trellis width and algorithms for supplied or bounded-width decompositions are classical structural tools; see Kashyap https://doi.org/10.1137/070691152 and https://doi.org/10.1109/TIT.2017.2740283. We do not claim novelty for these ingredients or for all possible covering-radius dynamic programs. The overlap bounds, explicit reconstruction, noisy recognition conditions, and profile identity are proved here in the precise forms used. Independent literature-priority and proof review remain outstanding.

**Status:** two exact positive algorithms for structured abstract dictionaries. They work on substantial high-rate examples that invalidate the earlier low-degree or small-interface assumptions. They do not solve general high-rate remote-set construction, do not establish a universal structural disjunction, and do not establish an improved end-to-end exponent for the preceding project. No real key recovery or deployed cryptographic experiments are included.
