# Overlapping block-ball recognition and flag-lifted pseudo-arc certificates

27 November 2026. Abstract finite-geometry follow-up to PRs [#28](https://github.com/davidiach/elliptic-relation-counting/pull/28), [#31](https://github.com/davidiach/elliptic-relation-counting/pull/31), [#37](https://github.com/davidiach/elliptic-relation-counting/pull/37), and [#40](https://github.com/davidiach/elliptic-relation-counting/pull/40), reviewed at main `13d2803e9f56e70d83d036a9ee426eba0d1f4221` and #40 head `763c581afd1b6340da7785bdbc5853d0a45285f4`. The complementary simultaneous-shadow work #36 and affine-core peeling #39 were checked at the level of their stated results. This note studies recognition and certification of explicitly supplied finite-field vector configurations. It does not claim a universal sparse-span or additive-avoidance algorithm, a new discrete-logarithm complexity bound, or an implementation on cryptographic inputs. All experiments are synthetic mathematical fixtures. Independent proof and priority review remain outstanding.

## 1. Recovering overlapping groups from a hidden block ball

Let q be prime, let V_1,...,V_L have common dimension a>=1 over F_q, and let W=direct_sum V_i. Put Q=q^a. The block weight bw(x) counts the nonzero components of x in W. For integer t>=1 let

    B_t={x in W:bw(x)<=t},   K=|B_t|=sum_{i=0}^t binom(L,i)(Q-1)^i.

Assume L>=2t+1 and let A:W -> F_q^d be linear and injective on the sum of any at most 4t of the formal blocks (or on all W if 4t>=L). No global injectivity is assumed. The hidden linear dictionary is S=A(B_t), consisting of K distinct vectors; an affine version is c+S for an unknown c. The complete vectors, including all scalar multiples, are supplied explicitly; replacing projective inputs by complete scalar lines costs up to a factor q and must be charged. The block labels, A, and the original coordinates are not supplied.

This model is substantially denser than the one-bridge-per-pair model in #40: for t>=2, every line joining core points in distinct blocks is contained in the two-block ball. Rich-line connected components can merge all blocks. Here the signal comes from pairwise overlaps rather than line density.

### 1.1 Exact intersection and gap identities

For v in W of weight j let C_j=|B_t intersect (B_t+v)|. This is independent of the particular v because each coordinate is an additive group of cardinality Q and the ball only tests zero versus nonzero. Define V(n,Q,h)=sum_{i=0}^{min(n,h)} binom(n,i)(Q-1)^i for h>=0, and zero for h<0. Then

    C_j = sum_{r+s+u=j} [j!/(r!s!u!)] (Q-2)^u
                    V(L-j,Q,min(t-j+r,t-j+s)).                (1)

Indeed, among the j changed coordinates, categorize x_i as zero, the value that cancels the shift, or any other value. If the respective counts are r,s,u, the weights within those coordinates in x and x+v are j-r and j-s. Outside the support the weights agree. Counting the allowable outside choices proves (1), using 0^0=1. In particular C_0=K and C_j=0 for j>2t.

For 0<=j<L, compare shifts of weights j and j+1 that differ only in one additional coordinate. Fix all other coordinates of x. To meet either ball constraint, the allowed set in the last coordinate is either all Q values, one specific value, or empty. The intersection size changes when the shift in that last coordinate goes from zero to nonzero if and only if both allowed sets were singletons. In that case it decreases by one. Both singletons occur exactly when the two weights in the other L-1 coordinates are each t. Thus

    C_j-C_(j+1) = sum_{r=0}^{floor(j/2)}
          [j!/(r!r!(j-2r)!)] (Q-2)^(j-2r)
          binom(L-1-j,t-j+r)(Q-1)^(t-j+r),                   (2)

where a binomial term outside its nonnegative integer range is omitted. Formula (2) is an exact shell-intersection count, in particular nonnegative. If Q>=3 and L>=2t+1, for each 0<=j<=2t take r=floor(j/2). Then 0<=t-j+r<=L-1-j, and its contribution is strictly positive. Consequently

    C_0>C_1>...>C_(2t)>C_(2t+1)=0.                           (3)

No probabilistic or asymptotic approximation is used. For instance

    C_1=Q E+Q^2 V(L-2,Q,t-2),
    C_2=2 E+Q^2 V(L-2,Q,t-2),
    E=binom(L-2,t-1)(Q-1)^(t-1).

Their gap is (Q-2)E. Q=2 is exceptional: (2) gives C_(2j-1)=C_(2j) whenever both indices are in range, and the positive even-to-odd gaps persist. Section 3 handles that ambiguity separately.

These are elementary Hamming-ball intersection identities; see Junnila, Laihonen, Lehtila and Devraj, [The Size of the Intersection of q-ary Hamming Balls](https://arxiv.org/abs/2606.09158), for broader intersection theory. No novelty is claimed for the intersection count in isolation.

### 1.2 Automatic recognition and decomposition

For a shift h in the observed ambient space compute D_S(h)=|S intersect (S+h)| by tabulating all ordered pair differences of the explicit vectors. A nonzero observed difference is A(v) for some v of block weight <=2t. Such a v is unique: two possible preimages have difference supported on at most 4t blocks. Moreover all differences of elements in B_t with this shift are precisely the pairs counted by C_j, so D_S(A(v))=C_j for j=bw(v). Every v of weight <=2t is a difference of two ball points. The identity remains true for an affine translate, since translating a set does not change its difference counts.

**Theorem 1.** For Q>=3, the complete affine dictionary determines the individual image subspaces A(V_i) uniquely, up to order, under the stated model. They and all block decompositions relative to the center can be recovered in deterministic time polynomial in the explicit input length and d,log q. Neither L,a,t nor the hidden map need be supplied in the noiseless case.

Proof. The set R of nonzero shifts attaining the maximum nonzero overlap is exactly union_i (A(V_i) minus {0}) by (3). For a nonzero x in R, its group is

    {y in R: x+y in R union {0}}.

If x,y come from the same formal block this follows from closure; if they come from different blocks, equality to a core vector or zero contradicts local injectivity on at most three blocks (4t>=4). Thus the displayed relation splits R into its complete nonzero subspaces. Compute their bases by elimination. Their ranks determine a, and their number determines L. The unique t in [1,floor((L-1)/2)] with V(L,Q,t)=K is determined by the explicit cardinality; a failure is reported honestly. All these computations use polynomial work in the observed list and its difference table. This also proves uniqueness for other descriptions satisfying the same hypotheses.

For a point x of the centered dictionary, if x is zero stop. If x belongs to R, take it as its sole remaining component. Otherwise, among all v in R choose one maximizing D_S(x-v), where an unobserved difference has count zero. In true formal coordinates, subtracting a core vector can decrease the weight by one precisely when that core vector is the full nonzero component in its block; other choices preserve or increase the weight. By (3), all maximizing choices reduce the weight. The current weight is at most t and the comparison shifts have weight at most t+1<=2t, so the local-injectivity argument applies. Repeat at most t times. This proves exact decomposition without enumerating support subsets. Each step uses |R| table lookups and vector operations. In the affine case the center can be recovered first by section 2.2. QED.

On arbitrary data, the difference statistic merely proposes groups. One may check that the proposed core components are complete subspaces by rank and size, and verify every returned vector sum explicitly. If the algorithm fails a guard it has not proved absence of a useful factorization. Section 4 gives a separate algebraic certificate of the local-independence property that does not assume the overlap inference was sound. If independence on 2t blocks is certified, the decompositions of all supplied distinct points are valid and their number equals K, their equality with the full block ball follows by cardinality. No need to infer it from an unverified statistical model.

## 2. Robust recognition, thinning, and affine centers

### 2.1 Arbitrary edits and independent thinning

Changing a set by at most e insertions/deletions (symmetric-difference size at most e) changes its overlap count at every shift by at most 2e: each changed point participates in at most one ordered pair at either end for a fixed shift. Suppose the original parameters L,Q,t are known in the noisy case. Put J=max(1,t-1) and Delta_j=C_j-C_(j+1) for j=1,...,J. If

    4e < min_{1<=j<=J} Delta_j,

thresholding the observed counts at (C_1+C_2)/2 recovers exactly R. Moreover each of the necessary nonzero-weight comparisons for the greedy reduction is correct, since the counts of weight j and every higher weight differ by at least Delta_j. Direct detection of a zero residual or a residual in R handles the final step. Every retained inlier factors correctly in its original blocks. For an arbitrary outlier, verification of the returned decomposition is essential; an actual point outside the modeled ball cannot pass the t-component reconstruction. The theorem does not guarantee that the heuristic factors every point of an arbitrary nonconforming input. The explicit bound concerns edits, not a constant-fraction corruption tolerance uniformly for all parameters.

There is a stronger sampling result. Retain each original vector independently with probability theta in (0,1], then make at most e arbitrary edits to the sample, even adversarially after observing it. Write S' for the supplied result. If Q>=3 and

    8e <= theta^2 min_{1<=j<=J} Delta_j,

then thresholding at theta^2(C_1+C_2)/2 and performing the same greedy reduction recovers the groups and correctly factors all centered retained inliers with probability at least

    1 - 6 J q^d exp(-theta^2 min_{1<=j<=J}(Delta_j^2/C_j)/640).    (4)

Interpret a negative probability lower bound as vacuous. The retention probability is supplied for this theorem; parameter estimation is not part of its guarantee.

Proof. Fix a nonzero ambient shift h with a positive full count C. The underlying undirected pairs counted by D_S(h) form a graph of maximum degree two, since each point has at most one neighbor in either shift direction. For odd q the directed pairs correspond to distinct undirected edges, of weight one. For q=2 each undirected edge accounts for two ordered pairs, of weight two; the graph is a matching in this case. Color the undirected edges with at most three matchings. Within any matching the indicators that both endpoints survive are independent Bernoulli variables of mean theta^2. Each contribution has weight at most two. For a matching with expected weighted count nu, the standard Bernstein bound (variance <=2nu, bound on each centered term <=2) gives two-sided failure at deviation u at most 2 exp(-u^2/(4nu+4u/3)). Union over at most three matchings, using u=A0/3, yields for the total sampled count with mean mu=theta^2 C

    Pr(|D_sample(h)-mu| >= A0) <= 6 exp(-A0^2/(36mu+4A0)).

For a fixed j let T_j=theta^2(C_j+C_(j+1))/2. After edits, if C>=C_j and observed count is <=T_j, the pre-edit count deviated below mu by at least A0=mu-T_j-2e. At the smallest possible mu=theta^2 C_j this is at least theta^2 Delta_j/4. The expression (mu-U)^2/(40mu-4U) with U=T_j+2e is increasing for mu>=U>=0; hence the exponent is at least theta^2 Delta_j^2/(640 C_j). If C<=C_(j+1) and observed count is >T_j, the analogous upper-tail deviation is at least theta^2 Delta_j/4 and its denominator is at most 40theta^2 C_j. Shifts with full count zero have pre-edit count zero and cannot pass the threshold under the edit guard. Union over at most q^d-1 nonzero shifts and J thresholds gives (4). This simultaneously separates all weight classes needed in the proof of Theorem 1. The proof is analytical; no algorithm checks all q^d shifts. QED.

**Thin-sample scale.** For fixed t and L asymptotic to Q, as Q tends to infinity, (1) and (2) give C_j and Delta_j both Theta_t(Q^(2t-j)) for j<=2t. To see this, a term of (1) with counts r,s,u and outside bound h=t-j+min(r,s) has order at most Q^(u+2h)=Q^(2t-j-|r-s|), and the r=s=floor(j/2) term in (2) attains the top exponent. In particular for t>=2 the minimum ratio in (4) has order Q^(t+1); for t=1 it has order Q. If d log q grows polylogarithmically in Q, t is fixed, and theta=Q^(-beta) with fixed beta<(t+1)/2 (beta<1/2 when t=1), the probability bound tends to one with e=0. This is substantial thinning, not merely the deterministic small-edit bound. One family meeting the local-injectivity premise has L=Q+1, blocks obtained by reducing a Vandermonde arrangement over F_Q to F_q, and d=sa with Q=q^a and 4t<=s=polylog(Q); section 4 certifies such arrangements. The original ball size is Theta_t(Q^(2t)); its retained sample size is Theta_t(theta Q^(2t)) with high probability. Thus the full ball size, if needed for affine recovery below, is polynomial in the retained input size in this regime. This is a structured example family, not an inference that arbitrary explicit dictionaries have such a structure.

### 2.2 Unknown affine centers, including a deleted center

Suppose the groups have been recovered and let B denote the image of their full block ball before translation. Enumerate B only under an explicit cap on K; its enumeration has at most K distinct outputs and uses polynomial work in K and the group basis descriptions. In the noiseless and bounded-edit setting K is comparable to the supplied input length. In the thin-sample setting charge K separately, or require K<=|S'|^C for a fixed C. This avoids materializing a huge unobserved ball on a tiny input.

For a full or bounded-edit affine dictionary, a candidate center c' can be tested by exact overlap with c'+B. The true center c achieves symmetric difference at most e. Every distinct translate of B intersects B in at most C_1 in the nonbinary case: a shift with no representation as a difference gives zero overlap, and every other nonzero shift is governed by (3). Thus the symmetric difference between two translates is at least 2(K-C_1). If e<K-C_1, at most one center can pass the e-edit guard. Enumerate cross-differences x-b for x in S', b in B and count multiplicities; the count is exactly the intersection size for that candidate center. The true center occurs among the cross-differences when at least one inlier survives. This costs O(|S'|K) pair operations with explicit sorting/hashing; there is no ambient scan. The verifier uses an equivalent small-anchor candidate method with at most e+1 anchors and explicit intersection checks. The bound is polynomial in K and |S'|. The actual center need not be present in S'.

For independently thinned data, the same cross-difference count maximization recovers the original center with high probability when the original ball is affordable to enumerate. Indeed for any other center, at least K-C_1 of the original ball points are absent from its translated ball. With probability at least 1-exp(-theta(K-C_1)/8), at least theta(K-C_1)/2 of those points are retained (multiplicative Chernoff). Union over q^d-1 competing centers. If 4e<theta(K-C_1), the true center strictly wins the overlap comparison even after e arbitrary edits, since one edit changes the difference of the two overlap scores by at most two (in fact one suffices). The failure probability adds at most q^d exp(-theta(K-C_1)/8) to (4). The same argument ensures at least one retained inlier on the stated success event. Decompose each centered input point only after selecting this center, verifying every reconstruction. No affine-invariance of a sparse linear span is asserted: a center is a separate quantity and must not be silently dropped in other applications.

## 3. The binary ambiguity is geometric, not merely a failure of counts

For Q=2, a=1. Assume t>=2, L>=2t+1 and the same 4t local injectivity. Embed F_2^L into the even-weight space of F_2^(L+1) by x -> (x,parity(x)). The formal weight-one and weight-two vectors correspond precisely to the edges {u,v} of a complete graph on L+1 vertices, represented by e_u+e_v. By (2), C_1=C_2>C_3. Hence the maximum nonzero shift set R in the complete dictionary is precisely this edge set after the hidden linear map.

Two distinct observed edge vectors have sum in R exactly when the corresponding edges meet. In the formal even-weight space this is elementary symmetric difference. Any erroneous equality after projection would give a nontrivial dependence involving at most six original coordinate blocks, impossible by 4t>=8. Thus the observed adjacency graph on R is the line graph of a complete graph, with no vertex labels supplied.

Recover a star without enumerating cliques: pick adjacent edges e,f. Their common neighbors, together with e,f, form the star of their common endpoint plus the edge between their two other endpoints. The extra edge has degree exactly two inside that induced graph; each star edge has degree at least L-1>=4. Delete the degree-two vertex. The resulting L edge vectors are a star and hence a coordinate basis in the formal even-weight space. In the projected space they may have global dependencies, as permitted. All of this costs polynomial work in |R| and the ambient dimension.

With even t, the original ball corresponds to all even-weight vectors of length L+1 and weight <=t, independent of which vertex is the extra parity coordinate. Thus any recovered star gives a valid block-ball description with the same affine center. With odd t, use the odd-parity embedding (x,1+parity(x)); the ball corresponds to all odd-weight vectors of length L+1 and weight <=t. Switching the extra vertex changes the coordinate description by translation by the edge joining the old and new extra vertex. Consequently an arbitrary recovered star gives an affine description, possibly with a different center. The original star is not in general identifiable from an affine dictionary, and a unique linear decomposition is not claimed for odd radii without accounting for the center.

To recover a center in the chosen star, use the finite cross-difference search of section 2.2. The new star ball has exactly K distinct images: a relation involving at most 2t new star vectors involves at most 2t+1 original blocks, and 2t+1<=4t. It is a translate of the original ball after projection. Thus the same center uniqueness proof applies with maximum overlap C_1=C_2. The method also works under 4e<C_2-C_3 edits for a supplied binary L,t: threshold at (C_2+C_3)/2 recovers R, then use the same star and center routine. The verifier checks complete binary dictionaries and a radius-two fixture with one deletion. In that fixture the midpoint recovers all 45 edge shifts, of which 28 attain the observed maximum. The t=1 case is an affine simplex model and is not covered by the line-graph proof here (the local injectivity guard there would not justify every adjacency test).

This explains why enforcing a strict maximum at weight one would wrongly reject a natural class. The ambiguity is inherent to the Hamming geometry. It does not imply impossibility of recognition up to a valid alternative layout.

## 4. Certifying all maximal-group independences without enumerating them

In sections 1-3 local injectivity is a model assumption. On arbitrary data, simply observing high overlaps is not proof that every small set of inferred groups is independent. Checking all sets of s groups can cost binom(L,s) rank tests. We now give a separate polynomial-time certificate for an important broader class of subspace configurations. It includes extension-field rational normal curves (field-reduced arcs), commutative operator-ring arcs, and noncommutative flag-lifted arcs with arbitrary off-diagonal perturbations. No hidden extension field or flag basis is needed in the certificate.

Let U_1,...,U_g be explicit a-dimensional subspaces of F_q^(sa) for s>=2, g>=s. We seek to certify that any s distinct groups form a direct sum (a subspace arc, also called a pseudo-arc; equivalently a block MDS condition). If the ambient has redundant coordinates, first restrict to the actual span of the groups and require its dimension to be sa. The test is sufficient on arbitrary inputs and complete for the flag-lifted class below. It is not a decision procedure for the general pseudo-arc condition.

### 4.1 Normalize an arbitrary frame

Use the first s groups as coordinate blocks, checking their combined matrix invertible. If g=s this suffices. If g>=s+1, project the next group onto each reference block and require each projection is an isomorphism. By separate changes of basis in each reference block, make this anchor the diagonal group {(u,...,u)}. When g=s+1 this already certifies any s groups independent, since omitting any one reference coordinate, the anchor's projection to it is an isomorphism. If any of these guards fails, reject this branch.

For each further nonreference group, also require all its coordinate projections to be isomorphisms. In normalized coordinates each nonreference group is now the graph

    {(u,T_(j,2)u,...,T_(j,s)u): u in F_q^a},

where T_(j,i) is invertible. Include the anchor j=0 with all T_(0,i)=I, and define Z_(j,i)=T_(j,i)^(-1), Z_(j,1)=I. This is all computed by ordinary elimination. For a block basis matrix C with blocks C_i, T_i=C_i C_1^(-1).

Pick any further group j=* different from the anchor, which exists when g>=s+2. Define parameter operators

    beta_j = Z_(j,2),
    alpha_2=0,
    alpha_i=(I-Z_(*,i))^(-1)(beta_*-Z_(*,i)), i=3,...,s.       (5)

Reject if an inverse does not exist. We will interpret reference group 1 as infinity, reference group i>=2 as alpha_i, and the nonreference groups as beta_j. The order of multiplication in (5) is important in a noncommutative algebra.

### 4.2 The exact flag certificate

Let A0 be the unital F_q-algebra of matrices generated by all the graph operators Z_(j,i). Let J0 be the two-sided ideal in A0 generated by commutators of pairs of generators. Require the following three tests:

1. J0^a=0.
2. Every difference of two distinct-index members of the list of finite parameter operators (alpha_2,...,alpha_s,beta_0,...,beta_(g-s-1)) is invertible, including differences from the anchor beta_0=I.
3. For every nonreference j and i>=2, the residual

       R_(j,i)=Z_(j,i)(I-alpha_i)-(beta_j-alpha_i)

   satisfies R_(j,i)^a=0.

All tests use exact matrix operations. In particular the third test does NOT require R_(j,i)=0. Noncommutative nilpotent off-diagonal parts are allowed. Every operator in (5) lies in A0, since the inverse of an invertible matrix in a unital algebra is a polynomial in it by Cayley-Hamilton.

**Theorem 2.** If these tests pass, any s distinct groups form a direct sum. The certificate can be discovered and checked deterministically in time polynomial in g,s,a,log q, with absolute exponents. There is no enumeration of binom(g,s) subsets or construction of an algebraic closure.

Proof. First, J0^a=0 is equivalent to simultaneous triangularizability of A0 over the algebraic closure. The forward implication has an elementary proof: after extending scalars, a nilpotent ideal J has a nonzero common annihilator on every nonzero module. For example, choose the largest k with J^k U nonzero for the current module U and take a nonzero vector from J^k U. The common annihilator is A0-invariant since J is an ideal. On it, the action factors through the commutative algebra A0/J0; a commuting family of matrices over an algebraically closed field has a common eigenvector. Take that invariant line and induct on the quotient to obtain a complete invariant flag. Conversely, if all matrices are upper triangular in a common basis, all their commutators and the ideal generated by them are strictly upper triangular in that basis; their a-fold products vanish. This equivalence is classical (see also [ECCC TR13-103](https://eccc.weizmann.ac.il/report/2013/103/), Lemma 20 and Corollary 21).

In such a common flag, all the graph and parameter matrices are triangular. Invertibility of parameter differences guarantees that in every diagonal layer the corresponding scalar parameters are distinct. Nilpotence of R_(j,i) guarantees its diagonal entries are zero. In any fixed layer, denote the scalar parameters by alpha_i,beta_j. The graph vector for group j is then

    (1, (1-alpha_2)/(beta_j-alpha_2),...,(1-alpha_s)/(beta_j-alpha_s)).

Here alpha_2=0 and beta_0=1. Multiply it by the nonzero product prod_{i=2}^s(beta_j-alpha_i). The resulting coordinates are the values at beta_j of the degree <=s-1 polynomials

    L_1(z)=prod_{i=2}^s(z-alpha_i),
    L_i(z)=(1-alpha_i) prod_{2<=h<=s,h!=i}(z-alpha_h),  i>=2.

These form a basis of the polynomial space, since evaluation at the s-1 finite reference nodes isolates L_i for i>=2 and the leading coefficient isolates L_1. Their homogeneous degree-(s-1) values at infinity give the first reference coordinate. Thus every selected scalar set of s projective columns is nonsingular by the Vandermonde determinant (or interpolation at s projective nodes). This holds in each diagonal layer. For any selection of s original groups, reorder the rows and columns of its block matrix by flag coordinate, obtaining an upper block triangular matrix with those s-by-s scalar matrices as diagonal blocks. Its determinant is nonzero over the closure, hence over F_q. QED.

**Resource accounting.** Flatten matrices as a^2-dimensional vectors. Construct A0 by adjoining I and the generators, repeatedly closing their linear span under right multiplication by generators until stable. At most a^2 basis additions are needed. Similarly form the span of generator commutators and close it under left/right multiplication by generators to obtain J0. Repeated products of its basis and elimination compute J0^2,...,J0^a; each intermediate span has dimension at most a^2. The number of graph generators is at most gs. Computing pairwise parameter differences uses at most binom(g-1,2) invertibility tests; residual powers use repeated squaring. Ordinary field arithmetic and elimination therefore give a polynomial-time test with no exponential dependence on s. Neither explicitly triangularizing matrices nor factoring polynomials over an extension field is necessary. If A0 happens to be commutative, optional Frobenius diagnostics distinguish a field from other commutative algebras: injective Frobenius x->x^q means the finite algebra is reduced, hence a product of finite fields; its fixed-space dimension is the number of factors. Dimension one of that fixed space and injectivity together certify a field. These classical facts are not needed for Theorem 2.

### 4.3 Completeness for hidden flag-lifted rational normal curves

The theorem is not merely a sufficient test with no useful completeness guarantee. Suppose that after an unknown ambient linear transformation and group-basis choices (over the algebraic closure if needed), the g groups are represented by s-block columns of a-by-a upper triangular matrices. In each of the a diagonal layers, the corresponding g scalar projective points lie on some rational normal curve in projective (s-1)-space, with distinct parameters within that layer. Each layer may use its own curve and its own parameter assignment to the group indices. The strictly upper triangular entries of all blocks may be arbitrary. Call these flag-lifted curves.

**Claim.** Every such arrangement passes Theorem 2, using any s+2 distinct groups for the reference frame when g>=s+2. To verify this, first use the latent flag basis. Any s group block matrix is invertible by the same layerwise Vandermonde argument; after reordering coordinates by flag layer, its inverse has blocks that are triangular in the same flag. Projections of the anchor onto all reference blocks have nonzero scalar diagonal entries, so they too are invertible. The frame normalization and graph maps are therefore upper triangular in one common flag. Arbitrary bases of the input groups cancel during normalization except for a common conjugation on the anchor coordinate space: if group j's basis is postmultiplied by S_j, its normalized graph maps become conjugate by S_0, the anchor's basis change. An ambient change cancels with the reference frame inverse. Thus the actual computed graph algebra is simultaneously triangularizable over the closure.

For clarity, the scalar normalization in each layer is elementary. On the rational normal curve move the first reference parameter to infinity, the second to zero, and the anchor to one by a fractional linear change of parameter. Such a change acts invertibly on homogeneous degree-(s-1) polynomials. Normalize the coordinates using the s reference points and anchor as above. Then the projective coordinates at any remaining finite parameter beta are exactly the ratios from the L_i polynomials in the preceding proof: each reference finite point isolates one basis polynomial and the point at infinity isolates the leading coefficient, and the anchor sets the relative scalings. The ratios satisfy z_i(beta)=(beta-alpha_i)/(1-alpha_i). At the extra reference beta_*, z_i !=1 for i>=2 because beta_*!=1. Formula (5) thus recovers the correct scalar reference parameters on each diagonal layer. Differences are nonzero on every diagonal, and every ratio residual has zero diagonal, so the matrix tests pass. This proves completeness, even for arbitrary off-diagonal perturbations of this model. It does not assert completeness for every block MDS code or every pseudo-arc.

An extension-field rational normal curve, viewed over F_q by representing F_(q^a) multiplications as a-by-a matrices, is a special case: over the closure the multiplication algebra diagonalizes with one layer per field embedding. The certificate recovers enough of its common operator structure without assuming the extension degree or an extension-field representation is supplied. Commutative product rings with unit-separated parameters and nonreduced rings (e.g. dual numbers) also fit. More importantly, arbitrary noncommuting triangular lifts fit even when the residuals R_(j,i) are nonzero matrices; merely testing commutation or exact scalar ratio equalities would reject such examples.

Related finite-geometry background: Pavese and Santonastaso, [On pseudo-arcs from normal rational curve and additive MDS codes](https://arxiv.org/abs/2602.23130), 2026, discusses both field-reduced and genuinely non-Desarguesian constructions. The present result is an explicit polynomial certificate for the stated flag-lifted class, not a classification of all pseudo-arcs or an independent novelty claim about constructing additive MDS codes. Simultaneous triangularization and Vandermonde interpolation are classical; the combined normalization, nilpotent-ratio test and overlap recognition are the contributions developed here.

### 4.4 Combining discovery and independent verification

When d=sa and an input really is a hidden complete block ball of a projected flag-lifted arrangement with s>=4t, Theorem 1 recovers the blocks and all representations, and Theorem 2 validates their maximal-subset independence from the recovered bases. No enumeration of all s-group combinations is necessary. The same certificate can be used on groups proposed by other recognition methods; its soundness is independent of the statistics that proposed them. For example, the output records explicit group bases and component decompositions and verifies each vector sum. If the verified arrangement gives independence on at least 2t groups and the complete cardinality guard holds, the input is exactly the proposed block ball. For sparse samples this is only a certified representation of the retained points, not a certification of unobserved points or an assumption that all points of an arbitrary input belong to the model.

The regime g large relative to s is real: take Q=q^a, g=Q+1 and distinct extended-field projective parameters, with 4t<=s<=g. Then the formal block domain has dimension ga, its image dimension sa and its kernel dimension (g-s)a. For fixed t and s=polylog(Q), the ratio of these dimensions tends to zero while the input ball has Theta_t(Q^(2t)) vectors. The algebraic test has polynomial cost in g,s,a and log q, regardless of the enormous number of maximal-group subsets when s grows. This is a high-overlap structured family, not a universal structural dichotomy.

## 5. Reproduction, caveats and open frontier

Run `python -B verify.py` with Python 3.10+ and NumPy. The recorded run used NumPy 2.4.4. The certificate uses Python integers in NumPy arrays for exact modular multiplication and elimination. Binary differences in the finite fixtures use packed XOR. The final probability examples display rounded logarithms of analytical bounds. `verification.txt` records the complete output.

The verifier checks 41 overlap counts, 28 strict gaps and two positive binary ties. Across 26 hidden arc layouts it independently checks the ranks of 2,212 maximal subsets, including eight noncommutative layouts with 24 nonzero nilpotent ratio residuals. Alongside the original 21 negative controls, separate fixtures check rejection of a singular anchor, a nonnilpotent commutator ideal and a failed ratio residual. Over each of the primes 2^61-1 and 2^89-1 it rejects an opposite pair of columns and certifies an independent five-group arc.

Across five complete configurations it uses the inferred radii for 2,304 decompositions and checks all 341,290 shift counts against their formal block supports. Exchanging counts between shifts is detected even when the histogram is preserved. A separate fixture infers radius three and decomposes a point that needs all three blocks. The bounded-edit trials verify 3,201 inlier decompositions, including a radius-one example with gap five and an edit budget of one. The original 15 affine trials cover 2,052 input points and nine binary star layouts. An additional binary trial recovers the full ball after deleting one point.

The four centered thinning controls succeed at retention 3/4 and 1/2 and fail at 1/4 and 1/8. In a separate translated 3/4 sample, the planted center is absent. Cross-difference maximization recovers it with a score of 143 against a runner-up of 67, and all 143 retained points decompose after centering. A singleton sample exercises the tied-center rejection. These seeded fixtures illustrate finite outcomes. Three of the four larger analytical examples have nonvacuous probability lower bounds, and their balls are never enumerated. The frame-only case with s+1 groups and the projected binary radius-three case illustrate the separate frame and local-injectivity hypotheses.

**Frontier.** We now have an automatic recognition theorem for complete and suitably thinned block balls, including overlapping mixed blocks, and an independent polynomial-time maximal-group independence certificate beyond a single extension-field model. The certificate is complete for flag-lifted curves but not for arbitrary pseudo-arcs. The unresolved issue is to derive useful groups for arbitrary explicit configurations that do not contain the needed overlap information, and to certify broader group arrangements without an exponential search. No universal structural dichotomy, general avoidance routine, or new end-to-end cryptographic complexity bound is established. Do not change the project's main theorem or its assumptions on the basis of this note. Independent novelty and proof review remain outstanding.
