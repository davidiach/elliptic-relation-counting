# Hidden determinantal geometry from thin samples and affine rank-ball avoidance

27 November 2026. Follow-up to elliptic-relation-counting PRs [#25](https://github.com/davidiach/elliptic-relation-counting/pull/25), [#28](https://github.com/davidiach/elliptic-relation-counting/pull/28) and [#30](https://github.com/davidiach/elliptic-relation-counting/pull/30). Main inspected at `f38f02b71160c0859e012339fb9cdeabc2dcc402`; PR #28 at `e6677bbc48dac6ade5f63fd8cbd0e4ce94368e82`; PR #30 at `60309a8e34c8dc7ba00e26d7acea82716d658481`. This note concerns abstract finite-field matrix geometry, additive-set separation, and sparse-span avoidance. No general avoidance algorithm, improved end-to-end ECDLP bound, or cryptographic-scale experiment is claimed. The prior complete-ball overlap method and hidden Hankel frame method are not claimed new here. Independent mathematical and literature-priority reviews remain necessary.

## 1. Reconstruction theorem and the precise sample condition

Let F=F_q, where q is an odd prime, a,b>=2, d=ab and 1<=k<min(a,b). Let t=k+1. Assume q>t; when reconstructing a rank-one point via the derived Lie algebra, additionally assume q divides neither a nor b. In the following theorem the unknown affine map has the form

    X -> c + T vec(X),    T in GL_d(F), c in F^d.

Let the explicit sample A be contained in the image of the rank-at-most-k determinantal variety. Put

    D = binom(d+t,t),
    B = binom(a,t) binom(b,t),
    E = D-B.

The sample-completeness hypothesis is that the nullspace of its evaluation matrix for all monomials of degree <=t has dimension B. This is directly checkable; it does not require the hidden coordinates. Under the containment promise, this means that the sample has exactly the degree-at-most-t equations of the entire affine image. On arbitrary inputs, passing this dimension test alone is not a proof of containment; validate the recovered coordinate system and all needed ranks before returning a certificate.

**Theorem 1.** Under these conditions, there is a deterministic reconstruction algorithm that recovers the affine center c and an invertible linear coordinate system S in which every vector has the same matrix rank as in the original hidden geometry (up to possibly transposing the matrix when a=b). No rank-one reference point or supplied frame is necessary. If a rank-one centered sample is available, the characteristic restriction q not dividing a or b is unnecessary for the tangent construction. For general samples, a conservative arithmetic-cost bound is

    O(|A| D^2) + poly(D,d,k,q) + q^k poly(d,q),

where the second term includes explicit construction of the linear symmetry algebra, its commutator and centroid, and basis operations. The polynomial exponents are absolute. Field operations and integer indices have their ordinary polynomial bit overhead. The explicit representation of the input and all evaluation-space dimensions are charged. The algorithm can retain a streaming evaluation basis, of size at most E, instead of the entire sample evaluation matrix.

In particular if q<=d^C for a fixed C, the q^k term is polynomial in D: k+1<=sqrt(d) and D>=(d/(k+1))^(k+1). Thus the full procedure is polynomial in the explicit sample length and D for that field-size regime (with a dependence on the fixed C). It is not polynomial in d,k alone when D is large. One may enumerate all factorizations d=ab and rank caps k satisfying a specified cost cap, with polynomial overhead in d. This is a guarded structured recognition theorem, not a promise that arbitrary dictionaries have useful geometry.

### 1.1 Degree-(k+1) finite-field ideal: an elementary proof

In ordinary matrix coordinates, consider a polynomial f of total degree <=t vanishing on every matrix of rank <=k. Restrict to matrices with at most k nonzero rows; all such matrices have rank <=k. Every individual variable degree is <q, so repeated univariate interpolation forces all coefficients supported on at most k distinct rows to vanish. The same argument with columns forces all remaining monomials to use at least t distinct columns. Therefore any surviving monomial has degree exactly t, using t distinct rows and columns once each.

Fix such row and column sets I,J of size t, and restrict to that square submatrix. The surviving terms are permutation monomials. Set two rows equal. The rank is <=k, and the resulting polynomial is zero on the whole grid. Individual degrees are <=2<q, so each coefficient vanishes. The permutation monomials become paired under swapping the two row assignments; their coefficients must sum to zero. Repeating for every row pair forces the alternating coefficients of a determinant. Consequently f is a linear combination of the t-by-t minors. Those minors are independent because different I,J have disjoint monomial supports. There are exactly B of them, and no nonzero polynomial of degree <=k vanishes on the rank ball.

Affine change of variables is an invertible transformation of the space of degree <=t polynomials. Thus the full vanishing space of the affine image is precisely the B-dimensional span of the translated minors. The rank test in the theorem proves equality of the sample ideal with this full space under the stated containment promise. This proof does not rely on Zariski density over a finite field or on a bound for arbitrary high-degree polynomials. Classical determinantal ideal theory gives a stronger algebraic statement; the finite-grid argument here establishes exactly the degree slice needed by the algorithm.

### 1.2 Thin random samples, including rank-exact samples

Assume q>=4t. Draw matrices as U V, where U in F^(a x k), V in F^(k x b) are independent uniform matrices, and apply the same unknown affine map. For any polynomial outside the full vanishing space, substitution of U V yields a nonzero formal polynomial of total degree <=2t. Nonzeroness follows from the preceding characterization and the fact every rank <=k matrix factors this way. Schwartz-Zippel bounds its vanishing probability by 2t/q<=1/2.

There are fewer than q^E nonzero cosets of polynomials modulo the full vanishing space. For L independent samples the failure probability is at most q^E 2^(-L), by a union bound. The same conclusion, with a weaker numerical constant, holds for samples uniformly distributed among rank-EXACTLY-k matrices. To see this, U and V are both full rank except with probability at most 2/(q(q-1)): a union bound over one-dimensional kernel or left-kernel lines bounds either failure by (q^k-1)/(q-1)*q^(-min(a,b)) <=1/(q(q-1)). Conditioned on both being full rank, every rank-k matrix has the same number |GL_k(F)| of factorizations. Thus the conditional vanishing probability is at most (2t/q)/(1-2/(q(q-1)))<=2/3 under the guard. The failure probability for L rank-exact samples is therefore <=q^E (2/3)^L. Taking

    L >= ceil((E log q + log(1/epsilon))/log(3/2))

suffices for success probability >=1-epsilon. Repeated samples are harmless: they are independent draws for this bound and can subsequently be deduplicated without changing the evaluation kernel. This is a probability theorem about the input ensemble, not randomized computation required by the recognition algorithm. The probability bound does not assert that the actual projected additive lists in the repository are determinantal samples.

## 2. Recover the exact affine center and descend to quadrics

Compute a basis f of the degree-at-most-t sample ideal. Write f_t for its top homogeneous part and f_(t-1) for the next part. Since f(v)=g(v-c) for a homogeneous polynomial g of degree t in the hidden linear coordinates, the top-degree relation is

    f_(t-1) = - D_c f_t,

where D_c is the formal directional derivative. These are explicit linear equations in the d coordinates of c. They are consistent under the promise. The common derivative kernel of all t-minors is zero: for a matrix entry ij choose a t-by-t minor containing it, and consider a monomial in the derivative on t-1 distinct complementary rows and columns. Its coefficient isolates the derivative direction's ij entry (with sign). This gives uniqueness of c in the ordinary coordinates and hence in the hidden coordinates. No factor-two loss from subtracting an arbitrary sample center is necessary.

Take the span Q of all (t-2)-fold partial derivatives of the top homogeneous forms, successively differentiating bases and retaining a basis after each step. In ordinary coordinates a derivative of a t-minor is either zero or a smaller minor up to sign. Every 2-minor is obtained by adjoining complementary rows and columns and differentiating a matching of distinct complementary entries. Linear coordinate transformations commute with forming the full span of directional derivatives. Hence Q is exactly the space of 2-by-2 minors in hidden linear coordinates, of dimension binom(a,2)binom(b,2). The argument does not divide by a high factorial. In characteristic odd, represent each quadratic by a symmetric matrix M with q_M(v)=v^T M v, dividing off-diagonal coefficients by two.

The common zero set of Q is precisely the cone of hidden rank <=1 matrices. This is an exact membership predicate computed from a basis of quadrics. It is not yet an algorithm for finding a nonzero member from a rank-k sample.

## 3. Recover the two matrix-action algebras

Let V=F^d. Solve the following linear equations for X in End(V):

    X^T M + M X in Q    for every M in Q,

where Q is now the space of symmetric matrices representing the quadrics. To implement them, take a basis of the annihilator of the quadratic coefficient space. Represent its elements as symmetric matrices H with entries equal to the dual monomial coefficients (including both off-diagonal entries). The constraints are Tr(H(X^T M+M X))=0, or equivalently the linear equations <M H,X>=0 after cancelling the common factor 2. The result is the Lie algebra g of linear infinitesimal preservers of Q. This approach is standard in polynomial equivalence algorithms; we need its exact form and a deterministic recovery procedure in the present finite-sample setting.

**Lemma 2.** In the original tensor coordinates V=F^a tensor F^b,

    g = { L tensor I_b + I_a tensor R : L in M_a(F), R in M_b(F) }.

Proof. The displayed operators preserve the quadratic ideal infinitesimally by the left and right rank-preserving actions. Conversely, at every nonzero rank-one vector u tensor v, X(u tensor v) must lie in its tangent space u tensor F^b + F^a tensor v. In particular X(e_i tensor f_j) is supported only on its row and column. Applying the tangent condition at e_i tensor (f_j+f_l) forces the off-diagonal left-action coefficients to be independent of the column. Applying it at (e_i+e_h) tensor f_j forces the right-action coefficients to be independent of the row. Subtract these two off-diagonal actions. The remaining operator is diagonal in the e_i tensor f_j basis, with coefficients c_ij. At (e_i+e_h) tensor (f_j+f_l), tangency gives c_ij+c_hl=c_il+c_hj. Hence c_ij=s_i+t_j, the sum of diagonal left and right actions. QED.

This proves dim g=a^2+b^2-1. Compute the commutator span g'=[g,g] from pairwise brackets of basis operators. Under q odd and q not dividing a,b,

    g' = (sl_a tensor I_b) direct-sum (I_a tensor sl_b).

Indeed [gl_n,gl_n]=sl_n; the two factors commute and their trace-zero scalar intersection is zero. Both simple ideals are absolutely simple under these conditions. One elementary proof extends scalars to an algebraic closure, takes a generic trace-zero diagonal matrix whose off-diagonal root eigenvalues are nonzero and distinct, and uses polynomials in its adjoint action to isolate a matrix unit in any nonzero ideal with an off-diagonal entry. A nonzero diagonal element in an ideal cannot be scalar because q does not divide n, and brackets with matrix units produce an off-diagonal entry. Matrix-unit brackets and characteristic not two then generate the whole sl_n. The root eigenvalue forms are distinct since their differences have sum zero, and the only vector proportional to the all-ones vector with sum zero is zero; q odd prevents equality of opposite roots. This also covers n=2.

Compute the centroid in the form of the commutant of all adjoint operators on g'. In a computed basis this is a linear system. The centroid is exactly F direct-sum F: the two cross maps vanish because each simple ideal is perfect and the other acts trivially on it, and on each absolutely simple adjoint representation the commutant is scalar (Schur's lemma after extending scalars, followed by descent). Choose a nonscalar element of the two-dimensional centroid. Its eigenspaces are the two simple ideals; the two distinct eigenvalues lie in F, so a scan of q field elements suffices. Adjoin the identity operator to each ideal. The recovered commuting associative action algebras are exactly

    A_left = M_a(F) tensor I_b,  A_right = I_a tensor M_b(F),

up to exchange. Their dimensions determine a,b up to order. All computations use explicit linear systems and bracket products. No randomized conjugacy search or finite-algebra isomorphism oracle is required. The restrictions on characteristic matter; the proof does not assume a split of simple ideals when those guards fail.

## 4. Find a rank-one direction in O(q^k) projective trials

Let x be any nonzero centered sample. In hidden coordinates its rank is some h<=k. Compute the action-orbit spans

    U=A_left x, V=A_right x.

If C and W are the column and row spaces of the hidden matrix, these are F^a tensor W and C tensor F^b. Their dimensions are a h,b h, and their intersection is C tensor W, dimension h^2. Thus h itself is recovered by ordinary linear algebra, without a tensor layout. In the h=1 case we already have a rank-one vector.

For h>=2 select any (h+1)-dimensional linear subspace J of U intersect V. In the hidden h-by-h matrix coordinates on C tensor W, projection onto one row maps J to an h-dimensional vector space. Its kernel contains a nonzero matrix, and that matrix has rank <h. Enumerate the projective directions in J, normalizing the first nonzero coefficient to one, and test rank decrease by the dimension of the left orbit. The guarantee finds a nonzero lower-rank point in at most (q^(h+1)-1)/(q-1) trials. Repeat. The sum of these bounds over h=2,...,k is O(q^k) with an absolute constant for q>=2, because rank strictly decreases at each step. Each trial is polynomial-dimensional linear algebra. This improves over direct enumeration of the entire h^2-dimensional corner, and avoids having to split an abstract matrix algebra into primitive idempotents. It does not remove the separate cost D of learning the vanishing ideal.

## 5. Recover the tensor layout using only quadrics and tangent spaces

Write B_Q(x,y)=q(x+y)-q(x)-q(y) for the vector of polarized quadratic forms. At a nonzero rank-one point x=u tensor v, the kernel of B_Q(x,.) is its tangent space T_x=u tensor F^b+F^a tensor v. All computations below are in the observed vector space.

Choose a complement W of the line F x in T_x. On this quotient the quadratic equations identify two ruling components, with coordinates (s,t) in F^(a-1) direct-sum F^(b-1), and the restrictions of the minors span all bilinear products s_i t_j. (Choose hidden bases with x the (1,1) unit; minors in rows 1,i and columns 1,j restrict to minus these products; other minors vanish.) The nonzero rank-one directions within the quotient have s=0 or t=0.

For a nonzero y in W, if all quadrics vanish we have such a direction. Otherwise both components s,t are nonzero. The kernel in W of the polar map B_Q(y,.) is a line, since the equation s tensor t'+s' tensor t=0 forces (s',t')=(alpha s,-alpha t). For a nonzero generator z of that line, the vectors y+lambda z become the two nonzero pure-component directions at the two opposite reciprocal values of lambda. Scan F to find them. This argument uses odd characteristic. With a rank-one y in W, the intersection T_x intersect T_y is the ruling through x,y. Pick a vector in W outside that ruling and repeat the pure-component split to find a point on the other ruling. We obtain the two ruling subspaces through x without enumerating their points or assuming they appear in the input sample.

Orient the two rulings as L=u tensor F^b and R=F^a tensor v (either orientation is fine in the square case). Extend x to bases a_j=u tensor v_j of L and b_i=u_i tensor v of R, taking a_1=b_1=x. For i,j>1, compute the other ruling at a_j and at b_i by the same algorithm. Their intersection is the line F(u_i tensor v_j). For a nonzero point z on this line, select the unique scalar lambda for which x+a_j+b_i+lambda z has rank one. The hidden two-by-two block is [[1,1],[1,c]], of rank one exactly at c=1, so the unique scalar exists. This gives the basis u_i tensor v_j in the observed coordinates, hence an explicit invertible map S preserving all matrix ranks. Verify the dimensions, uniqueness and invertibility at each step. This completes the reconstruction theorem.

The construction can be viewed as a tangent-space version of reconstructing the two Segre rulings. Its input is a basis of the quadratic ideal, not a complete rank-one dictionary with all pairwise addition queries. The additional sample, derivative, and Lie-algebra steps above remove the need for a supplied rank-one anchor.

## 6. What this changes for structured additive certificates

After reconstructing and verifying a coordinate system, let A_1,...,A_m be explicit lists with the same recovered center c, and compute actual rank bounds k_i=max rank(S^-1(x-c)) for x in A_i. If sum k_i<min(a,b), take Z of full rectangular rank and output m c + S vec(Z). This is outside the sumset by rank subadditivity. The local rank validations suffice on ANY input; no inferred geometric promise is used for correctness. If they fail, this branch honestly fails. Different centers are also allowed if supplied or separately learned and verified; replace m c by the sum of the centers.

For relevance to small fields, consider RECTANGULAR examples rather than requiring q to exceed a rank cap growing like sqrt(d) in a square shape. Let w grow and choose k>=1, q>=4(k+1) an odd prime polynomially bounded in w,d. Choose a approximately 4wk, b much larger than a, q dividing neither dimension, and d=ab. Assume a^2=o(d) and k log(d+k)=o(d log q /w). For example w,k can grow polylogarithmically in d, with a suitable sequence of rectangular dimensions and q. Then log_q D=O(k log(d+k)/log q)=o(d/w), and the rank-exact sampling threshold is q^o(d/w). In contrast the full rank ball has size q^Theta(d/w), since its rank-k population is given by the standard full-rank factor count. Thus a polynomially complete rank-exact sample can be an extremely small fraction of the ball and still permit recognition. For dictionaries of size q^Theta(d/w) the charged reconstruction cost is polynomial in that explicit size. Even allowing w arbitrary-coefficient low-rank atoms has rank cap wk<a. The usual sparse-span volume bound is meaningful: log_q |B_k|<=k(a+b), so w log_q(1+q|B_k|)<=d/4+o(d) at these parameters. This illustrates compatibility with the abstract high-rate, small-field regime. It is not a theorem that the actual digit-lifted dictionaries have such a rectangular layout. Nor does the square-rank asymptotic comparison silently require q>k when it does not hold.

## 7. A stronger affine rank-ball theorem: avoid a whole line of offsets

The exact affine center is useful beyond ordinary additive sums. In an arbitrary-coefficient w-sparse span of translated rank-at-most-k matrices, the offset coefficient can be ANY scalar. Prior rank-metric probe arguments (#28) supply several probes under stronger radius guards, or enumerate centers. A single affine offset permits a sharper elementary conclusion.

**Theorem 3 (full-rank pencil).** For any field F=F_q and any a-by-b matrix C, with min(a,b)>=2, a matrix Z can be constructed deterministically in poly(a,b,q) arithmetic work such that

    rank(Z-lambda C)=min(a,b) for every lambda in F.

Consequently, for every R<min(a,b), Z is outside F C + B_R. This does not require 2R<min(a,b), a rank-metric code, enumeration of all possible sparse representations, or even enumeration of the q centers in the construction. A certificate can additionally check the q ranks explicitly. The theorem applies in particular to a single affine center plus a sum or arbitrary-coefficient sparse span of low-rank deviations.

Proof. Transpose if necessary to assume a<=b. By invertible row and column operations put C into canonical form [diag(I_s,0) | 0] where s=rank C. Work in the first a columns and fill the extra columns with zeros. If s<a, choose Z_0 an a-cycle permutation matrix. Its inverse times diag(I_s,0) is nilpotent: the full cyclic basis shift is broken by at least one deleted link. Thus I-lambda Z_0^-1 diag(I_s,0) is invertible for all lambda, and so is Z_0-lambda diag(I_s,0).

If s=a, choose Z_0 with no eigenvalue in F. Every a>=2 is a sum of twos and threes. Use a block diagonal of companion matrices of irreducible quadratics and (when needed) one irreducible cubic. Over a finite field these exist; a deterministic scan of the q^2 or q^3 monic coefficient lists, rejecting polynomials with a root, is polynomial in q. Degrees two and three have no proper factor without a root, so the scan is exact. The companion matrices have no eigenvalues in F. Hence Z_0-lambda I is invertible for every lambda. Restore the original row and column coordinates. QED.

Notice that this theorem does not claim a rank-(R+1) ball centered at every point of an arbitrary higher-dimensional center subspace can always be avoided. The offset being one-dimensional is load-bearing. For s=1 and square size one the statement could fail, explaining the minimum-dimension guard.

### 7.1 Common-kernel/common-image center subspaces

**Theorem 4.** Let W be any explicitly specified linear subspace of a-by-b matrices, a<=b, and let

    U=sum_(C in W) im(C),   K=intersection_(C in W) ker(C),
    t=dim U.

If dim K>=t, a full-row-rank Z can be constructed by ordinary linear algebra with rank(Z-C)=a for every C in W. For a>b apply the condition and construction to transposes. The condition is checkable using a basis of W; neither its q^dim(W) elements nor any center combinations need be enumerated.

Proof. Choose a t-dimensional subspace K0 of K, and extend its basis to an a-dimensional domain subspace L and then a full domain basis. Choose an isomorphism Z|_L onto F^a that maps K0 onto U, and extend Z to the whole domain arbitrarily. If (Z-C)x=0 for x in L, then Zx=Cx in U, hence x in K0. But C kills K0, so Zx=0 and x=0. Thus Z-C restricted to L is injective, and has full row rank. QED.

This yields a certificate for some arbitrary-exception models: let verified low-rank deviations have total rank cap R<a, and let every remaining center contribution belong to W. Then S is contained in W+B_R and Z is certified absent. It can handle many distinct centers when their common geometry satisfies the dimension guard, even if q^dim(W) is too large to enumerate. The common-kernel guard is not implied by low dimension of W; in particular a line of full-rank centers fails it, and is instead handled by Theorem 3. These are separate sufficient tests, not a universal characterization of affine subspaces avoiding a determinantal variety.

## 8. Verification and limitations

From this directory run `python -B verify.py` and `python -B center_checks.py` (Python 3.10+, NumPy, assertions enabled). The scripts use exact modular integer arithmetic with small finite field primes and bounded matrix sizes, avoiding numerical rank. They do not implement a large-parameter bit-complexity solver. `verification.txt` records seven hidden layouts with 4,076 training vectors and 700 independent rank tests. Four fixtures have rank exactly two at every training vector, so the rank-one anchor is not supplied. The affine centers match the actual hidden offsets exactly. The rank-reduction search took 37 projective trials in those fixtures; this is a finite observation, not an asymptotic improvement over the proven worst-case bound. The sample sizes are twice the monomial count and their completeness is checked directly; some rank-one fixtures use primes below the sufficient random-sampling probability guard. The insufficient-sample control fails the ideal dimension test and the high-rank outlier control violates the validated local rank cap. `center_verification.txt` records 1,170 affine center matrices, including exhaustive 2-by-2 matrices over F_2,F_3,F_5, and 6,200 exact scalar shifts. It also records 48 center-space fixtures with 5,392 exact center combinations, plus one failed common-kernel guard. The line-pencil theorem works even in characteristic two, independent of the odd-characteristic reconstruction theorem. All finite outputs are illustrative abstract matrix controls.

Relevant prior research: determinantal ideals and finite-field evaluation codes (Beelen and Ghorpade et al., *Linear codes associated to determinantal varieties*, [DOI](https://doi.org/10.1016/j.disc.2015.03.009)); polynomial symmetries and Lie-algebra recovery (*Linear preservers of secant varieties and other varieties of tensors*, [DOI](https://doi.org/10.1016/j.jsc.2025.102449)); determinant equivalence and full matrix algebra isomorphism (Ivanyos et al., *Determinant Equivalence Test over Finite Fields and over Q*, [ICALP 2019](https://doi.org/10.4230/LIPIcs.ICALP.2019.62)); constructive finite-algebra structure (Rónyai, [DOI](https://doi.org/10.1016/S0747-7171(08)80017-X)). The present use of Lie-algebra linear equations, centroids, rank-metric dimension arguments, tangent rulings, and matrix pencils builds on classical mathematics. No claim is made that these components individually are novel. Independent priority review of the precise thin-sample reconstruction, complexity bounds and affine-center combinations is outstanding.

**Research verdict.** This replaces complete-ball and supplied-frame recognition by a different checkable sample condition, handles arbitrary affine centers exactly, and removes the factor-two radius restriction for one-dimensional center avoidance. The results still require a useful structured matrix layout and, for recognition, enough determinantal equations in the sample. There is no theorem that an arbitrary high-growth additive dictionary admits such a layout or other compact certificate. Therefore uniform additive avoidance and the project's main end-to-end exponent remain open. The actual elliptic relation predicate, contractor and coverage assumptions are not used by these abstract statements and have not been independently reverified in this note.
