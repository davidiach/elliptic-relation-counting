# Normalized chord residuals: exact box cancellation, dense odd support, and the first informative sketch

**Research checkpoint, 24 November 2026.** This is an original mathematical continuation of `average_case_detector_note.md`, Section 6, and the uploaded v3 preprint. It establishes exact statements about a particular modular residual representation and affine mask sampler. It does **not** improve the proposed general-prime ECDLP coefficient or the logarithmic saving power. It does not establish a lower bound on arbitrary implicit parity algorithms, nor a novelty or priority claim.

## 1. A normalized residual and complete-fiber cancellation

Let `E: y^2=f(x)=x^3+a x+c` be smooth over F_p, p>3. Set

    h(u,v)=(u+v)(uv+a)+2c,
    N(u,v,w)=(u-v)^2 w-h(u,v),
    g(u,v,w)=N(u,v,w)^2/(4 f(u) f(v))

when u!=v and f(u)f(v)!=0; set g=0 otherwise. On the nondegenerate domain the third summation polynomial satisfies

    (u-v)^2 S3(u,v,w)=N(u,v,w)^2-4 f(u)f(v).

Thus g=1 is the guarded chord equation. Given rational non-two-torsion target T, impose g(X_(j-1),x_j,X_j)=1 for j=2,...,m, with X_1=x_1, X_m=x(T). The internal abscissas range over the entire field, and the leaves over the full field-affine masks, including nonliftable abscissas. A g=1 chain has non-two-torsion intermediates (every one is checked as an input of the next edge). Backward rationality gives rational lifts of all leaves, as in the previous abscissa-chain note. For fixed signed target, a nondegenerate edge has exactly one ordered signed input pair with the two given distinct abscissas: the two possible output abscissas differ because f(u)f(v)!=0 and p!=2. The output sign determines the simultaneous input signs. Hence the exact g=1 count equals the guarded signed relation count, excluding equal-abscissa steps and two-torsion. We do not assert equality with the unguarded count.

For every function phi:F_p->F_2, even on degenerate rows,

    sum_{w in F_p} phi(g(u,v,w)) = phi(0)  (mod 2).             (1)

Indeed the nondegenerate map is a nonzero multiple of the square of an affine bijection. Every nonzero value has an even number of preimages and zero has exactly one. A degenerate row is constant zero on an odd number of inputs. This proves an exact cancellation identity in the actual elliptic arithmetic.

Let r=m-1. Encode each residual in b=ceil(log2 p) canonical binary bits (p odd). Let c(q) in F_2 be the parity of the number of leaf/internal-state assignments with residual vector q in F_p^r, extended by zero outside the valid binary encodings. Denote the genuine residual by tau=(1,...,1). A moment is the parity sum of a product of selected residual bits. Its block support is the subset of chord edges having at least one selected bit. If a selected edge is immediately followed by an unselected edge, the shared internal abscissa is free, and (1) kills the moment: the bit monomial in that edge vanishes at zero. Therefore any nonempty moment with nonzero value has suffix block support. This holds for arbitrary mask cardinalities. A degree-k sketch with k<r only sees a suffix of at most k edges, up to the known parity of earlier mask sizes. For fixed odd mask sizes and the actual affine seed law, it is a function of the suffix seeds and target alone; the prior unread-slot Q/R marginal estimate then applies. The stronger even-width results below do not require that marginal estimate.

## 2. Exact even-margin code and its box decomposition

Use deterministic even mask cardinalities w_i>=2, and enumerate all w_i abscissas, not only liftable ones. This is permissible for the normalized guarded chain just proved. If an entire residual block j is summed out, the corresponding leaf x_j appears in no remaining predicate (edge indices j=2,...,m). Summation over x_j supplies the factor w_j=0 in F_2. Hence c has zero sum on every coordinate line:

    sum_{q_j in F_p} c(q_2,...,q_m)=0

with every other coordinate fixed. Equivalently c is in the tensor product of r binary even-weight codes on the p-element alphabet. This proof of the margins needs only the chain locality and the even mask sizes; (1) supplied the stronger odd-size suffix claim.

For q having all coordinates nonzero, let V_q be the r-dimensional binary subspace with independently chosen q_j or 0 in each disjoint bit block. Then there is the exact, unique anchored-box decomposition

    c = XOR_{q in (F_p^*)^r : c(q)=1} 1_{V_q}.              (2)

Both sides have zero line sums, and agree on all full-support q. Their difference vanishes by induction on the number of zero coordinates: the line-sum identity writes any entry with one chosen zero coordinate as the sum of entries with that coordinate nonzero. Equivalently {e_a+e_0 : a!=0} is a basis of the even-weight code. In particular every nonzero c has at least 2^r nonzero entries, by induction: a nonzero slice along one axis is accompanied by another nonzero slice, and each satisfies the remaining even-margin conditions.

The true residual tau, when its count is odd, contributes the anchored box V_tau. It has 2^r entries, including tau itself and patterns with one or more zero sentinels. Thus the desired certificate participates in the cancellation structure, rather than sitting outside it as an unpaired singleton.

For any binary k-by-rb matrix A, consider the complete parity table

    Z_A(y)=sum_{q: A q=y} c(q),    y in F_2^k.

The contribution of each box V_q is zero identically if rank(A|V_q)<r; otherwise it is the indicator of the r-dimensional image A V_q. Consequently:

* k<r: Z_A is identically zero for every A, even on positive inputs.
* k=r: Z_A is constant over all y, equal to the XOR over active boxes of [rank(A|V_q)=r]. In characteristic two this is a determinant of the columns A_j q_j. The genuine syndrome is not distinguished from any other syndrome.
* k>=r: Z_A is a multilinear Boolean polynomial of degree <=k-r as a function of y. This follows because the indicator of a codimension k-r linear subspace is a product of k-r affine-linear equations. If Z_A is not zero, its support has at least 2^r entries, by the standard elementary degree-weight inequality for Boolean polynomials (split on the last variable and induct using the derivative of degree one less).

These are exact assertions for every seed and arbitrary A, not just independent random A. They rule out preserving the true odd count while cancelling all false parity by fewer than r checks. They also rule out a sketch of any width whose projected parity table consists of just the one genuine syndrome. They do NOT exclude prediction using one distinguished syndrome, when k>=r, based on a statistical difference between seed laws. Nor do they exclude a nonlinear transformation of the residual coordinates before the projection.

## 3. A dense-support theorem for the actual normalized elliptic residual map

We now prove that the canonical box representation is almost as large as the entire candidate assignment domain at suitably thinned widths. This is uniform at every fixed non-two-torsion target, including the polynomial-cofactor class when the target is eligible.

Use independent field-affine masks a_i+t_i S_i, with a_i uniform in F_p, t_i uniform in F_p^*, and each support having the same even size w. Their exact one-point retention is rho=w/p and distinct-pair retention is w(w-1)/(p(p-1))<=rho^2. Set

    lambda_0=w^m/p,
    B_0=w^m p^(m-2)=lambda_0 p^r.

B_0 is the number of all candidate paths at each seed, before any residual test. Let I count paths having a zero residual in some block. Let C count ordered pairs of DISTINCT paths with identical full-support residual q. Let U count full-support residuals with exactly one path, and B the number of active full-support boxes in (2). Deterministically

    B_0-I-C <= U <= B <= B_0.                              (3)

Indeed a fiber of multiplicity t>=2 accounts for t paths and t(t-1)>=t ordered colliding pairs.

**Theorem.** Uniformly in the fixed eligible target,

    E I / B_0 <= 9m/p,
    E C / B_0 <= 4^r lambda_0 (1+1/w)^m + 4m 2^(m-2)/p.    (4)

Thus if epsilon is the sum of the displayed upper bounds, Markov's inequality gives B>= (1-sqrt(epsilon)) B_0 with probability at least 1-sqrt(epsilon) (for epsilon<=1). For growing m<=sqrt(b) and w=2 floor(p^(1/m)/32), lambda_0=16^(-m)(1+o(1)), m/w=o(1), and epsilon=o(1). Therefore

    B=(1-o(1)) lambda_0 p^(m-1)                             (5)

with probability 1-o(1), and log2 B=(m-1)b-O(m) up to the replacement of b by log2 p and the small relative error. The thinning changes the inverse nominal relation density by only 2^O(m), acceptable in the earlier square-root recovery accounting; this does not assert a detector has sufficient advantage at those widths.

### 3.1. Invalid residuals and bounded fibers

Count first paths on the full field grid of p^(2m-2) leaf/internal choices. The exceptional cases u=v or f(u)f(v)=0 cost at most 7/p per edge. At an interior edge N=0 costs at most 1/p after those cases, since w has nonzero linear coefficient. At the final edge, w=x(T) and f(w)!=0, the polynomial N(u,v,w) in v is nonzero of degree at most two. This also follows from the invertible matrix in 3.3 below. Thus each edge is bad on at most 9/p of full-grid choices. Every first path's mask-inclusion probability is rho^m, giving the first inequality.

For any prescribed nonzero q and two of u,v,w fixed, there are at most four possibilities for the third. In w there are at most two, since the nondegenerate equation is quadratic. In u or v the equation is of the form N^2=4q f(u)f(v), a polynomial equation of degree at most four. It is nonzero: for example, with v fixed and f(v)!=0 the right side has degree three in u, whereas a nonzero polynomial square has even degree. All degenerate solutions are excluded by q!=0.

### 3.2. Pairs with at least one different leaf

Fix a first full-support path and its residual vector. Let J be the exact set of coinciding leaf abscissas of a second path, with s=|J|<m. Choose an unmatched leaf j0. Choose the other unmatched leaves freely, at most p^(m-s-1) possibilities, and fix the matching leaves to the first path's values. If j0>=2, propagate the second path forward from its x1 until edge j0 and backward from the fixed endpoint through edge j0+1. The bounded fiber result determines at most 4^r possible paths, solving the unmatched leaf at the bridge with its two endpoints now specified. If j0=1, simply propagate backward to determine x1. This counts every second path. Joint mask inclusion probability is at most rho^(2m-s). Summing over first paths and proper J gives

    E C_proper <= p^(2m-2) 4^r sum_{s=0}^{m-1} binom(m,s) p^(m-s-1) rho^(2m-s)
               <= B_0 4^r lambda_0 (1+1/w)^m.

### 3.3. Same-leaf pairs: fixed-target reconvergence

A second path with the identical leaf tuple and identical full-support residual vector is specified by at most two choices at each intermediate state, once the preceding state is fixed. Since both paths have the same starting and ending abscissa, every distinct second path has an edge where the preceding states u,u' differ but the outgoing state w agrees. At this edge the common leaf v must satisfy

    f(u') N(u,v,w)^2 - f(u) N(u',v,w)^2=0.                (6)

This is a nonzero polynomial in v of degree <=4, even with w fixed. To prove it, the three coefficients of N as a polynomial in v are the matrix product

    [ 0   -1    w ] [u^2]
    [-1  -2w   -a ] [ u ]
    [ w   -a   -2c] [ 1 ].

The determinant of the matrix is 2f(w). On a full-support chain, all shared internal w have f(w)!=0, since they enter the following nonzero residual; the prescribed terminal w also has f(w)!=0. If (6) were identically zero, the two N polynomials would be proportional over the algebraic closure, by taking square roots of the nonzero scalar ratio. Invertibility would make (u^2,u,1) and (u'^2,u',1) proportional, hence u=u', contradiction.

For a candidate reconvergence edge j, fix all first-path variables other than its leaf v. Earlier second-path branches depend only on earlier first-path variables, so can be fixed independently of v (use a fixed order on the at-most-two choices). There are at most four possible v by (6). The number of branch patterns before and after the reconvergence is at most 2^(m-2) altogether; imposing reconvergence fixes that outgoing state if it is internal. Sum over at most m edges to bound the number of ordered same-leaf collision pairs by 4m 2^(m-2) p^(2m-3). Their inclusion probability is rho^m, proving the remaining term of (4). No averaging over targets, genericity assumption, or independence of the residual values is used.

## 4. Even after random projection, there are too many active boxes to enumerate

Fix a seed whose active set Bset has B~lambda_0 p^r as in (5). For an independent uniformly random k-by-rb binary matrix A, a false active box q!=tau survives at the genuine syndrome A tau exactly if A|V_q has rank r and A tau belongs to A V_q. Since tau is not in V_q, the restriction and A tau are independent. Conditional on A tau=y!=0, the probability is

    pi=prod_{i=0}^{r-1}(1-2^(i-k)) * (2^r-1)/(2^k-1),  k>=r.  (7)

The first factor is bounded below by an absolute positive constant, so pi is Theta(2^(r-k)). Write Y for the number of surviving false active boxes (ordinary integer number, not parity). Its conditional mean is (B-1_{tau active}) pi.

There is also concentration when r<=k=o(b), m<=sqrt(b) grows and lambda_0=16^(-m)(1+o(1)). For two distinct q,q' not sharing any coordinate, V_q intersects V_q' trivially. If additionally tau is not in V_q+V_q', the restrictions to these two subspaces are independent conditional on A tau=y. To bound exceptional ordered pairs, count at most B r p^(r-1) for a shared coordinate and at most B r p^(r-1) for q having some coordinate equal to one. For any remaining q, tau in V_q+V_q' forces every coordinate q'_j to be either 1 or (the binary vector) q_j XOR 1, hence at most 2^r choices. The total exceptional pairs is at most B(2r p^(r-1)+2^r). Bounding each covariance by one and the diagonal variances by pi gives

    Var(Y | A tau=y) <= B pi + B(2r p^(r-1)+2^r).

Relative to mean squared this is o(1): pi^-2=2^o(b), lambda_0^-1=2^O(m), and r/(lambda_0 p)=2^(-b+o(b)). The probability A tau=0 is 2^-k=o(1). Therefore with probability 1-o(1), over both seeds and A,

    Y=Theta(lambda_0 p^r 2^(r-k)),
    log2 Y=(m-1)log2 p-k-O(m).

This excludes enumeration of surviving boxes at subexponential widths. It is **not a parity-bias theorem**: an enormous number of surviving boxes can still have even parity, and pairwise concentration is insufficient to bound the parity of Y. In particular, it does not rule out a genuinely implicit evaluation of their XOR.

## 5. Exact bias of two random parity checks

This separate identity applies to any binary residual histogram c(v) in F_2^n, after translating the desired residual to zero if necessary. Put

    t0=sum_v c(v),  t1=sum_v c(v)v,  M2=sum_v c(v) vv^T.

All sums are in F_2. For two independent uniform binary rows a,b, the parity of the residuals satisfying both checks is

    Z=t0+a.t1+b.t1+a^T M2 b.

**Exact identity:** writing R=rank_F2(M2),

    E_{a,b} (-1)^Z = (-1)^(t0+R) 2^(-R).                 (8)

Proof: average over b. The conditional bias is zero unless M2 a=t1. Since t1=diag(M2), it is orthogonal to the kernel (x^T M2 x=x.t1=0 for every kernel vector). Thus the equation is solvable and holds for 2^(-R) of a. On these solutions a.t1 is constant and equals R modulo two. For the last assertion, pass to the quotient by the radical of the symmetric bilinear form. Its characteristic vector a represents x->B(x,x). Split off any norm-one vector orthogonally and induct. A remaining nondegenerate alternating form has even dimension and characteristic vector zero. Each split one-dimensional component contributes one both to rank and characteristic norm. This proves (8).

With more than two independent checks, condition on the first k-2 rows, filter the residual histogram, and apply (8) to the remaining two. This gives a rank-based exact average over the final two random rows rather than a Monte Carlo estimate over them; computing the filtered moments still has a cost.

For the three-slot even-width normalized chain, all moments missing either residual block vanish, so t0=t1=0 and M2=[0 C; C^T 0], with C the b-by-b cross-moment matrix of the two blocks. Its rank is 2 rank(C). Thus

    P_A(Z_A(A tau)=1 | seed) = (1-4^(-rank(C)))/2.          (9)

The syndrome translation is irrelevant at k=r=2, as already proved. This permits exact finite evaluation of the mean sketch response from a b-by-b rank, without sampling A. It is not a subexponential algorithm for computing C at growing bitlength: the implementation forms it by summing over all intermediate abscissas.

A second tested score uses the complement of the zero-sentinel test, accepting all nondegenerate genuine residuals. With even mask sizes, every proper-subset product of local checks has even sum; hence the parity of the product of the nonzero tests equals that of the product of the zero tests. Away from degeneracies g=0 is the cubic vertex equation (u-v)^2 w=h(u,v). This provides a cheaper-looking alternative local equation, but degeneracies, encoding multiplicities and the detector advantage must still be charged. No new numerical coefficient follows from this observation.

## 6. Finite checks and statistical experiments

Run the two Python programs with NumPy. `modular_cancellation_audit.py` tests 180 arbitrary binary histograms and all 163,800 ordered random-check pairs in their dimensions; 12,028 complete field rows of (1); 432 guarded-chain counts against direct group addition; and 2,497 parity-histogram, syndrome, and moment identities. `box_collision_audit.py` tests 96 exact box decompositions and 1,872 random projections, 1,827 bounded bridge fibers and 27,968 fixed-output reconvergence comparisons. It enumerates complete candidate grids at four small configurations, including 6,290,016 paths over F_65521 (12 seeds, arity three, width two). The observed active-box fraction there is 0.99978029, consistent with (4). These are not executions of a subexponential relation evaluator.

The finite statistical probe uses Q and R exactly as in the affine-seed and averaged-detector notes, on prime-order curves at fields 101,257,509,1021. Each planted tuple is uniform conditioned on sum to a fixed nonidentity target, and its offsets and slopes are independently sampled; R uses independent uniform affine point marks. Both have full abscissa masks of fixed even width. There are 60,000 independent draws per law in total (16,000 per law on each of the first three fields, 12,000 on the fourth), with three slots. The two-check rate is averaged exactly over the random check rows using (9), and the vertex score is evaluated by a separate sparse-per-row transition. These are fixed, untrained scores; the raw averaged gap does not rule out using A as part of a more elaborate decoder.

| p | width | Q-R two-check gap | Q-R vertex gap |
|---:|---:|---:|---:|
| 101 | 4 | -0.0000296860 | +0.000437 |
| 257 | 4 | +0.0000311787 | -0.002125 |
| 509 | 6 | -0.0000001913 | +0.003875 |
| 1021 | 8 | -0.0000000421 | -0.000083 |

All pointwise 95% Hoeffding intervals contain zero. For n draws per law the two-check score is in [0,1/2], so its difference half-width is sqrt(log(40)/n)/2; the vertex bit has half-width sqrt(log(40)/n). These conservative intervals do not establish vanishing asymptotic advantage. Sample rank histograms and sample standard errors are in the transcript. No evidence of sufficient advantage is supplied.

## 7. A separate gate-level degree limit

For a Boolean encoding with N free candidate and gate bits z, suppose the base weight B has degree d0 and every residual bit has degree <=d in z. For any fixed seed and any matrix A chosen for that seed, the relaxed parity sum of B times k affine linear residual checks vanishes identically when d0+kd<N. The product has Boolean degree at most d0+kd, and summing every monomial over the full cube annihilates it unless it contains all N variables. If A is fixed independently of seeds and the corresponding *joint* degrees in (z,seed bits) obey the same bounds, the output polynomial has seed-degree at most d0+kd-N, by the same full-cube coefficient argument. In particular, a quadratic gate-consistency encoding on a full cube with constant base weight cannot keep nonzero parity with fewer than ceil(N/2) linear checks; at exact equality its output is seed-independent under the joint-degree hypothesis. This does not apply to arbitrary nonlinear constraints, unrestricted base weights or non-cubical candidate domains.

## 8. Status and next gate

The normalized residual gives genuine structured cancellation in the elliptic chain. But the true relation is part of an even-margin tensor code, and at thinned widths the anchored-box expansion is essentially maximally dense. Random low-width sketches cannot isolate a sparse list of boxes, and the threshold two-check and cubic vertex scores did not produce a statistically established planted advantage. A new approach must evaluate the XOR of exponentially many boxes **implicitly**, or use a substantially different, nonlocal residual representation. One precise remaining target is a paid algorithm for a higher-order full-support moment/filtered rank statistic with a provable planted-versus-matched advantage. The present result does not meet that target. The earlier candidate coefficient 12.39481703 remains conditional and unchanged; the logarithmic saving power 1/3 and the square-root objective remain open.

Sources within the project: uploaded `paper (1).pdf` v3, Sections 2--5 and 8; `abscissa_chain_note.md`; `affine_digit_sampler_note.md`; `average_case_detector_note.md`, Sections 2 and 6. The proofs in this checkpoint are mathematical derivations, not claims taken from those sources. The finite verifier's prime-order curve and affine-mark helpers follow the preceding checkpoint, while the residual, box, collision and rank checks are new here.
