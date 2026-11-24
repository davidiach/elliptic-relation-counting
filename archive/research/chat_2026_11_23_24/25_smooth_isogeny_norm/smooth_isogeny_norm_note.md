# Smooth auxiliary isogeny towers: exact scalar norms and rational traces

Research continuation, 24 November 2026. Previous checkpoint: `auxiliary_torsion_note.md`. **No improved general-prime ECDLP exponent is established.** The new results replace the auxiliary division-polynomial membership oracle with a succinct circuit for the exact rational kernel polynomial, prove a width-independent scalar norm and rational trace elimination, and identify a geometric obstruction to keeping bounded-degree polynomials or rational functions during simultaneous isogeny descent. This checkpoint does not re-audit the earlier conditional numerical coefficient 11.31486814 or the complete relation-reuse reduction. No priority claim is made.

## 1. Smooth-order auxiliary support for all sufficiently large p

Let the target curve E/F_p have order hr, r odd prime and h=poly(b), b=ceil(log_2 p). Take m tending to infinity, m<=sqrt(b). Set L=7^k to be the largest power of seven at most r^(1/m)/8, and w=(L-1)/2. All asymptotic statements exclude bounded exceptional inputs (in particular p<=7). Then log_2 w=b/m+O(log(b)/m+1), and w^m/r=2^{-Theta(m)}. Hasse gives r=p/poly(b) up to constant factors.

Sample auxiliary nonsingular short Weierstrass curves E' over F_p with uniformly random coefficients. This induces Howe's weighted isomorphism-class law, because each class has (p-1)/|Aut(E')| models. Instead of requiring that the rational 7-primary subgroup be cyclic, accept a curve with L^2 dividing #E'(F_p). Howe, *On the group orders of elliptic curves over finite fields*, Theorem 1.1, gives probability at least L^{-2}-O(L^2/sqrt(p)): the main term for any prime power N is >=1/N and the multiplicative error factors are bounded for the fixed prime 7. Since L^4=o(sqrt(p)), this is Omega(L^{-2}). Schoof's point-counting algorithm is polynomial in b, so the expected setup cost is 2^{O(b/m)} poly(b). This bound does not rely on useful factors of p-1 or p+1. Source: https://arxiv.org/abs/math/0110262 (Howe's original Compositio Math. paper, 1993); Theorem 1.1 also appears in https://www.researchgate.net/publication/242380492_On_the_group_orders_of_elliptic_curves_over_finite_fields_Compositio_Mathematica_85_1993_229-247.

The 7-primary subgroup of E'(F_p) is a product of at most two cyclic 7-groups. If its size is at least 7^{2k}, its exponent is at least 7^k. Project a uniform point onto the 7-primary subgroup by multiplying by the prime-to-seven part of the order. In a largest cyclic factor, at least 6/7 of the elements have maximal order. Repeated multiplication by 7 determines the order, and multiplication by the appropriate power of 7 yields a rational point H of exact order L. This works even when p=1 mod 7.

Let S={x([j]H):1<=j<=(L-1)/2}. There are exactly w distinct abscissas and no points of order two. One may enumerate S and choose disjoint affine offsets within the preceding setup bound. The shared-translation average coverage proof of Section 2 of `auxiliary_torsion_note.md` applies to *any* public support of cardinality w: greedily choose disjoint translates, avoid O((mw)^2) bad common multipliers for branch-locus intersections, and take a uniform common translation A. The first moment is (1+o(1))w^m/(r-1) by Hasse–Weil on the multiquadratic curve; the factorial second moment uses the degree h^2 2^(d-1) sum c_j^2 of the projected group-sum difference for d distinct offsets. This gives average inverse success at most 2^{O(m)} poly(b). Neither a new targetwise Poisson statement nor a full relation-reuse proof is claimed. Enumerating the bases, avoiding bad multipliers and selecting the auxiliary curve are all charged.

## 2. A division-free circuit for the exact kernel polynomial

Build a tower of normalized degree-seven isogenies

    E'=E_0 -> E_1 -> ... -> E_k=E'/<H>,

where the kernel at step j is the order-seven subgroup of the image of H in E_j. Only three nonidentity abscissas up to sign are needed at each step. The image of H loses a factor of seven in order. Normalized Vélu formulas (or elementary pairwise group addition) give the Kummer x-map

    X -> a_j(X)/b_j(X)^2,

with monic a_j,b_j of degrees 7,3. In detail, for the three abscissas u of the local kernel, with f_j(X)=X^3+a X+c the domain cubic,

    b_j(X)=prod_u(X-u),
    a_j(X)=b_j(X)^2 [ X+sum_u(2 f_j'(u)/(X-u)+4 f_j(u)/(X-u)^2) ].

The displayed expression is a polynomial, with canceled denominators evaluated as polynomial products; no divisions are used inside the circuit. This follows by adding the two x-coordinates of P±Q and subtracting 2x(Q) in Vélu's x-map. Codomain coefficients use the usual Vélu formulas. The constant-degree stage setup and evaluation cost is polynomial in b and k. For a reference on normalized isogenies and the squared odd-kernel denominator, see https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html and https://doi.org/10.1007/s40993-026-00754-x (standard Kummer denominator discussion; the composition below is elementary).

Denote by a_j^h,b_j^h the separate homogenizations to their declared degrees. For arbitrary ring elements U,V, initialize X_0=U,Z_0=V,K_0=1 and compute

    B_j=b_j^h(X_j,Z_j),
    X_(j+1)=a_j^h(X_j,Z_j),
    Z_(j+1)=Z_j B_j^2,
    K_(j+1)=K_j B_j.

Then Z_j=V K_j^2 identically. X_j,Z_j are homogeneous degree 7^j, K_j degree (7^j-1)/2, and the pair is the homogeneous Kummer map of the composite isogeny. Since compositions of finite morphisms of the projective line have no base points, the pullback of infinity has exactly the multiplicity-two finite preimages at the nonidentity odd-kernel abscissas, and a simple pole at infinity. The kernel of the composite is exactly <H>. Monicity follows inductively by specializing V=0: X_j(U,0)=U^(7^j) and K_j(U,0)=U^((7^j-1)/2). Therefore

    K_k(U,V) = prod_{u in S}(U-uV).                 (1)

There are O(k) arithmetic gates (degree seven is fixed), rather than O(w) coefficients. Equation (1) is an exact rational-support polynomial, not the full division polynomial. The original prime-order auxiliary support required retaining U^p-U to exclude nonrational roots in a zero test. This tower has no such extraneous roots and computes scalar values, not just a vanishing predicate. It can be evaluated inside any commutative coefficient algebra, including an algebra with zero divisors.

### Prime-power lifts

The construction works modulo p^P, with polynomial overhead in P and b. Lift the chosen nonsingular coefficients to Z/p^P. Since L is prime to p, the L-torsion is finite etale; H has a unique torsion lift to every such precision. Concretely x(H) is a simple root of the odd division polynomial psi_L; use the logarithmic-index recurrence from the previous checkpoint and Hensel lifting (first derivatives by automatic differentiation). The ordinate lifts by 2y(H) a unit. This constructs the lift in poly(P,b) time, not by enumerating L torsion points. All local-kernel points and their images in the tower then lift by prime-to-p group operations. Apply the same polynomial Vélu formulas over the ring. Distinct kernel abscissas have differences that are units because their reductions are distinct. Equation (1) remains true over the lifted ring: prove it over the unramified p-adic integer ring by the characteristic-zero morphism argument and reduce, or by uniqueness of the lifted finite etale kernel. The supports at higher precision are the Hensel-lifted kernel abscissas, not an arbitrarily selected set of integer lifts. Previous mixed-characteristic readouts only require the reductions, so this choice is allowed. An exact product for arbitrary chosen integer lifts would need a different circuit.

## 3. Width-independent exact norm elimination

Let B denote a coefficient algebra, e.g. a product of field lanes on the remaining grid or a product of unramified p-power rings. For F(U) monic of degree d, calculate K_k(U,1) in B[U]/(F), using (1)'s circuit. Each gate costs poly(d) underlying B-gates, so the full reduction costs poly(d,b) B-gates. Resultant reciprocity gives

    prod_{s in S} F(s) = (-1)^(wd) det(m_{K_k mod F}),              (2)

where m_q is multiplication by q in B[U]/(F). Compute its d-by-d determinant by a division-free polynomial-size algorithm, e.g. Berkowitz. The finite verifier uses a small determinant subset DP instead; the polynomial-time claim rests on the standard division-free determinant algorithm, not the finite implementation.

For a declared degree d but possibly dropping leading coefficient, use the projective chart (U,V)=(alpha X+Y,X), with determinant -1. On a lane with lambda=F^h(alpha,1) a unit, put g(X)=F^h(alpha X+1,X)/lambda, which is monic of degree d even if the original affine degree drops. Evaluate K_k(alpha X+1,X) modulo g using the homogeneous circuit. Homogeneous resultant covariance and swapping the two forms both introduce (-1)^(wd); they cancel. The result is the exact identity

    prod_{s in S} F(s) = lambda^w det(m_{K_k(alpha X+1,X) mod g}).   (3)

The declared degree of the transformed kernel is w even if its leading coefficient vanishes. The exponentiation lambda^w takes O(log w) ring gates. Unit inversion is performed by exponentiation over a finite field or Newton lifting in the p-power ring. Over N product lanes that reduce to nonzero F polynomials modulo p, choosing alpha in an extension F_{p^e} costs failure probability at most dN/p^e by the root bound; pick e=2 or larger as necessary. The needed extension and inverses must be charged. Zero-polynomial lanes in characteristic p have norm zero for w>0; in precision P, if all coefficients are divisible by p, the norm is zero for w>=P. They can be masked in a product ring using the reduction-to-zero idempotents (or just excluded under the earlier nonzero-lane hypotheses). If w<P, direct multiplication of the w factors is polynomial in P. The norm claim is for a coefficient algebra operation count; it is *not* a statement that a B operation costs one base-field operation.

This improves the previous interval norm cost ~O(d+sqrt(d w)) B-operations and the previous prime-torsion zero detector. It yields exact scalar norms modulo p^P, compatible with the mixed-characteristic coordinate recovery from earlier checkpoints. For d=2^{O(m)} poly(b,h), the symbolic cost is 2^{O(m)} poly(b,h,P) with no exponential width term. It is meaningful for m much smaller than sqrt(b), where the earlier sqrt(w) factor dominated. It does not itself solve the trace over B.

## 4. Direct rational traces, including repeated poles

Let Q=K_k(U,1), with roots S in the coefficient ring and pairwise unit differences. Let F monic of degree d, P(U) of degree <d, and assume F(s) is a unit for all s in S. In B[U]/(F), Q is a unit: the resultant is the product of these units up to sign. Compute q=Q mod F and q'=Q' mod F by evaluating the short homogeneous circuit and propagating its first derivative via dual numbers. Invert q in the d-dimensional algebra by the adjugate/determinant or other polynomial-size ring arithmetic; the determinant is a unit. Set

    g = P q' q^{-1} mod F.

Then

    sum_{s in S} P(s)/F(s) = - [U^(d-1)] g.                      (4)

Proof: because P Q'-g Q is divisible by F, write it as F B_0. We have deg B_0<=w-1. At each root s, B_0(s)=P(s)Q'(s)/F(s). Interpolation at the split simple roots gives [U^(w-1)] B_0 = sum_s B_0(s)/Q'(s). Compare the coefficient of U^(d+w-1) in P Q'-g Q=F B_0. The first term has degree at most d+w-2, so [U^(w-1)] B_0 = -[U^(d-1)]g. This argument needs no assumption that F itself is squarefree, and works over commutative rings with unit inverses as stated. The cost is poly(d,b) B-gates, again not a paid scalar runtime over the full grid.

A polynomial part of a rational function can also be summed without expanding Q. Evaluate its reverse polynomial Q*(t)=K_k(1,t)=prod_s(1-st) in a truncated ring B[t]/(t^(jmax+1)). The formal identity -t(Q*)'/Q* = sum_{j>=1}(sum_s s^j)t^j yields all desired low-order power sums. The constant term Q*(0)=1 is a unit. Both procedures work at prime-power precision, so (4) applies to the valuation-aware rational filters with unit denominator or in a fraction ring after a suitable parameter shift. If F has zeros on S then (4) cannot invert Q; one can use a formal parameter t with F+t or the separate norm/valuation construction, not apply a nonexistent inverse.

## 5. Why simultaneous descent does not preserve small degree

Here is an exact geometric obstruction for the generic isogeny descent, as distinct from a runtime lower bound for a circuit.

Fix an algebraic closure k of F_p. On the auxiliary product E_0^m with independent variable points R_i, set x_i=A_i+t_i x(R_i) (t_i nonzero) and y_i^2=f(x_i) for the target cubic f. Choose affine constants generically so that the three roots of f, transformed by the affine map, avoid the three finite branch abscissas of E_0 -> P^1_x. Then each f(x_i) has six distinct simple zeros on the auxiliary curve and an even-order pole at O. Denote its branch set by B_i (six points). The constants can be chosen over the algebraic closure, or for sufficiently large p avoiding finitely many forbidden choices. For the shared-translation scheme this is a generic-parameter statement; it does not require every actual seed to meet these conditions.

Let K=k(E_0^m). Adjoining all y_i gives a multiquadratic extension with group (Z/2)^m: the odd valuations at the disjoint vertical branch divisors B_i prove independence of the square classes. Put X=x([h] sum_i (x_i,y_i)). For m>=2, this function has stabilizer exactly the simultaneous sign reversal (including the identity). Indeed equality of abscissas of two signed sums implies equality or opposite equality of the group points on the integral product cover. Any nontrivial such identity would impose a nonzero multiple (with coefficient +/-2h) of a nonconstant coordinate map equal to a function of the other independent variables. This is impossible because multiplication by a nonzero integer on the elliptic curve is nonconstant. Thus the monic minimal polynomial F_0(X) over K is irreducible of degree d_X=2^(m-1), with splitting field the quotient of the multiquadratic extension by global sign reversal. Every individual sign flip acts nontrivially in this quotient. The branch locus of the splitting field over K is exactly the union of the vertical divisors coming from the six points B_i in each coordinate.

Let J be a subset of the coordinates, and translate each R_i for i in J by the rational subgroup H of order L=7^k. Write G=H^J. The field extension K/K^G is Galois, realized by the product isogeny quotient in the coordinates of J. If a translation g in G stabilizes the monic polynomial F_0(X), it stabilizes the splitting-field branch locus. Its ith component must preserve the six-element set B_i. Translation by a nonzero element of the 7-group has all orbits of cardinality a positive power of seven, which cannot partition a nonempty set of six points. Hence the stabilizer is trivial.

It follows that all |G|=L^|J| translates of F_0 are distinct monic irreducibles in K[X]. Their product

    N_J(X)=prod_{g in G} g(F_0(X))

is irreducible over K^G: any monic factor over the fixed field, viewed over K, would select a G-invariant subset of the irreducible factors; transitivity leaves only the empty or full subset. It has exact degree

    deg_X N_J=2^(m-1) L^|J|.                                    (5)

In particular it is squarefree over K (the original minimal polynomial is separable). The corresponding trace of a logarithmic derivative is N_J'/N_J, whose reduced denominator has exactly this degree, since all roots are simple with residue one. Thus even if a univariate norm along a curve isogeny preserves divisor degree in its *own* quotient coordinate, the norm of the multi-coordinate elliptic relation has exponentially growing target degree; it cannot be kept as a bounded-degree polynomial or a low-degree rational log derivative in X throughout generic simultaneous descent by removing repetitions. The use of the fixed prime seven makes the stabilizer argument uniform. In a target curve isogenous to the auxiliary curve with a compatible addition map, the divisor can behave differently; the six-point branch-locus hypothesis is important.

**Scope:** (5) is for generic auxiliary-point parameters, not an assertion of irreducibility or maximal degree after the specialization to the kernel fiber over O. That specialization has x(P)=x(-P) duplicates, possible infinity factors and further accidental collisions. It also does not establish a lower bound against computing one value by a succinct circuit. A polynomial of huge degree may have a short circuit. It specifically rejects the degree-preserving dense-polynomial or reduced-rational-function descent as a generic invariant, not all implicit algorithms.

## 6. Charged end-to-end status

The exact support/norm circuit removes the slot-width dependence in one-coordinate elimination. If the other m-1 slots are materialized, the evaluation costs about w^(m-1) times 2^{O(m)} poly(b,h,P) scalar work, with leading log-time b-b/m+O(m+log b). This is below the previous explicit one-coordinate interval-block estimate b-b/(2m)+O(m+log b) in the relevant degree regime but is still worse than the balanced half-list evaluator at b/2+o(b). There is no gain to the best full relation search. The full conditional contraction exponent b/(log b)^(1/4) remains unbeaten. A fast collective norm/trace evaluator specialized to the actual inputs is missing. Equation (5) shows why a proof that simply iterates the one-slot norm routine while asserting that all coefficient degrees stay small would be invalid.

The most productive next target would need to be a circuit-aware algorithm insensitive to that generic degree, or a genuinely different outer relation-generation method. More improvements in the local scalar readout do not establish the desired exponent.

## 7. Reproducibility

Run `python verify.py` (Python standard library). Tests use eight small auxiliary curves, including a subgroup of order 343; precision p,p^2,p^3; exact kernel evaluation, two projective charts, small quotient reductions, planted vanishing factors, rational traces including repeated poles, reversed-kernel power sums, and group-homomorphism checks. Actual torsion lifts are constructed with the odd division-polynomial recurrence and Hensel lifting. The test fixture search uses literal small point enumeration, not Schoof or the asymptotic L^2 divisibility setup; that setup is justified analytically in Section 1. The code's unit inverse routine is local Gaussian elimination and small determinant computation is subset DP; the general product-algebra circuit uses charged adjugate/determinant algorithms instead. No large-field timing or collective trace benchmark is claimed. The generic irreducibility and degree growth in Section 5 is proved algebraically, not asserted from a numerical test. `verification.txt` records the executed checks.
