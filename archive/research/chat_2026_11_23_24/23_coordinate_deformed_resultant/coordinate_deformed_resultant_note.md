# Coordinate-deformed elliptic resultants: signed recovery and a prime-power evaluator

Research checkpoint continuing `mixed_characteristic_norm_note.md` and `shared_norm_aggregate_note.md`. **The general-prime asymptotic ECDLP exponent is not improved.** The new result is a stronger evaluator interface: keep the original elliptic resultant and perturb its input coordinates, instead of evaluating the norm of an arbitrary perturbation of the multivariate summation polynomial. This yields exact recovery of signed points, even for several algebraic hits. A homogeneous resultant/remainder-tree construction evaluates the scalar at the old half-list scale, also over the required prime-power rings. The earlier candidate coefficient 11.31486814 and the sampler-to-ECDLP reduction are not re-audited here. No novelty or research-priority claim is made.

## 1. Setup and declared-degree resultants

Let E/F_p be the nonsingular short Weierstrass curve y²=f(x)=x³+ax+c, p>3. Let h be a positive integer prime to p (in particular a polynomially bounded cofactor for sufficiently large p). Let T=(x_T,y_T) be rational with y_T nonzero. For finite public sets S_i ⊂ F_p, put x_i(s)=A_i+t_i s, t_i nonzero. Assume f(x_i(s)) is nonzero for every address. We may use independent A_i or the shared-translation sampler; this checkpoint is an algebraic result for each fixed grid, not a new sampling theorem. Let W=∏|S_i| and let H count **oriented algebraic signed tuples** with [h]∑P_i=T, x(P_i)=x_i(s_i). The y_i belong to F_{p²}, and a signed tuple with projected sum -T gives the same global-sign orbit as one with sum T. Each such orbit has exactly one orientation to T, since T≠-T. H can exceed the number of distinct abscissa addresses.

Use the homogeneous degree-two-per-coordinate version of

    S3(u,v,X)=(u-v)^2 X² -2((u+v)(uv+a)+2c)X
              +(uv-a)^2-4c(u+v).

For fixed x(P),x(Q) in P¹, S3 as a binary form in X has the two roots x(P+Q), x(P-Q), counted with multiplicity. The assertion includes O, with x(O)=∞. It is crucial that this binary form is nonzero at every fixed pair of coordinates on the smooth curve, including the special cases. For finite u≠v its leading coefficient is nonzero; for u=v with f(u)≠0 the linear coefficient is -4f(u); for u=v with f(u)=0 the constant coefficient is f'(u)^2≠0. At one infinity its form is (X-u)^2 in homogeneous notation, and at two infinities its only root is infinity, twice. The usual chord and tangent formulas verify the roots in all cases.

For a single tuple of leaf abscissas, recursively form the binary polynomial F_I for any block I using leaf forms X_0-x_i X_1 and homogeneous resultant joins with S3. Always use **declared homogeneous degrees**, not a degree that changes after specialization. Its degree in X is d_I=2^(|I|-1). For projection by [h], write its induced Kummer-line map as (U_h:V_h), homogeneous degree h², with no common zero on the special fiber. Standard division-polynomial recurrences construct this map in poly(h) field/ring operations; it is the quotient of the multiplication morphism by negation. For a finite target coordinate X, the additional resultant against U_h-X V_h has degree d_I in X. One may equivalently perform projection inside the final pair join. Work with raw homogeneous resultants rather than canceling a factor that might vanish modulo p.

**Local factorization lemma.** Lift the smooth curve, leaf abscissas and target coordinate to a complete unramified DVR of residue characteristic p containing the ordinates in its residue field. The nonzero leaf ordinates lift uniquely after a sign is specified (Hensel, since 2y_i is a unit). Represent the x-coordinate of any resulting group sum by a primitive projective pair (a:b), meaning at least one of a,b is a unit. The raw recursively formed F_I is, up to a unit, the product of the binary linear forms bX_0-aX_1 of its 2^(|I|-1) signed orbits. Projection sends each of these roots to x([h]P), retaining multiplicities.

Proof: at the fraction-field level S3 is a scalar times the two appropriate linear forms. Both the S3 form and the product of primitive forms are primitive modulo p, by the cases above. Hence the scalar is a unit. Homogeneous resultant multiplicativity and evaluation against primitive linear factors then prove the assertion inductively. Equivalently one can use the polynomial identities first on generic inputs and specialize. The projection form has the same property because (U_h,V_h) has no base point. For two blocks, one relative sign ± is chosen for each pair of block sign orbits, giving exactly the parent sign orbits even when two sums happen to have the same x-coordinate. This reasoning does not discard partial sums at infinity or torsion points. ∎

Denote by G the product of the projected tuple forms, evaluated at (X,1). It is the ordinary collective resultant of the earlier checkpoint (a balanced rather than sequential join gives the same local product with a unit normalization). In all results below arbitrary unit factors are allowed: at the leading p-adic valuation only their reductions matter, and for first derivatives at a zero their differentiated part vanishes. No assumption about an explicitly factorized representation of G is made for a hypothetical fast evaluator.

## 2. Exact first derivative recovery, over the original field

**Theorem 1.** If H=1, the partial derivatives of G in the slot parameters and target recover the unique oriented signed tuple:

    s_i = G_{t_i}/G_{A_i},
    y_i = -h y_T G_X/G_{A_i}.                       (1)

All derivatives are evaluated at the base instance. Every denominator is nonzero. Moreover a unique algebraic oriented hit is rational: Frobenius preserves its address tuple and target, so by uniqueness it fixes every ordinate.

Proof. Over the residue-field extension, each leaf x-coordinate is a local parameter because its ordinate is nonzero. The addition and multiplication morphisms satisfy the invariant differential identity

    ([h]∑P_i)^*(dx/(2y)) = h ∑ dx_i/(2y_i).

On the oriented witness branch this gives ∂X_sum/∂x_i=h y_T/y_i. The unique vanishing factor in G can locally be normalized to X-X_sum. Every other factor is a unit at the point, including the factors corresponding to sums at infinity. Thus there is a common nonzero scalar C with

    G_X=C,  G_{A_i}=-C h y_T/y_i,
    G_{t_i}=-C s_i h y_T/y_i.

The claimed identities follow. The product-unit assertion also shows G vanishes exactly on the algebraic incidence event, and when H>1 all its first derivatives vanish. ∎

**Circuit transfer.** A division-free polynomial-identity circuit for the ordinary collective resultant as a function of (A_i,t_i,X), with S_i fixed, gives all the derivatives by a single reverse differentiation pass at constant-factor arithmetic overhead (plus charged input/output access). The reverse tape can require space proportional to the circuit size; no polynomial-space claim follows automatically. An algorithm that merely agrees on scalar base-field inputs is not automatically such a circuit. Unlike the previous checkpoint, we need no arbitrary perturbation of the internal summation polynomial coefficients. Independent infinitesimal shifts can be introduced even when A_i are equal at the base instance. Differentiating only a common A instead of allowing independent A_i is not sufficient for (1).

## 3. One-call mixed-characteristic recovery of all signed hits

This is the main new scalar readout. It preserves every leaf factor and every homogeneous resultant join.

Choose K≥1 as a cap on H. Work with an unramified coefficient ring modulo p^(K+1) whose residue field has a tower

    k0 ⊃ F_{p²},
    k=k0(v), [k:k0]=2m+1,
    L=k(θ), [L:k]=K+1.

For the stated cost bounds k0=F_{p²} suffices, or one may use a slightly larger fixed even extension for projective charts. Lift the defining monic irreducibles to construct the tower over Z/p^(K+1)Z, as in the preceding mixed-characteristic checkpoint. Choose canonical integer lifts of base parameters and addresses, and deform each leaf and the target by

    x_i^+(s)=A_i+t_i s+p(v^(2i-1)+s v^(2i)),
    X^+=x_T+p θ.                                    (2)

All the original grid polynomials remain products of linear factors in their own coordinates. The deformed grid is still a Cartesian grid of affine slot parameterizations. In particular, this does not add a term to F, or square F, which would generally change the nested resultant representation.

**Theorem 2.** Evaluate the single scalar G^+ at (2) over this tower ring. If G^+=0 at the chosen precision then H>K. Otherwise the p-adic valuation is H≤K and

    (G^+/p^H) mod p = C ∏_{(s,y):[h]∑P_i=T}(θ-e_{s,y}),  C≠0,
    e_{s,y}=-c_{s,y}+∑ (h y_T/y_i)(v^(2i-1)+s_i v^(2i)),  (3)

where c_{s,y}∈k0 is an unknown base-lift discrepancy and C is a unit reduction independent of v,θ. The roots e_{s,y} are distinct. Read the θ coefficients of (3), normalize by its leading coefficient, factor the degree-H monic polynomial over k, and recover all signed tuples exactly. In the v basis, if u_i,z_i are the coefficients of v^(2i-1),v^(2i), then

    s_i=z_i/u_i,  y_i=h y_T/u_i.                     (4)

Proof. Work first over an unramified DVR lifting the tower; reduce to the stated precision at the end. For a positive oriented branch, its baseline projected point Q^0 reduces to T. The difference between its baseline x-coordinate and the chosen lift of x_T is divisible by p. Write c_{s,y}=(x_T^lift-x(Q^0))/p mod p. Both baseline lifts are defined over the unramified lift of k0, so c_{s,y} lies in k0. The first differential of the group morphism gives

    (X^+-x(Q^+))/p mod p
      = θ+c_{s,y}-∑ (h y_T/y_i)(v^(2i-1)+s_i v^(2i)).

This is θ-e_{s,y}; it is a unit because θ is outside k. A factor corresponding to a negative branch, including a branch at infinity, is a unit without division. The local factorization lemma then proves the valuation and (3). All baseline unit reductions are independent of v,θ; the p-adic perturbations vanish modulo p. The coefficients u_i are nonzero by p∤h and the hypotheses on the target and leaves. Formula (4) proves injectivity even for two distinct sign choices at the same address, hence the roots in (3) are distinct. Since H≤K, neither the θ-polynomial nor the linear-in-basis encoding of its individual roots loses information under extension reduction. If H≥K+1, the product vanishes modulo p^(K+1), and conversely if H≤K its normalized product is nonzero. ∎

The ring has O(b m(K+1)²) bits per scalar for a fixed degree k0. Irreducible generation, ordinary field factoring of a split squarefree polynomial and decoding take expected polynomial time in b,m,K. Verify recovered rational signed points and projected sums; the encoding also correctly returns nonrational algebraic hits, which must not be treated as rational relations. This checkpoint does not establish a new bound on their frequency. The preceding shared-translation notes have separate coverage/regularity assumptions which have not been re-audited here.

**Singleton simplification.** Take one extension of degree 2m+2 over F_p modulo p². In (2) replace θ by v^(2m+1). On H=1 the unique hit is rational, so the baseline discrepancy and all derivative coefficients lie in F_p. Divide G^+ by p, reduce mod p and normalize by the highest v coefficient. The coefficient of v^(2i-1) is -h y_T/y_i, and the coefficient of v^(2i) is s_i times that. Recover addresses by their ratio and ordinates by -h y_T divided by the odd-index coefficient. H=0 is a unit and H≥2 gives zero. No factoring, random projection or arbitrary polynomial-norm input is needed. This is one scalar call to the original resultant over an O(bm)-bit scalar ring. Formal equal-characteristic jets can be used instead of p-adic shifts for a polynomial-identity circuit; their interface remains distinct from a base-field black-box evaluator.

## 4. An exact scalar evaluator with charged half-list costs

The preceding mixed-characteristic norm of f²+p·tag was not preserved by elliptic resultant composition. The input perturbation (2) is preserved. We give a scalar algorithm at the already-known half-list time scale; it does not improve the full ECDLP asymptotic bound. It is useful because it computes the *same* structured resultant over the coefficient ring needed for signed recovery and does not need to evaluate arbitrary perturbed multiplication polynomials.

Split the slots into I,J. Put d_I=2^(|I|-1), d_J=2^(|J|-1), W_I=∏_{i∈I}w_i, W_J=∏_{i∈J}w_i. For each tuple in J form its unprojected homogeneous relation polynomial F_{s_J}(Z) of degree d_J. Define Q(Z)=∏ F_{s_J}(Z), declared degree D=d_J W_J. Form for each tuple s_I the small binary polynomial

    B_{s_I}(Z)=Res_Y(F_{s_I}(Y), C_h(Y,Z,X^+)),
    C_h(Y,Z,X)=Res_V(S3(Y,Z,V), U_h(V)-X V_h(V)).

Use C_1=S3 directly (up to harmless sign). The degree of B_{s_I} is declared e_I=2h²d_I. The B forms are primitive modulo p: fixing any partial point Y, a given finite projected target restricts the other point to a finite set, since [h] and translation are finite/nonconstant morphisms. Thus the projected pair form is not identically zero as a form in Z; the factorization lemma and resultant multiplicativity give the assertion for B. The scalar

    G_split=∏_{s_I} Res_Z(Q,B_{s_I})                  (5)

has precisely the factorization of §1, including multiplicities and units, so it is an admissible ordinary collective resultant for Theorems 1–2. This remains true at all prescribed deformations. We use (5) as the definition of the evaluated scalar; equality to another normalization of G is not needed.

The impediment to simultaneous reductions Q mod B is that the specialized leading coefficient may be a nonunit (including when partial sums at infinity are present). A public projective chart fixes it:

    (Z_0,Z_1)=(αz+1,z).

Choose α in k0 uniformly until all B_{s_I}(α,1) are nonzero modulo p. If |k0|≥2 e_I W_I, the union of their root sets occupies at most half of k0, so expected two attempts suffice. Test this condition explicitly; choose a larger residue extension if necessary. At the balanced widths and m≤sqrt(b), polynomial h, a fixed extension (e.g. F_{p²} for sufficiently large b) suffices. The perturbations vanish mod p, so chart tests depend only on the base instance. Every B' in the new chart has leading coefficient λ a unit; divide by λ to get a monic b. Form Q' by multiplying the transformed small F_{s_J} forms (no expensive large-degree change of coordinates is required). Compute all r_{s_I}=Q' mod b_{s_I} by a monic subproduct/remainder tree. Pairwise coprimality is not required, and monic division works over the prime-power extension ring with zero divisors.

The exact homogeneous resultant identity is

    Res_Z(Q,B) = λ^D Res_z(b,r).                    (6)

Indeed under the chart matrix of determinant -1 the resultant gains the factor (-1)^(D e_I); swapping Q',B' gains the same sign, so they cancel. Multiplying b by its leading coefficient contributes λ^D. Finally the root product identity for a monic b permits reduction of Q' modulo b, including degree drops in the remainder. Formula (6) is a polynomial identity after the allowed unit normalization. All signs and scalar powers are accounted for; multiplying the right sides computes (5) exactly. One can also use their unit-associate version for witness readout.

**Cost.** Each tuple polynomial has separate degrees bounded by 2^{O(m)}poly(h). Compute its coefficients by small division-free determinants (Berkowitz works over arbitrary commutative rings) and coefficient arithmetic or unit-separated interpolation. For arbitrary finite ring size, extend the residue field if necessary to obtain enough distinct interpolation residues. This charges 2^{O(m)}poly(b,h,K) bit work per tuple. Use fast polynomial multiplication and monic division for the subproduct and remainder trees. Monic division can be implemented via the reciprocal of a reversed monic polynomial's unit-constant power series. Fast integer convolution on the lifted integer coefficient arrays implements fast polynomial multiplication in these finite free prime-power rings with polynomial overhead in the extension degree, precision, and log degree. The trees take soft-O(D+e_I W_I) coefficient operations. The last small resultants can be evaluated division-free at cost polynomial in e_I and the coefficient size per left tuple. Chart search, unit inversions, λ^D powers and output readout are charged. Therefore

    time ≤ (W_I+W_J) 2^{O(m)} poly(b,h,K,log W)

up to logarithmic factors, including coefficient-ring construction. Space at this scale suffices (we make no stronger space claim). Choosing the split by log widths gives

    log_2 time ≤ (1/2) log_2 W + (1/2) max_i log_2 w_i
                 + O(m+log b+log h+log(K+1))

with suppressed logarithmic factors absorbed. At balanced widths with log W=b+o(b), m growing, m≤sqrt(b), polynomial h and log K=O(m+log b), this is b/2+o(b). It is a genuine improvement to the cost of evaluating the **coordinate-deformed scalar** compared to the explicit perturbed-product calculation in the previous checkpoint. Ordinary meet-in-the-middle relation search already has this scale; it is not a new ECDLP complexity bound. The remaining need is to evaluate (5) without materializing either half's family of polynomials or an equivalently large representation.

## 5. Verification and reproduction

Run `python -B verify.py` here. The two previous verifier files are copied from the preceding checkpoint for standard-library ring arithmetic and small polynomial helpers. Their previous test results are not counted as new results. The new script independently enumerates signed points in F_{p²}, computes actual group sums and compares them with declared-degree homogeneous resultants over finite rings. Nonrational branches are included. The projection factor tested is h=2 (not assumed to be an actual cofactor on each toy curve); the proof covers p∤h. The tests deliberately include empty, singleton, multi-hit and overflow fixtures, selected by the exact truth counts; these are not random-sampler performance estimates.

The final transcript reports 35 elliptic grids, 435 addresses and 3,164 oriented sign choices; 34 sums at infinity before projection and 48 after projection; 14 projected grids. It verifies 35 first-jet calculations including 60 singleton address/ordinate comparisons, 35 one-extension scalar calculations including 12 exact singleton recoveries, and 21 two-hit-cap tower calculations, with 16 recovered signed tuples (including two nonrational tuples) and one overflow. A separate 14 small-grid and eight four-slot balanced-grid remainder-tree suite verifies 59 monic remainders against direct division and 217 homogeneous pairwise resultants against the exact scalar product. Five of these grids have nonunit leading coefficients before chart selection. Counts refer to actual executed tests, not independent asymptotic evidence. All assertions passed.

The tests use schoolbook polynomial arithmetic and a small division-free subset-determinant implementation. The asymptotic cost proof uses fast polynomial arithmetic and polynomial-cost division-free determinants; we do not claim benchmark validation of these fast algorithms or the complete ECDLP reduction. The previous mixed-characteristic verifier's input-ring construction is reused but its previous tests are not rerun by import.

References for standard primitives: the previous `shared_norm_aggregate_note.md` derives the elliptic collective resultant from the classical Semaev polynomial; `mixed_characteristic_norm_note.md` derives the scalar readout for arbitrary norm perturbations. For division-free determinants over commutative rings, see S. J. Berkowitz, *On computing the determinant in small parallel time using a small number of processors* (1984), https://cs.uwaterloo.ca/~r5olivei/courses/2021-winter-cs487/berkowitz-paper.pdf. For multipoint evaluation over general rings with monic subproduct trees, see Bostan–Gaudry–Schost, https://mathexp.eu/bostan/publications/BoGaSc04.pdf, and the generic subproduct-tree discussion at https://fredrikj.net/fdoc/html/algorithms/polynomials/subproduct-tree.html. Our proofs use the monic remainder identity directly; no assumption about fast evaluation in the original exponentially large grid algebra is smuggled into the local cost.

## 6. Frontier

The earlier generic mixed-characteristic readout is now specialized to the ordinary elliptic resultant. It recovers signed hits without arbitrary perturbations of F, either by derivatives or a single scalar coordinate-deformation call. The half-list scalar evaluator gives a concrete upper bound for that primitive at prime-power precision. A smaller asymptotic evaluator is still missing: construction of the half families dominates the charged time. In particular b/2+o(b) does not establish o(b/(log b)^(1/4)). No improved general-prime coefficient or asymptotic exponent is claimed. The most relevant new gate is to evaluate the ordinary nested resultant on these deformed inputs while bypassing the explicit half-family materialization, or exploit the derivative ratios without evaluating its entire scalar value. Low p-adic precision alone has not supplied such a method.
