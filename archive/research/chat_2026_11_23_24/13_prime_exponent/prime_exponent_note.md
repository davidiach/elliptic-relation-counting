# Direct elliptic evaluation-code route: transfer target and an oracle obstruction

Research checkpoint, 24 November 2026. Continuation of `paper (1).pdf` (v3), the section-norm compiler, and the exterior-packet checkpoint. This note does **not** beat the general-prime exponent. It supplies a precise direct-solver target, an actual elliptic planted-list rank-oracle lower bound, and a stronger obstruction to determinant zero tests of public linear combinations of rows. The planted-list model below is a new uniform signed-list experiment, **not** the original dyadic affine hash law or the later field-affine digit sampler. The arguments do not rule out algorithms reading numerical row entries, nonzero determinant values, or the actual hash seeds. No novelty claim against the literature is made.

## 1. Direct solver interface

Let p>3, b=ceil(log2 p), and let E/F_p be a smooth short Weierstrass elliptic curve with known odd prime subgroup order r and polynomial cofactor, as in the original v3 paper. For an eligible target T!=O, let V=L((m+1)O) have its usual Riemann--Roch monomial basis with the constant first. Subtract the evaluation at -T from the evaluation at P and remove the constant coordinate, obtaining v(P) in F_p^m. Assume sampled slots are regular: no abscissa overlaps and no retained point has the target abscissa. A transversal of m points has determinant zero iff its sum is T, by the divisor criterion. This criterion and the section-space formulation already occur in the original project; the reduction here is a change of computational interface, not a new elliptic identity.

If a direct algorithm for this particular dependent-transversal problem has expected bit time

    poly(b,m) w^C exp_2(h(m))

for a fixed C independent of m, where w is the nominal slot width, and returns a verified witness with probability at least 2^-O(m) at the paid sampler widths, then the inherited one-pilot/training/sparse-solve reduction gives log expected time O(b/m + m + h(m)), assuming its hypotheses hold at the chosen arity. Generation, repetitions, and linear algebra are included at their inherited 2^O(b/m+m) cost. Use the balanced field-affine sampler with constant targetwise coverage from the earlier affine-digit note, or another sampler with the same guarantees. Regularity can be checked after enumeration. The direct solver is not supplied by this checkpoint.

Examples: h(m)=O(m) gives m~sqrt(b) and O(sqrt(b)) log time; h(m)=O(m^2) gives m~b^(1/3) and O(b^(2/3)); h(m)<=2^(C' m) gives m=c log b for a sufficiently small c and O(b/log b). All are strict improvements over b/(log b)^(1/4). These are transfer calculations, not established algorithms. The generic colourful dependent-transversal problem is distinct from the polynomial-time independent-transversal problem (Bárány and Onn, *Colourful Linear Programming and Its Relatives*, Mathematics of Operations Research 22(3), 1997, Theorems 5.1 and 5.3). Their hardness theorem is over Q and does not settle this parameterized finite-field elliptic subclass.

## 2. Actual elliptic planted-list experiment

For this section restrict to a prime-order curve group G=E(F_p) of odd prime order r (hence log2 r=b+O(1)), and fixed T!=O. Let m>=3 and w>=1. There are n=mw underlying labeled points A_ij, j=1,...,w in each slot i. Each underlying point supplies two labeled rows for +A_ij and -A_ij. An address alpha selects exactly one labeled signed row per slot, so W=(2w)^m possible addresses.

Choose alpha uniformly. Sample A uniformly conditioned on the signed sum selected by alpha being T. Equivalently draw n-1 independent uniform group points and force the final selected point. This is a simple planted distribution of uniform signed lists; it is not the affine sampler. Nonidentity points have the expected rational evaluations v(P) defined in Section 1. For the proof, call the instance regular if all 2n signed rows correspond to distinct affine points avoiding ±T.

Each elementary irregularity (A_j=O, A_j=±T, A_j=±A_k for j!=k) is a linear constraint independent of the plant in the abstract cyclic prime-order group, because the plant uses m>=3 distinct underlying labels. Its conditional probability is 1/r. Thus

    Pr(not regular) <= (3n+2 binom(n,2))/r.

Every m-element subset of the 2n signed row labels defines a signed linear form c.A in the underlying group variables with c_j in {-1,0,1}, even when it includes both signs of one point. If it is not the planted subset, its coefficient vector cannot equal the planted vector (the latter already uses m nonzero coordinates). If c is proportional to the plant, it is either zero or the negative, and then its sum cannot be T since r is odd and T!=O. Otherwise its form is independent of the plant, so the conditional probability of c.A=T is 1/r. Union bound gives

    Pr(an extra m-subset sums to T) <= (binom(2n,m)-1)/r.

Set epsilon to the sum of these two bounds (capped at one). On the complementary good event, the actual evaluation matroid on the signed rows is the uniform rank-m matroid with exactly one basis removed: the planted subset. Indeed, for k<m distinct points avoiding -T, the evaluation conditions on L((m+1)O) at those k points and at -T are independent. Riemann--Roch gives dimension m-k for the kernel, since the residual divisor has positive degree m-k on a genus-one curve. Thus the projected rows have rank k. For k=m, dependence is exactly the sum-T criterion. Any larger subset contains an independent m-subset, since only one m-subset is dependent.

For any randomized adaptive algorithm given only the labels, target, and an oracle returning the rank of *any* subset, follow the answers from the uniform matroid. Its at most J rank queries on that path and its final output (if a proposed transversal) form a fixed list of at most J+1 candidate addresses once its coins are fixed. On the good event, an answer can differ from the uniform answer only when the oracle is asked for exactly the planted m-subset. Therefore its probability of recovering a relation under the planted experiment is at most

    (J+1)/W + epsilon.

This union bound is unconditional; it does not require the hidden address to remain uniform after conditioning on regularity and uniqueness. Fixing coins and averaging proves the randomized bound.

At w=floor(r^(1/m)/16), with m->infinity, m=O(sqrt(b)), we have epsilon=o(1), using binom(2mw,m)/r <= (2e/16)^m and n^2/r=o(1). Also log2 W=log2 r-O(m), since w->infinity. A constant-success rank-only solver therefore needs 2^(b-O(m)) queries. This applies to actual elliptic evaluation matroids in this planted uniform-list model, including the ability to query ranks of nontransversal subsets; it is not an oracle lower bound on all ordinary arithmetic access to the curve.

## 3. Stronger obstruction: determinant zero tests of linear combinations

Retain the uniform signed-list experiment. An algorithm can specify arbitrary public matrices C in F_p^(k x 2n) (and choose k columns of the m-coordinate rows, for 1<=k<=m), and ask whether the corresponding determinant of C times the selected-column evaluation matrix is zero. Coefficients can depend on E, T, labels, earlier zero-test answers and private coins, but cannot depend on numerical row values or other access to the hidden points. Include the case k=m, the full determinant of arbitrary linear combinations. This is substantially more permissive than testing single transversals, but less permissive than reading nonzero determinants.

Write d=m+1, D=2d. For a fixed test, its determinant g is a rational function on E^n over F_p. By Cauchy--Binet, g is a linear combination of minors of the underlying evaluation rows. Each term uses distinct signed row labels, so a single underlying variable occurs at most twice in any term. Each evaluation coordinate has pole order at most d at O; replacing P by -P does not change the pole location. Thus g has pole order at most D along each coordinate identity divisor, and no other poles. This also applies after selecting fewer columns and rows. If g is identically zero on E^n it provides the zero answer on every regular input and can be ignored.

For a nonzero g, consider the algebraic planted divisors

    D_alpha = {A in E^n : sum_i epsilon_i A_(i,j_i) = T}.

Each is irreducible and isomorphic to E^(n-1), via solving for any selected coordinate. For distinct addresses they are distinct: when they both involve a given coordinate, their graphs as functions of the other coordinates have different constants (±T) or different nonzero coefficient maps (±identity); 2T!=O and characteristic >3 suffice. For a generic fiber varying coordinate j and fixing all others, every divisor containing j cuts out a single point, and these points are distinct for distinct addresses involving j. If g vanishes identically on a planted divisor then its restriction to that generic coordinate fiber vanishes at that point. A nonzero rational function with sole pole of order <=D has at most D zeros counting multiplicity. Hence at most D exceptional planted divisors involve j. Sum over n coordinates and count each exceptional divisor m times:

    #{alpha : g|D_alpha identically zero} <= nD/m = wD.

For a nonexceptional divisor, parameterize D_alpha by n-1 independent group points. The restriction of g is nonzero. In a univariate slice with other coordinates fixed away from wholly forbidden fibers, there are at most two possible pole locations (the variable's own identity and the forced coordinate's identity), with total pole order <=2D. Excluding the corresponding undefined points and zeros costs at most 2D+2 points for a nonzero slice.

For clarity, the elementary product-grid bound used here is as follows. Induct on the number of independent coordinates. Choose a generic algebraic specialization of the final coordinate where the remaining rational function is nonzero and all identically forbidden specializations are avoided. Any fixed remaining tuple for which the final-coordinate slice is identically zero or everywhere forbidden is bad in this specialization. Induction bounds its probability by (k-1)s/r. For every other tuple the univariate zero-plus-forbidden bound is s/r. This works over an algebraic extension containing the specialization: the sampled coordinates still lie in the original r rational points, which are distinct in that extension. At each stage the two possible pole loci come from the identities of an unforced point and the forced sum point, and a generic specialization prevents a globally constant forbidden identity. When a previously fixed coordinate is O, its fiber is counted as forbidden, not evaluated. With k=n-1 and s=2D+2, the bound is (n-1)(2D+2)/r.

Consequently for every fixed nontrivial zero test,

    Pr_plant(g=0 or undefined) <= wD/W + (n-1)(2D+2)/r.    (1)

Under the fresh independent uniform-list law, the analogous bound is n(D+1)/r, with only one identity pole per coordinate. The bounds can be capped at one; they are not designed to be tight on tiny fields.

For an adaptive J-query algorithm with only this zero-test interface, define the generic path by zero for identically zero tests and nonzero for all other tests (this path is a theoretical object, not a routine the algorithm needs to compute). Fix its coins. A departure from the generic path on a regular input requires a nontrivial test along that fixed path to vanish. By (1) and the union bound the probability is at most J[wD/W+(n-1)(2D+2)/r]. On the generic path the proposed address is fixed, so its chance to equal the uniformly planted address is at most 1/W, unconditionally. Recovering any different valid relation or encountering an irregular input is already accounted for by epsilon. Thus success is at most

    1/W + J [wD/W + (n-1)(2D+2)/r] + epsilon.             (2)

Averaging coins preserves (2). At the thinned widths above, and m->infinity, m=O(sqrt(b)), the bracket is 2^(-b+b/m+O(m+log m+log b)). In particular J=2^o(b) gives success o(1). This rules out a generic-path sequence of determinant equality tests as the desired constant-success direct solver in this model. It does not bound full-value finite-field arithmetic, algorithms with input-dependent numerical coefficients, approximate real readouts, or the original affine seed law. The proof concerns the underlying point variables and their elliptic sum fibers, not the earlier fixed affine-seed polynomial obstruction.

## 4. Reproduction

Run `python -B verify.py` (standard library). The code checks the prime group order by explicit enumeration on the toy curve y^2=x^3+44x+1 over F_65521, draws uniform planted signed lists, compares determinants with independent group addition for every m-subset of the regular test lists, and independently verifies Cauchy--Binet for public dense linear-combination matrices. See `verification.txt` for the actual transcript. These finite checks exercise the algebraic interfaces; the asymptotic bounds are the proofs in Sections 2 and 3. No fast witness finder or full ECDLP benchmark is claimed.

## 5. Status

The previous conditional numerical candidate 11.31486814 is unchanged and not independently audited. The strict general-prime exponent target remains open. The direct route must use numerical arithmetic or richer joint access to the actual elliptic rows; ordinary matroid intersection solves the opposite problem (finding an independent transversal), and generic determinant zero tests do not supply the missing relation evaluator. The promising next mathematical target is a uniform direct solver using *numerical* row structure at time w^O(1) exp_2(exp(O(m))), or any weaker bound that still yields o(b/(log b)^(1/4)) log time after optimizing arity. No existence claim for such a solver is made.

Sources: `paper (1).pdf`, Sections 2 and 5; `affine_digit_sampler_note.md`; project section-space derivations; I. Bárány and S. Onn, *Colourful Linear Programming and Its Relatives*, MOR 22(3) (1997), 513--768, https://pubsonline.informs.org/doi/10.1287/moor.22.3.550. Only the distinction between generic independent and dependent transversal algorithms is attributed to that external reference.
