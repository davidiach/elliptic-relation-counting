# Global elliptic section norms and pairwise compression of dyadic modes

24 November 2026. Continuation of the v3 paper and the preceding auxiliary-polynomial, abscissa-chain, short-projector, affine-digit-sampler, and contractor-degree notes.

**Result.** A new regular-support relation compiler has leading Boolean word budget `(20/3+o(1)) m b`. Under the same quantitative finite-packet and outer-reduction estimates used in the preceding numerical candidates, its time coefficient is 11.31486813704766 and its space coefficient is 5.65743406852383, both in `b/(log2 b)^(1/4)`. The earlier 8mb budget gives 12.39481702771492. This is a candidate improvement of the numerical coefficient, not an improved power of log b. The second result is an exact pairwise-mode tensorization of the public Boolean contractor's crossing sum. The straightforward two-stage packing analysis cancels its potential extra logarithmic saving; no stronger generic contraction theorem is claimed.

The finite regression below verifies the norm/divisor equivalence and the projected polynomial residual identities, but does not implement the full contractor or certify a cryptographic-scale timing. Novelty against all unpublished literature has not been established.

## 1. Reused polynomial-ring interface

Let `b=ceil(log2 p)`, `L=log2 b`, and work on the original smooth elliptic curve `y^2=x^3+a x+c`, `p>3`. Use the preceding polynomial-ring construction with `R=Z[Z]/(Z^n-lambda)`, `lambda in {1,2}`, `n=q1 q2=Theta(b/L^2)`, two distinct balanced odd primes, and a field evaluation `phi: Z -> gamma` into F_p. Write `d=n` in the binomial case and `d=phi_Euler(n)` in the cyclotomic case, and set `k=b/d=Theta(L^2)`. Thus `n/d=1+o(1)`. In the first case `X^n-2` is irreducible over Q and vanishes at gamma modulo p; in the second case use `Phi_n` of degree d and gamma of exact order n. These are the irreducible polynomials used for the resultant estimates below. The dimension choices cover arbitrary field primes as in the previous notes.

Retain the preceding strictly diagonally dominant polynomial M in ker(phi), with multiplication matrix A_M, determinant absolute value D=p c0 with p not dividing c0, and a short public projector E satisfying

    phi(F)=0  iff  E F is in (M).

The previous construction has full field representative coefficient bound `S=2^(k+O(log b))`, projector coefficients `2^O(log b)`, and `||A_M^(-1)||_infinity=O(1/M_0)`, where M_0 is the positive constant coefficient. In particular the quotient Q of EF by M, when it exists, has an explicit uniform magnitude bound. Reused preprocessing costs `2^O(n log n)`. See `abscissa_chain_note.md` and `large_arity_research_note.md` for these earlier arguments. Nothing in the new norm compiler requires canonical representatives of its private field coefficients: bounded representatives suffice for soundness, and a surjective box suffices for completeness. The projector is needed for arbitrary representatives.

### 1.1 A short-ratio lemma and a short nonsquare

For each nonzero `v in F_p`, there exist U,V of degree <d with

    phi(U)=v phi(V),   phi(V) != 0,
    log2 max(1,||U||_infinity,||V||_infinity) <= k/2+O(log b).

Indeed the coefficient lattice of pairs with `phi(U)-v phi(V)=0` has rank 2d and index p (the coefficient of the constant term of U is 1). Minkowski's theorem or a pigeonhole argument gives a nonzero pair of the indicated size, with an additional polynomial-in-d factor allowed. A nonzero polynomial g of degree <d of this size cannot evaluate to zero: its nonzero integer resultant with the irreducible polynomial just specified has absolute value bounded by `(2d ||g||_infinity)^d = p^(1/2+o(1)) < p`. The resultant would be divisible by p if `phi(g)=0`. Therefore V and U both have nonzero images. A sufficiently short pair can be found by exact shortest-vector search in `2^O(d log d)` bit time (standard exact lattice algorithms, with polynomial factors in the input bit size); the finite code uses checked LLL, not an assumption that LLL gives the asymptotic bound.

Choose a nonsquare v by ordinary random sampling and Euler testing, expected polynomial time. The above ratio forces one of U,V to have nonsquare image. This produces a short representative of a nonsquare without assuming a small integer quadratic nonresidue or GRH. Hence either square class has a representative with coefficient bound `2^(k/2+O(log b))` (use 1 for the square class). This argument also supplies a direct injectivity certificate for thin digit cubes.

## 2. Shrink the three weights of the elliptic norm

If `ac != 0`, change only the abscissa by `x=u X`, where `u=c/a`. The curve is

    y^2 = f_u(X)= A3 X^3 + C(X+1),
    A3=u^3, C=c.

For `a=0` or `c=0`, use `u=1`; the cubic then still has at most two distinct nonzero coefficients. The curve is smooth, so the nonzero coefficients needed below exist. Denote its two nonzero weight values by c1,c2 (allow equality). Apply the ratio lemma to c2/c1 to find short U,V, and set `h=phi(V)/c1`. The coefficients of `h f_u` have short representatives: c1 maps to V, c2 to U, zero to zero. This uses an exact field-affine change of abscissa and does not replace the rational group with its twist.

Choose a short D0 with image in the square class of h and choose nonzero `s` with `h s^2=phi(D0)`. All weights in the norm of the section `s A(X)+y B(X)` can now be represented by coefficient polynomials of bit size `k/2+O(log b)`:

    h (s^2 A(X)^2 - f_u(X) B(X)^2)
      = phi(D0) A(X)^2 - (h f_u(X)) B(X)^2.

The square-root computation is in F_p, expected polynomial bit time. No product of the full-size scalar h with a short polynomial is transmitted. Instead each coefficient of the scaled cubic and h s^2 is directly represented by its already short polynomial.

## 3. The section-norm relation predicate

Let m be the number of slots, `d_s=m+1` and assume the retained abscissa sets are disjoint and avoid the target abscissa. This regularity is checkable before contraction by enumerating and sorting the factor base. For a proposed tuple of abscissas X_1,...,X_m and target T=(x_T,y_T), set X_T=x_T/u. Use

    deg A <= floor(d_s/2),
    deg B <= floor((d_s-3)/2),

and fix the coefficient of the unique basis monomial with pole order d_s to one (in A when d_s is even and in B when odd). There are d_s-1=m remaining field coefficients. Let c_lead be `h s^2` in the even case and `-h A3` in the odd case, where A3 is the nonzero cubic coefficient. Impose the polynomial equality in F_p[X]

    h s^2 A(X)^2 - h f_u(X) B(X)^2
        = c_lead (X-X_T) product_i (X-X_i),

plus the one linear field equality `s A(X_T)-y_T B(X_T)=0`.

**Soundness.** The right-hand side has d_s distinct roots and nonzero leading coefficient. At every such root B is nonzero, since common roots of A and B give a multiple root of the norm. Thus `y=-s A(X)/B(X)` is a rational ordinate and the point is on the original curve after scaling x back by u. The linear equality fixes the target point of the divisor to -T. The section has exactly pole order d_s at the identity, since the chosen leading coefficient is nonzero and the x-only and y-times-x monomials have distinct pole orders. All its d_s affine zeros have just been exhibited; the elliptic divisor criterion gives sum zero. Therefore the slot points sum to T. This reasoning handles rational two-torsion without extra pins and excludes any non-rational lift on a regular grid.

**Completeness.** A rational tuple summing to T gives an effective divisor of degree d_s whose sum is zero after adjoining -T. Its difference from d_s times the identity divisor is principal; by the usual elliptic-curve Riemann-Roch basis the corresponding function is a section of the indicated form. Normalize the highest-pole coefficient to one. Its norm has the d_s distinct specified roots and leading coefficient c_lead, and the target equality holds. Thus the test is exact for existence on regular inputs. It need not count distinct signed tuples when arbitrary coefficient representatives are allowed. The graph counts certificates, and a nonzero certificate count is all that witness extraction uses.

## 4. Injective thin support with an adjustable owner fraction

Take the balanced two-stage transform of length n as before. Choose a fixed fraction `theta` of the first-stage owners, `0<theta<=1`, and use their coefficient positions among 0,...,d-1. There are `(theta+o(1)) n` such positions. Put `w=ceil(p^(1/m))` and distribute `ell=ceil(log2 w)` free bits as equally as possible across those positions; set all other coefficients to zero. Each variable digit is a nonnegative coefficient of at most

    k/(theta m)+O(1)+o(k/m)

bits. Use the first w words of this cube. For m tending to infinity, the resultant argument of Section 1.1 shows that the evaluation map on the full cube is injective: the difference polynomial has degree <d and resultant absolute value at most `(2d*2^ceil(ell/number_positions))^d = p^(1/(theta m)+o(1)) < p`. The first w words therefore have exactly w distinct field images, without searching for a large-image block. An owner-ordered prefix comparison uses O(number of owners) constant-size state wires, a lower-order charge. The enumeration cost is `2^(b/m+O(1)) poly(b,m)`.

Independently in each slot choose uniform `a_i in F_p`, `t_i in F_p^*`, and accept `X_i=a_i+t_i phi(U_i)` for the thin variable polynomial U_i. Retain both rational ordinates. The abscissa map X->x=uX preserves affine two-point transitivity. The previous affine sampler's first- and second-moment estimates apply verbatim: one-point retention is w/p, two distinct abscissas have joint retention `w(w-1)/(p(p-1)) <= (w/p)^2`. With w as above, `lambda=N^(m-1)(w/p)^m=1+o(1)` for the growing polylog arities, and targetwise coverage is at least `1/2-o(1)` by the stated second-moment estimate (target 2T != O). Nonregular masks have probability at most `binom(m,2)w^2/p + m w/p=o(1)` uniformly in T. Reject such inputs deterministically, retaining constant coverage. No twist-error term is necessary for this norm predicate.

## 5. Public coefficient expansion and exact Fourier ledger

For fixed seeds, precompute the scalar coefficients of every public polynomial

    P_I(X) = c_lead (X-X_T) product_{j not in I}(X-a_j)
                     product_{i in I}(-t_i),    I subseteq [m].

For each coefficient choose a short full field representative in R. The polynomial product on the right-hand side of the norm identity has a lift obtained by summing these public coefficient polynomials times `product_{i in I} U_i`. In particular the high degree in all the **public** seed words does not increase the auxiliary Fourier precision; all products of those public words were already reduced in F_p. The highest coefficient may use the identical short leading-norm representative on both sides. This preprocessing costs `2^m poly(b,m)`, polynomial in b at the optimized arity, and does not require enumerating any private tuples.

For each formal-X coefficient of the norm residual the left-hand side in R is a sum of quadratic products of full private section coefficients, multiplied by short public weights from Section 2. Its raw coefficient width is

    (5/2) k + O(log b+log m).

The lifted right-hand side has coefficient width at most

    (1+1/theta+o(1)) k + O(m log b).

Indeed its public coefficient has k+O(log b) bits and all m private thin factors together have at most `(1/theta+o(1)) k+O(m)` bits of magnitude exponent, plus O(m log n) for convolutions and summing over subsets. Here `m=Theta(L^(1/4))` and `m log b=o(k)`. Let

    W = (max(5/2, 1+1/theta)+o(1)) k.

Multiplying the residual by E adds only O(log b) bits. The norm bound for `A_M^-1` and the construction of M imply that the quotient Q can be bounded so that both sides of `EF=M Q` and their difference have coefficient width W (with explicit lower-order safety margins). In particular `||M Q||_infinity <= poly(n)||EF||_infinity` for all accepted bounded Q. Choose the auxiliary prime power modulus larger than twice a safe coefficient bound for the difference. Its bit length is W. Use the inherited construction with auxiliary prime `q=-1 mod 8n`, real twist beta, norm-one primitive n-th root, and paired conjugate frequencies; multiplication in R is pointwise after the twist and invertible transform. The finite ring operations then verify exact integer equality of each residual, not a random fingerprint. The separate target condition is linear and costs one additional equation and quotient array.

The two-stage transform has n first-stage boundary words per complete array, counting two components for every conjugate pair. On the middle processors compute every formal-X coefficient of the norm and the public-expanded product in the frequency ring and verify the projected equation against the quotient. No output arrays are transmitted. Private section coefficient arrays number m, quotient arrays m+O(1), and the variable U_i arrays use only the theta-fraction active rows (the other rows have public zero input). Therefore the leading word-bit budget is

    a(theta) m b = (2+theta) max(5/2,1+1/theta) m b.

Auxiliary prefix checks and relays are lower order; public tables and spectra are polynomial-size. The graph is bipartite, with polynomially many vertices, maximum degree `O(m sqrt(n) W)` and parallel multiplicity `O(m W)` up to fixed factors and lower-order routing. These are square-root-in-b times polylog b and polylog b, respectively. Local predicates invert their row transforms to check digit bounds, canonical ring word encodings and any bit pins. The middle predicates compute only polynomially many field/ring operations at this arity, including the 2^m public expansion terms. An exact Boolean contraction over Z/2^(M_edges+1)Z can count all certificates without wraparound. Standard bitwise pinning recovers a certificate with polynomially many count queries; the soundness proof recovers the rational tuple. All stated graph conditions require the inherited detailed numerical transfer audit, but no new numerical contractor constant is introduced.

For theta<=2/3, the expression `(2+theta)(1+1/theta)=3+theta+2/theta` decreases; for theta>=2/3, the expression `(5/2)(2+theta)` increases. The optimum is theta=2/3 and

    a=20/3,   W=(5/2+o(1)) k.

The scalar ratio approximation and the choice of active owner fraction are both essential. Using an ordinary three-term cubic with full-size public coefficients restores the larger width; spreading thin variables across all owners pays 7.5mb rather than 20mb/3.

## 6. End-to-end numerical and order comparison

Let `K=2*pi*sqrt(14*log2(90))/(9 ln 2)` as in the prior manuscript. Under the inherited quantitative contraction estimate, the log relation-query cost is `(10K/3+o(1)) m b/sqrt(L)`. With the original sparse solver and one-pilot transfer, the log-time terms are `b/m + (10K/3) m b/sqrt(L)` and `2b/m`. The constant coverage and regularity check have only lower-order costs. Auxiliary lattice construction (including all short pairs) costs `2^O(b/L)`, and factor enumeration has log cost `b/m+O(log b+log m)`. Thus optimize `m=(1+o(1)) L^(1/4)/sqrt(10K/3)` to get the claimed time coefficient `2 sqrt(10K/3)=11.31486813704766` and the inherited space coefficient `sqrt(10K/3)=5.65743406852383`. This is an 8.7129% reduction relative to the previous 12.39481702771492 candidate; the scale is unchanged.

The big-O order can alternatively use the public Boolean contraction theorem in the October 25, 2026 preprint *The Exponential Time Hypothesis is False*, which gives `2^O(M/sqrt(log M))` for the relevant degree, description, entry and precision bounds. See `https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-reduction.tex` and the proof in `quantitative-contraction.tex`. This route uses the stated public theorem, not the separate numerical finite-packet assumption; it does not imply the numerical leading coefficient and is not an independently verified proof of that preprint. Ordinary polynomial-time linear algebra in the factor-base size suffices for the big-O order. The one-pilot size-bias proof in the earlier large-arity note applies: choose targets as independent uniform `[u]G+[v]Q`, excluding identity, run the same verified witness finder in training and testing, collect more rows than labeled factor-base columns, and compute a dependence without using u,v. Conditional on the targets and relation rows the v's remain uniform, allowing verification of the extracted log. The expected query overhead is bounded by the average cost divided by average success via the size-biased pilot. Neither route improves the logarithmic saving power beyond 1/4.

## 7. Pairwise tensorization of the contractor's crossing sum

This is a separate exact observation about the public contractor. At an interpolation node its mode expansion has the form

    F = sum_{sum_Z w_Z beta_Z < P} product_g Z_g(beta_incident(g)) mod 2^P.

Every mode has one source and one target group, and the child contraction Z_g depends only on its incident mode exponents. This locality follows from the public packet coefficient extraction identity. Group modes into blocks by their unordered endpoint pair (direction remains part of the label). Let the total number of nonempty blocks be e<=r(r-1)/2 for r groups. Enumerate nonnegative exact block budgets h_e with sum<P, using at most P^e budget vectors. In a block of weights w_1,...,w_v, define

    N(h) = #{beta>=0 : sum_j w_j beta_j=h}.

For each fixed budget vector encode the admissible block exponent list by a binary index of length `ceil(log2 N(h))`. Omit blocks of size 1 and skip empty blocks. A polynomial-size dynamic program over weights and budget computes counts and lexicographic unranking exactly (P and the number of modes are polynomial in the original size). Pad invalid binary codes with zero at either endpoint. Now the sum over all mode exponents for this budget vector is precisely a contraction on the r group vertices, with one bundle of binary edges per nontrivial block. A group tensor evaluates its ordinary child contraction with the decoded incident exponents. If bipartite form is required, subdivide each binary edge using an equality tensor. This costs a factor two in wires and preserves the count. All creation-time precision and ring arithmetic must still be paid in evaluating the group entries.

For each budget vector the product of its block cardinalities is bounded by the number of global feasible mode exponents. The public partition bound for M original wires, u retained valuation levels and its tapered precision profile gives `exp(O(M/sqrt(u)))` such exponent vectors. Hence before subdivision the compressed graph has at most

    S = O(M/sqrt(u)+r^2)

binary edges. This is an exact compression, not yet a fast algorithm: the tensor entries have the cost of ordinary child contractions and the r main vertices may be very high degree. The budget enumeration has log cost `O(r^2 log P)`, which is lower order when r=M^o(1). The code independently verifies all padded-index and exact-budget identities on 80 arbitrary finite mode networks.

### 7.1 Why the naive second application does not improve the power

The public packing proof uses `a=O(log M/log r)` label coordinates for r groups: its array size is `r^(2a)` and the squared array size times maximum degree must be bounded by M in the stated sufficient packing condition. Thus its proven number of useful levels is `u=O(log M/log r)` when r grows. Using the optimized level count in this bound gives a first-stage certified size bound of order `M sqrt(log r/log M)` (add the lower-order r^2). This is an upper-bound accounting, not an entropy lower bound on actual nonzero contributions.

The new compressed graph has only r main vertices; even after edge subdivision its maximum degree is at least `2S/r` up to a constant. In applying the same variance-based packing proof to this graph, the sufficient packing inequality for a constant number of second-stage groups allows at most `O(log(S/Delta))=O(log r)` retained coordinate levels. The corresponding optimistic second-stage estimate would have exponent of order

    [M sqrt(log r/log M)] / sqrt(log r) = M/sqrt(log M).

The two logarithmic gains cancel. This already ignores the extra inherited precision and charged child-entry routines, which are real further requirements. The formula is not a lower bound on a different packing proof, specialized mode graphs, multilinear coefficient-pairing algorithms, or ECDLP. It only explains why composing the *existing* estimates in the most straightforward way is not a proof of a stronger exponent. Large groups require fewer usable first-stage interpolation levels; small groups leave a dense high-degree mode network.

### 7.2 Multidimensional deformation diagnosis

A product of scalar Cauchy kernels over several independent coordinate chunks does permit grouping by all their active coordinates. Its degree in each scalar parameter is still proportional to the full original wire count times that chunk width. Cross interactions that differ in only one grouped coordinate remain: the tensor product retains the other coordinate factors, which can have full rank. One cannot charge every crossing term the product of the dyadic suppressions of all coordinates. A tensor-product Hermite readout also has the separate coordinate precision losses to pay. Neither the previous XOR variant nor this direct product construction establishes the missing improved mode spectrum. A different geometry that provides many independent grouping labels per unit of interpolation loss, or a method for contracting the high-degree compressed mode network without the generic packing loss, remains the concrete next target. This is a diagnostic, not a multivariate impossibility theorem.

## 8. Exact finite regressions

Run from this directory with Python and SymPy:

    python -B section_norm_audit.py
    python -B mode_budget_audit.py

The section test independently enumerates all normalized sections on five smooth toy curves over F_7,F_11,F_13,F_17,F_19 at arities 2,3,4 and compares the distinct reduced divisor support with direct elliptic group sums. It also reconstructs every target sign, tests thin digit injectivity on two large-field fixtures, verifies the short ratio and square-class normalization, constructs genuine norm equations from random rational point tuples, checks all projected quotient equalities, rejects the mutated norm systems, and checks separate target-sign lifts. The 96-bit binomial and 128-bit cyclotomic finite fixtures reuse the prior prime-degree n=5 principal-ideal helper; they are not numerical tests of the asymptotic balanced composite-degree layout or its wire coefficient. Exact SVP is the asymptotic preprocessing bound, not the LLL used in these finite fixtures. The mode audit uses arbitrary deterministic local integer tables with weighted exponent blocks and compares literal enumeration with binary-index contraction, including invalid padding labels. Both transcripts are retained.

**Status:** a smaller conditional leading coefficient and a new exact mode compression, with a specific accounting reason the obvious second use of the same contractor does not change the asymptotic power. The desired general-prime `2^O(b/(log b)^(1/3))` and square-root-in-b bounds remain unproved.
