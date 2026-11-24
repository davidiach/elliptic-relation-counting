# Field-affine thin-digit samplers and a planted polynomial-readout obstruction

24 November 2026. Independent research checkpoint following `large_arity_research_note.md` and the uploaded v3 paper, *Counting relations in independently hashed elliptic-curve slots*. The new results are targetwise coverage, controlled thinning, and an algebraic obstruction on seed-polynomial zero tests. They do not establish an improved ECDLP exponent or an end-to-end implementation. The proposed previous numerical coefficient 12.39481703 remains conditional on the previous compiler and finite-packet assumptions.

## 1. Thin digit blocks with exact support-size control

Assume the previous auxiliary polynomial ring R=Z[X]/(X^n-lambda), lambda in {1,2}, evaluation phi:R -> F_p, and short projector E and diagonally dominant polynomial M. Write k for a common coefficient bit width of a signed dyadic box whose image contains F_p, A=n k=b+o(b), with k/log b -> infinity. The balanced Fourier owner partition has V rows, each owning n/V coefficient positions up to rounding. All estimates in this section are conditional only on the previous representation construction and graph-transfer assumptions, not on a new heuristic. For the leading constant use n=Theta(b/(log b)^2) and balanced transform dimensions. For the qualitative large-arity version use the preceding prime-degree/padded Bluestein construction.

Let r grow and take m=r^2, which can approximate any growing asymptotic arity up to a relative o(1). Partition the V owners into r nearly equal groups, and partition each coefficient's k bits into r nearly equal consecutive slices. The Cartesian product gives m disjoint blocks partitioning all A bits. Within each block the bit slice begins at the same bit position lo. Its contribution to phi is 2^lo times the image of a small nonnegative coefficient polynomial U supported on that owner group. Every coefficient of U is less than 2^ceil(k/r). The signed offset of the whole box is public and has no effect on surjectivity of the sum of block images.

Let S_j be the image of a block after removal of the common factor 2^lo. Surjectivity gives product_j |S_j| >= p. Hence a largest block has |S| >= p^(1/m). Find it by enumerating the m blocks, deduplicating their images. For balanced n, V/r and k/r tending to infinity, every block has at most (1+o(1))b/m bits, including rounding. At m=Theta(sqrt(b)) in the prime-degree version, the number of owners in a group is O(1), k/r=O(b^(1/4)), and block length is O(sqrt(b)); the qualitative budget suffices there.

For any desired integer w<=|S|, scan the selected block's words in owner-major lexicographic order, maintaining their image set. Stop at the first prefix with w distinct images. Every increment adds at most one image, so every integer cardinality is available. The prefix threshold can be checked by one initialized comparison state passed between consecutive active owner rows. It uses O(V/r+1) constant-size messages, not a new full transform array. Enumeration and prefix construction cost 2^(O(b/m)) poly(b,m) at the relevant tunings; at the leading-constant tuning the log cost is (1+o(1))b/m. A graph may retain repeated representations of a point, but a deduplicated factor base and verified witness interpretation do not require globally canonical representatives. Exact counts of distinct point tuples must not be identified with counts of arbitrary representation assignments.

## 2. Exact random affine masks

For each slot independently draw a_i uniformly in F_p and t_i uniformly in F_p^*. Accept every rational affine curve point whose abscissa belongs to a_i+t_i S, where S now denotes the chosen support of cardinality w. Write rho=w/p. For every x, and any x!=y,

    P(x retained)=rho,
    P(x,y retained)=w(w-1)/(p(p-1)) <= rho^2.

Proof: for each s there are p-1 affine maps taking s to x. For any ordered distinct pair (s,u), there is exactly one field-affine map taking it to (x,y). Count these maps out of p(p-1). This argument requires no additive-combinatoric property of S and holds for arbitrary finite support.

To implement the mask, precompute short public representatives A_i,T_i of a_i,t_i. The variable abscissa polynomial is A_i+T_i U_i in R. Absorb the omitted factor 2^lo into the uniformly random multiplier when necessary. The middle Fourier processors compute its transform from the public spectra and the variable U_i spectrum. Only the O(V/r+1) active first-transform owners transmit U_i words. No full variable abscissa transform or quotient is required for this membership check. In coefficient space, multiplication by T_i raises the size bound to at most 2n ||T_i||_infty 2^ceil(k/r), so the bit width is k+o(k) when k >> log b and r->infinity. All short public field constants have k+O(1)-bit coefficients. Every term in the chord quartic, counting its public curve coefficients as well as variable operands, has coefficient size 2^((4+o(1))k) up to polynomial(n) factors. The preceding short projector adds O(log b) bits and diagonal dominance bounds the integral quotient. Thus the quartic exact equality still needs W=(4+o(1))k bits per Fourier word at the balanced asymptotic tuning. The additional variable leaf transforms cost O(m b/r)+o(mb)=o(mb), and the owner comparison messages are also lower order. The 8mb leading word budget of the prior abscissa-only compiler is unchanged. This is a graph-ledger result, not a verification of the downstream finite-packet realization.

## 3. Uniform targetwise second moment and sparse isolation

Let G=E(F_p), N=|G|, kappa=|G[2]|<=4. Let H count ordered rational affine tuples in the independently sampled slots summing to fixed T with 2T!=O. Define nu=N rho and lambda0=N^(m-1) rho^m. More general unequal supports give the same formulas with nu_i=N rho_i.

The expected count is exactly

    mu = ((N-1)^m-(-1)^m) rho^m/N.

The count of affine tuples is obtained by convolving the indicator 1_G-delta_O or by Fourier diagonalization on the finite abelian group. Its relative difference from lambda0 is O(m/N) for m=o(N).

For two distinct tuples, let I be the positions whose abscissas coincide. If I is proper, choose the first tuple in at most N^(m-1) ways. Choose one of at most two signs at each position in I and all but one unmatched point arbitrarily, the last being forced by the target. The number is at most 2^|I| N^(2m-|I|-2). Pairwise mask probabilities give the contribution bounded by lambda0^2 product_{i in I}(2/nu_i). If every abscissa coincides, a nonempty proper subset of sign flips must sum to G[2]. The all-sign flip is impossible because 2T!=O. For each such subset there are at most kappa N^(m-2) tuples, forcing one point on each side. The retention probability is product_i rho_i. Consequently

    E H(H-1) <= lambda0^2 product_i (1+2/nu_i)
                 + lambda0 kappa(2^m-2)/N.

This is the paper's isolation bound, now recovered with the field-affine thin-digit law. No average over targets is taken. If X additionally requires nonidentity prefix sums, at most O(m N^(m-2)) tuples are removed, by forcing the last point of a prefix and its complementary suffix. Hence E X >=(1-O(m/N))lambda0. Genuine doubling is allowed in the relaxed abscissa chain; the previous mixed-twist error bound and verification apply. Paley-Zygmund gives

    P(X>0) >= (1-O(m/N))^2 lambda0 /
      (1 + lambda0 product_i(1+2/nu_i)+kappa(2^m-2)/N).

For balanced widths w=ceil(p^(1/m)), m=o(b), nu>>m, and m/ sqrt(p)->0, lambda0=1+o(1) if m/w->0. Thus P(X>0)>=1/2-o(1) uniformly at every eligible T. The previous average-shift coverage could cost 2^O(m+log R_S) retries; that term is no longer necessary. For sparse lambda*=1/log_2 b, choose w=floor((p lambda*)^(1/m)) using the prefix construction. At m=Theta(sqrt(b)), or the earlier polylog arities, lambda0=(1+o(1))lambda*, 2^m/N=o(lambda*), and nu>>m. The same bound yields P(H=1)=lambda0+O(lambda0^2), P(H>=2)=O(lambda0^2), uniformly in T. X has the corresponding estimates. The mixed-twist probability is still at most 4(2^m-2) max(N_+,N_-)^(m-2) rho^m, exponentially smaller than the genuine positive probability. A returned witness is always verified in E(F_p).

## 4. A fixed-target obstruction to polynomial zero-test decoders

This is a new algebraic limitation for the new seed law, not a lower bound on arbitrary algorithms, Boolean circuits, character tests, or high-degree finite-field calculations. Fix T with 2T!=O, m>=2, and nonempty supports S_1,...,S_m. Write W_prod=product_i |S_i|. Let P be the fresh independent affine-seed law. Define Q by sampling a uniformly random rational affine tuple (P_1,...,P_m) of sum T, offsets s_i uniform in S_i, and t_i uniform in F_p^*, then setting a_i=x(P_i)-t_i s_i. Every fresh seed containing an ordered tuple has exactly one compatible offset for that tuple at each slot, hence dQ/dP=H/mu, exactly as in the uploaded paper.

Let W_T(x_1,...,x_m) be an absolutely irreducible equation for the abscissa image of the fixed-target sum fiber. The existing project's chord-geometry supplement proves degree 2^(m-1) separately in each x_i and total degree d_T=m 2^(m-1) for 2T!=O. Only irreducibility and positive degree in each x_i are needed for the following argument. For each offset tuple s define

    W_s(a,t)=W_T(a_1+t_1 s_1,...,a_m+t_m s_m).

The substitution is part of an invertible linear change of the 2m affine variables, so W_s is absolutely irreducible and has total degree d_T. Distinct offset tuples give nonassociate polynomials: at t=0 all specialize to W_T(a), so an associate scalar must be one. If s and s' differ in position i, compose with the inverse change for s; equality would imply invariance of W_T under translation of x_i by an independent indeterminate times a nonzero scalar. A nonconstant polynomial in x_i cannot satisfy that identity, including in positive characteristic. Thus all W_prod irreducible factors are distinct.

Let g(a,t) be a nonzero polynomial of total degree D. It can be divisible by at most D/d_T of the W_s. For every other offset tuple, the substitution a_i=x(P_i)-t_i s_i into g is a nonzero rational function on E^(m-1) x A^m, where P_m=T-sum_{i<m}P_i. Indeed the abscissa image is dense in W_T by the definition of that hypersurface, and the t coordinates are free. A polynomial vanishing identically on that product would be divisible by W_s, contradicting the choice of offset tuple.

For such a nonzero rational function, separately varying one P_i gives a rational function on E with pole divisor of degree at most 2D: x(P_i) and x(P_m) have a double pole each (possibly the same point), and the sum of exponents of any monomial is at most D. The same total bound holds for a sum, using the common effective divisor with multiplicities 2 max exponent at each point; a conservative uniform bound 4D is sufficient and is used below, because two separate monomials can attain maximal poles at different points. There are at most two forbidden points (P_i=O or P_m=O). Separately varying one t_i gives a polynomial of degree at most D. The following product-grid zero bound follows by induction: choose a generic algebraic specialization of one variable for which the remaining rational function is nonzero. Any fiber identically zero or wholly forbidden is bad also at that specialization. Apply induction to these fibers. Every other fiber has at most 4D zeros and two forbidden points on E, or D zeros on F_p^*. The specialization can be taken away from identities and from t=0; after specialization the same separate degree and forbidden-pole bounds hold. If all curve variables have been specialized, choose the generic values with P_m!=O. Therefore under independent uniform P_1,...,P_{m-1} in G and t_i in F_p^*, the probability of zero or forbidden input is at most

    (m-1)(4D+2)/N + m D/(p-1).

Conditioning all points to be affine (including the forced P_m) costs a factor at most 1/(1-m/N). Thus

    Q[g=0] <= D/(d_T W_prod)
          + ((m-1)(4D+2)/N + m D/(p-1))/(1-m/N).

The fresh law also satisfies P[g=0]<=D/(p-1) by Schwartz-Zippel on F_p^m x (F_p^*)^m. The bounds may of course be capped at one. In particular when m=O(sqrt(b)), W_prod=Theta(p/polylog b) or larger, and D=2^o(b), both probabilities are 2^(-b+o(b)). The same is true for a union of 2^o(b) nontrivial tests of such degrees. An adaptive equality-only decision tree with 2^o(b) tests along its default path has the same conclusion: follow the formal default path (zero for identically zero polynomials, nonzero otherwise), and union-bound the exceptional outcomes along that path. Hence the output is asymptotically constant under both laws. Tests selected from actual numerical field values, or arbitrary tests on the field values, are outside this statement.

For sparse isolation P(H>0)=lambda0+O(lambda0^2). For an event A with Q(A)=o(1), P(H>0,A)<=mu Q(A)=o(lambda0). A low-degree equality detector reporting positive only on exceptional zero-tests consequently misses essentially all positive mass. If it reports positive on the default path instead, the fresh false-positive rate tends to one. The proof applies to arbitrary polynomial combinations of the seeds, not merely a fixed finite monomial fit. It complements the original paper's separate bounds concerning target pole order, rank and moment spectrum.

## 5. Application to determinant moments and explicit limitations

The following observation applies to the twist-inclusive moments with both signed roots in F_{p^2}, counting a zero ordinate twice, or to formal signed abscissa tuples. This convention is not the ordinary rational-point moment of the v3 paper: the latter includes the quadratic-residue selection, which is high degree in the field-affine seeds. The mixed-twist probability bound concerns the existence of genuine nondegenerate elliptic sums, not all determinant zeros with repeated points. Do not conflate the two.

For the formal full-sign moments, in each slot an integrand Delta^j is a polynomial in x_i,y_i of pole weight at most j d, d=m+1; summing both signs cancels odd y_i powers, and replace y_i^(2u) by (x_i^3+a x_i+c)^u. The resulting polynomial has x_i-degree at most floor(j d/2). With fixed support S_i and x_i=a_i+t_i s_i, the moment Z_j is therefore a polynomial in all 2m seed variables of total degree <= m floor(j d/2). For a polynomial equality readout of total degree h in moments through order k0, the substituted seed polynomial has degree <=h m k0 d/2. If that is subexponential in b, the preceding obstruction applies unless the substituted polynomial is identically zero (then the equality supplies no information). The theorem covers arbitrary exact low-degree polynomial combinations, not only low-degree linear fits, for these formal moments. Character tests, rational-point-only moments, nonlinear value lookup, and high-degree powering (notably Fermat's indicator of degree p-1) are not ruled out.

## 6. Remaining exponent and candidate cost accounting

At m=Theta((log b)^(1/4)) the preceding 8mb compiler still gives the proposed conditional coefficient 12.39481703, with all the earlier finite-packet and graph/reduction assumptions. The new sampler changes coverage and corrects the affine-mask overhead, not the leading exponent. At m=Theta(sqrt(b)), all factor-base, affine sampling, coverage and polynomial-ring preprocessing costs are within 2^O(sqrt(b) polylog b). The missing relation evaluator still has O(mb) transmitted bits. Applying the public generic contractor gives log cost O(mb/sqrt(log(mb))), which does not meet the desired goal. An improvement to O(M/(log M)^alpha) with alpha>1/2 on the relevant network class would improve the outer exponent to O(b/(log b)^(alpha/2)) by balancing b/m and m b/(log b)^alpha. No such stronger contraction bound is proved here.

The public generic contractor's weighted-mode sum can be factored by grouping modes according to unordered source/target group, with exact budget allocations. This is an alternative computation of its already-known group join, not by itself an exponent improvement. If r groups grow, generic balanced packing requires label-pair arrays r^(2a), giving u=O(log M/log r). A second generic contraction on grouped mode labels would, at best under a high-degree extension of the contractor, have a logarithmic scale controlled by log r; the two savings cancel in this accounting. This is a heuristic cost audit, not a general obstruction theorem. The supplied compressed-relation manuscript already uses a nontrivial nine-group matrix join, so a simpler three-way join must not be claimed as an improvement over its coefficient.

## 7. Finite checks

Run `python affine_digit_sampler_audit.py` (NumPy). The verifier checks exact pair retention over every affine seed on fields 13,17,37; four toy rectangular digit-block decompositions; exact first and second moments on 98,320 ordered seed pairs over two curves, at 36 eligible targets; and full planted monomial evaluation ranks in the following finite settings. These are component tests. The general proofs above do not follow just from those ranks, and no cryptographic-scale relation oracle was run.

- F_13, m=2, degree <=8: 495/495.
- F_17, m=2, degree <=8: 495/495.
- F_13, m=3, degree <=4: 210/210.

Sources: uploaded v3 paper, Sections 2–5 and 8; previous `large_arity_research_note.md` and `abscissa_chain_note.md`; project supplement `derivations/group_equations/chord_geometry.md` at github.com/davidiach/elliptic-relation-counting; public OpenAI quantitative Boolean tensor-contraction theorem, October 25, 2026. This note is an AI-assisted research checkpoint; novelty and external review are not established.
