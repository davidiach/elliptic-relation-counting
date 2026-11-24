# Average-case weak relation detection: recovery, robust slice rank, and character readouts

Research checkpoint, 24 November 2026. Continuation of `seed_local_slice_rank_note.md`, `contractor_degree_note.md`, `affine_digit_sampler_note.md`, and the uploaded v3 preprint. The constructive result is an average-target weak-score-to-verified-witness reduction. The negative results concern specified real-valued and finite-field readout models. No improved general-prime ECDLP exponent, independently verified contractor, or cryptographic-scale timing is claimed. The previous proposed coefficient 12.39481703 remains conditional and is not re-certified here.

## 1. Setup and exact laws

Let p>3, b=ceil(log2 p), G=E(F_p), N=|G|, and let T range over nonidentity elements of the supplied odd prime-order subgroup. For slot i choose a public nonempty proper support S_i of size w_i, and independently choose a field-affine seed (a_i,t_i) uniform in F_p x F_p^*. Keep all affine rational points with x-coordinate in a_i+t_i S_i. Write rho_i=w_i/p, nu_i=N rho_i, W=product_i w_i, w_max=max_i w_i, lambda=N^(m-1) product_i rho_i. H_T is the ordinary ordered signed relation count, and D_s=product_i |F_i|. All statements below concern these masks, not the original modulo-2^b affine hashes unless explicitly noted.

P is the independent fresh law. K_{i,x} is the law of a seed conditioned on including x: sample t uniformly nonzero and s uniformly in S_i, then take a=x-ts. R_i averages K_{i,x(P)} over a uniform nonidentity rational point P, and R=product_i R_i. Q_T samples an ordered uniform rational affine tuple of sum T and then applies these kernels independently. The exact Radon-Nikodym derivatives are

    dQ_T/dP = H_T/mu,     dR/dP = D_s/D_0,
    mu = (((N-1)^m-(-1)^m)/N) product_i rho_i,
    D_0 = (N-1)^m product_i rho_i.

Here mu is independent of the chosen nonidentity target. See the preceding affine digit sampler note for the exact pair-retention and factorial moment bounds. In particular E H_T^2=O(mu), uniformly, at sparse widths with lambda<=1, sum_i 1/nu_i=o(1), and 2^m/N=o(1). The conclusions below allow lambda=2^{-O(b/m)}, rather than requiring exactly 1/log b, provided those conditions hold. For balanced growing arity up to O(sqrt b), the earlier sampler and its coverage estimates apply. The factor base has at most 2 sum_i w_i labeled points.

For reference, the one-slot abscissa incidence covariance is exactly c_i (I-J/p) restricted to the liftable abscissas, where

    c_i = rho_i (1-rho_i) p/(p-1) <= rho_i.

This follows by counting the affine maps of each ordered pair of distinct abscissas. The covariance on group points, duplicating coordinates for opposite signs, has operator norm at most 2c_i. In particular Var(|F_i|)<=2 E|F_i|, so E_P[(D_s/D_0)^2]<=product_i(1+2/E|F_i|)=O(1).

## 2. Constructive theorem: only average-target score advantage is needed

Let f(T,s) in [-1,1] be a possibly randomized evaluator. Its expectation over private coins is denoted by the same letter. Put

    Delta_T = E_{Q_T} f(T,s) - E_R f(T,s),
    delta = E_T |Delta_T| > 0,
    Cbar = E_{T,P,coins} cost(f).

Assume a known lower bound on delta and Cbar, and 2^m/N=o(delta). Targets are independent and uniform in the nonidentity subgroup elements. The latter small-error assumption is automatic at the proposed subexponential inverse-advantage tunings. Then a verified witness attempt can be constructed with cost per success bounded by

    (1+sum_i w_i)(1+Cbar) poly(b,m,lambda^{-1},2^m,delta^{-1}).       (2.1)

The statement includes conditional score simulations, root-slot enumeration, sign reconstruction, and retries. Public support setup, where needed, is paid separately. It is a sufficient condition, not the construction of a score meeting it.

### 2.1. Average-time cutoff

If necessary, cap each score call at C Cbar/(mu delta^2), returning zero on timeout. Markov's inequality bounds the fresh timeout probability by mu delta^2/C. Cauchy-Schwarz, the two likelihood identities in Section 1 and E H^2=O(mu) bound its planted and independently marked probabilities by O(delta/sqrt C). This holds after averaging over targets and private coins. Taking a large absolute C preserves at least half the absolute average advantage, because | |Delta_T|-|Delta'_T| | is at most the sum of the two timeout probabilities. The capped score is uniformly evaluable even when called on the conditional seeds required below. Rename the preserved lower bound delta.

### 2.2. Proper marginal and mixed score

For any fixed T!=O, the distribution of the sum of j independent uniform affine points is [1+(-1)^j(N 1_{u=O}-1)/(N-1)^j]/N. Consequently, observing k<m marks from the planted tuple differs from independent marks by total variation O((N-1)^{-(m-k)}). To see the average rather than pointwise estimate, condition on the observed sum and use that its probability of equaling T is O(1/N). Passing marks through seed kernels does not increase total variation.

Set B_{i,x}=(K_{i,x}-R_i)/2. For 0<=k<=m let

    g_k(T,x_1,...,x_k;s_{k+1},...,s_m)
       = (tensor_{i<=k} B_{i,x_i}) f(T,s_{>k}).

Write Gamma_T=E_{planted marks} g_m. Expanding the product of differences and applying the proper-marginal bound gives, uniformly,

    |2^m Gamma_T - Delta_T| = O(2^m/N).

Thus E_T |Gamma_T| >= delta/(2*2^m) eventually. Neither the sign nor the size of Gamma_T must be known for an individual target.

### 2.3. Full-path survival, even for rare useful targets

Take h=delta/(16*2^m), decreasing the absolute constant if necessary to absorb the cutoff. Draw T and a planted tuple and root seeds. Independently for each slot choose a fair sign sigma_i. On plus use its root seed as the working seed; on minus draw an independent R_i seed. Let

    U = f(T,s_work) product_i sigma_i.

For fixed T and planted marks, after revealing the working suffix and its signs backwards, the conditional expectation of U when k prefix slots are hidden is (product_{i>k} sigma_i) g_k. In the analysis only, multiply U by sign(Gamma_T), taking either sign if Gamma_T=0. Its mean is E_T |Gamma_T|>=8h and its absolute value is at most one. Stop the backwards revelation the first time |g_k|<4h at a k in {m,...,1}. The contribution of all stopped paths to the mean of the signed U has absolute value at most 4h, by conditioning at the stopping time (a finite martingale argument). Hence the probability that the entire path clears every barrier is at least 4h.

Repeat the working-vector experiment with fresh roots under P and let J_good count the signed tuples of sum T whose prefixes clear the barriers. The planting likelihood gives E J_good>=4 mu h, averaged over T, roots and working coins. Since J_good<=H_T and E H_T^2=O(mu), Cauchy-Schwarz gives

    Pr(J_good>0) >= c mu h^2 =: zeta.

When the sampler additionally satisfies E H(H-1)=o(mu h), one can improve the success bound to Omega(mu h) directly by subtracting the factorial moment from E J_good. This optional stronger bound requires the score to retain its advantage at that thinner width; it is not a general way to transfer an existing score to a new width.

### 2.4. Charged search and finite cap

For any real one-slot function v, the covariance identity implies

    sum_x |K_{i,x}v-P_i v|^2 <= rho_i^{-1} E_{P_i}|v|^2.

Indeed the vector of differences is the adjoint of the centered incidence map divided by rho_i. Centering instead at R_i v minimizes the squared deviation with weights pi_x=|{P:x(P)=x}|/(N-1). The ratio between the positive weights is at most two. Including the factor 1/2 in B gives

    sum_x |B_{i,x}v|^2 <= (2 rho_i)^{-1} E_{P_i}|v|^2.

For odd prime-order full groups the weights are uniform and the stronger (4 rho_i)^(-1) holds. Tensorize this inequality to obtain, for every fixed suffix and target,

    sum_{x_1,...,x_k} |g_k|^2 <= product_{i<=k}(2 rho_i)^(-1).

The fresh root prefix is independent of the working suffix, and each prefix is present with probability product_{i<=k} rho_i. Thus the expected number of present prefixes with |g_k|>=h is at most 2^(-k) h^(-2), and the expected total is at most h^(-2).

Enumerate the root slot abscissas once. Extend retained prefixes one slot at a time and estimate each g_k by sampling the prefix signs and conditional K or matched R seeds independently, holding the original working suffix fixed. Each estimate is an average of bounded signed score evaluations. Keep prefixes with estimated |g_k|>=2h. Abort if the total retained count exceeds ceil(C/(zeta h^2)). With accurate estimates, Markov bounds this abort probability by at most zeta/10 for a sufficiently large C. There are at most (1+cap) sum_i w_i prefix tests. Hoeffding with O(h^-2 log((2+tests)/zeta)) samples per test controls every adaptive test's failure probability in total by zeta/10. Every full path above 4h is retained on that event. Enumerate up to 2^m sign combinations at each surviving full prefix, and verify the point sum exactly. No classifier output is trusted as a relation. Conditional kernels are sampled as in Section 1; uniform affine point sampling by rejection and square roots costs expected polynomial time. This proves (2.1), with a further global runtime cap if required, and no use of sign(Gamma_T) in the actual search.

### 2.5. Transfer to the outer reduction

The earlier one-pilot argument accepts a seed with probability proportional to its conditional verified-witness success on a fresh uniform target. Conditional expected work per subsequent row therefore averages to average attempt cost / average success. B+1 successful labeled rows for a factor base of size B suffice for a nonzero dependence, and uniformly random hidden target coefficients make the resulting logarithm determination verifiable. This is the `large_arity_research_note.md` reduction, and needs no targetwise guarantee. At m=Theta(L^(1/3)) for L=log2 b, (2.1) and the factor-base setup/sparse solve give the stronger O(b/L^(1/3)) log-time objective if log Cbar, log delta^{-1}, log lambda^{-1}=O(b/m). At m=Theta(sqrt b), the analogous tilde-O(sqrt b) bounds suffice, including the previous auxiliary-ring preprocessing if used. These are implications; the missing cheap score is not supplied by this theorem.

A finite control has nonzero advantage only at T=+G and -G, with opposite signs, so the signed average is exactly zero while the absolute average is positive. The exact existence score used for this control is deliberately expensive, not a candidate faster evaluator.

## 3. Robust slice-rank bound for bounded real scores

This is distinct from the previous exact characteristic-two slice-rank theorem. All functions and ranks in this section are over the real numbers with the product probability measure P. A bound on real multilinear degree is not interchangeable with one on degree modulo two. In particular, this result does not rule out low-degree *modular* contractions followed by parity extraction.

Let f(s,T) be a score with |f|<=1 and real slice rank q on the m seed axes for each fixed target. It may be any real linear combination of slice terms and may depend arbitrarily on the target. Then, for the same sparse parameters,

    |E_{Q_T} f-E_R f|
       <= C 2^(m/2) sqrt((q+m)/(lambda N)) + O((N-1)^(-m)).       (3.1)

The constant is absolute once mu/lambda is bounded away from zero and m/N is sufficiently small. Thus an advantage delta=2^{-o(b)}, with m=o(b), lambda=2^{-o(b)}, requires q>=2^{b-o(b)} (when the bound is nontrivial). The same bound applies to a randomized mixture of uniformly bounded scores of rank at most q, regardless of the number of choices of private coins.

### 3.1. Orthogonal centered decomposition

For each slot write I_i(P)=1_{P in F_i}, set z_i(P)=I_i(P)-rho_i 1_{P!=O}, and regard P->z_i(P) as a map from group coordinates to centered L^2 functions of the seed. Its operator norm is at most sqrt(2 c_i). The count difference Z=H_T-D_s/N is the nontrivial-character convolution part. Expand by the subset I of centered factors, using that convolution by 1_G-delta_O acts as minus the identity on the zero-sum subspace. Explicitly,

    Z_I = (-1)^(m-|I|) (product_{j notin I}rho_j)
          [(*_{i in I} z_i)(T) - (1/N) product_{i in I} sum_P z_i(P)].

The empty convolution is delta_O. The terms Z_I are orthogonal under P because every centered factor has mean zero and the seeds are independent.

For the top-order term, the core tensor on G^m is C_T=1_{sum P_i=T}-1/N. Flattening on any one mode has Gram matrix N^(m-2)(I-J/N), since the indicator rows have disjoint supports of size N^(m-2). After applying all m incidence maps, the squared operator norm of every top-order flattening is at most

    A=N^(m-2) product_i (2 c_i) <= 2^m lambda/N.                 (3.2)

For a nonempty I of size k, the centered count core on k modes has squared Frobenius norm <=N^(k-1). Hence

    ||Z_I||_2^2 <= lambda product_{j notin I}(rho_j/N)
                             product_{i in I}(2 c_i/rho_i).

The empty term satisfies this upper bound too by direct inspection at T!=O. Summing all proper subsets and using orthogonality, c_i<=rho_i, gives

    ||Z-Z_[m]||_2^2 <= lambda 2^m (product_i(1+rho_i/(2N))-1)
                    = O(lambda 2^m m/N).                       (3.3)

### 3.2. Slice-space projection

A real slice decomposition with q terms can be collected into subspaces V_i of dimension q_i, sum q_i<=q, on the individual axes: f belongs to the sum of the subspaces V_i tensored with all other axes. Their orthogonal projectors P_i commute since each acts on a different axis. The projector onto their sum is P_union=I-product_i(I-P_i), satisfying P_union<=sum_i P_i. By (3.2),

    ||P_union Z_[m]||_2^2 <= sum_i ||P_i Z_[m]||_2^2 <= q A.

Therefore |<f,Z_[m]>|<=sqrt(q A) ||f||_2. Combine (3.3), ||f||_2<=1 and the two likelihood identities. Replacing D_s/D_0 by D_s/(N mu) causes error at most |D_0/(N mu)-1|=O((N-1)^(-m)). Divide the count correlation by mu~lambda to obtain (3.1).

### 3.3. Consequence for real-degree local compilers and direct certificate counts

The combinatorial monomial-to-slice upper bound in Section 3 of `seed_local_slice_rank_note.md` holds over the reals as well. If a bounded real contraction score is a linear combination of K=2^{o(b)} seed-local networks, all with M D=o(m^2 b^2) (D is now real multilinear degree), its slice rank is 2^{o(b)}. The same conclusion for each random seed of an evaluator suffices. One cannot obtain a subexponential inverse planted advantage by that route. To be precise, if all MD <= m^2 B^2/16, each network has rank <=m 2^B by the prior entropy argument. Taking B=2 sqrt(MD)*O(1)/m yields the stated o(b) logarithm. This is not a theorem about arbitrary nonlinear combinations, nor about degree computed in a dyadic ring.

One may also thin the support to guarantee D_s/N <= 1/4 for every seed, e.g. balanced w with (2w)^m <= N/4. At growing m<=O(sqrt b) the extra inverse mean is only 2^{O(m)}, within the recovery budget when lambda and widths are chosen accordingly. Let Z_cert>=0 be any real contraction output supported on true positive seeds (not necessarily their exact count). Then <Z_cert,H-D_s/N> >= (3/4) E Z_cert. Repeating the projection proof with ||Z_cert||_2 gives

    (3/4) E Z_cert/||Z_cert||_2
        <= sqrt(2^m lambda/N)(sqrt(q)+O(sqrt m)).

If success mass is at least lambda 2^{-o(b)} and nonzero multiplicities lie between 1 and 2^{o(b)}, then E Z_cert/||Z_cert||_2 >= sqrt(success mass)/2^{o(b)}, forcing q>=2^{b-o(b)}. This extends the real-rank limitation to selective, one-sided, verified certificate counts with controlled redundancy. It does not constrain a modular readout with exponentially larger real degree or unbounded cancellations.

## 4. Polynomial-character and bit-readout obstruction

The earlier affine sampler note ruled out polynomial equality tests. The new result addresses low-degree **quadratic-character** and **additive-character** readouts, which can be much more nonlinear than equality tests. This applies uniformly to every eligible fixed target and any fixed supports, not just a random-support heuristic.

Let g(a_1,t_1,...,a_m,t_m) be a polynomial over F_p of total degree at most D, with D<p. Let chi be the quadratic character, chi(0)=0, and psi a nontrivial additive character of F_p. Under N~p, N>=2m and the stated model, both f=chi(g) and f=psi(g) satisfy

    |E_Q f-E_R f| <= C [ (D+1)/sqrt(p) + m(D+1)^2/p
                     + m(D+1)^2 w_max/(d_T W) + 1/W + 1/(N-1) ],   (4.1)

where d_T=m 2^(m-1) is the degree of the fixed-target abscissa hypersurface. The constant is absolute for sufficiently large p; cap the bound by 2 when it is vacuous. For multiplicative rational functions use numerator times denominator inside chi and add the pole-zero probability, which is of the same order for bounded total degrees. No analogous rational additive extension is required for the results below.

For the intended growing arities, W=2^{b-o(b)}, w_max=2^{o(b)}, D=2^{o(b)}, m=o(b). The right side is 2^{-b/2+o(b)}. This is too small for the subexponential inverse-advantage criterion in Section 2.

### 4.1. Geometry and exceptional polynomial zeros

Let W_T(x_1,...,x_m) be the absolutely irreducible polynomial defining the abscissa image of the sum-T fiber. The prior chord geometry proof gives total degree d_T and positive degree in every x_i: the sum fiber is E^(m-1), the last abscissa generically has 2^(m-1) distinct signed completions, and separately homogenizing shows the product of top powers has nonzero coefficient for T!=O and 2T!=O. For each offset vector s the transformed equation

    W_s(a,t)=W_T(a_1+t_1 s_1,...,a_m+t_m s_m)

is absolutely irreducible. Distinct offsets give nonassociate polynomials: specialize all t_i=0 to normalize the unit scalar, and vary an independent t_i where two offsets differ to contradict the positive x_i-degree. This is precisely the family in the earlier polynomial zero-test proof.

A nonzero polynomial B of degree d, evaluated after substituting a_i=x(P_i)-t_i s_i, is identically zero on the sum fiber times the t-space only if W_s divides it. Otherwise its conditional zero probability on Q at that fixed s is O(m(d+1)/p). For completeness, parametrize P_m=T-sum_{i<m}P_i on E^(m-1). Separately varying each P_i gives a rational function with pole divisor of degree at most 4d, at P_i=O and P_m=O, because each abscissa has pole order two. Each t variable has degree at most d. A product-grid zero bound follows by induction: choose a generic specialization of one variable where the remaining rational function is nonzero; fibers identically zero are zero at this specialization, and every other fiber has at most its pole-degree many zeros (or polynomial-degree many in a t variable). Include the at-most-two excluded points on each elliptic fiber and condition on all points being affine; this costs factor at most (1-m/N)^(-1). One may choose the generic specializations away from forbidden poles for the finitely many remaining grid choices. Thus the stated bound also holds when one t variable is absent. Under R, independent rational abscissas have atoms at most 2/(N-1), so weighted Schwartz-Zippel after the invertible change of variables gives O(d/p), uniformly in s. Averaging Q over offsets adds at most d/(d_T W) for the exceptional offset vectors. These claims reproduce the needed zero estimate of the prior affine digit sampler note rather than assuming independent marks under Q.

### 4.2. A good direction except on very few offsets

First consider a polynomial g that depends on every slot. For a fixed offset s_i change variables a_i=z_i-s_i t_i. A direction is bad if the resulting polynomial does not depend on t_i. Since D<p, this is equivalent to the formal derivative (partial_{t_i}-s_i partial_{a_i})g being zero. There is at most one bad s_i: if two were bad, both partials would be zero and, by the degree condition, g would not depend on slot i. Thus the probability all m directions are bad is at most 1/W. Choose the first good i by the offset tuple alone.

In the additive case take B to be the nonzero leading coefficient of the resulting polynomial in t_i over the polynomial ring in (z_i,a_{-i},t_{-i}). Its degree is <=D. In the multiplicative case reduce g first to c A^2 G with G a squarefree polynomial (factorization over F_p suffices), and apply the direction test to G. For a good direction G is squarefree over the fraction field of the other variables: by Gauss's lemma distinct nonconstant irreducibles remain nonassociate and separable, as their degrees are <p. Factors independent of t_i are units in that field. Take B to be leading coefficient times the resultant of G and its t_i derivative. It is a nonzero polynomial of degree O(D^2); its nonvanishing guarantees nonconstant squarefree specialization. The resultant degree follows immediately from its Sylvester determinant of size <=2D-1.

### 4.3. Count bad planted fibers using distinct hypersurfaces

Fix i and a good s_i. View B in the variables (z_i,a_{-i},t_{-i}). As s_{-i} varies, the equations W_T(z_i,a_j+s_j t_j for j!=i) are again irreducible and pairwise nonassociate by the same independent translation argument. Hence B can vanish identically on at most deg(B)/d_T of those hypersurfaces. The probability that the chosen i is good but has such an exceptional offset tuple, summed over possible i, is at most

    sum_i deg(B_i) w_i/(d_T W) = O(m D^2 w_max/(d_T W)).

This is an unconditional bound over offsets; it does not assume that the other offsets remain uniform conditional on the 'first good i' choice. For all remaining offset tuples, the zero estimate of 4.1 bounds the probability B=0 on the rational marked-point and other-multiplier draws by O(m(D+1)^2/p). Under R no hypersurface exception is needed.

Conditional on these data and B!=0, the remaining t_i is independent uniform in F_p^*. Its univariate polynomial is nonconstant of degree <p in the additive case, and nonconstant squarefree in the multiplicative case. The usual univariate Weil bounds give absolute character average <=((D-1)sqrt p+1)/(p-1), including removal of t_i=0. For reference see the standard polynomial additive bound (degree not divisible by p) and multiplicative bound (not a constant times a square) in the literature; e.g. https://math.ucr.edu/~mcc/paper/DubProc.pdf (multiplicative theorem 1.5), and https://riemann.tamnd.com/learn/analytic-number-theory/ (both stated). The argument here uses only these standard one-variable estimates, not a growing-dimension Lang-Weil bound. On bad fibers bound the score by one. The difference of Q and R averages is bounded by the sum of their absolute averages. This proves the full-slot case of (4.1).

If the relevant polynomial depends on a proper subset of slots, use the proper-marginal estimate in 2.2 instead. For chi(g), removing square factors changes the score only at zeros of A, which by 4.1 have probability O(D/(d_T W)+m(D+1)/p) under Q and O(D/p) under R. Constant or identically zero polynomials are covered as well. This completes the proof without assuming irreducibility of g itself.

### 4.4. Bounded transcript depth, not just individual tests

A depth-k ternary decision tree of Legendre-symbol tests of polynomials of degree <=D has output gap at most 6^k times the bound (4.1) at degree 2kD (provided 2kD<p). Indeed each branch indicator is one of (chi^2+chi)/2, (chi^2-chi)/2, or 1-chi^2. Expand each leaf product, use chi(g_1)chi(g_2)=chi(g_1 g_2) including zeros, and bound the coefficient l1 norm by 2^k per leaf. There are at most 3^k leaves. Branch polynomials can be selected adaptively from earlier symbols and public parameters; the proof is on each fixed transcript. A randomized mixture obeys the same bound. Thus k=o(b) such tests, with D=2^{o(b)}, cannot supply a subexponential inverse gap.

For additive readouts, every interval indicator on canonical field representatives has additive Fourier coefficient l1 norm O(log p). For a canonical bit j, first expand the square wave of period 2^(j+1) into its finite Fourier series, whose coefficient l1 norm is O(log p). Each restricted exponential exp(2 pi i h x/2^(j+1)) has an F_p Fourier expansion of l1 norm O(log p), by the geometric-series estimate min(1,1/(2p||h/2^(j+1)-k/p||)). Thus each bit indicator and its complement have l1 norm O((log p)^2). A depth-k binary tree of bit or interval tests on polynomials of degree <=D has gap at most (C (log p)^2)^k times (4.1) at degree D. To see this, expand the product of indicators on a fixed leaf in additive characters; the product of characters of polynomials is the character of their field-linear combination, still degree <=D. Sum over at most 2^k leaves. In particular k=o(b/log b) is insufficient for subexponential inverse gap.

These statements permit arbitrary field constants and target-dependent polynomial choices. They do NOT bound arbitrary computation using all the field values: a single full canonical residue exposes b bits and an arbitrary decoder may have exponentially large Fourier l1 norm. Nor does a low arithmetic circuit size imply small polynomial degree (repeated squaring is a counterexample). No general circuit, learning, or ECDLP lower bound is claimed.

### 4.5. Implication for the determinant-moment route

For formal full-sign moments Z_j, the sign sum cancels odd ordinate powers and y^2=x^3+ax+c gives a polynomial of total seed degree <=m floor(j(m+1)/2). Thus character tests and short bit/interval transcripts of subexponential-degree polynomial combinations of such moments fall under the theorem. The ordinary rational-only moments include a high-degree quadratic-residue selection and are not covered by this direct seed-degree argument. This explains why a successful nonlinear decoder must exploit something beyond the previous low-degree moment fits, even after switching from exact parity to a weak score.

## 5. Full-value moment experiment (outside the short-transcript theorem)

The held-out experiment tests the entire scalar Z2, not just its quadratic character. Choose the first prime-order smooth toy curve found by a deterministic small-coefficient search at p=101,257,1009,4093, with m,w as in the transcript. Supports are intervals {0,...,w-1}. Use the independently marked R law, not the fresh P law, as the control. At each field, generate independent training and held-out sets under Q and R. Train the full-field lookup by reporting positive whenever the empirical Q frequency at the measured value exceeds the R frequency. Test three scalars: formal full-sign Z2 (including non-rational signed roots), ordinary rational Z2, and rational Z2 divided by the product of rational slot sizes in F_p. Since Q and R are marked, the slots are nonempty and this normalization is well-defined for these toy parameters.

The second moment is evaluated by the mixed determinant: subtract the fixed target column from every variable column to remove the constant row. Form the Gram sum for each of the m slots on the remaining m coordinates. The coefficient of the product of all slot variables in the determinant of their weighted sum equals the sum of squared determinants, by Cauchy-Binet. Inclusion-exclusion on slot subsets extracts the coefficient since total degree is m. This avoids enumeration of all signed tuples. A separate literal check on 100 toy roots confirms the rational moment identity.

For each held-out difference of acceptance proportions, the displayed pointwise 95% half-width is sqrt(log(40)/n), where n is the number of held-out draws per law. This follows from Hoeffding for the sum of 2n independent bounded variables with ranges 1/n, conditional on the trained rule. These are pointwise, not simultaneous, intervals. Also record the unbiased collision estimator of p times squared L2 distance between the two scalar value distributions:

    p [sum c_Q(c_Q-1)/(n(n-1)) + sum c_R(c_R-1)/(n(n-1))
          - 2 sum c_Q c_R/n^2].

It may be negative from noise and is not a proved upper bound on total variation. The exact asymptotic full-value decoder question remains open. In these runs all twelve held-out confidence intervals include zero. The experiment supplies no positive advantage certificate for any of the tested scalar lookup rules; it cannot exclude weaker gaps or a better decoder. Storing/training a full p-entry table is also not an acceptable asymptotic algorithm at the desired exponent.

See `average_case_detector_audit_output.txt` for all numeric results. The experiment is intentionally restricted to small fields, and does not use a deployed cryptographic curve or perform end-to-end discrete-log solving.

## 6. A remaining constructive modular route: parity cancellation, not pointwise approximation

The robust real-rank bound does not apply to a low-degree contraction over F_2 with a parity readout. The exact F_2 slice-rank theorem from the previous checkpoint does not apply to an average-case detector. Here is a specific algebraic gate for that route, rather than treating an arbitrary pointwise probabilistic polynomial as an approximate parity oracle.

Suppose a candidate network has an efficiently represented Boolean predicate B_s(x) and a binary residual map r_s(x) in F_2^n, with zero residual corresponding to a verified certificate. Let c_s(v)=sum_{x:r_s(x)=v} B_s(x) mod 2 be its residual parity histogram. Its exact desired certificate parity is c_s(0), which need not equal ordinary count parity if certificate multiplicities vary. For a k x n binary matrix A use the relaxed parity score

    Z_A(s) = sum_x B_s(x) product_{j<=k}(1+(A r_s(x))_j) mod 2.

This is exact equality after k linear checks. Define the false residual moments

    M_J(s) = sum_{v!=0} c_s(v) product_{i in J} v_i mod 2,
            |J|<=k,

including the empty moment. Expanding the product of linear factors shows that Z_A(s)-c_s(0) is exactly a known linear combination of these moments, with coefficients depending on A. Repeated variables are idempotent on binary vectors. In particular:

* If the false residual parity histogram is a sum of indicators of affine flats of dimension >k that avoid zero, every such moment vanishes. Indeed restricting a degree<=k multilinear polynomial to an affine flat gives degree<=k on >k free binary variables, whose full sum is zero modulo two. Thus all linear-check sketches give exact certificate parity despite arbitrarily many false assignments.
* More generally, if these moments have predictable background values with high probability under both planted and independently marked laws, subtracting the corresponding background yields an approximate parity oracle. On a sparse genuine singleton, if certificate parity is odd with sufficient probability, this gives a planted score gap. All conditions, including the cost and degree of computing the relaxed score and the multiplicity of valid certificates, must be proved in an actual elliptic compiler.
* Merely approximating every local predicate pointwise does not establish this. Parity sums can turn many tiny individual errors into a large output error. For instance in a model of K independent uniform nonzero residuals, the conditional bias of the false parity given A is (1-2q_A)^K where q_A=(2^(n-rank A)-1)/(2^n-1). If k<n and K is much larger than 2^k, this generally gives negligible bias, not cancellation. This independent-residual example is a model calculation, not a claim about actual elliptic residuals.

The exact expansion and affine-flat cancellation lemma are verified on 200 finite fixtures. No suitable low-degree, high-advantage elliptic residual construction has been found. The next concrete mathematical gate is to establish such **predictable low-degree residual parity moments** for a paid compiler, or to produce a score from genuinely joint high-degree field computations that evades Section 4. This is a precise sufficient condition for further work, not another claimed speedup.

## 7. Reproduction and status

Requirements: Python 3, NumPy, SymPy. The finite verifier uses floating linear algebra for Hilbert norm/projection checks with explicit tolerances, exact modular arithmetic for the polynomial/character/moment checks, and fixed-seed independent Monte Carlo samples for the held-out probes. The analytic proofs, not the numerical SVD or Monte Carlo, establish the stated theorems.

    OPENBLAS_NUM_THREADS=1 python -B average_case_detector_audit.py
    python -B average_case_detector_audit.py 101 3 3 30000
    python -B average_case_detector_audit.py 257 4 3 30000
    python -B average_case_detector_audit.py 1009 4 4 50000
    python -B average_case_detector_audit.py 4093 5 4 60000

The four spectral fixtures check 2,380 covariance entries, 11 single-mode flattenings, 240 slice-space projection comparisons and 320 centered-expansion comparisons. The character fixture checks one exceptional factor among five offsets, 101 complete quadratic sums, 202 symbolic directional tests, and numerical bit Fourier norms. The average-gap control and literal moment comparison and residual-parity tests are included in the transcript. These finite checks are not exhaustive validation of the complete recovery or contraction algorithm.

Primary project sources: uploaded v3 preprint, Sections 3–5 and 8; project weak-detection derivation at https://github.com/davidiach/elliptic-relation-counting/blob/8171080efb0d9024f49eb47e1e311914f0f7a4c8/derivations/strided_slots/weak_detection.md ; fixed-target chord geometry at the same revision; previous checkpoint notes named above. All theorem transfers retain the relevant hypotheses. The 1/3 logarithmic saving and square-root-in-bitlength goals remain open. The results are AI-assisted research for independent scrutiny, not a novelty/priority or peer-review claim.
