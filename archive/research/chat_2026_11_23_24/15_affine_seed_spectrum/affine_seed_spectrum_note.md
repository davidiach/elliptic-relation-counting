# Exact spectral normal form and higher-order incidence for field-affine elliptic samplers

Research checkpoint, 24 November 2026. Continues `affine_digit_sampler_note.md`, `numerical_sketch_note.md`, and the v3 paper *Counting relations in independently hashed elliptic-curve slots*. This is a new analysis of the actual **field-affine** sampler, not independent uniform lists or the paper's original modulo-2^b sampler. The strict general-prime ECDLP exponent is **not improved**. The earlier numerical coefficient 11.31486814 is conditional and is not re-audited here. No claim of literature novelty or peer review is made.

## 1. Model

Let E/F_p be a smooth elliptic curve in characteristic p>3, with G=E(F_p) of odd prime order r. Write b=ceil(log_2 p). Let S_i be a public subset of F_p of size 1<=w_i<p, W=product_i w_i, and m>=2. A fresh seed consists of independent a_i uniform in F_p and t_i uniform in F_p^*, retaining the rational affine points whose abscissas belong to a_i+t_i S_i. This is the field-affine sampler of the preceding checkpoint; the theorem does not require thin-digit support or any special arithmetic structure of S_i.

Use Q_T for the planted seed law: choose an ordered tuple of uniformly random affine points conditional on sum T, choose independent offsets s_i uniform in S_i and multipliers t_i uniform nonzero, then set a_i=x(P_i)-t_i s_i. Use R for the matched law in which the marks are independent uniform affine points and seeds are conditioned on the respective marks in the same way. Q_T and R are the relevant comparison, not Q_T against fresh uniform seeds (the latter also measures slot-size bias).

For several exact Fourier calculations first allow the identity as a mark, assigning it a public sentinel x(O) unused by the affine abscissas. Such a sentinel exists for all sufficiently large p by Hasse, and is checked in the finite fixtures. Write Q_T^+, R^+ for this extension; the superscript is not a twist designation. Under Q_T^+, every proper subset of the marked coordinates is independently uniform on G, because the remaining coordinate can complete any partial sum. Under either extended law, the probability that some mark is O is at most m/r. Conditioning on no identity recovers Q_T and R exactly. Thus for every score of modulus <=1, the difference of the two planted-minus-matched gaps (actual and extended) is at most 4m/r. This convention is only a mathematical device, not a change to the actual relation solver. The support-equivalence theorem below holds without it.

All normalized L2 norms and orthogonal projections on seeds refer to the *fresh uniform* seed measure. A group Fourier transform uses the normalized convention fhat(chi)=E_{P in G} f(P) overline(chi(P)); for additive phase coefficients we use the positive sign and state signs explicitly.

## 2. Exact incidence isometry; no support-shape improvement at second order

For one slot let L_x(a,t)=(p/w) 1_{x in a+tS}. It is the density of the seed conditioned on retaining x, relative to the fresh law. The sharply two-transitive affine action gives

    E L_x = 1,
    E L_x L_y = p/w                    if x=y,
              = p(w-1)/(w(p-1))        if x!=y.

Consequently, writing c=p(w-1)/(w(p-1)) and d=p(p-w)/(w(p-1)), the Gram matrix of the p functions is exactly c J+d I. In particular on zero-sum coefficients, the map from coefficient vectors v to sum_x v_x L_x is a scaled isometry of squared norm d. The Gram matrix is positive definite when w<p.

**Support-equivalence theorem.** For any two supports S,S' of the same width, there exists a real orthogonal transformation of the fresh seed function space mapping L_x^S to L_x^{S'} simultaneously for all x. It fixes the constant function. Indeed, their Gram matrices coincide, so the prescribed linear map is an isometry on the incidence span; extend it by orthonormal bases of the complements. The constant belongs to that span since sum_x L_x=p at every seed. Taking a tensor product of these transformations across slots maps both likelihood densities Q_T and R to the corresponding densities for the replacement supports, and hence maps their difference. The proof uses only their expansions as sums of products of incidence functions and applies to the actual affine-only marking law, any target admitting an affine tuple, and arbitrary choices of group-point weights independent of the seeds. It also applies to the artificial identity extension.

Thus **all singular values across all cuts, all L2 approximation costs defined by cut rank, and the real tensor's orthogonally invariant properties are exactly independent of support shape** at fixed widths. This does not assert equivalence of total variation, computational cost, positivity, nonlinear functions, or Boolean/modular tensors. The orthogonal transformations may themselves be extremely expensive, and are not proposed as part of an algorithm.

## 3. Exact translation and multiplier Fourier spectrum

Use the full-group marking extension. For additive character e_p(z)=exp(2 pi i z/p), define

    d_T(xi) = E_{sum P_i=T} e_p(sum_i xi_i x(P_i))
                 - E_{independent P_i} e_p(sum_i xi_i x(P_i)).

This is zero whenever at least one xi_i=0, by the proper marginal identity. Let

    B_i(v) = (1/w_i) sum_{s in S_i} e_p(v s).

Conditional on any fixed vector of nonzero multipliers t, the positive-sign Fourier coefficient in the translations of the density difference D_{T,t} is

    Dhat_{T,t}(xi) = d_T(xi) product_i B_i(-xi_i t_i).       (1)

The coefficient identity follows simply by substituting a_i=x(P_i)-t_i s_i and averaging the independent offsets. The multiplier marginal is identical and uniform under Q_T^+ and R^+; D_{T,t} is the difference of the conditional translation densities relative to uniform translations.

There is also an exact support statement for a convenient analytic extension of the multiplier law to all of F_p. At t_i=0 the conditional seed construction still draws the uniform offset, so incidence multiplicities must be counted (the actual algorithm never uses these multipliers). The positive-sign Fourier coefficient in (a,t) with frequencies (xi,eta) is

    Dhat_T(xi,eta) = W^(-1) d_T(xi)
         if every xi_i!=0 and eta_i/xi_i in S_i;
       = 0 otherwise.                                    (2)

This follows from the additive character identity E_t e_p(t(eta-xi s))=1_{eta=xi s}. The Fourier support is therefore contained in a product of explicit unions of lines. The extension differs from the nonzero-multiplier law by at most m/p in TV under either marked distribution; there is no claim that sparsity of a score persists when it is extended in this way. Equation (1), used for the main results, works directly with the nonzero multipliers.

### 3.1 A uniform hybrid character bound

For xi!=0 put f_xi(P)=e_p(xi x(P)) on G, with the specified sentinel value at O. Kohel and Shparlinski, *On Exponential Sums and Group Generators for Elliptic Curves over Finite Fields*, Theorem 1 and its single-pole specialization, https://www.i2m.univ-amu.fr/perso/david.kohel/pub/character.pdf , establish the mixed additive/group-character estimate. Applying it to the degree-two abscissa with its sole pole at O yields

    |fhat_xi(chi)| <= alpha := min(1,(3 sqrt(p)+1)/r)       (3)

uniformly over every group character chi; the +1 restores the value at O. This is a classical cited theorem, not a new claimed estimate. Without that theorem, a weaker alpha=O(p^-1/4) follows by squaring the Fourier coefficient, averaging the correlations e_p(xi(x(P+D)-x(P))), and applying the Bombieri curve sum bound to the two order-two poles for D!=O. The stronger (3) is used below.

Character orthogonality on G gives d_T as a sum over the nontrivial characters of products of fhat_xi_i, times a unit character at T (depending on Fourier sign convention). Parseval gives sum_chi |fhat_xi(chi)|^2=1. Hence for all xi with nonzero coordinates,

    |d_T(xi)| <= alpha^(m-2),                           (4)
    E_{T uniform in G} |d_T(xi)|^2 <= alpha^(2m-2).       (5)

For (4), bound m-2 factors by alpha and the sum of the product of the other two by Cauchy-Schwarz. For (5), use Parseval in the target and bound m-1 factors by alpha. These formulas require m>=2 and remain valid when the right sides are trivial for small fields.

### 3.2 Slope-aware sparse Fourier readouts

Suppose a score f(T,t,a) of modulus <=1 is, for every T and t, a trigonometric polynomial in the translations with at most K additive Fourier frequencies. It may choose the frequencies and coefficients arbitrarily using all multipliers, the target, and public or private randomness. From (1), (4), and Parseval, the absolute fixed-target advantage on the actual affine-only law is at most

    sqrt(K) alpha^(m-2) + 4m/r.                          (6)

Indeed the L1 norm of the score's Fourier coefficients is <=sqrt(K) and |B_i|<=1. Thus an advantage 2^-o(b) at growing m=O(sqrt(b)) requires log_2 K >= (m-2) log_2 p - O(m)-o(b). This is a restriction on the explicit Fourier support of the score, NOT on the number of arithmetic gates used to compute it. A nonlinear operation can generate huge Fourier support cheaply. The theorem allows arbitrary dependence on all the multipliers, even if that dependence is not separated across slots.

## 4. Stronger results for joint linear compression

### 4.1 Arbitrary target-independent, multiplier-aware projections

Fix nonzero t. Let a matrix L(t) over F_p of rank q<=m-1 be chosen arbitrarily using the multipliers and public coins, but independently of the target and translations. The decoder receives T,t,L(t)a, and may process these values arbitrarily nonlinearly. The rowspace V of L has an injective coordinate projection onto some q pivot coordinates J. Projection of a function of translations onto functions of L a retains exactly the additive Fourier frequencies in V. Moreover

    sum_{xi in V} product_i |B_i(-xi_i t_i)|^2
        <= product_{j in J} sum_{xi_j in F_p}|B_j(-xi_j t_j)|^2
        = product_{j in J} (p/w_j).                     (7)

The inequality drops the other factors (which have magnitude <=1) and uses the injective pivot map. The equality is Parseval on F_p; t_j!=0 is essential. Combining (5) and (7) bounds the expected squared L2 norm of the projected density difference, averaged over T, by alpha^(2m-2) product_{j in J}(p/w_j). This is a bound on the density, so arbitrary bounded nonlinear decoding is allowed. Jensen gives

    E_T | E_{Q_T} f - E_R f |
        <= alpha^(m-1) (p/w_min)^(q/2) + 4m/r,           (8)

where w_min=min w_i. For the mean over nonzero T, multiply the first term by sqrt(r/(r-1)); the same conditioning error works for each T. If w_min=p^(1/m) 2^(-o(b/m)) and m grows with m=O(sqrt(b)), even q=m-1 gives an upper bound 2^(-b/2+o(b)). Thus losing a single translation field word in this model kills the advantage needed by the weak-recovery theorem. The matrix may be different for every multiplier vector and is not required to have any local structure. Target-dependent choice of its rowspace is not covered by the use of target Parseval in (5).

For fixed-target matrices, using (4) instead of (5) gives alpha^(m-2)(p/w_min)^(q/2)+4m/r; it is already 2^(-b/2+o(b)) at q<=m-2. In particular this theorem is not merely the earlier independent-list summary obstruction: it permits a single joint linear combination of all translations, and uses the short actual seed description.

### 4.2 Random projection: exact energy loss, including target-dependent scores

For a uniformly random q-dimensional rowspace V chosen independently of the translations (and uniformly conditional on T,t), every nonzero frequency belongs to V with probability (p^q-1)/(p^m-1). Equation (1) has zero constant coefficient. It follows exactly that

    E_V || E[D_{T,t} | L a] ||_2^2
       = ((p^q-1)/(p^m-1)) ||D_{T,t}||_2^2.             (9)

Averaging over t and using the exact norm in Section 5 bounds any bounded score of the projected translations, for this random compression, by sqrt((p^q-1)/(p^m-1)) ||D_T||_2 +4m/r. This permits a target-aware readout and is pointwise in nonzero T. At the sparse tuning lambda=Theta(1/log b), dropping one field word gives O((p lambda)^(-1/2))+O(m/r). Equation (9) is verified by exact sums over all one-dimensional kernels in small fields.

### 4.3 Adaptive linear queries

For completeness, after the slopes have been revealed an adaptive linear sketch of the translations has an exact affine-subspace description on each fixed transcript. If the sequence of queries has rank d on that transcript, the leaf is a coset with uniform translation mass p^-d and rowspace V of dimension d. The relative Q_T^+-minus-R^+ mass on the leaf, relative to the fresh uniform mass, is bounded by the sum of the magnitudes of the Fourier coefficients in V. By (4) and the pivot argument, now using the L1 rather than L2 form of Parseval,

    |Q_T^+(leaf|t)-R^+(leaf|t)|
       <= p^-d alpha^(m-2) product_{j in J}(p/sqrt(w_j)).

Summing over disjoint leaves and using d<=q bounds any bounded transcript score's fixed-target gap by alpha^(m-2)(p/sqrt(w_min))^q+4m/r. This rules out sufficiently short adaptive linear-query algorithms, e.g. q<=floor((m-2)/2) at the stated widths. It is weaker than (8), and does not cover arbitrary nonlinear queries of the full translations.

## 5. Exact all-cut spectrum and point-to-point communication

Let s=(r-1)/2 and choose one representative chi_j from each inverse pair of nontrivial group characters. For the full-group marking extension define

    g_{i,j}(seed_i)=(1/r) sum_{P in G} chi_j(P)L_{x(P)}(seed_i).

These are real since x(P)=x(-P); they have mean zero. Their exact Gram matrix is

    <g_{i,j},g_{i,k}> = (d_i/r) (1_{j=k}-1/r),
    d_i=p(p-w_i)/(w_i(p-1)).

Proof: the constant part of the incidence Gram drops out since each nontrivial character sums to zero. Group the other part by abscissa. A fiber is a pair P,-P, except for the singleton sentinel O. The sum of products of the fiber character sums equals r 1_{j=k}-1, by group character orthogonality and subtracting the one double-counted identity. This proves the formula.

Writing c_j=chi_j(T)+chi_j(-T), the density difference on the full seed vector is

    D_T=sum_{j=1}^s c_j product_i g_{i,j}.               (10)

This is character orthogonality on the relation constraint and cancellation of the trivial character by the matched density. Put A=product_i (d_i/r). For a subset I of k slots, the normalized Gram of the functions product_{i in I} g_{i,j}, divided by product_{i in I}(d_i/r), is

    B_k=[(1-1/r)^k-(-1/r)^k] I_s + (-1/r)^k J_s.

It is the k-fold entrywise product of B_1=I-J/r. Its eigenvalues are the coefficient of I, and that coefficient plus s times the coefficient of J. In particular 0<B_k<=I; the upper bound also follows inductively from the Schur product theorem, as B_1<=I. For k>=2 both eigenvalues are 1+O(k/r) when k=o(r), while for k=1 they lie in [1/2,1]. These matrices completely determine the singular spectra across every bipartition in (10).

For T!=O, sum_j c_j=-1 and sum_j c_j^2=r-2 by character orthogonality and 2T!=O. Hence the squared L2 norm is exactly

    ||D_T||_2^2 = A [(r-2)(1-1/r)^m+(3-r)(-1/r)^m].     (11)

At the relevant widths, m=O(sqrt(b)), lambda=r^(m-1)W/p^m, and each w_i/p small, (11) equals (1+o(1))/lambda. In particular a uniformly small Fourier coefficient or a small projection advantage does NOT mean the complete seeds are information-theoretically uninformative. The total energy is large at sparse tuning. The real rank across any nontrivial cut is exactly s for nonzero T, because both Gram matrices are positive definite and c_j!=0 for odd r. The operator norm of the flattening is at most 2 sqrt(A) by the bound on its two Gram matrices. Thus every score with modulus <=1 and real separated rank K across a cut has advantage at most

    2 sqrt(K A) + 4m/r.                                (12)

This follows from the Frobenius norm <=1 of the score kernel and the nuclear norm <=sqrt(K) of any rank-K such kernel. A randomized mixture is handled by fixing its coins and averaging the gap. At balance A=(1+o(1))/W. This is an exact spectrum for the field-affine law's full-group extension and a stable bound on the actual law; the uploaded paper already has a separated-approximation result for its original dyadic sampler. We do not claim that all rank-based restrictions are new.

**Communication consequence.** Two coalitions on opposite sides of any cut can process their seeds using arbitrary computations and an interactive binary protocol. Suppose at most B bits cross the cut (including all control information), and a designated party returns a score. Fix private/public coins. Each transcript defines a rectangle of possible inputs, and on that rectangle the output is a function of just the designated party's inputs. Thus the bounded score has separated rank <=2^B (or <=2^(B+1) if one uses a crude count of short binary strings). By (12), a planted advantage delta>4m/r requires

    B >= log_2(1/A)+2 log_2(delta-4m/r)-O(1).            (13)

In particular B>=b-o(b) at balanced widths and delta=2^-o(b). This holds even for adaptive rounds and unlimited local computation, provided the usual rectangle property of the charged transcript holds. It is a communication bound, not a RAM-time lower bound.

In a fixed-schedule private point-to-point m-party protocol with one local seed per party, the total communicated bits C satisfy sum_i B_i=2C, where B_i is the actual schedule's bit count crossing the singleton cut i. Applying (13) to each cut gives C>=(1/2-o(1))mb. One must not sum unrelated worst-case per-cut maxima for an adaptive schedule. A variable-schedule version is possible under the following explicit model: for each cut the other parties can simulate their internal conversation freely, a complete cut transcript is an ordinary binary two-party transcript with routing/termination charged, and the lengths B_i(x) sum to at most 2C+O(m) on every input. Partition the score by the least i whose cut uses <=floor(2C/m)+O(1) bits; such an i exists on every input. For any fixed short cut transcript, the complement can reconstruct the full communication history (every message has an endpoint in the complement), and decide the chosen-i event. Consequently the masked score for cut i still has rank at most 2^(floor(2C/m)+O(1)); the output, if at i, is multiplied by a complement-only indicator on that rectangle. Sum the m bounds (12), applying the sentinel conditioning error only once to the complete score. This yields

    delta <= 2m sqrt(2^(2C/m+O(1)) A)+4m/r,

and the same leading lower bound. This variable-schedule result is stated only for protocols with these cut simulations and no free timing, free broadcast, or hidden routing information. Ordinary RAM access to all seeds, tensor contraction as a single centralized computation, and arbitrary approximate modular certificates are not constrained in this way. In particular (13) must NOT be advertised as an ECDLP time lower bound.

## 6. A constructive higher-order primitive: interval incidence in polynomial time

The support-equivalence result is limited to second order. Higher-order incidence contains actual support geometry, and one can compute it very cheaply for a special but completely general-field sampler. Take S={0,...,w-1}, with 2(w-1)^2<p and w>=2. This condition holds at balanced widths when m grows. This support is a valid two-transitive affine sampler with the same previous coverage guarantee; it is not automatically compatible with the previous thin-digit Fourier ledger, which is not needed for this incidence primitive.

For k prescribed abscissas x_1,...,x_k, discard repeats. If there is one distinct value, probability of retention is w/p (zero prescribed points: probability one). If at least two, fix distinct x_1,x_2, and put tau_j=(x_j-x_1)/(x_2-x_1). For an arbitrary S, the exact probability that the mask contains all these values is N_S(tau)/(p(p-1)), where

    N_S(tau) = #{s_1,s_2 in S, s_1!=s_2:
                    s_1+tau_j(s_2-s_1) in S for every j}.

Each ordered pair determines a unique nonzero affine map; hence no multiplicity error. For three distinct points this is the familiar triple ratio kernel. Unlike the pair Gram matrix, it depends on the support.

**Interval theorem.** For S=[0,w-1] under the size condition above, compute every tau_j as a rational reconstruction A_j/B_j with gcd(A_j,B_j)=1, |A_j|<=w-1, and 1<=B_j<=w-1. If any reconstruction fails, N_S=0. Include tau_1=0=0/1 and tau_2=1=1/1. Let B=lcm_j B_j and scale each numerator a_j=A_j(B/B_j). If B>=w, N_S=0. Otherwise let R=max_j a_j-min_j a_j. If R>=w, N_S=0. In every remaining case put K=floor((w-1)/R). Then

    N_S(tau) = 2 K w - R K(K+1).                        (14)

The denominator lcm and numerator span are ordinary integers, not field residues. Formula (14) takes polynomial time in b and k, independently of w.

Proof: any accepted pair has integer difference d=s_2-s_1 with 0<|d|<=w-1. For each j, e_j=s_j-s_1 also has absolute value <=w-1, and e_j/d equals tau_j mod p. Reduction to lowest terms gives a bounded rational reconstruction. Such a reconstruction is unique: the absolute value of the cross-multiplied difference of two bounded candidates is at most 2(w-1)^2<p. Thus their congruence implies equality as rational numbers. We must have B|d, and writing d=B v for a nonzero integer v, every s_j=s_1+v a_j. Conversely these values give a legal pair and all its other images when they lie in the interval. The number of legal s_1 for a given v is max(0,w-|v| R). Summing over the two signs and v=1,...,K proves (14). R>=B>=1 because 0 and B are among the a_j.

For efficient reconstruction of tau=t mod p, write t B-A = h p with integer h and gcd(h,B)=1 (since gcd(A,B)=1 and B<p). Then |t/p-h/B|=|A|/(p B)<1/(2B^2); Legendre's continued fraction criterion shows h/B is a convergent of t/p. Enumerate its O(log p) convergents, recover A=tB-hp, and verify the bounds and coprimality. Handle zero directly. All operations are on polynomial-bit integers; stop lcm construction once B>=w. This supplies a genuine algorithm rather than enumeration of w(w-1) ordered pairs.

For example, for three distinct abscissas in affine arithmetic progression with tau=2, R=2 and N_S is about w^2/2. A generic triple with no bounded rational ratio has zero joint retention. More generally every nonzero higher-order incidence is described by an integer direction of small span, and the same closed form applies however many prescribed points there are. This is the first constructive support-specific primitive in this checkpoint. It prices each term of a higher-order likelihood moment without expanding an interval or computing elliptic discrete logs. It does **not** sum the enormous family of relation tuples occurring in those moments; doing that collectively is the remaining algebraic task.

## 7. Why this does not yet beat the exponent

The old independent-list smoothing result could not be transferred to the actual affine sampler. Equations (1)-(14) establish a different, exact analysis on the actual field-affine masks. Arbitrary support choices are equivalent at second order, target-independent linear compression loses the useful signal, and the joint likelihood is genuinely high-dimensional. The new interval formula gives a precise place where nonlinear support geometry differs. None of these statements restricts a centralized nonlinear field algorithm reading all 2m seed words, or a modular contraction with large real rank.

For orientation, the missing direct evaluator can be formulated as an affine-grid product: let F_{m,T}(x) be the fixed-target summation hypersurface polynomial, whose existence and degree are established in the project's chord-geometry derivation. Given seeds, evaluate product_{s_i in S_i} F_{m,T}(a_1+t_1 s_1,...,a_m+t_m s_m) at one target, including restricted subgrids for witness extraction. Its zero is a relation over the algebraic closure. For 2T!=O, rational abscissas give rational or twist points; the same mixed-twist estimate as the preceding abscissa compiler makes false witnesses exponentially rare on balanced independent affine masks. Exact verification is still required. This is a useful *target operation*, not a new fast algorithm: literally evaluating the product enumerates W choices. The dyadic block identity for interval supports merely expresses the product as two smaller products at different translations; it does not make them the same subcomputation. No claimed saving comes from that identity. Under the inherited outer reduction, an evaluator with log bit time O(b/m) at m~(log b)^(1/3), with charged construction and subgrid access, would beat the exponent. I do not have such an evaluator.

A natural next target is a nonlinear joint evaluation that exploits the small-span rational structure (14) without constructing a full high-order relation tensor. The higher-order incidence result is not proof such an evaluation exists. The earlier candidate coefficient 11.31486814 and the saving power 1/4 remain unchanged. The strict exponent goal is open.

## 8. Verification and scope

Run `OPENBLAS_NUM_THREADS=1 python -B verify.py` in a Python environment with NumPy. The verifier enumerates actual toy prime-order curves, constructs all incidence Gram matrices by integer counts, and forms the actual seed-likelihood numerator tensors by summing over group tuples. It checks the additive Fourier coefficients through integer phase histograms (equality modulo the prime cyclotomic polynomial), the conditional-translation formula by integer convolution, and the exact squared norm by rational arithmetic. It checks the sentinel conditioning TV error and the random projection identity by exact integer coset sums. Numerical singular value comparisons are explicitly floating-point checks, not exact certificates; the Gram and spectrum theorems have the algebraic proof above. It then exhaustively compares bounded rational reconstruction against literal rational fraction enumeration on seven fields and checks higher-order interval incidence by literal pair enumeration, including repeated abscissas, planted masks, and random negatives. See verification.txt for the actual successful transcript. The finite checks are not end-to-end ECDLP tests or evidence of a new asymptotic runtime.
