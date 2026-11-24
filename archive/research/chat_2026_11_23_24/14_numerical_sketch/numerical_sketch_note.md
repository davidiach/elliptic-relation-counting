# Adaptive numerical sketches and full-value moment smoothing on elliptic planted lists

Research checkpoint, 24 November 2026. Continues `prime_exponent_note.md` and the original v3 preprint `paper (1).pdf`. The strict general-prime exponent is **not improved**. The results below apply to the previous checkpoint's uniform signed-list planted experiment, not the original dyadic affine hash or the later field-affine digit sampler. They provide information/query-model obstructions, not general circuit or time lower bounds. The numerical coefficient 11.31486814 remains an earlier conditional candidate; it is not re-audited here. No novelty claim against the literature is made.

## 1. Setting and scalar character estimate

Let E/F_p be a smooth short Weierstrass curve with G=E(F_p) of odd prime order r, b=ceil(log2 p), and p>3. Take m slots, each with w underlying labeled points A_ij. There are n=mw independent uniform points in G under the fresh law P. Each underlying point supplies both signs. Under Q_alpha, for a signed address alpha selecting one label in each slot, condition the signed sum of those m points to be a prescribed target T. Under Q choose alpha uniformly among W=(2w)^m addresses and then use Q_alpha. Every address has the same conditioning probability 1/r. Every proper subset of the selected coordinates is still independently uniform, since an unobserved selected coordinate can complete the sum uniquely. Unselected points are independent of the selected tuple under Q_alpha. We use the whole group, including O, throughout the probability arguments.

A feature dictionary phi=(f_1,...,f_l) consists of rational functions in L(DO), linearly independent modulo constants. Extend their values at O by zero. For the elliptic evaluation code, use the m nonconstant monomials in L((m+1)O), so l=m and D=m+1. The evaluation row v(P) is phi(P)-phi(-T). For the signed row v(-P), each monomial in x,y differs from its value at P by its public sign, plus the same constant translation. With the chosen extension at O this remains an algebraic identity; the extended rows are artificial at O and agree with the intended rows on regular inputs.

Let psi(z)=exp(2 pi i z/p). For a nonzero coefficient vector c in F_p^l, the nonconstant rational function c.phi has a sole pole at O, of positive order <=D<p. It is not geometrically of the Artin--Schreier form h^p-h+constant, since a pole of such a nonconstant function has order divisible by p. The Bombieri--Weil curve-sum bound (genus one, a single pole) yields

    |E_A psi(c.phi(A))| <= eta := ((D+1) sqrt(p)+1)/r.

The +1 accounts for the artificial value at O. By Hasse, r~p. In particular, if D=p^o(1), then eta<=p^(-1/3) for sufficiently large p. We use this conservative inequality below. See Bombieri, *On exponential sums in finite fields*, American Journal of Mathematics 88 (1966), 71--105, Theorem 5; a statement with the needed single-pole degree bound is given in Section 2.2 of `Tracing AG Codes: Toward Meeting the Gilbert-Varshamov Bound`, ICALP 2026, https://drops.dagstuhl.de/storage/00lipics/lipics-vol374-icalp2026/html/LIPIcs.ICALP.2026.68/LIPIcs.ICALP.2026.68.html .

## 2. Adaptive full-value linear measurements

A scalar measurement may be an arbitrary field-linear combination of all feature values of all underlying points, plus a known constant. An algorithm can choose its coefficients using the target, public parameters, private coins and all preceding numerical answers. It may process the full field residues by an arbitrary nonlinear algorithm. It has no other access to the hidden points. Constants can be absorbed by translating answers. Pad early termination with zero queries, so there are q scalar field responses. On every fixed transcript y, the coefficient choices specify a fixed simultaneous linear sketch

    Y = sum_j A_j^* phi(A_j) in F_p^q,

where A_j^* is a public q-by-l coefficient block on that transcript (not the underlying point A_j). The event that the adaptive transcript is y equals the event that this fixed sketch has value y, with constants subtracted along the path. No efficiency claim about enumerating the transcripts is required.

**Theorem 1.** Set K=l(m+12),

    theta = min(1, (K q/(mw))^m),
    epsilon = (1+(p-1) eta^12)^q - 1.

Assume eta<1 and epsilon<1. Then the total variation distance between the complete fresh and planted adaptive transcript laws is at most

    1 - (1-theta)(1-epsilon)/(1+epsilon)
    <= theta + 2 epsilon/(1+epsilon).

The probability of recovering the exact planted address is at most

    theta + 2 epsilon/(1+epsilon)
        + (1+epsilon)/((1-epsilon) W).

Both bounds may be capped at one. The proof holds for arbitrary randomized algorithms by fixing and averaging their private coins, and for arbitrary deterministic or randomized readouts. No complexity restriction is imposed on the decoder or on its selection of the linear coefficients. Its data access, however, is restricted to these scalar measurements.

### 2.1 A linear-matroid peeling lemma

For a fixed sketch, regard all scalar columns of all blocks A_j^* as individually labeled vectors in F_p^q. Let their span be V of dimension t. If K disjoint bases of V do not exist among these columns, Edmonds' matroid base-packing theorem gives a set of columns of span U such that

    number of columns outside that set < K(t-dim U).

In particular the number outside the subspace U is at most this count. Remove every *whole block* containing a column outside U, putting its label into B. At most K(t-dim U) blocks are removed, and the remaining span is contained in U. Repeat with the remaining span. The dimension drops at each failure, and summing the drops proves

    |B| <= K q.

The process terminates with zero span (vacuously packed) or a remaining span V with K disjoint scalar-column bases. The proof is existential, sufficient for a query lower bound. It does not charge the hidden-input solver for constructing B. The packing criterion used here is the standard consequence of the rank formula for matroid union; e.g. MIT 18.438, Lecture 13, Theorem 4, https://math.mit.edu/~goemans/18438F09/lec13.pdf .

### 2.2 Relative, not just absolute, smoothing

After peeling, condition on the removed blocks and on all m selected points. At most m underlying blocks of the remaining sketch are now fixed by the plant. Each block has at most l scalar columns, and packed bases are column-disjoint. Consequently at most ml bases meet selected blocks. At least 12l disjoint bases remain, all formed entirely from independent uniform unselected points. They span the same V. Under the fresh law one can use the same argument without the extra conditioning.

For any linear functional lambda in V^*, the character of the independent residual sum factors over its blocks. If h(lambda) is the number of residual blocks with at least one column pairing nontrivially with lambda, the magnitude of its Fourier coefficient is at most eta^h(lambda). For each of the k'=12l surviving bases let h_i(lambda) count the active scalar columns in that basis. Across all bases a block contains at most l used columns, so sum_i h_i<=l h. Thus

    sum_(lambda!=0) |hat mu(lambda)|
        <= sum_(lambda!=0) product_(i=1..k') eta^(h_i(lambda)/l)
        <= product_i [sum_(lambda!=0) eta^(k' h_i(lambda)/l)]^(1/k')
        = (1+(p-1)eta^12)^t-1
        <= epsilon.

The second inequality is Holder's inequality. In the equality, each basis identifies V^* bijectively with F_p^t; h_i is ordinary Hamming weight. Fourier inversion gives for every z in V the pointwise bounds

    (1-epsilon)/p^t <= mu(z) <= (1+epsilon)/p^t.

This is the reason the result extends to adaptive transcript probabilities. A mere TV estimate on each fixed sketch would not justify summing over all possible adaptive paths. Fixed selected contributions outside B only translate the distribution within V.

### 2.3 Compare adaptive transcript probabilities

On a fixed transcript path y, let B=B(y) be the peeled set constructed from its coefficient matrix. Define the reference distribution U_y as the sum of the removed blocks with independent uniform inputs, plus an independent uniform vector in the surviving V. The subscript indicates the fixed query sequence, not a probability distribution on adaptive paths. Let u_y be this reference distribution's mass at the adjusted answer y.

Under P the fixed sketch mass p_y lies in [(1-epsilon)u_y,(1+epsilon)u_y]. For Q_alpha the same pointwise interval holds whenever alpha does not have all its selected underlying labels in B. Indeed, B then omits at least one selected point and its input distribution is exactly the fresh product law. Given its inputs and all selected points, residual smoothing applies, so average over the conditioning. We call these addresses good for this transcript. This designation depends on the transcript's coefficients and alpha, and need not define a single fixed subset of addresses across transcripts.

For any fixed B, the fraction of addresses entirely in B is

    product_i (|B intersect slot_i|/w) <= (|B|/(mw))^m <= theta,

by AM--GM (the 2^m choices of signs cancel). Hence for every y the total planted mass contributed by good addresses is at least

    Q_good(y) >= (1-theta) (1-epsilon)/(1+epsilon) P(y).

In particular this is a pointwise lower bound for Q(y). Summing the negative part of Q-P proves the TV inequality. No independence between alpha and the actual transcript, and no union bound over its exponentially many values, is used. Summing the same lower bound proves that the total planted probability of the covered-address event (bad for its own transcript) is at most 1-(1-theta)(1-epsilon)/(1+epsilon). On a good transcript with a proposed address alpha_hat, its joint planted mass with the correct hidden address is at most P(y)(1+epsilon)/[W(1-epsilon)], by the upper bound on Q_alpha for good addresses. Sum over y and add the covered-address mass. This proves recovery. If a solver is judged by returning *any* relation rather than the planted address, add the preceding checkpoint's probability of an extra relation or irregularity. That probability tends to zero at the thinned widths there.

### 2.4 Consequences and scope

For the ordinary evaluation dictionary l=m, D=m+1, K=m(m+12). At m growing with m=O(sqrt b), w=floor(r^(1/m)/16), and q<=w^(1-gamma) for a fixed gamma>0,

    theta <= 2^(-gamma b+O(m log m)) = 2^(-gamma b+o(b)),
    epsilon <= (1+p^(-3))^q-1 = O(q/p^3) = 2^(-3b+o(b)).

Any bounded score based on such adaptive numerical transcripts has gap at most this TV (or twice TV for scores in [-1,1]). To achieve inverse advantage 2^o(b) requires q>=w^(1-o(1)) in this asymptotic sense. For fixed full numerical determinants of linear combinations of the signed rows, each k-by-k matrix (k<=m) can be exposed using at most m^2 scalar measurements. Thus any algorithm using J adaptive full determinant values is simulated by q<=m^2 J scalar measurements with exactly the same outputs. The lower bound transfers with that factor. The previous zero-only interface had a much stronger 2^(b-o(b)) query obstruction; numerical values do carry more information. Querying each raw row individually is allowed by the direct-solver budget. This result therefore does not eliminate a direct arithmetic algorithm reading the full input and using joint nonlinear operations.

## 3. Entire vectors of determinant moments can be uninformative

Let d=m+1 and let k0 be the highest determinant-moment order. On a regular signed-list instance, for each slot separately the integrand Delta^j has pole order <=jd as a function of its point. Summing the two signs cancels odd powers of y; replacing y^2 by x^3+ax+c leaves a polynomial in x of degree at most floor(jd/2). This is the usual pole-weight argument from the uploaded paper, Section 5, and the affine-digit note. It is valid for the rational points in these signed uniform lists; no quadratic-residue filtering is needed because the underlying samples are actual group points. Put s=floor(k0 d/2). All moments Z_0,...,Z_k0, jointly, are deterministic functions of T and the m vectors

    U_i = (sum_(j=1..w) x(A_ij)^h)_(h=1..s).

The constant moment of each slot is known (w underlying points, both signs). At O set x(O)=0 and extend the determinant moments by the same polynomial in the U_i; the extension is artificial and agrees with the intended moments on regular inputs. Conditioning on regularity at the end changes TV by at most the irregularity probabilities.

For one underlying random point, every nonzero additive character of (x,x^2,...,x^s) has magnitude <=eta_s=((2s+1)sqrt(p)+1)/r, provided 2s<p. This is Section 1 with pole order 2s. Given the planted address and *all its selected point values*, each slot still has w-1 independent uniform points. Hence by Fourier inversion its U_i has relative pointwise deviation from uniform on F_p^s at most

    epsilon_s=(p^s-1) eta_s^(w-1).

The selected point translates the law and does not change the bound. The unselected random points are independent across slots even under this conditioning. Thus the joint law of the vectors is within m epsilon_s/2 in TV of the uniform law on (F_p^s)^m, by the product TV inequality. The fresh law has the same bound using w independent points. In particular:

* The TV between the *whole* extended moment vector under P and Q_T is at most m epsilon_s, for every fixed T.
* For any decoder of the whole vector, the probability of recovering alpha is <=1/W+m epsilon_s/2, because the joint law of alpha and the sketch is within that distance of an independent uniform address and uniform sketch.
* The TV of the actual moments conditioned on regular inputs is at most m epsilon_s + P(not regular)+Q_T(not regular). The preceding checkpoint gives O(n^2/r) bounds for these irregularities in both laws. For any-relation recovery also add its extra-relation probability.

For the paper's inexpensive moment prefix k0=floor(2^(b/m^2)), at growing arity m=O(sqrt b) and the thinned widths, s=o(w). More generally the argument is effective whenever s<=w/12 eventually and s=p^o(1). Then eta_s<=p^(-1/3) for large p and epsilon_s<=p^(-w/4+1/3). At the thinned widths, the regularity cost is 2^(-b+o(b)). So in this uniform-list model even an unbounded, perfectly trained nonlinear decoder of the *complete affordable moment vector* cannot get a gap 2^-o(b). This conclusion is stronger than a character or scalar-moment test and relies on independent unselected list entries. It does NOT follow for the affine sampler in the actual proposed ECDLP algorithm.

## 4. Arbitrary separate slot summaries: an averaged information bound

The next result needs no algebraic features, character bounds or field-degree restrictions. Let Z_i be an arbitrary randomized summary of slot i, computed independently of the other slots and of the target, with at most s_i bits of output. Shared public randomness is allowed by conditioning on it. Private randomization is independent per slot. Use the fresh lists, mark an independently uniform signed address, and call its sum T. Then T is uniform in G. This joint law has conditional distribution Q_T on lists when the target is fixed. Compare with the fresh lists and an independent uniform T.

**Theorem 2.** The average fixed-target TV distance of the summary laws is at most

    (1/2) product_i min(2, sqrt(2 ln(2) s_i/w)).

The same bound applies to the advantage of any target-aware bounded [0,1] decoder of the separate summaries, averaged in absolute value over targets. For a [-1,1] decoder multiply by two.

Proof: write U for uniform on G and M_i for the independently marked signed point in slot i. Before applying the independent sign let J be its uniform underlying index. By independence of list entries and the mutual-information chain rule (all entropies in bits),

    I(A_J;Z_i) <= I(A_J;Z_i,J)
                 = (1/w) sum_j I(A_ij;Z_i)
                 <= H(Z_i)/w <= s_i/w.

The penultimate inequality follows by subadditivity of conditional entropy and additivity of unconditional entropy of the independent points. Independent signing can only reduce mutual information. For Q_i(z), the conditional distribution of M_i given Z_i=z, Pinsker and Jensen yield

    E_z ||Q_i(z)-U||_1 <= sqrt(2 ln(2) s_i/w),

with the trivial upper bound 2. Conditioned on all summaries the marks remain independent; their sum has law Q_1*...*Q_m. Each delta_i=Q_i-U has total mass zero, so all proper terms of the expansion vanish upon convolution with U. Therefore

    Q_1*...*Q_m - U = delta_1*...*delta_m.

Use the l1 convolution bound and independence of the summaries to take expectations. Half the resulting quantity is the TV between (T,Z) under the marked-sum law and under the independent-target law. It equals the mean over T of the fixed-target TV because the target marginal is uniform in both. Conditioning away T=O costs at most factor r/(r-1). The score claim follows from the defining variational bound on TV. This proof does not apply when a summary is allowed to depend on T in advance or to inspect other slots jointly. In particular it does not constrain all interactive or joint nonlinear computations.

At s_i=2^o(b/m) uniformly, w as before, the right side is 2^(-b/2+o(b)). Thus even arbitrary nonlinear, target-independent summaries that are substantially smaller than the input lists cannot retain a useful average planted advantage in this uniform-list model.

## 5. Why these statements do not solve the affine-sampler problem

The actual field-affine sampler has only two field words of seed entropy per slot. Its entries are correlated, not independent. In fact, for a public formal support S with known size w<p, raw full-support power sums can reveal the entire affine mask from only two or three moments. This is a useful counterpoint, not a relation solver. Write mean sbar and central moments c_j of S, and use the corresponding raw mask moments C_j of a+t S. Then the mask mean is a+t sbar and C_j=t^j c_j. If c_2 and c_3 are nonzero, recover t=C_3 c_2/(C_2 c_3) and then a. If S is a symmetric interval of length 1<w<p with c_2=(w^2-1)/12 nonzero in F_p (exclude exceptional w such as p-1), recover t up to sign from C_2; reflection about the mean produces the same mask. These claims concern the **formal unfiltered affine support**, not the subset with rational elliptic lifts. The original algorithm already receives the seed words. Recovering them supplies no relation finder. It illustrates why iid-list smoothing and the information bound cannot be silently transferred to the structured sampler.

The remaining direct-solver goal is unchanged: joint, nonlinear numerical arithmetic on actual structured seeds at expected cost w^O(1) 2^(2^O(m)), or another fully charged evaluator meeting o(b/(log b)^(1/4)) after the inherited reduction. The present results neither construct such an algorithm nor disprove its existence. In particular they do not give an unconditional ECDLP lower bound, a new time coefficient, or a fast numerical determinant implementation.

## 6. Finite verification

Run `python -B verify.py` (standard library only). The script tests the matroid packing/peeling lemma by exhaustive flat searches and disjoint-basis searches on tiny linear matroids, and exact convolution of numerical sketch distributions on the order-19 elliptic curve over F_13. It calculates an adaptive two-query transcript law using three different second queries selected by the first response. It verifies pointwise relative smoothing on every good signed type and every possible answer, and tests the corresponding transcript and covered-address bounds using exact integers and fractions. It also checks all marked-point translations of the two-power-sum vector and the arbitrary local-summary convolution bound. Large decimal displays are rounded only for readability. The first test run used w=64 as a stopping point and an excessively tight 10^-9 threshold on the toy full-moment distribution: its exact relative error was about 1.90049e-9. The final verifier retains the w=64 value, extends the test to w=128, and passes the unchanged 10^-9 threshold there. No theoretical inequality was changed. See `verification.txt` for the successful transcript. The finite checks are not validation of a general character-sum theorem, a high-dimensional algorithm, or an end-to-end ECDLP runtime.
