# Direct group incidence and the cost of sampling its Fourier formula

Fix the mixed sampled slots $F_i\subset E(\mathbb F_p)$, the signed address $a_i(P)$ within slot $i$, and $D=\prod_i |F_i|$. Each legal point has a distinct address within its slot. For a binary challenge $r=(r_1,\ldots,r_m)$, let

$$
h_T(r)=\sum_{P_1+\cdots+P_m=T}\prod_i(-1)^{r_i\cdot a_i(P_i)},\qquad P_i\in F_i .
$$

On a target-regular grid the determinant zero criterion and Fermat give $M_0(r)-M_{p-1}(r)=h_T(r)\bmod p$. If there is a single witness with address $a_*$, the integer $h_T(r)=(-1)^{r\cdot a_*}$. Here we study an exact formula that bypasses the determinant recurrence. The error bounds below concern specific ways to approximate this formula. [Determinant criterion and regularity](../nonlinear_moments/signed_spectrum.md)

## Exact formula and uniform sampling

Write $G=E(\mathbb F_p)$, $N=|G|$, and $\widehat G$ for its group of complex characters. Character orthogonality gives

$$
h_T(r)=\mathbb E_{\gamma\in\widehat G}Y_\gamma(r),\qquad
Y_\gamma(r)=\overline{\gamma(T)}\prod_i A_i(r_i,\gamma),\qquad
A_i(r_i,\gamma)=\sum_{P\in F_i}(-1)^{r_i\cdot a_i(P)}\gamma(P).
$$

This identity works on every grid for the group-incidence quantity. Grant access to every character and its slot sum at unit cost. The full average still uses $N$ terms. On *every* fixed singleton, sampling one uniform character has complex mean squared error exactly $D-1$, averaged over the fresh challenge and the draw. Indeed, for each fixed $\gamma$, binary character orthogonality gives $\mathbb E_{r_i}|A_i|^2=|F_i|$. The independent challenge blocks give $\mathbb E_r|Y_\gamma|^2=D$, while $|h_T|=1$.

## Four-term relations in the sampled masks

**Translate lemma.** Suppose $E:y^2=x^3+Ax+B$ is nonsingular over $\mathbb F_p$, $p>5$. Fix four offsets $d_j\in\mathbb F_p$ whose multiset cannot be partitioned into two equal-value pairs. At most 16 values of $C\in\mathbb F_p$ admit four geometric points $P_j$ with $x(P_j)=C+d_j$ and $\sum_jP_j=O$.

The fourth summation polynomial $S_4$ has degree at most four in each variable, and its vanishing characterizes such a sum over the algebraic closure. Its restriction $S_4(C+d_1,\ldots,C+d_4)$ has degree at most 16. To see it is a nonzero polynomial, suppose the relation holds at the generic $C$. Substitute $C=t^{-2}$ into a Laurent series field. For each offset the two lifts have opposite $y$; choose the signs $\epsilon_j$ giving the identity. Pull back the invariant differential $dx/(2y)$, which is additive under the group law. Divide out the common nonzero factor $dC\,t^3/2$. The vanishing sum becomes

$$
\sum_{j=1}^4\epsilon_j\left[1-\tfrac32 d_jt^2+
  \left(\tfrac{15}{8}d_j^2-\tfrac12 A\right)t^4+O(t^6)\right]=0 .
$$

The constant forces two plus signs and two minus signs. The next two coefficients equate the sums and sums of squares of the two pairs of offsets. Since $2,3,5$ are invertible, the pairs have the same sum and product and are the same multiset. This contradiction proves the lemma. [Semaev, *Summation polynomials and the discrete logarithm problem on elliptic curves*](https://eprint.iacr.org/2004/031)

Let $M=2^b$, $p\le M<2p$, and let a slot have $w$ raw abscissa addresses before discarding out-of-field values or abscissas without points. Call a slot bad if some four signed points sum to $O$ with an unpaired abscissa multiset. In a fixed-bit slot the common fixed bits give an integer $C$ uniform on $M/w$ possibilities. The $w$ free-bit patterns give fixed offsets $d$, with canonical abscissa $C+d$ when it is below $p$. For each ordered four raw addresses that could give an unpaired multiset, at most 32 fixed seeds work by the lemma, since an element of the field has at most two representatives in $[0,M)$. If different raw addresses have equal offsets modulo $p$, they cannot both be canonical at the same seed. Thus the bad probability is at most $32w^5/M$.

For the odd affine anchor, condition on its odd multiplier $\alpha$. Inverting the mask gives $x=(C+(\alpha^{-1}u\bmod M))\bmod M$ for each retained raw hash address $u$, with $C$ uniform on $[0,M)$. For any four addresses these abscissas are simultaneously ordinary integer translates on at most five intervals in $C$, according to the four wrap points. On each interval either the offsets pair modulo $p$ and the relation is irrelevant, or at most 32 choices work. This gives bad probability at most $160w^4/M$. The union bound for one affine and $m-1$ fixed-bit slots is $O(mw^5/M)$.

For a good slot of size $s_i=|F_i|$ its additive energy satisfies

$$
\#\{(P,Q,R,S)\in F_i^4:P+Q=R+S\}\le12s_i^2 .
$$

Indeed the four abscissas pair, since changing the signs of $R,S$ preserves their abscissas. For each of three pairings choose one signed point from each pair, then at most two signs for each partner. This counts torsion points as well.

## Even an optimized independent importance sampler has large mean squared error

Condition on any fixed good singleton and a fixed challenge. Let $L(r)=N^{-1}\sum_\gamma |Y_\gamma(r)|$. For a draw from a proposal $\pi$ the usual unbiased term is $Y_\gamma/(N\pi_\gamma)$, with positive sampling probability wherever $Y_\gamma\ne0$. Its conditional second moment is at least $L(r)^2$ by Cauchy–Schwarz, with equality at $\pi_\gamma\propto |Y_\gamma|$. The proposal can use all the public seeds, the target, and the actual challenge.

To bound $L(r)$, set $k_i(\gamma)=\mathbb E_{r_i}|A_i|^4/s_i^2$. The second-moment identity gives $k_i\ge1$. Expanding the fourth power and averaging over $\gamma,r_i$ counts ordered quadruples with equal pair sums and an additional binary address constraint. The energy estimate gives $\mathbb E_\gamma k_i\le12$. Interpolation between first and fourth moments gives $\mathbb E_{r_i}|A_i|\ge\sqrt{s_i/k_i(\gamma)}$. Independence of the challenge blocks and Jensen's inequality now give

$$
\mathbb E_rL(r)\ge\sqrt D\ \mathbb E_\gamma
       \exp\!\left(-\tfrac12\sum_i\log k_i(\gamma)\right)
\ge\sqrt D\ \exp\!\left(-\tfrac12\sum_i\mathbb E_\gamma\log k_i\right)
\ge\sqrt{D/12^m}.
$$

Thus an average of $Q$ independent unbiased terms has mean squared error at least $(D/12^m-1)/Q$, even if determining and sampling an optimal proposal is free. At the thinned widths with $m\to\infty$, $m=O(\sqrt b)$, and $\log_2 w=b/m+O(1)$, the probability of any bad slot is $2^{-b+o(b)}$. Under the sampler's free-bit run condition $m2^s/w=o(1)$, conditioning on the fixed-target guarded singleton $U$ costs a factor at most $2^{O(m+\log m)}$. Conditional on $U$, existing size and regularity estimates give $D\ge 2^{b-O(m\log m)}$ on regular grids with probability $1-O(m^{-2})$. Consequently constant mean squared error on these grids requires $Q\ge 2^{b-O(m\log m)}$ for this estimator. A sign decoder with rare, very large errors could have good prediction probability despite high mean squared error. Correlated quadrature, nonlinear decoding and algorithms using another representation also fall outside this result. [Sampling widths](../sampled_slots/isolation.md) · [Conditional size bound](../nonlinear_moments/signed_spectrum.md)

## Bounded scores from a common feature space

There is also an exact restriction across targets. For one fixed seed write $H_T=\#\{a:\sum_iP_i=T\}$. Orthogonality of distinct address characters and disjointness of target fibers give $\mathbb E_r h_T(r)h_S(r)=\mathbf1_{T=S}H_T$. With uniform challenge weight, the nonzero singular values of the target-by-challenge matrix are $\sqrt{H_T}$.

Suppose $n$ eligible targets are uniform and independent of the seed. For each seed allow an arbitrary common $k$-dimensional subspace $V$ of square-integrable challenge functions. For each target allow any $f_T\in V$ with $\|f_T\|_2\le1$. Distinct singleton targets have distinct witness address characters $\chi_{a_*}$. Thus

$$
\mathbb E_{\sigma,T}\big[{\bf1}_{U_T}|\langle f_T,\chi_{a_*}\rangle|^2\big]\le k/n.
$$

For a fixed seed bound each term by $\|\operatorname{Proj}_V\chi_{a_*}\|_2^2$, then sum over the orthonormal witness characters before averaging over seeds. For singleton probability $q$ this gives conditional squared correlation at most $k/(nq)$. Target-dependent subspaces, unrestricted nonlinear transformations of features and models where the seed depends on the target need separate arguments.

## Sign decoding from a single Fourier phase

One natural nonlinear decoder takes the sign of the real part of a single Fourier term. For a target-independent choice of character, even the **prediction correlation** of this decoder admits an upper bound. More generally, let the seed and shared randomness determine $K$ orderings of the $2^d$ challenges, independently of the target. For each target the predictor can choose any one ordering and any union of $J$ intervals in it, returning $+1$ inside the union and $-1$ outside. Choices of the order and interval endpoints can depend arbitrarily on the target and seed. Then

$$
\mathbb E_{\sigma,T}\big[|\langle f_T,\chi_{a_*}\rangle|^2\mid U_T\big]
\le \frac{K(16J(d+1)^2+2)}{n q}.
$$

For one order, partition its positions into dyadic blocks at each of the $d+1$ scales. Any union of $J$ intervals decomposes into at most $2J(d+1)$ disjoint blocks. For each block $I$, Bessel's inequality over the distinct singleton characters bounds the sum of squared coefficients of $1_I$ by its squared norm $|I|/2^d$. Summing this bound over all blocks gives $d+1$. For each target separately, Cauchy–Schwarz on its selected blocks contributes a factor at most $2J(d+1)$. Finally $f_T=2\,1_{\rm union}-1$, so $|2a-b|^2\le8|a|^2+2|b|^2$ and Bessel for the constant function give the displayed numerator for one order. Sum over the $K$ orders, divide by $n$, and condition on $U_T$. Shared randomization preserves the bound by averaging.

For a single character, order the challenges by the argument of $\prod_i A_i(r_i,\gamma)$, placing zero products in a separate block. Multiplication by $\overline{\gamma(T)}$ rotates the decision semicircle. Its preimage is a union of at most three intervals including a possible uniform choice at zero. The same argument allows a target-specific choice among a menu of $K$ seed-dependent characters. For $n=2^{b+O(1)}$, $q\ge2^{-O(m+\log m)}$, $d=O(b+m)$, and $KJ=2^{o(b)}$ the conditional mean squared correlation is $2^{-b+o(b)}$. This is an average over independent uniform targets and seeds conditioned on success. Combining several Fourier terms *before* taking the sign, allowing an unrestricted number of phase intervals, or choosing from all $N$ orders needs another argument.

## Thresholding a sum of several common characters

Suppose a seed determines $Q$ complex Fourier features on challenges. A target can use any affine halfspace on the $2Q$ real and imaginary parts, including the sign of a weighted sum of these features. The resulting challenge decisions for one seed have VC dimension at most $v=2Q+1$. More generally, any seed-dependent family of binary decisions with VC dimension at most $v\ge1$ satisfies, for $0<\tau\le1$,

$$
\Pr_{\sigma,T}\!\left(|\langle f_T,\chi_{a_*}\rangle|\ge\tau\mid U_T\right)
\le \frac{4(20e/\tau)^{2v}}{\tau^2nq}.
$$

Fix the seed and take a maximal set of decisions pairwise separated by probability of disagreement at least $\delta=\tau/4$. If its size $P$ exceeds one, take $h=\lceil3\log(P)/\delta\rceil$ uniform challenges. With positive probability they distinguish every pair, by the union bound. For $h<v$, $P\le2^v$. Otherwise Sauer's lemma gives $P\le(eh/v)^v$. Since $h\le5\log(P)/\delta$, writing $t=\log(P)/v$ yields $e^t\le5et/\delta$. Using $t\le e^{t/2}$ gives $P\le(5e/\delta)^{2v}$. This bound also covers $P=1$. [Sauer, *On the density of families of sets*](https://www.sciencedirect.com/science/article/pii/0097316572900192)

Each target predictor disagrees on at most $\delta$ with one center, so correlation of magnitude at least $\tau$ with a witness character implies magnitude at least $\tau/2$ for that center. At most $4/\tau^2$ distinct witness characters can meet that condition per center by Parseval. Count singleton targets at each seed, average, then divide by $nq$. Thus achieving constant correlation on a constant fraction of the successful seed-target pairs requires $Q=\Omega(b)$ common characters in the large-target regime above. Inverse-polynomial correlation on such a fraction still requires $Q=\Omega(b/\log b)$. These bounds permit polynomial-sized common feature sets. A choice of frequencies depending on the target can fall outside the common-feature model.

The finite checks enumerate translate roots over two primes, the full affine and fixed-bit slot laws at width four over $\mathbb F_{101}$, and all signed characters for one regular singleton over $\mathbb F_{1031}$. The latter has 384 legal tuples, 346 represented targets and 512 challenges. Uniform sampling gives one-draw mean squared error 383 and the optimal challenge-dependent proposal gives 109.559381. The sum of best squared correlations from one arbitrary interval in a fixed phase ordering, over 308 singleton targets, was 3.582672. [Check](../../checks/check_direct_fiber.py) · [Recorded output](../../evidence/check_direct_fiber.txt)
