#set document(title: "Counting relations in independently hashed elliptic-curve slots", author: "davidiach research repository", date: datetime(year: 2026, month: 11, day: 20))
#set page(paper: "a4", margin: (x: 23mm, top: 22mm, bottom: 22mm),
  footer: context [#set text(size: 8.5pt, fill: rgb("#425450"))
    Hashed elliptic-curve relations #h(1fr) #counter(page).display("1")])
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.66em)
#set heading(numbering: "1.")
#show heading.where(level: 1): h => block(above: 0.3em, below: 0.75em)[#text(size: 18pt, fill: rgb("#17483e"))[#h]]
#show heading.where(level: 2): h => block(above: 0.85em, below: 0.25em)[#text(size: 12pt)[#h.body]]
#let label(x) = text(weight: "bold", fill: rgb("#17483e"), x)
#let grid(cols, ..cells) = block[#set par(justify: false); #table(columns: cols, inset: (x: 6pt, y: 4pt), stroke: none,
  fill: (_, y) => if y == 0 {rgb("#e2ece8")} else if calc.even(y) {rgb("#f4f6f5")}, ..cells)]
#let public = "https://github.com/bochum-crypto/subexponential-ecdl/blob/1ee44e562b1b9a50ea339f3b8f81112b5a24b96f/subexponential-ecdl.pdf"
#let project = "https://github.com/davidiach/elliptic-relation-counting"

#text(size: 9pt, tracking: 0.8pt, fill: rgb("#31584f"))[RESEARCH PREPRINT · 20 NOVEMBER 2026]
#v(8pt)
#block[#set par(justify: false); #text(size: 26pt, weight: "bold")[Counting relations in independently hashed elliptic-curve slots]]
#v(4pt)
#text(size: 14pt)[Conditional complexity, planted distributions, and certificates from determinant moments]
#v(4pt)
Published by #link(project)[davidiach] · AI-assisted mathematical research · Version 1
#v(9pt)
#label[Abstract.] For a fixed point on an elliptic curve, we study the number of ordered representations of that point as a sum from independently hashed slots. A short inverse-hash circuit gives a conditional improvement to the leading time coefficient in a supplied relation-search framework, from 27.01388 to 26.29338. The power of the logarithmic saving is $1/4$ in both bounds. Separately, we give a finite second-moment bound for isolation and an exact likelihood formula for planted relations. At the proposed sparse widths, a predictor that leaves one whole slot unread has error at least $lambda-O(lambda^2)$ even with adaptive access to the other slots. Determinant power sums can be evaluated jointly on a smaller interpolation grid. For three slots of sizes 2, 2 and 4, nine pooled moments certify 145,776 exact counts among 168,480 eligible seed triples on one small curve. Tests on larger pools locate the outstanding difficulty in the decoding step. A square-root-in-bitlength bound remains open.

= The question and its baselines
Let $b=ceil(log_2 p)$ and $L=log_2 b$. Akshima, Bacho, Huang, Kiltz, Nishaburi, Nooripoor, Tittelbach and Wiedijk derive a conditional expected-time bound $2^(O(b/L^(1/4)))$ on prime-order curves over prime fields, with an extension to fixed cofactors. Their version 0.2 assumes a Boolean tensor-contraction result. It gives no numerical leading coefficient for this comparison. A separate, supplied manuscript develops an explicit finite-packet model and reports a coefficient of 27.01388. We compare numerical coefficients only within that framework on the overlapping input class.

#grid((2.65fr,1fr,1fr), [*Source of bound*], [*Time*], [*Space*],
  [Supplied simplex manuscript, v6], [32.15695], [16.07848],
  [Supplied compressed-relation manuscript], [27.01388], [13.50694],
  [Short inverse hash, this project], [26.29338], [13.14669])
#text(size: 9.5pt)[Table 1. Entries multiply $b/L^(1/4)$ in the base-two logarithm. They are asymptotic upper-bound coefficients. None is a measured speedup.]

#pagebreak()
= Where the coefficient comes from
Use a smooth short Weierstrass curve $E$ over $bb(F)_p$, $p>3$, with supplied group order $N=h r$, $r$ an odd prime, $gcd(h,r)=1$, and polynomially bounded cofactor in the conditional comparison. A slot accepts an affine point when its canonical abscissa $x$ satisfies $(alpha x+beta) mod M < tau$, where $M=2^b$, $alpha$ is uniform odd, $beta$ is uniform modulo $M$, and $tau=2^ell$. Different slots use independent seeds. Both signs of a point are retained.

== A shorter operand in the hash circuit
Write $A=alpha^(-1) mod M$ and $C=-A beta mod M$. For a fixed accepted $x$, there is exactly one index $0<=t<2^ell$ with $x=(A t+C) mod M$. At paid arity $m$, $ell=b/m+O(log b+log m)$. Split the long public multiplier $A$ into blocks, reuse the forward transform of the short index $t$, and send one inverse product per block. The hash distribution and the multiplicity of each accepted assignment are preserved.

More precisely, take radix $2^k$ with $k=Theta((log b)^2)$, $n=ceil(b/k)$ and $s=ceil(ell/k)$. Choose transform length $n' asymp n/sqrt(m)$ and an owner-aligned block width $H_0$ just below $n'-s+1$. Then $H_0/n' -> 1$. Each product fits one transform, and $J=ceil(n/H_0)$ blocks need $(J+1)n'$ transmitted words, including the common input transform. With word width $(2+o(1))k$, the leading cost per slot is $(2+o(1))b$. Owner alignment enables local coefficient addition and carry checks. The full accounting of maximum degree, parallel edges and lower-order carries is given in the supplementary short-hash derivation.

The earlier alternating-chord graph, called C38 for its leading word cost, pays $38 m b$. It assigns $4 m b$ to hashing and $34 m b$ to the other arithmetic. The replacement changes the total to $(36+o(1))m b$. Let
$ K = (2 pi sqrt(14 log_2 90))/(9 ln 2). $
Under the finite-packet assumption and the inherited contraction estimate, the log cost of a relation query becomes $(18K+o(1))m b/sqrt(L)$. The outer reduction pays a polynomial multiple of $B C_m+m B^2$, with $log_2 B=b/m+O(log b+log m)$ and $C_m$ the query cost. Balancing gives
$ m tilde.eq L^(1/4)/sqrt(18K), quad log_2 bb(E) T <= (2 sqrt(18K)+o(1)) b/L^(1/4). $
The ratio to the compressed-relation coefficient is $sqrt(18/19)$, a 2.6671% decrease. Relative to the earlier simplex coefficient it is an 18.2342% decrease.

#label[Dependency.] Assumption A.1 of the supplied compressed-relation manuscript requires a uniform finite dyadic packet realization over unramified coefficient algebras. Its graph has square-root maximum degree up to polylogarithmic factors, bounded parallel multiplicity, charged entry access, and a finite parent-to-child identity with the prescribed precision and packing bounds. The transfer also uses the earlier graph encoding, contraction and relation-reduction estimates. The checks in this release verify finite arithmetic and counting identities. The complete packet algorithm and a practical crossover remain unverified.

#pagebreak()
= Ordinary counts, guarded counts and isolation
For fixed $T$, let $H$ count ordered tuples $(P_1,...,P_m)$ with $P_i in F_i$ and $sum_i P_i=T$. The arithmetic graph excludes intermediate identities and equal-abscissa chord steps. Write $X$ for the number of tuples passing these guards.

== Recovering an ordinary count by rounding
Insert a free point $D_0$ before a tuple and $-D_0$ after it. When $T=O$, insert a fixed nonidentity point before the closing step and shift the terminal target accordingly. Let $Z$ count the guarded extended paths. An ordinary tuple has at most $2m+7$ bad offsets. An interior step with addend $A$ forbids at most the two offsets for which the previous partial point equals $plus.minus A$. The closing step can introduce at most four solutions of a doubling equation. Hence
$ 0 <= N H-Z <= (2m+7)H. $
Fixing all but one slot leaves at most one completing point. If $(2m+7) min_i product_(j!=i) abs(F_j)<N$, then $H=ceil(Z/N)$, including at the identity. In the counting graph the offset is a variable summed internally. The finite reference experiment enumerates it independently to test the rounding.

== A fixed-target second moment
For slot densities $rho_i=tau_i/M$ put $nu_i=N rho_i$, $lambda=(product_i nu_i)/N$, and $kappa=abs(G[2])<=4$ for $G=E(bb(F)_p)$. For $2T != O$,
$ bb(E)[H(H-1)] <= lambda^2 product_i (1+2/nu_i) + lambda kappa (2^m-2)/N. $
To see the first term, group pairs of tuples by the slots where their abscissas coincide. At a matching abscissa there are at most two signs. At a pair of distinct abscissas the joint retention probability is at most $rho_i^2$: an odd multiplier preserves the 2-adic valuation of their difference, and a dyadic window admits at most the uniform fraction at that valuation. For a proper match set the prescribed target forces the last unmatched point. When every abscissa matches, a nonempty proper set of sign flips must sum to a point of order dividing two. Counting these flip sets gives the second term. Flipping every sign would require $2T=O$.

At balanced widths with $nu_i$ large compared with $m$, $2^m/N=o(lambda)$ and $lambda$ of order $1/L$, the chance of one usable relation is $lambda+O(lambda^2)$ and the chance of two or more is $O(lambda^2)$. In particular a fixed-target existence classifier must resolve most positive inputs to attain absolute error $O(L^(-2))$. This statement averages over fresh multipliers and offsets at the given target. It requires no distributional assumption on the target.

#pagebreak()
= What planting reveals about the seeds
Let $P$ denote the fresh seed law. To draw from the planted law $Q$, sample a uniform ordered affine tuple of sum $T$, then independently draw each hash conditional on retaining its marked point. The last point is forced by the first $m-1$, with rejection if it equals $O$. For a given marked abscissa, sample $alpha$ as usual, draw $u$ uniformly below $tau$, and set $beta=u-alpha x mod M$. Counting the ways to plant the same seeds gives
$ (d Q)/(d P) = H/mu, quad mu=bb(E)_P H=((N-1)^m-(-1)^m)/N product_i rho_i. $
The numerator counts affine tuples by convolution of the indicator of the nonidentity elements. For any acceptance event $A$, the score $mu+P(A)-2mu Q(A)$ differs from its fresh existence error by at most $bb(E)_P (H-1)_+$. This permits calibration using planted examples and unlabeled fresh seeds at a prescribed target.

== An unread slot forces first-order error
Let $R$ plant *independent* uniform affine marks in the slots, with no sum condition. If $D=product_i abs(F_i)$, then $(d R)/(d P)=D/(bb(E)_P D)$. Under $R$, conditionally on the hashes, the marks are uniform independently in the observed slots. Enforcing the target sum gives $Q$. For a nonempty observed subset $I$ of size $s$, write $t=m-s$, $D_I=product_(i in I) abs(F_i)$ and $H_I$ for its count at $T$. The exact marginal likelihood is
$ (d Q_I)/(d R_I) =
  (1+(-1)^t (N H_I/D_I-1)/(N-1)^t)/(1-(-1)^m/(N-1)^m). $
Indeed $t$ independent affine marks sum to $u$ with probability
$ frac(1,N) (1+ frac((-1)^t (N bold(1)_(u=O)-1),(N-1)^t)). $
Average this identity over the observed marks. Under $R$, $W=H_I/D_I$ lies in $[0,1]$ with mean $a_s=(1-(-1)^s/(N-1)^s)/N$. The inequality $bb(E) abs(W-a_s)<=2a_s (1-a_s)$ now gives total variation $O((N-1)^(-t))$ between these marginal laws.

An adaptive algorithm reading at most $m-1$ slots may be partitioned by the least index it leaves unread. Its decision on the part corresponding to index $i$ is measurable on all slots except $i$. Summing gives $abs(Q(A)-R(A))=O(m/N)$. Meanwhile $op("Var")(abs(F_i))<=2bb(E) abs(F_i)$ by the two-abscissa bound. Independence gives relative variance of $D$ at most $product_i (1+2/(bb(E) abs(F_i)))-1$, exponentially small in $b/m$ at the balanced widths. Thus $Q(A)<=P(A)+o(lambda)$. As $P(H>0,A)<=mu Q(A)$, fresh error is at least $lambda-O(lambda^2)$. The statement permits arbitrary computation on every slot read. A statistic supplied from all slots, including a single pooled moment, lies outside this access model.

#pagebreak()
= Evaluating the joint moments
Write $d=m+1$. In the space of functions on $E$ with at most a pole of order $d$ at $O$, choose the monomials $x^i$ with $2i<=d$ and $y x^i$ with $2i+3<=d$. Let $v(P)$ be their evaluation vector. The determinant
$ Delta=det[v(P_1),dots,v(P_m),v(-T)] $
vanishes on a tuple of distinct points precisely when their sum is $T$. A vanishing determinant supplies a section through the $d$ points. Its degree-zero divisor is principal exactly when their group sum is $O$, the elliptic-code criterion used by Cheng. Call the masks regular if their abscissas are pairwise disjoint and avoid $x(T)$. Failure has probability at most $N sum_(i<j) rho_i rho_j+sum_i rho_i=o(lambda^2)$ in the growing-arity tuning.

== Interpolation on a common set of curve points
Let $Z_j=sum_"tuples" Delta^j$ be a field element. Fix $k>=1$ and $n=k d$ with $n+1$ at most half the affine population. Choose $n+1$ distinct affine nodes. A nonzero function of pole order at most $n$ cannot vanish at all of them. Thus the evaluation matrix of the $n$ basis monomials has rank $n$. Select $n$ independent rows. For each slot solve for weights $w_(i,a)$ so that its sum functional agrees on those nodes for every function of pole order at most $n$. For each $j<=k$ the integrand $Delta^j$ belongs to that space in each slot separately. Applying the identity successively gives
$ Z_j=sum_(a_1,...,a_m) (product_i w_(i,a_i)) Delta(B_(a_1),...,B_(a_m),-T)^j. $
The cost, up to polynomial factors in $b,m,log(k+1)$, is
$ n sum_i tau_i + n^3 + m n^2 + k n^m . $
The first term includes lifting slots and accumulating monomials, the next two solve for the weights, and the last evaluates all moments on the common grid. For $k=floor(2^(b/m^2))$ at the proposed growing arity, the logarithm of the expected bit cost is $(1+o(1))b/m$. The moment vector therefore fits the proposed budget. Decoding the multiplicity of zero from that vector remains the additional step.

== A restriction on reconstruction
On the regular prime-order family, collision bounds on the determinant values give a large typical spectrum at the proposed widths. A generic recurrence that first represents all distinct determinant values can require $2^(b-o(b))$ terms. Independently, low-pole-order zero tests whose coefficients are fixed before choosing the target can change a constant decision only on as many targets as their total pole order. These statements have specific access models. They allow the possibility of a nonlinear, fixed-target decoder trained on all joint values.

#pagebreak()
= A zero-count certificate for two sign orbits
Suppose $p>7$, $m=3$, and the regular slots consist of one, one and two abscissa pairs in any order. Fixing an abscissa in each slot leaves eight sign choices. The basis $1,x,x^2,y$ uses the $y$ row once, so each orbit has determinant $c_sigma+sum_(i=1)^3 a_(sigma i) epsilon_i$, with $epsilon_i=plus.minus 1$. There are two equally sized orbits and the centers are unobserved.

Write $c=Z_1/16$ and $c_sigma=c+sigma h$. Compute the normalized pooled central moments $M_j$ about $c$. Starting with $A_0=1$, $B_0=0$, the two local even central moments at order $2q$ are $A_q plus.minus B_q/h$ for $q=1,2,3$ whenever $h!=0$. Binomial expansion gives
$ A_q=M_(2q)-sum_(r=0)^(q-1) binom(2q,2r) h^(2(q-r)) A_r, $
$ (2q+1) B_q=M_(2q+1)-sum_(r=0)^(q-1) binom(2q+1,2r) h^(2(q-r)) B_r . $
Scale the centered sign values by $h$ to clear denominators. Denote the local means of their second, fourth and sixth powers by $u_2,u_4,u_6$. The sums $s_j$ of the $2j$th powers of the three scaled sign coefficients satisfy
$ s_1=u_2, quad s_2=(3s_1^2-u_4)/2, quad s_3=(u_6-15u_4s_1+30s_1^3)/16. $
Newton's identities give the elementary symmetric functions $e_2,e_3$ of their squares. For each orbit, its centered scaled sign values satisfy the octic
$ K_sigma (z)=((z^2-s_1)^2-4e_2)^2-64z^2e_3 . $
It follows either by multiplying the eight linear sign factors or by pairing opposite signs. The octic predicts higher moments by a linear recurrence. Equating the predicted pooled eighth and ninth moments to the observations gives polynomial residuals in $h$. Intersect their common factors with $h^(p-1)-1$ by modular powering. The result contains every possible nonzero separation. If the observed odd centered moments all vanish, we allow equal centers and abstain.

At a surviving separation the zero multiplicity is the order at $z=0$ of $K_+(z-h c-h^2)K_-(z-h c+h^2)$. Only coefficients of orders zero through four are needed. Each of the four ordered pairs in the two smaller slots admits at most one completing point. Successive gcds with these five coefficients record the multiplicities compatible with the moments. Agreement certifies the count, and agreement on positivity certifies existence. Spurious separations can cause abstention. For fixed prefix length the readout uses bounded-degree arithmetic and $O(log p)$ field operations. It neither visits the 16 determinants nor reconstructs their spectrum at inference time.

#pagebreak()
= What the finite experiments measure
The small-curve reference enumerates masks and computes determinants independently of the decoder. On $y^2=x^3+4x+1$ over $bb(F)_13$, $N=19$, $T$ is the first affine point in lexicographic order. At widths $(2,2,2)$ modulo 16, the regular inputs with slot sizes $(2,2,4)$ up to permutation account for 168,480 seed triples with multiplicity. Exactly 133,008 have a hit.

#grid((0.9fr,1.15fr,1.15fr,1.2fr), [*Moments*], [*Count certified*], [*Bit certified*], [*Misses on default*],
 [8], [123,312], [123,792], [36,096], [9], [145,776], [146,256], [13,632],
 [10], [151,008], [151,008], [10,272], [12], [151,008], [151,008], [10,272])
#text(size: 9.5pt)[Table 2. The default reports absence on unresolved inputs. Every issued certificate agreed with direct counting. The denominator conditions on this slot shape.]

Keep the nine-moment certificates and bin the unresolved cases by their candidate count sets. Use 60,000 independently generated plants and 60,000 unlabeled fresh seeds to fit the test $2mu Q(v)>P(v)$ for each bin $v$. Across the entire eligible population, existence mistakes fall from 13,632 to 7,680, including 7,200 false positives. For integer counts, on training plants only, check the four possible completions and use $P(H=j,v)=(mu/j)Q(H=j,v)$ for $j>0$. Pick the most probable admissible count. This makes 10,320 mistakes, versus 14,112 for the smallest compatible count, using zero on abstention.

== Independent fresh-seed tests at larger widths
We drew 12,000 plants per field to fit equations, then 16,000 independent fresh seeds per field to score the resulting rules. One family contained the 11 weighted monomials of degree at most eight in three normalized even central moments. Another contained every monomial of ordinary degree at most two in the nine raw field moments, including the constant. On regular masks containing at least three orbits, the planted evaluation matrices attained full rank 11 and 55 in their respective eligible populations on both curves. Every nonzero equation in either family therefore fails on a genuine plant.

#grid((0.85fr,1.25fr,1.1fr,1.2fr,1.1fr), [*Field*], [*Fresh hits*], [*Two-orbit counts certified*], [*Larger-pool hits missed*], [*All errors*],
 [101], [3,669], [1,715 / 1,732], [2,299 / 2,299], [3,363],
 [1,031], [3,274], [44 / 44], [2,929 / 2,929], [3,272])
#text(size: 9.5pt)[Table 3. Widths were $(4,4,4)$ and $(8,16,16)$, with orders 103 and 1,019. The fitted rule returned absence where it had no test. Both trials had zero false positives.]

The larger-field rule captured two of 3,274 fresh positives. These experiments identify a failure of the tested families on typical larger pools. Two finite fields do not estimate an asymptotic error rate for unrestricted decoders.

#pagebreak()
= Other approaches and the open exponent
*Generic packet span.* For a complete generic binary packet of width $w=2^a>=2$, source assignments of Hamming weight one or two and arbitrary targets span the entire rational polynomial space of degree at most $w(2w-1)$. The proof constructs a spanning family of marginals one degree layer at a time. This closes the proposal to exploit a smaller linear span for every such complete packet. Restrictions of a particular arithmetic graph may have a different span. Earlier degree-aware Hermite readouts did save precision on finite fixtures, including losses of 27 instead of 45 bits and 102 instead of 128.

*Separated approximation.* In the odd prime-order case $r=N$, split the independent slots into two nonempty groups and approximate $H-D/r$ by sums of $q$ products of one real function per side. For $m/r+sum_i rho_i=o(1)$, set $n=(r-1)/2$ and suppose $q/n -> x$ in $[0,1]$. The optimal residual mean-square error divided by $lambda$ tends to $1-x-sin(pi x)/pi$. The formula comes from the tail of the squared cosine spectrum after controlling the Gram matrices of the actual hash masks. Rank $o(r)$ preserves first-order residual error in this model.

*Digits, resultants and field units.* Repeating a slot collapses powers in characteristic two but also shrinks the support of detectable targets. The second base-prime digit has an endpoint formula for runs along a common group step. Independent abscissa windows do not supply sufficiently short common runs at the tuning. A scalar norm and its derivatives detect simple zeros. In the paired, odd prime-order setting, a balanced two-block evaluation of the single-target unit and its moments takes $2^((1/2+o(1))b)$ time and space. A three-block evaluation of the normalized scalar resultant and derivatives trades a time exponent of $2/3$ for a space exponent of $1/3$. Normalized section products admit a group convolution with two integer moments. The exceptional field-unit correction can even be projected into one integer count given a primitive field element, but the charged circuit is too large for a better exponent. These are costs for particular representations. Detailed hypotheses appear in the supplement.

*The missing algorithm.* At $m asymp L^(1/4)sqrt(g)$, $g=log_2 L$, a fresh predictor with absolute error $O(L^(-2))$ at each eligible target and log bit cost $O(b/m)$, with setup charged, would give an outer exponent of order $b/(L^(1/4)sqrt(g))$. The fixed-width reduction rerandomizes a slot to read a dyadic branch and lists the sparse set of possibly wrong branch decisions. At $m asymp sqrt(b)$, an exact restricted-window parity oracle of cost $2^(tilde(O)(sqrt(b)))$, uniform over the adaptive transcript, would reach the square-root-in-bitlength objective under the outer reduction. Both evaluator hypotheses remain open. In the first regime a typical seed has on the order of $N lambda/2^m=2^(b-o(b))$ sign orbits. The two-orbit certificate gives no inexpensive way to combine them.

#pagebreak()
= Provenance and reproducibility
*Public baseline.* Akshima, Renas Bacho, Yu-Hsuan Huang, Eike Kiltz, Aysan Nishaburi, Samin Nooripoor, Noah Tittelbach and Emiel Wiedijk. _A Subexponential Time Algorithm for Elliptic Curve Discrete Logarithm_, version 0.2, 29 October 2026, especially Theorems 1, 2 and 4. #link(public)[Authors’ manuscript at the pinned revision.]

*Elliptic evaluation codes.* Qi Cheng. _Hard Problems of Algebraic Geometry Codes_, 2005. The divisor criterion on an elliptic curve supplies the geometric input for the determinant. #link("https://arxiv.org/abs/cs/0507026")[arXiv:cs/0507026.]

*Supplied working manuscripts.* _Simplex Alphabets and Backward Secants for Elliptic-Curve Discrete Logarithms_, v6.0, 31 October 2026, and _Compressed Relation Search and Reusable Preprocessing for Elliptic-Curve Discrete Logarithms_, 31 October 2026, especially Assumption A.1. These are separate sources supplied to the project. The latter lists no named author. Their numerical baselines and finite-packet premises are attributed to those documents. Archive hashes are recorded in the repository README. The authors of the public baseline are not attributed these later constants.

*Supplement.* The repository contains the project’s mathematical derivations for the inverse-hash accounting, circuit joins, target isolation, packet span, determinant moments, planted certificates, and the norm, section-product and field-unit alternatives. The topic index maps them to the sections of this paper. These are working proofs, including scope conditions and failed approaches, available for independent review.

*Reproduction.* The release includes four command-line checks and their recorded output. The accounting check recomputes the coefficients, 8,192 small inverse-hash cases, 160 larger products, rounding at 103 targets, and 24 comparisons of interpolated to enumerated moments. A second script exhausts the small pooled inputs and checks the marginal planting law on 756 partial mask patterns. A third runs the larger fresh-seed trials using fixed pseudorandom seeds and separate training and evaluation streams. A fourth verifies finite modular ranks behind the generic packet-span theorem. The small-curve reference arithmetic is included. No source archive is required to run these checks. #link(project)[Repository, paper source, derivations, and experimental evidence.]

*Status.* This is an AI-assisted working paper prepared for public scrutiny. The algebraic propositions have direct arguments and finite checks under their stated hypotheses. The end-to-end numerical bound remains conditional on the named upstream ingredients. No full-scale timing or independent peer review is claimed.
