#set document(title: "Counting relations in independently hashed elliptic-curve slots", author: "davidiach research repository", date: datetime(year: 2026, month: 11, day: 22))
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
#let note(path, body) = link(project + "/blob/8171080efb0d9024f49eb47e1e311914f0f7a4c8/derivations/" + path, body)

#text(size: 9pt, tracking: 0.8pt, fill: rgb("#31584f"))[RESEARCH PREPRINT · 22 NOVEMBER 2026]
#v(8pt)
#block[#set par(justify: false); #text(size: 26pt, weight: "bold")[Counting relations in independently hashed elliptic-curve slots]]
#v(4pt)
#text(size: 14pt)[Conditional complexity, planted distributions, and certificates from determinant moments]
#v(4pt)
Published by #link(project)[davidiach] · AI-assisted mathematical research · Version 3
#v(9pt)
#label[Abstract.] We count ordered representations of a fixed elliptic-curve point $T$ as a sum from independently sampled slots. With affine hashes at the proposed sparse widths and $2T != O$, the probability of exactly one representation is $lambda+O(lambda^2)$, where $lambda$ is the product of the nominal slot sizes divided by the group order. Every adaptive existence predictor that leaves a whole slot unread has error at least $lambda-O(lambda^2)$. Joint determinant moments recover exact counts on finite examples. For regular slots of sizes 2, 2 and 4 over $bb(F)_13$, nine moments certify 145,776 of 168,480 counts and the twelfth moment determines them all. Weighted measurements also decide existence on dense three-slot inputs. Efficient decoding at growing arity remains open. Under the supplied finite-packet assumption and inherited reduction estimates, the leading time coefficient for general field primes is 26.29338 with affine hashes throughout, or 24.00246 when the algorithm may choose a mixed sampler. The latter coefficient falls to 19.59793 for primes with bounded signed radix support and to 13.85783 for primes $p=2^b-c$ with positive, polynomially bounded $c$.

#v(5pt)
= The question and its baselines
Let $b=ceil(log_2 p)$ and $L=log_2 b$. Akshima, Bacho, Huang, Kiltz, Nishaburi, Nooripoor, Tittelbach and Wiedijk derive a conditional expected-time bound $2^(O(b/L^(1/4)))$ on prime-order curves over prime fields, with an extension to fixed cofactors. Their version 0.2 assumes a Boolean tensor-contraction result. It gives no numerical leading coefficient for this comparison. A separate, supplied manuscript develops an explicit finite-packet model and reports a coefficient of 27.01388. We compare numerical coefficients only within that framework on the overlapping input class.

Write $Z_b=b/L^(1/4)$. The supplied simplex and compressed-relation manuscripts give time/space coefficient pairs of 32.15695/16.07848 and 27.01388/13.50694, respectively. Table 1 gives the bounds derived here.

#grid((1.8fr,2.5fr,0.9fr,0.9fr), [*Field prime $p$*], [*Sampler and arithmetic*], [*Time*], [*Space*],
  [General], [All affine, short inverse hash], [26.29338], [13.14669],
  [General], [Mixed, fixed abscissa bits], [24.00246], [12.00123],
  [Bounded signed radix support], [Mixed, local quotient digits], [19.59793], [9.79896],
  [$p=2^b-c$, $1<=c<=b^C$], [Mixed, weighted cyclic transform], [13.85783], [6.92891])
#text(size: 9.5pt)[Table 1. Entries are rounded leading coefficients. For an exact coefficient $a$ from Section 2, the upper bound is $(a+o(1))Z_b$ on base-two log expected bit time or log space. “General” allows every field prime in the inherited curve and group-order class of Section 2. The last two rows require the indicated modulus structure, with constants independent of $b$.]

#pagebreak()
= Conditional bounds by sampler and field prime
The conditional comparison uses a smooth short Weierstrass curve $E$ over $bb(F)_p$, $p>3$, with supplied group order $N=h r$, odd prime subgroup order $r$, $gcd(h,r)=1$, and polynomially bounded cofactor $h$. The discrete logarithm inputs are a generator of that subgroup and a point in it. The prime restrictions in Table 1 apply to the field characteristic $p$. A slot accepts an affine point when its canonical abscissa $x$ satisfies $(alpha x+beta) mod M < tau$, where $M=2^b$, $alpha$ is uniform odd, $beta$ is uniform modulo $M$, and $tau=2^ell$. Different slots use independent seeds. Both signs of a point are retained.

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
== General primes with a choice of sampler
Keep the odd affine hash in slot 1. In each other slot choose $ell$ public free bit positions and fill the remaining positions independently with uniform bits. Each legal point survives with probability $rho=w/2^b$, where $w=2^ell$. This changes the seed law. Sections 3 and 4 use the original all-affine law for their isolation and unread-slot bounds.

To bound the target coverage, write $S$ for the words supported on the free bits of one slot. If those positions form $s$ runs, then $abs(S-S)<=2^s w$. For $D!=O$ the function $x(P+D)-x(P)$ has degree four, so each admissible abscissa difference has at most four preimages. In a fixed-bit slot the sum over $P$ of the joint retention probabilities for $P,P+D$ is therefore at most $4rho 2^s w$. In slot 1 two distinct abscissas survive together with probability at most $rho^2$. Applying these estimates to pairs of tuples and using Cauchy–Schwarz gives inverse success probability $2^(O(m))$ for a guarded relation at each nonidentity target, at the paid widths. Here $m$ is the number of slots and $w$ is the smallest dyadic width with $N rho >= max(2N^(1/m),20m)$ and $m asymp L^(1/4)$. These retries add only $O(m)$ to log time.

Choose the free bits in order of the Fourier rows that own them. Only an $O(1/m+1/V)$ fraction of the $V$ rows then contain a free bit, and the number of runs satisfies $m 2^s/w=o(1)$. The other rows have public input symbols. Substituting them and removing the redundant hash outputs preserves each accepted assignment with its multiplicity. The graph carries $(30+o(1))m b$ leading word bits. On two finite compiler fixtures this specialization reduced the edge counts from 4,688 to 3,988 and from 8,122 to 7,202. #note("strided_slots/results.md")[Coverage proof and compiler accounting.]

== Why removing quotient inputs ties the general bound
For a canonical quadratic residual $D$ with $abs(D)<8p^2$, put $E=D+8p^2$. In radix $B_0=2^k$, choose a padded length $n$ and a sign so that $R=B_0^n plus.minus 1>16p^2$ is coprime to $p$. Let $z$ be the least nonnegative residue of $p^(-1) E$ modulo $R$. If $p$ divides $D$, then $z=E/p<16p$. Conversely, when $z<16p$, the congruent integers $p z$ and $E$ both lie in $[0,R)$ and hence are equal. A public multiplication and a cyclic carry test implement this comparison. They eliminate the transmitted quotient. Recovering the larger product coefficients widens a word from $2k+O(log b)$ to $3k+O(log b)$, bringing the total back to $(30+o(1))m b$. The coefficient for general primes stays at 24.00246. #note("group_equations/results.md")[Residual bounds and finite comparisons.]

#pagebreak()
== Two prime families with lower coefficients
*Bounded signed radix support.* Use the compiler radix $B_0=2^k$, with $k=Theta((log b)^2)$. Suppose $p=sum_(d in S) p_d B_0^d$, where $0 in S$, $abs(p_d)<B_0$, $abs(S)=O(1)$, and $max S<=ceil(b/k)+O(1)$. For a quadratic residual shifted to a nonnegative integer $E$, write $E_i$ for its raw radix coefficients. At position $i$, the quotient digits $q_i$ and carries $a_i$ obey
$ E_i+a_i-sum_(d in S) p_d q_(i-d)=B_0 a_(i+1). $
The constant digit $p_0$ is odd. Given the incoming carry and earlier digits, reduction modulo $B_0$ determines $q_i$ uniquely. Initialization, terminal carry zero, and the quotient range enforce the integer equality $E=p q$. Bounded support lets the processor groups exchange the few shared digits using $o(b)$ bits per equation. The leading graph budget falls to $(20+o(1))m b$. #note("group_equations/sparse.md")[Digit and routing proof.]

*Primes just below a power of two.* Assume $p=2^b-c$ with $1<=c<=b^C$ for a fixed $C$. The identity $2^b equiv c mod p$ folds the high part of a quadratic product into its low part. Using place exponents $ceil(i b/n)$ lets a length-$n$ weighted cyclic transform do this even when $n$ does not divide $b$. The folded coefficients fit in $2k+O(log b)$ bits since $c$ is polynomially bounded. Pairing conjugate frequencies accounts for both components of each extension value. The total per array is $(2+o(1))b$ bits and the graph budget is $(10+o(1))m b$. More generally the same argument applies when a known odd positive $v$ and positive $c$ are polynomially bounded in $b$ and $v p=2^s-c$ for $s=b+O(log b)$: divisibility of $v D$ by $v p$ tests divisibility of $D$ by $p$. #note("group_equations/weighted.md")[Weighted transform, carries and small multiples.]

*From word budgets to time and space.* If the leading word budget is $(a+o(1))m b$, put $u=a K/2$. The inherited contraction, enumeration and sparse solve give the leading log-time terms $b/m+u m b/sqrt(L)$ and $2b/m$. Balancing at $m tilde.eq L^(1/4)/sqrt(u)$ gives coefficients $2sqrt(u)$ for time and $sqrt(u)$ for space. Thus the mixed-sampler budgets $a=30,20,10$ give the last three rows of Table 1. All retain the same $Z_b$ scale. Even on the weighted family, the coordinates and slopes alone cost $Theta(m b)$ in this graph. Taking $m/L^(1/4) -> infinity$ would reduce the outer term $2b/m$ to $o(Z_b)$ but would raise the inherited contraction estimate above $Z_b$.

*Limits of two further simplifications.* In the alternating chord chain with $e=2m-2+floor((m-1)/2)$ residuals, checking at most $e-3$ independent uniform linear combinations admits a false assignment with probability tending to one on the negative-input family of the batching proof when $2^m/p -> 0$. Exact algebraic elimination faces a different constraint: over the algebraic closure in characteristic above three and with $m>=2$ and $2T!=O$, the generic fixed-target abscissa hypersurface has total degree $m 2^(m-1)$. Bézout’s bound forces the product of the defining degrees of any exact lift on a dense open to be at least this large. In particular a quadratic lift needs at least $m-1+ceil(log_2 m)$ scalar equations on that open set. #note("group_equations/batching.md")[Batching argument.] #note("group_equations/chord_geometry.md")[Degree argument.]

#pagebreak()
= Ordinary counts, guarded counts and isolation
Return to independent affine hashes in every slot. For fixed $T$, let $H$ count ordered tuples $(P_1,...,P_m)$ with $P_i in F_i$ and $sum_i P_i=T$. The arithmetic graph excludes intermediate identities and equal-abscissa chord steps. Write $X$ for the number of tuples passing these guards.

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
Fix a nonidentity target $T$ for the likelihood formulas below. For the sparse-error bound we also assume $2T != O$, so that the fixed-target second-moment bound of Section 3 applies.

Let $P$ be the fresh law. For the planted law $Q$, draw a uniform ordered affine tuple of sum $T$ by sampling the first $m-1$ points and rejecting an identity completion. At each marked abscissa $x$, draw $alpha$ as usual and $u$ uniformly below $tau$, then set $beta=u-alpha x mod M$. Counting the possible plants gives
$ (d Q)/(d P) = H/mu, quad mu=bb(E)_P H=((N-1)^m-(-1)^m)/N product_i rho_i. $
Convolving the indicator of the nonidentity elements gives $((N-1)^m-(-1)^m)/N$ ordered affine tuples at each nonidentity target. For any acceptance event $A$, the score $mu+P(A)-2mu Q(A)$ differs from its fresh existence error by at most $bb(E)_P (H-1)_+$. This permits calibration using planted examples and unlabeled fresh seeds at a prescribed target.

== An unread slot forces first-order error
Let $R$ plant *independent* uniform affine marks in the slots, with no sum condition. If $D=product_i abs(F_i)$, then $(d R)/(d P)=D/(bb(E)_P D)$. Under $R$, conditionally on the hashes, the marks are uniform independently in the observed slots. Enforcing the target sum gives $Q$. For a nonempty observed subset $I$ of size $s$, write $t=m-s$, $D_I=product_(i in I) abs(F_i)$ and $H_I$ for its count at $T$. Then
$ (d Q_I)/(d R_I) =
  (1+(-1)^t (N H_I/D_I-1)/(N-1)^t)/(1-(-1)^m/(N-1)^m). $
Indeed $t$ independent affine marks sum to $u$ with probability
$ frac(1,N) (1+ frac((-1)^t (N bold(1)_(u=O)-1),(N-1)^t)). $
Average this identity over the observed marks. Under $R$, $W=H_I/D_I$ lies in $[0,1]$ with mean $a_s=(1-(-1)^s/(N-1)^s)/N$. The inequality $bb(E) abs(W-a_s)<=2a_s (1-a_s)$ now gives total variation $O((N-1)^(-t))$ between these marginal laws.

An adaptive algorithm reading at most $m-1$ slots may be partitioned by the least index it leaves unread. Its decision on the part corresponding to index $i$ is measurable on all slots except $i$. Summing gives $abs(Q(A)-R(A))=O(m/N)$. Meanwhile $op("Var")(abs(F_i))<=2bb(E) abs(F_i)$ by the two-abscissa bound. Independence gives relative variance of $D$ at most $product_i (1+2/(bb(E) abs(F_i)))-1$, exponentially small in $b/m$ at the balanced widths. Thus $Q(A)<=P(A)+o(lambda)$. Since $P(H>0,A)<=mu Q(A)$, fresh error is at least $P(H>0)+(1-2mu)P(A)-o(lambda^2)$. At the proposed widths $mu$ tends to zero, so the coefficient of $P(A)$ is eventually nonnegative and the lower bound is $lambda-O(lambda^2)$. This access bound permits arbitrary computation on the observed slots.

== A sufficient weak-detection criterion
At the balanced all-affine widths, let a score $f in [0,1]$ cost $D_m$ and have a known targetwise gap $abs(bb(E)_Q f-bb(E)_R f)>=delta$ with $m 2^m/N=o(delta)$. Assume the second-moment bound and $bb(E)(H-X)/mu=O(m/N)$. Then resampling slots and verifying completions gives inverse success polynomial in $lambda^(-1),2^m,delta^(-1)$ and cost $(1+sum_i tau_i)(1+D_m) op("poly")(b,m,lambda^(-1),2^m,delta^(-1))$. Under the inherited reduction, $m asymp sqrt(b)$ and $D_m,delta^(-1)<=2^(O(sqrt(b)))$ would suffice for the square-root-in-bitlength goal. Section 5 tests a candidate score. #note("strided_slots/weak_detection.md")[Recovery proof.]

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

== Testing a joint moment under planting
Section 4 compares planted seeds with independently marked seeds because both laws favor larger slots. Consider the score that tests whether $Z_2(T)$ is a nonzero square. Over $bb(F)_65521$, one million draws per law gave planted-minus-matched gaps of $+0.00090$ at four slots and $-0.00005$ at six. Both individual 95% intervals, of half-width $0.00192$, contain zero. A uniform advantage of the size needed for recovery remains open. #note("strided_slots/weak_detection.md")[Finite tests and score definition.]

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

#grid((1.5fr,1fr,1fr,1.2fr), [*Readout and moments*], [*Count certified*], [*Bit certified*], [*Misses on default*],
 [Octic, 1–8], [123,312], [123,792], [36,096], [Octic, 1–9], [145,776], [146,256], [13,632],
 [Octic, 1–10], [151,008], [151,008], [10,272], [Octic, 1–12], [151,008], [151,008], [10,272],
 [Fermat, 12 alone], [168,480], [168,480], [0])
#text(size: 9.5pt)[Table 2. The default reports absence on unresolved inputs. In $bb(F)_13$, each nonzero determinant has twelfth power one, so the raw moment satisfies $Z_12 equiv 16-H mod 13$. Each of the four pairs from the smaller slots admits at most one completion. Thus $0 <= H <= 4$ and the residue determines the integer count.]

Keep the nine-moment certificates and bin the unresolved cases by their candidate count sets. Use 60,000 independently generated plants and 60,000 unlabeled fresh seeds to fit the test $2mu Q(v)>P(v)$ for each bin $v$. Across the entire eligible population, existence mistakes fall from 13,632 to 7,680, including 7,200 false positives. For integer counts, on training plants only, check the four possible completions and use $P(H=j,v)=(mu/j)Q(H=j,v)$ for $j>0$. Pick the most probable admissible count. This makes 10,320 mistakes, versus 14,112 for the smallest compatible count, using zero on abstention.

== Independent fresh-seed tests at larger widths
We drew 12,000 plants per field to fit equations, then 16,000 independent fresh seeds per field to score the resulting rules. One family contained the 11 weighted monomials of degree at most eight in three normalized even central moments. Another contained every monomial of ordinary degree at most two in the nine raw field moments, including the constant. On regular masks containing at least three orbits, the planted evaluation matrices attained full rank 11 and 55 in their respective eligible populations on both curves. Every nonzero equation in either family therefore fails on a genuine plant.

#grid((auto,1fr,1fr,1fr,1fr,1fr,1fr), [*Field*], [*Fresh hits*], [*Irregular*], [*1 orbit*], [*2 orbits*], [*3+ orbits*], [*Misses*],
 [101], [3,669], [1,064], [45], [261], [2,299], [3,363],
 [1,031], [3,274], [343], [0], [2], [2,929], [3,272])
#text(size: 9.5pt)[Table 3. Fresh hits are partitioned by mask regularity and, for regular masks, by orbit count. All misses came from irregular masks or regular masks with at least three orbits. Both trials had zero false positives. Widths were $(4,4,4)$ and $(8,16,16)$, with group orders 103 and 1,019, respectively.]

Across all regular two-orbit seeds, including negatives, nine moments certified the exact count on 1,715 of 1,732 seeds over $bb(F)_101$ and all 44 over $bb(F)_1031$. On the latter field the fitted rule detected two of 3,274 fresh positives. The behavior of unrestricted decoders at growing field size remains open.

#pagebreak()
== What target derivatives and weighted measurements recover
*Varying the target.* Suppose there are three slots, $p>5$, and every abscissa occurs with both nonzero ordinates. Write $f(X)=X^3+a X+c$ and view the second determinant moment as a function of the target abscissa $X$. Expanding in the ordinate column and summing over the eight sign choices cancels every mixed term. The resulting $Q(X)=Z_2(X)$ has degree at most four. Its five Taylor coefficients capture its values at every target abscissa. This permits exact finite tests of how much information the moment carries. It also exposes a loss of information. Two regular seed triples over $bb(F)_37$ have identical unweighted slot Gram matrices and all target derivatives, yet at $T=(30,5)$ one has a unique witness and the other has none. #note("joint_parity/target_jets.md")[Sign identity and Taylor-coefficient tests.] #note("joint_parity/quartic_decoder.md")[Collision data.]

*A certificate from six field values.* Suppose two slots each have one abscissa, $u$ and $v$, and the third has $1<=k<=4$ distinct abscissas. Assume regularity and full sign pairs. Along with the five coefficients of $Q$, retain $G$, the Gram coefficient of $y(T)^2$ before the curve equation is substituted. Both $u$ and $v$ are roots of $Q(X)-2G f(X)$. For each candidate pair, matching coefficients gives equations for the first four power sums of the third slot. Newton’s identities give candidate support polynomials. A signed pair-sum test decides existence for each surviving candidate, and agreement certifies the bit. An unresolved coefficient system or disagreement among surviving candidates leads to abstention.

In 20,000 fresh draws over $bb(F)_101$ at width four, 3,284 positives had eligible masks absent from training. Of these, 599 had the supported shape. The algebraic certificate found all 599 with zero false positives on unseen masks. A lookup table on the same six values found zero unseen positives. On the $bb(F)_13$ population of Table 2 the six-value rule certifies the existence bit on all 168,480 seeds. The twelfth-moment rule there recovers the integer count. #note("joint_parity/quartic_decoder.md")[Holdout and exhaustive populations.]

*Dense slots.* Assume $p>5$, an odd prime-order rational group, full sign pairs and an affine target. Let $A_i$ be the retained abscissas and $k_i=abs(A_i)$. Weighting one slot by $x^j$ and tagging the other two yields sums $sum_(x in A_i) h_i(x)x^j$ for $0<=j<2k_i$. A recurrence recovers the polynomial with roots at the visible abscissas. Independent retagging recovers every root in a constant expected number of rounds on nonempty disjoint slots. Translating the three support polynomials then decides existence. Overlapping slots can be completed with two additional sums per slot and an original membership test for repeated singletons. For $k_i<=n$, the expected bit cost after enumerating the slots is $tilde(O)(n^2 op("poly")(log p))$. In 3,000 fresh draws over $bb(F)_65521$, weighted measurements decided all 2,935 regular dense cases, including all 666 positives. Directly enumerating pairs has the same quadratic exponent at three slots. #note("joint_parity/dense_weighted.md")[Reconstruction, overlaps and cost.]

*What the root test detects at larger arity.* With $m>=2$ nonempty sign-paired slots disjoint in abscissa and an odd prime-order group, give each signed point its own formal weight. Sum the products of the two determinants at $T$ and $-T$, each multiplied by its tuple’s weight. This polynomial in the weights vanishes identically precisely when one slot forces $x(T)$ as its only abscissa. At the balanced widths of Section 3, the probability that even one slot has at most eight points is $O(m/min_i nu_i)=o(lambda^2)$. A rule restricted to these shapes misses first-order positive mass. #note("joint_parity/quartic_decoder.md")[Opposite targets and the sparse-shape bound.]

#pagebreak()
= Other approaches and the open exponent
*Generic packet span.* For a complete generic binary packet of width $w=2^a>=2$, source assignments of Hamming weight one or two and arbitrary targets span the entire rational polynomial space of degree at most $w(2w-1)$. The proof constructs a spanning family of marginals one degree layer at a time. This closes the proposal to exploit a smaller linear span for every such complete packet. Restrictions of a particular arithmetic graph may have a different span. Earlier degree-aware Hermite readouts did save precision on finite fixtures, including losses of 27 instead of 45 bits and 102 instead of 128.

*Separated approximation.* For a nonidentity target in the odd prime-order case $r=N$, split the independent slots into two nonempty groups and approximate $H-D/r$ by sums of $q$ products of one real function per side. For $m/r+sum_i rho_i=o(1)$, set $n=(r-1)/2$ and suppose $q/n -> x$ in $[0,1]$. The optimal residual mean-square error divided by $lambda$ tends to $1-x-sin(pi x)/pi$. The formula comes from the tail of the squared cosine spectrum after controlling the Gram matrices of the actual hash masks. Rank $o(r)$ preserves first-order residual error in this model.

The supplement extends this calculation to cofactor groups when the target has odd prime order $r -> infinity$. Removing the at most four real-character terms leaves the same limiting cosine spectrum at balanced widths. Any real rank $q=o(N)$ across the chosen split leaves mean-square error at least $(1-o(1))lambda$ for the existence bit. If a score in $[-1,1]$ of rank $2^(O(b/m))$ gives the correct sign with margin $gamma$ outside fresh probability $O(lambda^2)$, polynomial approximation to its threshold forces $gamma=O(log(1/lambda)/m)$. This limits classifiers with a fixed positive margin. #note("joint_parity/dense_weighted.md")[Cofactor and margin proofs.]

*Digits, resultants and field units.* Repeating a slot collapses powers in characteristic two but also shrinks the support of detectable targets. The second base-prime digit has an endpoint formula for runs along a common group step. Independent abscissa windows do not supply sufficiently short common runs at the tuning. A scalar norm and its derivatives detect simple zeros. In the paired, odd prime-order setting, a balanced two-block evaluation of the single-target unit and its moments takes $2^((1/2+o(1))b)$ time and space. A three-block evaluation of the normalized scalar resultant and derivatives trades a time exponent of $2/3$ for a space exponent of $1/3$. Normalized section products admit a group convolution with two integer moments. The exceptional field-unit correction can even be projected into one integer count given a primitive field element, but the charged circuit is too large for a better exponent. These are costs for particular representations. Detailed hypotheses appear in the supplement.

*The missing algorithm.* At $m asymp L^(1/4)sqrt(g)$, $g=log_2 L$, a fresh predictor with absolute error $O(L^(-2))$ at each eligible target and log bit cost $O(b/m)$, with setup charged, would give an outer exponent of order $b/(L^(1/4)sqrt(g))$. The fixed-width reduction rerandomizes a slot to read a dyadic branch and lists the sparse set of possibly wrong branch decisions. At $m asymp sqrt(b)$, an exact restricted-window parity oracle of cost $2^(tilde(O)(sqrt(b)))$, uniform over the adaptive transcript, would reach the square-root-in-bitlength objective under the outer reduction. Both evaluator hypotheses remain open. In the first regime a typical seed has on the order of $N lambda/2^m=2^(b-o(b))$ sign orbits. The two-orbit certificate gives no inexpensive way to combine them.

#pagebreak()
= Provenance and reproducibility
*Public baseline.* Akshima, Renas Bacho, Yu-Hsuan Huang, Eike Kiltz, Aysan Nishaburi, Samin Nooripoor, Noah Tittelbach and Emiel Wiedijk. _A Subexponential Time Algorithm for Elliptic Curve Discrete Logarithm_, version 0.2, 29 October 2026, especially Theorems 1, 2 and 4. #link(public)[Authors’ manuscript at the pinned revision.]

*Elliptic evaluation codes.* Qi Cheng. _Hard Problems of Algebraic Geometry Codes_, 2005. The divisor criterion on an elliptic curve supplies the geometric input for the determinant. #link("https://arxiv.org/abs/cs/0507026")[arXiv:cs/0507026.]

*Supplied working manuscripts.* _Simplex Alphabets and Backward Secants for Elliptic-Curve Discrete Logarithms_, v6.0, 31 October 2026, and _Compressed Relation Search and Reusable Preprocessing for Elliptic-Curve Discrete Logarithms_, 31 October 2026, especially Assumption A.1. These are separate sources supplied to the project. The latter lists no named author. Their numerical baselines and finite-packet premises are attributed to those documents. Archive hashes are recorded in the repository README. The authors of the public baseline are not attributed these later constants.

*Supplement.* The repository contains the project’s mathematical derivations for the inverse-hash accounting, circuit joins, target isolation, packet span, determinant moments, planted certificates, and the norm, section-product and field-unit alternatives. The topic index maps them to the sections of this paper. These are working proofs, including scope conditions and failed approaches, available for independent review.

*Reproduction.* The release contains 52 standalone checks and four C38 integration checks, each with a recorded transcript. Run `python reproduce.py` for the standalone checks. With the two original source ZIPs available locally, run `python reproduce.py --source-archives DIR` for all 56. The integration checks exercise the fixed-bit compiler, the general-prime reciprocal, sparse quotient digits and weighted transforms on finite inputs. Other checks compare the counting and decoding results with direct elliptic addition or determinant enumeration. The 29 historical replay entries have separate inputs and instructions. #link(project)[Reproduction instructions and recorded evidence.]

*Status.* This is an AI-assisted working paper prepared for public scrutiny. The algebraic propositions have direct arguments and finite checks under their stated hypotheses. The end-to-end numerical bound remains conditional on the named upstream ingredients. No full-scale timing or independent peer review is claimed.
