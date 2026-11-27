#set document(title: "A polynomial certificate for deterministic elliptic-curve relation systems", author: "elliptic-relation-counting research project", date: datetime(year: 2026, month: 11, day: 27))
#set page(paper: "a4", margin: (x: 23mm, top: 22mm, bottom: 20mm), footer: context [
  #set text(size: 8.5pt, fill: rgb("#46565b"))
  Fixed-shift relation systems #h(1fr) #counter(page).display()
])
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: none)
#show heading.where(level: 1): it => block(above: 0.7em, below: 0.45em, text(size: 17pt, fill: rgb("#184856"), it.body))
#let small(body) = text(size: 9pt, fill: rgb("#46565b"), body)
#small[RESEARCH NOTE  ·  27 NOVEMBER 2026]
#v(7pt)
#block[#set par(justify: false); #text(size: 25pt, weight: "bold")[A polynomial certificate for deterministic elliptic-curve relation systems]]
#v(3pt)
#small[elliptic-relation-counting research project  ·  AI-assisted working manuscript]
#v(11pt)
*Abstract.* A universally decomposing factor base determines the elliptic discrete logarithm through a family of fixed-shift linear equations. We encode all these equations in a matrix polynomial of full generic rank. Each matrix entry can be evaluated with at most two slot pins per contraction and small public weights on address bits. Assuming the closed Boolean tensor-contraction theorem used by Liu, Zhang and Ren, one trial takes $2^(O(b/(log b)^(1/4)))$ bit time on every eligible input. A single shared challenge of $O(b)$ random bits gives verified success with probability $1-2^(-Omega(b))$. Uniform deterministic selection of the challenge remains open. We also give a deterministic improvement to their partial expansion by varying the number of vertex groups at the root.

= 1. What the certificate guarantees

Let $p>3$ be prime, $b=ceil(log_2 p)$, $L=log_2 b$, and let the rational points of a nonsingular elliptic curve form a group of prime order $r$. Given a generator $P$ and a target $Q$, write $lambda$ for the table of logarithms of the factor-base columns and two distinguished columns containing $P,Q$. All complexity estimates charge bit operations on a uniform multitape machine.

*Theorem 1.* Assuming LZR Theorem 2.2, one can construct an implicit matrix $B$ and vector $h$ over a polynomial ring over $bb(F)_r$ such that $det B$ is nonzero and $h=B lambda$. At any specified point, evaluation, a rank test and, on full rank, solution and verification take $2^(O(b/L^(1/4)))$ bit time and at most that many bits of space. Independent trials with fair bits give a Las Vegas algorithm on this time scale using $O(b)$ expected random bits. A deterministic algorithm with $O(b)$ bits of advice per input length has the same worst-case time.

The hypothesis is the closed Boolean contraction bound. The result uses neither the growing-alphabet finite-packet assumption of the repository's v4 manuscript nor its sparse randomized solver. It concerns the whole prime-order group. The primary uniform deterministic time target would require an efficiently computable hitting list.

#pagebreak()
= 2. Full generic rank from ordered addresses

Use the $m asymp L^(1/4)$ blocks of #link("https://eprint.iacr.org/2026/3020.pdf#page=11")[LZR Theorem 3.1]. Each block has a compact membership test and an injective $k$-bit address, with $k=Theta(b/m)$. Every group element has an ordered representation with two entries from each block. Let $N$ count all labeled block elements and the two distinguished columns. Enumerate the blocks once. For each of the $N$ columns, including $P$ and $Q$, use the theorem's deterministic decomposition of its point. The proof fixes the tuple bits in order using $O(m b)$ prefix-pinned contractions per target (#link("https://eprint.iacr.org/2026/3020.pdf#page=14")[LZR, p. 14]). Equate each column to its tuple sum and normalize the column for $P$ to one, giving $C=N+1$ explicit rows.

For every pair of representations $(X,Y)$ with $"sum"(Y)="sum"(X)+P$, add conceptually a row with coefficient vector $a="cnt"(Y)-"cnt"(X)$ and right side one. Repeated slots count with multiplicity. Any solution agrees on all representations of one point: they have a common successor. The induced function increases by one under addition of $P$ and equals one at $P$. The explicit consistency rows therefore force every true label. The resulting matrix $A$ has column rank $N$ over $bb(F)_r$.

Give explicit rows distinct keys $0,...,C-1$. A shift pair has key $C+v$, where $v$ is the integer obtained by concatenating its $beta=4m k$ address bits. For row $i$ write $e_i$ for its key and $c_i$ for its right side. Define
$ B(Z)=sum_i Z^(e_i) a_i a_i^t, quad h(Z)=sum_i Z^(e_i)c_i a_i. $
Cauchy–Binet expresses the determinant as a sum of squared basis minors times $Z^(sum_(i in I)e_i)$. The greedy basis that takes the first independent rows in increasing key order reaches every prefix rank first. Its ordered keys are componentwise no larger than those of any other basis, with strict inequality somewhere. Its total weight is uniquely minimal. The coefficient at that degree is a nonzero square. This gives full generic rank even in finite characteristic. Termwise consistency gives the identity for $h$.

= 3. A challenge in the prime field

Set $d=max(1,floor(b/4))$, $D=2^d$. Choose $s$ so that $D^s$ covers all addresses and exceeds $C$, taking the least such positive integer. Write $M_v (X)=product_j X_j^(v_j)$ for the base-$D$ digits of $v$. Use $M_e$ on explicit rows and $M_C M_v$ on pair rows. Substituting $X_j=Z^(D^j)$ recovers the isolated basis term. The total determinant degree is at most $2 N s(D-1)$.

At the balanced arity $s=O(1)$ and $log N=O(b/m+log b)$. Hasse's bound gives $D^2<r$ for large $b$. Draw each coordinate using exactly $2d$ bits from $[0,D^2)$. The Schwartz–Zippel bound gives failure at most $2 N s(D-1)/D^2=2^(-b/4+O(b/m+log b))$. Detect failure by elimination, and retry independently. All inequalities can be tested directly; finitely many smaller lengths use exhaustive arithmetic.

#pagebreak()
= 4. From modular weights to charged bit time

Each component $a_j$ is a signed sum of at most $4m$ single-slot equality indicators. Its product with $a_k$ requires at most $16m^2$ two-slot pins. Keep signs outside the network and use only single pins for the right-hand side. A contradictory pair of pins contributes zero. Each group addition has a unique internal extension for valid inputs, including identity operands, inverse pairs and doubling. It uses $O(b)$ edges and has maximum vertex degree $O(b^(3/4))$ (#link("https://eprint.iacr.org/2026/3020.pdf#page=10")[LZR, Lemma 2.5]). Thus the $O(m)$ additions preserve ordered-pair multiplicities.

Partition the address within each public base-$D$ digit and each slot into contiguous pieces of length at most $ceil(sqrt(b))$. A piece with value $u$ starting at bit $j$ in digit $i$ contributes the canonical residue of $X_i^(2^j u)$, computed by repeated squaring. There are $q=O(sqrt(b)+m+s)$ pieces. Each contributes one small-arity vertex on dedicated address wires. With subdivision this gives $E=O(m b)$ edges, degree $O(b^(3/4)+sqrt(b))$, and polynomial description and entry time. The scalar $M_C$ is public.

Choose $S=Theta(m b)$ with $S>=E+b+q+1$. Take $t=floor(S/(q+1))>=1$, so $q t<=S$, and split each canonical field weight into $J=ceil(ceil(log_2 r)/t)$ nonnegative base-$2^t$ digits. For each of $J^q$ digit choices contract the network with the selected local digits. Its exact integer answer is below $2^(E+q t+1)$. Boolean edges carry bits, and vertex entries may lie in $ZZ slash 2^(P_0) ZZ$ at truncation length $ell=1$ (#link("https://eprint.iacr.org/2026/3020.pdf#page=7")[LZR, Theorem 2.2, pp. 7–8]). Hence the digit weights are admissible. At precision $P_0=E+q t+1=O(S)$ the least residue equals the integer answer. Multiply by the digit shifts modulo $r$ and sum. This respects the linear precision cap and uses no homomorphism between the two coefficient rings.

The deterministic seed requires polynomially many counts of size $O(m b)$. Enumeration and constructing all explicit rows add a polynomial factor in $N,b$. The digit enumeration adds $O(q log(b+2))$ to the log cost; all $N^2$ matrix entries, the dense solve and scalar verification add $O(log N+log b)$. Altogether, for growing polylogarithmic $m$,
$ log_2 T_"trial" = O(b/m+(m b)/sqrt(log(m b))+(sqrt(b)+m)log(b+2)). $
Balancing the first two terms proves Theorem 1. Theorem 2.2 bounds the space for each sequential contraction by a polynomial in $S$ and the network description length. The matrix and lists use $op("poly")(N,b)$ bits.

At a given length there are at most $2^(O(b))$ canonical curve inputs. A sufficiently large constant number of independent trial tapes succeeds on all inputs with positive probability by a union bound. Storing these tapes proves the advice assertion. No uniform procedure for finding them in the desired budget follows.

*One-candidate version.* Sum $(a_i^t y-1)M_(v_i)$ over shift pairs. Its distinct monomials prevent symbolic cancellation for any wrong candidate satisfying the explicit equations. After a nonzero specialization, the values for the two next prefix branches sum to the current value. One can therefore extract a violated relation deterministically using only single-point pins and public address weights.

#pagebreak()
= 5. A shorter deterministic partial expansion

Write $S=Theta(m b)$ and $ell_0=ceil(log_2 S)$. In LZR Lemma 5.3 the total internal mass at depth $j$ is at most $2^j S$. With $u_j asymp ell_0/2^j$ the cost of distributing siblings adds up to $O(S 2^(3D_0/2)/sqrt(ell_0))$ for $D_0=ceil(log_32(4m))$. The extra power of $m$ comes from enumerating complete products.

*Theorem 2.* Take $m=Theta((log b)^(1/4))$ and $beta=Theta(b)$ terminal bits on distinct hosts. Under the explicit Appendix D construction, the pointwise expansion admits $O(m)$ stopped tensors per term and
$ log_2 abs(cal(I)) = O(S sqrt(log(2m)/ell_0)+ell_0). $
Each tensor has $O(b/m)$ terminal bits and a pinned evaluation costs $2^(O(S/sqrt(ell_0)))$ times a polynomial in the description length. At the original arity the expansion loses a factor $sqrt(log log b)$ over the target exponent.

*Proof.* Use $R$ groups at the root only, taking the smallest power of two above or equal to $max(64,4m)$. Set $u=floor(frac(epsilon ell_0, log_2 R))$, $a=2u$, $W=2^a$, with $epsilon=2^(-17)$. At this arity $u>=1$ for sufficiently large $b$. Handle the finite initial range directly. The square of the pair-label array size is $R^(4a)<=2^(8epsilon ell_0)<=(2S)^(8epsilon)$. The packing guard then holds at an active root with $M>S/sqrt(ell_0)$ internal edges and maximum degree at most $ceil(S^(9/10))$. Conditional expectation balances the squared edge-class deviations and the terminal-bin deviations. Distinct hosts give at most $beta/R+sqrt((a+2)beta)$ terminals per bin. This is $O(beta/m)$ because $R>=4m$ and $m^2 (a+2)=O(beta)$ at the stated arity. The padded common class size is at most $3M/R^(2a)$, hence at most $3M/W$ packets and at most $6M/R$ real child edges per group.

Within any packet and any active coordinate there are still only two possible symbols. At a Hermite node of valuation $i$, $D_t=W/2^i$ and there are at most $6M/D_t$ crossing modes per weight. At the active root $P_0=O(S)=O(M sqrt(ell_0))$. At this arity $2^u/u$ grows faster than any fixed power of $ell_0$. Thus $D_t>=2^u$ gives $p_t=P_0+8ceil(3M/u)D_t+8M(a+1)=O(M D_t/u)$. The weighted partition estimate gives the claimed log cardinality after one step. Keep all $R$ children as canonical stopped tensors. Pinned evaluation continues with 64 groups and the original scalar threshold, starting at absolute depth one with $a_j=a+2sum_(i=1)^j floor(epsilon ell_0/2^i)$. The inherited input precision is at most $S ell_0^(O(1))2^(a_(j-1))$. The next tile scale is at least $2^(a_(j-1)+u_j)$, so the mode bound persists. The per-path exponents sum geometrically; creation-precision representatives and the original Hermite guard suffice. Early small calls can be directly summed. Full bounds, including digit reconstruction, appear in the branch supplement.

The residual graph retains one degree-$Theta(b/m)$ lookup per slot. Theorem 2 alone therefore leaves deterministic end-to-end time at the original double-logarithmic saving.

#pagebreak()
= 6. Verification and the remaining decision

On a cyclic fixture of order five, 307 of 625 test points yield a nonsingular six-column matrix. The all-ones point has rank five. Direct outer products agree with two-pin evaluations throughout. All 490 nonzero evaluations for a chosen incorrect candidate produced violated rows. On $y^2=x^3+4x+1$ over $bb(F)_13$, direct group addition gave 379 ordered shift pairs. Each of 17 nonsingular univariate specializations over $bb(F)_19$ recovered all eight labeled columns. These finite checks exercise cancellation, repetitions and witness recovery.

The explicit single-point rows alone can admit incorrect labels. Take two blocks each containing all nonzero elements of the cyclic group of order eleven. Represent $5x$ by four copies of $4x$, two per block. Multiplication by four has two orbits, independently scalable in a candidate table. The resulting 20-column matrix has nullity two. Scaling one orbit preserves normalization at one and passes the single-point rows when $x$ is represented by four copies of $x/4$. This tests the abstract matrix criterion. Additional restrictions from a particular sparse seed could help.

The univariate key ordering proves nonvanishing but yields exponential degree in the total boundary width. Evaluating one supplied challenge is inexpensive, whereas enumerating the resulting interpolation grid exceeds the target budget. It would suffice to construct a short list on which the matrices, stacked vertically, have full column rank. An explicit list with total construction and test cost $2^(O(b/L^(1/4)))$ would finish the uniform deterministic reduction.

= Sources and reproducibility

#block[#set par(justify: false)
Shuchang Liu, Bingsheng Zhang and Kui Ren. _A Deterministic Subexponential-time Algorithm for the Discrete Logarithm Problem._ 23 November 2026. Sections 2–5. #link("https://eprint.iacr.org/2026/3020")[Preprint 2026/3020.]

OpenAI. _The Exponential Time Hypothesis is False._ 25 October 2026. Appendix D, quantitative contraction and its uniform resource hypotheses. #link("https://github.com/openai/math/blob/3344da59991ed94f9ba711217b29c1d760badfa6/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex")[Pinned contraction source.]

H. Narayanan, H. Saran and V. V. Vazirani. _Randomized parallel algorithms for matroid union and intersection, with applications to arborescences and edge-disjoint spanning trees._ Symbolic determinant isolation is classical. #link("https://ics.uci.edu/~vazirani/NSV.pdf")[Authors’ paper.]

_Counting relations in independently hashed elliptic-curve slots._ Repository version 4, 24 November 2026. The expected-time coefficient there relies on finite-packet Assumption A.1. #link("https://github.com/davidiach/elliptic-relation-counting/blob/fc3a8019f1763bd15250b34cc9c15cfabda8131c/paper.pdf")[Prior manuscript.]
]

#v(3pt)
#small[The accompanying audit records the public Lean statements and their premises. The proofs in this note have not been kernel-checked. The finite arithmetic checks, transcript and complete parameter derivations are included in the research branch.]

#pagebreak()
= Addendum. Deterministic evaluation lists

Let $E$ be the $C=N+1$ explicit rows. For $0<=z<C$, evaluate the original matrix at $(z,z^D,...,z^(D^(s-1)))$. Every shift incidence vector has sum zero. Therefore summing all the rows of the evaluated matrix gives $sum_(e=0)^(C-1) z^e (bold(1)^t a_e)a_e^t$. Its scalars on the explicit rows are one or $1-2m$. For $C<=r$ and $r$ not dividing $2m-1$, Vandermonde interpolation puts every explicit row inside the evaluation stack.

*A second list.* Set $K=C$, $H=2^(floor(b/3))$ and $J=1+(s-1)binom(K,2)$. Construct a primitive root $gamma$ modulo the group order $r$ using LZR Proposition 3.3 under the closed theorem. For $1<=t<=J$, $0<=h<K$, take $(gamma^h,gamma^(h t),...,gamma^(h t^(s-1)))$. If a direction $u$ survives this stack, every nonzero polynomial combination of the entries of $B(X)u$ with individual degrees at most $H$ has more than $K$ terms. Indeed, differences of exponent vectors exclude at most $s-1$ integer slopes per pair. At a remaining slope the projected exponents are distinct and below $r-1$, since $H sum_(j=0)^(s-1) J^j<r-1$ eventually. The geometric samples form a Vandermonde matrix on their distinct frequencies. Both lists together fit the target bit budget, including testing the ranks. Uniform full column rank on all inputs is open.

*Density of witnesses.* Let $D_0$ be the group values of one-per-block tuples, $d=abs(D_0)>7r/8$, and $q=2d-r$. Each group element has at least $q$ full tuples by intersection of two translates. Every direction beyond the block constants on factor-base columns violates at least $q^2/2$ ordered shifts. To see this, a fiber with two different label sums yields at least $q$ violations by comparison with every successor. If some half-tuples of equal value have different sums, appending a common half-tuple gives at least $d$ mixed fibers. Otherwise half-tuple sums define $h$ on $D_0$. If $h(x+t)-h(x)$ takes two values at $a_0,a_1$, compare each eligible $x$ to an $a_i$ with unequal difference. The fiber at $x+t+a_i$ is mixed, with multiplicity at most two in this count. There are at least $q/2$ such fibers. If all differences are constant, intersections of three translates yield their additivity. The labels then have the form $kappa log_P A+c_i$ on each block. A nonzero slope violates every homogeneous shift. For an affine direction satisfying the homogeneous explicit rows, consistency gives $c_i=2 sum_j c_j$ and normalization at $P$ gives $kappa+2 sum_j c_j=0$. Hence $c_i=-kappa$ and $(2m-1)kappa=0$ in $bb(F)_r$. When $r$ does not divide $2m-1$, the explicit rows alone force all labels, including the two distinguished coordinates, to zero.

Large witness counts do not control modular cancellation. For any prescribed $L$ points, omitting at most $(r-1)L(N(N+1)/2+N)$ shift rows can make the remaining shifts' matrix and right side vanish at all the points. In the target budget, the remaining shifts still span the original shift space by the density bound. The proof uses the zero-sum subsequence bound for elementary abelian groups. The full elliptic family includes every ordered edge between neighboring fibers; controlling its exact fiber sum is the remaining problem. Complete proofs are in the deterministic-points supplement.

#pagebreak()
= Exact neighboring fibers

*Full-grid theorem.* On every eligible input the stack over the whole base-field grid has full column rank. For $r>2$ the individual degrees of $B(X)u$ are at most $2(D-1)<r$ by Hasse's bound and the definition of $D$. A constant direction killed everywhere therefore kills each monomial coefficient. Each ordered neighboring pair occurs exactly once in its fiber, with a distinct monomial and coefficient $a a^t$. The explicit and shift rows together force $u=0$. For the published arity the short-list guard fails whenever $m<=45$: at $b>=40$ coverage gives $s>=16$, and $J>=(7/4)(s-1)r^(2/m)$. Consequently $H J^(s-1)>r-1$ throughout this range. Smaller bit lengths give the same conclusion directly.

*The order-two case.* On $y^2=x^3+2x$ over $bb(F)_5$ choose $P=(0,0)$. The two blocks are the whole group and the canonical tuple for $A$ is $(A,O,O,O)$. The unguarded formulas yield only the all-zero and all-one points. At all ones every tuple fiber has mass eight and every first-moment coordinate is zero modulo two. The entire neighboring-fiber Gram contribution vanishes. Two equal explicit rows cancel, leaving rank five of six, with kernel generated by $(1,1,0,0,0,0)^t$ in the column order $(O_1,P_1,O_2,P_2,P_*,Q_*)$. At the grid point $(1,0,...,0)$ the consistency row for $O_1$ detects this vector. Thus the prescribed fallback succeeds. Both targets and all 256 grid points were checked directly.

*Displacements.* Let $G_T$ be the unweighted shift matrix at displacement $T$, with block sizes $n_i$ and $M=product_i n_i^2$. Summing the neighboring fibers over all $T$ counts each ordered pair of independent tuples once. On the base columns their covariance gives
$ sum_T G_T = 4 M^2 ⨁_i (I/n_i - bold(1)bold(1)^t/n_i^2). $
When $0<n_i<r$ and $r$ is odd this has rank $N_0-m$. Prefix pinning the elliptic coordinates of $T$ selects at most $N_0-m$ displacements attaining this joint rank within the target budget under closed contraction. Their right sides involve the unknown logarithms of the selected points.

*Fixed coefficient.* At a torus point the weighted displacement counts are the coefficients of a product in $bb(F)_r[U]/((U-1)^r)$. A slot on $n$ different group elements has vanishing order at one at most $n-1$, by the binomial Vandermonde matrix on their group exponents. In the guarded regime the sum of these orders is at most $4(N_0-m)<r$, even with up to two slots pinned. The product survives, but its coefficient at $P$ can vanish. Controlling the Gram sum of those fixed coefficients across the guarded list remains open. Proofs and finite controls appear in the neighboring-fiber supplement.

#pagebreak()
= Cancellation at the prescribed generator

*Theorem.* In the guarded regime, partition the terms of the all-ones shift matrix by the source group value and by the point in one specified slot on each side. The slot may be the first occurrence of a largest block. The resulting matrices have joint rank $N_0-m$ before adjoining the explicit rows. The same rank holds when one sums every successor separately for each ordered source.

*Pin counts.* Write $d$ for the one-per-block count on the group and $cal(E)=sum_g d(g)^2$. If a block of size $n_i$ is omitted, the resulting count $d_(-i)$ has squared norm at most $cal(E)/n_i$, by retaining the diagonal terms in the square of the sum of its $n_i$ translates. Cauchy–Schwarz now bounds each one-slot pin count in a full tuple fiber by $cal(E)/sqrt(n_i)$. The actual Stage I potential gives $cal(E)<321r$: for guarded $m>=46$ the tuning gives $r^(1/m)>2^20 m$, so the two ceilings in (3.2) yield $mu<=16r exp(3m/(16r)^(1/m))<17r$. The largest block has $n_i>(7r/8)^(1/m)>321^2$. Every nonempty pinned class therefore has size between one and $r-1$. Some pin occurs at every group value. Twice its count is a nonzero coordinate of the full fiber's incidence moment $v_g$, even if its mass vanishes modulo $r$.

*Summing both sides.* A class $B$ specifies the tuple value and the selected pin. Let $z_B,V_B$ be the mean and covariance of its incidence vector in the field. If $B,D$ lie at $g,g+P$, their full Gram sum divided by the product of their sizes is
$ V_B+V_D+(z_D-z_B)(z_D-z_B)^t. $
The true label vector satisfies $lambda^t V_B=0$ and $lambda^t (z_D-z_B)=1$. Thus any common kernel direction $u$ gives equal mean scores on adjacent classes. The displayed identity then gives $V_B u=-V_D u$. All classes are connected, and one class at each group value forms an odd cycle. Hence all mean scores are equal and every class covariance kills $u$. Expanding the unconditional covariance by classes shows that it too kills $u$. Independence of the two slots per block gives this covariance as $2 ⨁_i (I/n_i-bold(1)bold(1)^t/n_i^2)$, with zero dedicated columns. Its kernel consists of the block constants and the free dedicated coordinates, which annihilate every individual shift.

*Summing the successors.* For a fixed source $X$ let $K_X$ sum all successor Gram terms at $g=op("val")(X)+P$, and let $H_g$ be the successor second moment. The identity $(I+c_X lambda^t)K_X=H_g-v_g c_X^t$ forces any two source tuples of the same value to have equal scores under a common kernel direction, since $v_g!=0$. Calling the score $f(g-P)$ gives $H_g u=v_g f(g)$ and hence $f(g)=f(g-P)$. The score is constant and varying an individual slot proves the same rank as above.

Each specified matrix fits the closed contraction budget by pinning slots and testing group values. Retaining the classes for every group value requires at least $r$ matrices. In the order-19 arithmetic check the class matrices have rank eight together and their sum has rank seven. This instance is outside the guard. The rank of the prescribed short list, which couples the source and successor evaluations, remains open.
