#let title = [A polynomial certificate for deterministic elliptic-curve relation systems]
#let project-credit = "David Iach, project lead and repository maintainer"
#set document(
  title: "A polynomial certificate for deterministic elliptic-curve relation systems",
  description: "Research report of the Elliptic Relation Counting Project. " + project-credit + ".",
  date: datetime(year: 2026, month: 11, day: 27),
)
#set page(paper: "a4", margin: (x: 29mm, top: 23mm, bottom: 23mm), numbering: "1",
  header: context if counter(page).get().first() > 1 [
    #set text(size: 8.5pt, fill: luma(85))
    Fixed-shift relation systems #h(1fr) Elliptic Relation Counting Project
  ])
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.64em, spacing: 0.78em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 14pt)
#show heading.where(level: 2): set text(size: 11pt)
#show heading: set block(above: 1.25em, below: 0.48em)
#set math.equation(numbering: none)
#set figure(gap: 0.6em)
#show figure.caption: set text(size: 9.5pt)
#let repository = "https://github.com/davidiach/elliptic-relation-counting/"
#let revision = "d9488fa24a53344b05cbe4f701f52f8070c3ef3d"
#let artifact(path, body) = {
  let route = if path.ends-with(".md") { "blob" } else { "tree" }
  link(repository + route + "/" + revision + "/" + path, body)
}

#align(center)[
  #set par(justify: false)
  #text(size: 21pt, weight: "semibold")[#title]
  #v(0.5em)
  #text(size: 11pt)[Technical report of the Elliptic Relation Counting Project] \
  #text(size: 10pt)[27 November 2026]
]
#v(0.55em)
#block(inset: (x: 6mm))[
  #set text(size: 10pt)
  *Abstract.* Given a prime-order elliptic-curve group over a prime field, the relations between ordered representations of neighboring points determine a table of discrete logarithms. Directly summing their outer products can lose rank in the scalar field. We assign each ordered pair a monomial weight and prove that the resulting matrix polynomial has nonzero determinant. Its entries require counts with at most two prescribed slots. Under the closed Boolean tensor-contraction bound of Liu, Zhang and Ren, a single random specialization shared by all entries yields a Las Vegas algorithm whose expected bit time has binary logarithm $O(b slash (log_2 b)^(1 slash 4))$ and whose expected random-bit cost is $O(b)$, where $b$ is the field bit length. The full field grid suffices deterministically. A uniform choice of evaluations with the same time bound remains open.
]

#block(inset: (x: 6mm))[
  #set text(size: 10pt)
  *Project credit.* #project-credit. \
  *AI use and verification.* AI agents helped develop the mathematical arguments and experiments. OpenAI Codex edited the manuscript and reproduced the finite counts and ranks for the fixed-shift matrices, evaluation lists, and neighboring fibers, as well as every entry in @pinned. @evidence gives the data and programs. The project has no recorded independent human verification of the general proof of Theorem 1.
]

= Introduction

A fixed displacement gives a linear equation for each pair of representations of neighboring group elements. Liu, Zhang and Ren (LZR) show that a factor base covering the group supplies enough such equations to determine every label @lzr. There can be far more equations than unknowns. One way to compress the system is to sum the outer products of its rows. Over a finite field this can destroy a direction even when the original rows span the entire space. On a six-column cyclic example, the ordinary Gram matrix has rank five.

We attach a weight to each relation using the addresses of its slots. The ordered addresses distinguish a basis term in the determinant, while the digits of each address factor into weights that the counting algorithm can evaluate. This yields a certificate for the entire table at one specialization. It uses the closed contraction bound in LZR, Theorem 2.2. LZR's uniform deterministic end-to-end bound, $2^(O(b slash sqrt(log_2 log_2 b)))$, also invokes the terminal-preserving construction from Appendix D of the upstream work @upstream. The earlier project report treats independently hashed slots under a different finite-packet assumption @prior.

*Theorem 1 (fixed-shift certificate).* Let $p>3$ be prime and let $G=E(bb(F)_p)$ be the group of rational points of a nonsingular elliptic curve, of prime order $r$. Given a generator $P$ and a target $Q in G$, put $b=ceil(log_2 p)$ and $L=log_2 b$. Assume LZR's uniform closed Boolean tensor-contraction bound. There is an implicit square matrix polynomial $B$ over $bb(F)_r$ and a polynomial vector $h$ such that $det B != 0$ and $h=B lambda$, where $lambda$ is the table of base-$P$ logarithms, including the target column. At any specified field point, evaluation, rank testing and solution at full rank take $2^(O(b slash L^(1 slash 4)))$ bit time and at most that many bits of space. Independent trials give a verified answer in the same expected time using $O(b)$ expected fair bits. With $O(b)$ bits of advice per input length, the time bound holds deterministically in the worst case.

Throughout, logarithms in resource bounds are binary and costs count bit operations on a multitape machine. A _pin_ fixes the value in a specified slot of an ordered tuple. A column represents an occurrence in a labeled block, so a point shared by two blocks has two columns.

= Relations and the weighted Gram matrix <relations>

*The explicit equations.* LZR, Theorem 3.1, supplies $m asymp L^(1 slash 4)$ blocks with injective $k$-bit addresses, $k=Theta(b/m)$. Every $g in G$ has an ordered representation $X$ with two slots from each block. Write $g(X)$ for the group sum and $c_X$ for the column vector of multiplicities among the labeled block elements. Extend $c_X$ by zeros on two distinguished columns for $P$ and $Q$. Write $N_0$ for the total block size and $N=N_0+2$. For each column $j$, the Stage I procedure supplies a representation $U_j$ of its group element by fixing tuple bits with $O(m b)$ prefix counts. The $C=N+1$ explicit equations on a label vector $y$ are
$ y_j=c_(U_j)^t y quad (1 <= j <= N), quad y_(j_P)=1. $
Here $j_P$ denotes the distinguished generator column. The normalized true table satisfies these equations even when several columns refer to the same group element.

*The shift equations.* For every ordered pair $(X,Y)$ with $g(Y)=g(X)+P$, take the row $a=c_Y-c_X$ with right side one. Let $A$ contain the explicit and shift rows, written $a_i^t$, and let $c_i$ be their right sides. Any solution assigns equal scores $c_X^t y$ to all representations of the same value, because they share a successor. The resulting function $f$ on $G$ satisfies $f(g+P)=f(g)+1$. Consistency and normalization give $f(P)=1$, so $f(j P)=j$ throughout the group. The explicit equations now fix every column. Hence $A$ has column rank $N$.

*Isolation.* Give the explicit rows keys $0,...,C-1$. Concatenating the $4m$ slot addresses of a shift pair gives an integer $v$ of width $beta=4m k$. Assign that pair the key $C+v$. With $e_i$ the key of row $i$, define
$ B(Z)=sum_i Z^(e_i) a_i a_i^t, quad h(Z)=sum_i Z^(e_i)c_i a_i. $
Consistency gives $h=B lambda$ term by term. By Cauchy–Binet,
$ det B(Z)=sum_(I: abs(I)=N) det(A_I)^2 Z^(sum_(i in I)e_i). $
Scan rows by increasing key, keeping a row whenever it increases rank. At every prefix, the number kept equals its rank. Thus the $j$th retained row occurs no later than the $j$th row of any competing basis, for every $j$. A different basis has a strict inequality at some position. The minimum sum of keys belongs to just one basis, whose squared minor is nonzero in $bb(F)_r$. This use of determinant isolation for represented matroids has a classical precedent @matroid. The ordered addresses allow it to be applied to the implicit elliptic relation system.

= Specialization and bit complexity <specialization>

*Sampling in the scalar field.* Put $d=max(1,floor(b/4))$, $D=2^d$, and take the least positive $s$ such that $D^s>=2^beta$ and $D^s>C$. For $v<D^s$, with base-$D$ digits $v_j$, set $M_v(X)=product_(j=0)^(s-1) X_j^(v_j)$. Use $M_e$ for an explicit row and $M_C M_v$ for a shift pair. The common factor $M_C$ lets the private address retain its original digits. Substituting $X_j=Z^(D^j)$ recovers the univariate construction. Consequently $det B(X)$ is nonzero and its total degree is at most $2 N s(D-1)$.

Here $s=O(1)$ and $log N=O(b/m+log b)$. For large $b$, Hasse's bound gives $D^2<r$. Choose each coordinate independently from $0,...,D^2-1$ using exactly $2d$ fair bits. The Schwartz–Zippel bound @schwartz @zippel gives
$ Pr[det B(X)=0] <= frac(2 N s(D-1),D^2) = 2^(-b/4+O(b/m+log b)). $
A rank test detects failure. On success, solving gives $lambda$ and the proposed target logarithm can be checked by scalar multiplication. Independent retries therefore have zero error. Finitely many smaller bit lengths can use exhaustive group arithmetic.

*Closed counts with two pins.* Each component of a shift row is a signed sum of at most $4m$ slot indicators. Expanding $a_j a_k$ gives at most $16m^2$ weighted counts with two pins. Contradictory pins contribute zero. A single pin suffices for the right side. The group-law gadgets of LZR, Lemma 2.5, have a unique internal extension for each valid input, including identities, inverses and repeated points. Each addition has $O(b)$ edges and maximum degree $O(b^(3/4))$. These properties preserve multiplicity when counting ordered tuples.

Cut each digit at slot boundaries and into pieces of length at most $ceil(sqrt(b))$. A piece of value $u$ starting at bit $j$ in digit $i$ contributes the field weight $X_i^(2^j u)$, with $0^0=1$. Binary powering computes its canonical residue in polynomial bit time. The $q=O(sqrt(b)+m+s)$ pieces attach to address wires. After subdivision the network has $E_0=O(m b)$ Boolean edges, degree $O(b^(3/4)+sqrt(b))$, and polynomial description length.

The contraction returns residues modulo powers of two. To recover the desired scalar modulo $r$, choose $S=Theta(m b)$ with $S>=E_0+b+q+1$ and let $t=floor(S/(q+1))$. Split each of the $q$ nonnegative field representatives into $J=ceil(ceil(log_2 r)/t)$ digits of $t$ bits. For each of the $J^q$ choices of digit positions, the corresponding network has an integer count below $2^(E_0+q t+1)$. Apply the closed theorem at truncation length one and precision $E_0+q t+1=O(S)$. The least dyadic residue is the exact count. Multiply it by the appropriate power of $2^t$, sum modulo $r$, and apply the public factor $M_C$. This charges the precision required to pass between the two rings.

*Time, space and advice.* Polynomial factors in $N,b$ cover the explicit rows, all $N^2$ entries, elimination and verification. Each closed contraction costs $2^(O(S slash sqrt(log S)))$ up to polynomial factors. Digit choices add $O(q log(b+2))$ to the log cost. Thus for growing polylogarithmic $m$ in the stated regime,
$ log_2 T_"trial" = O(b/m + (m b)/sqrt(log(m b)) + (sqrt(b)+m)log(b+2)). $
Balancing the first two terms gives $m asymp L^(1 slash 4)$. A contraction uses polynomial space in its input size. Sequentially computing entries leaves the lists and the dense matrix as additional storage, polynomial in $N,b$.

At each length there are at most $2^(O(b))$ canonically encoded inputs. Pad trial tapes to a common length of $O(b)$. Since a fixed input fails a trial with probability $2^(-Omega(b))$, a sufficiently large constant number of independent tapes has positive probability of succeeding on every input of that length by a union bound. Such tapes give the advice in Theorem 1. The argument counts tapes without constructing them uniformly.

*Checking a single candidate.* For an incorrect $y$ satisfying the explicit rows, the sum of $(a_i^t y-1)M_(v_i)$ over shift pairs is nonzero because its monomials are distinct. At a specialization with nonzero value, the two extensions of a given prefix sum to the parent value. Following a nonzero extension recovers a violated shift equation. This variant requires only single-slot pins.

= Deterministic evaluation and cancellation <deterministic>

For a list $cal(H) subset.eq bb(F)_r^s$, stack the matrices $B(x)$ vertically for $x in cal(H)$. The right sides are consistent. The stack determines the table whenever the intersection of the kernels is zero, even if every member is singular.

== The full grid and a shorter list

*Proposition 2 (full grid).* The stack over $bb(F)_r^s$ has full column rank on every eligible input. For $r>2$, each entry of $B(X)$ has individual degrees at most $2(D-1)<r$. If $B(x)u=0$ everywhere, interpolation in one variable at a time gives $B(X)u=0$. Distinct row keys give distinct formal monomials, as the substitution in @specialization shows. The coefficient at row $i$ is $a_i(a_i^t u)$, so every row annihilates $u$. Full rank of $A$ forces $u=0$.

For $r=2$, Hasse's bound forces $p=5$, and the prescribed blocks are both the whole group. The canonical representation of a point $A$ is $(A,O,O,O)$, with $O$ the identity. On the Boolean grid use $X_j^2=X_j$. Each shift monomial then has support containing the binary support of $C$, and its support interpreted as an integer is at least $C$. The explicit rows keep their distinct coefficients at keys below $C$. They already have full rank: the canonical tuple for the identity has zero incidence vector modulo two, so both identity labels vanish. The canonical tuple for the generator then has score equal to its first-block label. Consistency equates this score to the other generator columns, and normalization fixes it. The same argument in the homogeneous system forces every component to zero.

The full grid has $r^s$ points. A smaller list starts with $(z,z^D,...,z^(D^(s-1)))$ for $0<=z<C$. Multiplying $B$ on the left by the all-ones row removes every shift term. The remaining row sums are $1$ or $1-2m$. If $C<=r$ and $r$ does not divide $2m-1$, Vandermonde interpolation recovers every explicit row from this stack.

Put $K=C$, $H=2^(floor(b/3))$, and $J_0=1+(s-1)binom(K,2)$. With a primitive element $gamma$ modulo $r$ from LZR, Proposition 3.3, add the points with coordinates $X_j=gamma^(h t^j)$ for $1<=t<=J_0$, $0<=h<K$, $0<=j<s$. Suppose
$ H sum_(j=0)^(s-1) J_0^j < r-1. $
For any nonzero polynomial with at most $K$ monomials and individual degree at most $H$, each pair of exponent vectors rules out at most $s-1$ slopes. At a remaining slope their projected exponents are distinct integers below $r-1$. The geometric samples give a Vandermonde matrix on distinct frequencies, as in sparse interpolation @sparse. Thus a surviving direction $u$ annihilates the explicit rows, and every nonzero polynomial in the ideal generated by the coordinates of $B(X)u$ with those degree bounds has more than $K$ terms. Construction and testing of the short list fit the time of Theorem 1. Whether every nonzero residual direction admits such a sparse consequence is open.

The parameter tests hold eventually at the published tuning. The no-wrap inequality above fails for $m<=45$. For $b>=40$, coverage gives $s>=16$ and $J_0>=(7/4)(s-1)r^(2/m)$, which already contradicts the inequality in that range of $m$. The finite smaller values select the same full-grid fallback. The #artifact("NEIGHBORING_FIBERS.md")[neighboring-fiber supplement] gives the bounds for both ranges.

== Why many violated equations can still cancel

Let $D_0$ be the set of sums using one slot per block. Stage I gives $d_0=abs(D_0)>7r/8$. Write $q_0=2d_0-r$. Intersecting $D_0$ with a translate of $-D_0$ gives at least $q_0$ full tuples at each group value. Every direction beyond the block constants on the base columns violates at least $q_0^2/2$ ordered shifts.

Indeed, any fiber with two unequal tuple scores gives at least $q_0$ violations by comparison with its successors. If half-tuples of equal value have unequal scores, appending half-tuples of each value in $D_0$ gives at least $d_0$ mixed fibers. Otherwise the half-tuple scores define a function $f$ on $D_0$. If the difference $f(x+t)-f(x)$ takes different values at two points $a_0,a_1$, compare each eligible $x$ with an $a_i$ having unequal difference. The fiber at $x+t+a_i$ is mixed, and at most two choices of $x$ give any one fiber. This yields at least $q_0/2$ mixed fibers. If every difference is constant, a triple intersection of translates of $D_0$ proves that the differences are additive. Varying one slot shows that base labels have the form $kappa log_P A+c_i$ in block $i$. Nonzero $kappa$ violates every homogeneous shift. For a direction satisfying all homogeneous explicit rows, consistency gives $c_i=2sum_j c_j$ for every $i$, and normalization gives $kappa+2sum_j c_j=0$. When $r$ does not divide $2m-1$, these equations force the direction to vanish.

Density alone gives little control over evaluations. Given $l$ points, record each shift's upper-triangular matrix and right-side contributions as a vector of dimension $d_l=l(N(N+1)/2+N)$. Any sequence of more than $(r-1)d_l$ such vectors contains a nonempty zero-sum subsequence. In the group algebra of the additive vector space the augmentation ideal has nilpotency index $(r-1)d_l+1$. The product of the factors $1-[v]$ therefore vanishes, whereas its identity coefficient would be one in the absence of a zero-sum subsequence. Successively removing zero-sum subsequences leaves at most $(r-1)d_l$ rows. The union of the removed subsequences vanishes at all the listed points. When $(r-1)d_l<q_0^2/2$, that union still has the full original shift span. This obstruction concerns subfamilies of rows. The elliptic matrix sums every edge between neighboring fibers.

== Rank after conditioning on a source and successor

Write $G_T$ for the unweighted shift Gram matrix at displacement $T$, let $n_i$ be the block sizes and $M=product_i n_i^2$. Summing over all $T$ counts each ordered pair of independent tuples once. On the factor-base columns this gives
$ sum_T G_T = 4 M^2 ⨁_i (I/n_i - bold(1)bold(1)^t/n_i^2). $
For odd $r$ and $0<n_i<r$ this has rank $N_0-m$. Prefix counts can select at most $N_0-m$ displacements with that joint rank in the stated budget. Their logarithms, which enter the right sides, are unknown. For a prescribed displacement there is another obstruction: the weighted counts are coefficients of a product in $bb(F)_r[U] slash ((U-1)^r)$. A slot of size $n$ has vanishing order at one at most $n-1$ by binomial Vandermonde interpolation. In the guarded regime the total order is at most $4(N_0-m)<r$, even with two pins. The polynomial product survives, yet its coefficient at $P$ can vanish.

*Proposition 3 (conditioned rank).* Under the short-list parameter conditions, partition the unweighted shift sum at $P$ by the source group value and the value in one specified slot of a largest block on each side. The resulting matrices have joint rank $N_0-m$ before the explicit rows are added. The same rank holds if each source tuple is kept separately and all its successors are summed.

*Proof.* Write $d(g)$ for the half-tuple count and $cal(E)=sum_g d(g)^2$. If $d_(-i)$ omits block $i$, the diagonal terms in the square of the sum of its translates give $sum_g d_(-i)(g)^2<=cal(E)/n_i$. Cauchy–Schwarz bounds each full-tuple count pinned in this block by $cal(E)/sqrt(n_i)$. The Stage I potential gives $cal(E)<321r$: at guarded $m>=46$, $r^(1/m)>2^20 m$, and the Stage I expected half-tuple mass satisfies $mu<=16r exp(3m/(16r)^(1/m))<17r$. A largest block has $n_i>(7r/8)^(1/m)>321^2$. Every nonempty pinned class thus has invertible size modulo $r$. At least one pin occurs in each fiber. Twice its count is a nonzero coordinate of that fiber's first incidence moment $v_g$, even if the total fiber size vanishes in the field.

For adjacent pinned classes $B_1,B_2$, let $z_j,V_j$ denote their mean incidence vector and covariance. Divide the class-pair Gram sum by the product of the class sizes to obtain
$ V_1+V_2+(z_2-z_1)(z_2-z_1)^t. $
The log table annihilates each covariance on the left and takes value one on $z_2-z_1$. For a common kernel direction $u$, the two mean scores are equal and $V_1u=-V_2u$. All classes are connected, and choosing one class at each group value produces an odd cycle. All covariance vectors therefore vanish and all class means have the same score. The unconditional covariance also annihilates $u$. Independence of slots gives this covariance as $2 ⨁_i (I/n_i-bold(1)bold(1)^t/n_i^2)$. Its kernel consists of the block constants and free distinguished columns, which already annihilate every shift.

For the source-tuple partition let $K_X$ be the sum over successors of $X$ at $g=g(X)+P$, with second incidence moment $H_g$. Direct multiplication gives $(I+c_X lambda^t)K_X=H_g-v_g c_X^t$. Because $v_g!=0$, every two sources at $g-P$ have the same score under a common kernel direction. Calling the score $f(g-P)$ gives $H_g u=v_g f(g)=v_g f(g-P)$. Thus the score is constant. Varying one slot gives the same kernel as above. $square$

Evaluating one such class fits the closed-count budget. Keeping every source group value requires at least $r$ matrices. The short list combines these classes with linked weights, so Proposition 3 leaves its rank unresolved.

= Finite checks and the uniform question <evidence>

The cyclic group of order five gives 51 ordered shifts on six columns. Exactly 307 of 625 field points yield full rank, and the all-ones evaluation has rank five. On a separate elliptic fixture over $bb(F)_13$ of order 19, direct group addition gives 379 shifts and each of the 17 nonsingular univariate evaluations recovers all eight columns. The single-candidate check extracts a violated row at each of 490 nonzero evaluations. On the order-two fallback, the unguarded zero and one points have joint rank five on six columns. The additional point $(1,0,...,0)$ detects the missing direction. These exact counts test algebraic identities on the specified inputs.

Two examples locate different sources of rank loss. In the cyclic group of order eleven, take two blocks containing all nonzero elements. Both the rule representing $5x$ by four copies of $4x$ and the rule representing $x$ by four copies of $x/4$ admit separate rescalings of the two orbits of multiplication by four. A normalized incorrect label table can therefore pass these consistency checks. In the order-19 neighboring-fiber example, retaining pinned classes gives joint shift rank eight, while summing the classes at all ones gives rank seven. This example lies outside the short-list guard.

The sparse-target experiments ask a related question with a different matrix. For one target $t$ form $A_t$ with one incidence row per ordered tuple of sum $t$. Row $(i,x)$ of $M_t$ sums the incidence vectors with the first slot in block $i$ pinned at $x$. Swapping the two occurrences of each block gives $A_t^t A_t=2M_t$. Consecutive multiples of the prescribed generator can all lose rank in this sum (@pinned). In the first three examples, the vertical stack recovers the forced nullity $m-1$ from block contrasts. Two six-target runs have an additional common direction. In the order-271 example, every pin count is nonzero modulo 271 at every target in the run. The shared direction comes from two translated copies of a two-point block. The original parameter floor forces every block size to grow, so these examples leave the asymptotic question open.

#figure(
  table(columns: (auto, auto, auto, auto, auto, auto), align: center, stroke: none, inset: (x: 7pt, y: 4pt),
    table.hline(stroke: 0.7pt),
    table.header([Order $r$], [Blocks $m$], [Targets $k P$], [Rank of $A_t$], [Rank of $M_t$], [Stack nullity]),
    table.hline(stroke: 0.4pt),
    [223], [5], [$0<=k<=3$], [31], [30], [4],
    [67], [4], [$0<=k<=4$], [28], [27], [3],
    [109], [5], [$0<=k<=4$], [25], [24], [4],
    table.hline(stroke: 0.25pt),
    [67], [4], [$0<=k<=5$], [28], [27], [3],
    [109], [6], [$0<=k<=5$], [24], [23], [6],
    [271], [7], [$0<=k<=5$], [27], [26], [7],
    table.hline(stroke: 0.7pt)),
  caption: [Ranks for six reduced-floor affine seeds. Each row uses its own generator and seed. Each rank holds at every target in the indicated run. The last column uses the vertical stack of the pinned matrices.],
) <pinned>

*Exact moments.* For abstract integer row populations, the degrees required to distinguish modular ranks and to recover every row multiplicity are different. For example, take rows in $ZZ^n$ with entries in $\{-1,0,1\}$, coordinate sum zero and sum of positive entries at most $w$. If $n>=3w$ and $w>=1$, rational moments through total degree $3w$ determine all multiplicities, and this degree is sharp. The #artifact("research/exact_moment_oracle_barriers/moment_rank_boundary/balanced_incidence_moments/RESEARCH.md")[moment supplement] proves the bound by interpolation and alternating examples. Applying it to an implicit elliptic family would also require computing those moments within budget.

*Alternative reductions.* Shorter fixed-shift relations use one slot per block on each side. Conditional-mean equations remove unknown displacement logarithms, but selecting enough displacements within budget remains open. A second reduction maps the residual rank problem to additive avoidance in $bb(F)_r^2$. A constant block of common translates gives full group coverage with tuple product below $r^2$ under the closed contraction bound. If the product of the planar list-size caps is at most $r^(1+sigma)$ for fixed $sigma<1$, there exists a universal probe list of size linear in the sum of the caps. Each admissible sumset misses some probe.

There are constructive results for restricted inputs. Two lists of size at most $n$ admit an explicit probe set of size $n^(1+o(1))$ as $n$ tends to infinity with $n=r^(o(1))$. A carry-aware digit lift gives a related sparse-span problem. Hidden bounded-rank matrix dictionaries and dictionaries with a supplied tree of small linear interfaces admit structured solutions. The general growing-arity setting still requires either an efficient uniform probe construction or a different argument for full rank. A short stack of evaluations in the original certificate would also suffice.

*Data and verification records.* The #artifact("checks")[verification programs] and #artifact("evidence")[recorded transcripts] give the finite counts and ranks. The #artifact("research/sparse_target_runs")[sparse-target data] include the affine seeds and an independent enumeration of the pinned entries. The #artifact("research")[research supplement] gives the proofs for the alternative reductions. The #artifact("COMPARATIVE_AUDIT.md")[formal-source audit] identifies the public theorem statements inspected. The complexity estimates assume the closed contraction bound of LZR.

#heading(numbering: none)[Appendix. A smaller partial expansion]

This appendix concerns an alternative to the weighted Gram construction. In the terminal-preserving contraction of LZR, a partial expansion expresses the network as a sum of products of tensors whose terminal assignments can be fixed later. Write $S=Theta(m b)$ and $ell_0=ceil(log_2 S)$. In LZR Lemma 5.3 the total internal mass at depth $j$ is at most $2^j S$. With $u_j asymp ell_0/2^j$ the cost of distributing siblings adds up to $O(S 2^(3D_0/2)/sqrt(ell_0))$ for $D_0=ceil(log_32(4m))$. The extra power of $m$ comes from enumerating complete products.

*Proposition 4 (root expansion).* Take $m=Theta((log b)^(1/4))$ and $beta=Theta(b)$ terminal bits on distinct hosts. Under the explicit Appendix D construction, the set $cal(I)$ of expansion terms satisfies the bound below. The expansion has $O(m)$ stopped tensors per term, with
$ log_2 abs(cal(I)) = O(S sqrt(log(2m)/ell_0)+ell_0). $
Each tensor has $O(b/m)$ terminal bits and a pinned evaluation costs $2^(O(S/sqrt(ell_0)))$ times a polynomial in the description length. At the original arity the expansion loses a factor $sqrt(log log b)$ over the target exponent.

*Proof.* Use $R$ groups at the root only, taking the smallest power of two above or equal to $max(64,4m)$. Set $u=floor(frac(epsilon ell_0, log_2 R))$, $a=2u$, $W=2^a$, with $epsilon=2^(-17)$. At this arity $u>=1$ for sufficiently large $b$. Handle the finite initial range directly. The square of the pair-label array size is $R^(4a)<=2^(8epsilon ell_0)<=(2S)^(8epsilon)$. The packing guard then holds at an active root with $M>S/sqrt(ell_0)$ internal edges and maximum degree at most $ceil(S^(9/10))$. Conditional expectation balances the squared edge-class deviations and the terminal-bin deviations. Distinct hosts give at most $beta/R+sqrt((a+2)beta)$ terminals per bin. This is $O(beta/m)$ because $R>=4m$ and $m^2 (a+2)=O(beta)$ at the stated arity. The padded common class size is at most $3M/R^(2a)$, hence at most $3M/W$ packets and at most $6M/R$ real child edges per group.

Within any packet and any active coordinate there are still only two possible symbols. At a Hermite node of valuation $i$, $D_t=W/2^i$ and there are at most $6M/D_t$ crossing modes per weight. At the active root $P_0=O(S)=O(M sqrt(ell_0))$. At this arity $2^u/u$ grows faster than any fixed power of $ell_0$. Thus $D_t>=2^u$ gives $p_t=P_0+8ceil(3M/u)D_t+8M(a+1)=O(M D_t/u)$. The weighted partition estimate gives the claimed log cardinality after one step. Keep all $R$ children as canonical stopped tensors. Pinned evaluation continues with 64 groups and the original scalar threshold, starting at absolute depth one with $a_j=a+2sum_(i=1)^j floor(epsilon ell_0/2^i)$. The inherited input precision is at most $S ell_0^(O(1))2^(a_(j-1))$. The next tile scale is at least $2^(a_(j-1)+u_j)$, so the mode bound persists. The per-path exponents sum geometrically. The creation-precision representatives and the original Hermite guard suffice. Early small calls can be directly summed. The #artifact("PARTIAL_EXPANSION.md")[partial-expansion supplement] gives the inherited precision, terminal reconstruction and digit bounds.

The residual graph retains one lookup of degree $Theta(b/m)$ per slot. Contracting that graph gives the uniform deterministic bound $2^(O(b/sqrt(log log b)))$ under the Appendix D assumptions. Reducing the number of expansion terms alone leaves this bound unchanged.

#show bibliography: set text(size: 10pt)
#bibliography("references.bib", style: "ieee", title: [References])
