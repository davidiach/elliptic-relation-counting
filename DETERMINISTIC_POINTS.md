# Deterministic points for the fixed-shift matrix

The matrix, digit map, row keys and field are exactly those of MAIN_RESULT.md. The proposed list can be constructed and its rank tested within the target bit budget. On the guarded branch its stack contains all explicit rows and detects every direction admitting a sparse consequence within the stated degree bound. Only the closed contraction theorem is used for the computational bounds.

[The neighboring-fiber analysis](NEIGHBORING_FIBERS.md) proves full rank on the prescribed full-grid branch, including the order-two exception to the polynomial degree argument. The two displayed point families lose rank on an eligible order-two curve if used without their guard. At the prescribed generator, separately retaining the source tuple or the source group value and two slot pins preserves full shift rank. Whether the coupled evaluations on the guarded short list always have full column rank remains open.

## 1. Interpolating the explicit rows inside the actual stack

Write $E$ for the $C=N+1$ explicit rows, $a_e^\top$ for row $e$, $w=2m$ and $D=2^{\max(1,\lfloor b/4\rfloor)}$. Use $s$ digit variables as in MAIN_RESULT.md. For $z=0,\ldots,C-1$ take the points

$$x(z)=(z,z^D,\ldots,z^{D^{s-1}})\quad\text{in }\mathbb F_r^s.$$

**Proposition.** If $C\le r$ and $r\nmid w-1$, the row space of the stacked $B(x(z))$ contains the entire row space of $E$. Both parameter conditions hold at the asymptotic tuning. In particular the construction gives a full-rank stack whenever $E$ has rank $N$.

Let $\mathbf1$ be the all-ones vector. Every shift row has sum zero. The normalization row has sum one and every consistency row has sum $1-w$. Hence

$$\mathbf1^\top B(x(z))=\sum_{e=0}^{C-1}z^e(\mathbf1^\top a_e)a_e^\top.$$

Inverting the ordinary $C$-point Vandermonde matrix recovers every scaled explicit row. This argument works also when the private address digits of $C$ have carries, since $M_CM_v$ pulls back to $z^{C+v}$. Each projected vector is already in the literal stack. No row is appended by hand. For the finitely many small lengths where the test on parameters fails, this proposition gives no interpolation guarantee.

## 2. A list for sparse polynomial consequences of the kernel equations

Put $K=\max(2,C)$ and $H=2^{\lfloor b/3\rfloor}$. For a nonzero column direction $u$ let $I_u$ be the ideal generated in $\mathbb F_r[X_0,\ldots,X_{s-1}]$ by the coordinates of $B(X)u$. Membership here means a polynomial linear combination, with arbitrary polynomial multipliers. They are used only in the proof. The list requires neither their coefficients nor the direction $u$.

**Theorem.** Construct a primitive element $\gamma$ of $\mathbb F_r^\times$. Put $J=1+(s-1)\binom K2$. Append to the first list all points

$$x(t,h)=(\gamma^{h},\gamma^{ht},\ldots,\gamma^{h t^{s-1}}),\qquad
1\le t\le J,\quad 0\le h<K.$$

Suppose $H\sum_{j=0}^{s-1}J^j<r-1$. Every surviving nonzero $u$ lies in $\ker E$, and every nonzero polynomial in $I_u$ of individual degrees at most $H$ has **more than $K$ monomials**. In particular the stack has full column rank if each nonzero $u\in\ker E$ admits even one such polynomial with at most $K$ terms. The latter is a sufficient structural condition. We have no general proof that it holds for LZR blocks.

For the proof, let $f$ be a nonzero polynomial with at most $K$ terms satisfying the degree bound. For distinct exponent vectors $v,v'$, the integer polynomial $\sum_j(v_j-v'_j)t^j$ is nonzero and has at most $s-1$ integer roots. One of the $J$ slopes separates every pair of terms of $f$. Each projected exponent is between zero and $H\sum J^j<r-1$. The powers of $\gamma$ at these exponents are distinct. The first $|\operatorname{supp}f|$ successive evaluations in $h$ are a nonsingular Vandermonde matrix applied to its nonzero coefficient vector. At least one evaluation is nonzero. If $B(x)u=0$ throughout the list, every element of $I_u$ vanishes throughout, proving the assertion. This geometric progression test uses the classical sparse-interpolation argument. [Ben-Or and Tiwari](https://doi.org/10.1145/62212.62241)

Dense inputs can have short polynomial consequences. For instance $\prod_{j=1}^2(1+\cdots+X_j^{n-1})$ has $n^2$ terms and multiplication by $(1-X_1)(1-X_2)$ leaves four. Applying a sparse theorem to an arbitrary entry without establishing such a consequence gives no rank guarantee.

**Bit budget and output.** LZR Proposition 3.3, applied with the prime $r$, constructs $\gamma$ deterministically from the closed theorem. Its proof factors $r-1$ via Proposition 2.6, enumerates a generating union from multiplicative Stage I, extracts each prime-power component and multiplies the components. Hasse gives $\log r=b+O(1)$. Since $s=O(1)$ and $\log K=O(b/m+\log b)$, the exponent of the no-wrap bound has logarithm at most $b/3+O(b/m+\log b)<\log(r-1)$ eventually. The explicit list size is at most $C+JK=O(C+sK^3)$. Generation per point costs polynomial bit time in $b$ at this tuning. Evaluate the original matrix and right side by the charged two-pin contractor. Reduce the rows incrementally over $\mathbb F_r$, retaining a point only when the rank increases. At most $N$ retained points attain the tested stack's rank. If it reaches $N$, solve and verify. Otherwise return the attained rank and a nonzero vector in the common kernel of the tested matrices. Total time and space are each at most $2^{O(b/(\log b)^{1/4})}$, including root construction, all evaluations and elimination. For the bounded branch where these parameter tests fail, test the entire field grid. Its stack always has full rank by Proposition 1 in NEIGHBORING_FIBERS.md, including the order-two case.

## 3. How many fixed shifts detect a direction?

Here $G=\langle P\rangle$ has prime order $r$. Denote by $D_0$ the sums of tuples with **one** entry in each block and by $d$ its cardinality. The actual Stage I guarantee is $d>7r/8$ in the large-input regime. Each full tuple is a pair of one-per-block tuples. Set $q=2d-r$. Every $g\in G$ has at least $q$ ordered full tuples, since $|D_0\cap(g-D_0)|\ge q$ and we can choose one representative for each half. Let $S$ denote the matrix of all ordered fixed-shift rows.

**Density theorem.** Its kernel on factor-base columns consists exactly of tables constant on each block. Dedicated columns are free. If $u\notin\ker S$, at least $q^2/2>9r^2/32$ different shift rows have $a^\top u\ne0$. If all full tuples of each group value have the same label sum, then the factor-base labels have the form

$$u_{i,A}=\kappa\log_P A+c_i.$$

In particular a nonzero direction in $\ker E$ cannot have this last property when $r\nmid(2m-1)$.

**Proof.** If every shift annihilates $u$, two full tuples of equal value have equal label sum by comparison with a common successor. The resulting function on $G$ is invariant under $+P$, hence constant. Changing one slot gives equality of all labels within that block. The converse follows from the fixed block multiplicities.

For the quantitative claim, call a fiber mixed if two full tuples with its group value have unequal label sums. Each mixed fiber at $g$ gives at least $q$ distinct violations. Fix two witnesses in it and compare both to each successor tuple of value $g+P$. At least one comparison violates, and successor tuples distinguish the comparisons.

If two one-per-block tuples of the same value already have unequal sums, append to both a representative of each element of $D_0$. This gives $d$ distinct mixed fibers. Otherwise label sums of half tuples define $h:D_0\to\mathbb F_r$. Write $D_t=\{x:x,x+t\in D_0\}$, with $|D_t|\ge q$, and $\delta_t(x)=h(x+t)-h(x)$. If some $\delta_t$ takes different values at $a_0,a_1$, each $x\in D_t$ has $\delta_t(x)\ne\delta_t(a_i)$ for at least one $i$. The two half-value pairs $(x+t,a_i)$ and $(x,a_i+t)$ give unequal sums in fiber $x+t+a_i$. Each such fiber has at most two preimages among the chosen $x$, giving at least $q/2$ mixed fibers.

In the remaining case all $\delta_t$ are constant. For $t,v\in G$, choose $x$ with $x,x+t,x+t+v\in D_0$ by $3d-2r>0$. Telescoping proves that the constant differences define a homomorphism $G\to\mathbb F_r$, necessarily $\kappa\log_P$. Thus $h$ is affine on $D_0$. Fix the other half-tuple slots and vary one entry to obtain the displayed block-label form. If $\kappa\ne0$ every shift violates and there are at least $rq^2$ shifts. If $\kappa=0$ every shift annihilates $u$. This proves the bound and the classification. For $u\in\ker E$ the explicit equations on base columns give $c_i=2\sum_j c_j$ for every $i$, whence $c_i=0$ if $r\nmid(2m-1)$. Consistency at $P$ and normalization give $\kappa=0$, and consistency fixes the two dedicated coordinates.

A consequence for the actual polynomial is that at least one coordinate of $B(X)u$ has at least $q^2/(2N)$ shift monomials for each nonzero $u\in\ker E$. An offending row contributes a nonzero entry to $a(a^\top u)$. Different ordered addresses give distinct digit monomials, also distinct from all explicit keys by Kronecker substitution. Density gives a lower bound on the number of witnesses; a modular sum of them may still vanish.

## 4. A limitation of arguments that retain only row coverage and density

Fix **any** $L$ base-field evaluation points on one eligible input, including points chosen using the input. Let $d_L=L(N(N+1)/2+N)$. There exists a subset of the actual shift rows, omitting at most $(r-1)d_L$ of them, whose contributions to **both** $B$ and the right side vanish at every listed point. If $(r-1)d_L<q^2/2$, the retained subset has the entire original shift row span. This inequality holds asymptotically for every $L=2^{O(b/(\log b)^{1/4})}$.

For each shift record its $d_L$ upper-triangular Gram and right-side coordinates at the listed points as a vector in $\mathbb F_r^{d_L}$. Any sequence of more than $(r-1)d_L$ vectors has a nonempty zero-sum subsequence. To see this, write the group algebra as $\mathbb F_r[T_1,\ldots,T_{d_L}]/(T_1^r,\ldots,T_{d_L}^r)$ via $T_j=[e_j]-1$. Every $1-[v]$ is in the positive-degree ideal. A product of more than $(r-1)d_L$ such factors vanishes. Its identity coefficient would be one if no nonempty subset summed to zero. Repeatedly remove zero-sum subsequences until at most $(r-1)d_L$ vectors remain. The union of the removed subsequences is the desired retained row family. If its span were smaller, a direction annihilating it but not the original shifts would have all its violations among those omitted, contradicting the density theorem.

With the explicit rows restored, this modified family retains full symbolic column rank. At the specified points its matrices have rows in the span of $E$. Thus even omitting a vanishing fraction of *valid* rows can defeat a proof based solely on spanning and abundant violations whenever $E$ is deficient. The actual contractor sums **every** valid row. The omitted subset need not have an efficient description, and this argument supplies neither a counterexample to the full elliptic certificate nor a lower bound against input-dependent algorithms on that certificate. The zero-sum bound is the elementary abelian case of the Davenport constant. [Olson, 1969](https://www.sciencedirect.com/science/article/pii/0022314X69900213)

## 5. What remains to prove

The exact structure discarded by the obstruction can be written as a fiber sum. Split each address monomial into $W_-(X)W_+(Y)$ on source and successor tuples. This factors even through a radix digit crossing the cut, because its two exponent contributions add. For each side and group element $g$, let $n_\pm(g),v_\pm(g),H_\pm(g)$ be the weighted sums of $1,c,cc^\top$, respectively, with $c$ the incidence vector of a tuple of value $g$. After the common offset monomial is factored out, the shift part of $B$ equals

$$\sum_{g\in G}\big[n_-(g)H_+(g+P)+n_+(g+P)H_-(g)
 -v_-(g)v_+(g+P)^\top-v_+(g+P)v_-(g)^\top\big].$$

This follows by expanding $(c_Y-c_X)(c_Y-c_X)^\top$ on each complete bipartite pair of fibers. Explicit enumeration of the $r$ values of $g$ exceeds the target. The closed contractor already evaluates the sum at a supplied point. Controlling its cancellation across points is the missing step. Across the same natural cut, the unweighted feasibility table itself contains an $r$-by-$r$ permutation submatrix, by retaining one tuple per group element on each side. A factorization passing through a single linear space at this cut needs dimension at least $r$.

The first interpolation removes all explicit constraints from the unknown part of the problem. Every remaining direction has many mixed-fiber witnesses. To finish via the second list one could use the exact all-pairs convolution of neighboring fibers to exhibit a short polynomial consequence in $I_u$ for each remaining direction. Neither the number of violations nor the total number of generic basis minors proves the existence of that consequence. A different deterministic point construction could also work. The checks in checks/check_fixed_shift_points.py evaluate the original matrices on a small elliptic curve, including an example where two singular evaluations suffice together.

A second order-19 test uses parity-affine masks on the actual curve with threshold three and block sizes four and six. Its half sums cover 18 group elements. The explicit rank is ten and two retained evaluations reach rank twelve. The fiber formula at scalar two was checked against the direct sum over every ordered shift. These are finite examples of indexed blocks, with the second using the mask formula. The LZR small-group fallback is separate.
