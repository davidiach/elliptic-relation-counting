# What the complete neighboring-fiber sum implies for the proposed stack

The full-field-grid branch of DETERMINISTIC_POINTS.md always has full column rank. A genuine prime-order curve loses rank if the two displayed lists are used outside their guard. The same curve has full rank under the prescribed fallback. For the guarded short list, the calculation below proves a rank identity over all displacements and isolates the difficulty with fixing a single generator. It does not decide the uniform short-list question.

## The full-grid branch

Keep the original $b,D,s,C,N$ and the address convention $M_CM_v$. Let $A$ be the matrix of the explicit rows and every ordered shift row. For clarity, the coefficient of $M_CM_v$ in the shift sum is precisely $aa^\top$, where $a=c_Y-c_X$. Indeed, $(X,Y)$ occurs once in the neighboring-fiber sum, at $g=\mathrm{val}(X)$ in the product of the fibers at $g$ and $g+P$. Expanding $(c_Y-c_X)(c_Y-c_X)^\top$ gives its four moment terms. Distinct valid address pairs give distinct exponents and they cannot collide formally with the explicit rows, as substituting $X_j=Z^{D^j}$ gives disjoint integer keys.

**Proposition 1.** For every eligible curve the stack over $\mathbb F_r^s$ has full column rank. If $r>2$, any Cartesian grid with $2D-1$ distinct values in each coordinate already suffices.

For $r>2$, each entry of $B$ has individual degrees at most $2(D-1)<r$. For $3\le b\le7$ use $D=2$. For $b\ge8$ Hasse gives $r\ge(\sqrt p-1)^2>p/2>2^{b-2}\ge2(D-1)$. If $B(x)u=0$ on the specified grid, repeated univariate interpolation makes $B(X)u$ the zero polynomial. The coefficient attached to each row now gives $a(a^\top u)=0$. Thus $Au=0$, including any zero rows harmlessly. The complete shift equations make the label sums equal for tuples of the same value by comparison with a common successor. The homogeneous shift increments are zero, so the sums are constant on the group. The explicit equations at $P$ set this constant to zero and the other consistency rows give $u=0$.

If $r=2$, Hasse and $p>3$ force $p=5$. Then $b=3,m=2,D=2$, both blocks are $G$ and the prescribed decomposition of $A$ is $(A,O,O,O)$. The explicit rows already have full rank over $\mathbb F_2$, as can also be seen in the example below. On the Boolean grid the unique multilinear representatives use $X_j^2=X_j$. Each shift monomial reduces to the bitwise union of the masks $C$ and $v$. Its integer key is at least $C$. The explicit monomials, keyed $0,\ldots,C-1$, therefore retain their distinct outer-product coefficients. If the stack kills $u$, every explicit row kills it. This proof also covers either target and any injective choice of addresses for the two elements.

**Which bit lengths does this settle?** With the published tuning $m=\max(2,\lceil\lceil\log_2(b+2)\rceil^{1/4}\rceil)$ the short-list no-wrap guard necessarily fails whenever $m\le45$. In particular the full-grid proposition covers all $b+2\le 2^{4\,100\,625}$. This is a bounded range of bit lengths. The full-grid cost on it fits the target big-O with an absolute constant, since $s=O(1)$ and $(\log b)^{1/4}\le45$. This observation gives no asymptotic improvement by enumerating grids at unbounded $m$.

To check the cutoff write $d_0$ for the size of the one-entry-per-block sumset, so $d_0\ge7r/8$, and $n_i=|F_i|$. Then $2^\beta\ge(\prod_i n_i)^4\ge d_0^4$, and $C\ge3+m d_0^{1/m}$ by AM–GM. For $b\ge40$, Hasse gives $r>2^{b-2}$, and hence $4\log_2 d_0>15b/4$. As $D\le2^{b/4}$, necessarily $s\ge16$. The short-list number of slopes obeys

$$J=1+(s-1)\binom C2\ \ge\ (s-1)(m^2/2)(7/8)^{2/m}r^{2/m}\ \ge\ (7/4)(s-1)r^{2/m}.$$

Also $H=2^{\lfloor b/3\rfloor}>2^{-4/3}r^{1/3}$, using $r<2^{b+1}$. For $m\le45$, $1/3+2(s-1)/m\ge1$ and the constant from $J^{s-1}$ more than absorbs $2^{-4/3}$. So $H\sum_{j=0}^{s-1}J^j>r-1$. For $b<40$, $m=2$ and $s\ge2$ (for $b\le7$ use $d_0\ge2$, and otherwise use the same Hasse bound). Already $J>r-1$. The proof in DETERMINISTIC_POINTS.md that the guard passes eventually at growing $m$ still applies.

## A boundary example checked at every point

Take $E:y^2=x^3+2x$ over $\mathbb F_5$. It is nonsingular, and its only affine point is $P=(0,0)$, so $G=\{O,P\}$ has prime order two. Either $Q=O$ or $Q=P$ works. Use column order $(O_1,P_1,O_2,P_2,P_*,Q_*)$, normalization key zero, then the six consistency keys in column order. Take the natural one-bit addresses $O=0,P=1$ for the exhaustive calculation. The formulas for the two unguarded lists yield only the all-zero and all-one points.

At all ones each unweighted tuple fiber has eight elements. Each of the two slots of a given block, when fixed to a specified element, occurs four times in either fiber. Thus the mass and first incidence moment of both fibers vanish over $\mathbb F_2$. Every term $n_gH_{g+P}+n_{g+P}H_g-v_gv_{g+P}^\top-v_{g+P}v_g^\top$ is zero. At the origin the offset monomial kills the shifts. Among the explicit rows the consistency rows for $O_1$ and $P_1$ are both $e_{O_1}$, so their outer products cancel at all ones. The five remaining rows are independent: the normalization row fixes $P_*$, the row for $O_2$ fixes $O_2$, the row for $P_*$ then fixes $O_1+P_1$, and the rows for $P_2,Q_*$ fix their own coordinates. Their Gram has rank five because their transpose is injective. The origin supplies only normalization. The two-point stack consequently has rank five, with nonzero kernel vector $(1,1,0,0,0,0)^\top$.

On the prescribed fallback grid the point $(1,0,\ldots,0)$ retains just keys zero and one. It detects that kernel vector. Together with all ones it gives rank six. The finite verifier independently forms the neighboring-fiber and direct ordered-pair matrices at all 256 points for both targets. This example disproves unguarded use of the two formulas. It is not a counterexample on either branch of the proposed piecewise algorithm.

## Summing the exact fibers over displacements

Write $G_T$ for the **shift part alone** at the all-one evaluation for displacement $T\in G$. Keep all original blocks. Let $M=\prod_i n_i^2$ be the number of full ordered tuples. In the LZR general regime $1\le n_i<r$ since the blocks exclude $O$, and $r$ is odd. Then, on the factor-base columns,

$$\boxed{\quad\sum_{T\in G}G_T
=4M^2\bigoplus_i\left(n_i^{-1}I_{n_i}-n_i^{-2}{\bf1}{\bf1}^\top\right).\quad}$$

The two dedicated columns contribute zeros. Summing the exact neighboring fibers counts every ordered tuple pair $(X,Y)$ once. Dividing by $M^2$ amounts to two independent uniform tuples in the field. Their incidence difference has mean zero. A block of one tuple is the sum of two independent categorical observations of covariance $I/n_i-{\bf1}{\bf1}^\top/n_i^2$. Subtracting the other independent tuple doubles this covariance again, proving the formula entirely over $\mathbb F_r$. Each block has kernel exactly its constant vectors. Thus the aggregate has rank $N_0-m$.

**Constructive consequence.** There is a deterministic selection of at most $N_0-m$ displacements whose unweighted shift matrices jointly have this rank, in the target total bit budget. Start with their row space $W=0$. If the displayed aggregate has a row $j$ outside $W$, find $u$ annihilating $W$ with the $j$th aggregate row evaluating to one. For a bit prefix on the canonical elliptic coordinates of $T$, compute the $j$th row of the partial displacement sum applied to $u$. For each column expand the two incidence factors as signed slot pins, as in MAIN_RESULT.md. The arithmetic core computes $T=\mathrm{val}(Y)-\mathrm{val}(X)$ once with the full group law and tests its prefix. Apply coefficients of $u$ outside the graph. This uses $\operatorname{poly}(N,b)$ closed calls on $O(mb)$ edges and linear precision, with exact nonnegative counts reconstructed before field reduction. Starting at one, the two child sums add to the parent. A nonzero branch reaches a displacement $T$ where row $j$ of $G_T$ increases the span. All these rows have zero sum in each block, so the process stops at exactly $N_0-m$. Elimination, $O(b)$ prefix decisions per selected point and all storage fit $2^{O(b/(\log b)^{1/4})}$ bit time and space.

This selects group **displacements**, whereas PR #3 requires evaluation points for the **single prescribed displacement $P$**. An arbitrary selected $T$ has right-side scalar $\log_P T$, which the construction does not give. Appending the original explicit rows restores full homogeneous column rank but does not supply those missing right sides.

**Why a short shift recurrence does not transfer the result.** For any direction $u$ that varies within some factor-base block, the aggregate has a coordinate $j$ with nonzero value on $u$. Put $f(k)=(G_{[k]P}u)_j$ for $k\in\mathbb F_r$. Its cyclic sequence has scalar linear recurrence order exactly $r$ over $\mathbb F_r$. Indeed $F(U)=\sum_{k=0}^{r-1}f(k)U^k$ satisfies $F(1)\ne0$, so it is a unit in $\mathbb F_r[U]/((U-1)^r)=\mathbb F_r[U]/(U^r-1)$. A recurrence valid at every origin would annihilate this unit by a polynomial in translation, and its annihilator is zero in that quotient. Equivalently the $r$ translates of $f$ form a basis. The unique translation-invariant scalar filter giving the nonzero full sum at all origins uses all $r$ offsets, since $1+U+\cdots+U^{r-1}$ does so and multiplication by $F$ is invertible. This rules out one proposed way to move from the aggregate to a fixed shift. It places no lower bound on tests that vary address weights at fixed $P$ or use the elliptic encoding.

## A nonvanishing result at every torus point of the short list

One possible exact cancellation mechanism can be excluded in the guarded regime. For a listed torus point $x$ and any consistent choice of up to two valid slot pins, let $c_x(k)$ count weighted tuples of displacement $[k]P$ modulo $r$. Include the address monomials, omitting the common nonzero public offset. Form the group-algebra element

$$F_x(U)=\sum_{k=0}^{r-1} c_x(k)U^k
=\prod_{\ell=1}^{4m}\left(\sum_{A\in F_{i(\ell)}^{\,\mathrm{pin}}}
 w_\ell(A)U^{\sigma_\ell\log_P A}\right)
\quad\bmod (U-1)^r.$$

Here a pinned slot is a singleton and every $w_\ell(A)$ is nonzero. The equality follows by grouping the very same ordered pairs by their displacement; a radix digit split across slots only adds their exponent contributions. Discrete logarithms in this identity serve as coordinates for the proof.

**Proposition 2.** This product never vanishes as a whole. Its cyclic coefficient sequence has recurrence order at least $r-4(N_0-m)$. Both conclusions hold at every torus point of the guarded list.

For one factor with $n$ allowed elements, reduce the group exponents to distinct integers $d_i\in[0,r-1]$. The first $n$ coefficients after setting $U=1+V$ are obtained from the weights by the matrix $\binom{d_i}{j}$, $0\le j<n$. Its determinant up to sign is $\prod_{i<k}(d_k-d_i)/\prod_{j=0}^{n-1}j!$, nonzero modulo $r$. Thus this factor's order of vanishing at one is at most $n-1$. The orders add before reduction, and their sum $\nu$ is at most $4\sum_i(n_i-1)=4(N_0-m)$. In the short-list regime $s\ge2$ and $r>J\ge\binom C2>4(N_0-m)$, using $C=N_0+3$. Reduction modulo $V^r$ therefore preserves the product. Its annihilator is generated by $V^{r-\nu}$, giving recurrence order exactly $r-\nu$. Pins only decrease the bound on $\nu$.

This excludes a counterexample in which nonzero torus weights annihilate an entire elementary displacement convolution. The specific coefficient $[U^1]F_x$, which is the one required at $P$, can still vanish. Distinct pinned terms can also cancel when assembled into a Gram entry. No inference from whole-product nonvanishing to either of these more specific statements is justified.

The verifier illustrates the distinction on $y^2=x^3+20x+1$ over $\mathbb F_{251}$ with generator $(0,1)$ of order 257. The two arithmetic masks use $x$ and $x+23$ modulo 251, threshold 23, and the parity offset $23(y\bmod 2)$. Their sizes are 24 and 28, and the half sums cover 242 group elements. With six-bit addresses, radix four, $\gamma=3$ and $(t,h)=(4,48)$, the elementary sum at displacement $P$ is zero while 254 displacement coefficients are nonzero. All four slots on each side contribute, and the test uses their exact fiber distributions. These masks exercise the algebraic identity; this input size falls in LZR's small-group regime and the masks are not its prescribed fallback blocks. This is a control for a single coefficient, not a stack counterexample.

For the guarded regime the unsolved statement is still $\bigcap_{x\in\mathcal L}\ker B_P(x)=0$ for the actual prescribed list $\mathcal L$. The full-grid theorem uses points outside that list, and the displacement construction changes $P$. Neither step resolves that statement.
