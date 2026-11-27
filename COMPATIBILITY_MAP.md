# Interfaces and the weighted section norm

| Producer | Consumer | Contract and transfer |
|---|---|---|
| LZR Stage I | Ordinary group-law core | Enumerated, parity-specific field-affine blocks, possibly overlapping. Canonical coordinates and block-specific injective addresses. All targets have a $2m$-slot decomposition. Deterministic. |
| v4 sampler | v4 normalized section | Two signs, disjoint raw masks avoiding the target. Thin lifts are scattered across Fourier owners. Random masks and inherited A.1 estimates give an efficient existence query. |
| Ordinary closed core | Gram evaluation | At most two labeled point pins and public geometric address weights. Digit splitting transfers these weights to closed dyadic scalar contractions with exact integer reconstruction. Deterministic at a given point. |
| LZR open core | Arbitrary candidate weights | Stop with canonical entries at their inherited precisions and attach lookups. Digit reconstruction costs a residual contraction in addition to the partial expansion. |
| Gram evaluation | Linear solve | If nonsingular, the field solution is the entire true table; if singular, reject the trial. One random vector for all rows; a worst-case deterministic choice is open. |
| v4 pilot and sparse solve | Masked targets | A consistent sparse solution suffices when fresh rows lie in the trained span and exact answer checks succeed. This does not certify all coordinates on all targets. |

## Exact geometric transfer, including multiplicities

For the large prime-order curve case let a fixed shift pair comprise $Y$ and $X$. Form the effective divisor with points from $Y$, the negatives of the points in $X$, and $-P$, retaining repetitions. Its degree is $d=4m+1$. The group shift condition is equivalent to this divisor minus $d[O]$ being principal. An element of $L(d[O])$ vanishing to the indicated multiplicity at every point exists exactly in that case. A nonzero vanishing section has pole order exactly $d$; otherwise its zero divisor could not contain $d$ affine points. Normalization of the coefficient of the pole-$d$ monomial is unique, by the one-dimensionality of the space for a principal degree-zero divisor. This is over $\mathbb F_p$ by the rational divisor criterion, or equivalently by the rank of the rational linear conditions. For nonzero ordinates $x-x_0$ is a local parameter and the coefficients of the branch $y(x)$ are solved successively using the unit $2y_0$. Prescribing its first $j$ section coefficients enforces multiplicity $j$, including tangent constraints. Opposite points are different local expansions. There are no nonidentity rational two-torsion points for odd prime group order. Small groups and $P=O$ outside the nontrivial generator case use the direct branch. This supplies a multiplicity-aware field predicate; the v4 simple-root compiler by itself cannot enforce these jets.

Let $\mathcal A$ be the ordered address grid, use an idempotent $I(a)\in\{0,1\}$ for this geometric condition and valid membership, and take $R=(\mathbb Z/2^{P_0})[z]/(z^\ell)$. Put $W(a)=\prod_s(1+z)^{\sigma_s y_{s,a_s}}$ for valid coordinates and zero on padded coordinates. Over the split finite free algebra $R^{\mathcal A}$,

$$
       \sum_a I(a)W(a)
       =\operatorname{Tr}_{R^{\mathcal A}/R}(IW)
       =[t]\operatorname{Norm}_{R^{\mathcal A}/R}(1+t IW).
$$

Indeed multiplication is diagonal in the address idempotents, and the coefficient of $t$ in its determinant is the sum of its diagonal entries. This holds with nilpotents and negative exponents of the unit $1+z$. Apply the Newton functional and the LZR precision cap to recover the exact nonnegative integer. The idempotents are canonical zero/one decisions, so this argument does not presume a unital homomorphism $\mathbb F_p\to R$. A unique normalized field section does not make every *lifted* certificate unique.

The trace is one formal norm coefficient, but the split module has rank $\prod_s|\mathcal A_s|$ (or $2^\beta$ padded). Direct determinant or trace arithmetic on it is outside the desired bound. For a product expansion into $H$ explicitly enumerable pin rectangles with at most $j$ pinned slots, ordinary pinned counts transfer termwise by linearity. More generally arbitrary weights on $j$ slots can be expanded into at most $2^{jk}$ pins per summand, with their actual ring cost. Taking $j=4m$ loses the desired saving. MAIN_RESULT.md achieves $j\le2$ for incidence weights and separates the remaining geometric factors across short address pieces. Its matrix certificate replaces the positive integer separator, so it has different output guarantees.

## Lookup obstruction

For an arbitrary label table on $k=2a$ bits, set $y(u,v)=1[u=v]$ with two $a$-bit halves. Modulo $(2,z^2)$ the positive lookup is $1+zy$. Across this cut the $z$ coefficient is the identity of order $2^a$. A network factoring through a single bond with $D$ states gives a sum of $D$ products. The $z$ coefficient has rank at most $2D$ over $\mathbb F_2$, so at least $2^{a-1}$ states must cross, or at least $a-1$ Boolean wires in a multiwire cut. Entry tables and enumeration of these states must be charged. This example assumes all addresses legal; it is a statement about a universal lookup gadget.

A second counting obstruction allows *arbitrary* topology rather than a single cut. With $O(b)$ Boolean edges, delete isolated scalar vertices, combining their product into one scalar. Labeled bounded-size wiring and port order use $O(b\log b)$ bits. Over a finite ring with $\mathrm{poly}(b)$ bits per entry, the truth tables of all degree-$d$ vertices together use at most $\mathrm{poly}(b)2^d$ bits, even when entry programs have huge descriptions. Thus all such gadgets represent at most $2^{O(b\log b)+\mathrm{poly}(b)2^d}$ different functions. Already binary labels on a legal set of $T$ addresses yield $2^T$ distinct first coefficients. When $T=2^{\Theta(b/m)}$ and $d\le b\,2^{-c\sqrt{\log b}}$, the former count is strictly smaller. This excludes a single universal exact low-degree, linear-edge replacement over this ring. It does not bound representations restricted to actual surviving candidates, sums of many graphs, or a joint contraction that never builds the lookup.

A mere count with fixed masks loses label information even for a single tuple. Interpolating arbitrary labels as a univariate field function on $T$ distinct abscissas in general needs degree at least $T-1$ by dimension. The preprocessing to store a table can fit a trial bound; a joint evaluator on the Cartesian grid remains a different cost.

## Candidate separation and extraction with one checked specialization

For any candidate $y$ passing the explicit equations, form the scalar polynomial

$$F_y(X)=\sum_{\text{valid shifts }i}(a_i^ty-1)M_{v_i}(X)\quad\text{over }\mathbb F_r.$$

Its distinct ordered addresses give distinct exponent vectors, so it is identically zero exactly when every fixed-shift row passes. This statement separates an arbitrary wrong table symbolically. Its total degree is at most $s(D-1)$, independent of the number of rows. At a given point expand $a_i^ty$ as a sum over one pinned slot and one column, multiply by the public $y_j$ outside the contraction, and subtract the unpinned count. The closed-core digit method in MAIN_RESULT.md applies even under an existing prefix. All terms take time polynomial in $N,m,b$ times one count. If the scalar value is nonzero, the two values after pinning the next coordinate bit sum to that value. Select a nonzero child. After $O(mb)$ pins the remaining parsimonious tuple violates the congruence. No further randomness or faith in a zero test is needed for extraction. A random trial misses a fixed wrong table with probability at most $s(D-1)/D^2$. A zero specialization alone cannot certify correctness. The matrix construction makes one specialization work for all candidates together when its rank is full. This supplies an alternative to integer Newton filtering under precisely stated conditions, without an efficiently evaluated arbitrary weighted elliptic norm.
