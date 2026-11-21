# Exact product equality by a balanced chord calculation

For a regular sampled grid this algorithm decides whether two lifted determinant products agree modulo $p^2$ without constructing either one. It handles empty, singleton and multiple target fibers. For a single address-bit speed, inequality on a singleton recovers the bit. The complexity is essentially that of building the two balanced lists, even in the empty case where ordinary group matching alone fails to decide product equality.

## Reduce the modulus before computing the answer

Use a smooth short Weierstrass curve over $\mathbb F_p$, $p>3$. Include the endpoint $-T$ as a final slot of size one. Regularity means that distinct slots have disjoint abscissas, including the endpoint. Write $d=m+1$ for the number of rows, $\Delta$ for their evaluation determinant in $L(dO)$, and $D$ for the invariant derivation $D(x)=2y,\ D(y)=3x^2+A$. Let $s_i(P)$ be the **difference** of the two tangent speeds at the evaluated point. Thus moving the target rather than the endpoint reverses that speed. The two calls have the same initial lifts. Put $k(a)=\sum_i s_i(P_i)$ and let $h$ count the tuples whose row points sum to $O$. If any slot is empty both products are one.

The first-order product rule gives three cases. For $h\ge2$ the two products vanish modulo $p^2$. For $h=1$ the products agree exactly when $k(a_*)=0\bmod p$, because the common invariant derivative of the one vanishing factor is nonzero and all its other factors are units. For $h=0$, divide the product difference by $p$ and by the nonzero residue of the baseline product *in the proof*. Equality is equivalent to

$$
L=\sum_a\sum_i s_i(P_i)\,\frac{D_i\Delta(a)}{\Delta(a)}=0\quad\text{in }\mathbb F_p.
$$

No product or lifted scalar is needed by the algorithm. An empty fiber can give either answer. [Precision-two identity and regularity](../direct_parity/deformation.md)

## A finite-field formula for the logarithmic derivative

For affine $P,Q$ of distinct abscissas let $l(P,Q)=(y_P-y_Q)/(x_P-x_Q)$. For $P=Q$ with $y_P\ne0$ use the tangent slope, and set $l=0$ if either input is $O$ or if $Q=-P$. In a binary addition tree with row points as leaves, let $C(a)$ be the sum of $l(S_1,S_2)$ at all internal nodes whose two children have group sums $S_1,S_2$. For original row points define $K(P,Q)=(y_P+y_Q)/(x_P-x_Q)$. For every nonsingular tuple determinant on a regular grid,

$$
\frac{D_i\Delta}{\Delta}=C(a)+\sum_{j\ne i} K(P_i,P_j). \tag{1}
$$

Here is an algebraic proof which also handles exceptional partial sums. For two nonidentity points that are not opposites, use the chord or tangent line divided by the vertical line through their sum. This rational function has divisor $[P]+[Q]-[P+Q]-[O]$. For opposites use just the vertical line, and for an identity input use the constant one. In the parameter $t=-x/y$ at $O$ the ordinary line quotient is $-t^{-1}(1+l(P,Q)t+O(t^2))$. The vertical line is $t^{-2}(1+O(t^2))$. Thus the first coefficient after its leading monomial, normalized by the leading coefficient, is always the indicated slope. Products over trees for the same list have the same divisor, hence are constant multiples of each other. Their normalized first coefficients agree, so $C$ is independent of the tree.

Fix all rows except $i$, with sum $S$ of those fixed rows. The determinant as a function of row $i$ is a scalar multiple of the section with affine zeros at those fixed points and at $-S$ when $S\ne O$, counted with multiplicity. This follows from the degree-one residual divisor in the elliptic linear system. If $S=O$ the residual zero is at infinity and the affine function has a pole of order $d-1$ there. The tree for the fixed points followed by the residual gives the same divisor as this function and has normalized first coefficient $C(\text{fixed rows})$, since $l(S,-S)=0$.

For each finite zero $Z$ of that function, $K(X,Z)\,dx/(2y)$ has residue one at $Z$ and minus one at $O$, with no other pole. For a point of order two this follows using $y$ as local parameter. At $-Z\ne Z$ the apparent singularity is removable. Subtract the sum of the $K(X,Z)$ from the logarithmic derivative $Df/f$. The difference has no poles and is constant: the local logarithmic derivative has at most a simple pole, also at repeated zeros. At $O$, $D=(1+O(t^4))\partial_t$, and $K(X,Z)=-t^{-1}+O(t)$. The constant is exactly $C(\text{fixed rows})$. Finally the residual term $K(P_i,-S)$ equals $l(P_i,S)$ with the tangent convention, or is omitted when $S=O$. Since $\Delta\ne0$, $P_i=-S$ cannot occur. Adding $P_i$ to the tree now gives (1). The formula is related to the classical elliptic alternant identity; this derivation takes place in the base field. [Frobenius–Stickelberger formula](https://arxiv.org/abs/2603.27466)

## Sum the chords without enumerating cross pairs

Antisymmetry of $K$ in (1) gives

$$
L=\sum_a k(a)C(a)+
\sum_{i<j}\left(\prod_{s\ne i,j}|F_s|\right)
 \sum_{P\in F_i,Q\in F_j}(s_i(P)-s_j(Q))K(P,Q).
$$

Split the varying slots into balanced halves, attaching the endpoint to one side. For each half tuple store its group sum, its speed sum and its internal chord sum. Write the two tables as $(P,u,c)$ and $(Q,v,e)$. An ordinary join $P=-Q$ gives the first two cases above by accumulating hit counts and speeds. In the empty case $kC=(u+v)(c+e+l(P,Q))$. The contribution $(u+v)(c+e)$ factors into separate sums over the two tables. It remains to evaluate the single weighted cross sum $\sum_{P,Q}(u+v)l(P,Q)$.

After grouping right points by abscissa, compute four Cauchy transforms $\sum_z w_z/(X-z)$ for the right-side weights $1,y_Q,v,v y_Q$, aggregated at each $z=x_Q$. Omit identity points in the cross sum. At a left abscissa $x$ evaluate the four transforms with the diagonal $z=x$ excluded. Denote their values by $F,F_y,G,G_y$. The off-diagonal contribution of $(P,u,c)$ is $u(y_P F-F_y)+y_P G-G_y$. At equal abscissas only actual doublings with $y_P\ne0$ contribute; their exact tangent slopes can be added from a dictionary keyed by the signed point. The regular row-pair terms involving $K(P,Q)=l(P,-Q)$ use the same computation with the right point and weight negated.

For clarity, write a transform as $A(X)/B(X)$ where $B=\prod_z(X-z)$ over distinct abscissas. At an abscissa shared with the left list the finite part is

$$
\frac{A'(x)}{B'(x)}-w_x\frac{B''(x)}{2B'(x)}.
$$

All other abscissas use the ordinary quotient. A product tree constructs the rational function by pairwise addition; a remainder tree evaluates the needed polynomials at all left abscissas. Thus the cross sum takes quasi-linear field work in the two list lengths. Repeated group points, equal abscissas, two-torsion and internal identity sums do not require expanding colliding pairs. The univariate products in these evaluation trees contain abscissa labels; neither original determinant product is formed. [Fast multipoint evaluation](https://cs.uwaterloo.ca/~eschost/publications/evalint.pdf)

With $w$ the common raw width, bit time and space are bounded by $(2w)^{\lceil m/2\rceil}\operatorname{poly}(m,b)$ on nonempty regular grids. Slot expansion and square roots over the base field are included in expected bit time. At the project's balanced widths this is log cost at most $b/2+O(b/m+m+\log b)$. The method improves on literal full-grid product evaluation. The isolated bit alone was already obtainable at this half-list scale by a group join. The new calculation covers equality on the empty fibers too. Typical successful inputs still have $2^{b/2-o(b)}$ half entries, so this algorithm does not beat the project's conditional general-prime exponent. [Sizes on the actual sampler](../direct_parity/scalar_join.md)

Once the geometry of a seed and target is fixed the computed zero test in each case is a linear form in the input speeds. In the empty case the arithmetic for this form is linear, because the denominators, subtree counts and all chord sums are independent of the speeds. Transposing that calculation (reverse accumulation from one output) gives every coefficient at the same asymptotic setup cost. Subsequent equality queries need just a dot product. On a singleton the normalized form is supported on the witness, with one coefficient equal to one in each slot. Reusing the form does not remove the cost of a fresh seed and target.

## Reading the bit as a residue at the target

More generally, vary the endpoint $Q$ while leaving its speed zero and keeping the slot speeds fixed. Write $S(a)$ for the sum of just the slot points and $k(a)$ for their speed sum. Away from row collisions and target relations, (1) gives the normalized product difference as a rational function

$$
L_c(Q)=C_c+\sum_{a:S(a)\ne O} k(a) K(Q,-S(a))
                 -\sum_{a,i}c_i(P_i(a))K(Q,P_i(a)).
$$

Here $C_c$ is independent of $Q$: it collects the internal slot slopes and row-pair terms in the tree with the endpoint added last. The last chord is $l(S,Q)=K(Q,-S)$, continued to the tangent where appropriate. This identity is an identity of rational functions in $Q$. As proved above, $K(Q,Z)\omega_Q$ has residue one at an affine $Z$ and no other finite pole. Thus at a regular endpoint $Q=-T$ the residue of $L_c(Q)\omega_Q$ is **exactly** $\sum_{a:S(a)=T}k(a)$. On a regular singleton it recovers any bit, including the ordinate sign, or the whole packed speed sum. The additive constant and every row-pair term drop out of this readout. The unresolved operation is computing one such residue from the individual seed descriptions without representing the distribution of sums.

In ordinary abscissa coordinates set $A_x=\sum_{a:x(S(a))=x}k(a)$ and $B_x=\sum_{a:x(S(a))=x}k(a)y_{S(a)}$, excluding the identity. For $2T\ne O$ the one residue is $(A_{x_T}+B_{x_T}/y_T)/2$: each tuple with sum $T$ contributes $k(a)$ and each with sum $-T$ contributes zero. This separates the two target signs even for speeds that depend on the ordinate sign of a slot point. It suffices to obtain the two entries at one abscissa.

Take only one moving slot $i$ and set its speed to a chosen raw abscissa bit. Write $\mathcal A_1$ for the full tuples for which this bit is one, omitting the endpoint from each tuple, and $S(a)$ for their group sum. Both the original slots and $\mathcal A_1$ are invariant under simultaneous negation. Let $E=-T$, assume $2T\ne O$, and define over $\mathbb F_p$

$$
R_i(X)=\sum_{a\in\mathcal A_1,\ S(a)\ne O}\frac1{X-x(S(a))}
      -\sum_{a\in\mathcal A_1}\frac1{X-x(P_i(a))}.
$$

On an empty regular grid $L=y_E R_i(x_E)$. To see this, put the endpoint last in the chord tree. Under simultaneous negation of the varying points all their internal slopes and their mutual $K$ terms change sign and cancel in the sum. The average of $l(S,E)$ and $l(-S,E)$ is $y_E/(x_E-x(S))$ for $S\ne O$. Neither $S=E$ nor $S=-E$ occurs on a symmetric empty grid. The average of $K(P_i,E)$ and $K(-P_i,E)$ gives the negative slot term after division by $y_E$. This also covers rational two-torsion among the variable points.

On a regular singleton the coefficient of $(X-x_E)^{-1}$ in $R_i$ is exactly twice the witness bit: a selected sum $T$ has exactly one negated partner of sum $-T$, the targets are distinct, and no slot point has abscissa $x_E$. Hence the **absence of a pole** tests equality on a singleton even if the value there happens to be nonzero. This gives a seed-defined bit-specific function with no lifted quantities. The two different tests on empty and singleton roots are essential.

There is a cost to explicitly constructing this function for every target. At the thinned widths there are fewer than $p$ legal tuples. Every distinct represented finite sum abscissa of $\mathcal A_1$ is therefore a genuine pole, except possibly at abscissas of slot $i$ where the subtraction can cancel it. For a fixed $T$, regular success with bit one implies the selected anchor is nonempty (use its full set if $i$ is another slot). If the selected full sumset has fewer than $B$ **group elements**, its tail sumset has fewer than $B$ by translation by any selected anchor point. Freeze the tails. In this event fewer than $B$ possible anchor points can complete $T$, each surviving with probability at most $\rho$. Thus $\Pr(U_{\rm reg},\text{bit}=1,\ |\mathrm{sumset}(\mathcal A_1)|<B)\le B\rho$. With $B=2^{(1-\eta)b}$ for any fixed $0<\eta<1$, this event conditioned on $U_{\rm reg}$ has probability $2^{-\eta b+o(b)}$ as an upper bound. Outside that exception on bit-one roots, an explicit reduced denominator has degree at least $(B-1)/2-2w$. This bound prices a full rational representation; evaluation at a single target need not construct it.

The checker evaluates literal products over $\mathbb Z/p^2\mathbb Z$ only for its reference answers. Its alternate evaluator uses base-field group additions, joins and polynomial remainder trees. [Finite check](../../checks/check_chord_difference.py) · [Recorded output](../../evidence/check_chord_difference.txt)
