# Pooling half sums by abscissa and canceling all adjacent poles

Fix a target $T$ with $2T\ne O$, $t=x(T)$, on a smooth short Weierstrass curve $y^2=x^3+ax+c$ over $\mathbb F_p$, $p>3$. Consider a challenge consisting of one **raw abscissa bit**. Restrict the half containing that slot to addresses with bit one. Let $A_P$ and $B_Q$ be the integer multiplicities of group sums in this selected half and the other half, respectively, *before attaching the endpoint*. Both are invariant under negation because the slots retain both signs and the restriction only reads the abscissa. A restriction on an ordinate sign does not have this symmetry.

On an empty regular grid, the nonlocal part of the normalized lifted-product difference, after dividing out the endpoint ordinate, is

$$
\sum_{P,Q} A_P B_Q h_t(P+Q),\qquad h_t(S)=1/(t-x(S)),\quad h_t(O)=0.
$$

The remaining slot correction is a sum over individual points. The following calculation identifies exactly which right-half coefficients survive after cancellation at this **fixed target**. [Reduction to the chord sum](../product_differences/chords.md)

## An explicit translation difference

Write $\bar A_s$ and $\bar B_u$ for the total multiplicities over the inversion orbits labeled by $s=x(P),u=x(Q)$, with infinity denoting the identity. Let

$$
H_t(s,u)=h_t(P+Q)+h_t(P-Q).
$$

The double sum equals $\frac12\sum_{s,u}\bar A_s\bar B_u H_t(s,u)$. Equal signs or an input of order two require no changes of weight. For finite inputs the kernel has the rational expression $H_t=\partial_t F_3/F_3$, since the two projective roots of

$$
F_3(s,u,t)=(s-u)^2t^2-2((s+u)(su+a)+2c)t+(su-a)^2-4c(s+u)
$$

are $x(P+Q)$ and $x(P-Q)$ with multiplicities. A root at infinity contributes zero to the logarithmic derivative. On the diagonal the rational continuation, or the two projective roots, gives the same expression whenever neither finite root is $t$. [Summation polynomials](https://eprint.iacr.org/2004/031)

After summing **all** terms on the right, define the rational function $R_B(x(P))=\sum_Q B_Q h_t(P+Q)$, extending around its poles in the function field. It descends to the abscissa line by inversion symmetry. Write $[R]$ for a pair of opposite points. Its complete partial fractions are

$$
R_B(s)=C_B+\frac1{y_T}\sum_{[R]:\,2R\ne O}
 \frac{y_R(B_{R+T}-B_{R-T})}{s-x_R},\qquad
C_B=\sum_{Q\notin\{O,T,-T\}}\frac{B_Q}{t-x_Q}
                    +\frac{B_T(3t^2+a)}{2y_T^2}.
$$

All coefficients are in $\mathbb F_p$. Indeed at $P=R$ the only possibly singular summands have $Q=-R\pm T$. Under invariant motion of $P$, $x(P+Q)$ has derivative $\pm2y_T$. Dividing by that derivative and changing the variable from the invariant parameter to $s=x(P)$ gives the displayed coefficient. Replacing $R$ by $-R$ changes both numerator signs and preserves it. An inversion-fixed $R$ cannot support a pole: each summand has at most a simple pole on the curve, whereas an abscissa function has even pole order there. To find the constant at infinity use the odd parameter $z=-x(P)/y(P)$ at $O$. The two terms with $Q=\pm T$ have $x(P+Q)=t\pm2y_Tz+(3t^2+a)z^2+O(z^3)$. Their simple poles cancel and their constants sum to the term in $C_B$.

Thus the reduced denominator has one linear factor per finite nonfixed orbit with $B_{R+T}\ne B_{R-T}\bmod p$. Whole stretches of constant multiplicity along addition by $2T$ disappear. For example, a $2T$-periodic right distribution gives a constant function. Given the differing coefficients and the constant, the formula only uses those surviving terms. Finding those coefficients cheaply from independent hash seeds is a separate operation.

On a promised singleton under the thinned widths, all counts here are below $p$. For a nonfixed left orbit of positive multiplicity the number of selected witnesses it contributes is $A_P(B_{T-P}+B_{T+P})$. The integer singleton promise means that at most one summand is positive. The orbit contributes a selected witness precisely when $x_P$ is a pole of $R_B$. Left points with $2P=O$ instead require checking $A_P B_{T-P}$ at at most four points. In fact a witness with fixed left sum must have every individual left point fixed: otherwise simultaneous reversal of the left points gives a second witness. A regular singleton with at least four left slots thus **cannot** have a fixed left sum, since only three distinct affine two-torsion abscissas exist. At smaller arity the correction is real. Over $\mathbb F_{29}$ on $y^2=x^3+5x+1$, the two left abscissas $15,21$ and the right abscissas $0,1$ give a unique relation for $(11,16)$ whose left sum is fixed. Thus an intersection with the reduced denominator's roots, with the exceptional checks at small arity, exactly determines the raw bit on a singleton. The pole is in the *left half abscissa variable*. An empty grid still requires the scalar value together with its local correction.

## Size after cancellation on the mixed sampler

Use a fixed balanced split with the marker in the left block. Let $k$ be the number of right slots and $B_Q$ their unweighted multiplicity. For the actual independent masks and any fixed nonzero displacement $D$,

$$
\mathbb E\sum_Q B_Q B_{Q+D}
 \le k (6\,2^{s_0}w^2/M)(2w)^{2k-2},
$$

where $s_0$ bounds the number of free-bit runs and $M=2^b$. To see this, two block tuples with difference $D$ have at least one nonzero slot displacement. Mark one such slot. Its displacement is forced by the others and has expected pair multiplicity at most $6\,2^{s_0}w^2/M$ from the individual hash laws. For each remaining slot the sum of these multiplicities is at most $(2w)^2$. This is the same argument as for the zero-displacement off-diagonal count $\sum_Q B_Q(B_Q-1)$. [Pointwise displacement bound](../sampled_slots/block_support.md)

At the growing thinned widths with $m=O(\sqrt b)$, $\log_2w=b/m+O(1)$, $s_0=o(b)$ and the isolation run condition, both upper bounds are $2^{o(b)}$ for $k\le\lceil m/2\rceil$. Conditioning on the fixed target's guarded singleton $U$ costs at most $2^{O(m+\log m)}$. Markov gives the same $2^{o(b)}$ upper scale with failure probability $O(m^{-2})$ conditional on $U$. On the same high-probability event there are $2^{b/2-o(b)}$ right tuples, because every slot has at least $\nu/m^3$ points, $\log_2\nu=b/m+O(1)$. All but $2^{o(b)}$ of these tuples have multiplicity one and a sum $Q$ with **neither $Q-2T$ nor $Q+2T$ represented**, applying the displacement bound at $D=2T\ne O$. Remove the at most constantly many $Q$ with $Q\pm T\in G[2]$. Each remaining $Q$ gives a genuine pole at $x(Q-T)$; a pole is obtained at most twice in this count. Thus the fully canceled, target-specific right transform already has $2^{b/2-o(b)}$ distinct poles with probability $1-O(m^{-2})$ conditional on $U$. This estimates the explicit reduced denominator without requiring any conjecture about the group order. [Widths](../sampled_slots/isolation.md) · [Conditional slot sizes](../nonlinear_moments/signed_spectrum.md)

Moving the translate difference through the convolution to a single original slot does not make its explicit list sparse either. With conditional probability $1-2^{-b+o(b)}$ every slot $F$ satisfies the [four-point translate lemma](../direct_parity/fiber.md): no four points summing to $O$ have an unpaired abscissa multiset. This implies $|F\cap(F+D)|\le4$ **for every nonzero $D$**, including a shift chosen after reading the seeds. Indeed any two pairs $(P,P+D)$ and $(Q,Q+D)$ have the same displacement. If one pair has distinct abscissas, every other pair must have those same two abscissas in some order for the four abscissas to pair. This allows at most four possible starting points. If all pairs have equal abscissas, each satisfies $2P=-D$, also permitting at most four starts. The difference of two translates of a slot indicator a nonzero distance apart therefore has at least $2|F|-8$ nonzero coefficients. This prices listing the local difference. Its implicit mask description may still be short.

## Rank on the occupied abscissas

For a generic right column $[Q]$, its two poles are $x(Q-T)$ and $x(Q+T)$. They are distinct and simple as long as $2Q\ne O$ and $Q\pm T\notin G[2]$. The function has no other poles, including infinity. Two such columns can share a pole only if their orbits coincide or differ by one step $[Q]\mapsto[Q\pm2T]$. In any collection of $n$ eligible columns one can consequently choose at least $n/3$ with pairwise disjoint poles by greedy selection. Linear combinations of $j$ such columns have at most $2j$ poles, and a nonzero combination is nonzero as a rational function because it has a pole unique to a participating column. Evaluation at more than $2j$ different ordinary left abscissas off their poles preserves linear independence.

On an empty grid every occupied left abscissa is off the occupied right poles. On a singleton remove the at most one left orbit taking part in a match to $T$ or $-T$: the negated full relation gives the same orbit. On the good event above, **if the bit equals one**, the selected left half has $2^{b/2-o(b)}$ distinct sums. Fix the marked point of the witness and vary the other left slots. This gives at least $(\nu/m^3)^{\lfloor m/2\rfloor-1}$ tuples, of which at most $2^{o(b)}$ have a repeated sum even within the unrestricted block, by the zero-displacement bound. Passing to inverse pairs loses at most a factor of two. The right supply has the same scale. Hence on these bit-one roots the matrix of the pooled kernel, restricted entirely to ordinary occupied off-pole labels, has rank $2^{b/2-o(b)}$ over $\mathbb F_p$. This rules out a short exact sum of separated left and right label functions for that matrix. It does not bound the time for a joint arithmetic circuit to evaluate one weighted scalar or test one pole from the seeds. Neither the boundary formula nor the rank estimate supplies an improved end-to-end exponent. [Finite check](../../checks/check_implicit_chord_fiber.py) · [Recorded output](../../evidence/check_implicit_chord_fiber.txt)
