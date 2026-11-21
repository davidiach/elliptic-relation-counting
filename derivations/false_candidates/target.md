# A signed-address filter using only abscissas inside the tree

Take a smooth short Weierstrass curve over $\mathbb F_p$, $p>3$, with rational group $G$ of order $N$ and fixed target $2T\ne O$. Write $\pi(P)=x(P)$ for the orbit under $P\mapsto-P$, with $\pi(O)=\infty$. On the project's independent mixed slots there is a deterministic relaxation of the signed chord test once a public offset is fixed. Every true signed tuple passes. At the original thinned widths the expected number of wrong accepted **signed addresses**, conditioned on a guarded singleton at a fixed target, is at most $(1+o(1))/(2m^2)$. The same bound controls surplus internal assignments. The construction still carries full-size abscissas and gives no improved end-to-end exponent.

## Moving the leaves before removing signs

Fix public rational offsets $R_i$ independently of the tuple. Put the orbit of $P_i+R_i$ at leaf $i$ and the orbit of $T+\sum_i R_i$ at the root. At each binary vertex, with children $\pi(U),\pi(V)$, allow either $\pi(U+V)$ or $\pi(U-V)$ at its parent. Internal vertices need only an abscissa or the identity marker. The leaf uses the **signed** point encoded by its address. Pairing the two signs of a leaf before translating it would change the test.

The local relation is the homogeneous three-point summation polynomial. For finite child abscissas $u,v$, the quadratic in the parent coordinate $z$ is

$$
(u-v)^2 z^2-2((u+v)(uv+A)+2B)z +(uv-A)^2-4B(u+v)
$$

for $y^2=x^3+Ax+B$. Its two projective roots are the permitted parent orbits, by the chord formulas. At $u=v$ and $u^3+Au+B\ne0$, one root is infinity and the other is the tangent double. At a two-torsion abscissa both roots are infinity since its derivative $3u^2+A$ is nonzero. Homogenizing in each coordinate gives a repeated root at the other child when one child is infinity, including when both are infinity. Thus rational inputs never produce an additional nonrational orbit. Equivalently these cases may be checked by the group law. [Existing summation polynomial and conventions](../implicit_chords/fiber.md)

For an assignment accepted by the tree, start with the **signed** lift $T+\sum R_i$ at the root. The local relation permits lifts of both children summing to that chosen parent, since a simultaneous change of sign reverses the sum. Continue to the leaves. Conversely, signed lifts define every internal abscissa uniquely. Thus acceptance is exactly

$$
\sum_i\epsilon_i(P_i+R_i)=T+\sum_i R_i
\quad\text{for some }\epsilon_i\in\{1,-1\}.
$$

The all-positive choice retains every exact witness. Every other accepting assignment comes from at least one sign pattern with a negative sign. Bounding the number of such patterns with multiplicity bounds both wrong addresses and surplus internal assignments, even if some shifted leaves are fixed by negation.

For a finite offset $R=(r,h)$ with $h\ne0$ and $x=x(P)\ne r$, a leaf can read the proposed $v=x(P+R)$ and infer the ordinate:

$$2h\,y(P)=f(x)+f(r)-(x-r)^2(v+x+r),\qquad f(x)=x^3+Ax+B.$$

The curve equation and canonical sign then verify the signed address. At $x=r$ handle $P=\pm R$ explicitly. When $R=T$, regular roots already avoid that abscissa. Computing this check, or looking up the shifted orbit of each allowed signed address, incurs its own cost.

## Independent offsets give a bound for each fixed seed

First suppose the $R_i$ are independent and uniform on $G$. For a wrong tuple and a nonempty set $J$ of flipped signs, acceptance for that sign choice implies

$$[2]\sum_{i\in J}R_i=\sum_i\epsilon_i P_i-T.$$

The left side is uniform on $2G$. A single branch therefore has probability at most $|G[2]|/N\le4/N$. Summing over all $2^m-1$ nontrivial patterns and all $D=\prod_i|F_i|$ tuples bounds the expected surplus assignments by $4D(2^m-1)/N$ on any fixed singleton. This includes extra branches of the correct address.

At the thinned widths, $M=2^b$, $D_0=2M/N$, $\rho=w/M$ and $\lambda=N^{m-1}\rho^m\le(4m^2D_0^{2m-2})^{-1}$. Hence $D/N\le(2w)^m/N=\lambda D_0^m\le D_0^{2-m}/(4m^2)$. Hasse gives $D_0\ge2(1-O(p^{-1/2}))$. For growing $m=O(\sqrt b)$ the bound above is at most $(4+o(1))/m^2$. Preparing uniform public offsets takes expected polynomial work per offset by sampling an abscissa and a sign with appropriate rejection of infinity and two-torsion multiplicities. They can also be sampled uniformly by rejection from signed encodings including one encoding for infinity.

## A fixed public offset at the original widths

Set every $R_i=R$. Require $Q_j=T+[2j]R\notin\{T,-T\}$ for $1\le j\le m$. One may take $R=T$ if $\operatorname{ord}(T)>2m+2$. For other targets all but at most $K=2\sum_{j=1}^m(2j)^2$ group elements work, since the fibers of multiplication by $2j$ contain at most $(2j)^2$ points. One can pick and check such an offset before drawing the seeds, with expected polynomial cost. If a rational point of order greater than $K$ is supplied, one can also test its first $K+1$ distinct multiples deterministically. The bounds hold for **each** eligible fixed offset, including when two different $j$ give the same $Q_j$.

Let $H_Q$ count all signed legal tuples of sum $Q$, and $U_T$ be the event of one ordinary hit, which also passes the guards. For a particular $J$ of size $j$, changing the signs of the selected leaf points permutes the legal addresses of the grid. The number of tuples satisfying that branch is therefore exactly $H_{Q_j}$. These counts may be correlated. It suffices to bound their sum.

For $Q\ne\pm T$ write $\nu=N\rho$ and let $s$ bound the runs of free bits in every fixed-bit slot. Then

$$
\mathbb E H_T H_Q
\le \rho^2(2w\nu)^{m-1}
   \left(1+\frac{4(m-1)2^s}{\nu}\right)
=\lambda^2D_0^{m-1}\left(1+\frac{4(m-1)2^s}{\nu}\right).
$$

To see this, fix both tails (slots 2 through $m$). If their forced affine-anchor abscissas are distinct, joint retention has probability at most $\rho^2$. Each remaining slot satisfies $\mathbb E|F_i|^2\le 2w\,\mathbb E|F_i|\le2w\nu$. This stronger second-moment estimate for the tail pays for the later sum over all signs.

If the anchor abscissas agree, the two anchors are equal or opposite, with retention probability at most $\rho$. In the equal case the tail displacements sum to $Q-T$; in the opposite case negate the first tail and their displacements sum to $Q+T$. Both are nonzero. Mark one slot with a nonzero displacement. After fixing the others this displacement is forced, and its expected pair multiplicity is at most $4\rho\,2^s w$. All other slots contribute their total pair masses, at most $2w\nu$. Summing over the two cases and the marked position proves the claim. [Mixed-slot displacement bound and guarded singleton probability](../sampled_slots/isolation.md)

Put $q_T=\Pr(U_T)$. Since $1_{U_T}\le H_T$, the conditional expected number of extra sign patterns across all addresses is at most

$$
\frac{(2^m-1)\lambda^2 D_0^{m-1}}{q_T}
\left(1+\frac{4(m-1)2^s}{\nu}\right)
\ \le\ \frac{1+o(1)}{2m^2}.
$$

Here $q_T=(1-o(1))\lambda$, $m2^s/w=o(1)$, and $D_0\ge2(1-O(p^{-1/2}))$. Conditioning further on regularity costs $1+o(1)$ under the mixed sampler. Markov gives one intended accepting assignment and no others with conditional probability $1-O(m^{-2})$. This estimate uses the singleton **seed law at each fixed target**. Unlike the independent-offset estimate it need not hold on each fixed seed. It also explains why deliberately using $R=O$ fails: the all-minus branch at $-T$ is the inverse of the actual witness on every such sign-symmetric grid.

## Cost in the inherited representation

A binary tree on $m$ leaves still uses $m-2$ unknown internal abscissas excluding the root. Each full canonical label uses $b$ bits and possibly an infinity flag. A direct implementation carrying these labels has $\Theta(mb)$ internal bits before any polynomial or leaf arithmetic. A membership check on a supplied address can instead enumerate its $2^m$ sign patterns with polynomial work per pattern, while searching all legal addresses still costs up to $(2w)^m$ such checks. Exhaustively checking the exact signed sum on a supplied address already takes polynomial time. The gain proved here is a one-sided, sign-free internal relation with a quantitative bound on wrong addresses; it does not give a cheaper search or counting exponent. [Cost model and previous conditional coefficient](../sampled_slots/README.md)

The finite check compares the homogeneous relation with all projective roots on two curves, computes cross-target moments over a complete mixed slot law, and counts surplus signs and wrong addresses separately on fixed-target singleton seeds. [Checker](../../checks/check_candidate_targets.py) · [Recorded output](../../evidence/check_candidate_targets.txt)
