# Testing the product differences in one homogeneous degree

On a regular singleton the equality $V_C=V_{jE}\pmod {p^2}$ is an ordinary ideal-membership test over $\mathbb F_p$. The same graded quotient counts multiple witnesses, and a pairing of linear forms separates empty from singleton inputs. These constructions use polynomials in section coordinates and do not form either lifted grid product. They test equality or recover the ratio; the absolute nonzero difference digit is not needed for either task.

## The degree that records the witnesses

Let all $m$ sampled signed slots be nonempty. Through the fixed endpoint $-T$, a section is a point in $\mathbb P^{m-1}$. Each retained point $P$ imposes a linear equation $\ell_P(z)=0$. Put $n_i=|F_i|$, $Q_i=\prod_{P\in F_i}\ell_P$, $B=k[z_1,\ldots,z_m]/(Q_1,\ldots,Q_m)$ and $D=\sum_i(n_i-1)$, with $k=\mathbb F_p$. Assume regularity, meaning disjoint abscissas between slots and avoidance of $x(T)$. There are $h$ common projective points, in bijection with ordinary relation tuples. Each such section is rational and reduced, and exactly one factor of each $Q_i$ vanishes there. [Section argument](../direct_parity/algebra.md)

Then

$$
\dim_k B_D=\max(1,h).
$$

For $h>0$, evaluation on the $h$ projective representatives, with any consistent choice of scaling, identifies $B_D$ with $k^h$.

To prove this, first omit $Q_m$ to get the partial complete intersection $R$. Choose a linear chart coordinate $u$ nonzero at every partial section, extending scalars for this proof if necessary. Its finite affine algebra $A=(R[u^{-1}])_0$ has length $L=\prod_{i<m}n_i$. Write $s=\sum_{i<m}(n_i-1)$. Sending $K\in R_j$ to $K/u^j$ embeds $R_j$ in $A$. The partial Hilbert series gives $R_s=A$, and $R_{s-1}$ has codimension one there, taking $R_{-1}=0$ when $s=0$. The annihilator of $Q_m/u^{n_m}$ in $A$ is the $h$-dimensional span of the idempotents supported at the individual witnesses. At every other local factor $Q_m$ is a unit, even when that factor is nonreduced.

None of those individual idempotents lies in $R_{s-1}$. Indeed a homogeneous representative for one would be annihilated by the $m-1$ independent linear forms defining its point. Modulo $u$ their images generate the homogeneous maximal ideal, and the Artinian complete intersection $R/(u)$ has socle exactly in degree $s$. Any representative in smaller degree would be divisible by $u$ and could be divided repeatedly, an impossibility. This proof applies separately to each reduced witness even if other witnesses exist. Thus for $h>0$ that $h$-dimensional span meets the hyperplane $R_{s-1}$ in dimension $h-1$. The map $R_{s-1}\xrightarrow{Q_m}R_D$ has source dimension $L-1$, target dimension $L$ and kernel dimension $\max(0,h-1)$. Its cokernel is $B_D$, giving the formula. Surjectivity on the witnesses follows because $D\ge s$ and $R_D$ represents all of $A$. Evaluation commutes with field extension, so the conclusion also holds over $k$.

## The actual product equality

For any local tangent velocities $c_i(P)\in k$, in each slot set

$$
A_i=\sum_{P\in F_i} c_i(P)\prod_{Q\in F_i\setminus\{P\}}\ell_Q,
\qquad H_i=\sum_{P\in F_i}\prod_{Q\in F_i\setminus\{P\}}\ell_Q.
$$

Their short arithmetic descriptions can be built with a product tree and its derivative, or with prefix and suffix products, in $O(\sum_i n_i)$ operations on linear factors. At a witness section $z_*$, $H_i(z_*)\ne0$ and $A_i(z_*)/H_i(z_*)=c_i(P_i^*)$. At degree $D$ form

$$
G=\sum_i A_i\prod_{j\ne i}H_j,\qquad H=\prod_i H_i.
$$

On a singleton, $V_C=V_{jE}\pmod {p^2}$ if and only if $G-jH\in(Q_1,\ldots,Q_m)$. Indeed the tangent product rule identifies the equality with $\sum_i c_i(P_i^*)=j$, and evaluation identifies the degree-$D$ quotient with the one-dimensional space at the witness. The ratio $[G]/[H]$ of the two classes is the packed value when the speeds encode address digits. Ordinary homogeneous membership suffices, including on a positive root. On several witnesses this membership would require the same velocity sum at *all* witnesses, whereas precision-two product equality is then automatic. [Meaning of the lifted equality](../direct_parity/deformation.md)

## Compact probes of the top-degree pairing

When $\dim B_D=1$, let $\Lambda$ be any nonzero scalar coordinate on $B_D$, used consistently across queries. In the empty case the $m$ forms have height $m$, so they form an Artinian complete intersection with socle its degree-$D$ piece. The pairing $B_1\times B_{D-1}\to B_D$ is perfect: every nonzero homogeneous element generates an ideal meeting the socle, by taking an element of maximum degree in its ideal; symmetry gives the same statement on both sides. On a singleton this pairing instead has rank one, since evaluation expresses every entry as a multiple of $r(z_*)f(z_*)$. Assume $\dim B_1\ge2$ and $1\le D<p$, which hold when all growing sampled slots have the stated typical sizes.

Choose four independent uniform linear forms $r,s,u,v$. Test the determinant of

$$
\begin{pmatrix}
\Lambda(ru^{D-1}) & \Lambda(rv^{D-1})\\
\Lambda(su^{D-1}) & \Lambda(sv^{D-1})
\end{pmatrix}.
$$

It is always zero on singletons. On empty inputs this is a nonzero polynomial of total degree $2D$ in the probe coefficients. Indeed pure $(D-1)$st powers span the homogeneous polynomials in that degree, since all multinomial coefficients are nonzero and a polynomial of degree below $p$ vanishing at all field points is zero. The perfect pairing has rank at least two, so some choices give independent rows and columns. Random probes therefore give a false singleton with probability at most $2D/p$. One can also test the full pairing deterministically. The probes have only $O(m)$ field coefficients even when their degrees are high.

On a promised singleton, the vector $(\Lambda(z_j u^{D-1}))_{j=1}^m$ gives all projective coordinates of the witness section with probability at least $1-1/p$ for uniform $u$. Its single unknown multiplier is the same across coordinates. Factoring the section and checking its zeros against each slot recovers every address digit in expected $\operatorname{poly}(m,b)$ additional bit time. An algorithm that returned these normalized scalar evaluations directly could skip even the explicit velocity markers. [Recovering a tuple from section coordinates](../direct_parity/coefficient.md)

## Cost

A direct Macaulay calculation reduces the degree-$D$ probes against the multiples of the $Q_i$ in the same degree. In the partial quotient the presentation is a map from $L-1$ to $L$ dimensions. Representing those whole spaces costs at least $L$ entries. On almost all guarded singleton roots, even choosing which slot to omit gives $\log_2 L\ge(1-1/m)b-O(m\log m)$ in the thinned regime. Thus this product-free equality test gives no better end-to-end exponent. The short input to each pairing probe and the one-dimensional output leave open implicit computation of the probe from the seeds. [Conditional size estimate](../direct_parity/algebra.md) · [Executable checks](../../checks/check_difference_algebra.py)
