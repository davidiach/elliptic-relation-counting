# Degree constraints on algebraic checks of an elliptic sum

Let $k$ be algebraically closed of characteristic greater than three, let $E:y^2=f(x)=x^3+ax+c$ be smooth, and write $O$ for its identity. Consider exact polynomial descriptions evaluated on a Zariski-dense open part of a relation variety. A lift may introduce scalar variables and project them away. Open guards may discard components and lower-dimensional exceptional sets. Any equality used to select a lower-dimensional subvariety must be charged as an equation.

## The degree at a fixed target

Fix $2T\ne O$ and $m\ge2$. In $(\mathbb P^1)^m$ let $W_{m,T}$ be the image, under the coordinatewise abscissa map, of the complete variety

$$
V_{m,T}=\{(P_1,\ldots,P_m)\in E^m:\ \sum P_i=T\}.
$$

The projective map is finite, and $V_{m,T}\cong E^{m-1}$ by forgetting the last point. Its image is a closed irreducible hypersurface. The irreducible defining polynomial on the affine chart of finite abscissas has **degree $2^{m-1}$ in each variable and total degree $m2^{m-1}$**. This is the fixed-target form of the familiar summation-polynomial degree phenomenon. [Semaev's original construction](https://eprint.iacr.org/2004/031)

To prove the separate degree, take generic independent points $P_1,\ldots,P_{m-1}$ and fix their abscissas. Their function field over $k(x_1,\ldots,x_{m-1})$ is Galois of degree $2^{m-1}$, from the independent separable double covers $E\to\mathbb P^1$. Its sign changes send the final abscissa to

$$
x\left(T-\sum_{i<m}\epsilon_iP_i\right),\qquad\epsilon_i\in\{1,-1\}.
$$

These values are generically distinct. Equality for two different sign vectors requires the corresponding points to be equal or opposite. In the first case a nonempty signed subset has double equal to $O$. In the second case the double of the signed sum over unchanged signs is $2T$. Each equation involving a nonempty subset cuts out a proper locus since varying one of its points makes the doubled sum nonconstant. For an empty unchanged subset the equation would be $2T=O$. Thus the orbit has size $2^{m-1}$, giving the minimal polynomial degree in the last variable. Symmetry gives the same degree in each variable.

For the total degree separately homogenize the affine equation in each coordinate pair. This defines its closure in $(\mathbb P^1)^m$. Indeed the affine polynomial is irreducible, and separate homogenization adds no divisor supported on a coordinate at infinity, since its leading coefficient in each individual variable is nonzero. The image of the complete sum fiber omits $(\infty,\ldots,\infty)$ because its only possible preimage is $(O,\ldots,O)$. Thus the separately homogenized equation has nonzero value at that point. Equivalently the monomial $\prod_i x_i^{2^{m-1}}$ has nonzero coefficient. This proves the claimed total degree even after specializing the target.

If a dense open of this hypersurface is exactly the image of an algebraic lift specified by equations of degrees $d_1,\ldots,d_r\ge1$ and open guards, then

$$
\prod_{j=1}^r d_j\ \ge\ m2^{m-1}.
$$

The affine Bézout inequality bounds the degree of every component of the equality locus by the product above. Open guards may remove whole components, and the closure of the projection of at least one surviving component must equal the desired irreducible hypersurface. Linear projection cannot increase degree. In particular an exact generic quadratic lift needs at least $m-1+\lceil\log_2m\rceil$ scalar equations, regardless of how many auxiliary scalars are allowed. This counts defining equations. Shared word computations, algorithms working solely on finite masks, randomized tests and bit-operation bounds require separate analysis. [Affine Bézout and linear projection, Krick–Pardo–Sombra, §1](https://mate.dm.uba.ar/~krick/KrPaSo01.pdf)

On the order-1019 curve $y^2=x^3+9x+1$ over $\mathbb F_{1031}$, the explicit choices $P_i=[2^{i-1}]G,\ T=[250]G,\ G=(0,1)$ give $4,8,16,32,64,128$ distinct final abscissas for $m=3,\ldots,8$. Direct enumeration of the signs reproduces these numbers with the existing group fixture. The proof establishes the generic degrees for every arity. [Finite curve fixture](../../checks/fixtures.py)

## Why two quadrics with one auxiliary scalar cannot encode one generic chord

Here all three abscissas $u,v,w$ vary. The chord hypersurface is the irreducible quartic

$$
(u-v)^2w^2-2\bigl((u+v)(uv+a)+2c\bigr)w+(uv-a)^2-4c(u+v)=0.
$$

Indeed with $P=(u,y),Q=(v,z)$ and $w=x(P+Q)$, the addition formula gives

$$
2yz=f(u)+f(v)-(w+u+v)(u-v)^2 .
$$

Consequently its function field is $k(u,v,yz)$, the fixed subfield of $k(E\times E)$ under simultaneous negation. The remaining degree-two extension of $k(u,v)$ is genuine since $f(u)f(v)$ is not a square there. The displayed quartic is the quadratic minimal equation in $w$ up to a scalar in $k(u,v)$, and its coefficients in $k[u,v]$ are primitive. Any common factor would divide $u-v$, while the coefficient of $w$ at $u=v$ is $-4f(u)$. Thus its projective closure in $\mathbb P^3$ is an integral quartic. Its smooth projective model is birational to the Kummer surface of $E\times E$. In odd characteristic, resolving the quotient by simultaneous negation gives a K3 surface. That surface has $h^2(\mathcal O)=1$ by its trivial canonical bundle and Serre duality. [Shioda, introduction, p. 718](https://www.jstage.jst.go.jp/article/pjab1945/50/9/50_9_718/_pdf)

Suppose two quadrics in coordinates $u,v,w,z'$, allowing the auxiliary scalar $z'$ and open guards, had exactly this generic image. If the quadrics share a linear factor, the complete zero set is a hyperplane together with a linear intersection. Neither component projects onto an irreducible quartic. Proportional quadrics give only quadric or linear components and also cannot have that image. Otherwise the projective intersection in $\mathbb P^4$ has total degree four. A component projecting dominantly onto the quartic surface must itself have degree at least four, and therefore exhausts the intersection with multiplicity one. As a complete intersection it has no embedded components, hence it is integral. The projection is generically finite. Its generic field extension is separable, since the forgotten scalar satisfies a polynomial of degree at most two over the image field and the characteristic exceeds three. At least one such polynomial is nonzero there, since a surface cannot have one-dimensional generic fibers over a surface. Two general linear sections of the quartic pull back to at least four times the extension degree in isolated intersection points on the source. Bézout forces degree one. The source would thus be birational to the Kummer surface.

For a projective complete intersection $X$ of two quadrics in $\mathbb P^4$ the Koszul resolution has terms $\mathcal O(-4),\mathcal O(-2)^{\oplus 2},\mathcal O$, and yields $H^2(X,\mathcal O_X)=0$. In particular $H^4(\mathbb P^4,\mathcal O(-4))=0$. This argument includes singular integral intersections. On passing to the normalization, $H^2(\mathcal O_X)$ surjects onto $H^2(\mathcal O_{X^\nu})$ since the cokernel of the map of structure sheaves is supported in dimension at most one. On resolving the normal surface, $H^2(\mathcal O_{X^\nu})$ surjects onto $H^2(\mathcal O_{\widetilde X})$ by Leray. The higher direct image in degree one is supported on finitely many points and the higher direct images in degrees at least two vanish by the fiber dimension bound. The smooth model therefore has $h^2(\mathcal O)=0$. This contradicts the K3 value, which is birationally invariant for smooth projective surfaces. [Cohomology of projective space](https://stacks.math.columbia.edu/tag/01XS) · [Proper fiber dimension and higher direct images](https://stacks.math.columbia.edu/tag/02V7)

This excludes the proposed two-quadrics, one-auxiliary, linear-projection encoding on the generic chord surface. It allows a projective singularity at infinity and makes no ordinary-curve assumption. Extra variables, coupled equations for several chords and arithmetic tests specialized to the rational points of one finite field are outside this obstruction.

As a separate finite check of the chord equation, take every rational endpoint $(w,t)$ with $t\ne0$ and every ordered pair of distinct field abscissas $u,v$. For each such input enumerate both sets of rational ordinates and count the signed pairs whose sum is the endpoint. The direct count equaled one at every quartic zero and zero everywhere else.

| $(p,a,c)$ | Group order | Quartic zeros | Other inputs |
|---|---:|---:|---:|
| $13,4,1$ | 19 | 288 | 2,520 |
| $17,2,2$ | 19 | 288 | 4,608 |
| $19,1,0$ | 20 | 308 | 5,848 |
| $37,3,1$ | 36 | 1,124 | 44,164 |

The degree examples and all four signed-pair counts can be reproduced together. [Finite verification](../../checks/check_geometry.py)
