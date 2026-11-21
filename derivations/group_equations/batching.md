# How many curve equations can a random linear batch discard?

For the alternating chord chain with $m\ge3$ leaves there are $e=2m-2+\lfloor(m-1)/2\rfloor$ residuals over $\mathbb F_p$. Replacing them by $s\le e-3$ independent uniform linear combinations gives false acceptance with probability approaching one on the negative inputs described below, when $2^m/p\to0$. This concerns a verifier that publicly chooses the combinations and then searches for simultaneous solutions. [Source of the original chord equations](../../README.md#attribution-and-provenance)

## The residual map and its fibers

Take an odd prime $p>3$ and a nonsingular short Weierstrass curve with linear coefficient $a$. Fix a nonidentity affine target, and fix one legal liftable abscissa $x_j$ in each slot. Suppose $x_1\ne x_2$ and, for even $m$, $x_m\ne x(T)$. A negative restricted instance has no guarded signed tuple on those leaves summing to $T$. The same obstruction applies to a negative full slot instance whenever it permits such a choice of leaves.

Write $X_1=x_1$, $(X_m,Y_m)=T$, and let $d=(m-2)+(m-1)+\lfloor(m-1)/2\rfloor=e-1$. The $d$ free field coordinates are the internal $X_j$, the slopes $\lambda_j$, and the internal even $Y_j$. Use canonical field representatives. At step $j$, write $u=X_{j-1}$, $v=x_j$, $w=X_j$, and let

$$
E=\lambda^2-u-v-w,\quad
L=2\lambda Y_{j-1}-2u^2-uv-uw+vw-a,\quad
R=2\lambda Y_j+2w^2+uw+vw-uv+a .
$$

The residual vector records $E$ at every step, $L$ at odd steps and $R$ at even steps, and the usual links for odd $j\ge3$. Retain the original guard $u\ne v$ and further restrict the input coordinates to $w-v\ne0$ at even steps and $w-v-2E\ne0$ at odd steps. Under independent uniform coordinates, any nonautomatic guard or even leading-coefficient exclusion removes at most $1/p$ of assignments. For odd steps the leading-coefficient exclusion is quadratic in the free $\lambda_j$, with nonzero quadratic coefficient, so it removes at most $2/p$. At least

$$
M_0=(1-3(m-1)/p)\,p^d
$$

assignments survive. Only use this estimate when its prefactor is positive.

**Fiber bound.** Each prescribed residual vector is attained by at most $2^{m-1}$ surviving assignments. Define every missing odd ordinate, including $Y_1$, by $Y_i=\lambda_{i+1}(X_i-X_{i+1})-Y_{i+1}$. Then $C_j=\lambda_j(u-w)-Y_{j-1}-Y_j$ is zero at even steps and equals the recorded link at odd steps. This includes the terminal link if $m$ is odd.

Start at the public target and reconstruct backwards. Given the current $w,Y_j$, the known leaf $v$, and the prescribed step residuals $\epsilon,\eta$, we have $u=\lambda^2-v-w-\epsilon$. At an even step, $R=\eta$ becomes a quadratic in $\lambda$ with leading coefficient $w-v$. At an odd step use

$$
L+R=2(u-w)E-2\lambda C_j .
$$

Since $L=\eta$ and $C_j$ is prescribed, this gives a quadratic with leading coefficient $w-v-2\epsilon$. These coefficients are nonzero on the restricted domain. Each step thus has at most two slopes, and every slope determines both $u$ and $Y_{j-1}$. Branching over $m-1$ steps proves the bound. The first public abscissa and any additional guards can only decrease the fiber size.

## False acceptance after projection

On a negative instance the zero residual is absent by correctness of the full chord system. The image of the restricted assignments contains at least $M_0/2^{m-1}$ nonzero residual vectors. Since each one-dimensional subspace has $p-1$ nonzero elements, these vectors occupy at least $M_0/((p-1)2^{m-1})$ projective directions.

For a uniform $s$-by-$e$ matrix over $\mathbb F_p$, a fixed direction is killed with probability $p^{-s}$. The indicators for any two distinct directions are independent, because the two representatives are linearly independent. If $Z$ counts attained directions killed by the matrix, then $\operatorname{Var} Z\le\mathbb E Z$. Consequently

$$
\Pr(Z=0)\le\frac{1}{\mathbb E Z}
\le\frac{2^{m-1}(p-1)p^{s-e+1}}{1-3(m-1)/p}.
$$

Whenever $Z>0$, a guarded canonical assignment with a nonzero true residual satisfies every batched equation. At $s=e-3$, the upper bound for avoiding this false acceptance is $O(2^m/p)$. The same conclusion holds for fewer rows. This defeats the batched existential test itself. It does not estimate the cost of a different procedure that lists projected candidates, verifies them independently, or arranges fresh challenges after an irrevocable witness commitment.

## A degree requirement for universal polynomial checks

Suppose $F_1,\ldots,F_s\in\mathbb F_p[R_1,\ldots,R_e]$ vanish simultaneously exactly at the zero vector of the whole residual space. Then $\sum_j\deg F_j\ge e$. Otherwise the polynomial $\prod_j(1-F_j^{p-1})$ has degree below $e(p-1)$. In every monomial some coordinate has exponent below $p-1$, and summing over that coordinate gives zero modulo $p$, including when the exponent is zero. But summing this polynomial on the full space counts its one common zero and gives one. This is the elementary Chevalley–Warning argument. It allows an anisotropic quadratic norm to replace a pair of congruences and rules out $o(e)$ universal tests of bounded degree. A special residual image, a compact high-degree circuit, or additional certificates require separate analysis. [Warning's original paper](https://doi.org/10.1007/BF02940715)

## Finite verification

The regression checks existence directly by signed curve addition. It exhausts the free field coordinates, measures the residual fibers, and independently solves the backward quadratics on sampled full fibers. It then samples 64 independent uniform matrices at each indicated row count. [Reproduction script](../../checks/check_batch.py)

| Curve and fixed data | Surviving assignments | Distinct residual vectors | Largest fiber | Projective directions | Sampled false acceptances |
|---|---:|---:|---:|---:|---|
| $y^2=x^3+x+1$ over $\mathbb F_{13}$, leaves $(0,1,4)$, target $(0,1)$ | 22,308 | 20,201 | 4 of bound 4 | 15,221 | 64/64, 64/64, 33/64 at 2, 3, 4 of 5 rows |
| $y^2=x^3+1$ over $\mathbb F_7$, leaves $(0,1,0,0)$, target $(2,3)$ | 61,397 | 51,006 | 6 of bound 8 | 43,422 | 64/64, 62/64, 18/64 at 4, 5, 6 of 7 rows |

Both fixtures are negative. The distribution of the sampled matrices concerns only these two fixed leaf tuples. The asymptotic probability estimate above follows from the fiber argument.
