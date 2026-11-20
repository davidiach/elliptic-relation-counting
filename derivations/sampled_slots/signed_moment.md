# Why a cheap joint determinant score loses the address

For the sampled slots, there is a natural public-seed calculation with an affordable cost per query. At growing arity its correlation with the isolated relation is exponentially too small for the parity-recovery argument. This conclusion holds over general prime fields in the stated sampling regime.

## The score

Let $v(P)$ be the binary address of a signed point in its slot. It consists of its dyadic abscissa address and a bit for the ordinate. For a challenge $r=(r_1,\ldots,r_m)$, define a row for each slot by

$$
A_i(r_i)=\sum_{P\in F_i}(-1)^{r_i\cdot v(P)}\,e(P).
$$

Here $e$ is the usual dimension-$(m+1)$ elliptic evaluation row with columns $1,x,\ldots,x^{\lfloor(m+1)/2\rfloor}$ and $y,xy,\ldots,x^{\lfloor m/2\rfloor-1}y$. Append the row $e(-T)$. Denote the determinant of the resulting matrix over $\mathbb F_p$ by $z$. Multilinearity makes $z$ the Walsh-weighted first determinant moment on the entire product of slots. The bounded real score is $f_\sigma(r)=\left(\frac{z}{p}\right)$, with value zero at $z=0$.

A single score query can enumerate at most $2mw$ signed points and then use a polynomial-size field determinant. This costs $w\,\operatorname{poly}(b,m)$ expected bit operations with the inherited field routines, including constructing the slots. The finite experiments compute the *entire* Walsh transform to measure truth, which is much more expensive. A predictor only needs to answer its challenge.

## The sign bits split the matrix

Use $m\ge4$. Temporarily omit affine two-torsion points. Each pair $P,-P$ has the same abscissa address and opposite sign bits. Among the columns of $e$, let $o=\lfloor m/2\rfloor$ contain $y$ and let the other $m+1-o$ be even in $y$. If the sign-challenge bit in a slot is zero, its row vanishes on all odd columns. If it is one, the row vanishes on all even columns.

Fix the set $I$ of slots whose sign-challenge bit is one. For a nonzero determinant, $|I|$ must be $o$ or $o-1$. When $|I|=o$, the odd columns are filled by the rows of $I$, and the even columns are filled by the remaining slots and the endpoint row. When $|I|=o-1$, the endpoint row completes the odd block instead. This follows by expanding the endpoint row into its odd and even parts. In the surviving term a row and column permutation makes the matrix block diagonal. Both groups of slots are nonempty for $m\ge4$.

On each of these fixed sign patterns, the quadratic character of the determinant is therefore a constant sign times a product of two bounded functions, one on the seeds and challenge bits of $I$ and one on those of its complement. Multiply each function by the indicator of its prescribed local sign pattern to extend it by zero. In the norm averaging over independent seeds and uniform challenge bits, the product of the two norms is at most $2^{-m/2}$. Apply the two-block inequality separately on each slice. For $U=\{H=X=1\}$, let $v_*$ denote the signed address on $U$, and write $N=\#E(\mathbb F_p)$ and $\lambda=N^{m-1}\rho^m$. The result for the original slots is

$$
\mathbb E[1_U|\widehat f_\sigma(v_*)|]\ \le\
\left[\binom{m}{o}+\binom{m}{o-1}\right]2^{-m/2}\sqrt{\lambda/N}
+3m\rho .
$$

Indeed, at most three affine points have $y=0$ on a smooth short Weierstrass curve in characteristic greater than three. The chance any slot contains one is at most $3m\rho$. On every other seed the original score and the omitted-point score coincide. On the exceptional seeds the original absolute coefficient is at most one. The two-block proof bounds the omitted-point score against the *original* event $U$ by summing over all original legal target tuples, which avoids any independence assumption after conditioning. [Two-block inequality](address_parity.md)

At the sparse widths, $\Pr(U)\asymp\lambda=2^{-O(m+\log m)}$ and $\ell=b/m+O(1)$. For growing $m=O(\sqrt b)$ in the isolating regime, dividing the displayed bound by $\Pr(U)$ gives at most $2^{-b/2+O(m+\log m)}$. The polarity of the score may even be chosen separately on each seed in this comparison. At $m\asymp\sqrt b$ the recovery criterion asks for correlation at least $2^{-O(\sqrt b)}$. The score meets the per-query time budget and falls short on correlation.

The same proof works for any unit-group character extended by zero at zero, using complex absolute values and Parseval. More generally, for a bounded field readout $\phi$ write $\phi(z)=\phi(0)+\sum_\chi c_\chi\chi(z)$ by Fourier expansion of $\phi-\phi(0)$ on $\mathbb F_p^\times$. Each nonconstant-in-the-challenge contribution satisfies the same slice estimate on two-torsion-free roots. The constant score itself is a single separated product. Thus the numerator is at most $(|\phi(0)|+2^{m/2}\sum_\chi |c_\chi|)\sqrt{\lambda/N}+3m\rho$. This also rules out readouts with total absolute coefficient weight $2^{o(b)}$ at the stated growing arities. Arbitrary lookups may have much larger weight. Combining several moments nonlinearly is open.

A single prescribed nonzero field value gives an indicator with total absolute multiplicative Fourier weight one: all $p-1$ coefficients have magnitude $1/(p-1)$. A set of $k$ nonzero values has weight at most $k$ by addition. Consequently direct tests for individual determinant values and for buckets of $2^{o(b)}$ values also have insufficient correlation. The zero indicator is covered by taking the complement of the principal character extended by zero.

## Finite comparison

The table gives mean absolute Walsh coefficient **on the same unique guarded roots**. The controls choose a uniform address either from the full binary address space or from the tuples of retained signed points on that particular root. A legal control tuple need not hit the target or satisfy the addition guards. Here $G$ is the fixture's enumerated generator.

| Field $p$ | Target | Width | Conditional sample | True address | Uniform binary address | Uniform legal address |
|---:|---|---:|---|---:|---:|---:|
| 13 | $G$ | 2 | All 940 guarded singletons among 8,192 seeds | 0.056516 | 0.028840 | 0.060805 |
| 13 | $2G$ | 2 | All 1,054 guarded singletons among 8,192 seeds | 0.071928 | 0.032197 | 0.070855 |
| 101 | $G$ | 2 | 1,500 sampled singletons | 0.028979 | 0.035275 | 0.029559 |
| 1,031 | $G$ | 4 | 1,000 sampled singletons | 0.064668 | 0.010685 | 0.064587 |

On the larger fields the test plants an ordinary signed relation, then keeps draws with $H=X=1$. Planting weights each seed by $H$, which is identically one on this event, so the accepted draws have precisely the intended conditional law. At $p=1031$ the marked address has about six times the mean magnitude of a uniform binary address. A uniform legal address gives almost the same sample mean as the marked one. These finite comparisons use $m=3$, where an additional all-even slice is possible. The growing-arity bound above is a separate argument. [Probe and sign-slice checks](../../checks/probe_signed_moment.py) · [Recorded output](../../evidence/probe_signed_moment.txt)
