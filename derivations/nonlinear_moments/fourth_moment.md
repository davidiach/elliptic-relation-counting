# A fourth-moment test for a fixed nonlinear character score

Fix a seed, let $r$ range over the $D=2^\ell$ binary challenges on an $\ell$-bit full address, and let $T$ vary on the curve. Put $d=m+1$. Consider either of the two field functions

$$
F_r(T)=M_2(r,T)+M_1(r,T)^2,\qquad F_r(T)=M_2(r,T).
$$

The score is $f_T(r)=\chi(F_r(T))$, with the quadratic character extended by zero. Both field functions belong to $L(2dO)$ since each tuple determinant has pole order at most $d$ as a function of the endpoint. For each seed there is an explicit count that controls the ability of *any* signed address to stand out over uniformly masked targets.

For independent uniform $r_1,r_2,r_3$ put $r_4=r_1+r_2+r_3$, with binary addition, and form the rational function $G=\prod_{j=1}^4F_{r_j}$. Let $\eta_\sigma$ be the fraction of these quadruples for which $G$ is nonzero as a function and has **even order at every geometric point**. An identically zero product contributes zero to the score, including when only one factor is identically zero. The test includes even divisors of functions that need not be squares. This choice makes the result applicable after a fixed-cofactor pullback.

For a uniform nonidentity target in the large prime-order subgroup, on the project's fixed coprime-cofactor class, we obtain

$$
\mathbb E_T\sum_v|\widehat f_T(v)|^4
\le \eta_\sigma + O(h^2 d/\sqrt p),
$$

where $h$ is the cofactor. The sum over $v$ includes the full binary address space. In particular, if $U$ denotes a unique guarded tuple and $q=\Pr_{\sigma,T}(U)$, its conditional mean absolute coefficient satisfies

$$
\mathbb E_{\sigma,T}[|\widehat f_T(v_*)|\mid U]
\le \left((\mathbb E_\sigma\eta_\sigma+O(h^2d/\sqrt p))/q\right)^{1/4}.
$$

The seeds before conditioning must be independent of the masked target. Exponential decay of $\mathbb E_\sigma\eta_\sigma$ in $b$ would imply exponentially small correlation at the isolating widths. Establishing that decay for the sampled law at growing arity remains open.

## Proof and price of the exceptional quadruples

Expand the fourth power of the normalized Walsh transform and sum over addresses. Character orthogonality imposes the single condition $r_1+r_2+r_3+r_4=0$, with net normalization $D^{-3}$. The result equals the mean over these quadruples of $\chi(G(T))$. Nonzero even-divisor products cost at most one each. All other nonzero products have an odd-order zero or pole. Their character sums admit a uniform estimate.

Indeed, $G$ has at most $8d$ poles counted with multiplicities and the same bound on its zeros. For $p>h$, the map $[h]$ on the elliptic curve is separable and unramified of degree $h^2$. An odd-order geometric place of $G$ gives odd-order places of $G\circ[h]$. Hence adjoining a square root gives a geometrically connected degree-two cover of the elliptic curve, tamely ramified at the odd-order places. Riemann–Hurwitz gives genus $1+t/2$, where $t=O(h^2d)$ counts them. Away from zeros and poles, the number of rational points in each fiber is $1+\chi(G([h]P))$. The Hasse–Weil bounds for the two curves give an absolute character sum $O(h^2d\sqrt p)$, including $O(h^2d)$ possible endpoint corrections at zeros and poles. Multiplication by $h$ maps $E(\mathbb F_p)$ uniformly onto the subgroup of order $N/h$ for $N=hn$ and $\gcd(h,n)=1$. Dropping the finitely many preimages of the identity and averaging therefore gives $O(h^2d/\sqrt p)$. For a prime-order group one can take $h=1$. [Poonen, *Lectures on rational points on curves*](https://math.mit.edu/~poonen/papers/curves.pdf)

The resulting bound on the full fourth-power sum also dominates the fourth power of the coefficient of any selected legal address on $U$. Dividing by $q$ and applying the power-mean inequality proves the assertion. The method allows a large value on a few specially selected targets.

## Two exact finite comparisons

Over $\mathbb F_{1031}$, use one mixed-mask fixture for each of three and four slots, at width two. To get nonempty comparisons, the fixtures select disjoint slots with all four signed points present. They are finite examples rather than independent draws from the original ensemble. For each of the two scores and each seed we compute the parity of the divisor of every $F_r$ over the algebraic closure. Two extra targets check interpolation against direct signed sums. We count quadruples exactly by grouping pairs $(r,r+a)$ by their odd-divisor symmetric difference for each $a$. A nonzero quadruple has even divisor precisely when the pair fingerprints agree. This diagnostic uses $O(D^2)$ pair operations after constructing the fingerprints. Its computational cost is separate from the cost of one score query.

| Slots | Score | Nonzero products | Even divisors | All quadruples | Mean fourth-power sum over targets |
|---:|---|---:|---:|---:|---:|
| 3 | $\chi(M_2)$ | 154,112 | 9,296 | 262,144 | 0.035176 |
| 3 | $\chi(M_2+M_1^2)$ | 154,112 | 9,296 | 262,144 | 0.036127 |
| 4 | $\chi(M_2)$ | 16,777,216 | 196,096 | 16,777,216 | 0.011747 |
| 4 | $\chi(M_2+M_1^2)$ | 16,777,216 | 196,096 | 16,777,216 | 0.011633 |

**Every even divisor in these four rows came from pairwise repeated masks.** With $a$ masks having nonzero field functions, those unavoidable quadruples number $3a^2-2a$. Here $a=56$ and $256$ respectively.

The finite divisor calculation factors a target section as $g(x)(A(x)+B(x)y)$ with coprime $A,B$. An odd factor of $g$ toggles both sheets away from the branch locus of $x$. At a branch point it has even order on the curve. For the primitive part, factor $A^2-(x^3+ax+c)B^2$. An odd factor gives the unique odd sheet $y=-A/B$ at an unramified abscissa, or the sole ramified place above a branch point. Add the identity when the pole order is odd. Encoding these places by Frobenius orbits preserves symmetric differences. The checker also verifies that multiplying representative functions combines these fingerprints by symmetric difference. [Checker](../../checks/check_character_fourth.py) · [Recorded output](../../evidence/check_character_fourth.txt)
