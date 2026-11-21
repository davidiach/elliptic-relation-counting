# Moving the sampled points in their elliptic tangent directions

For an odd prime $p>3$, one quadratic-character bit of a deformed determinant product gives the Walsh parity of a unique witness. The unknown field factors from every nonwitness disappear under this character. With the full leading coefficient, one deformation instead recovers the entire address. The cost of evaluating either version from seeds is the outstanding step.

## First-order calculation

Write the curve as $y^2=x^3+Ax+B$. Normalize its translation-invariant tangent vector by $D(x)=2y$ and $D(y)=3x^2+A$, the vector dual to $dx/(2y)$. In a formal parameter $t$, move each retained signed point $P$ in slot $i$ with speed $c_i(P)$ along $D$ and move the target $T$ with speed $v$. These speeds determine the first derivatives of every evaluation row. The row at $-T$ moves with speed $-v$. Smoothness permits formal arcs with these first derivatives. Put $F(t)=\prod_a\Delta_a(t)$, over all legal signed tuples, using the same elliptic evaluation determinant as in the [section calculation](algebra.md).

On a regular grid the $m+1$ evaluated points are pairwise distinct. At a relation the partial derivative of the determinant with respect to the endpoint, with all slots fixed, is nonzero. Indeed the $m$ slot points impose independent conditions on $L((m+1)O)$ by Riemann–Roch. Their unique section up to scale has as its last zero the endpoint, distinct from the first $m$, so the endpoint is a simple zero. Locally the equation that the $m+1$ points sum to $O$ is smooth, since the endpoint maps to the sum by translation. It is a graph over the other $m$ points, and the determinant vanishes identically on this graph by the divisor criterion. In invariant tangent coordinates the differential of addition is the sum of the tangent coordinates. Consequently the derivative at that relation is

$$
u_a\biggl(\sum_i c_i(P_i)-v\biggr),\qquad u_a\in\mathbb F_p^\times.
$$

All $m+1$ row derivatives have the same coefficient $u_a$ when parametrized by the speed of their *evaluated* point. No group index or discrete logarithm is required to form a derivative row. The formula also applies in the presence of rational two-torsion if a retained point there uses its single legal address.

Let $K=\mathbb F_p(\sqrt d)$ for a quadratic nonresidue $d$, take all $c_i(P)$ in $\mathbb F_p$, and take $v=a_0+\sqrt d$. Each relation contributes one simple zero of $F$, since its velocity sum lies in the base field. Every other tuple contributes its nonzero base-field determinant at order zero. Thus the order is precisely the number $h$ of ordinary relations, and the first nonzero coefficient is

$$
[t^h]F(t)=U\prod_{a\,\mathrm{a\ relation}}\biggl(\sum_i c_i(P_i)-v\biggr),\qquad U\in\mathbb F_p^\times.
$$

For exactly one relation, write this coefficient as $X+Y\sqrt d$. Then $Y\ne0$ and $a_0-X/Y$ is its velocity sum. Choosing $c_i(P)$ to be its positionally weighted address digit returns the packed address in a single readout whenever the address alphabet product is at most $p$. The project achieves $(2w)^m<p$ at its thinned widths. Arbitrary labels in each slot can be used, including the bit for a Walsh challenge.

## One quadratic-character bit suffices

Choose the local speeds to be the Walsh bits in $\{0,1\}$, and put $k=\sum_i c_i(P_i)\in[0,m]$. For $m<p$, choose $a_0$ so that

$$
\chi_p((a_0-k)^2-d)=(-1)^k\quad(0\le k\le m).
$$

The same $a_0,d$ work for every seed, challenge and root at that arity. For fixed $d$, the number of such shifts differs from $p/2^{m+1}$ by at most $(m+1)\sqrt p$. Indeed the indicator is $2^{-(m+1)}\prod_{k=0}^m(1+(-1)^k\chi_p((a-k)^2-d))$. A product of $j$ distinct translated quadratics is squarefree of degree $2j$: their roots are $k\pm\sqrt d$ and all these roots are distinct. Weil's character sum bound limits its absolute complete sum to $(2j-1)\sqrt p$. Averaging these upper bounds over the subsets gives the stated error. For $(m+1)2^{m+1}=o(\sqrt p)$, as holds in the project's growing regimes, rejection sampling needs $O(2^m)$ expected shifts and $O(m\,2^m\operatorname{poly}(b))$ expected bit operations. A nonresidue is also found by uniform trials. [André Weil, *On Some Exponential Sums*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1079093/)

The quadratic character of $K$ equals $\chi_p$ of the norm and is one on every base-field unit. Applying it to $[t^h]F$ therefore gives $\prod_{a\,\mathrm{a\ relation}}(-1)^{r\cdot\operatorname{addr}(a)}$. For a singleton this is exactly $M_0-M_{p-1}$. For several witnesses it is the Walsh sign of the bitwise XOR of their addresses. For an empty grid of witnesses it is one. In the latter two cases the original residual is the *sum* of the witness signs. A constant local marker of one in one slot and zero elsewhere gives $(-1)^h$ in the same computation.

## Charged evaluation

Direct row differentiation uses polynomial work per tuple and one product of the leading coefficients, at cost $O((2w)^m\operatorname{poly}(b,m))$. Dual numbers suffice on a singleton. With more hits a dual-number product would be zero. Adding the valuations and multiplying the leading coefficients factor by factor handles all hits. The computation accepts the public sampled slots but enumerates their grid. A balanced group-sum join can recover the same output in $(2w)^{\lceil m/2\rceil}\operatorname{poly}(b,m)$ by directly matching halves. Neither estimate improves $b/(\log_2 b)^{1/4}$ in the logarithm of end-to-end time. Charging a quadratic character after first constructing an expensive scalar would still pay for that scalar. The open operation in this version is the character of the initial coefficient computed jointly from the slot seeds. [Balanced scalar join](scalar_join.md) · [Finite verification](../../checks/check_section_deformation.py)

## Three ordinary product evaluations at precision two

There is a related reduction which needs no field extension, fitted shift or differentiation by the product evaluator. Use the ring $R=\mathbb Z/p^2\mathbb Z$, a fixed lift of the Weierstrass coefficients, and arbitrary curve-point lifts of every retained signed point and of the endpoint $-T$. These lifts need not preserve addition. Smoothness supplies each lift by solving one linear equation over $\mathbb F_p$. Membership and digit labels always refer to the original sampled points. In all calls keep the same initial choice of lifts.

Let $V$ be the ordinary product of all the evaluation determinants on these lifts. Compute $V_E$ a second time after moving the **evaluated endpoint** $-T$ by $pD$, and $V_C$ a third time after moving every slot point by $p c_i(P)D$ from its original lift, keeping the endpoint fixed in this third call. Explicitly the displacement with speed $c$ is $p c(2y,3x^2+A)$ on the coordinates, which preserves the lifted curve equation modulo $p^2$. All products use the original ordered tuples.

* If $V$ is nonzero modulo $p$, there are no relations on this regular grid.
* If $V$ vanishes modulo $p$ and $V_E=V$ in $R$, there are at least two relations.
* Otherwise there is one relation. With $k=\sum_i c_i(P_i^*)\bmod p$ at its witness,

$$
k=\left((V_C-V)/p\right)\left((V_E-V)/p\right)^{-1}\pmod p.
$$

Here division by $p$ maps the ideal $pR$ to $\mathbb F_p$, so both quotients are well defined. For the proof, with two or more zero determinants each product contains at least two multiples of $p$, hence vanishes in $R$. With a single zero determinant, its baseline lift may be any multiple of $p$. Subtracting the baseline product cancels that choice. The first-order product rule modulo $p^2$ kills the contributions from perturbing nonwitness factors. The two differences are $pU$ and $pUk$ for the same $U\in\mathbb F_p^\times$, by the tangent calculation above. This proves both the classifier and the ratio even if $V=0$ in $R$ on a singleton.

The first two evaluations suffice for the three-way count classification. With positionally weighted digit speeds, the third gives the entire packed address at the project's thinned widths. On a known singleton even an equality-only oracle for two lifted products suffices: move just one slot with speed equal to one address bit there, and compare the result with $V$. The bit is zero exactly when the values agree. Repeating for the $m(\ell+1)$ address positions recovers the word. For a full challenge one can instead compare the combined-speed product with endpoint shifts $pjD$ for $0\le j\le m$; exactly one shift matches when $m<p$, and its index has the desired parity.

The packed construction charges **three scalar norm calls over $R$**, whose concrete cost is still exponential by direct expansion. An efficient joint routine for those three scalar values on the lifted sampled inputs would immediately give verified witness extraction on regular singletons. The checker uses ordinary modular determinants, two different choices of lifts and ground truth from the original curve. [Numerical results](../../evidence/check_section_deformation.txt)
