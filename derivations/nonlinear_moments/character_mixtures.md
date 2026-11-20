# What random field coefficients can add to a nonlinear score

Consider a nonlinear field expression $A_\sigma(r,T)$ in the signed determinant moments, followed by the quadratic character $\chi$ over an odd prime field, with $\chi(0)=0$. Randomizing the constants in this arithmetic is inexpensive. For two natural families, the probability of obtaining a useful witness signal from this randomization is exponentially small on the masked-target experiment. One of the families includes the experimentally tested expression $M_2+\alpha M_1^2$.

## Shifting an arbitrary joint expression

Let $T$ be uniform on $n$ nonidentity targets, independently of the seeds. Use $U$ for the unique guarded event, $q=\Pr(U)$, and $v_*$ for its signed address. Suppose that for each seed and challenge $A$ is a function of the target with no poles outside $O$ and pole order at most $D$. For example $A=M_2+M_1^2$ has $D\le2(m+1)$, even though its quadratic character can have much higher degree. Draw $c\in\mathbb F_p$ uniformly once and hold it fixed across the queries of an attempt. Then

$$
\mathbb E_{\sigma,T,c}[1_U|\widehat{\chi(A+c)}(v_*)|^2]\le (1+D)/n .
$$

For $x,y\in\mathbb F_p$, the exact normalized autocorrelation is

$$
\mathbb E_c\chi(c+x)\chi(c+y)=1_{\{x=y\}}-1/p .
$$

At $x=y$ the character square is one at every nonzero argument. Otherwise the map $c\mapsto(c+x)/(c+y)$ bijects the complement of $-y$ with the complement of $1$. Away from the excluded argument, the product equals the character of that ratio. The character sum over the whole field is zero. The identity follows. Expanding the squared Walsh coefficient therefore gives the collision energy of the raw $A$ at $v_*$ minus a nonnegative term. Distinct symbolic values of $A$ collide at no more than $D$ targets. The [two-challenge bound](collisions.md) finishes the proof.

## Varying the relative coefficient with no added constant

Write $d=m+1$, $A=M_2$ and $B=M_1^2$. At each challenge split the score into

$$
\chi(A+\alpha B)=
 1_{\{B=0\}}\chi(A)
 + g_\alpha,\qquad
 g_\alpha=1_{\{B\ne0\}}\chi(\alpha+A/B).
$$

Here $\chi(B)=1$ whenever $B\ne0$. For uniform $\alpha$ we have

$$
\mathbb E_{\sigma,T,\alpha}[1_U|\widehat g_\alpha(v_*)|^2]
\ \le\ [1+8(m+1)]/n.
$$

To prove it, give each active challenge ($B\ne0$) the label $A/B$. Inactive challenges contribute zero. Apply the autocorrelation identity to a pair of active challenges. Its $-1/p$ term contributes the negative of a squared coefficient of the active-set indicator divided by $p$. Thus the averaged square is at most the sum of the squared Walsh coefficients of the active label indicators.

For fixed seeds, group challenges with symbolically nonzero $B$ by their ratio as a function on the curve. The generic indicators are disjoint, so their squared coefficient sums over all addresses total at most one. For a pair of challenges whose denominators are symbolically nonzero, their specialized equal-label indicator agrees with the generic one except possibly at zeros of $B$, $B'$, and the cross product $AB'-A'B$ when it is not identically zero. Their pole bounds are respectively $2d$, $2d$ and $4d$. A symbolically inactive challenge is always inactive. Summing the resulting exceptional-target bound for the pairs gives the displayed estimate, as in the collision proof.

If one draws $Q$ uniform coefficients independently of the seeds and target and then chooses the **best coefficient for each witnessed instance**, the conditional mean of the largest absolute contribution from $g_\alpha$ is at most $\sqrt{Q(1+8d)/(nq)}$. This follows by bounding the square of the maximum by the sum of squares. At the isolating widths, $Q=2^{o(b)}$ still gives $2^{-b/2+o(b)}$ on large subgroups. The parameter-independent term $1_{\{M_1=0\}}\chi(M_2)$ survives every trial. The theorem does not constrain that term, specially chosen fixed coefficients, or coefficients computed algebraically from the input.

The finite check verifies the character identity on all field-value pairs over $\mathbb F_{13}$. It also averages the coefficients directly over every shift and every relative coefficient over $\mathbb F_{101}$ on the 102-target fixture and compares the results to the predicted collision energies. [Check and recorded output](../../checks/check_nonlinear_readout.py) · [Finite results](../../evidence/check_nonlinear_readout.txt)
