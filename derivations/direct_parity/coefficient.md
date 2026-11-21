# A single coefficient for signed witness parity

Fix an odd prime $p$, an affine target $T$ and sampled slots $F_1,\ldots,F_m$. Give each retained signed point its local Walsh weight $\chi_i(P)=(-1)^{\langle r_i,\operatorname{addr}_i(P)\rangle}$. Use the project's $m+1$ elliptic evaluation functions with constant first column. Subtract the row at $-T$ from each point row and discard the constant column to get $z(P)\in\mathbb F_p^m$. Up to a fixed sign the determinant $\Delta$ of a tuple is the determinant of its $m$ resulting rows.

Write $L_P(t)=z(P)\cdot t$. From each slot alone form the homogeneous polynomial

$$
S_i(t)=\sum_{P\in F_i}\chi_i(P)L_P(t)^{p-1}.
$$

For an arbitrary grid, including repeated rows, the exact field identity is

$$
M_0(r)-M_{p-1}(r)
=\prod_i\sum_{P\in F_i}\chi_i(P)
-[t_1^{p-1}\cdots t_m^{p-1}]\prod_i S_i(t).
$$

On a regular singleton this is $+1$ for an even witness parity and $-1$ for an odd one. The determinant moments of intermediate orders never enter the formula.

**Proof.** For a fixed tuple sum $\prod_i(1-L_{P_i}(t)^{p-1})$ over $t\in\mathbb F_p^m$. It counts kernel vectors and is congruent to one for full rank and zero for deficient rank. A univariate monomial sums to zero in the field unless its exponent is a positive multiple of $p-1$, when it sums to $-1$. The total degree is at most $m(p-1)$, so only the indicated top monomial survives. Its two signs $(-1)^m$ cancel, giving the top coefficient of $\prod_i L_{P_i}^{p-1}$. By Fermat this equals $\Delta^{p-1}$. Sum over tuples using the factored weights.

**What an explicit evaluation costs.** Each field evaluation of the product uses at most $2mw$ retained rows and polynomial overhead in the bit length and arity. Summing the product over projective directions with weight $(-1)^{m+1}$ extracts the coefficient: every nonzero vector occurs in an orbit of $p-1$ nonzero multiples with identical value. This takes $(p^m-1)/(p-1)$ directions. Alternatively a dense coefficient tensor with each exponent capped at $p-1$ has $p^m$ positions. Neither calculation meets the exponent budget. These costs leave room for a symbolic algorithm exploiting the restrictions on elliptic rows and on how each seed chooses them.

**Can one just sample fewer directions?** For universal projective evaluation formulas there is an exact obstruction. To isolate the direction spanned by $e_m$, use

$$
B(t)=t_m^{p-1}\prod_{j<m}(t_m^{p-1}-t_j^{p-1}).
$$

This is a product of allowed signed sums of powers of arbitrary linear forms, each with at most two terms. Its projective values are one at $e_m$ and zero elsewhere. A linear change of coordinates places its sole nonzero at any direction. Its coefficient is $(-1)^{m+1}$ in every case by the projective sum. Thus an exact linear formula using point values of the whole product needs every projective direction with that weight. An adaptive algorithm restricted to whole-product point queries has probability at most $k(p-1)/(p^m-1)$ of querying a uniformly hidden nonzero in $k$ queries, by following its all-zero transcript. These examples concern unrestricted row directions. They prove no such sampling bound for algorithms reading the actual elliptic row vectors, the hash seeds, or individual slot evaluations.

On the actual regular singleton there is a related exact sparsity statement. Let $H_i(t)=\sum_{P\in F_i}\chi_i(P)1_{\{L_P(t)=0\}}$. Then $\prod_i H_i(t)$ is supported at the unique projective kernel direction of the witness, where its value is the witness sign. Indeed, expanding the product selects precisely those tuples killed by $t$. The witness has corank one by independence of any $m$ distinct evaluation points in the degree-$(m+1)$ linear system on the elliptic curve. A direct proof and a way to retrieve this direction using equations for the slots appear in the [section calculation](algebra.md). For the classical grid coefficient formula used above, see [Schauz, *Quantitative Combinatorial Nullstellensatz*](https://arxiv.org/abs/1212.5694).

## Pricing a section query from the seeds

There is a cheaper way to evaluate all the local incidences at a *specified* direction. Write the section as $A(x)+yB(x)$, placing the constant so that it vanishes at $-T$. For the curve $y^2=g(x)$, list the field roots of $A^2-gB^2$, of degree at most $m+1$ by the pole orders of $x,y$. The polynomial is nonzero because $g$ is a square-free cubic and the section is nonzero. At a root with $B(x)\ne0$, the ordinate must be $-A(x)/B(x)$. When $A(x)=B(x)=0$, retain every rational lift at this abscissa. Discard the point at infinity throughout, since slots contain affine points. Using modular powering to take the gcd with $X^p-X$ first removes multiplicities and nonrational roots, and finite-field factorization lists the remaining roots in expected $\operatorname{poly}(m,b)$ bit time. Test the at most $m+1$ affine points against each seed, calculate their signed addresses, and add their weights. This computes every $H_i(t)$ and their product in expected polynomial bit time, **independent of the slot width**. On a regular singleton, testing the right direction even returns the entire tuple at this cost. [Cantor and Zassenhaus, finite-field factorization](https://doi.org/10.1090/S0025-5718-1981-0606517-5)

How expensive is choosing that direction? Write $U_{\rm reg}$ for the regular guarded singleton event, $q_{\rm reg}=\Pr(U_{\rm reg})\asymp N^{m-1}\rho^m$, $N=\#E(\mathbb F_p)$, $\rho=w/2^b$. Fix the target and a challenge independent of the seeds. Give a searcher any specified $m-k$ full slot seeds, with $1\le k\le m$. Suppose the other $k$ seeds enter only through queried values of the whole product $\prod_i H_i(t)$. On $U_{\rm reg}$ every query before the witness section answers zero. For $J$ such adaptive queries,

$$
\Pr(\text{queries the witness section}\mid U_{\rm reg})
\le \min\!\left(1,\frac{J\,k!\,\rho^k}{q_{\rm reg}}\right).
$$

To see this, freeze the known seeds and coins and follow the all-zero answers. Each tested section already vanishes at $-T$ and so has at most $m$ remaining distinct rational zeros. A regular witness uses all $m$. It must have exactly one of them in each visible slot, since abscissas of distinct slots are disjoint on this event. If the visible tests pass, at most $k!$ assignments of the remaining zeros to hidden slots are possible. Each assignment survives with probability $\rho^k$ by independent slot sampling. Union over the fixed path before conditioning. The case $k=m$ bounds the hit probability by $O(Jm!/N^{m-1})$.

Even if only **one** seed is hidden and its answers expose the entire list of retained points on each tested section, seeing a first hidden retained point other than $-T$ on a regular singleton has probability at most $Jm\rho/q_{\rm reg}$. Follow the empty-list transcript, excluding $-T$, which is absent on regular roots. The at most $Jm$ candidate points are fixed after revealing all other seeds.

**Parity prediction, including guesses on misses.** Let $f_\sigma(r)$ be any such adaptive score with modulus at most one on every transcript, including the empty transcript. With uniform fresh challenge $r$ independent of the masks,

$$
\mathbb E_{\sigma\mid U_{\rm reg}}
\left|\mathbb E_r f_\sigma(r)(-1)^{r\cdot a_*}\right|
\le \sqrt{\rho/q_{\rm reg}}+2Jm\rho/q_{\rm reg}.
$$

Fix the visible seeds $V$ and use the empty-list responses to define the baseline score $g_V(r)$. For each full address $a$, its visible digits determine every visible point or are illegal. A witness with this address forces one particular missing group point, which survives its independent mask with probability at most $\rho$. Thus $p_a=\Pr(U_{\rm reg},a_*=a\mid V)\le\rho$. By Parseval, $\sum_a|\widehat g_V(a)|^2\le1$. Cauchy–Schwarz bounds $\sum_a p_a|\widehat g_V(a)|$ by $\sqrt{\rho\Pr(U_{\rm reg}\mid V)}$. Average over $V$ and divide by $q_{\rm reg}$. Finally, for each fixed $r$ and $V$, a departure on $U_{\rm reg}$ from the empty-list path requires a first retained point on one of its $J$ predetermined sections. The preceding union bound applies even though the sections can depend on $r$. Each departure changes the score by at most two. Independent coins can be fixed and then averaged. For a score in $\{-1,+1\}$ the signed correlation is twice parity prediction accuracy minus one.

At the thinned widths, constant conditional correlation requires $J\ge 2^{(1-1/m)b-O(m+\log m)}$. Even conditional correlation $2^{-o(b)}$ needs $J\ge 2^{(1-1/m)b-o(b)}$ in the growing-arity regime. This bounds a parity predictor directly, including one that can guess without finding a section. A predictor with simultaneous arithmetic access to all seed words is outside this query model.

There is an intermediate case. Give the predictor a label of at most $L$ possible values prepared from **all** the seeds and $T$ before the challenge. The same correlation bound becomes $\sqrt{L\rho/q_{\rm reg}}+2JmL\rho/q_{\rm reg}$. Indeed, for each label $z$ also impose that label in the definition of $p_a$; the individual probabilities are still at most $\rho$. The Cauchy–Schwarz argument summed over labels now costs $\sqrt L$. For the first hidden response, union over the $L$ candidate empty-list paths before conditioning. Thus for conditional correlation $2^{-o(b)}$, the advice and query budget must satisfy $\log_2 L+\log_2(J+1)\ge (1-1/m)b-o(b)$. This permits arbitrary target-specific preprocessing and charges its output for one hidden slot. Feeding all raw seed words to a joint challenge-dependent circuit is still permitted by the original problem.

[Executable coefficient and projective tests](../../checks/check_projective_coefficient.py)
