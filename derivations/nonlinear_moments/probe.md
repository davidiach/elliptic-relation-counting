# Nonlinear readouts on isolated signed tuples

We tested the signed first and second determinant moments on prime fields through $p=65521$ and arities through six. Every non-witness control comes from the same sampled slots as the witness. A fitted nonlinear readout distinguishes the witness on a small selected population over $\mathbb F_{13}$. Its ranking advantage does not persist in the larger held-out populations.

## Computing two moments from the slots

Let $D(v)$ be the elliptic evaluation determinant of a legal signed tuple of address $v$, with endpoint row $e(-T)$. Equation values are in $\mathbb F_p$. Give each slot point the integer weight $\epsilon_i(P)=(-1)^{r_i\cdot v_i(P)}$ under a uniform binary challenge $r$. Define

$$
M_j(r)=\sum_{v\ {\rm legal}}(-1)^{r\cdot v}D(v)^j,\qquad j=0,1,2,
$$

where $M_0$ includes every legal tuple, even at a determinant zero. The constant column of $e$ is one. Subtract $e(-T)$ from each row and delete the constant column to obtain the $m$-vector $\eta(P)$. For each slot form $a_i=\sum_P\epsilon_i(P)\eta(P)$ and $G_i=\sum_P\epsilon_i(P)\eta(P)\eta(P)^\mathsf{T}$. Then $M_0=\prod_i\sum_P\epsilon_i(P)$, $M_1=(-1)^m\det(a_1,\ldots,a_m)$ with the $a_i$ as rows, and

$$
M_2=\sum_{I\subseteq [m]}(-1)^{m-|I|}\det\!\left(\sum_{i\in I} G_i\right).
$$

Indeed, expand $\det(\sum_i t_iG_i)$ by Cauchy–Binet. The coefficient of $t_1\cdots t_m$ chooses one point per slot and squares its determinant, multiplying by the signed point weights. Inclusion–exclusion extracts this coefficient since the polynomial has total degree $m$. This argument works in the base field for every prime in the project's class. At width $w$, one challenge costs $(w+2^m)\operatorname{poly}(b,m)$ expected bit operations with the inherited point lifting. The full-grid calculation in the checker serves as reference truth and is cross-checked at individual challenges against these independent formulas.

The fixed nonlinear scores are the quadratic characters of $M_2+M_1^2$ and $M_0M_2-M_1^2$, using zero at zero. There are also single-moment controls and two learned readouts. The learned feature jointly records the four quadratic characters of $M_1,M_2,M_2+M_1^2,M_0M_2-M_1^2$, giving at most 81 cells. Once fitted, each query needs a constant-size lookup after computing the moments. The finite calibration procedures enumerate complete challenges and use labeled examples. We have no uniformly cheap asymptotic calibration bound.

## The variance score reduces to an abscissa-only test

For $m\ge4$ there is a precise restriction on the score $f(r)=\chi(M_0M_2-M_1^2)$. Split the challenge into its abscissa bits $r_x$ and $m$ sign bits $t$. Suppose first that every retained affine point has a different negative. If any sign bit is one, the paired point weights in that slot cancel in $M_0$, giving $M_0=0$. If all sign bits are zero, the averaged first-moment rows vanish on every odd-in-$y$ column. At least two such columns exist for $m\ge4$, while the endpoint supplies only one additional row, giving $M_1=0$. Consequently

$$
 f(r_x,t)=\chi(-1)\,1_{\{M_1(r_x,t)\ne0\}}+
 1_{\{t=0\}}\,h(r_x),\qquad
 h(r_x)=\chi\!\left(M_0(r_x,0)M_2(r_x,0)\right).
$$

The first term is a scalar test on the first signed moment and obeys the earlier sign-slice bound. The second has coefficient $2^{-m}\widehat h(a)$ at *every* signed address $(a,s)$ above abscissa address $a$. In particular, with $U$ the unique guarded event and $(a_*,s_*)$ the true signed address, the triangle inequality gives on the mixed sampling law

$$
\left|\mathbb E[1_U|\widehat f(a_*,s_*)|]
-2^{-m}\mathbb E[1_U|\widehat h(a_*)|]\right|
\ \le\ 2^{m/2}\sqrt{\lambda/N}+12m\rho .
$$

Here we used the previous bound for the principal multiplicative character of $M_1$ including its $3m\rho$ allowance. An affine two-torsion point appears in some slot with probability at most $3m\rho$, and the residual error in the displayed decomposition has absolute value at most three on those seeds. The first term thus has negligible witness correlation at the isolating widths. Any larger advantage of the variance score must come from the abscissa-only second-moment test. This also gives a positive simplification: a sufficiently heavy abscissa coefficient can be recovered by parity listing on the shorter address, followed by verification of at most $2^m$ signs. We obtained no growing-arity lower bound on the correlation of $h$. [First-moment sign-slice bound](../sampled_slots/signed_moment.md) · [Parity listing and costs](../sampled_slots/address_parity.md)

There is an exact obstruction specific to the geometry of the *retained binary addresses*. Suppose the legal addresses form an affine subspace $a_0+L$ and write $M_0(r)=\sum_{a\in a_0+L}(-1)^{r\cdot a}$. Over any odd prime field, $M_0$ is nonzero exactly on $L^\perp$, by character orthogonality. For any bounded function of the moment values, challenge and public seeds that vanishes whenever $M_0=0$, the Walsh coefficients at any two legal addresses are **identical**. Indeed the two characters agree everywhere on $L^\perp$, which contains the entire support of the score. More generally, if translation by a legal-address difference $d$ preserves the whole address set, pairing $a$ with $a+d$ shows $M_0(r)=0$ when $r\cdot d=1$, and the same conclusion holds for each pair separated by $d$.

With both affine signs admitted and no two-torsion, apply this to the abscissa-only challenges. Each abscissa has the same number of signed lifts, so the factor of $2^m$ in $M_0$ is harmless. At width two every nonempty set of liftable abscissa addresses in one slot is a singleton or the entire one-bit space. Their product is affine. Thus $h=\chi(M_0M_2)$ gives exactly the same coefficient at every legal abscissa tuple on *every* such width-two root. This also applies at larger widths when the liftable address sets are affine, and to other nonlinear readouts gated on nonzero $M_0$. At large widths arbitrary liftable sets need not have a nonzero additive period.

## Paired controls and fits

For a bounded score $f(r)$ compute its normalized Walsh coefficients $\widehat f(v)$. On every evaluated seed the guarded witness is unique even among ordinary signed tuples. Measure the absolute coefficient at that witness $u$ and its average over **all other legal signed addresses** on that seed. The reported rank is the probability of beating a uniform wrong legal address, assigning half credit to ties. Every singleton seed receives its original conditional probability in the exact case and equal weight in sampled cases. The checker also reports the signed coefficient gap. Rank is an analysis of a full spectrum and does not supply an efficient address search.

The first table fit chooses each cell's sign to maximize the empirical signed witness-minus-control coefficient. The second trains on the difference of squared coefficients, allowing the witness correlation's sign to change between seeds. That objective is a quadratic form in the cell weights. We use coordinate ascent from eight fixed-seed starting points. Training the squared objective uses up to eight uniformly selected controls per seed. Evaluation of both trained rules uses every wrong legal tuple. Each field and arity has its own calibration.

For the larger fields, plant a uniform ordered affine tuple summing to $T$, then sample independent masks conditioned on retaining that tuple. This weights a fresh root by its ordinary count $H$. Keeping $H=X=1$ gives exactly the fresh law conditional on the unique guarded event. Training and evaluation use independent random streams. We also report the subpopulation where the witness is the only *determinant* zero among legal addresses. Repeated evaluation rows can otherwise give additional zeros.

## Results

The table shows values for the squared-coefficient fit. SE is the sample standard error of the paired absolute-coefficient gap. At $p=13$ all 8,192 equally weighted mixed seeds were included in evaluation before conditioning, so the conditional means for the fixed trained rule are exact up to display rounding.

| $p$ | $m$ | $w$ | Seeds | Witness | Other legal | Gap (SE) | Rank |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 | 3 | 2 | All 940 singletons | 0.09684 | 0.08862 | +0.00823 (exact) | 0.4997 |
| 13 | 3 | 2 | All 276 with one determinant zero | 0.15625 | 0.10211 | +0.05414 (exact) | 0.5801 |
| 101 | 3 | 2 | 1,200 singletons | 0.08810 | 0.08704 | +0.00106 (0.00246) | 0.5053 |
| 1,031 | 4 | 4 | 1,200 singletons | 0.04408 | 0.04508 | −0.00100 (0.00122) | 0.4949 |
| 1,031 | 5 | 2 | 1,200 singletons | 0.06551 | 0.06681 | −0.00130 (0.00155) | 0.4933 |
| 65,521 | 6 | 2 | 1,200 singletons | 0.02109 | 0.02183 | −0.00075 (0.00054) | 0.4846 |

At $p=1031$, the two fixed nonlinear scores had gaps $-0.00074$ and $+0.00004$ at arity four, and $+0.00040$ and $+0.00072$ at arity five. Almost all held-out roots had exactly one legal determinant zero (1,146 and 1,166 respectively out of 1,200). The trained power scores on those restricted cohorts gave ranks $0.4973$ and $0.4948$. The numerical advantage on the selected small field is insufficient to establish a correlation guarantee as $b$ and $m$ grow. A candidate end-to-end improvement still requires such a guarantee, the charged calibration or a uniform readout, and verified witness recovery. [Executable probe](../../checks/probe_nonlinear_moments.py)

At $p=65521$, the fixed score $\chi(M_2+M_1^2)$ gave rank $0.4783$ in the initial 1,200 six-slot roots, suggesting a possible preference in the reverse direction. We then selected that one score for independent repetition, using 5,000 new conditional roots for each target. At $T=G$ its rank was $0.5003$ (SE $0.0041$), and the mean absolute-coefficient gap was $+0.00019$ (SE $0.00026$). At $T=2G$ these were $0.4963$ (SE $0.0041$) and $-0.00037$ (SE $0.00026$). The reverse preference from the first run did not reproduce at its original size. These are two fixed targets on one finite curve. The repeat imposes no bound on other readouts.

The abscissa-only test $h$ gave exact ties at every legal abscissa tuple on each of the 4,920 and 4,909 corresponding roots with alternative abscissas, as predicted by the width-two obstruction. [Recorded output](../../evidence/probe_nonlinear_moments.txt)
