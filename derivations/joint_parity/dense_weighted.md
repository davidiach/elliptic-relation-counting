# What weighted moments recover on dense elliptic-curve slots

We obtained an exact three-slot readout on regular dense masks, using a number of field values proportional to the sum of the three slot sizes. The readout reconstructs one polynomial per slot and decides incidence using polynomial translation. Its balanced cost is quadratic up to logarithmic factors in the slot size. Directly visiting pairs already has this exponent. The asymptotic ECDLP bound therefore stays at the project's previous conditional value.

## Three slots of arbitrary size

Let $A_i$ be the distinct retained abscissas in slot $i$, with $k_i=|A_i|$, and let both signs of each point be present. Assume an odd prime-order rational group, characteristic $p>5$, and an affine target of abscissa $t$. Throughout the proof, all group operations are used for their mathematical meaning. The decoder uses no table of discrete-log coordinates. Such coordinates are used only by the independent small-curve test oracle.

For $d(u,v,w)=(u-v)^2(u-w)^2(v-w)^2$, the odd Gram block from the [quartic decoder](quartic_decoder.md) gives eight times $\sum d(u,v,w)$ over abscissa triples. Weighted versions give the same sum with a product of known slot weights. In channel $i$, weight the $i$th slot by $x^j$ and the other two by point tags, constant on sign pairs. After dividing out eight the measurements are

$$
 S_{i,j}=\sum_{x\in A_i} h_i(x)x^j,\qquad 0\le j<2k_i .
$$

The polynomial $h_i$ has degree at most four. For tags $e_v$ in one of the other slots form $a_r(X)=\sum_v e_v v^r(v-X)^2$ for $r=0,1,2$, and form $b_r(X)$ from the remaining slot. Then $h_i=a_2b_0-2a_1b_1+a_0b_2$. Thus obtaining all the measurements costs $O(\sum_i k_i^2+1)$ field operations by direct slot power sums, without visiting the tuple grid.

The monic polynomial of least degree that annihilates the sequence has roots exactly at $x\in A_i$ with $h_i(x)\ne0$. Indeed its Hankel matrix factors as a Vandermonde matrix, the diagonal of the visible amplitudes, and the transposed Vandermonde. Two samples per potential root suffice, including a root at zero. Berlekamp–Massey returns the recurrence in quadratic time. Accept the locator only when its degree equals the known $k_i$. Thus a vanished amplitude causes an explicit abstention, never a wrong support. [Massey's original paper](https://www.isiweb.ee.ethz.ch/archive/massey_pub/pdf/BI411.pdf)

On nonempty disjoint slots, every abscissa has a nonzero Vandermonde term. Consequently $h_i(x)$ as a polynomial in the tags of the other two slots is a nonzero bilinear polynomial. Uniform independent nonzero field tags make any such amplitude zero with probability at most $2/(p-1)$. The implementation starts with unit tags and takes the least common multiple of the visible locators across independent retaggings. After $s$ random rounds, the probability of missing any potentially visible abscissa is at most $(\sum_i k_i)(2/(p-1))^s$. Since $\sum_i k_i\le3p$, the expected number of rounds is bounded by an absolute constant for $p>5$.

**Completion on overlapping masks.** For $u\in A_i$, a Vandermonde term through $u$ exists exactly when the other two slots have distinct representatives different from $u$. Except when those two slots are the same singleton, at most two values of $u$ can fail this test. Indeed, a failure means one of the other sets is the singleton $\{u\}$, or their union after deleting $u$ has at most one member. When their union has more than two members only the singleton obstruction remains. When their union has one or two members, all obstructing $u$'s lie in it unless the two sets are the same singleton. Thus two extra known sums $\sum x,\sum x^2$ per slot suffice to complete the locator whenever at most two roots are missing. Subtract the visible roots' sums using their polynomial coefficients. For missing sums $s_1,s_2$, the missing quadratic has coefficients $1,-s_1,(s_1^2-s_2)/2$. A single missing root is $s_1$. Every completion is sound regardless of why the roots were invisible.

If two slots force the same abscissa $u=x(P)$, their signed sum is one of $O,2P,-2P$. The exact decision then needs at most three membership queries to the remaining original hash, at $x(T)$ and the finite abscissas of $T\pm2P$. These are rational by construction. The complete procedure checks empty slots and this exceptional case first, then uses the moment sequences, the two extra sums, and certified retries. It terminates with the exact bit in constant expected rounds on every legal three-slot input. This is a fixed-arity result.

## A decision from the three locator polynomials

For an abscissa $s$ and a monic locator $B$ define

$$
 D_{B,s}(X)=\prod_{B(w)=0} S_3(w,s,X)
 \quad\text{where}\quad
 S_3(w,s,X)=(s-X)^2w^2-2((s+X)(sX+a)+2c)w+(sX-a)^2-4c(s+X).
$$

The curve equation is $y^2=x^3+ax+c$. The third summation polynomial gives the finite abscissas of signed pair sums. If the two abscissas coincide, one sign choice gives the identity and the finite sum is given by doubling. [Semaev's original paper](https://eprint.iacr.org/2004/031)

Order the locators $C,B,A$ by nondecreasing degree. If $C(t)=0$ and $\gcd(A,B)$ has positive degree, both pairs can give the identity and the answer is positive. For finite matches, factor $D_{C,t}$. For each distinct rational root $s$, test whether $\gcd(A,D_{B,s})$ has positive degree. A match is equivalent to choosing signs in the three slots that sum to the specified target. All support points lift rationally. The case with target of opposite sign is equivalent by simultaneous sign reversal. This proof allows overlapping input sets.

Each translation can be built directly from the coefficients of $B$. Write the quadratic in $w$ as $Lw^2+Mw+C_0$ and work in the free rank-two polynomial algebra with $J^2=-MJ-LC_0$. For $k=\deg B$, form $U=L^k B(J/L)$ by homogeneous divide-and-conquer polynomial evaluation. Then

$$
 D_{B,s}(X)=\operatorname{Norm}(U)/L^k,\qquad
 \operatorname{Norm}(u+vJ)=u^2-Muv+LC_0v^2 .
$$

The equality follows first over the rational function field from the two roots of the quadratic and the monicity of $B$. The final division is exact as a polynomial identity, including at $X=s$. The implementation pads the coefficient vector for its balanced evaluation tree and accounts for the additional power of $L$. Degrees within a subtree are linear in its size. Fast polynomial arithmetic gives quasilinear field work per translation up to logarithmic factors. The decoder factors only the translated smallest slot, of degree at most $2k_{\min}$; it represents one other translated slot at a time. For all three sizes bounded by $n$, the complete method uses $\widetilde O(n^2\operatorname{poly}(\log p))$ bit work per successful attempt, in addition to enumerating inverse-hash windows. This expression includes finite-field root finding. A direct signed pair join has the same quadratic exponent and starts with the original inputs.

## Fresh tests with no training

Each row uses 3,000 independent fresh hash-seed triples at $T=(0,1)$. A regular input has no abscissa shared by slots and avoids the target abscissa. A dense input in the table is regular and has at least two distinct abscissas in every slot. The reference compares signed group sums independently. To measure what the weighted measurements achieve alone, each input gets one unit-tag attempt and up to seven independent retaggings before allowing the two extra sums or original hash queries.

| Field and hash widths | Dense draws | Dense negatives certified | Dense positives certified | Dense single-witness cases certified | Abstentions on all draws | True positives among abstentions |
|---|---:|---:|---:|---:|---:|---:|
| 103, (8,8,8) | 829 | 101 / 101 | 728 / 728 | 203 / 203 | 111 | 71 |
| 1031, (8,16,16) | 1,515 | 1,099 / 1,099 | 416 / 416 | 325 / 325 | 15 | 4 |
| 65521, (16,32,32) | 2,935 | 2,269 / 2,269 | 666 / 666 | 578 / 578 | 0 | 0 |

All returned decisions agreed with the reference. All 126 abstentions were irregular inputs. The two extra sums resolved 123, including every positive among them. The three remaining inputs each had a duplicated singleton and were resolved by the original hash membership test. The largest fixture uses a group of order 65,419 and $\lambda=N^2\prod_i(\tau_i/65536)$ within the project's target interval $1/(2\log_2 b)<\lambda\le1/\log_2 b$ at $b=16$. The implementation also checks the two previously colliding dense inputs, both the visible and invisible locator cases, symbolic polynomial products for degrees 1 through 17, and the complete decision on 1,980 additional support/target combinations permitting overlaps.

## Limits of replacing the support with a short random sketch

There is a distinction between weighting one slot by a known power, as above, and hiding all the slot data behind random weights. Give each abscissa an independent uniform field tag, use the same tag on both signs, and reveal only the resulting measurements. Put $d=m+1$. All single-slot Gram entries for the weights $\eta_x x^j$, $0\le j\le J$, are determined by $(\sum_x\eta_xx^\ell)_{0\le\ell\le J+d}$. This follows from the even/odd split in the curve basis. Products of two basis functions have pole order at most $2d$ and their even parts have abscissa degree at most $d$.

If each slot has at least $J+d+1$ abscissas, the relevant Vandermonde map is surjective and each displayed vector is **exactly uniform**, conditionally on the complete slot. Independent tag rounds preserve this property. Hence any joint statistic computed solely from those fully tagged vectors adds zero information about the seeds beyond separately supplied seed-only side information, on this stratum. This observation requires hidden uniform tags. Published tag values, deterministic weights, and measurements mixing tagged and untagged slots can behave differently. The exact decoder uses such mixed measurements.

## Which nonlinear real scores can reach the accuracy target?

The project already derived the [optimal prime-order real separated approximation](../mixed_coefficient/frontier.tex). Two extensions concern a cofactor and a nonlinear threshold with a margin. These are restricted-model results.

Take $N=|G|$, $T$ of odd prime order $r$, and $\kappa=|G[2]|\le4$. Split the independent seeds into two nonempty groups, in any predetermined way. A real rank-$R$ score is a sum of $R$ products of arbitrary real functions of the two halves. Subtract from the integer witness count $H$ the $\kappa$ real-character terms of its Fourier expansion. Their combined real rank is at most $\kappa$.

For the remaining $n=(N-\kappa)/2$ opposite character pairs, let $V_{x,j}=\sum_{P:x(P)=x}\chi_j(P)$, summing over affine points. These entries are real. Write $v_U=(\chi_j(U))_j$ for a nonidentity point of order two. Orthogonality and the singleton fibers give

$$
 V^\mathsf TV=NI-2\mathbf1\mathbf1^\mathsf T-\sum_{U\in G[2]\setminus\{O\}}v_U v_U^\mathsf T,\qquad V^\mathsf T\mathbf1=-\mathbf1 .
$$

The first correction removes the identity. Each additional term corrects a singleton fiber that full sign pairing would double count. Let $z_i$ be the slot indicator on abscissas, $\rho_i=\tau_i/2^b$, and let $B_i$ be the abscissa-to-residue indicator matrix modulo $\tau_i$. The [dyadic two-abscissa identity](../mixed_coefficient/mixed.tex) gives $\mathbb E z_i z_i^\mathsf T=\rho_i I+\rho_i^2(\mathbf1\mathbf1^\mathsf T-B_iB_i^\mathsf T)$. Consequently the Gram matrix $A_i=\mathbb E[V^\mathsf Tz_i z_i^\mathsf T V]/(N\rho_i)$ is a positive contraction with

$$
 \operatorname{tr}(I-A_i)\le (\kappa+1)n/N+\rho_i\,a_G,\qquad a_G=(N+\kappa-2)/2 .
$$

Here $a_G$ counts affine abscissas, including the singletons. The trace bound uses $\|V\|^2_{\rm op}\le N$ and $\|B_i\|_F^2=a_G$. For any nonempty half the Gram matrix of the products of its normalized features is the entrywise product of the slot Grams. It is a positive contraction, and its trace deficit is at most the sum of the slot deficits.

The remaining count kernel is $\sqrt{\lambda/N}\,U D V_0^*$ in probability-weighted row and column spaces, with $U,V_0$ contractions and diagonal entries $D_{jj}=2\operatorname{Re}\chi_j(T)$. The normalized squared singular-value tail after retaining $q$ terms is bounded above by the tail of the squared diagonal entries and below by that tail minus $4\sum_i\operatorname{tr}(I-A_i)$. Indeed each singular value is bounded by the corresponding absolute diagonal entry, while the total squared norm loses at most four times the trace deficits. The finite loss is $o(N)$ at the balanced tuning.

The characters of $G$ restrict equally often to each character of $\langle T\rangle$. Removing at most four real characters changes only $O(1)$ entries. For $r\to\infty$, the same cosine integral as in the prime-order proof applies. At the balanced sparse widths, any $q=o(N)$ rank approximation to the residual leaves $(1-o(1))\lambda$ squared error. The same conclusion applies to the full count, allowing a shift of at most $\kappa$ in the rank. It applies to the existence bit $Y$ because $\|Y-H\|_2^2\le 2\mathbb E[H(H-1)]=O(\lambda^2)$. The calculation was checked on actual noncyclic curve groups of orders 20 and 44, including all four real characters. For order 20 the probability-weighted count matrix was computed directly as a separate check.

**Margin consequence.** Suppose $|S|\le1$ everywhere and, outside fresh probability $O(\lambda^2)$, the score $S$ has the correct sign for $Y$ and absolute value at least $\gamma>0$. A polynomial in $S$, bounded in $[0,1]$, approximates $Y$ in mean square to $O(\lambda^2)$ with degree $D=O(\gamma^{-1}\log(1/\lambda))$.

To see this, write $S=\cos\theta$ and convolve the step $1_{\{\cos\theta\ge0\}}$ with the normalized Fejér kernel of order $h$. The result is an even trigonometric polynomial and thus a degree-$(h-1)$ polynomial in $S$, globally in $[0,1]$. At angular distance at least $\gamma$ from a jump its error is at most $\pi/(h\gamma)$, by integrating $K_h(u)\le\pi^2/(hu^2)$ against $du/(2\pi)$ on both tails. Taking $h\ge8\pi/\gamma$ gives error at most $1/8$ whenever $|S|\ge\gamma$. Apply the odd-binomial majority polynomial of degree $O(\log(1/\lambda))$ to reduce this error to $O(\lambda)$ while preserving the global bounds. The exceptional mass then contributes at most $O(\lambda^2)$ to squared error.

For a rank-$R$ score, that polynomial has real rank at most $\binom{R+D}{D}$, by expanding monomials in the separated terms. At $R\le2^{Cb/m}$ its logarithmic rank is $O(D b/m)$. Since $\log_2 N=b+O(1)$, $D=o(m)$ would contradict the squared-error lower bound. Quantitatively a successful such score must have $\gamma=O(\log(1/\lambda)/m)$, with a constant depending on $C$ and the accuracy constants. At the proposed growing arity this rules out a fixed positive margin. An exact finite-field gcd, a narrow-margin comparison, or a computation combining the full seeds jointly need not have this representation.

## What remains open

The dense example uses $2\sum_i k_i$ field measurements to recover the slot polynomials, information already present in the seeds. The improvement is a training-free readout that can distinguish a single witness in arbitrary-size regular three-slot pools. A faster group-sum decision at growing arity is still required. Enumerating translations on balanced blocks reproduces the growth of a partial-sum set. No estimate here meets the simultaneous growing-arity error and time budgets, so the best complete conditional bound is still $\log_2\mathbb E T\le(26.2933775156+o(1))b/(\log_2 b)^{1/4}$.

From the project root, with Python, NumPy and python-flint, run

~~~sh
python -B checks/check_dense_weighted.py
python -B checks/check_seed_rank.py
~~~
