# How many signed determinant values survive a Walsh challenge?

On a regular grid with exactly one target tuple, the signed moments do determine its address parity. A finite-field Hankel calculation recovers it using twice the number of distinct nonzero determinant values whose signed multiplicities survive. On the mixed sampling law, this number has conditional mean $2^{b-O(m\log m)}$ or larger even on almost all guarded singleton seeds. Thus explicitly identifying the minimal scalar recurrence remains expensive after signing the grid.

## A nonlinear decoder that is exactly correct

Use the [thinned mixed sampler](../sampled_slots/isolation.md) with one odd affine anchor, $M=2^b$, $\rho=w/M$, $N=\#E(\mathbb F_p)$, $\nu=N\rho$, $\lambda=\nu^m/N$, and an affine target $T$ with $2T\ne O$. Call a root regular if the slot abscissa sets are pairwise disjoint and avoid $x(T)$. Write $\Delta(a)$ for the determinant on the signed tuple at address $a$, with endpoint $-T$. On a regular grid it vanishes exactly when the points sum to $T$. Let $U$ be the event that there is one ordinary tuple and it passes the guards. [Elliptic determinant criterion](../joint_parity/target_incidence.tex)

For a binary challenge $r$ and field value $z$, write

$$
n_z=\#\{a:\Delta(a)=z\},\qquad
h_z(r)=\sum_{a:\Delta(a)=z}(-1)^{r\cdot a},\qquad
M_j(r)=\sum_z h_z(r)z^j\pmod p .
$$

The zeroth power here is constantly one. On regular roots in $U$, $h_0(r)=(-1)^{r\cdot a_*}$ and $h_0(r)=M_0(r)-M_{p-1}(r)$ in the field by Fermat. Let $t(r)$ count the nonzero $z$ at which $h_z(r)$ is nonzero in the field. The monic polynomial $Q(X)=\prod_{z\ne0,\ h_z(r)\ne0}(X-z)$ of degree $t=t(r)$ has $Q(0)\ne0$. Its coefficients $q_j$ give

$$
h_0(r)=Q(0)^{-1}\sum_{j=0}^t q_j M_j(r).
$$

Given an upper bound $L\ge t$, the $L$-square matrix with entries $M_{i+j+1}$, $0\le i,j<L$, factors as a Vandermonde matrix on the surviving nonzero values, a diagonal of $z h_z$, and the Vandermonde transpose. Its rank is exactly $t$, since the rectangular factors have full column and full row rank. With the rank known, the leading $t$-square matrix is invertible. Solving the ensuing monic recurrence gives $Q$ from $M_1,\ldots,M_{2t}$. For $t=0$, take $Q=1$. Conversely, multiplying any recurrence of the sequence starting at $M_1$ into its sum of exponentials and using Vandermonde inversion shows it vanishes at every one of these values. Its degree is at least $t$. This supplies an exact nonlinear decoder on spectra with small support. Knowing that the zero term alone is small imposes no such bound on $t$.

## An energy estimate for the actual mixed masks

Write $\mathcal C=\sum_z n_z^2$ and $d=m+1$. For the fixed target,

$$
\mathbb E[1_{\rm regular}\mathcal C]
\ \le\ 4\nu^m + d(N\rho^2+2\rho)(2w)^{2m-2}.
\tag{1}
$$

To prove this, freeze the non-anchor slots, discarding roots already irregular on these slots. Choose signed tails $A=(P_2,\ldots,P_m)$ and $B$ in them. The function $f_A(P)=\Delta(P,A,-T)$ is a nonzero section of $L(dO)$, vanishing at the $m$ distinct supplied points $P_2,\ldots,P_m,-T$. Indeed, the kernel of evaluation at these points is a degree-one divisor space on a genus-one curve and has dimension one by Riemann–Roch. Thus the $m$ fixed rows are independent and the determinant with its variable row is nonzero. It is nonconstant and has at most one additional zero counted with multiplicity. The same holds for $B$. We bound pairs of anchor points $P,Q$ that give equal determinant values, dropping the regularity condition involving the anchor to get an upper bound.

* For $x(P)\ne x(Q)$ there are at most $dN$ pairs, since each fixed $P$ allows at most $d$ values of $Q$. Distinct abscissas survive the affine anchor together with probability at most $\rho^2$.
* At equal abscissas, use $Q=\epsilon P$ for $\epsilon=\pm1$. If $f_A(P)-f_B(\epsilon P)$ is a nonzero function, it has at most $d$ zeros, each retained with probability at most $\rho$. Counting both signs also bounds two-torsion cases.
* If the functions are identical, the points $\epsilon B_j$ must all be zeros of $f_A$. Their abscissas belong to their respective slots, since the slots keep both signs. Slot $j$ can contain only $A_j$ among the prescribed zeros, by regularity. There is at most one additional zero in the whole divisor, and its abscissa can belong to at most one tail slot. Thus for fixed $A$ and $\epsilon$ there are at most **two** candidate tails $B$. Each identical pair costs at most $\mathbb E|F_1|\le\nu$ anchor choices. The expected number of tails $A$ is at most $\nu^{m-1}$ by independence.

There are at most $(2w)^{2m-2}$ ordered tail pairs. The three estimates give (1). Under the [isolating widths](../sampled_slots/isolation.md), $\lambda(2M/N)^{2m-2}\le1/(4m^2)$. Dividing the last term of (1) by $\nu^m$ gives exactly $d(1+2/\nu)\lambda(2M/N)^{2m-2}$. Hence the full expectation is $O(\nu^m)$.

## Conditioning on success and paying for signed cancellations

Use the growing regime $m=O(\sqrt b)$ of the isolating widths, with $m2^s/w=o(1)$ for the free-bit run bound $s$. Then $q=\Pr(U)\asymp\lambda$ and $\log_2(1/\lambda)=O(m+\log m)$. The chance of any irregularity is at most $m\rho+\binom m2 p\rho^2$ by one-point marginals and independence. This is $o(q/m^2)$. Conditional on $U$, with probability $1-O(m^{-2})$ all of the following hold.

1. The root is regular and $\mathcal C=O(m^2 N)$, by (1) and Markov's inequality applied before dividing by $q$.
2. All slots have at least $\nu/m^3$ signed points. Given any slot of size $\le h$, the conditional expected ordinary target count is at most $hN^{m-2}\rho^{m-1}=h\lambda/\nu$, since the last point of a completion is forced. Take $h=\nu/m^3$, sum over the slots and divide by $q$.

The number $L_0$ of legal signed tuples therefore satisfies

$$
L_0\ge(\nu/m^3)^m=N\lambda/m^{3m}.
$$

By Cauchy–Schwarz the number of occupied field values is at least $L_0^2/\mathcal C$. After excluding zero this leaves $R\ge 2^{b-O(m\log m)}$ distinct nonzero values. This statement concerns the actual curve and is conditional on a positive singleton. It does not assume that a determinant resembles a random field element.

For each such $z$, orthogonality of the binary characters gives $\mathbb E_r h_z(r)^2=n_z$ as an integer identity. Also $|h_z(r)|\le n_z\le L_0\le(2w)^m<p$ at the specified thinned widths for sufficiently large $b$. Thus every nonzero signed integer survives reduction in the field and $\Pr_r(h_z(r)\ne0)\ge1/n_z$. A second use of Cauchy–Schwarz yields

$$
\mathbb E_r t(r)\ge\sum_{z\ne0,\ n_z>0}\frac1{n_z}
\ge \frac{R^2}{L_0}\ge 2^{b-O(m\log m)} .
$$

This is a lower bound on **mean** order for a fresh uniform challenge at those roots. A decoder that constructs or outputs the full minimal recurrence therefore exceeds the targeted average budget. Other nonlinear readouts can try to isolate $h_0$ without finding the recurrence. Skipping challenges is a different method whose surviving correlation must be quantified.

## A uniform-in-challenge estimate under size bias by the legal grid

There is a stronger conclusion for a different sampling experiment. Keep $T$ fixed, draw one independent uniform affine marked point in each slot **without imposing a sum**, and draw each mask conditioned on its mark. This weights each seed by $L_0/\mathbb E L_0$. In this experiment most retained tuples have a determinant value occurring exactly once. The corresponding nonzero bins survive **every** sign challenge, including an adaptively chosen one.

To see this, count only off-diagonal pairs in (1). Outside the identical-section cases the old bounds still apply. If $f_A(P)=f_B(P)$ identically and $A\ne B$, some $B_j$ must be the single additional zero of $f_A$, since different slots have disjoint abscissas. In particular it is an affine zero distinct from the $m$ prescribed zeros. The pole at $O$ then has order $d$, and the divisor sum gives $\sum_{h=2}^m A_h+B_j=T$. With $m\ge3$, choose the two points $A_j,B_j$ from slot $j$, at most $(2w)^2$ choices. Conditional on that slot, the other $m-2$ slots have at most $N^{m-3}\rho^{m-2}$ completing tuples in expectation. Sum over $j$. If instead $f_A(P)=f_B(-P)$ identically, $T$ must also be a zero of $f_A$. Regularity and $2T\ne O$ make it the extra zero, hence $\sum_{h=2}^m A_h=O$ and $B_h=-A_h$ in every tail slot. There are at most $N^{m-2}\rho^{m-1}$ such tails in expectation. The independent anchor contributes at most $\nu$ point choices in either case.

Writing $D_0=2M/N$, division by $\nu^m$ gives the explicit bound

$$
\frac{\mathbb E[1_{\rm regular}(\mathcal C-L_0)]}{\nu^m}
\ \le\ d(1+2/\nu)\lambda D_0^{2m-2}+(m-1)D_0^2\rho+1/N .
$$

At the isolating widths this is $O(1/m)$. A bin of size $k\ge2$ contains $k$ nonunique tuples and contributes $k(k-1)\ge k$ ordered collisions. Since $\mathbb E L_0=((N-1)\rho)^m\sim\nu^m$, the expected fraction of nonunique tuples on regular size-biased roots is $O(1/m)$ (assign zero to this statistic on irregular roots). Irregularity itself has probability at most $(2w)^m/\mathbb E L_0$ times its fresh probability, which is $o(1/m)$ in the stated growing regime. A size-biased slot has probability at most $(\nu/m^3)/((N-1)\rho)$ of having at most $\nu/m^3$ points. By Markov and a union bound, with probability $1-O(m^{-1/2})$ the grid has $L_0\ge N\lambda/m^{3m}$ legal addresses and all challenges obey $t(r)\ge(1-m^{-1/2})L_0-1$. A bin occurring just once contributes a single sign and cannot disappear modulo an odd prime. Conditioning the marks to sum to $T$ changes their law; the earlier singleton statement controls that case in mean over the challenge.

## Finite checks

On the order-19 curve over $\mathbb F_{13}$ we enumerated 32,768 weighted width-one mixed seeds, checked (1), and compared all 1,792 candidate tail-section pairs at all affine points. On regular guarded singletons the mean nonzero support was $6.0769$ and the mean signed recurrence order was $5.6154$. On 256 conditional width-two singleton samples over $\mathbb F_{101}$ the corresponding means were $15.4570$ and $14.8159$. The checks verify signed Parseval and the order bound on each measured root. Fermat and independent Hankel reconstruction both returned the correct address parity in 1,760 challenge checks. [Reproduction script](../../checks/check_signed_spectrum.py)
