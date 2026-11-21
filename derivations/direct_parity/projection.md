# Target-specific sparse group-character projections

The direct group-incidence formula gives a separable parity feature for each complex group character. Even choosing a different set of features for every seed and target has a limitation. On the project's mixed slots, any normalized linear combination of at most $2^{(1/2-\varepsilon)b}$ such features has exponentially small witness correlation for every fixed $0<\varepsilon<1/2$. This applies to arbitrary curve order even when choosing features and coefficients is free. The choice must precede the binary challenge.

## Statement

Let $E/\mathbb F_p$ be a smooth short Weierstrass curve, $G=E(\mathbb F_p)$ and $N=|G|$. Set $b=\lceil\log_2 p\rceil$, $M=2^b$, $w=2^\ell$, and $\rho=w/M$. Use the independent mixed slots $F_i$, injective signed addresses $a_i(P)$, and thinned widths of the [isolation theorem](../sampled_slots/isolation.md). Assume $m\to\infty$, $m=O(\sqrt b)$, $\ell=b/m+O(1)$, $m2^s/w=o(1)$ and $s\log b=o(b)$, where $s$ bounds the number of runs of free bits. These assumptions cover both contiguous slots and the project's strided layout. For a fixed target with $2T\ne O$, the guarded singleton event $U$ has probability $q\ge 2^{-O(m+\log m)}$. On $U$ write $a_*$ for the concatenated witness address.

For a group character $\xi$, put

$$
Z_\xi(r)=\prod_i\sum_{P\in F_i}(-1)^{r_i\cdot a_i(P)}\xi(P).
$$

Fix $0<\theta<1/2$. For each seed and target choose an arbitrary set $A$ of at most $2^{\theta b}$ characters, and coefficients $c_\xi\in\mathbb C$, independently of the fresh uniform challenge $r$. Set $f(r)=\sum_{\xi\in A}c_\xi Z_\xi(r)$ and require $\mathbb E_r|f(r)|^2\le1$ on the singleton seeds. In particular, bounded linear scores satisfy this condition. Then for some $c_\theta>0$, uniformly in the target,

$$
\mathbb E_{\sigma\mid U}
\left|\mathbb E_r f(r)(-1)^{r\cdot a_*}\right|
\le 2^{-c_\theta b+o(b)}.
$$

The [exact signed fiber](fiber.md) is the full character average with target phases. On regular singleton roots this fiber is the signed lift of $M_0-M_{p-1}$. The number of group modes in the restricted linear score and the cost of computing an arbitrary nonlinear decoder are different quantities. For instance, clipping a large linear score to its sign changes the model.

## Exact evaluation on one realized grid

When $N$ is prime, for a fixed singleton the minimum number of modes needed for an exact linear score on **all** challenges equals $K=|F_1+\cdots+F_m|$. Indeed, Walsh inversion in the signed addresses identifies the coefficient at each legal address as $B(P_1+\cdots+P_m)$, where $B(S)=\sum_\xi c_\xi\xi(S)$. Equality to witness parity forces $B(T)=1$ and zeros at all $K-1$ other represented sums. On a cyclic group of prime order $N$, a nonzero complex function with $k$ Fourier modes has at most $k-1$ zeros by the prime-order uncertainty principle. Thus $k\ge K$. Conversely let $\xi_0$ be a faithful character. Interpolate a polynomial of degree at most $K-1$ on the $K$ distinct complex numbers $\xi_0(S)$, taking value one at $T$ and zero at the other sums. Its powers use at most $K$ consecutive modes and yield the exact score. [Tao, *An uncertainty principle for cyclic groups of prime order*](https://arxiv.org/abs/math/0308286)

For every fixed target, $\Pr(U,K<B)\le\rho B$. In fact on $U$ the anchor is nonempty, so the tail sumset has cardinality at most $K$ by translating it by any anchor point. Freeze the tail. When its sumset has fewer than $B$ elements, fewer than $B$ possible anchor points can complete the target. Their combined acceptance probability is at most $\rho B$. Setting $B=2^{(1-\eta)b}$ for fixed $\eta>0$, exact evaluation requires at least $B$ modes outside conditional probability $2^{-\eta b+o(b)}$. The interpolation construction can require knowledge of the entire sumset.

## An even moment for one actual fixed-bit slot

Write $E:y^2=g(x)$. For each fixed integer $h\ge1$ choose a constant $B_h$ greater than every degree bound for $\Phi$ in the proof below. Take any group character $\omega$ whose order exceeds $B_h$. Let $t(x)=\sum_{P:x(P)=x}\omega(P)$. Opposite lifts are conjugate, so $t$ is real and $|t|\le2$. The uniform fixed-bit slot satisfies

$$
\mathbb E\left|\sum_{P\in F_i}\omega(P)\right|^{2h}
\ll_h w^h + \frac{w^{2h+1}(Cb)^{s+1}}{\sqrt p}
\ll_h w^h . \tag{1}
$$

The constants may depend on $h$, which is fixed before $p$ grows. Here are the complete-sum and completion details.

**Complete shifted trace sum.** For $2h$ offsets $d_j\in\mathbb F_p$, if some offset has odd multiplicity, then for all $k\in\mathbb F_p$,

$$
\left|\sum_{x\in\mathbb F_p}e^{2\pi ikx/p}\prod_{j=1}^{2h}t(x+d_j)\right|
\ll_h \sqrt p. \tag{2}
$$

Let $u_1,\ldots,u_l$ be the distinct offsets. For sufficiently large $p$ depending only on $h$, the square classes of $g(X+u_j)$ are independent over $\overline{\mathbb F}_p(X)$. Otherwise the product over a nonempty subset of shifts has even multiplicity at each finite root. Partition the three distinct geometric cubic roots by their additive $\mathbb F_p$-cosets and select a nonempty coset. Within that coset the vanishing of all multiplicities modulo two says that the indicators of two nonempty subsets $B,A\subset\mathbb F_p$, of sizes at most 3 and $2h$, have zero cyclic convolution over $\mathbb F_2$. Here $B$ represents the roots in the coset after a common translation, and $A$ represents the negated shifts. To rule this out, list the at most $2h+3$ residues $v_j$ in these sets. Partition $[0,p)^{2h+3}$ into at most $16^{2h+3}$ boxes, using intervals of length at most $\lceil p/16\rceil$. Two of the $p$ vectors $(av_j\bmod p)_j$ lie in one box. Their distinct multipliers differ by some $d\ne0$, so all $dv_j$ have signed representatives of absolute value $<p/8$. Every sum of one representative from each set lies in $(-p/4,p/4)$. Now cyclic convolution agrees with ordinary Laurent polynomial multiplication on these representatives, and the latter has no zero divisors over $\mathbb F_2$.

Consequently the smooth projective normalization $C$ of $y_j^2=g(x+u_j)$ is geometrically connected. The $x$-map has degree $2^l$, is tamely ramified only at the at most $3l+1$ roots and infinity, and has genus bounded in terms of $h$ by Riemann–Hurwitz. At unbranched finite $x$, expanding the product of $t$ gives a sum over $C(\mathbb F_p)$ and at most $2^{2h-l}$ relative sign choices for occurrences of repeated offsets. Each sign choice gives a character $\omega(\Phi)$, where $P_j=(x+u_j,y_j)$ and $\Phi=\sum_j e_jP_j:C\to E$. In an odd-multiplicity block $e_j$ is a nonzero integer of magnitude at most $2h$. Changing the sign of $y_j$ while fixing all other coordinates changes $\Phi$ by $2e_jP_j$. Hence $\Phi$ is nonconstant. Its degree is bounded in terms of $h$: each $P_j$ has degree $2^{l-1}$, multiplication by $e_j$ scales degree by $e_j^2$, and the elliptic degree parallelogram law bounds the degree of the sum. All maps extend over the smooth projective curves. Branch fibers omitted in this expansion contribute $O_h(1)$.

The Lang isogeny $\mathrm{Frob}_p-1:E\to E$ is finite étale of degree $N$ with geometric deck group $G$. Its quotient by $\ker\omega$ is a cyclic étale cover of degree $\operatorname{ord}(\omega)$, and the induced faithful character supplies an unramified rank-one sheaf with trace $\omega(P)$, using the inverse character if required by the trace convention. The pullback by $\Phi$ is geometrically nontrivial. A trivial pullback would split this cyclic cover on the geometric curve, giving a lift of $\Phi$ through the cover. Its degree exceeds $\deg\Phi>0$, ruling out a lift. For $k\ne0$ also tensor by the Artin–Schreier sheaf for $kx$. Its poles on $C$ have positive orders at most $2^l<p$ with total order $2^l$. Hence the twist has wild ramification and bounded total Swan conductor, while the Lang factor is unramified. In either case the trace formula for a geometrically nontrivial rank-one sheaf and the Riemann hypothesis for curves bound the sum by $O_h(\sqrt p)$, since the genus, removed points and total conductor are bounded in terms of $h$. This proves (2). [Farashahi–Fouque–Shparlinski–Tibouchi–Voloch, §4](https://www.normalesup.org/~tibouchi/papers/welldistributed.pdf) · [Deligne, Weil II](https://www.numdam.org/articles/10.1007/BF02684780/)

**Completion of the actual bit pattern.** Let $B$ be the $w$ nonnegative words supported on the free bit positions, and $z$ a uniform word among the $M/w$ words with zeros in those positions. The accepted curve abscissas are $z+d<p$, $d\in B$. In the $2h$-th moment expansion, the $O_h(w^h)$ tuples with all multiplicities even cost $O_h(w^h)$ in total. For any remaining tuple with $\max d_j\ge p$ the term vanishes identically. Otherwise the distinct offsets remain distinct modulo $p$ and the truncation is the single interval $0\le z<p-\max d_j$. For every run with bit positions $v,\ldots,v+t-1$, the mask on $z$ is $z\bmod 2^{v+t}<2^v$. The normalized Fourier coefficients of an interval in this cyclic group have absolute sum at most $Cb$ by the geometric-series estimate. Multiplying the $s$ expansions costs at most $(Cb)^s$ in absolute coefficient sum. Each resulting phase $e^{2\pi i\beta z}$ on the truncated interval has normalized Fourier coefficient sum at most $Cb$ when completed on $\mathbb F_p$, uniformly in real $\beta$. This follows by summing $\min(1,[2p\|\beta-j/p\|]^{-1})$ over $j$ with the natural value $1$ at a vanishing distance. Apply (2), then divide by the *actual* number of fixed words $M/w$. Each unpaired tuple contributes at most $O_h(w(Cb)^{s+1}/\sqrt p)$. Summing all terms proves the first bound of (1). The second follows from $w^{h+1}(Cb)^{s+1}=2^{o(b)}$ and $p=2^{b+O(1)}$.

## From slot moments to witness correlation

Let $D=\prod_i|F_i|$. For distinct group characters the off-diagonal entry of their challenge Gram matrix is

$$
\mathbb E_r Z_\xi(r)\overline{Z_\zeta(r)}
=\prod_i\sum_{P\in F_i}(\xi\bar\zeta)(P).
$$

Its diagonal is $D$. The Walsh signs impose equality of each *entire* address including its sign. This is why the group-character sum in each factor of the Gram entry is unsigned, even though the features themselves are signed.

Choose positive $\eta,\gamma$ with $\theta<\gamma$ and $\eta+\gamma<1/2$. Then choose fixed $h$ satisfying $h(1-2\eta-2\gamma)>1$. For each character of order exceeding $B_h$, independence, (1), and the trivial $2w$ anchor bound give a $2h$-th moment at most

$$
(2w)^{2h}(C_h w^h)^{m-1}=2^{hb+o(b)}
$$

for its product across slots. Markov and a union bound over at most $N$ characters show that the probability any such product exceeds $2^{(1-\eta-\gamma)b}$ is at most $2^{[1-h+2h(\eta+\gamma)]b+o(b)}$. Conditioning on $U$ multiplies this by at most $q^{-1}=2^{o(b)}$.

For arbitrary curve order we must also handle low-order ratios. Set $J=\mathrm{lcm}(1,\ldots,\lceil B_h\rceil)$, let $H=\widehat G[J]$, and let $C_T=T+H^\perp$ be the target's coset in $G$. The elliptic group has at most two invariant factors, so $|H|\le J^2=O_h(1)$. Every character outside $H$ has order greater than $B_h$. In the address-Walsh expansion of $f$ retain only tuples with sum in $C_T$. This orthogonal projection contracts the norm and preserves the witness coefficient. On these addresses the modes $\xi$ and $\xi\alpha$ for $\alpha\in H$ are proportional by the constant $\alpha(T)$, so combine their coefficients and retain at most $k=|A|$ representatives of distinct cosets modulo $H$. Let $D_T$ count the retained tuples. Its value is the diagonal of the restricted Gram matrix. For any pair of distinct representatives of ratio $\omega$, its off-diagonal is the following average by character orthogonality,

$$
\frac1{|H|}\sum_{\alpha\in H}\overline{\alpha(T)}
             \prod_i\sum_{P\in F_i}(\omega\alpha)(P).
$$

All $\omega\alpha$ lie outside $H$, so the preceding simultaneous estimate bounds this average in absolute value.

There is a separate targetwise estimate for a small $D_T$. Write $\bar T$ for the target class in $G/H^\perp$. If $U$ occurs, the witness anchor lies in a class $j$ containing at least one anchor point. Thus its complementary tail class $\bar T-j$ contains at most $D_T$ tail tuples. Freeze the tails. For each $j$ whose complementary class has fewer than $B$ tail tuples, fewer than $B$ possible anchor points can give a target hit with that class. Since there are $|H|$ classes, $\Pr(U,D_T<B)\le |H|\rho B$. Taking $B=2^{(1-\eta)b}$ gives $2^{-\eta b+o(b)}$. This proof uses only the anchor marginal.

Outside the two exceptional events, each restricted off-diagonal has modulus at most $D_T2^{-\gamma b}$. The smallest eigenvalue on the retained representatives is at least $D_T(1-(k-1)2^{-\gamma b})$ by the Hermitian row-sum bound. Writing $c'_\xi$ for the combined coefficients, the contracted norm bound therefore gives

$$
\left|\mathbb E_r f(r)(-1)^{r\cdot a_*}\right|
=\left|\sum_{\text{representatives }\xi}c'_\xi\xi(T)\right|
\le \sqrt{\frac{k}{D_T(1-(k-1)2^{-\gamma b})}} .
$$

This is exponentially small since $\theta<1-\eta$. On exceptional singleton roots Cauchy–Schwarz bounds the coefficient by one. Combining the estimates proves the statement. The argument gives a restriction on the explicit group-character representation over arbitrary sufficiently large prime fields. It leaves open a joint arithmetic evaluator, challenge-dependent mixing, and nonlinear readouts. [Numerical identity checks](../../checks/check_direct_projection.py) · [Output](../../evidence/check_direct_projection.txt)
