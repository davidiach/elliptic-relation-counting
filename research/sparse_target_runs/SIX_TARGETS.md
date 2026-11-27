# A nonconstant vector in six consecutive pinned kernels

On $y^2=x^3+x+46$ over $\mathbb F_{101}$, the group order is $r=109$.
Let $B=(2,37)$ and $P=20B=(4,66)$. A reduced-floor seed with $m=6$ has
block sizes $2,2,4,6,8,8$ and threshold $T=4$. Its 6,144 half-tuples
represent every group element. The collision count is 348,112, giving

$$
\Phi=\frac{2674241382816338354696941832713}{2640492748413577099820018958336}
<\frac{64}{57}.
$$

Each of $M_0,M_P,\ldots,M_{5P}$ has nullity 7. On 30 tagged columns,
each raw tuple-incidence matrix has rank 24 and each pinned matrix has
rank 23. The six stacked matrices have nullity **6**, one more than the
five dimensions of block constants whose block values sum to zero.
At $6P$ the individual nullity is still 7, and at $7P$ it is 6.
[Exact affine masks](witnesses.json) · [Verification transcript](verification.txt)

## The shared vector

Before translating the first block, the first two tagged blocks are the
same set $F=\{a,b\}$. Put values $(1,-1)$ on the first copy, $(-1,1)$
on the second, and zero elsewhere. Call this vector $v$. Let $C(u)$ count
the completions to sum $u$ after removing two of the four slots in these
copies. The remaining slots consist of two occurrences of $F$ and two
occurrences of each other block.

For a row pinned in any other block the two contributions cancel by
exchanging the copies. On the first two blocks put $K_{xy}=C(t-x-y)$ for
$x,y\in F$. The first diagonal block of $M_t$ is
$\operatorname{diag}(K\mathbf1)+K$ and the cross block is $2K$.
On opposite vectors on the two copies the action therefore reduces to

$$
\operatorname{diag}(K\mathbf1)-K
= C(t-a-b)\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
$$

Consequently $M_t v=2C(t-a-b)v$ over the integers. After translating
the first block by $H$, carry the labels of $v$ with it and replace $t$
by $t-2H$ on the right. This computation applies at every target.

In the example $F=\{-B,B\}$ and $H=25B$. The two tagged blocks after
translation are $\{24B,26B\}$ and $\{-B,B\}$. Assigning $(1,-1)$ and
$(-1,1)$ in those respective orders gives the shared vector. The six
completion counts, in target order, are

$$
109\, (793,792,793,793,792,793).
$$

The scalar vanishes modulo 109 throughout the run. The verifier also
checks the integer eigenvector identity at all 109 targets and counts
the entries of the six pinned matrices by direct enumeration.

On $y^2=x^3+x+4$ over $\mathbb F_{251}$, the same calculation gives a
second shared witness with group order 271, $P=(40,87)$, and $\Phi =
1.0541007592\ldots$. All seven blocks have sizes between 2 and 8.
The six individual nullities are 8 and the stacked nullity is 7, one
above the forced dimension 6. All pin counts are nonzero modulo 271
at all six targets. The completion counts are 699451, 714085, 701348,
701348, 714085, and 699451, each divisible by 271. At $6P$ the
individual nullity returns to 7.

## What the comparison says

On $y^2=x^3+x+13$ over $\mathbb F_{59}$, a different accepted seed has
$m=4$, group order 67, and generator $P=(18,9)$. The six individual
nullities are all 5. All other targets have nullity 4, and all 2,211
unordered pairs of distinct targets have joint nullity 3. Every pin count
at the six exceptional targets is nonzero modulo 67. Thus individual
run length alone does not determine the rank of the stack.

The two-point construction also has a limit under LZR's original parameters.
For the full prime-order elliptic group set $q=(r-1)T/p$ and $\mu=q^m$.
Acceptance implies $W>(1-1/\sqrt{57})\mu$ from the nonnegative
collision term and the penalty $64(W-\mu)^2/\mu^2$ in $\Phi$.
Each parity class supplies at most $T$ elements to a block because the
abscissa is injective within the class and the affine multiplier is nonzero.
Since $W=\prod_i |F_i|$, every block obeys

$$
|F_i| > (1-1/\sqrt{57})\,q\left(\frac{r-1}{2p}\right)^{m-1}.
$$

In the original general regime $q\ge(16r)^{1/m}$. Hasse's bound gives
$r/p\to1$, and the published arity satisfies $m^2=o(\log r)$.
The logarithm of the displayed lower bound is at least
$(\log r)/m-O(m)$ and tends to infinity. Hence a two-element block is
excluded on all sufficiently large full-parameter curves. Whether a
shared mode survives on blocks of the required growing size remains
open. [Liu, Zhang and Ren, §3, equations (3.1)–(3.2)](https://eprint.iacr.org/2026/3020)

For two identical larger tagged blocks $F$ the same derivation gives
$L_t=\operatorname{diag}(K_t\mathbf1)-K_t$ on opposite block values,
with $(K_t)_{xy}=C(t-x-y)$ for $x,y\in F$. The concrete question is
whether several such $L_t$ can share a nonconstant null vector when
the blocks grow at the required rate and come from the canonical seed.
Increasing the length of the reduced-floor two-point run leaves that
question open.
