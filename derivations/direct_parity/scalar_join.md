# One scalar equality, with a high-rank intermediate

At a fixed target the signed Fermat residual can be read from one equality between two scalars. At the thinned widths a random projection preserves every legal equality simultaneously with probability at least $1-(2w)^m/(p+1)$. This eliminates the determinant values and their recurrence. A tempting strategy would keep low-rank factors while raising the scalar difference to the Fermat power. We show why that particular strategy costs too much even on actual successful seeds.

## Project the matching partial sums

Fix a regular root: the slot abscissas are pairwise disjoint and avoid the target abscissa. On this grid $\Delta(v)=0$ exactly when the signed group sum equals $T$. Split the slots into $I,J$ of sizes $\lfloor m/2\rfloor,\lceil m/2\rceil$. The left group point of a half tuple is $P=\sum_{i\in I}P_i$, and the right point is $Q=T-\sum_{j\in J}P_j$. Use distinct affine coordinates $\iota(P)=(x(P),y(P))$ and an off-curve coordinate pair for $O$; one can choose the latter by trying three ordinates at $x=0$. This encodes every curve point injectively.

Choose a uniform direction in $\mathbb P^1(\mathbb F_p)$, giving one of the linear forms $L(x,y)=x+cy$ or $L(x,y)=y$. Two distinct coordinate pairs agree in exactly one of the $p+1$ directions. There are at most $D_{\rm grid}=(2w)^m$ pairs of half tuples. With probability at least $1-D_{\rm grid}/(p+1)$ there is no false cross-half match. This is a per-root, per-target bound that applies simultaneously to **all** challenges and address pins. The thinned sampler gives $D_{\rm grid}/p=O(m^{-2})$. Full group equality can also be verified directly when the halves have been listed. [Widths and isolation](../sampled_slots/isolation.md)

Let $a=L(\iota(P))$, $b'=L(\iota(Q))$ and aggregate the signed half-tuple weights by scalar label into $u_a$ and $v_{b'}$. On a good direction we have the field identity

$$
M_0(r)-M_{p-1}(r)=\sum_{a,b'}u_a v_{b'}
\bigl[1-(a-b')^{p-1}\bigr]=\sum_a u_a v_a .
$$

On a guarded singleton this equals the sign of its address parity. With the halves explicitly aggregated, a table join computes it in $(2w)^{\lceil m/2\rceil}\operatorname{poly}(b,m)$ bit time by sorting. This gives log cost at most $b/2+O(b/m+m+\log b)$, and involves no recurrence or enumeration of determinant values. It still materializes potentially large tables. The indicator matrix indexed by distinct labels has **rank one** when there is just one true matching pair. Its complementary Fermat-power matrix then has rank at most two.

## The middle powers cannot all have small separation rank

Let $A,B\subset\mathbb F_p$ be the surviving distinct labels, $s=\min(|A|,|B|)$. For any univariate polynomial $g$ of exact degree $0\le j<s$ (hence $j<p$),

$$
\operatorname{rank}_{\mathbb F_p}\bigl(g(a-b')\bigr)_{A\times B}=j+1.
$$

Factor the matrix as two evaluation matrices on the monomials $1,X,\ldots,X^j$, with the coefficient matrix of $g(X-Y)$ between them. Both evaluation matrices have full column rank, by univariate interpolation. The middle matrix is zero beyond its antidiagonal, with nonzero antidiagonal entries from the leading coefficient times the signed binomial coefficients $\binom{j}{k}$. Its rank is $j+1$.

Consequently every forward addition chain for raising $a-b'$ to $p-1$ has an intermediate of rank at least $(s-1)/2$ for $s\ge4$. More generally, consider a straight-line computation of **any polynomial in the scalar difference** using binary field additions, subtractions and multiplications and arbitrary constants shared across the grid. Assume its output exactly detects equality on these lists and there is one matching label. On the column at that label its output is one in the matching row and zero at all the other $|A|-1$ differences, forcing degree at least $|A|-1$. Follow the first gate whose univariate polynomial has degree at least $\lceil(s-1)/2\rceil$. Its inputs have smaller degree, so its output still has degree $<s$. The rank formula applies. The same argument using the row holds for $|B|-1$. Seed- or target-specific scalar coefficients do not change the proof.

An attempted shortcut using the entrywise inverse with zero sent to zero also has a large intermediate. Delete the common label from one list; the remaining rectangular reciprocal matrix $1/(a-b')$ on two disjoint lists has rank at least $s-1$. In a square submatrix of size $t$, a vanishing linear combination of columns gives a rational function with $t$ distinct zeros, denominator of degree $t$ and numerator of degree at most $t-1$. The numerator vanishes identically and the simple-pole coefficients must all be zero. These are limits on explicitly expanding factors of those intermediate matrices. An implicit formula for a matrix, a proof system, a computation of its weighted sum alone, or a circuit accessing the separate slot seeds need not expand them.

## Size on the mixed sampler, including signed challenges

Write $N=\#E(\mathbb F_p)$, $b=\lceil\log_2p\rceil$, $M=2^b$, $\nu=Nw/M$ and $q=\Pr(U)$ for the unique guarded event at the prescribed $T$ with $2T\ne O$. At the thinned widths $\log_2\nu=b/m+O(1)$ and $q^{-1}=2^{O(m+\log m)}$. Take growing $m=O(\sqrt b)$ and a layout with $s_0=o(b)$ runs of free bits in each fixed-bit slot satisfying the isolation conditions. Conditional on $U$, all slot sizes exceed $\nu/m^3$ except with probability $O(m^{-2})$: each slot on the event its size is $\le h$ permits at most $h\lambda/\nu$ expected completions, and $q\asymp\lambda$. Regularity fails with conditional probability $o(m^{-2})$ by the pointwise marginals.

For either balanced block of size $k$ let $H_k(P)$ count its tuples of group sum $P$. The off-diagonal part of the displacement count gives

$$
\mathbb E\sum_P H_k(P)(H_k(P)-1)
\ \le\ C^k k\,2^{s_0}\nu^{2k}/N=2^{o(b)}
$$

in the sense of an upper bound, since $k\le(m+1)/2$. In that estimate the diagonal is exactly the all-zero displacement vector. Dividing by $q$ and using Markov leaves at most $2^{o(b)}$ ordered collisions in each block except on $O(m^{-2})$ of the successful roots. Each block contains at least $(\nu/m^3)^k$ tuples there. Every tuple in a cell of size at least two can be charged to an ordered collision, so there are $2^{b/2-o(b)}$ tuples with unique group sum in each block. [Displacement energy calculation](../sampled_slots/prefix_dictionary.md)

For each fixed root the expected number of ordered within-block tuple pairs of **different** group sums that collide after the projection is at most $(2w)^{2k}/(p+1)=2^{o(b)}$ as an upper bound. A second Markov estimate loses only another $O(m^{-2})$ in the probability over directions. Thus on almost all regular singleton roots and directions there are at least $2^{b/2-o(b)}$ scalar labels unique within each half. The signed aggregate at each such label is $+1$ or $-1$ **for every challenge**, so cancellations do not defeat the rank bound by removing most of the support. Restricting to nonzero signed aggregates keeps the true intersection on a singleton. Multiplication of rows and columns by nonzero weights preserves rank.

The expensive intermediate therefore occurs in this explicit separated-power approach even though the requested final bit is one sign and the equality matrix has rank one. A different way of extracting the common label from the seeds remains open. [Executable finite checks](../../checks/check_scalar_join.py)
