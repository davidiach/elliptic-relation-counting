# A unary chord function and its isogeny trace

The chord statistic along an addition tree can be replaced by one function of its final sum and one term per original point. The function takes $O(\log N)$ field operations per supplied point when its ambient finite group has known order $N$ prime to $p$. Its trace over a whole isogeny-kernel coset takes just one evaluation on the quotient. We give a fresh-seed estimate for the chance that any sampled slot admits even one such compression.

## Removing the internal tree

Use a short Weierstrass curve over $\mathbb F_p$, $p>3$, and write $l(P,Q)$ for the chord slope, with the tangent at a doubling and zero for opposite points or identity inputs. Tree independence of the summed slopes gives

$$
l(P,Q)+l(P+Q,R)=l(Q,R)+l(P,Q+R).
$$

This includes exceptional partial sums. On a finite subgroup $G$ of order $N$ prime to $p$, summing this identity over $R\in G$ gives

$$
g(P)=N^{-1}\sum_{R\in G}l(P,R),\qquad
l(P,Q)=g(P)+g(Q)-g(P+Q).
$$

The function $g$ is odd and $g(O)=0$. Its normalization is unique, since the difference between two solutions would be a homomorphism from $G$ to the additive field. Thus the sum $C$ of the slopes on *any* addition tree for $P_1,\ldots,P_t$ is $\sum_i g(P_i)-g(\sum_iP_i)$. [Tree identity](../product_differences/chords.md)

There is no need to sum over $N$ points. For any known positive annihilator $n$ of $P$ prime to $p$, follow a binary addition chain for $[n]P$, maintaining the sum of slopes $c_j$ for the subtree on $j$ copies of $P$. Combining $i,j$ gives $c_{i+j}=c_i+c_j+l([i]P,[j]P)$, with $c_1=0$. The telescoping formula gives $c_n=n g(P)$, so $g(P)$ costs $O(\log n)$ field operations. One may take $n=N$ on this subgroup. The construction is the first normalized Laurent coefficient of a Miller function. The prime-to-$p$ condition is essential to this division. [Miller, *Short Programs for Functions on Curves*](https://crypto.stanford.edu/miller/miller.pdf)

Include the endpoint in the slots and let $s_i(P)$ be its specified speed or the specified speed of a sampled point. With $S=\sum_iP_i$ and $k=\sum_i s_i(P_i)$, the tree contribution to the empty-fiber digit is now

$$
\sum_{(P_i)}k\sum_i g(P_i)\; -\;\sum_{(P_i)} k\,g(S).
$$

The first term factors into sums over at most two individual slots by expanding $k$. The separate pair-of-rows correction involving $K(P,Q)$ is already computable by the small-slot multipoint calculation. The sole remaining many-slot term is a marked group convolution paired with the unary function $g$. Explicitly, convolve the group-ring elements $\sum_{P\in F_i}(1+\epsilon s_i(P))[P]$, with $\epsilon^2=0$, and pair their $\epsilon$ coefficients with $g$. In a full group table the unmarked coefficient at $O$ also gives the hit count modulo $p$, which is exact at the thinned widths; the marked coefficient there gives the singleton speed. For example, successive small-slot convolutions take $O(mwN\operatorname{poly}(b,m))$ bit time and $O(N\operatorname{poly}(b,m))$ space using group points as keys. This eliminates the half lists but costs more than the balanced method. No $g$ table is needed for individual values.

## Exact compression of a complete coset

Let $\phi:E\to E'$ be a separable isogeny between short Weierstrass curves. Suppose its geometric kernel $H$, of size $e$, lies in the chosen prime-to-$p$ finite subgroup $G$. Write $\omega=dx/(2y)$, $\omega'=dx'/(2y')$, and $\phi^*\omega'=\alpha\omega$. With $g'$ computed on $\phi G$, one has

$$
\sum_{h\in H}g(P+h)=\alpha\,g'(\phi P).
$$

To prove it, first use the rational chord on both curves. As a function of $P$ for generic $Q$, $l(P,Q)\omega_P$ has simple poles at $O,-Q$ of residues $-1,+1$. Thus the following two rational functions have equal principal parts in $P$:

$$
\sum_{h,h'\in H} l(P+h,Q+h'),\qquad e\alpha\, l'(\phi P,\phi Q).
$$

For the sum, each pole over $O$ or $-\phi Q$ has residue $-e$ or $e$ after multiplying by $\omega$. For the right side, separability and the displayed differential pullback give the same residues at the same simple poles. Their difference is independent of $P$. By symmetry it is independent of $Q$, and oddness under simultaneous negation makes that constant zero. This proves a rational identity. On points with $\phi P,\phi Q,\phi(P+Q)$ all nonidentity, it holds with the slope and tangent conventions above. Applying the coboundary identities on both sides shows that the difference between $P\mapsto\sum_h g(P+h)$ and $P\mapsto\alpha g'(\phi P)$ is additive on these pairs. For the exceptional pairs additivity also holds, since that difference is odd, constant on $H$-cosets and zero on $H$. A homomorphism from $G$ to the additive field is zero.

For multiplication by $d$ prime to $p$, if $G$ contains all $E[d]$, the $d^2$-term trace is $d\,g([d]P)$ and can be computed in $O(\log d+\log N)$ field operations given $P,d,N$. For a supplied general isogeny the cost includes evaluating $\phi$. In the marked convolution, a full slot coset on which its speed is constant can be represented by one quotient entry: for each fixed selection from the other slots the weighted $g$ sum uses the trace. For the separate hit count, a quotient match for a chosen slot coset and an ordered selection of the other labels has exactly one lift from that coset summing to the original target. Labels from other slots must retain their multiplicities and individual speeds. Local corrections in the original curve still count their original points. This is a genuine saving for cosets given compactly, and applies to full cosets split out of a larger slot by linearity.

## How often do sampled slots contain one?

Let $G=E(\mathbb F_p)$, $N=|G|$, $s$ the maximum number of free-bit runs and $\nu=Nw/2^b$. For a fixed nonzero displacement $W$ and one mixed slot $F_i$, the known displacement estimate gives $\sum_P\Pr(P,P+W\in F_i)\le C 2^s\nu^2/N$. For any chosen subgroup of prime order $\ell$, each full coset supplies $\ell$ distinct starting points for such pairs at a fixed $W\ne O$ in the subgroup. The number of order-$\ell$ subgroups of $G$ is at most $\ell+1$, from the rank-at-most-two $\ell$-torsion. The number of distinct primes dividing $N$ is $O(b)$. Every nontrivial finite subgroup contains a prime-order subgroup. A union bound therefore gives

$$
\Pr(\text{any sampled slot contains any nontrivial subgroup coset})
\ \le\ O(mb\,2^s\nu^2/N).
$$

This bound permits choosing the isogeny and coset *after* reading all the seeds. At the balanced widths its exponent is $-b+2b/m+o(b)$, even after conditioning on the guarded singleton event of probability $2^{-o(b)}$. [Displacement bound](../sampled_slots/prefix_dictionary.md) · [Singleton sampling law](../sampled_slots/isolation.md)

An estimate for partial cosets also prices subtracting missing points from a full trace. For $0<\eta\le1$, a coset of a subgroup $H$ of order $e$ containing $t\ge\max(2,\eta e)$ retained points supplies $t(t-1)\ge\eta^2e^2/2$ ordered pairs at nonzero displacements in $H$. For fixed $H$ the expected number of such pairs in a slot is at most $(e-1)C2^s\nu^2/N$. The number of subgroups of order $e$ is at most $\sum_{d\mid e}d\le e\tau(e)$, where $\tau$ counts divisors. Indeed taking annihilators gives distinct index-$e$ subgroups of the character group, generated by at most two elements; taking preimages in $\mathbb Z^2$ and counting triangular bases gives the divisor sum. Summing over $e\mid N$, the chance that *any* original slot contains such a coset intersection is at most $O(m\,2^s\nu^2\tau(N)^2/(\eta^2N))$. Since $\tau(N)=N^{o(1)}$, even $\eta=2^{-o(b)}$ preserves the exponential rarity at the balanced widths.

## Coset density on the half sums themselves

For a balanced block $I$ of $k$ slots put $H_I(U)=\#\{\text{tuples in }I\text{ summing to }U\}$, and let $S_I$ be its support. For each fixed $W\ne O$ the same displacement argument gives

$$
\mathbb E\sum_U H_I(U)H_I(U+W)\le C^k k\,2^s\nu^{2k}/N=:B_I.
$$

Indeed in the sum over slot differences adding to $W$ at least one difference is nonzero. Mark its slot; the other slot differences force its value. Use the pointwise bound there and sum the total masses $\mathbb E|F_i|^2\le(2w)^2$ at the other slots. At balanced thinned widths $B_I\le2^{o(b)}$, including with an extra fixed endpoint translation.

For a fixed subgroup $H$ of size $e$, the expected number of **distinct half sums** in cosets meeting $S_I$ in $t\ge\max(2,\eta e)$ points is at most $2B_I/\eta$. To see this, each such cell contributes $t\le 2t(t-1)/(\eta e)$ points and its pairs are counted, at least once, among the $e-1$ possible nonzero displacements. For any menu of $J$ subgroups, the expected number of points covered by all such cells is at most $2J B_I/\eta$, and conditioning on a guarded singleton costs at most $1/\Pr(U)$. If $J,1/\eta=2^{o(b)}$, Markov's inequality bounds the conditional probability of covering more than $2^{\varepsilon b}$ points by $2^{-\varepsilon b+o(b)}$ for each fixed $\varepsilon>0$. Typical singleton half supports have $2^{b/2-o(b)}$ distinct sums. This permits adaptive selection from the whole menu. It includes all subgroups of subexponential order by the divisor estimate above, and on cyclic rational groups it includes *all* subgroups, since their number is $\tau(N)$. On prime-order rational groups any isogeny is injective on that group or kills it. These estimates do not charge algorithms computing partial traces by a different symbolic identity. [Typical half supports](../direct_parity/scalar_join.md)

The trace identity removes all internal tree corrections and gives an exact speedup on cosets. On the general-prime sampler it supplies no improved end-to-end exponent. The weighted value of $g$ at the summed points still has to be aggregated from the actual seeds. A separate split-block argument treats kernels selected after inspecting the half sums. [Adaptive kernels](adaptive_cosets.md) · [Finite check](../../checks/check_implicit_chord_algebra.py)
