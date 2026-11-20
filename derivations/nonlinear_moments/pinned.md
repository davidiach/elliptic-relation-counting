# Guessing all ordinate signs before a nonlinear readout

Guess one ordinate sign $\eta_i$ uniformly and independently for every slot, keeping points with that sign at every occupied abscissa. On an original guarded singleton the full sign vector is correct with probability $2^{-m}$. On a match, one needs to recover only the concatenated abscissa address. This prevents the cancellation between opposite points that split the first determinant into even- and odd-column blocks. It costs $O(m)$ in the logarithm of the expected number of attempts. [Original sign-slice argument](../sampled_slots/signed_moment.md)

Let $\xi$ be the challenge on abscissa addresses, $s$ the challenge on ordinate signs, and $Z_j(\xi,s)$ the original signed power moment. The projected moment $M_j$ on the pinned points satisfies

$$
 M_j^\eta(\xi)=2^{-m}\sum_{s\in\{0,1\}^m}(-1)^{s\cdot\eta}Z_j(\xi,s)\quad \text{in }\mathbb F_p.
$$

Orthogonality of binary characters selects exactly those points with the prescribed signs. A nonlinear score applied after this field addition thus uses a joint combination of original signed measurements. It can be computed directly on the restricted slots. For orders 0, 1 and 2 the product of signed slot masses, determinant of weighted row sums, and Gram subset formula respectively give per-query cost $(mw+2^m)\operatorname{poly}(b,m)$ including constructing the slots. A successful parity correlation $\delta\ge2^{-O(b/m)}$ on matched, isolated signs would still meet the project's targeted cost through $m=O(\sqrt b)$ with the inherited extraction and outer reduction; it remains necessary to establish that correlation at growing widths and arities. [Moment evaluation](../joint_parity/target_incidence.tex) · [Address recovery](../sampled_slots/address_parity.md)

## An affine symmetry of the score

Fix a nonempty pinned grid. In each slot choose an occupied binary address $a_i$ and let $V_i$ be the binary linear span of differences (XORs) of occupied addresses. Set $a=(a_i)_i$ and $V=\bigoplus_i V_i$. Every legal address lies in $a+V$. For any $k\in V^\perp$, changing $\xi$ to $\xi+k$ multiplies every moment by the same sign $(-1)^{k\cdot a}$, since this factor is constant on all summands.

Suppose $a\notin V$. Then some such $k$ has $k\cdot a=1$. If a bounded readout satisfies $\phi(-z)=\phi(z)$ on simultaneous negation of the vector of field moments, its score at $\xi+k$ is unchanged, while the Walsh sign at **every** legal frequency flips. Its coefficient is therefore exactly zero at every legal address. For an arbitrary readout only its odd part $(\phi(z)-\phi(-z))/2$ can contribute at those frequencies. The discriminant character $\chi(M_0 M_2-M_1^2)$ is always even. The individual characters of $M_1,M_2$ are even when $p\equiv1\bmod4$. The mixed character $\chi(M_2+M_1^2)$ has both parts.

There is a cheap remedy for this particular symmetry. Choose one known legal address $c$ and replace every address $v$ by $v+c$. The translated occupied set has affine hull $V$ and includes zero. Each moment acquires the same known challenge sign $(-1)^{\xi\cdot c}$. In the experiment $c$ was chosen uniformly among the legal addresses on the root, independently of which one is marked, and the same translation was applied to the witness and every control. The symmetry is a statement about the chosen coordinates; it places no general restriction on all encodings. It also applies to an unpinned grid when its full signed addresses obey the same affine-hull condition.

Write $f_e(\xi)=(\phi(M(\xi))+\phi(-M(\xi)))/2$ and $f_o(\xi)=(\phi(M(\xi))-\phi(-M(\xi)))/2$. The recentered coefficient at the new address of an old tuple $v$ is exactly $\widehat f_e(v+c)+\widehat f_o(v)$, by linearity and the Walsh frequency-shift identity. Its magnitude must be taken **after** addition. If the occupied set itself is an affine coset and the readout is purely even, the distribution of this magnitude over a uniform occupied origin $c$ is identical for every legal $v$: $v+c$ is uniform on $V$. For every nonempty width-two pinned grid the occupied set is such a coset. For arbitrary occupied sets the finite test also averages the actual absolute value over every legal origin on each seed.

## Trials with real alternatives

The trials plant a signed relation, then require exactly one guarded relation in the **original** signed grid, distinct abscissas across slots, and no abscissa shared with the endpoint. Ordinary planting weights a seed by the number of witnesses, which equals one on these selected roots. The matching $\eta$ used for evaluation therefore samples the correct law conditional on a regular singleton and a successful independent sign guess. The sign-guess success factor $2^{-m}$ must be charged in an application. At each root we evaluated the full abscissa Walsh spectrum, computed the absolute coefficient at the correct address, and subtracted the mean at every *other legal address on that pinned root*. A tied rank receives half credit. The three trials have independent sampling seeds. The two field-65521 trials also use distinct targets.

| Field and arity | Target | Roots with alternatives | Median legal pinned tuples | Mixed-character gap (paired SE) | Rank against wrong legal tuple |
|---|---|---:|---:|---:|---:|
| $p=1031$, $m=4$, $w=4$ | $G$ | 599 | 18 | $-0.01241$ (0.00610) | 0.5020 |
| $p=65521$, $m=6$, $w=4$ | $G$ | 600 | 144 | $+0.00088$ (0.00175) | 0.4969 |
| $p=65521$, $m=6$, $w=4$ | $2G$ | 600 | 144 | $+0.00132$ (0.00217) | 0.4985 |

The mixed character is $\chi(M_2+M_1^2)$. Requiring at least two occupied abscissas **in every slot** left 258 and 227 roots in the six-slot trials, with median sizes 288 and 256. Their mixed-character gaps were $-0.00094$ (0.00146) and $+0.00007$ (0.00150). In those original coordinates 540 and 546 of all 600 sampled roots had an occupied affine hull excluding zero. The checker confirms the exact annihilation for even scores on every such root.

With one random legal origin per root, the six-slot mixed-character gaps were $+0.00386$ (0.00255) and $+0.00394$ (0.00324). The recentered discriminant gap was $+0.01147$ (0.00540) at $G$ and $-0.00585$ (0.00360) at $2G$. Exact averaging over all legal origins on the same roots gave discriminant gaps $-0.00023$ (0.00025) and $+0.00016$ (0.00017), and mixed-character gaps $+0.00134$ (0.00128) and $+0.00170$ (0.00166). These standard errors are over independently sampled roots. The tests provide no replicated distinction strong enough to extrapolate a recovery cost. A small paired gap alone gives neither an asymptotic lower nor upper bound on the absolute witness coefficient. The exact projection and the symmetry hold independently of those empirical observations. [Executable probe](../../checks/probe_pinned_moments.py) · [Complete results](../../evidence/probe_pinned_moments.txt)
