# An average-case equality predictor and its singleton failure

At the thinned widths, the equality of two lifted products can be predicted with small error over an independent uniformly masked target, using the actual sampled slots and no half lists. This includes comparisons of binary slot speeds with any of the endpoint speeds $j=0,\ldots,m$. For $j>0$ a pole at infinity suffices. For $j=0$ a capped sum over the slots either certifies a pole or computes the answer exactly. This reflects how infrequently the masked endpoint is a zero of a particular rational function. On the ordinate-bit task the same predictor makes errors on about half the isolated instances. Comparing every endpoint shift misses the unique equality on almost all such instances.

## A certificate at infinity for positive endpoint shifts

Let $G=E(\mathbb F_p)$, $N=|G|$, and let $F_i$ be the $m$ varying slots, retaining both signs, with pairwise disjoint abscissas. The endpoint is $Q=-T$. Write $D=\prod_i|F_i|$. For each slot choose binary speeds $c_i(P)$ using the seeds and an optional public challenge, independently of the subsequent uniform target. A tuple $a$ has sum $S(a)$ and total speed $k(a)=\sum_i c_i(P_i)$. Compare its slot-motion product $V_C$ with the endpoint-motion product $V_j$ at speed $j\in[0,m]$. The two products are over $\mathbb Z/p^2\mathbb Z$ with the same initial lifts. Empty slots give $V_C=V_j=1$. Assume $D>0$ below.

On regular empty target fibers, equality means that a rational function $L_{c,j}(Q)$ is zero. The [chord identity](../product_differences/chords.md) gives, up to a constant independent of $Q$,

$$
L_{c,j}(Q)=
 \sum_{a:S(a)\ne O}(k(a)-j)K(Q,-S(a))
 -\sum_a\sum_i(c_i(P_i)+j)K(Q,P_i)+\mathrm{constant}.
$$

Indeed the fixed-slot terms for $j=0$ give the two sums with coefficients $k(a)$ and $-c_i(P_i)$. For the derivative of the evaluated endpoint, the determinant has zeros at the slot points and at $-S(a)$ when $S(a)\ne O$. Its logarithmic derivative contributes one copy of each corresponding $K$; the comparison subtracts $j$ times this derivative. The only poles of $K(Q,Z)$ are simple and are contained in $Z,O$. In the invariant differential normalization its residues there are $1,-1$. On the whole curve the pole-divisor degree of $L_{c,j}$ is at most

$$M_0=D+\sum_i |F_i|+1.$$

If $j\ge1$, the residue at infinity of $L_{c,j}(Q)\,d x/(2y)$ is

$$j(m+1)D+\sum_{a:S(a)=O}(k(a)-j).$$

Its ordinary integer representative in this formula lies between $jmD$ and $m(j+1)D\le m(m+1)D$. In particular, if $m(m+1)D<p$, the rational function has a genuine pole, for **every nonempty grid** and every such endpoint shift, including interior values of $j$. The bounded-work predictor returns *unequal* for those shifts. On a regular empty fiber it can be wrong at no more than $M_0$ geometric endpoints because a nonzero nonconstant rational function has as many zeros as poles counted with multiplicity.

## The unshifted product and a capped exact computation

For $j=0$ and all speeds zero the products are identical. Otherwise let $A$ be the union of the slot points with speed one. Choose one active slot $i$ and a selected point $P_*\in F_i$. If its complementary slots represent more than $1+|A|$ distinct sums, adding $P_*$ yields a sum $S$ outside $\{O\}\cup(-A)$. Every contributing tuple there has nonnegative weight, and at least one has positive weight. Its positive integer residue at $-S$ is at most $mD<p$. There is no slot subtraction at this point. This certifies a pole and allows the same prediction of inequality.

Compute the complementary sums in successive slots, stopping as soon as a support exceeds $1+|A|$. Its cardinality can never decrease on adding another nonempty slot, since a fixed point of the next slot translates the preceding support injectively. If the cap is never exceeded, finish the calculation exactly. At each group state keep the count $n$, the sum $u$ of the tuple speeds, the sum $c$ of the internal chord sums and the sum $v$ of their products. Appending one point of speed $s$ with chord slope $l$ contributes

$$
(n,\;u+ns,\;c+nl,\;v+sc+lu+nsl)
$$

to its new group state. Use ordinary integers for counts and field elements for the other entries. Fold in the omitted slot and the endpoint. The count at $O$ and its speed distinguish multiple and singleton fibers exactly. On an empty fiber add the original row-pair terms from the chord identity to the sum of the last accumulators. Each row-pair term costs at most the product of the two slot sizes and has no dependence on the other slots beyond their sizes. With $|F_i|\le2w$ and $|A|\le2mw$ all cases take $O(m^2w^2\,\operatorname{poly}(b,m))$ bit time, including generating the sampled slots. The large-support branch stops after at most $2mw+2$ distinct states. The small-support branch may create at most $(1+2mw)2w$ states when incorporating the omitted slot.

For one raw abscissa bit there is also an exact description of the zero rational function in the earlier [abscissa formula](../product_differences/chords.md). Write $A$ for its sign-symmetric selected subset in the moving slot and $W$ for the sumset of the other slots. For nonempty $A$ and tails the rational function $R_i(X)$ is identically zero precisely when

$$W\subseteq\operatorname{Stab}(A)=\{z\in G:A+z=A\}.$$

Indeed all abscissa multiplicities are smaller than $p$. Equality of the finite partial fractions is equality of the integer abscissa multisets, so total mass forces that no selected sum was $O$. Their support and sign symmetry imply $A+z\subseteq A$ for each $z\in W$. Each inclusion is an equality. Conversely translating $A$ by any tail sum then preserves its full point multiset, with the required multiplicities when tails repeat. Test at most $|A|$ candidate translations using $A-P_0$ for a fixed $P_0\in A$. To check $W\subseteq H$ for the resulting subgroup, select one representative per tail slot, check all within-slot differences against $H$, and test the sum of the representatives. This takes $O(mw^2\,\operatorname{poly}(b,m))$ bit time. Empty selected sets or empty tail slots give the zero rational function directly.

## What the masked error guarantees

Let the target be uniform on an independent set of $n$ eligible points, all satisfying $2T\ne O$. For fixed cross-slot regular seeds and speeds, at most $D$ targets have any relation, and at most $2\sum_i |F_i|$ have endpoint abscissa collisions. On the other targets any error on a certified-pole branch must be a zero of the rational function. Thus the failure probability per comparison, counting irregular endpoints as failures, is at most

$$\frac{2D+3\sum_i|F_i|+1}{n}.$$

For all $m+1$ shifts simultaneously sum the individual bounds. The denominator is the actual masking population $n$. A small target subgroup or a prespecified target does not inherit a bound with denominator $N$. Arbitrarily target-dependent choice of speeds is also outside this count.

For the independent thinned mixed roots, cross-slot irregularity costs at most $\binom m2p\rho^2$, with $M=2^b$, $\rho=w/M$. The exact width condition is $D_{\max}/N=\lambda D_0^m\le D_0^{2-m}/(4m^2)$ for $D_{\max}=(2w)^m$, $D_0=2M/N$. Since $M\ge p$ and $N\le p+1+2\sqrt p$, $D_0\ge 2(1-O(p^{-1/2}))$. For growing $m=O(\sqrt b)$ this proves $D_{\max}/N=O(2^{-m}/m^2)$ and $(m^2+2m)D_{\max}<p$ eventually. A per-comparison error bound is therefore

$$
O\!\left(\binom m2p\rho^2 + \frac{N}{n}\frac{2^{-m}}{m^2}
                      +\frac{mw+1}{n}\right).
$$

At a masking population $n=\Theta(N)$ and $m=\Theta((\log b)^{1/4})$ this is $o((\log b)^{-2})$, even after a factor $m+1$ for simultaneous shifts. The large-support branch never constructs half lists. Conditioning only on the much more common empty fibers preserves this order, since the fraction of positive targets is at most $D_{\max}/n$. [Mixed widths and singleton probabilities](../sampled_slots/isolation.md)

## Conditional error on the witness bit

For a concrete separation, move just one prescribed slot at speed equal to its ordinate sign. Use a negation-invariant eligible masking population, such as the nonidentity points of a large odd-order subgroup. On guarded singleton roots condition on $U$. For any fixed target, the chance of $U$ and at most $2w+1$ distinct complementary sums is at most $(2w+1)\rho$. Indeed each distinct complementary sum forces one point of the moving slot, retained with marginal probability at most $\rho$. This argument holds whether or not that slot is the affine anchor. Divide by the targetwise lower bound $q=\Pr(U)\ge2^{-O(m+\log m)}$. The resulting conditional probability is $2^{-b+o(b)}$. Irregularity also has conditional probability $2^{-b+o(b)}$ in the growing regime. A witness with prescribed moving point equal to rational affine two-torsion has probability $O(1/N)$ conditional on $U$: the unconditional expected number of completions over at most three such points is $3\lambda/N$, and $q\asymp\lambda$.

Outside these exceptions there is an active ordinate in the slot because both signs of the non-torsion witness point are retained. The $j=0$ query therefore certifies a pole; every $j\ge1$ query also certifies a pole. All answers are *unequal*. In truth, on each singleton exactly one comparison is an equality, at $j$ equal to the witness speed. Simultaneous failure conditional on $U$ is $1-o(1)$ even for each fixed eligible target. For the $j=0$ bit alone, negation pairs $T$ with $-T$ under the symmetric masking law, preserves regularity and the addition guards, and switches the ordinate bit of the unique point. The average conditional error for this bit is consequently $1/2+o(1)$. This last symmetry statement is over the masked population, not a fixed-target prediction probability.

## The full Walsh transcript on a fixed target

There is a separate bound for the absolute coefficient required by the [weak parity decoder](../sampled_slots/address_parity.md). Put $c_i(P)=r_i\cdot a_i(P)\bmod2$ at every slot. Allow a bounded score whose dependence on the fresh challenge is **solely through the predictor's equality answers for $j=0,\ldots,m$**. The postprocessing may depend on the fixed seed and target. Then for every fixed eligible target the conditional mean absolute witness coefficient is at most $2^{-b+o(b)}$.

To prove it, freeze all slots except slot $i$. If its complementary sumset has at most $1+2mw$ elements, the chance of a singleton hit is at most $(1+2mw)\rho$ by the forced-point argument. Union over all $m$ slots. Outside this small event and irregular roots, the capped calculation certifies a pole on **every** challenge with some active speed, because the complementary support for any chosen active slot is larger than the cap. The positive shifts always return unequal. Thus the transcript is constant except on challenges for which every local speed vanishes; on those the unshifted comparison returns equal.

Let $H$ be the direct sum of the binary spans of the legal addresses in each slot. The exceptional challenges form exactly $H^\perp$. We have $|H|\ge D$ by injectivity of the slot addresses, and $a_*\in H$. Any score with the two possible transcript values differs from a constant by a scalar of magnitude at most two times $1_{H^\perp}$. For nonzero $a_*$ this gives $|\widehat f(a_*)|\le2/|H|\le2/D$. A zero witness address on $U$ has probability at most $\rho$ before conditioning: freeze the tails of that one address and the moving slot must contain its forced point. For any $d\ge1$, $\Pr(U,D\le d)\le d\rho$ as well, since a hit with $D\le d$ has at most $d$ tail tuples when the nonempty moving slot is dropped. Summing this bound on dyadic ranges gives $\mathbb E[1_U/D]\le2\rho(1+\log_2 D_{\max})$, with the summand defined as zero when there is no hit. Divide all exceptional probabilities and bounds by $q\ge2^{-O(m+\log m)}$ to obtain the claimed fixed-target rate. A score using additional challenge-dependent arithmetic is outside this argument.

The fast average-case comparison therefore fails to improve the end-to-end exponent. The checker tests all endpoint shifts against both literal products at 96 regular targets on one curve, and measures the entire Walsh transcript on a second curve. [Checker](../../checks/check_implicit_chord_projection.py) · [Output](../../evidence/check_implicit_chord_projection.txt)
