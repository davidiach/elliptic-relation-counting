# How much can sampled partial sums be compressed?

The sampled slots have about $b$ bits of choices in total. This makes the *shape* of the arithmetic check important. An exact chain that passes a reusable running-sum state still needs $(1/2-o(1))mb$ bits on roots with appreciable target coverage. A balanced tree permits $O(b\log m)$ bits of state labels. However, building explicit tables to decode blocks introduces a separate exponential cost. The statements below distinguish a bound on the information in the labels from the cost of evaluating the local predicates that use them.

## Model

Let $G=E(\mathbb F_p)$, $N=|G|$, $M=2^b$, $w=2^\ell$, and $\rho=w/M$. Slot 1 uses a uniform odd affine permutation modulo $M$ and keeps addresses below $w$. Each other slot fixes a random assignment to the complement of $\ell$ public bits in the canonical abscissa, giving $s$ or fewer runs of free bits. All slots are independent and retain every affine point above an accepted abscissa. Write $F_i$ for the resulting sets of signed points. The statements use $|F_i|\le2w$, $1/3\le(N-1)/M$, and $N/M\le2$, which hold eventually in the intended regime. At the paid widths and growing polylogarithmic arities of interest, $\ell=b/m+O(1)$, $s=O(1+\sqrt b/\log b)$, and $s+m=o(\ell)$.

For a nonempty set $I$ of slots write $S_I=\sum_{i\in I}F_i$ for its set of distinct group sums, $H_I(y)$ for its number of signed tuples summing to $y$, and $T_I=\prod_{i\in I}|F_i|$. These are unrestricted block sums before testing a prescribed target or the addition guards.

## A coverage bound independent of the mask law

For any fixed root seed and any subset of $r$ slots,

$$
 |S_I|\ \ge\ |S_{[m]}|/(2w)^{m-r}.
$$

Indeed every full sum uses one block sum and one choice per remaining slot. If $|S_{[m]}|\ge 2^{b-\gamma}$, this gives $\log_2|S_I|\ge rb/m-\gamma-O(m)$ simultaneously for all subsets.

This has a consequence for a target-independent sequential representation that must combine a prefix with any legal suffix and return its sum, or answer exact equality for every target and supplied tuple. A prefix message must distinguish all attainable prefix sums. If two unequal sums had the same message, appending the same fixed suffix would give different totals. Querying one of those totals distinguishes them. Thus a fixed-width message at prefix length $j$ needs at least $\log_2|S_I|$ bits. If $\gamma=o(b)$ and $m\to\infty$, summing over $j=1,\ldots,m-1$ gives at least $(1/2-o(1))mb$ transmitted bits. This holds in any ordering, even an ordering selected from the seed. This statement counts the successive state interfaces. A single target-specific joint test can fall outside this model.

The same argument applies to guarded paths if the covered targets mean *guardedly* reachable targets. For a prefix occurring in a guarded full path, append that path's suffix and ask for its actual target. A prefix of different sum must be rejected for that query, whether its own concatenation satisfies the guards or fails them. Thus distinct extendible prefix sums need distinct messages. The guarded target count is at most their number times $(2w)^{m-r}$, giving the same inequality and chain bound under the guarded coverage premise.

Large-coverage roots matter even when many candidates are skipped. If a target is uniform on the nonidentity elements of a subgroup of order $n$, any root with $|S_{[m]}|<n/(4L)$ hits it with probability at most $n/(4L(n-1))$. Among at most $L$ attempts with fresh targets independent of their root conditional on the preceding history, the probability of any success on such a root is at most $n/(4(n-1))$. In particular on a group with $\log_2 n=b-O(1)$, a method succeeding with probability $1/2$ in $L=2^{o(b)}$ attempts must sometimes use roots with $\log_2 |S_{[m]}|\ge b-o(b)$.

## A simultaneous support estimate on the mixed sampling law

With probability at least $1/(2\cdot12^m)$ over one fresh root, every nonempty subset $I$ satisfies

$$
 |S_I|\ \ge\ 2^{-20m}\min(w^{|I|},M2^{-s}).
\tag{1}
$$

This deliberately modest probability establishes a common event for *all* choices of subsets, including trees or partitions selected after seeing the seed. It also gives expected explicit table-size lower bounds. It is not a high-probability assertion.

For $D\in G$, let $A_i(D)=\sum_P\Pr(P,P+D\in F_i)$. In a fixed-bit slot all retained abscissas are translates over the integers of a subset of one fixed free-bit set. A run of $a$ bits allows $2^{a+1}-1$ signed differences, so the abscissa difference set has at most $2^s w$ elements. For $D\ne O$, $x(P+D)-x(P)$ has distinct poles of order two at $O,-D$. Its degree is four, even when $D$ is nonzero two-torsion. Thus at most $4\cdot 2^s w$ point pairs can be retained at displacement $D$, each with probability at most $\rho$. In the affine slot, distinct abscissas have joint retention probability at most $\rho^2$: their hashed difference is uniform among residues of its fixed 2-adic valuation, which gives probability either zero or $\rho$ conditional on one accepted address. Equal abscissas with $D\ne O$ require $2P=-D$, with at most four solutions. Consequently all slots satisfy

$$
 \max_{D\ne O}A_i(D)\le 6\cdot 2^s w^2/M,\quad
 A_i(O)=(N-1)\rho,\quad \sum_D A_i(D)\le(2w)^2 .
$$

Let $C_I=\sum_y H_I(y)(H_I(y)-1)$. Two distinct ordered tuples with the same sum have at least one nonzero slot displacement, and their displacements sum to zero. Choose one nonzero slot, sum the other displacements freely and bound the forced displacement pointwise. Independence gives, with $r=|I|$,

$$
 \mathbb E C_I
 \le r(6\cdot2^s w^2/M)(2w)^{2r-2}
 \le 2r\,2^s(2w)^{2r}/M .
$$

Since $\mathbb E|F_i|\ge w/3$ and $|F_i|\le2w$, every slot has probability at least $1/12$ of containing $w/6$ or more signed points. By independence all $m$ do so with probability at least $p_0=12^{-m}$. By Markov's inequality and a union bound, outside probability $p_0/2$ every $C_I$ is bounded by $2^{m+1}/p_0$ times its displayed upper bound. Both properties hold together with probability at least $p_0/2$. On that event $T_I\ge(w/6)^r$. Cauchy–Schwarz gives

$$
 |S_I|\ \ge\ \frac{T_I^2}{T_I+C_I}
 \ \ge\ \tfrac12\min(T_I,T_I^2/C_I).
$$

Substitution gives (1), with the second term interpreted as unbounded when $C_I=0$.

For a proper block at the specified widths, $w^{m-1}\le M2^{-s}$ eventually. On the common event, its exact sum-state needs $r\ell-O(m)$ bits. For any binary tree, summing over the edges from proper subtrees thus gives a lower bound $\ell\sum_i d_i-O(m^2)$, where $d_i$ is the depth of leaf $i$. Each leaf contributes once for each edge between it and the root. Transmitting the signed leaf addresses in concatenated form along these edges uses at most $(\ell+1)\sum_i d_i$ bits, since $|F_i|\le2w$. In particular the balanced-tree label budget is $\Theta(b\log m)$ on this event. The root would receive the entire tuple if the addresses were simply concatenated. Evaluating the sum there does not meet the inherited low-degree local-predicate requirement. If a uniform admissible compiler achieved $O(b\log m)$ total Boolean-equivalent edges, including the join checks, with the inherited entry costs, then $m\asymp\sqrt{\log b}/\log\log b$ would give log time $O(b\log\log b/\sqrt{\log b})$ by the usual balance. The count of state bits alone gives no such compiler.

## What materializing blocks costs

Take a partition into at most $h$ nonempty blocks whose complete distinct sumsets are materialized. The largest block has at least $m/h$ slots. On (1) the largest table has at least $2^{-20m}\min(w^{m/h},M2^{-s})$ entries. The log of its unconditional expected size is at least $(1-o(1))b/h$ in the regime above if full materialization is performed on every root. The estimate holds for a partition picked after seeing the seed. A policy that discards roots before construction needs separate accounting.

There is a seed-law-free bound with early discards and reuse. Suppose a successful evaluation uses a partition into at most $h$ *complete* block tables, each with at most $K$ distinct entries. Then $|S_{[m]}|\le K^h$, since a full sum is obtained by choosing one value per block. Even a target-dependent partition obeys this constraint: all partitions of the same full slots have the same final sumset. Under a uniform nonidentity subgroup target independent of the root, the chance of success with such tables is at most $K^h/(n-1)$. This argument does not require any lower bound on the size of a typical slot.

In $L$ unit-or-higher-cost bit operations, at most $L$ fresh targets can be tried and each materialized table has at most $L$ entries if its construction is charged. Conditional on the history preceding each freshly randomized target, success within the budget requires a hit on a root with $|S_{[m]}|\le L^h$. A union bound gives

$$
 \Pr(\text{success by }L)\ \le\ L^{h+1}/(n-1).
\tag{2}
$$

Thus the expected cost, or the time budget for success probability $1/2$, is $\Omega((n-1)^{1/(h+1)})$. This charges work to build a saved table once. It permits repeated fresh queries using it. It assumes the root is fixed before the independent fresh target (as in the masked-query law), all preprocessing is charged, and a successful query actually uses the block tables. Lists populated only for a target-dependent restriction of the slot grid, implicit representations, and uncharged prior tables require a different statement.

For a large prime-order subgroup with $\log_2n=b-O(1)$, (2) rules out a time $2^{o(Z_b)}$ for this route unless $h/(\log_2 b)^{1/4}\to\infty$, where $Z_b=b/(\log_2 b)^{1/4}$. If the remaining $h$-way group arithmetic is then priced by the existing full-word packet estimate with order $hb/\sqrt{\log_2 b}$, that estimate already exceeds $Z_b$ in order. Thus explicit block tables followed by that estimate do not certify a better exponent. The computation-time lower bound in (2) applies to the explicit tables. The packet expression is a cost *estimate* for a particular remaining arithmetic construction, not a lower bound on all ways to combine the blocks.

## Reproduction

The finite check sums the actual mixed seed law in three toy cases. It verifies the displacement bounds, each block's exact second moment by an independent convolution of displacement masses, the support and collision inequalities, and the cheap-target bound for every admissible two-block partition. Run from the repository root.

~~~text
python -B checks/check_block_support.py
~~~

On the order-19 curve over $\mathbb F_{13}$, the three-slot and four-slot cases exhaust 2,048 and 65,536 weighted roots and check all 7 and 15 nonempty blocks. With four entries allowed per table, the respective numbers of eligible nonidentity root-target pairs are 1,192 and 103,480. On the noncyclic order-20 curve $y^2=x^3-x$ over $\mathbb F_{19}$, the three-slot case exhausts another 32,768 roots and gives 25,040 such pairs. Altogether the checks cover 100,352 weighted roots and 29 exact block collision identities. [Reproduction script](../../checks/check_block_support.py)
