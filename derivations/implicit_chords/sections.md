# Exact chord equality with bounded storage

Fix a regular grid of $m$ sampled slots, the endpoint $-T$, the original lifts, and any field-valued tangent speed at each row choice. Let $h$ be the number of ordinary hits. The [balanced evaluator](../product_differences/chords.md) returns a normalized field digit for the entire grid. A single smaller section can use that same evaluator. There is an exact rule for combining its answer with those of other sections, even when some sections have hits.

For a section store $(\min(h,2),d)$, where $d$ is its logarithmic derivative sum for no hits, its unique hit's sum of speeds for one hit, and zero for two or more hits. Empty Cartesian products of choices arising from a zero-size slot are treated as grids with no tuples and digit zero. Take two disjoint sections:

| Combined hit count | Combined digit |
|---|---|
| Zero | Add the two digits in $\mathbb F_p$. |
| One | Take the digit of the section with the hit. |
| Two or more | Set the digit to zero. |

The product rule proves the first row because all baseline factors are units. With one vanishing factor, all changes to the unit factors drop out of the product difference modulo $p^2$; the equality bit depends only on the summed speeds of the unique tuple. With two vanishing factors both full products vanish. Thus the summaries combine associatively, and the products are equal precisely when the final digit is zero. Comparing section equality bits or adding unnormalized digits across all three cases would lose information.

Fix $1\le k\le\lfloor(m-2)/2\rfloor$. Split the slots into two groups of $k$ free slots and $m-2k$ pinned slots. Enumerate assignments of the pinned slots depth first. Form a section by retaining each pin as a singleton row, with its own point and speed. Distribute the rows so each side of the chord calculation has $k$ free slots; the singleton endpoint may go on either side. Row permutations change the baseline and moved determinants by the same sign. Interleave original left and right slots when choosing the free groups. Only the two smaller lists for the current section and its field summary need to be kept. No original half list or lifted product is constructed.

Each slot has at most $2w$ points. There are at most $(2w)^{m-2k}$ sections, and each needs $(2w)^k\operatorname{poly}(b,m)$ time and space, including its row-pair terms. Slot expansion takes $O(mw\operatorname{poly}(b,m))$ expected bit time. Overall bounds, up to polynomial factors in $b,m$, are

$$
\mathrm{time}\ \le (2w)^{m-k}, \qquad \mathrm{space}\ \le (2w)^k .
$$

For even $m$ choosing $k=m/2-1$ removes a whole slot-width factor from the maximum table length at the expense of that factor in the worst-case arithmetic bound. At $\log_2 w=b/m+O(1)$ the log time is still $b/2+O(b/m+m+\log b)$, for either parity of $m$ with $k=\lfloor(m-2)/2\rfloor$. At $k=\lfloor m/3\rfloor$ (admissible for $m\ge5$), log time is $2b/3+O(b/m+m+\log b)$ and log space is $b/3+O(b/m+m+\log b)$. The [earlier three-block resultant](../implicit_resultant/resultant.tex) has the same latter leading tradeoff for a different scalar. Here all three fiber cases of the literal lifted chord comparison are covered, including the unrestricted speed of the endpoint.

These are worst-case upper bounds for a storage tradeoff. The section algorithm visits the surviving assignments and repeats arithmetic across pins. Even near the balanced split its largest-side work on typical successful seeds is $2^{b/2-o(b)}$: all but a vanishing fraction of those seeds have each slot larger than $\nu/m^3$, with $\log_2\nu=b/m+O(1)$. For $k=m/2-O(1)$, even one whole section thus has at least $(\nu/m^3)^k$ tuples in a stored side. This does not lower-bound methods working jointly on the seeds. [Slot-size law](../direct_parity/scalar_join.md)

The checker compares field digits with the unsliced calculation and final equality with literal determinant products. The fixtures exercise disjoint four- and six-slot bit masks, nonzero endpoint speeds, and sections with different hit counts within the same calculation. [Checker](../../checks/check_chord_sections.py) · [Transcript](../../evidence/check_chord_sections.txt)
