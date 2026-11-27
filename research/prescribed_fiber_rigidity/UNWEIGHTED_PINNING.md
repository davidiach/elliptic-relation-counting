# Reflection orbits inside a prescribed sum

**Positive result.** For any proper \(S\subset C_r\), \(r\) odd and \(3|S|>2r+1\), form the ordered three-sum fiber at any \(t\). Stack its cell Grams, fixing the same third entry \(z\) on both sides of each contrast. Their common kernel is exactly the affine functions on \(S\). This proof and its counterexample use ordinary unit weights.

In cell \(z\), the remaining entries occur in orbits of \(x\mapsto t-z-x\). Each unordered pair gives an incidence vector of disjoint support, and has multiplicity one or two. A fixed point gives \(2e_x\). If \(V\) has these columns, \(w\) is the vector of orbit sizes, and \(n\) is the number of ordered pairs, the Gram of all ordered differences is
\[
 2V\big(n\,\mathrm{diag}(w)-ww^\top\big)V^\top .
\]
Here \(2|S|-r\le n\le |S|<r\). Left injectivity and the unit weights force all orbit scores equal inside each cell. The prescribed-fiber theorem then forces a common affine score. This is a family of pinning results for all eligible one-block supports and every prescribed target. It uses a different relation family from PR #7's two-block fixed shift.

**General cell criterion.** In any number of slots, group equal incidences in one cell, writing \(p_j\) for its distinct incidence vectors, \(w_j\) for counts and \(n=\sum w_j\). If \(n\) and all \(w_j\) are units and the \(p_j\) are affinely independent modulo \(r\), all-pairs contrast Gram and raw contrasts have the same kernel. Indeed the coefficients \(w_j(n p_j^\top u-\sum_h w_hp_h^\top u)\) sum to zero. A Gram kernel vector gives an affine dependence with these coefficients, hence each coefficient vanishes. This suffices even if the columns are linearly dependent. It applies to common refinements with affinely independent orbit incidences, wherever their cost can be charged. The number of such orbits can exceed the affine dimension for growing arity.

For the fixed-shift Gram in PR #7 there is no subtraction of a cell reference. In that setting a matrix of linearly independent distinct shift rows and unit multiplicities suffices for equality with raw rank. PR #7 A3 supplies this via projection to the incidence columns of a forest. In PR #4's cross-mean system the defect is instead \(\dim((\operatorname{col}A+\langle1\rangle)\cap\operatorname{col}(H)^\perp)\) when a pin mass is a unit. These are different measurements.

**Failure without the pin.** On \(S=C_r\setminus\{0\}\), \(t=0\), the unpinned Gram has rank \((r-3)/2\) and the pinned stack has rank \(r-3\) for all primes \(r\ge5\). Single positional masses are \(r-2\) and the short law has chi-square \(1/(r-1)\). Neither the example nor this one-block positive theorem settles the several-block Stage I two-pin Gram. [Calculation](COUNTEREXAMPLES.md)
