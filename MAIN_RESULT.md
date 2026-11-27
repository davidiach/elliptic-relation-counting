# A lookup-free polynomial certificate for fixed-shift relations

Research branch, 27 November 2026. Logarithms in resource bounds are binary. The deterministic result is a symbolic rank certificate and an algorithm to evaluate it at a specified point. Random selection gives a Las Vegas theorem with a linear expected random-bit budget.

## Theorem

Let $p>3$ be prime. Let a nonsingular elliptic curve over $\mathbb F_p$ have prime rational group order $r$, and let $P$ generate that group. Put $b=\lceil\log p\rceil$ and $L=\log b$. Assume the uniform **closed** Boolean contraction theorem stated as Liu–Zhang–Ren (LZR) Theorem 2.2. There is a uniform deterministic construction of an implicit matrix polynomial $B$ and vector polynomial $h$ over $\mathbb F_r$ with these properties.

1. The determinant of $B$ is a nonzero polynomial. At every specialization $h=B\lambda$, where $\lambda$ is the true table of logarithms, including the target. Any nonsingular specialization determines this table uniquely.
2. Evaluating at a specified point and solving costs $2^{O(b/L^{1/4})}$ bit time and $2^{O(b/L^{1/4})}$ bits of space, including all setup. The network uses ordinary group arithmetic, at most two point pins per count, and public weights on small address pieces. There are a bounded number of specialization coordinates with $O(b)$ bits each.
3. A uniform fixed-length trial uses $O(b)$ independent fair bits, runs within those deterministic bounds, and fails with probability $2^{-\Omega(b)}$ for sufficiently large $b$. Independent retries give zero error, $O(b)$ expected random bits and expected time on the same scale. A constant list of trials works for every eligible input of a given bit length with $O(b)$ bits of **nonuniform advice**.

The theorem uses the whole prime-order group. No finite-packet A.1 assumption, cofactor extension or numerical leading constant is needed.

## Construction and rank

Use the deterministic seed, indexed blocks and canonical decompositions of LZR Theorem 3.1. Retain duplicate points across blocks as distinct columns. Set $N=N_0+2$, including distinct columns for $P,Q$. Form the explicit $N+1$ consistency and normalization rows from (4.1). Append conceptually every ordered fixed-shift pair from (4.2). Write $A$ for the resulting matrix of row vectors $a_i^t$ and $c_i$ for their right sides. Occurrences on opposite sides contribute opposite signs. The valid system contains $\lambda$. Two tuples for $g$ share a successor at $g+P$, so any solution assigns both the same label sum. Advancing by $P$ increments this function by one. The row at $P$ fixes the offset. The consistency rows then fix every column. Hence $\ker A=0$.

Concatenating the $4m$ fixed-length slot addresses gives a distinct $\beta$-bit integer $v$ for each ordered pair. Set $C=N+1$. Give the explicit rows keys $0,\ldots,C-1$, and the pair rows keys $C+v$. Define

$$
 B(Z)=\sum_i Z^{e_i}a_i a_i^t,\qquad h(Z)=\sum_i Z^{e_i}c_i a_i.
$$

Cauchy–Binet gives $\det B(Z)=\sum_{|I|=N}\det(A_I)^2 Z^{\sum_{i\in I}e_i}$. Scan rows by increasing key, retaining independent ones. After every prefix the number retained equals the prefix rank. Thus the ordered keys of the greedy basis are coordinatewise at most the ordered keys of any other basis, with a strict inequality in at least one coordinate. Its key sum is uniquely minimal and its coefficient is a nonzero square. Other degrees may cancel. This applies a classical linear matroid argument to a relation family too large to enumerate.

## Few-variable evaluation

Put $d=\max(1,\lfloor b/4\rfloor)$, $D=2^d$, and let $s$ be the least positive integer such that $D^s\ge 2^\beta$ and $D^s>C$. For an integer $v<D^s$ with base-$D$ digits $v_j$, set $M_v(X)=\prod_{j=0}^{s-1}X_j^{v_j}$. Replace $Z^e$ on an explicit row by $M_e$ and on a pair by **$M_C M_v$**. The latter choice avoids carries in the private address. Substitution $X_j=Z^{D^j}$ recovers the univariate matrix, proving nonvanishing. The determinant's total degree is at most $2Ns(D-1)$ and $h=B\lambda$ holds termwise. From Stage I, $N=2^{O(b/m+\log b)}$, $\beta=O(b)$ and $s=O(1)$ at the balanced arity.

Sample each $X_j$ as an integer in $[0,D^2)$ using exactly $2d$ fair bits. Hasse's bound ensures $D^2<r$ for large $b$. Schwartz–Zippel bounds singularity by $2Ns(D-1)/D^2=2^{-b/4+O(b/m+\log b)}$. Elimination detects singularity. A nonsingular solve recovers the exact table, which can also be checked by scalar multiplication. The inequalities can be tested on the actual parameters. A fixed finite initial range uses exhaustive group arithmetic.

For the advice statement use canonical field coordinates for the boundedly many curve coefficients and input points. There are at most $2^{cb+O(1)}$ inputs at each length. Pad trial tapes to a common $O(b)$ length. A sufficiently large constant number of independent tapes makes the probability that any input fails on every tape less than one. This proves existence only. Searching for such tapes across all inputs is outside the time bound.

## Bit cost of a single evaluation

Add the explicit rows directly. On the other rows $a_j$ is a signed sum of at most $4m$ indicators of a slot equaling column $j$. Expand $a_j a_k$ into at most $16m^2$ signed pairs of point pins. Incompatible pins on one slot give zero, and compatible identical pins give a single pin. For the right side use single pins since every shift row has $c_i=1$. Apply signs outside the counts. Each summand is the weighted count of ordinary shift pairs under at most two pins. Parsimony follows from the full group law, including repeats, inverses and identities, and the per-block address injection.

Within each base-$D$ digit cut the address into contiguous pieces of length at most $\lceil\sqrt b\rceil$, also splitting at slot boundaries. A piece of value $v'$ starting at bit $j$ of digit $i$ contributes the canonical field representative of $X_i^{2^j v'}$, using binary powering and $0^0=1$. There are $q=O(\sqrt b+m+1)$ pieces. Attach a tensor per piece to the dedicated wires and subdivide as usual. This graph has $E=O(mb)$ Boolean edges, polynomial description length, and maximum degree $O(b^{3/4}+\sqrt b)$. Entry time is polynomial in $b,m$.

The contractor works modulo powers of two. Choose $S=\Theta(mb)$ large enough to dominate $E+b$. Split each nonnegative field representative into digits of $t=\max(1,\lfloor S/(q+1)\rfloor)$ bits. There are $J=\lceil\lceil\log r\rceil/t\rceil$ positions. Enumerate the $J^q$ choices of one digit per weight. Every such scalar contraction is an integer less than $2^{E+qt+1}$. Apply Theorem 2.2 with scalar ring, precision $E+qt+1=O(S)$ and size parameter a sufficiently large constant times $S$. Take the canonical integer residue, multiply by $2^{t\sum j_v}\bmod r$, and sum modulo $r$. Apply the common public factor $M_C$ once. This proves exact recovery without a ring map from dyadic residues to $\mathbb F_r$. The logarithm of the digit enumeration overhead is $O(q\log(b+2))$.

For growing $m\le L^{O(1)}$, in the regime $b/m\gg\log b$, log bit time for evaluations, Stage I, explicit rows, elimination and verification is

$$
 O\left(b/m+mb/\sqrt{\log(mb)}+(\sqrt b+m)\log(b+2)\right).
$$

Each count streams in polynomial space; matrices, tables and sequential tape scans cost polynomial resources in $N,b$. Balance the first two terms at $m=\Theta(L^{1/4})$. No counterexample search or candidate iteration is needed.

## Uniform worst-case question

The follow-up in [Deterministic points](DETERMINISTIC_POINTS.md) gives an explicit target-budget list for the original $B$. Its stack contains all explicit rows and detects every direction with a bounded-degree sparse polynomial consequence. A separate density theorem locates many violations for each residual direction. Uniform full rank of this list remains unproved.

[Exact neighboring fibers](NEIGHBORING_FIBERS.md) settles the full-grid fallback for every eligible curve on that branch. It also gives the exact all-displacement rank and explains why that identity does not settle the prescribed shift on the short list.

The construction and each specified evaluation are deterministic. Enumerating the full interpolation grid costs $2^{O(b)}$ by this proof. An explicit hitting list constructible in $2^{O(b/L^{1/4})}$ time, of at most that cardinality and guaranteed to contain a nonsingular specialization for each eligible input, would give the target uniform worst-case bound. A list may depend on the input if its construction is charged. The original LZR separator remains a worst-case fallback.
