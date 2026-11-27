# Algebraic and Residual Barriers for the Uniform Deterministic Fixed-Shift Endgame

Research follow-up, 27 November 2026. Parent: `main` @ `5435264` (`research/fixed_shift_gram/manuscript.typ`, `DETERMINISTIC_POINTS.md`, `NEIGHBORING_FIBERS.md`, `PARTIAL_EXPANSION.md`, `COMPLEXITY_AUDIT.md`, `COMPATIBILITY_MAP.md`). These are AI-assisted proof candidates and exact finite checks, not an independent referee report or a completed uniform deterministic exponent improvement. The uniform deterministic quarter-power ECDLP time bound remains **OPEN**.

## Summary

The 27 November 2026 fixed-shift manuscript leaves two routes toward a uniform deterministic time bound $2^{O(b/(\log b)^c)}$ ($c > 0$) under Liu–Zhang–Ren (LZR, ePrint 2026/3020):

- **Lever A (Fixed-generator noncancellation):** Either establish full rank for $O(1)$-slot-pinned unweighted shift Gram matrices $K_{A_1,\dots,A_p}$ summed over all group values $g \in G$ (the $2^{O(b/(\log b)^{1/4})}$-budget relaxation directly below `NEIGHBORING_FIBERS.md` Proposition 4), or construct a deterministic hitting stack for decoupled polynomial evaluations $B_{\mathrm{shift}}(x, y)$.
- **Lever B (Residual lookup contraction after Theorem 2):** Break the $2^{O(b/\sqrt{\log\log b})}$ bottleneck left after Theorem 2's one-step $R = \Theta(m)$ partial expansion.

This note establishes:

1. **Pinning barriers (Theorems A1 & A2):** Unconditional $O(1)$-slot external pinning $K_{A_1,\dots,A_p}$ fails identically modulo $r$ via Lucas's Theorem at $m = r - 1$ for unweighted $F_i = \{1, 2\}$ ($p \le r - 7$), **and** fails identically modulo $r$ at strict $4m < r$ ($m \ge p + 3$, $r = m d + 1$, $d \ge 5$, $|D_0| = r = 100\%$ coverage) for positive unit binomial weights $F_i(U) = U(1+U)^d$. No proof using only $4m < r$, $|D_0| = r$, and unit group-algebra structure in $\mathbb{F}_r[U]/((U-1)^r)$ can settle $K_{A,B}$.
2. **Sharp positive & empirical boundary for unweighted 0/1 pinning (Propositions A3 & A4):** For $m = 2$, $2m$-slot (`cpb = 1`) shifts, $\{K_{A,B} : A, B \in F_1\}$ has the exact kernel of all individual shift rows for every prime $r$ and all proper subsets $F_1, F_2 \subsetneq \mathbb{Z}/r\mathbb{Z}$. For general unweighted 0/1 blocks with $2\cdot\mathrm{cpb}\cdot m < r$, $K_1$ and $K_2$ had **0 rank drops** across exhaustive sweeps ($r = 7, m = 3$: 969 triples; $r = 11, m = 2$: 185 identical pairs), Paley quadratic-residue blocks, and elliptic curve controls ($r = 19, 37, 257$).
3. **Residual estimates (B1–B3).** The core-only cost expression balances at $b/L^{1/6}$. With the lookups attached, their degree limits the saving certified by the published residual estimate. B2 bounds the expected packing potential for uniformly labeled clusters and accounts for both edge orientations. B3 proves a boundary of at least $(4/9-o(1))\beta$ bits for binary contraction trees on a balanced complete bipartite graph.

---

## 1. Lever A: Pinning Barriers, Forest Theorem, and Decoupled Evaluations

For $m$ blocks $F_1, \dots, F_m \subset G \setminus \{0\}$ (each appearing $\mathrm{cpb} \in \{1, 2\}$ times per tuple; $\mathrm{cpb} = 2$ gives LZR's $4m$-slot doubled pairs, $\mathrm{cpb} = 1$ gives $2m$-slot short pairs), write $K_{A_1, \dots, A_p}$ for the shift Gram matrix summing $(c_Y - c_X)(c_Y - c_X)^\top$ over all shift pairs $\mathrm{val}(Y) - \mathrm{val}(X) = P$ satisfying $p$ fixed slot pins, summed over all $g = \mathrm{val}(X) \in G$. Notice that $K_A$ ($p = 1$) and $K_{A, B}$ ($p = 2$) require $O(N^3)$ and $O(N^4)$ closed contractions respectively—both within the $2^{O(b/(\log b)^{1/4})}$ bit budget with zero polynomial weights or primitive-root setup.

### Theorem A1 (Unweighted Lucas Collapse at $m = r - 1$)
Let $r \ge 11$ be prime, set $m = r - 1$, and take $m$ identical blocks $F_1 = \dots = F_m = \{1, 2\} \subset \mathbb{Z}/r\mathbb{Z}$ with doubled tuples ($\mathrm{cpb} = 2$, $4m = 4r - 4$ slots). Then:
1. $|D_0| = |\{m, m+1, \dots, 2m\}| = r > 7r/8$, every block has size $2 < r$, and the full shift row space has rank $N_0 - m = r - 1$.
2. For every $0 \le p \le r - 7$ and every choice of $p$ pinned slots across $X$ and $Y$, the unweighted shift Gram matrix $K_{A_1, \dots, A_p}$ over $\mathbb{F}_r$ is identically zero.
3. There exists a valid explicit consistency+normalization matrix $E$ (using genuine $2m$-slot decompositions of every labeled column) with $\mathrm{rank}(E) = N - 2$, giving
   $$\mathrm{rank}\!\left(E \;\cup\; \bigcup_{p=0}^{r-7}\{K_{A_1,\dots,A_p}\}\right) = N - 2 < N.$$

**Proof.** Writing each element as $1 + x_\ell$ ($x_\ell \in \{0, 1\}$), any entry of $K_{A_1,\dots,A_p}$ is a signed sum of valid shift counts with $k = k_X + k_Y \le p + 2 \le r - 5$ total pinned slots, leaving $M = 4r - 4 - k = 3r + (r - 4 - k)$ free binary slots. The shift congruence $\mathrm{val}(Y) - \mathrm{val}(X) \equiv 1 \pmod r$ becomes
$$\sum_{\text{free } Y} y_\ell + \sum_{\text{free } X} (1 - x_\ell) \equiv 1 + (2r - 2 - k_X) + \sum_{\text{pin } X} x_\ell - \sum_{\text{pin } Y} y_\ell \pmod r,$$
so every valid Hamming weight $s$ satisfies $s \bmod r \in [r - 1 - k,\, r - 1]$. Because $k \le r - 5$, the base-$r$ units digit of $M$ is $M \bmod r = r - 4 - k < r - 1 - k \le s \bmod r$. Lucas's Theorem gives $\binom{M}{s} \equiv \binom{3}{\lfloor s/r \rfloor}\binom{r - 4 - k}{s \bmod r} \equiv 0 \pmod r$, killing every entry.

For (3), choose the decomposition for each indexed column as follows. Start all slots at label 1. A column of value $g$ needs $k_g=(g-2m)\bmod r\in\{3,4\}$ twos. Put them in the first $k_g$ slots of blocks 2 and 3 in block order. For each of the exceptional columns $(0,2),(1,2)$ put the first two-valued entry in its own block and the remaining $k_g-1$ in blocks 2 and 3. Take the dedicated columns to have values 1 and 2. In the homogeneous kernel every ordinary column of a given value has the same right side and thus the same coordinate. Normalization at the dedicated 1 sets the common coordinate for 1 to zero. Its consistency row then sets three times the common coordinate for 2 to zero. Since $r\ge11$, both vanish. Each exceptional coordinate occurs once on each side of its own row and nowhere else. These two free coordinates give rank $N-2$. This example allows the tuple choice to depend on the indexed column, including for repeated group values. $\square$

### Theorem A2 (Small-$m$ Dense Binomial-Weighted Collapse at Strict $4m < r$)
Let $p \ge 0$ be an integer, let $m \ge p + 3$, and let $r = m d + 1$ be any prime with $d \ge 5$ (so $4m = \frac{4}{d}(r-1) < r$ strictly, accommodating LZR's asymptotic regime $m = \Theta((\log b)^{1/4}) \ll r$). Take $m$ identical blocks $F_1 = \dots = F_m = \{1, 2, \dots, d + 1\}$ with doubled tuples ($\mathrm{cpb} = 2$) and positive binomial weights $w(j + 1) = \binom{d}{j} \in \mathbb{F}_r^\times$, so $F_i(U) = U(1 + U)^d$, $F_i(1) = 2^d \not\equiv 0 \pmod r$, and $|D_0| = m d + 1 = r$. Then:
1. For every $0 \le p \le m - 3$ and every choice of $p$ pinned slots, the weighted shift Gram matrix $K_{A_1, \dots, A_p}$ over $\mathbb{F}_r$ is identically zero.
2. There exists a valid explicit matrix $E$ with $\mathrm{rank}(E) < N$ such that
   $$\mathrm{rank}\!\left(E \;\cup\; \bigcup_{p=0}^{m-3}\{K_{A_1,\dots,A_p}\}\right) = \mathrm{rank}(E) < N.$$

**Proof.** Each free $Y$-slot contributes $(1 + U)^d$ and each free $X$-slot contributes $U^{-d}(1 + U)^d$. With $k = k_X + k_Y \le p + 2 \le m - 1$ pinned slots, $M = (4m - k)d = 4(r - 1) - k d = 3r + (r - 4 - k d)$ has $M \bmod r = r - 4 - k d \ge d - 3 \ge 2$. Meanwhile, $2m d = 2(r - 1) \equiv -2 \pmod r$ forces $s \bmod r \in [r - 1 - k d,\, r - 1] > r - 4 - k d = M \bmod r$, so $\binom{M}{s} \equiv 0 \pmod r$ by Lucas's Theorem.

For (2), take $P=1$, $Q=2$, and let $u$ be the indicator of the label $d+1$ on the indexed columns. The normalization row kills $u$. Start each of the $2m$ slots at 1. For a target $g=d+1$, set the first slot to $d+1$ and reserve it. For every other target reserve no slots. Write $h=1$ for the first case and $h=0$ for the second. Increase the remaining $2m-h$ slots by amounts in $[0,d-1]$ with total
$$ (g-2m-hd)\bmod r. $$
Every total from 0 to $(2m-h)(d-1)$ is attainable by filling slots in order. Even with a reserved slot, $(2m-1)(d-1)\ge md=r-1$ for $m\ge3,d\ge5$. Thus the tuple sums to $g$ and contains exactly $h$ entries with label $d+1$. Its score equals $u_g$, so every consistency row kills this nonzero vector. The construction uses a single tuple for each group value, including repeated labels. $\square$

*Scope boundary:* Theorem A2 rules out any proof of full rank for $\{K_{A_1,\dots,A_p}\}$ ($p \le m - 3$) that uses only $4m < r$, $|D_0| = r$, $n_i < r$, and unit group-algebra structure in $\mathbb{F}_r[U]/((U-1)^r)$. The family has nonconstant weights. The normalized one-per-block sum law is a translate of $\operatorname{Bin}(r-1,1/2)$ modulo $r$. Its $r$ possible integer totals occupy distinct residues, so its $\chi^2$ distance from uniform is $r\binom{2(r-1)}{r-1}/4^{r-1}-1=\Theta(\sqrt r)$.

### Proposition A3 (Exact Forest Theorem for $m = 2$, $\mathrm{cpb} = 1$)
For every prime $r$ and every pair of nonempty proper subsets $F_1, F_2 \subsetneq \mathbb{Z}/r\mathbb{Z}$, the $2m$-slot two-pin stack $\{K_{A,B} : A, B \in F_1\}$ satisfies
$$\ker\big(\{K_{A,B} : A, B \in F_1\}\big) = \ker(\text{all individual } 2m\text{-slot shift rows}).$$

**Proof.** Fix $A, B \in F_1$ and put $t = P + A - B$. If $t = 0$, then $Y_2 = X_2$ for every shift pair, so $K_{A,B} = |F_2|(e_B - e_A)(e_B - e_A)^\top$; since $0 < |F_2| < r$, $K_{A,B} u = 0$ forces $(e_B - e_A)^\top u = 0$. If $t \ne 0$, the directed shift graph $S_t = \{(x, y) \in F_2^2 : y - x = t\}$ sits inside a single $r$-cycle ($r$ prime) on $|F_2| < r$ vertices, hence is a forest of directed paths. Its edge-incidence vectors $\{b_e = e_y - e_x : e \in S_t\}$ (supported on $F_2$, disjoint from $F_1$) are linearly independent over $\mathbb{F}_r$. Restricting $K_{A,B} u = \sum_{e \in S_t} (a_e^\top u) a_e = 0$ to $F_2$-coordinates gives $\sum_{e \in S_t} (a_e^\top u) b_e = 0$, hence $a_e^\top u = 0$ for every individual shift row $e \in S_t$. $\square$

### Proposition A4 (Exact Reduction of $K_{A,B} u = 0$)
For fixed blocks $(i_1, i_2)$ and pins $A \in F_{i_1}$ (on $X$) and $B \in F_{i_2}$ (on $Y$), put $t = P + A - B$ and let $n(t), v(t), G'(t)$ be the 0th, 1st, and 2nd moments of $c_{Y'} - c_{X'}$ over unpinned complements with $\mathrm{val}(Y') - \mathrm{val}(X') = t$. Using $\lambda^\top(c_Y - c_X) = 1$,
$$K_{A,B} u = 0 \iff \begin{cases} v(t)^\top u + (u_{i_2, B} - u_{i_1, A})\, n(t) = 0, \\ G'(t) u + (u_{i_2, B} - u_{i_1, A})\, v(t) = 0. \end{cases}$$

**Proof.** Expand $K_{A,B} = n(t)(e_B - e_A)(e_B - e_A)^\top + (e_B - e_A)v(t)^\top + v(t)(e_B - e_A)^\top + G'(t)$, left-multiply $K_{A,B} u = 0$ by $\lambda^\top$ to obtain the scalar equation, and substitute back into $K_{A,B} u = 0$. $\square$

The number of pairs with the same $t=P+A-B$ depends on the differences in the blocks. An interval of $n$ elements has $n-1$ ordered pairs with difference one, even when $n\ll\sqrt r$. In the separate $r=7,m=3$ sweep, the scalar projection $\lambda^\top K_1$ fell below $N_0-m$ on 66 dense triples. The full one-pin stack retained the individual shift rank on all 969 triples.

A weighted example separates one pin from two. Over $\mathbb F_7$, take blocks $(1,2,3),(1,2,4),(1,3,5)$ with respective weights $(3,2,3),(1,5,4),(5,3,5)$ and one occurrence per block on each side. The source-only and target-only one-pin stacks each have rank 5 and annihilate the coordinate vector for label 2 in the first block. The stack with one pin on each side detects this vector and has full shift rank 6.

### Integer minors of finite two-pin stacks

The verifier enumerates 18 unweighted fixtures. Eleven use shift equality over the integers and seven use a shift modulo 7 or 11. In every case the outer products and the stacked Gram matrices have integer entries. Let $d=\sum_i(|F_i|-1)$. Each incidence difference sums to zero on each block. Hence the final row and column of each block in a pinned Gram are integer combinations of the others. Dropping one reference coordinate per block on each axis preserves the GCD of all $d$-by-$d$ minors of the stack.

Integer Euclidean row operations then give a triangular $d$-row matrix above zero rows. The absolute product of its pivots is the GCD. Fifteen fixtures give 1, two give 3, and one gives 27. Thus three of the eighteen integer stacks lose rank upon reduction modulo 3. For a fixture using a modular shift, the shift modulus determines which tuple pairs enter the fixed integer matrix. The prime used to test its rank can be chosen independently.

### Decoupled $B_{\mathrm{shift}}(x, y)$ Evaluations (Regime A-I)
On the order-19 curve over $\mathbb F_{13}$, the verifier compares the weighted shift Gram matrix with direct enumeration of all ordered shift pairs. It also evaluates $Bu$ and $u^\top B u$ from separate source and target fiber sums.

Splitting the single radix digit straddling the $X/Y$ cut yields $s_X + s_Y \le s + 1 = O(1)$ independent variables and drops $M_C$ (applied externally, with $E$ stacked directly). Against fiber-respecting dense $m = 2$ adversaries at $r = 7$, fixed coupled 2-point lists (`A1_coupled_L2`) drop to rank 4/6 on `[[1,2,3,4],[1,2,3,4]]` and 1-sided zero-projection 2-point lists (`A3_zero_L2`) drop to 5/6 on `[[1,3,5,6],[1,3,5,6]]`, while decoupled asymmetric 2-point lists (`A2_asym_L2`) hit full rank 6/6 on both. Constructing a worst-case deterministic hitting list of size $2^{O(b/L^{1/4})}$ remains open.

---

## 2. Residual cost estimates and balanced graph cuts

Anchor all parameters in `PARTIAL_EXPANSION.md` and `COMPLEXITY_AUDIT.md`: $L = \log_2 b$, $S = \Theta(mb)$, $\lambda = \lceil \log_2 S \rceil$, $\beta = 4mk = \Theta(b)$, $\epsilon = 2^{-17}$, $u = \lfloor \epsilon \lambda / \log_2 R \rfloor$, $a = 2u$, $W = 2^a$, $D = R^{2a}$, $H = \lceil S / \lfloor \sqrt{\lambda}\rfloor \rceil$, $\Delta = \lceil S^{9/10}\rceil$, and residual contraction parameter $\rho = \max(1, \lfloor \log_2(\beta / \max(\Delta_{\mathrm{res}}, \lceil \beta^{9/10}\rceil)) \rfloor)$ with log bit cost $O(\beta / \sqrt{\rho} + n\lambda + S/\sqrt{\lambda})$.

### Proposition B1 (Optimizing the published cost expressions)

1. At fixed total fan-out write $x_j=\log_2 R_j>0$ and $x=\sum_j x_j$. Budgeting each level at its frontier mass bound $2^j S$ and substituting the packing scale $u_j\asymp\lambda/x_j$ gives the comparison expression $(S/\sqrt\lambda)\sum_j 2^j\sqrt{x_j}$. For two or more levels, $\sum_j2^j\sqrt{x_j}>\sqrt x$, by squaring the sum of the positive square roots. A single level minimizes this expression. The parent expansion bounds its number of terms from above and supplies no matching count from below.
2. For the core alone, consider the maximum of the three budget terms $b/m$, $mb\sqrt{x/L}$, and $b/\sqrt{x}$, with $x=\log_2 R$. Their product is $b^3/\sqrt L$, so their maximum is at least $b/L^{1/6}$. Taking $m\asymp L^{1/6}$ and $x\asymp L^{1/3}$ attains this order in the expression. The verifier also evaluates the packing, degree and precision estimates at three finite parameter choices.
3. Reattaching the $4m$ lookups gives $\Delta_{\rm res}\ge k=\beta/(4m)$, hence $\rho\le\max(1,\lfloor\log_2(4m)\rfloor)$. In the published regime $m=L^{O(1)}$, the term $\beta/\sqrt\rho$ in the upper-bound formula is therefore $\Omega(b/\sqrt{\log L})$. The deterministic upper bound at the manuscript tuning remains $2^{O(b/\sqrt{\log L})}$. A running-time lower bound would require a separate argument.

### Proposition B2 (Expected edge-class deviations under tied labels)

Write $q=R^a$ and $D=q^2$. Orient the edges and assign uniform labels in $[R]^a$. A cluster of terminal hosts shares either its entire label or its first coordinate. All other coordinates are independent, except where two edges use the same vertex. Assume the two endpoints of each edge have independent labels. For two terminal edges from one cluster, assume the other endpoints are distinct vertices with labels independent of each other and of the cluster.

For ordered edge classes $C_e,C_f$, the sum of the same-class indicator covariances is $\Pr(C_e=C_f)-1/D$. With the cluster in the same position on both edges, a full tie needs the two other labels to agree, an event of probability $1/q$. A one-coordinate tie also needs the untied coordinates of the two hosts to agree, giving probability $R/D$. With the cluster in opposite positions the two other endpoints must each match the opposite host, with conditional probability $1/D$ in either tie model. The three covariance values are therefore $1/q-1/D$, $(R-1)/D$, and zero, respectively.

Let $C$ count the unordered within-cluster pairs with the cluster in the same position and independently labeled other endpoints. If $X_\gamma$ counts the edges of class $\gamma$ among $M$ edges, the expected edge potential $\mathbb E\sum_\gamma(X_\gamma-M/D)^2$ is at least $2C$ times the applicable same-position covariance. Indeed every edge is uniform on the $D$ classes, and every other pair contributes a nonnegative **summed** covariance: agreement in one coordinate imposes at most two equality constraints on independent uniform symbols, so the class collision probability is at least $R^{-2a}=1/D$.

For $4m$ slots partitioned into clusters of $s\mid k$ edges with a common orientation within each cluster and independent other endpoints, $C=4m(k/s)\binom{s}{2}$. With full-slot ties ($s=k\ge2$), $M=O(mb)$ and $\beta=4mk=\Theta(b)$, this contribution divided by $(M/D)^2$ is $\Omega(D^{3/2}/m^3)$ for a full-label tie, or $\Omega(RD/m^3)$ for one coordinate. Both grow for $R\ge4m$, $a=2u\ge2$ and growing $m$. In this model the initial expectation exceeds the scale needed for the packing certificate in Lemma 5.1. The minimum potential over deterministic labelings can still be zero. The verifier enumerates a four-edge example with mean potential 4 and minimum 0.

### Theorem B3 ($4/9$-Cutwidth Barrier on Balanced Residual $K_{R, 4m}$)
In a balanced Theorem 2 residual bipartite graph $K_{R, 4m}$ ($R, 4m \ge 4$) where each stopped child $j \in [R]$ shares $\frac{\beta}{4mR}(1 \pm o(1))$ bits with each lookup $\ell \in [4m]$, every binary contraction tree contains an intermediate tensor with boundary bit-width at least $(\frac{4}{9} - o(1))\beta = \Omega(b)$.

**Proof.** Assign leaf weights $\frac{1}{R}$ (left) and $\frac{1}{4m}$ (right); the root has weight $2$ and every leaf has weight $\le \frac{1}{4} < \frac{2}{3}$. Walking down the heavier child from the root reaches a subtree with left/right fractions $(\alpha, \gamma) \in [0, 1]^2$ satisfying $s = \alpha + \gamma \in [\frac{2}{3}, \frac{4}{3}]$. The cut bit-fraction is
$$\alpha(1 - \gamma) + \gamma(1 - \alpha) = s - 2\alpha\gamma \ge s - \frac{s^2}{2} \ge \min_{s \in [2/3, 4/3]}\left(s - \frac{s^2}{2}\right) = \frac{4}{9}. \qquad \square$$

B3 requires balanced edge multiplicities and measures the width of binary contraction trees. The running time of other evaluation methods remains open.

---

## 3. Reproducibility

The three scripts require only Python 3.10+ standard library (run from the repository root):

```sh
python3 -B research/deterministic_endgame_barriers/check_pinning_barriers.py
python3 -B research/deterministic_endgame_barriers/check_residual_barriers.py
python3 -B research/deterministic_endgame_barriers/check_integer_minor_gcd.py
```

The command outputs are recorded in `research/deterministic_endgame_barriers/verification.txt`. The scripts exercise finite-field examples, integer minors, finite samples of B1’s cost estimates, exact rational covariance and balanced cut fractions.

## 4. Open Frontier

Does the unweighted two-pin stack $\{K_{A,B}\}$ attain full shift rank under $4m<r$ and Stage I's density and flatness bounds? The weighted counterexample in A2 leaves this question open. Proposition A4 gives two equations for investigating the kernel.
