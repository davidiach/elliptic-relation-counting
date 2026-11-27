# Algebraic and Residual Barriers for the Uniform Deterministic Fixed-Shift Endgame

Research follow-up, 27 November 2026. Parent: `main` @ `5435264` (`research/fixed_shift_gram/manuscript.typ`, `DETERMINISTIC_POINTS.md`, `NEIGHBORING_FIBERS.md`, `PARTIAL_EXPANSION.md`, `COMPLEXITY_AUDIT.md`, `COMPATIBILITY_MAP.md`). These are AI-assisted proof candidates and exact finite checks, not an independent referee report or a completed uniform deterministic exponent improvement. The uniform deterministic quarter-power ECDLP time bound remains **OPEN**.

## Summary

The 27 November 2026 fixed-shift manuscript leaves two routes toward a uniform deterministic time bound $2^{O(b/(\log b)^c)}$ ($c > 0$) under Liu–Zhang–Ren (LZR, ePrint 2026/3020):

- **Lever A (Fixed-generator noncancellation):** Either establish full rank for $O(1)$-slot-pinned unweighted shift Gram matrices $K_{A_1,\dots,A_p}$ summed over all group values $g \in G$ (the $2^{O(b/(\log b)^{1/4})}$-budget relaxation directly below `NEIGHBORING_FIBERS.md` Proposition 4), or construct a deterministic hitting stack for decoupled polynomial evaluations $B_{\mathrm{shift}}(x, y)$.
- **Lever B (Residual lookup contraction after Theorem 2):** Break the $2^{O(b/\sqrt{\log\log b})}$ bottleneck left after Theorem 2's one-step $R = \Theta(m)$ partial expansion.

This note establishes:

1. **Pinning barriers (Theorems A1 & A2):** Unconditional $O(1)$-slot external pinning $K_{A_1,\dots,A_p}$ fails identically modulo $r$ via Lucas's Theorem at $m = r - 1$ for unweighted $F_i = \{1, 2\}$ ($p \le r - 7$), **and** fails identically modulo $r$ at strict $4m < r$ ($m \ge p + 3$, $r = m d + 1$, $d \ge 5$, $|D_0| = r = 100\%$ coverage) for positive unit binomial weights $F_i(U) = U(1+U)^d$. No proof using only $4m < r$, $|D_0| = r$, and unit group-algebra structure in $\mathbb{F}_r[U]/((U-1)^r)$ can settle $K_{A,B}$.
2. **Sharp positive & empirical boundary for unweighted 0/1 pinning (Propositions A3 & A4):** For $m = 2$, $2m$-slot (`cpb = 1`) shifts, $\{K_{A,B} : A, B \in F_1\}$ has the exact kernel of all individual shift rows for every prime $r$ and all proper subsets $F_1, F_2 \subsetneq \mathbb{Z}/r\mathbb{Z}$. For general unweighted 0/1 blocks with $2\cdot\mathrm{cpb}\cdot m < r$, $K_1$ and $K_2$ had **0 rank drops** across exhaustive sweeps ($r = 7, m = 3$: 969 triples; $r = 11, m = 2$: 185 identical pairs), Paley quadratic-residue blocks, and elliptic curve controls ($r = 19, 37, 257$). Any positive proof for unweighted $K_{A,B}$ must therefore exploit either 0/1 integrality or Stage I's $\chi^2$/no-wrap flatness.
3. **Unified Lever B residual barrier (Theorems B1, B2 & B3):** Inside the Theorem 2 + arbitrary-lookup architecture, (i) single-level core splitting ($d = 1$) strictly dominates all multi-level schedules $(R_0, \dots, R_{d-1})$ and already reaches a sharp core-only ceiling of $2^{O(b/L^{1/6})}$ at $m = L^{1/6}$, $\log_2 R = \lceil L^{1/3}\rceil$ (passing all 9 guards of `PARTIAL_EXPANSION.md`), isolating the $4m$ degree-$\Theta(b/m)$ lookups as the sole obstruction; (ii) co-locating each slot's $k = \Theta(b/m)$ terminal hosts into one root group incurs correlated edge variance $\Theta(\beta k / \sqrt{D})$ in Lemma 5.1, violating the $(M/D)^2$ packing guard by $\Omega(D^{3/2}/m^3) \ge 4096 m^3$ even at $u = 1$; (iii) any binary contraction tree on the balanced $K_{R, 4m}$ residual bipartite graph cuts at least $(\frac{4}{9} - o(1))\beta = \Omega(b)$ bits.

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
so every valid Hamming weight $s$ satisfies $s \bmod r \in [r - 1 - k,\, r - 1]$. Because $k \le r - 5$, the base-$r$ units digit of $M$ is $M \bmod r = r - 4 - k < r - 1 - k \le s \bmod r$. Lucas's Theorem gives $\binom{M}{s} \equiv \binom{3}{\lfloor s/r \rfloor}\binom{r - 4 - k}{s \bmod r} \equiv 0 \pmod r$, killing every entry. For (3), set $u_{(0,2)} = 1$, $u_{(1,2)} = -1$, and all other entries $0$; every element $g \in \{1, 2\}$ needs $k_g = (g - 2m) \bmod r \in \{3, 4\}$ twos, which blocks $\{2, 3\}$ supply with $u = 0$ (placing one `2` in block $0$ for $(0, 2)$ and one `2` in block $1$ for $(1, 2)$). $\square$

### Theorem A2 (Small-$m$ Dense Binomial-Weighted Collapse at Strict $4m < r$)
Let $p \ge 0$ be an integer, let $m \ge p + 3$, and let $r = m d + 1$ be any prime with $d \ge 5$ (so $4m = \frac{4}{d}(r-1) < r$ strictly, accommodating LZR's asymptotic regime $m = \Theta((\log b)^{1/4}) \ll r$). Take $m$ identical blocks $F_1 = \dots = F_m = \{1, 2, \dots, d + 1\}$ with doubled tuples ($\mathrm{cpb} = 2$) and positive binomial weights $w(j + 1) = \binom{d}{j} \in \mathbb{F}_r^\times$, so $F_i(U) = U(1 + U)^d$, $F_i(1) = 2^d \not\equiv 0 \pmod r$, and $|D_0| = m d + 1 = r$. Then:
1. For every $0 \le p \le m - 3$ and every choice of $p$ pinned slots, the weighted shift Gram matrix $K_{A_1, \dots, A_p}$ over $\mathbb{F}_r$ is identically zero.
2. There exists a valid explicit matrix $E$ with $\mathrm{rank}(E) < N$ such that
   $$\mathrm{rank}\!\left(E \;\cup\; \bigcup_{p=0}^{m-3}\{K_{A_1,\dots,A_p}\}\right) = \mathrm{rank}(E) < N.$$

**Proof.** Each free $Y$-slot contributes $(1 + U)^d$ and each free $X$-slot contributes $U^{-d}(1 + U)^d$. With $k = k_X + k_Y \le p + 2 \le m - 1$ pinned slots, $M = (4m - k)d = 4(r - 1) - k d = 3r + (r - 4 - k d)$ has $M \bmod r = r - 4 - k d \ge d - 3 \ge 2$. Meanwhile, $2m d = 2(r - 1) \equiv -2 \pmod r$ forces $s \bmod r \in [r - 1 - k d,\, r - 1] > r - 4 - k d = M \bmod r$, so $\binom{M}{s} \equiv 0 \pmod r$ by Lucas's Theorem. For (2), $u_{(0,2)} = 1$, $u_{(1,2)} = -1$ (all other coordinates $0$) is annihilated by genuine $2m$-slot decompositions whose non-anchored shifts are absorbed by blocks $\{2, \dots, m - 1\}$, which have capacity $2(m - 2)d \ge 2(r - 1 - 2d) > r - 1$ for $m \ge 5$. $\square$

*Scope boundary:* Theorem A2 rules out any proof of full rank for $\{K_{A_1,\dots,A_p}\}$ ($p \le m - 3$) that uses only $4m < r$, $|D_0| = r$, $n_i < r$, and unit group-algebra structure in $\mathbb{F}_r[U]/((U-1)^r)$. The binomial family violates both unweighted 0/1 integrality and Stage I's $\ell_2$ flatness ($\chi^2 = \Theta(\sqrt{r})$ and pair counts $\approx 4^{r-1}/r \gg r$).

### Proposition A3 (Exact Forest Theorem for $m = 2$, $\mathrm{cpb} = 1$)
For every prime $r$ and every pair of nonempty proper subsets $F_1, F_2 \subsetneq \mathbb{Z}/r\mathbb{Z}$, the $2m$-slot two-pin stack $\{K_{A,B} : A, B \in F_1\}$ satisfies
$$\ker\big(\{K_{A,B} : A, B \in F_1\}\big) = \ker(\text{all individual } 2m\text{-slot shift rows}).$$

**Proof.** Fix $A, B \in F_1$ and put $t = P + A - B$. If $t = 0$, then $Y_2 = X_2$ for every shift pair, so $K_{A,B} = |F_2|(e_B - e_A)(e_B - e_A)^\top$; since $0 < |F_2| < r$, $K_{A,B} u = 0$ forces $(e_B - e_A)^\top u = 0$. If $t \ne 0$, the directed shift graph $S_t = \{(x, y) \in F_2^2 : y - x = t\}$ sits inside a single $r$-cycle ($r$ prime) on $|F_2| < r$ vertices, hence is a forest of directed paths. Its edge-incidence vectors $\{b_e = e_y - e_x : e \in S_t\}$ (supported on $F_2$, disjoint from $F_1$) are linearly independent over $\mathbb{F}_r$. Restricting $K_{A,B} u = \sum_{e \in S_t} (a_e^\top u) a_e = 0$ to $F_2$-coordinates gives $\sum_{e \in S_t} (a_e^\top u) b_e = 0$, hence $a_e^\top u = 0$ for every individual shift row $e \in S_t$. $\square$

### Proposition A4 (Exact Reduction of $K_{A,B} u = 0$)
For fixed blocks $(i_1, i_2)$ and pins $A \in F_{i_1}$ (on $X$) and $B \in F_{i_2}$ (on $Y$), put $t = P + A - B$ and let $n(t), v(t), G'(t)$ be the 0th, 1st, and 2nd moments of $c_{Y'} - c_{X'}$ over unpinned complements with $\mathrm{val}(Y') - \mathrm{val}(X') = t$. Using $\lambda^\top(c_Y - c_X) = 1$,
$$K_{A,B} u = 0 \iff \begin{cases} v(t)^\top u + (u_{i_2, B} - u_{i_1, A})\, n(t) = 0, \\ G'(t) u + (u_{i_2, B} - u_{i_1, A})\, v(t) = 0. \end{cases}$$

**Proof.** Expand $K_{A,B} = n(t)(e_B - e_A)(e_B - e_A)^\top + (e_B - e_A)v(t)^\top + v(t)(e_B - e_A)^\top + G'(t)$, left-multiply $K_{A,B} u = 0$ by $\lambda^\top$ to obtain the scalar equation, and substitute back into $K_{A,B} u = 0$. $\square$

When $n_i = 2^{\Theta(b/m)} \ll \sqrt{r}$, blocks are nearly Sidon, so distinct pairs $(A, B)$ rarely share $t = P + A - B$; thus $\lambda^\top K_{A,B} u = 0$ alone cannot compare $u_{i, B} - u_{i, A}$ across pairs (confirmed empirically: 1-row $\lambda^\top K_1$ drops on 66 dense triples at $r = 7, m = 3$, whereas full $K_1$ and $K_2$ have 0/969 drops).

### Decoupled $B_{\mathrm{shift}}(x, y)$ Evaluations (Regime A-I)
Splitting the single radix digit straddling the $X/Y$ cut yields $s_X + s_Y \le s + 1 = O(1)$ independent variables and drops $M_C$ (applied externally, with $E$ stacked directly). Against fiber-respecting dense $m = 2$ adversaries at $r = 7$, fixed coupled 2-point lists (`A1_coupled_L2`) drop to rank 4/6 on `[[1,2,3,4],[1,2,3,4]]` and 1-sided zero-projection 2-point lists (`A3_zero_L2`) drop to 5/6 on `[[1,3,5,6],[1,3,5,6]]`, while decoupled asymmetric 2-point lists (`A2_asym_L2`) hit full rank 6/6 on both. Constructing a worst-case deterministic hitting list of size $2^{O(b/L^{1/4})}$ remains open.

---

## 2. Lever B: Unified Residual Barriers for Theorem 2 + Arbitrary Lookups

Anchor all parameters in `PARTIAL_EXPANSION.md` and `COMPLEXITY_AUDIT.md`: $L = \log_2 b$, $S = \Theta(mb)$, $\lambda = \lceil \log_2 S \rceil$, $\beta = 4mk = \Theta(b)$, $\epsilon = 2^{-17}$, $u = \lfloor \epsilon \lambda / \log_2 R \rfloor$, $a = 2u$, $W = 2^a$, $D = R^{2a}$, $H = \lceil S / \lfloor \sqrt{\lambda}\rfloor \rceil$, $\Delta = \lceil S^{9/10}\rceil$, and residual contraction parameter $\rho = \max(1, \lfloor \log_2(\beta / \max(\Delta_{\mathrm{res}}, \lceil \beta^{9/10}\rceil)) \rfloor)$ with log bit cost $O(\beta / \sqrt{\rho} + n\lambda + S/\sqrt{\lambda})$.

### Theorem B1 (Single-Level Optimality, Sharp Core-Only $1/6$ Ceiling, and Lookup Floor)
1. **Multi-level suboptimality:** For any $d$-level distributed core split with group counts $(R_0, \dots, R_{d-1})$, total fan-out $R_{\mathrm{eff}} = \prod_{j=0}^{d-1} R_j$, and packing-guard caps $u_j \le c_0 \lambda / \log_2 R_j$,
   $$\log_2 |\mathcal{I}| = \Theta\!\left(\sum_{j=0}^{d-1}\frac{2^j S}{\sqrt{u_j}}\right) \ge \Omega\!\left(\frac{S}{\sqrt{\lambda}}\sum_{j=0}^{d-1} 2^j \sqrt{\log_2 R_j}\right) \ge \Omega\!\left(S\sqrt{\frac{\log_2 R_{\mathrm{eff}}}{\lambda}}\right),$$
   with strict inequality for every nontrivial $d \ge 2$ (by $2^j \ge 1$ and strict subadditivity of $\sqrt{\cdot}$).
2. **Core-only $1/6$ ceiling:** Without the $4m$ reattached lookups, balancing factor base $b/m$, root split $m b \sqrt{\log_2 R / L}$, and core-only residual contraction $b / \sqrt{\log_2 R}$ attains the unique optimum $\Theta(b / L^{1/6})$ at $d = 1$, $m = \Theta(L^{1/6})$, $\log_2 R = \Theta(L^{1/3})$. This schedule satisfies every guard of `PARTIAL_EXPANSION.md` (including $R^a \le (2S)^{2\epsilon}$, $D^2 \le (2S)^{8\epsilon}$, $2^u/u \gg \sqrt{\lambda}$, $K\ell = O(R m \lambda) = o(\beta)$, and digit-tuple overhead $O(R\lambda) = o(b/L^{1/6})$).
3. **Lookup floor:** Reattaching the $4m$ degree-$k = \Theta(b/m)$ lookups forces $\Delta_{\mathrm{res}} \ge k = \Theta(b/m)$, hence $\rho \le \log_2(4m) + O(1)$. Any core contraction beating $2^{\epsilon b}$ forces $m = o(\sqrt{L})$, so $\rho \le \frac{1}{2}\log_2 L + O(1)$ and $\log_2 T_{\mathrm{res}} = \Omega(b/\sqrt{\log_2 L})$ for every core schedule $(R_0, \dots, R_{d-1})$ and every 2-level factor-base subdivision $m \to m' = mt$. $\square$

### Theorem B2 (Slot-Co-Location Variance Barrier in Lemma 5.1)
If Lemma 5.1's root labeling over $[R]^a$ ($a = 2u$, $D = R^{2a}$) constrains each slot's $k = \Theta(b/m)$ terminal hosts to be partitioned into $k/s$ clusters of size $s$ sharing a common $[R]^a$ label (with distinct clusters and internal vertices labeled independently and uniformly):
1. Two edges incident to the same tied cluster (with distinct other endpoints) have summed oriented-class covariance $\sum_{\gamma \in [R]^{2a}} \mathrm{Cov}(\mathbf{1}_{e\in\gamma}, \mathbf{1}_{e'\in\gamma}) = \frac{1}{\sqrt{D}} - \frac{1}{D}$.
2. Summing over all $4m$ slots gives correlated edge variance $\ge 4m \frac{k}{s}\binom{s}{2}(\frac{1}{\sqrt{D}} - \frac{1}{D}) = \Theta(\beta s / \sqrt{D})$.
3. Lemma 5.1's squared-deviation packing bound $\mathrm{Var} \le O(M^2/D^2)$ forces $s \le O(M^2 / (\beta D^{3/2})) = O(m^2 b / D^{3/2})$. Since $R \ge \max(64, 4m)$ and $u \ge 1$ imply $D = R^{4u} \ge (4m)^4$, full slot co-location ($s = k = \Theta(b/m)$)—which would disconnect the residual graph into $R$ independent $O(b/m)$-bit pieces—violates the packing guard by $\Omega(D^{3/2}/m^3) \ge 4096 m^3$ even at $u = 1$ (and by $S^{\Omega(\epsilon)}$ at Theorem 2's tuning). Tying only a single coordinate out of $a$ gives covariance $(R - 1)/D$ and still violates the guard by $\Omega(R D / m^3) \ge 1024 m^2$. Every legal cluster size leaves $k / s = \Omega(D^{3/2}/m^3) \gg R$ independent clusters per slot. $\square$

### Theorem B3 ($4/9$-Cutwidth Barrier on Balanced Residual $K_{R, 4m}$)
In a balanced Theorem 2 residual bipartite graph $K_{R, 4m}$ ($R, 4m \ge 4$) where each stopped child $j \in [R]$ shares $\frac{\beta}{4mR}(1 \pm o(1))$ bits with each lookup $\ell \in [4m]$, every binary contraction tree contains an intermediate tensor with boundary bit-width at least $(\frac{4}{9} - o(1))\beta = \Omega(b)$.

**Proof.** Assign leaf weights $\frac{1}{R}$ (left) and $\frac{1}{4m}$ (right); the root has weight $2$ and every leaf has weight $\le \frac{1}{4} < \frac{2}{3}$. Walking down the heavier child from the root reaches a subtree with left/right fractions $(\alpha, \gamma) \in [0, 1]^2$ satisfying $s = \alpha + \gamma \in [\frac{2}{3}, \frac{4}{3}]$. The cut bit-fraction is
$$\alpha(1 - \gamma) + \gamma(1 - \alpha) = s - 2\alpha\gamma \ge s - \frac{s^2}{2} \ge \min_{s \in [2/3, 4/3]}\left(s - \frac{s^2}{2}\right) = \frac{4}{9}. \qquad \square$$

Together with `COMPATIBILITY_MAP.md`'s single-graph counting lower bound ($E_{\mathrm{sub}} \ge 2^{(1-1/t)b/m}/\mathrm{poly}(b)$), Theorems B1–B3 close the Theorem 2 + universal arbitrary-lookup pipeline at $2^{O(b/\sqrt{\log\log b})}$.

---

## 3. Reproducibility

Both scripts require only Python 3.10+ standard library (run from the repository root):

```sh
python3 -B research/deterministic_endgame_barriers/check_pinning_barriers.py
python3 -B research/deterministic_endgame_barriers/check_residual_barriers.py
```

Compare outputs against `research/deterministic_endgame_barriers/verification.txt`. These scripts verify the exact finite-field algebra, Lucas collapse instances, Forest Theorem, rational covariances, bipartite cutwidth bounds, and asymptotic guard tables; they do not implement LZR's tensor contractor.

## 4. Open Frontier

The single highest-leverage open question remaining on `main` is whether the **unweighted 0/1** two-pin stack $\{K_{A,B}\}$ (or 1-pin $\{K_A\}$) always attains full shift rank $N_0 - m$ under $4m < r$ and Stage I's density + $\chi^2$/no-wrap flatness bounds. By Theorem A2, any proof must use 0/1 integrality or Stage I flatness (not mere unit group-algebra structure); by Proposition A4, it must use the vector equation $G'(P+A-B)u + (u_{i_2,B}-u_{i_1,A})v(P+A-B) = 0$.
