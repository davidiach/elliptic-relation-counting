# Single-target recovery from pinned target counts

27 November 2026. Research follow-up to the [reduced-floor examples](../sparse_target_runs/RESEARCH.md); baseline `main` at `13d2803`. This note formalizes the pinned single-target matrix $M_t$ used there. It proves that the explicit consistency rows plus **one** pinned matrix with the forced nullity determine every logarithm. It also gives an always-correct adaptive algorithm whose running time depends on the first such target, proves a counting bound for failures, and records two barriers. No uniform deterministic quarter-power bound is proved: the remaining step is a conjecture (§7). Resource claims assume LZR ePrint 2026/3020 Theorem 2.2 (closed Boolean contraction), its gadget lemmas and the Stage I seed of its Theorem 3.1, as in [the main theorem](../../MAIN_RESULT.md). The linear algebra is unconditional. The note was prepared with AI assistance. Its proofs have not been independently verified by a human.

## 1. Pinned rows

Let $G=E(\mathbb F_p)$ have prime order $r$ with generator $P$. Let $F_1,\dots,F_m\subseteq G\setminus\{O\}$ be the Stage I blocks. Columns are the tagged elements $(i,A)$, $A\in F_i$ ($N_0$ in total), plus $\star_P,\star_Q$, so $N=N_0+2$. Points within a block are distinct; duplicates across blocks are separate columns. A *representation* $\mathbf X$ chooses two slots $X_{i,1},X_{i,2}\in F_i$ in every block. Its incidence vector $\mathrm{cnt}(\mathbf X)\in\mathbb Z^N$ counts slot occurrences and is zero on $\star_P,\star_Q$. Let $\lambda$ be the true logarithm vector, with $\lambda_{\star_P}=1$ and $\lambda_{\star_Q}=\log_P Q$. Then $\mathrm{cnt}(\mathbf X)\cdot\lambda\equiv\log_P\mathrm{val}(\mathbf X)$.

*Explicit rows* are the normalization $e_{\star_P}$ with right side $1$, and for each column $j$ with point $g_j$ the row $e_j-\mathrm{cnt}(\mathbf Z_j)$ with right side $0$. Here $\mathbf Z_j$ is any representation of $g_j$ (a *decomposition rule*), which exists by LZR Theorem 3.1.

For a target $t\in\mathbb Z/r$ and a pin $h=(i,A)$ put

$$
R_{t,h}=\sum_{\mathrm{val}\,\mathbf X=[t]P,\ X_{i,1}=A}\mathrm{cnt}(\mathbf X),\qquad
n_h(t)=\#\{\mathbf X:\mathrm{val}\,\mathbf X=[t]P,\ X_{i,1}=A\}.
$$

Every such row is valid: $R_{t,h}\cdot\lambda\equiv t\,n_h(t)\pmod r$. Both sides are pinned counts, so no logarithm is needed to form them. Let $M_t$ be the $N_0\times N_0$ matrix of these rows on the block columns (target rows vanish on $\star_P,\star_Q$). This is the matrix of the [reduced-floor examples](../sparse_target_runs/RESEARCH.md), where $A_t^{\mathsf T}A_t=2M_t$. The verifier checks that the two constructions agree.

**Entry formula (analysis only).** Identify $G$ with $\mathbb Z/r$ through $\log_P$ and work in $\mathbb F_r[G]$. Put $\delta_i=\sum_{A\in F_i}[A]$ and $E_i=\delta_i\prod_{k\ne i}\delta_k^2$, the completions of the other $2m-1$ slots. Let $\varphi_{i,k}$ be $E_i$ with one factor $\delta_k$ removed. Counting slot by slot gives

$$
R_{t,(i,A)}[(k,B)]=[(k,B)=(i,A)]\,E_i(t-A)+c_k\,\varphi_{i,k}(t-A-B),\qquad c_i=1,\ c_{k\ne i}=2,
$$

where $f(g)$ denotes the coefficient of $[g]$, and $n_{(i,A)}(t)=E_i(t-A)$. The algorithm never uses this identification; it computes counts.

**Cost.** $R_{t,h}[(k,u)]$ is a sum over the two slots of block $k$ of counts of representations with value $[t]P$ and at most two address pins. This is the decomposition-counting network of LZR Theorem 3.1 with a fixed terminal. It has $O(mb)$ edges and degree $O(b^{3/4})$, like the pinned counts in the [bit-cost section](../../MAIN_RESULT.md#bit-cost-of-a-single-evaluation) with $2m$ slots and no weight pieces. Under LZR Theorem 2.2 one target costs $(2N_0N+N_0)\,2^{O(mb/\sqrt{\log b})}$. This inherits the LZR gadget bounds (Lemmas 2.3–2.5) and is not re-derived here.

## 2. All targets suffice

**Lemma 1.** An element of $\mathbb F_r[G]$ is a unit iff its augmentation is nonzero. In particular $[A]$, $\delta_i$ and $E_i$ are units.

*Proof.* $\mathbb F_r[\mathbb Z/r]\cong\mathbb F_r[U]/(U-1)^r$, and $U-1$ is nilpotent, so $x=\varepsilon(x)+n$ with $n$ nilpotent. The augmentations of $\delta_i$ and $E_i$ are products of the $|F_i|\in[1,r-1]$. ∎

**Theorem 2 (all-target completeness).** If $u\in\mathbb F_r^{N_0}$ satisfies $M_tu=0$ for all $t\in\mathbb Z/r$, then $u$ is constant on each block, with values $\kappa_i$ and $\sum_i\kappa_i=0$. That is, $\bigcap_t\ker M_t=\mathcal B_0$, the space of such vectors. Together with the explicit rows, the pinned rows at all targets determine $\lambda$.

*Proof.* Write $\tilde u_k=\sum_{B\in F_k}u(k,B)[B]$. Summing the value of each free slot over completions gives

$$
\sum_t (R_{t,(i,A)}\cdot u)[t]=[A]\,E_i\bigl(u(i,A)+\Psi_i\bigr),\qquad
\Psi_i=\tilde u_i\delta_i^{-1}+2\sum_{k\ne i}\tilde u_k\delta_k^{-1}.
$$

Here slot $(i,1)$ contributes $u(i,A)[A]E_i$, slot $(i,2)$ contributes $[A]\tilde u_iE_i\delta_i^{-1}$, and each of the two slots of block $k\ne i$ contributes $[A]\tilde u_kE_i\delta_k^{-1}$. If the left side vanishes, Lemma 1 gives $\Psi_i=-u(i,A)\cdot[O]$ for every $A\in F_i$. Since $\Psi_i$ does not depend on $A$, $u(i,\cdot)\equiv\kappa_i$. Then $\tilde u_i=\kappa_i\delta_i$ and $\Psi_i=(\kappa_i+2\sum_{k\ne i}\kappa_k)[O]$, so $-\kappa_i=\kappa_i+2\sum_{k\ne i}\kappa_k$. This means $2\sum\kappa=0$, so $\sum\kappa=0$ because $r$ is odd. Conversely $R_{t,h}\cdot\kappa=2n_h(t)\sum_i\kappa_i$ for block-constant $\kappa$, since each block has two slots.

For the second claim, let $w$ annihilate every explicit row homogeneously and every pinned row. Its block part is $\kappa\in\mathcal B_0$. The explicit row of a block column in $F_i$ gives $\kappa_i=\mathrm{cnt}(\mathbf Z_j)\cdot w=2\sum\kappa=0$. The normalization gives $w_{\star_P}=0$, and the row of $\star_Q$ gives $w_{\star_Q}=\mathrm{cnt}(\mathbf Z_{\star_Q})\cdot w=0$. ∎

Pins are needed. Without them the same computation gives $2D\sum_i\tilde u_i\delta_i^{-1}$ with $D=\prod_i\delta_i^2$. If $\delta_2=[a]\delta_1$, then $\tilde u_1=\theta$, $\tilde u_2=-[a]\theta$ is in the kernel for every $\theta$ supported on $F_1$.

## 3. One target: the forced kernel and recovery

**Lemma 3 (forced kernel).** For every $t$, $\ker M_t\supseteq V_t:=\mathcal B_0\oplus\langle u_t\rangle$, where $u_t=\lambda|_{\rm blocks}-t(2m)^{-1}\mathbf 1$. If some block has two points, $\dim V_t=m$, so $\operatorname{nullity}M_t\ge m$.

*Proof.* Every representation has $2m$ slots, so $R_{t,h}\cdot\mathbf 1=2m\,n_h(t)$. With validity, $R_{t,h}\cdot u_t=t\,n_h(t)-t\,n_h(t)=0$. $\mathcal B_0\subseteq\ker M_t$ was shown above. Finally $u_t$ is not block-constant, because $\lambda$ takes distinct values on distinct points of a block. ∎

**Definition.** $\mathrm H_1(t)$ holds if $\operatorname{nullity}_{\mathbb F_r}M_t=m$, i.e. $\ker M_t=V_t$.

**Theorem 4 (single-target recovery).** Let $r>2m$. If $\mathrm H_1(t)$ holds, then for **every** decomposition rule the explicit rows together with the $N_0$ rows $R_{t,h}$ have column rank $N$. So $\lambda$, and in particular $\log_PQ=\lambda_{\star_Q}$, is the unique solution of that linear system.

*Proof.* Let $w$ annihilate these rows homogeneously. The pinned rows vanish on $\star_P,\star_Q$, so $w|_{\rm blocks}\in\ker M_t=V_t$. Since $u_t$ and $\lambda|_{\rm blocks}$ differ by a block constant, $w|_{\rm blocks}=\kappa+a\lambda|_{\rm blocks}$ with $\kappa$ block-constant and $a\in\mathbb F_r$. Put $w'=w-a\lambda$. Each consistency row $e_j-\mathrm{cnt}(\mathbf Z_j)$ annihilates $\lambda$, hence $w'$, and $w'|_{\rm blocks}=\kappa$. For a block column $j\in F_i$ this gives $\kappa_i=\mathrm{cnt}(\mathbf Z_j)\cdot w'=2\sum_k\kappa_k$. So all $\kappa_i$ equal a common value $\kappa_0$ with $(2m-1)\kappa_0=0$, i.e. $\kappa_0=0$. The rows of $\star_P,\star_Q$ then give $w'_{\star_P}=w'_{\star_Q}=0$, so $w=a\lambda$. The homogeneous normalization gives $0=w_{\star_P}=a$. ∎

One target never satisfies the *hard form* $\bigcap_{t\in T}\ker M_t=\mathcal B_0$, because of $u_t$. The explicit rows remove exactly the extra direction $\lambda$ modulo block constants, so one target can still complete the *stated* system (explicit rows plus pinned rows). Recovery is direct: solve the explicit rows on $\ker M_t\oplus\langle\mathbf 1\rangle\oplus\langle e_{\star_P},e_{\star_Q}\rangle$, which has dimension $m+3$ and contains $\lambda$.

## 4. An adaptive, always-correct algorithm

**Theorem 5.** Assume LZR Theorem 2.2 and the cost statement of §1. Fix an explicit enumeration $t_1,t_2,\dots$ of $\mathbb Z/r$, for example $t_k=k-1$. Run Stage I; form the explicit rows; for $k=1,2,\dots$ append all $R_{t_k,h}$ and maintain an echelon form modulo $r$; stop at the first $K^\ast$ with rank $N$; solve. This algorithm is deterministic, always outputs $\log_PQ$, and satisfies

$$
K^\ast\le\min\{k:\mathrm H_1(t_k)\}\le X+1,\qquad X:=\#\{t\in\mathbb Z/r:\mathrm H_1(t)\text{ fails}\}.
$$

Its bit time is $K^\ast\cdot2^{O(b/(\log b)^{1/4})}$ at $m\asymp(\log b)^{1/4}$, and its space is $2^{O(b/(\log b)^{1/4})}$.

*Proof.* All rows are valid, so at full column rank the unique solution is $\lambda$. By Theorem 2 the rank reaches $N$ by $k=r$. By Theorem 4 it reaches $N$ as soon as some $t_k$ with $\mathrm H_1(t_k)$ has been added, and among $X+1$ distinct targets at least one satisfies $\mathrm H_1$. Stage I costs $N\,2^{O(b/m+mb/\sqrt{\log b})}$. Each target adds the cost of §1 and $O(N_0N^2)$ field operations. Since $\log N=O(b/m)$, balancing as in the main theorem gives the bound. ∎

The worst case $K^\ast=r$ is worse than baby-step giant-step. The algorithm meets the quarter-power bound exactly when $K^\ast=2^{O(b/m)}$. $\mathrm H_1(t)$ involves neither $Q$ nor the decompositions and is *checked* at run time; a failure costs time, not correctness.

## 5. A counting bound for failures

Let $W$ be a valid system containing the explicit rows. Its *residual* $\mathcal R$ is the block projection of its homogeneous kernel, and $d=\dim\mathcal R$. Let $S=\{\text{block constants}\}\oplus\langle\lambda|_{\rm blocks}\rangle$, so $\dim S=m+1$. Write $N_t=\ker M_t$, $\widehat N_t=(N_t+S)/S$ and $\hat e_t=\dim\widehat N_t$.

**Lemma 6.** (a) The block projection is injective on the homogeneous kernel of $W$, and $\mathcal R\cap S=0$. Hence $W$ together with the rows of a target $t$ has rank $N$ iff $\mathcal R\cap N_t=0$; in that case say that $t$ *completes* $W$. (b) If $n_h(t)\not\equiv0$ for some pin $h$, then $N_t\cap S=V_t$, so $\hat e_t=\operatorname{nullity}M_t-m$.

*Proof.* (a) Let $u$ be a homogeneous solution with $u_b=\kappa+a\lambda_b\in S$. The row of a block column $j\in F_i$ gives $\kappa_i+a\lambda_j=2\sum\kappa+a\lambda_j$, since $\mathrm{cnt}(\mathbf Z_j)\cdot\lambda_b=\lambda_j$. So $\kappa\equiv\kappa_0$ with $(2m-1)\kappa_0=0$, i.e. $\kappa=0$. The normalization gives $u_{\star_P}=0$, and the row of $\star_P$ gives $u_{\star_P}=\mathrm{cnt}(\mathbf Z_{\star_P})\cdot u_b=a$, so $a=0$. If $u_b=0$, the rows of $\star_P,\star_Q$ give $u_{\star_P}=u_{\star_Q}=0$. The homogeneous kernel of $W\cup\{R_{t,h}\}$ projects isomorphically onto $\mathcal R\cap N_t$. (b) For $v=\kappa+a\lambda_b$, $R_{t,h}\cdot v=(2\sum\kappa+at)\,n_h(t)$. So $N_t\cap S$ is the hyperplane $2\sum\kappa=-at$ of $S$, which has dimension $m$ and contains $V_t$. ∎

**Theorem 7 (failure inequality).** For a finite target set $T$, let $\delta(T)=\sum_{t\in T}\hat e_t-\dim\sum_{t\in T}\widehat N_t$ be the kernel dimension of the sum map $\Sigma:\bigoplus_{t\in T}\widehat N_t\to\mathbb F_r^{N_0}/S$. Then

$$
\#\{t\in T:\ t\text{ does not complete }W\}\le\min\bigl(\#\{t\in T:\hat e_t\ge1\},\ d+\delta(T)\bigr).
$$

*Proof.* Let $F\subseteq T$ be the failing targets and $s=|F|$. For $t\in F$ choose $0\ne v_t\in\mathcal R\cap N_t$. Since $\mathcal R\cap S=0$, the class $\hat v_t\in\widehat N_t$ is nonzero, which proves the first bound. The vectors $\iota_t(\hat v_t)$, $t\in F$, lie in distinct summands and span a space $Y$ of dimension $s$. Now $\Sigma(Y)\subseteq(\mathcal R+S)/S$, which has dimension $d$. So $s\le d+\dim(Y\cap\ker\Sigma)\le d+\delta(T)$. ∎

Three things must happen for a target to fail: $\mathrm H_1$ fails there, the residual is aligned with the excess kernel, and excess spaces of different targets are dependent. A related elementary fact: if the hard form fails for $T$, then for $v\in\bigcap_{t\in T}N_t\setminus\mathcal B_0$ the valid system $W=\mathrm{span}(\text{all valid rows})\cap v^\perp$ has residual $\langle v\rangle$ and is completed by no $t\in T$. So arguments that use only the residual need the hard form.

## 6. Barriers

**Theorem 8 (Chevalley–Warning barrier for fixed weightings).** Let $r$ be an odd prime and $w^{(1)},\dots,w^{(h)}\in\mathbb F_r^M$ arbitrary, with $M\ge N+2h$. Some $A\in\mathbb F_r^{M\times N}$ of rank $N$ has $\det\bigl(\sum_ew^{(k)}_ea_ea_e^{\mathsf T}\bigr)=0$ for every $k$.

*Proof.* Take a set $E'$ of $M-N+1>2h$ row indices. The $h$ diagonal quadratic forms $Q_k(\tau)=\sum_{e\in E'}w^{(k)}_e\tau_e^2$ have total degree $2h<|E'|$. By Chevalley–Warning the number of common zeros is divisible by $r$. Since $\tau=0$ is one, there is a common zero $\tau\ne0$. Give row $e\in E'$ the entry $\tau_e$ in column 1, and let the other $N-1$ rows be $e_2,\dots,e_N$. Then $\operatorname{rank}A=N$, while the first row of each weighted Gram matrix vanishes: its $(1,1)$ entry is $Q_k(\tau)$ and its other entries are zero. ∎

So a fixed set of evaluation points, a fixed stack of Gram weightings, or univariate projections tested as identities cannot certify rank uniformly over arrays. A derandomization must use the arithmetic of the actual count rows. This extends the remark in [adaptive rank](../short_relations_adaptive_rank/RESEARCH.md).

**Digit family (abstract model).** Work with logarithms: $G=\mathbb Z/r$ and arbitrary blocks. Let $D=2n-1$ and let $r$ be a prime with $D^m-\Delta\le r\le D^m$, $0\le\Delta\le D^m/4$. Put $F_i=s_i+D^i\{0,\dots,n-1\}$, with offsets such that $0\notin F_i$ and $\sigma:=2\sum_is_i\equiv1-\Delta$. A representation has value $\sigma+\Sigma$ with $\Sigma=\sum_iD^ie_i$ and $0\le e_i\le D-1$, so the $2m$-fold sums cover $\mathbb Z/r$.

**Theorem 9 (digit obstruction).** Column $(m-1,j)$ occurs in no representation of value $t$ for $1\le t\le D^{m-1}j-\Delta$. Consequently, for $1\le t\le D^{m-1}(n-1)-\Delta\approx r/2$, the coordinate vector of column $(m-1,n-1)$ lies in $\ker M_t$. So for $n\ge3$, $\mathrm H_1(t)$ fails at every such $t$, and the hard form fails for these consecutive targets.

*Proof.* Let $\tau=t-1+\Delta$. In the stated range $\Delta\le\tau<r$. Since $\Sigma\equiv\tau$ and $0\le\Sigma\le D^m-1\le r+\Delta-1<\tau+r$, we get $\Sigma=\tau$. If column $(m-1,j)$ fills a slot of the top block, then $e_{m-1}\ge j$, so $\tau\ge D^{m-1}j$, i.e. $t\ge D^{m-1}j+1-\Delta$. Every row $R_{t,h}$ then vanishes on column $(m-1,n-1)$, so that coordinate vector is in $\ker M_t$. It is not in $S$ when $n\ge3$, so the nullity exceeds $m$. ∎

**Lemma 10 (the digit family is not LZR-admissible).** Suppose the single-slot sumset $\Sigma_1=F_1+\dots+F_m$ covers more than $7r/8$ elements. Then for every pin $(i,A)$ the targets $t$ with $n_{(i,A)}(t)>0$ form more than $7r/8$ residues. In the digit family, $|\Sigma_1|\le n^m$ and $n^m/((2n-1)^m-\Delta)\le\frac43(n/(2n-1))^m\le\frac{16}{27}<\frac78$ for $m\ge2$, $n\ge2$. So the digit family violates the coverage guarantee.

*Proof.* Fix $a_k\in F_k$ for the slots $(k,2)$, $k\ne i$. Letting slot $(i,2)$ and the slots $(k,1)$, $k\ne i$, vary gives values $A+\sum_{k\ne i}a_k+\Sigma_1$, a translate of $\Sigma_1$. ∎

**Check needed.** Lemma 10 reads LZR's Stage I guarantee as "the sumset $F_1+\dots+F_m$ covers more than $7/8$ of $G$". The [LZR overview of Stage I and the proof of Theorem 3.1](https://eprint.iacr.org/2026/3020.pdf#page=11) derive this from $\Phi<64/57$ by Cauchy–Schwarz. This reading should be checked against the exact statement of LZR Theorem 3.1. The digit family is therefore an abstract-model obstruction only. It shows that no argument uniform over all covering block systems can bound $K^\ast$ for target progressions, but it says nothing about LZR seeds.

## 7. Conjecture: few exceptional targets

**Conjecture 11 (open).** For every seed produced by LZR Stage I ($\Phi<64/57$), $X\le2^{O(b/m)}$. A strong form is $X\le N_0$.

With Theorem 5, Conjecture 11 would give a uniform deterministic $2^{O(b/(\log b)^{1/4})}$ algorithm, conditional on LZR Theorem 2.2 and the Stage I imports. It is a counting statement with no explicit target set, so it needs no separate derandomization. Heuristically, $M_t$ with its forced kernel removed behaves like a random square matrix over $\mathbb F_r$, with extra nullity at about one target in $r$. **This is not a proof**, and symmetric block systems break the heuristic.

| Instance | $r$ | $m$ | $N_0$ | $X$ | Stated form at $t=0$ |
| --- | ---: | ---: | ---: | ---: | --- |
| Reduced-floor affine seed, $p=1009$, $\Phi=1.0447$ | 1013 | 2 | 257 | 2 (targets 642, 965) | rank $N$ at **every** target |
| [Sparse witness](../sparse_target_runs/witnesses.json), four-target run | 223 | 5 | 36 | 4 | rank $N$; explicit rank 34/38 |
| Sparse witness, five-target run | 67 | 4 | 32 | 5 | rank $N$; explicit rows already full |
| Sparse witness, five-target run | 109 | 5 | 30 | 5 | rank $N$; explicit rank 30/32 |
| Sparse witness, six-target run | 67 | 4 | 32 | 6 | rank $N$; explicit rows already full |
| Sparse witness, six-target shared vector | 109 | 6 | 30 | 12 | rank $N$; explicit rank 26/32 |
| Sparse witness, six-target shared vector | 271 | 7 | 34 | 12 | rank $N$; explicit rank 35/36 |
| Subgroup family $F_1=F_2=H$, $\lvert H\rvert=90$ (abstract) | 541 | 2 | 180 | **91** $\approx N_0/2$ | not tested |
| Digit family $n=17$ (abstract, not admissible) | 1087 | 2 | 34 | 1087 | not tested |

Stated-form results use the lexicographically first decomposition rule, $Q=[2]P$ and $Q=[r-3]P$. In each sparse witness $\mathrm H_1(0)$ fails, yet the explicit rows plus $M_0$ have rank $N$ and return the true logarithm. In the two order-67 witnesses the explicit rows alone already have rank $N$, so those cases are vacuous. In the other four, $M_0$ supplies the missing rank. Under Theorem 7, an excess at $t$ defeats a target only when the residual is aligned with it.

The subgroup family has blocks invariant under multiplication by $H\le\mathbb F_r^\times$, so $\operatorname{nullity}M_t=\operatorname{nullity}M_{ht}$. Its exceptional set is $\{0\}$ plus one orbit of size 90. The single-slot sumset $F_1+F_2$ covers all of $\mathbb Z/541$, and the energy ratio $rC/W^2=1.1065<64/57$. So in the abstract admissible class $X$ can be a constant fraction of $N_0$; the statement $X=O(1)$ is false there. Blocks from curves have no known multiplicative symmetry of logarithms, but that is not a proof. The six-target examples of [the reduced-floor note](../sparse_target_runs/SIX_TARGETS.md) use two translated copies of one two-point block, a symmetry of the same kind. Their own bound excludes two-point blocks under LZR's original floor.

## 8. Status

| Result | Status |
| --- | --- |
| Theorem 2, Lemma 3, Theorem 4, Lemma 6, Theorem 7, Theorem 9, Lemma 10 | Proved here; elementary linear algebra. Novelty not checked against the literature. |
| Cost statement (§1) and Theorem 5's time bound | Conditional on LZR Theorem 2.2, gadget Lemmas 2.3–2.5 and Theorem 3.1; not re-derived. |
| Theorem 8 | Proved; extends an existing remark in this repository. |
| Lemma 10's use of LZR coverage | Needs a check against LZR Theorem 3.1. |
| Conjecture 11 | **Open.** Numerical evidence at $r\le1087$ only. |

## Reproduction

Run from the repository root with NumPy and python-flint installed.

```sh
python -B -m research.single_target_recovery.check_single_target
```

The verifier checks the entry formula against direct enumeration of all tuples on two small curves at every target. It also checks the forced kernel and the all-target nullity $m-1$. It recomputes the exact potential of the [$p=1009$ seed](fixture_1009.json), scans all 1013 targets and tests single-target recovery at every target. For each sparse witness it confirms agreement with the construction in `sparse_target_runs`, counts $X$ over all targets and solves the stated system at $t=0$. It also covers the digit, subgroup and Chevalley–Warning examples. The run takes about 30 seconds. [Verifier](check_single_target.py) · [Recorded output](verification.txt)
