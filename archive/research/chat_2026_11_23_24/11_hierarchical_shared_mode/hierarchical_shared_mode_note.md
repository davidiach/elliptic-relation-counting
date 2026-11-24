# Hierarchical packet packing and shared-mode Newton contraction

24 November 2026 — continuation of `group_tensor_factorization_note.md`, `mode_contraction_note.md`, `contractor_degree_note.md`, and `section_norm_research_note.md`.

**Outcome.** We prove a shared-mode Newton interpolation identity that removes the requirement that each crossing mode have exactly one target group. It gives the same order of first-stage weighted enumeration as the ordinary Cauchy-packet contractor, with polynomial-time new packet entries rather than explicit coefficient extraction. This permits a larger class of hierarchical packet packings. For fixed labels in that class we give an exact minimum-padding algorithm, and exhibit packings with exponentially fewer packets than the flat unordered-label construction. On complete bipartite Fourier boundaries, however, a new entropy theorem prevents arbitrary reuse of group-label information across dyadic levels, even with bounded padding. A variable-precision lower bound and a matching taper (under its stated regularity conditions) then explain the cancellation of the two logarithmic savings for the ordinary two-stage accounting. A separate matrix-shift geometry produces many groups at one interpolation level but has full-rank first crossing coefficients. The strict general-prime ECDLP exponent improvement remains **unproved**. The previous 11.31486814 leading-coefficient candidate is unchanged and is not independently audited here.

Sources: the public quantitative Boolean tensor-contraction proof, *The Exponential Time Hypothesis is False*, 25 October 2026, `https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex`; its interface theorem in `quantitative-reduction.tex`; and the preceding project notes named above. We use the public Cauchy deformation, dyadic integrality and geometric expansion as stated there. The new Newton identity, packing and entropy arguments below are independent finite algebraic/combinatorial statements. Their algorithmic transfers retain the public contractor's precision, degree and entry-access assumptions. No claim of novelty against all unpublished work is made.

## 1. Shared crossing modes without target-dependent duplication

Let `A` be an unramified dyadic coefficient ring modulo `2^P`, possibly with an additional truncated polynomial variable. We can work in a larger unramified extension if necessary, with a residue field containing at least `P` distinct elements. Choose knots `xi_0,...,xi_(P-1)` with all pairwise differences units. All operations, including inverse differences, are polynomial-time in the explicit parameters. The original ring embeds in such an extension; the output of an exact identity with inputs in the original ring is the original element.

There are finitely many modes `Z_j` of positive integer weights `w_j`. Each has a source group `s(j)` and a source form `L_j(X)` that depends only on that group's bits. The same mode is allowed to enter arbitrarily many target groups. Fix Boolean row and column port assignments `X,Y`. Suppose the deformed row product is

```
A_(X,Y)(Z) = product_g A_g(Z; X_g,Y_g),
```

where `A_g` is a product of affine row factors with local part depending only on group `g` and crossing coefficients for `Z_j` divisible by `2^w_j`. In particular, the monomial coefficient of `Z^beta` in the product is divisible by `2^wt(beta)`, with `wt(beta)=sum_j w_j beta_j`. Thus in `A[Z]` the support is in the lower set

```
K = { beta>=0 : wt(beta)<P }.
```

Let `B_g(X_g,Y_g)` be the product of original tensors in group `g`, which can be entirely arbitrary subject to charged entry access. The desired contraction is the sum over `X,Y` of `product_g B_g` times `A_(X,Y)((L_j(X))_j)`.

### The identity

Define univariate monic Newton polynomials

```
N_b(Z) = product_(0<=h<b) (Z-xi_h),   N_0=1,
d_(b,a) = [product_(0<=h<=b, h!=a) (xi_a-xi_h)]^(-1),  0<=a<=b.
```

For any polynomial `A(Z)` supported on `K`, there is an exact identity

```
A(Z) = sum_(beta in K) sum_(0<=alpha<=beta)
         A((xi_(alpha_j))_j) product_j d_(beta_j,alpha_j) N_(beta_j)(Z_j).       (1)
```

The inequality is coordinatewise. To prove (1), the products of monic Newton polynomials with indices in `K` form a basis of the lower-set polynomial module: the monomial change of basis is unit triangular, and `K` is closed under decreasing any coordinate. The univariate divided-difference functional at the first `b+1` knots is `sum_(a=0)^b d_(b,a) A(xi_a)`. Applied to `N_c`, it is zero for `c<b` by the usual degree argument in Lagrange interpolation, one for `c=b` by its monic leading coefficient, and zero for `c>b` because all the evaluations vanish. Tensoring these functionals extracts the desired Newton coefficient for each `beta`. Every needed grid index `alpha<=beta` is in `K`. This proof uses only unit differences, so it is valid over the finite dyadic ring rather than only over a characteristic-zero field.

Substitute the actual source forms into (1). All factors then separate by groups after fixing `alpha,beta`. Define

```
Z_g(alpha,beta) = sum_(X_g,Y_g) B_g(X_g,Y_g)
          A_g((xi_(alpha_j))_j; X_g,Y_g)
          product_(j:s(j)=g) N_(beta_j)(L_j(X)).
```

The result is

```
F = sum_(beta in K) sum_(alpha<=beta)
        [product_j d_(beta_j,alpha_j)] product_g Z_g(alpha,beta).             (2)
```

Crucially, target groups use the same evaluation index `alpha_j`, and the source group uses `beta_j`. There is no separate exponent variable for each use of the mode. No hypothesis about the rank of the actual group tensor is needed. All choices of `X,Y` are included; (2) is an exact sum identity, not a heuristic on random inputs.

### Charged complexity

Suppose the number of modes of any one weight is at most `c`. There are exactly

```
J = sum_(beta in K) product_j (beta_j+1)
  = #{(alpha,gamma)>=0 : sum_j w_j(alpha_j+gamma_j)<P}
```

summands before the group contractions; set `gamma=beta-alpha`. This is an ordinary two-color weighted partition bound. For `theta>0`,

```
J <= exp(theta P) product_(k>=1) (1-exp(-theta k))^(-2c)
  <= exp(theta P + pi^2 c/(3 theta)).
```

The last inequality follows by expanding the logarithm and using `exp(theta h)-1>=theta h`. Choose `theta` proportional to `sqrt(c/P)` to get `log J=O(sqrt(cP))`. Polynomial factors and bounded exceptional parameters are harmless. Enumerate the pairs with polynomial delay by recursively choosing the two-color exponents subject to remaining weighted budget; all partial choices can be completed with zeros. The bound does **not** grow with the number of target groups. It changes the absolute partition constant (two colors instead of one).

In (2), after fixing the indices, each new packet-group tensor is evaluated by substituting public knots into its local row factors, and evaluating a product of univariate Newton factors on its source forms. Every such operation is polynomial in the explicit size and precision; no expansion of the packet row product into all mode monomials is needed. The original tensors keep their inherited entry-access costs and their original creation-precision representatives. The child graph is still formed from original vertices and packet-group vertices with one edge per retained original port; new vertex degree is at most the packet's row and column port count. Thus a first-level call has logarithmic cost bounded by `O(sqrt(cP))` plus the maximum paid child-contraction log cost and polynomial overhead. The exact reconstruction/readout and recursion invariants are still required.

### What is not yet compressed

This identity does not automatically give the pair-factor second-stage contractor from the previous checkpoint. The evaluation index of a shared mode may be incident to many target groups. Copying that index along ordinary equality wires charges its fanout, and could remove a second-stage saving. One can use a hyperedge description, but the generic Boolean contractor would then need a paid conversion or a separate hypergraph contraction result. Equation (2) removes target-dependent duplication in the *first-stage enumeration* and improves new packet-entry access; it does not imply an inexpensive second-stage join.

## 2. Hierarchical packet labels and the exact packing algorithm

Use a packet with `Lp=2^a` positions numbered `0,...,Lp-1`. A row port at position `x` is local to a column port at position `y` at an interpolation center `t` when `y-x=t` (within the actual finite interval). At valuation level `i=v2(t)<a`, choose row and column group labels for vertices. A hierarchical packet is allowed a different unordered pair of labels for each fixed low-bit prefix `x mod 2^i`. Within the two halves of that prefix, determined by bit `i`, the two orientations are reversed. The label at that level depends only on the first `i+1` low bits, not on higher bits. Unlike flat packetization, pairs at a later level can depend on the earlier prefix. Each dummy slot has virtual labels even if there is no graph vertex there.

For a full packet the condition that every allowed shift of valuation `i` localizes is equivalent to the hierarchical swap condition. Indeed, for a fixed residue modulo `2^i`, every row in one bit-`i` half must match every column in the other half (the difference has valuation `i`), and vice versa. So labels are constant on each half and swapped. Conversely, matching the low prefix and flipping bit `i` guarantees equal group labels. Even just the two shifts `+2^i` and `-2^i` suffice in a full finite packet: they make column labels periodic with step `2^(i+1)` by comparing the two shifts from an intermediate row, and every row has at least one valid shifted column. The virtual-label extension is a specified hypothesis for partial packets. Arbitrary holes without such an extension are not claimed to satisfy it.

The same integer Cauchy deformation from the public proof works at this geometry. For a tile of columns `y=y0+2^(i+1)j`, the source group is constant. If row `x` belongs to a different group, then `v2(t-y0+x)<=i`, since otherwise `y0-x` would have the matching low prefix and flipped bit `i`. At `T=t+2Lp*zeta`, the public geometric expansion is valid and has mode `L_k(X)=sum_j j^k X_y` with weight `k+1`. The same source mode may target several groups when their labels differ across prefixes, which is exactly the setting of (2). With `Mpad` padded slots in total, at most `(Mpad/Lp)*2^(i+1)=2Mpad/D_i` tiles exist, `D_i=Lp/2^i`. There is one mode per tile and weight, not per target group. Therefore the first-stage bound has `c<=2Mpad/D_i`. This statement is independent of how many labels occur in one packet at level `i`.

### Optimal padding for fixed labels

For consecutive retained levels `0,...,u-1`, each real edge has a sequence of ordered label pairs of length `u`. Fix those labels on the original vertices; we do not optimize the labels themselves. The minimum number of hierarchical packets of width `2^a`, `a>=u`, can be computed by a bottom-up recursion on the trie of these sequences.

At depth `u`, a subpacket has `2^(a-u)` slots. If the number of edges having that prefix is `n`, the minimum number of such subpackets is `ceil(n/2^(a-u))`.

At depth `i<u`, let `n_(x,y)` be the minimal number of child subpackets having next ordered pair `(x,y)`, computed recursively. For each distinct unordered pair `{x,y}`, the parent count required is `max(n_(x,y),n_(y,x))`; pad the deficient orientation by dummy child subpackets. For equal pair `(x,x)` the parent count is `ceil(n_(x,x)/2)`. Add these counts over the unordered pairs.

This recurrence is both necessary and sufficient. A parent subpacket has one unordered pair at this prefix, with exactly one child in each orientation (or two children with the same pair). Every child has the prescribed width. This gives the lower bound for each unordered pair independently; arbitrary pairing of the minimal child lists and adding dummy children achieves it. Induct from the leaves. The root count times `2^a` is the exact minimum padded slot count. The trie and counts have polynomial description length; explicit output costs time proportional to its size up to polynomial factors. A cap on padding can be enforced before materializing the output.

For flat packets each coordinate uses one unordered pair throughout. A synthetic width-`2^u` packet can have a distinct unordered pair at its final level for every one of its `2^(u-1)` low-bit prefixes. The hierarchical packing requires one packet, while the flat packing requires `2^(u-1)` packets. This demonstrates a real finite packing improvement for some label families. It is not an assertion that these synthetic patterns are dense elliptic Fourier boundaries. Also, the hierarchical algorithm guarantees localization but not low treewidth or pairwise mode locality; those are distinct issues.

## 3. Independence of localization labels on a product boundary

Consider a complete bipartite boundary between vertices `U` and `V`, with a uniform product edge distribution, or more generally multiplicities proportional to a product of endpoint weights. It can be a subset of a larger graph. Sample one real edge from this boundary. Let `A=(A_0,...,A_(u-1))` and `B=(B_0,...,B_(u-1))` be its row and column label sequences, with `r_i` possible labels at retained level `i`. Then `A,B` are independent. Distinct retained valuations may be sorted without assuming they are consecutive. We write the sequence order as `i=0,...,u-1` in the entropy argument.

Uniformly sample a virtual slot from the full padded packing. Its label sequences have distribution `Q`. The hierarchical condition implies that, conditional on any earlier label pairs, the law of `(A_i,B_i)` under `Q` is invariant under swapping its two entries: conditional on the packet and all lower position bits, bit `i` is uniform and swaps the pair, while earlier labels are determined. Averaging preserves this symmetry. The same condition holds for full cyclic or noncyclic packets that localize all shifts at each level, by Section 2.

### Full packing of a product boundary

If every slot belongs to the product boundary, its real distribution `P=Q` has the above symmetry. Under `P`, conditional on the earlier row prefix `a` and column prefix `b`, the pair law is a product `mu_a tensor nu_b`. Symmetry forces `mu_a=nu_b` (take marginals). Because the earlier prefixes are independent, every pair of prefixes of positive marginal probability occurs. Thus `mu_a` is constant in `a`, and `nu_b` is the same distribution independent of `b`. At every level, the label is independent of the earlier labels on its own side and has the same marginal law on the two sides.

Consequently the group labels at the different levels are jointly independent on each side. If `P(A_i=g)<=C_b/r_i` for every group, then

```
sum_i log2 r_i <= log2 |U| + u log2 C_b.                              (3)
```

Indeed `H(A_i)>=log2 r_i-log2 C_b` and `sum H(A_i)=H(A)<=log2|U|`. This is not an assumption of randomly chosen independent label words: it is a consequence of exact full-packet localization on a product boundary. In particular, hierarchical label correlations do not provide additional independently balanced groups at many levels for free.

### Robust padded version

Let the total number of virtual slots be at most `C` times the number of real edges in the specified product boundary, `C>=1`. The real label law `P` is dominated pointwise by `C Q` (sample the real boundary slots as a subset of virtual slots). Therefore `D_2(P||Q)<=log2 C`, with divergence in bits. This also covers padding containing other actual edges. Group balance is required on the real boundary, not on the padded population; a boundary occupying a fixed positive fraction of a balanced graph inherits such a bound up to a constant.

Let `D_i` in this paragraph denote the i-th conditional KL divergence term in the chain rule (not the packet width `D_t`). Its sum is at most `log2 C`. Under `P`, earlier prefixes `a,b` are independent and the next pair law is `mu_a tensor nu_b`. Under `Q` its conditional two marginals equal some `q_ab`. Pinsker's inequality and marginalization imply

```
TV(mu_a,nu_b) <= 2 TV(mu_a tensor nu_b,Q_i(.|a,b)),
E_(a,b) TV(mu_a,nu_b)^2 <= 2 ln(2) D_i.
```

Write `mu_bar=E_a mu_a` and `eps_i=E_a TV(mu_a,mu_bar)`. By convexity, the triangle inequality through an independent `nu_b`, and Cauchy--Schwarz,

```
eps_i <= E_(a,a') TV(mu_a,mu_a') <= 2 E_(a,b) TV(mu_a,nu_b),
sum_i eps_i^2 <= 8 ln(2) log2 C.
```

For distributions on `r_i` symbols, a maximal coupling with mismatch probability `d` gives the entropy difference bound `d log2(r_i-1)+h_2(d)` (condition on mismatch and the other variable). Thus

```
I(A_i;A_<i) <= eps_i log2 r_i + 1.
```

Set `R=sum_i log2 r_i` and `rmax=max r_i`. The chain rule for entropy gives

```
R <= log2|U| + u(1+log2 C_b)
       + sqrt(8 ln(2) log2 C * R * log2 rmax).
```

Using `sqrt(XR)<=(X+R)/2` yields the explicit bound

```
R <= 2[log2|U|+u(1+log2 C_b)] + 8 ln(2) log2 C log2 rmax.                (4)
```

This robust statement is exact for the described virtual-label packing model. In particular, when the packet sizes and label counts are polynomial in the original wire count `M`, the product boundary contains a fixed positive fraction of the real edges, `C_b=O(1)`, and `C=O(1)`, it gives `R=O(log M)`. The current two-stage Fourier compiler with merged source owners and arithmetic middle owners has full arrays joining these owner sets; the union of its full section-coefficient and quotient arrays supplies such a constant-mass product boundary (with endpoint-dependent word multiplicities if necessary). A changed compiler that splits these vertices or avoids that product boundary is outside the conclusion. We do not claim (4) for arbitrary incomplete packets without a compatible virtual-label extension, partial localization with large unmodeled exceptions, or a different deformation geometry.

A useful exact example: if the same balanced `r`-way partition is repeated at every one of `u` levels on both sides of a complete product boundary, with edge multiplicities divisible by `2^u`, Section 2 gives minimum padding factor

```
C = 1/r + (1-1/r) 2^(u-1).
```

Equal-label edges need no padding; an unequal-label sequence can only be paired with its reversed sequence at the first level, and each subsequent identical-orientation label forces a missing branch. Induction in the exact recurrence proves the formula. Repetition does not evade the entropy budget.

## 4. Nonuniform readout precision and the two-stage resource budget

Return to the consecutive integer grid `S={-Lp,...,Lp-1}`, `Lp=2^a`, local scale `alpha=2Lp`, and deformation degree bound `M(2Lp-1)`. For every retained valuation `i<a` let `p_i` be the maximum precision (in bits) requested at any node with that valuation. Unobserved levels have `p_i=0`. Coefficient orders and nodes may be chosen arbitrarily. Write `D_i=Lp/2^i` again, now for the packet width quotient. We prove a necessary *aggregate* precision bound, allowing nonuniform groups and node budgets.

**Theorem.** On the class of rational polynomials of degree at most `M(2Lp-1)` that are integral in every local expansion `f(t+alpha z)`, any readout from nonzero nodes of valuation `<a` which always recovers `f(0) mod 2` must satisfy

```
sum_(i:p_i>0) ceil(2 p_i/D_i) > M.                                  (5)
```

For adaptive readout, the condition applies to the precisions on the all-zero observation path. In particular `sum_i p_i/D_i > (M-u)/2` if at most `u` levels are observed. This holds for nonlinear readout from any specified jet coefficients, not merely a linear Hermite formula.

**Proof.** If the inequality fails, put `h_i=ceil(2p_i/D_i)` and `q_i=sum_(j>=i) h_j`, with `q_a=0`. At every nonzero grid point `s` of valuation `i` use root order `q_i`, and define

```
f(T)=product_(s in S\{0}) (1-T/s)^q_(v2(s)).
```

Then `f(0)=1`, and its degree is at most `M(2Lp-1)` because every root order is at most `sum h_i<=M`. At a nonzero node `t` of valuation `i`, the first nonzero local coefficient has valuation

```
V_t=(a+1-i) q_i + D_i sum_(j=1)^(a-i) (j/2^j)(q_i-q_(i+j)).            (6)
```

This follows directly by grouping the valuations of `(s-t)/s` at the full residue grid. Levels below `i` cancel; for each higher valuation `i+j` there are `D_i/2^j` points with that valuation and the same number of other points at valuation `i` whose difference from `t` has valuation `i+j`. At the root `s=t`, the leading scaled coefficient has valuation `a+1-i`. Every higher local coefficient has at least the same valuation, since `alpha/(s-t)` is integral when `s!=t`. At `t=0`, every coefficient is integral and the constant is one. All `q_i` are nonincreasing, so (6) is nonnegative at every center. For an observed level `i<a`, its first summand in the sum over j alone gives `V_t>=D_i h_i/2>=p_i`. Thus every requested coefficient is zero modulo its precision, just as for the zero polynomial, while `f(0)` is odd. This proves (5). The proof applies to the generic locally integral polynomial class; it is not a realization theorem for this polynomial as the deformation of fixed elliptic tensors.

### Matching upper taper up to constants

For a contiguous retained prefix `i=0,...,u-1`, take positive nonincreasing weights `ell_i` and let `R_ell=sum ell_i`. Set

```
h_i=ceil(2 M ell_i/R_ell),  q_i=sum_(j>=i) h_j.
```

The total number of Hermite coefficients is `sum_i D_i q_i>=Lp sum h_i>=2MLp`, enough for the degree bound. The public Hermite prefactor at a node of valuation `i` loses at most

```
D_i sum_(j>=1) j(q_i-q_(i+j))/2^j <= 6 D_i h_i,
```

since the slopes are nonincreasing and `sum j^2/2^j=6`. Its local evaluation at `-t/alpha` loses at most `(a+1-i)(q_i-1)` more powers of two. If `D_i>=u(a+1)` and `D_i>=R_ell`, and input precision is `O(M)`, these latter terms are `O(D_i h_i)` too. Thus precision `p_i=O(M D_i ell_i/R_ell + D_i)` is sufficient; with `R_ell=O(M)` the ceiling is absorbed. These conditions hold in the public large-packet schedule for the polylogarithmic weights in question. This is an actual nonuniform taper construction, not only a lower estimate.

### Consequence for the ordinary two-stage accounting

Take `ell_i=1+log2 r_i`. Consider the optimistic standard estimates for a second-stage join on `r_i` group vertices: mode enumeration/indexing has size scale `sqrt(Mpad p_i/D_i)`, and the degree-aware generic join buys a divisor `sqrt(ell_i)`. Ignore fanout, child entries, and all other costs, so the proposed log-time expression at node i is

```
C_i = sqrt(Mpad p_i/(D_i ell_i)).
```

Equation (5) implies

```
max_i C_i >= sqrt(Mpad (M-u)/(2 sum_i ell_i)).                        (7)
```

Indeed multiply the squared maximum by `sum ell_i`, and use the aggregate precision bound. The upper taper gives the matching scale when the `ell_i` are nonincreasing, under the specified regularity, and when the second-stage assumptions actually hold. Thus the total group-label resource `sum ell_i`, rather than an individual choice of a uniform r and u, is the relevant quantity for these estimates.

On the dense product boundary of Section 3, (4) bounds this resource by `O(log M)` at constant padding. Even allowing polynomial padding factor `C`, constant balanced-load parameters and polynomial packet sizes, it is at most `O(log M (1+log C))`. If the estimate charges padded slots `Mpad` (as the present shared-mode catalogue bound does), the factor `C/(1+log C)` is bounded below by a positive constant. Therefore (7) cannot certify a log-time scale smaller than `M/sqrt(log M)` by these changes alone. This is a barrier to the particular **precision and generic-join accounting**, not a running-time lower bound for actual tensor contractions; e.g. additional cancellations, tighter actual mode catalogues, better readouts on realizable tensor families, or a different join theorem could evade it. This conclusion is more robust than simply comparing the two uniform upper-bound formulas in earlier checkpoints.

For transfer to ECDLP, a genuine contraction bound `log T=O(M/(log M)^(2/3))` at `M=O(mb)`, plus the previous verified relation reduction and preprocessing, would balance at arity `m=Theta((log b)^(1/3))` and yield `log E T=O(b/(log b)^(1/3))`. The present estimates on the dense boundary still optimize at the old 1/4 power. No new coefficient or faster ECDLP algorithm is claimed.

## 5. Matrix-shift packet geometry (separate diagnostic)

A different packing supplies many group labels at a single interpolation level without independent coordinate label arrays. Let `U=V=Mat_(k x ell)(F_2)` with `ell>=k`. Retain edges of full-row-rank difference `Z=A+B`. For fixed such Z, the map `X -> XZ` on `Mat_(k x k)(F_2)` is injective. Partition row matrices into affine cosets of its image. One packet contains `q=2^(k*k)` pairs indexed by `X`:

```
row(X)=A0+XZ,    col(X)=A0+XZ+Z.
```

This partitions every retained edge exactly once. The discarded fraction is at most `(2^k-1)/2^ell`, by union over nonzero row vectors annihilating Z. It is not free to discard those edges in a full contraction.

If `I+T` is singular, choose `v!=0` with `v(I+T)=0`. At the XOR shift `X -> X+T`,

```
v row(X) = v col(X+T).
```

Thus grouping all vertices by `v A` into `2^ell` groups localizes every packet at this shift, including arbitrary affine cosets. The row and column marginals on the retained graph are uniform. This is an exact algebraic construction of many balanced groups at one shift.

It does not supply the useful crossing-mode spectrum. Identify the additive matrix index space with `F_q` and use unit-separated lifted interpolation knots for XOR shifts. The first crossing coefficient after division by the dyadic uniformizer at t is the Cauchy value `1/(t+x+y)` in the residue field. Inside a full packet, choose a nonzero linear functional of `vX` (equivalently of the group label restricted to that packet). Splitting the row and column groups by that bit, the off-diagonal block between halves has size `q/2`. Its row and column positions are cosets of the same additive hyperplane H; the denominators are `eta+h+h'` with `eta` outside H, by localization at t. This matrix is invertible: it is convolution on H by `c(h)=1/(eta+h)`, whose augmentation is `P_H'(eta)/P_H(eta)!=0`, with `P_H(X)=product_(h in H)(X-h)` an additive polynomial of nonzero constant derivative. The additivity and derivative follow by induction, adjoining a coset H+a and using the characteristic-two identity for the product. In the group algebra of the elementary abelian group H, the augmentation ideal is nilpotent (write the basis generators as `1+e_i` with `e_i^2=0`), so nonzero augmentation implies invertibility. This also recovers the previous unit-separated XOR rank obstruction for this more structured packet geometry.

There are therefore at least `q/2` first-order crossing directions per full retained packet in a normalized linear mode expansion. A generic interpolation readout on degree `M(q-1)` polynomials needs more than M bits at some nonzero knot: use the polynomial `product_(t!=0)(1-T/t)^M` with lifted nonzero residue knots, constant one, whose local jets at those knots are divisible by `2^M`. Thus a full occupied version of this geometry admits exponentially many feasible weight-one exponent subsets at the necessary precision. It does not yield the wanted subexponential mode catalogue through unrestricted enumeration. This limitation does not rule out different kernels or an implicit joint calculation.

## 6. Reproducibility and scope

Run `python -B packet_audit.py` and `python -B shared_mode_audit.py` in this directory (Python standard library only). The former includes:

- 197,136 exhaustive locality/prefix-swap comparisons at packet widths 2, 4, 8;
- 15,625 exact product-label-law comparisons, with 49 symmetric cases; these verify symmetry if and only if label independence and equal marginals in the tested laws;
- 250 independent minimum-packet comparisons against exhaustive covering by all 21 width-four hierarchical patterns;
- independent-label no-padding and repeated-label exact-padding examples, and hierarchical versus flat synthetic packet counts up to width 256;
- 150 padded entropy regressions, checking the real/synthetic distribution domination, exact conditional swap counts and numerical entropy inequalities;
- 40,800 exact valuation comparisons, 336 explicitly expanded rational local polynomials (10,016 coefficients), and 17,168 nonuniform upper-taper checks;
- 276 full-rank matrix differences, 3,660 coset packets, 41,316 localization identities, and 11 finite-field crossing-rank calculations.

The shared-mode verifier uses finite unramified rings modulo `2^P`, not rational approximations to dyadic arithmetic. It compares 40 random multi-target mode joins and 12 actual width-four hierarchical Cauchy packet joins against direct Boolean enumeration. It covers 3,939 Newton summands, target fanout three, and 32 nonzero output controls. Both transcripts report all assertions passed. These are finite component tests and do not benchmark a complete contractor, construct an unrestricted elliptic network, or prove a practical crossover.

**Next precise gate.** The shared-mode identity makes arbitrary fanout tolerable in the first-level sum, but the second-level compressed graph still needs a cheap *joint* contraction of source Newton factors and target evaluations. The dense product-boundary entropy theorem also shows why choosing correlated localization labels on the existing two-stage Fourier layout is insufficient. To beat the exponent, construct a relation compiler without that constant-mass product boundary at the same leading wire order, or a joint join/readout whose cost falls below the standard weighted-mode and generic-join expression even on that boundary. Actual full relation transfer remains open.
