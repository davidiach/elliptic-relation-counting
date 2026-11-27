# Common translation covers at the quarter-power scale

27 November 2026. Research follow-up to PR #18 (`9d3ae7ab3e9253391c392bf0a4a06e1e935a2631`) and merged PR #13 (`bef28551c5fdd352cc3e641cd25e7b3f6ec08eb8`). Baseline main checked at the PR #15 merge `18ce38e4b108a73497f0fa89be4883578255fead`. This note improves the **conditional reduction** to additive range avoidance, not the avoidance algorithm or the end-to-end deterministic ECDLP bound. All resource claims use LZR ePrint 2026/3020 Theorem 2.2 (closed exact Boolean contraction), its ordinary arithmetic compilation and point-counting setup, as in `MAIN_RESULT.md` and `COMPLEXITY_AUDIT.md`. The energy and covering identities are elementary and independent of that computational premise.

## 1. Why pointwise-flat short fibers are unnecessary for PR #18

The two-dimensional reduction in PR #18, section 2, needs universal tuple coverage, an affordable ordinary decomposition routine, and a tuple product smaller than r^2. It does not need enumeration of complete fibers, near-uniform multiplicity at each individual target, or ordered-pair pin counts below r. Its invocation of the batched high-degree masks therefore imposed a stronger seed property than necessary. Below we retain only global mean-square flatness and repair any uncovered group values by adding a small explicitly chosen block of translates.

## 2. Exact affine-mask autocorrelation identity

Let G be a finite abelian group of order r, with identity O. In each of m indexed blocks retain nonidentity elements with pairwise-independent indicator variables of common marginal theta; the blocks are independent. Put q=(r-1)theta. For a nonempty retained index set I of k blocks define d_I(g) as the number of one-per-block tuples with value g, W_I=prod_{i in I}|F_i| and C_I=sum_g d_I(g)^2. Then

    E[r C_I - W_I^2] = (r-1) [q-(r-2)theta^2]^k.             (1)

Proof. For the block indicator f_i, write H_i(t)=sum_x f_i(x) f_i(x+t). Its expected value at zero is q. At every nonzero t exactly r-2 ordered pairs (x,x+t) have both entries different from O. Pairwise independence gives expected value c=(r-2)theta^2 there. Thus E H_i=c 1_G + a delta_O with a=q-c in [0,q]. The autocorrelation of d_I is the convolution of the H_i. Independence of the blocks allows expectation to pass into each factor. Every expanded convolution containing a constant function is constant; the sole other term is a^k delta_O. Its value at zero is E C_I, and its total mass is E W_I^2. Subtracting r times the constant value from the total, in the order r C_I - W_I^2, leaves (r-1)a^k, as asserted. QED.

This also proves the bound E[(r C_I-W_I^2)/(r q^k)]<=1 without any division by a random block size. The numerator is always nonnegative by Cauchy--Schwarz. For the elliptic-curve application use independent field-affine masks ax+c for each of the two y parity classes in each block, a,c uniform in F_p, **including zero slopes**, and retain if the canonical mask value is below a threshold T_0. Points in one parity class have distinct x coordinates, so the affine masks are pairwise independent; opposite parity classes have independent seeds. This is the law of PR #13, not a claim that nonzero-slope masks have the same exact marginal. A size penalty below controls the constant-mask cases.

## 3. A deterministic exponentially small mean-square error, jointly over deletions

Let b=ceil(log2 p), L=log2 b, and assume the curve is nonsingular, p>3, and its whole rational group has prime order r, so log2 r=b+O(1). Let m grow polylogarithmically in b. Fix rational constants 0<=alpha<1/2 and 0<sigma<1 with

    kappa = (1-alpha)(1+sigma)-1 > 0.

Let a_0=m-floor(alpha m), and let J consist of all retained index subsets of size >=a_0. Write J0=|J|<=2^m. Set q_* = ceil(r^((1+sigma)/m)), T_0=ceil(p q_* /(r-1)), q=(r-1)T_0/p. This uses exact integer powering and binary search for rational sigma, not real exponentiation in the algorithm. Eventually q>=q_*, T_0<=p, and q^m=(1+o(1))r^(1+sigma). The initial bounded range of b is separate from the asymptotic construction.

Choose the affine mask seed by bitwise conditional expectation for the nonnegative rational potential

    Phi = sum_{I in J} (r C_I-W_I^2)/(r q^|I|)
          + J0 m^4 sum_i ((|F_i|-q)/q)^2.                   (2)

Equation (1) bounds the expectation of each energy term by one; pairwise independence gives Var(|F_i|)=q(1-theta). Consequently E Phi <= J0(1+m^5/q)<=2J0 eventually. Conditional expectation gives a deterministic seed with Phi<=2J0. Its size term ensures |n_i/q-1|<=sqrt(2)/m^2 for all n_i=|F_i|. In particular all blocks are nonempty and, uniformly over I, q^|I|/2<=W_I<=2 q^|I| for large m. Hence for the normalized tuple sum law

    chi_I^2 := (r C_I-W_I^2)/W_I^2
              <= 8 J0 r/q^|I| <= 8*2^m r^(-kappa).        (3)

The last bound uses |I|>=a_0>=(1-alpha)m. Its right side is r^(-kappa+o(1)). Every term in (2) is obtainable by exact ordinary counts: C_I uses two tuples with equal sums, W_I^2 two unrestricted tuples, and size moments one and two point memberships. When the affine coefficients are partially fixed, include the remaining coefficients as canonical field variables with the prefix restrictions. Divide the count by the number of allowed seed completions. The latter is a product of sizes of word intervals (or can be counted separately). Signs and rational scalings are external to the networks. Each network contains O(mb) Boolean edges with the inherited degree O(b^(3/4)) and exact count precision O(mb); unused coefficient variables can be omitted with their completion factors accounted for. There are O(2^m+m) count types per prefix and O(mb) bit decisions. No private arbitrary label lookup is used. This is the same counting interface as the empirical-mass second moment in PR #13, but its target is deliberately much more redundant.

## 4. A common small translating block, selected without reading the support

The next lemma is independent of the affine seed. Suppose finitely many nonempty tuple laws d_I have positive masses W_I. Fix P in G minus {O}. For integer s>=1 take a_1=P fixed, and a_2,...,a_s independent uniform elements of G minus {O}. Put

    f_I(g) = (r d_I(g)/W_I - 1)^2,
    Psi(a_1,...,a_s) = sum_{I in J} sum_{g in G} prod_{j=1}^s f_I(g-a_j).

Then

    E Psi <= sum_{I in J} r (r/(r-1))^(s-1) (chi_I^2)^s.  (4)

Indeed f_I>=0 and sum_g f_I(g)=r chi_I^2. For each fixed g a uniform nonidentity translate has E f_I(g-a_j)<=r chi_I^2/(r-1). Independence bounds the product of the s-1 random factors, and summing the fixed first factor proves (4). If Psi<1, for each I and every g at least one of the numbers f_I(g-a_j) is <1: otherwise that term alone contributes >=1. For the corresponding j,

    0 < d_I(g-a_j) < 2 W_I/r.

In particular a single common set T={a_1,...,a_s}, duplicates removed, satisfies

    T + sum_{i in I} F_i = G                  for every I in J.    (5)

There is no inference that the translated convolution is pointwise uniform or that its fibers are affordable to enumerate. The group-wide covering property, rather than complete fiber enumeration, is what PR #18 requires.

With (3) choose a fixed integer s>1/kappa. Bound (4) is at most r^(1-s kappa+o(1))<1 eventually, including the sum over J. When alpha=0 and sigma=2/3, s=2 suffices: just one new group point besides P is needed. For any fixed sigma>0 and alpha=0 take s>1/sigma. No requirement is imposed that the selected translates have known logarithms. Treat the small set T as one additional indexed block (P remains an available reference). This detail is essential: requiring a variable scalar multiplication by P inside the conditional expectation network would not be covered by the O(mb)-edge claim.

### Exact prefix implementation

For each I expand

    W_I^(2s) prod_j f_I(g-a_j)
       = prod_j (r d_I(g-a_j)-W_I)^2.

There are at most 3^s monomials. A term with exponents ell_j in {0,1,2} is evaluated by an ordinary nonnegative count with ell_j independent tuples of type I constrained to sum g-a_j. Include a common canonical group point g and the unselected translating points (on-curve nonidentity variables with the chosen prefixes) even in terms with exponent zero; otherwise the completion multiplicity would be wrong. The identity is allowed for g, with a canonical encoding. Apply the signed integer coefficients, divide by W_I^(2s), and normalize by the number of translating-point completions. The latter can be computed by ordinary point counts with their prefixes (the point domains are independent). The initial shift P is a fixed group coordinate, not an implicit unknown scalar. All computations are exact. At most O(sb) binary prefix choices fix all remaining translates, ignoring empty branches and choosing a nonempty child of no greater conditional expectation. The final potential is below one by (4).

There are O(J0 3^s) count terms per expectation. Each has at most 2s tuples, so the network has O(smb) edges, O(b^(3/4)) maximum degree, and O(smb) exact count precision. Canonical encodings and equality fanout trees for the group target and shift variables are charged, with no full scalar multiplication gadget on a private b-bit scalar. For fixed alpha,sigma,s the resulting log bit-time envelope for seed, translations, explicit block enumeration, and ordinary decompositions is

    O(b/m + mb/sqrt(log b) + m + log b).

One may enumerate the constant-size T choices externally and decompose the corresponding ordinary target; adding the block does not require an arbitrary lookup tensor. Ordinary group arithmetic verifies each returned tuple. With m asymp L^(1/4) this is 2^{O(b/L^(1/4))} bit time and space under the closed contractor. The size of the full tuple family is at most

    |T| W_full <= 2s q^m = r^(1+sigma+o(1)) < r^2.         (6)

Thus the required universal short coverage and strict product promise hold simultaneously, with **constant-degree affine masks**, at the quarter-power scale. There is no growing-degree membership polynomial. This is stronger for the PR #18 reduction but weaker than the globally pointwise/enumerable seed for applications that actually require small fibers.

## 5. End-to-end *implication*, with the avoidance problem still open

Apply the verified quotient reduction of PR #18 section 2 to the m+1 lists (one translating list of size at most s). A nonzero calibrated row comes from two ordinary decompositions at P and O, with known right side one. Given any consistent subsystem with nonzero residual dimension, project the quotient by its known zero-value row subspace surjectively to F_r^2. Project each explicitly listed block and remove duplicates. By (6) the product of the projected list sizes is less than r^2. A point supplied outside their Minkowski sum lifts to a coordinate vector whose group value can be decomposed by (5). Subtracting that tuple's incidence gives an independently checkable zero-value relation outside the current row span. Each such step increases rank. At most D steps are required, where D is the explicit number of nonreference columns, and all polynomial overhead in D,b is absorbed in the exponent. This is the previous reduction; the new theorem removes the more expensive seed cost. It does not supply a missing-point finder.

Suppose a uniform deterministic AdditiveAvoid(2) routine were proved for explicit lists with product <r^2 in bit time N^{O(m^theta)} poly(b) for fixed 0<=theta<1, where N is the total explicit input length up to polynomial factors in b. The m+O(1) lists, and the O(1) extra translating labels, preserve this envelope. Including up to D calls, the full conditional time exponent would now be

    O(b/m + mb/sqrt(L) + m^theta b/m).

At m asymp L^(1/(2(2-theta))) this yields

    delta = (1-theta)/(2(2-theta)).

The previous PR #18 result, using the stronger global pointwise seed, gave any saving below (1-theta)/(2(3-theta)). Compare:

| Hypothetical deterministic avoidance time | Previous resulting saving (approached) | New resulting saving |
| --- | --- | --- |
| N^{O(1)} | 1/6 | 1/4 |
| N^{O(m^{1/4})} | 3/22 | 3/14 |
| N^{O(sqrt(m))} | 1/10 | 1/6 |
| N^{O(m^{3/4})} | 1/18 | 1/10 |

These are sufficient implications only, not attained ECDLP exponents. The desired algorithm remains missing. For any fixed product promise r^(1+eta) with eta>0, choose 0<sigma<min(eta,1), alpha=0 and constant s>1/sigma. Then (6) fits that stronger promise eventually, with the same asymptotic complexity exponent. The constants depend on eta. No assumption about the actual residual output is used in this implication.

### Nearly linear product without the high-degree masks

If one insists on r^(1+o(1)) instead of the slack in (6), take alpha=0, sigma=1/ceil(log2 L) for sufficiently large L, and s=ceil(3/sigma). The exact integer threshold construction and size estimates above remain valid, and the product potential has expectation r^(-2+o(1)); the constant-factor bounds contribute only exp(O(s))=r^o(1). Then 3^s=poly(L), and the shift networks have O(smb) edges. The envelope becomes O(b/m+smb/sqrt(L)). Balancing m asymp L^(1/4)/sqrt(log L) gives log T=O(b sqrt(log L)/L^(1/4)), with product r^(1+o(1)). Hence a polynomial-time avoidance algorithm on this smaller-product regime would suffice for every fixed saving below 1/4. This version also does not claim pointwise small fibers. For growing s, one can test all s translate choices during a decomposition; their factor is charged.

## 6. Deletion transport at the same budget

The common cover can tolerate any fixed fraction alpha<1/2 of the main blocks being removed: choose sigma in (alpha/(1-alpha),1), and s>1/kappa as above. This is full coverage of the retained group values, not the higher-moment no-wrap guarantee of RESILIENT_SEEDS_AND_VARIANCE.md. The full product is still below r^2. To obtain full coverage after deletion of **any one block including the translating block**, use two identical indexed copies of T and choose any fixed alpha>0 with alpha<1/2. At sufficiently large m, removal of one ordinary block is allowed. With both T blocks retained, one copy gives coverage and choosing any fixed point in the other preserves it. Removal of a T block leaves the other copy. The product bound becomes s^2 W_full<r^2 and the time stays at quarter-power. PR #18's explicitly constructed short-fiber transport relations then apply to this family with O(D(m+2)b) ordinary/pinned decompositions. Their integer lattice contains r Z^D; all Smith factors are 1 or r. No claim of only one r factor follows. This transfers that constructive torsion-exponent guarantee to the original quarter-power resource scale without the globally pointwise seed. The repeated indexed copies are legal separate blocks; deleting a T block is not silently covered by (5) alone.

## 7. Exact tests, adversarial scope and literature

Run `python -B check_translation_cover.py` (standard library). `verification.txt` records the results. There are 8,512 complete independent Bernoulli mask arrays in four tests of (1), 32 exact affine autocorrelation/energy comparisons on two small prime-order elliptic curves, and 38 parent/child conditional-expectation identities for the actual affine-seed selection on the order-seven curve. That small seed does not meet the large asymptotic size/energy guards; it verifies the exact calculation, not the asymptotic runtime. Five independent toy two-block fixtures at orders 19,31,101,257,1009 deliberately have a hole before translation and no holes afterward; their shifted tuple products are all less than r^2. The order-19 fixture checks actual elliptic group additions. Three further cases test a common T across all retained subsets; some full toy products exceed r^2, so these cases verify the covering lemma only. Across both selector types the transcript records 80 exact translation prefix identities and 43 exact direct-vs-expanded count checks. These are finite mathematical fixtures with explicit small group coordinates; they do not implement the asymptotic contractor or any faster DLP solver.

The elementary idea that a dense set can be covered by a few translates is classical. The contribution here is the exact nonnegative potential expressed in the ordinary contraction interface, the simultaneous deletion guard, and the resulting improved reduction budget. For background see Juanjo Rue's lecture handout, *Covering by translates of the set*, Lemma 3.4 (https://web.mat.upc.edu/juan.jose.rue/BlockCourse-FUBerlin-Week1/Berlin_Handout2.pdf). This is not a claim of novelty for additive covering in general. Related range-avoidance results for bounded-locality Boolean circuits (Guruswami, Lyu, Wang, APPROX/RANDOM 2022; Kuntewar and Sarma, APPROX/RANDOM 2025; and ECCC TR25-049) do not directly give the required sublinear-in-arity exponent for arbitrary explicitly listed modular summands in this constant-stretch regime. A circuit with an output longer than its input is not by itself an efficiently derandomizable circuit. The literature scan is preliminary, not an exhaustive priority review or a hardness claim for AdditiveAvoid(2).

The main unsuccessful lines in this pass were generic low-depth circuit avoidance, Fourier small-bias estimates for the sparse sumset, and conditioning on progressively finer modular buckets. None supplied a uniform missing-vector constructor at the required budget. In particular a small per-bucket count does not make the next bucket count cheap, and computing arbitrary quotient labels inside a low-degree contraction would silently reintroduce the lookup bottleneck. No claimed progress depends on those lines.

**Status.** The full-coverage setup and improved reduction complexity are proved conditional component results. The universal deterministic positive-saving ECDLP theorem is still open. An affordable deterministic AdditiveAvoid(2) algorithm is a load-bearing missing lemma; the upstream closed contraction premise and the new proof-chain/priority review remain unverified independently. The mandatory end-to-end checklist therefore does not pass. No production cryptographic instances are used.
