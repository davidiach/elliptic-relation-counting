# Which directions disappear at a displacement?

This note studies the cross-multiplied pin means of [§2](RESEARCH.md), including the directions missed by a fixed list of displacements. Its statements concern the factor-base columns before adjoining the explicit consistency rows. The [single-anchor translation equations](ANCHORED_TRANSLATIONS.md) use different rows with additional known right sides. The per-displacement theorem below applies to the asymptotic punctured seed. The numerical examples are finite rank calculations.

## Two reasons a direction can disappear

Fix a displacement t with a nonempty set of ordered tuple pairs. Work over the group-order field. Let A have one signed incidence row per pair, and let H have one column for each source or target pin, with a one wherever that pin occurs. Thus n = Hᵀ1 contains the pin masses and HᵀA contains the pin incidence sums. Each column of A is a difference of two columns of H. Summing the indicators of all points in any one slot gives 1.

Write D for the matrix of raw contrasts: subtract one reference row of A from every row. Its kernel consists of the directions u for which Au is constant. Write R for the complete matrix of cross-mean rows at t. For any two pins h,k, their row is the sum of aᵢ − aⱼ over i satisfying h and j satisfying k. Consequently row(R) is contained in row(D), even when all pin masses vanish modulo the group order.

**Exact rank-loss identity.** Suppose at least one component of n is nonzero in the field. Put V = col(H) and L = col(A) + span(1). Then

    rank(D) − rank(R) = dim(L ∩ V⊥).

Indeed Ru = 0 exactly when there is a unique scalar c with HᵀAu = c n: solve for c using any nonzero pin mass. The linear map u ↦ Au − c1 sends ker(R) onto L ∩ V⊥. For surjectivity, write a vector of that intersection as Au − c1; orthogonality gives HᵀAu = c n. The kernel of the map is exactly ker(D). Rank-nullity gives the identity. Since L is contained in V, the intersection lies in the radical V ∩ V⊥. Over the reals that radical is zero; over a finite field its intersection with L can be nonzero.

This gives two distinct proof obligations. A short list can leave a nonaffine direction with constant raw increments at every selected displacement. Every pin comparison then misses it before modular cancellation enters. Even where the raw contrasts have the required rank, modular averaging can lose some of it.

**Reversing the displacement.** The bijection (X,Y) ↦ (Y,X) takes t to −t, exchanges the source and target pin names, preserves the corresponding masses, and negates every incidence. After the pin relabeling, the cross-mean rows change sign. Thus R(t) and R(−t) have the same row space. The same holds for D.

**A sufficient support condition for raw contrasts.** Suppose the short tuple sums cover the entire prime-order group, and t is nonzero. Then the kernel of D(t) is precisely the affine tables of §2. To see this, suppose the score increment on every pair is c. Two source tuples of value g share a successor of value g+t, so their scores agree. The induced function f on group values satisfies f(g+t) − f(g) = c everywhere. Since t generates the group, f is an additive function plus a constant. Varying one slot with the other slots held fixed gives one additive slope across all blocks, plus one constant per block. The converse follows by substitution. This specializes the common-successor argument in [the earlier fixed-source proof](../../NEIGHBORING_FIBERS.md). The seed energy estimate alone allows holes in the short-tuple sumset.

## Full raw rank at a prescribed displacement from punctured support

**Theorem.** Let the nonempty blocks F₁,…,Fₘ lie in a cyclic group G of prime order r. For each i let D₋ᵢ be the set of sums from choosing one element of each other block. Suppose |D₋ᵢ| ≥ (1−ε)r for every i, with 6mε < 1. For every prescribed displacement t, including t = 0, the raw contrasts D(t) have kernel exactly

    uᵢ(A) = λ(A) + cᵢ,

where λ is an additive map from G to Fᵣ and the cᵢ are arbitrary block constants. For at least one nonsingleton block this gives rank N₀−m−1. The theorem permits holes in the short-tuple sumset.

**Proof.** Denote the full sumset by S. It contains a translate of each D₋ᵢ, so its density exceeds one half. Fix t. There are pairs at displacement t by intersection of S and S−t. If u annihilates D(t), then every pair of tuples at displacement t has the same score increment. At every g in S with g+t in S, two source tuples at g have equal scores by comparison with the same successor. Let f(g) be this score and extend f arbitrarily to all other g.

Fix A,B in Fᵢ. For every x such that x and x+t both belong to D₋ᵢ, choose a remainder tuple of sum x. Adjoining A or B to the same remainder gives two scores whose difference is uᵢ(B)−uᵢ(A). Both f(x+A) and f(x+B) are defined by actual source scores because x+t also has a remainder representation. Hence

    f(g+B−A) − f(g) = uᵢ(B) − uᵢ(A)

on a set of starting points g of density at least 1−2ε. This includes the case t = 0.

Every group element s is a difference of elements of S since S has density greater than one half. Write that difference as a sum of m within-block differences. By intersecting the translated good sets of those steps, f(g+s)−f(g) agrees with the sum of the predicted increments on a set of density at least 1−2mε. Any two such paths for s agree on some starting point, so define λ(s) by their common increment. To compare λ(s), λ(v), and λ(s+v), intersect the domains of the direct path, the path for s, and the path for v translated to begin at g+s. Their intersection has density at least 1−6mε > 0. Telescoping gives λ(s+v) = λ(s)+λ(v). A single within-block difference is an m-step representation padded by zero differences; hence uᵢ(B)−uᵢ(A) = λ(B−A). This gives the claimed form. Conversely an affine table has constant increment λ(t) on all pairs, so it annihilates every raw contrast. ∎

**Consequence for the current seed.** For any fixed positive deletion fraction in [the deletion-resilient seed theorem, §1](RESILIENT_SEEDS_AND_VARIANCE.md), single-block deletions are eventually permitted and ε = 2^(−Ω(b/m)) satisfies 6mε < 1. Conditional on that seed theorem's construction, the raw rank is therefore maximal for **every** prescribed t. This removes the shared-raw-kernel obstruction throughout that regime. The proof counts no tuples algorithmically and gives no running-time improvement by itself.

At a shift where n ≠ 0, the exact remaining cross-mean deficiency is dim(L ∩ V⊥) by the identity above. When all ordered-pair pin masses are zero modulo r, all cross-mean rows are zero. The seed's no-wrap guarantee concerns single-tuple fibers and cannot supply this ordered-pair mass condition. Neither the finite examples nor the theorem bound the common kernel after averaging at an affordable list of shifts.

## Explicit witnesses on the order-19 curve

Use the curve y² = x³ + 2 over F₁₃ and P = (1,4). Coordinates of a block below are indices in the cycle generated by P, in the displayed order. Each vector u is written in concatenated block order. The increments are integer row scores modulo 19. For each witness, the complete list of invisible displacements is zero and the displayed shifts and their negatives. All other displacements detect that witness.

| Blocks | Shifts | u | Raw increments, in shift order | Individual ranks | Stacked / complete rank |
|---|---|---|---|---|---|
| (2,7,9,14,15), (4,6,9,11,18) | 7 | (0,0,0,0,1; 0,0,0,0,0) | 0 | 6 | 6 / 7 |
| (3,6,9,14,17), (5,8,11,12,16) | 3,8 | (0,0,0,0,0; 0,0,0,1,0) | 0,0 | 6,6 | 6 / 7 |
| (1,2,17,18), (4,6,13,15) | 1,2 | (0,0,0,0; 1,1,0,0) | 0,0 | 4,3 | 4 / 5 |
| (2,4,7,14), (12,14,15) | 2,3,9 | (1,0,18,0; 1,0,0) | −1,−1,0 | 2,2,3 | 3 / 4 |

The first two rows account for all six deficient nonidentity displacements in the original 68 accepted fixtures. The original eight sampled direction checks per fixture missed these extra kernel vectors. In the last two rows every integer one-pin count, for every displacement, lies between 1 and 6. Thus even full support of every pin with no wrap permits the common raw obstructions shown here. The short-tuple sumsets have sizes 14 and 10 respectively.

Both final examples come directly from coordinate masks. Keep a point when (ax+b) mod 13 < 3, excluding the identity. The third example uses (1,0) and (2,6) in both parities. The fourth uses the even/odd coefficient pairs ((7,3),(5,3)) for the first block and ((9,2),(6,11)) for the second. Its three shifts admit no choice of signs with signed sum zero modulo 19.

## Exhaustive finite family

On the same curve, take every nonempty distinct block cut out by the 156 maps with a ≠ 0, using the same map in both parities and threshold 3. There are 60 such blocks. Of the 1,830 unordered pairs with replacement, 996 meet the original fixture filters: the short sumset occupies over half the group, and 12 times the maximum number of modular zeros of an individual pin is strictly less than 19. The complete stack has rank N₀−3 on every accepted pair. There are nine nonzero displacement classes up to sign.

| Fixed displacement list | Lists failing on at least one of the 996 fixtures | Lists attaining full affine-quotient rank on all 996 |
|---|---:|---:|
| Two distinct classes | 36 / 36 | 0 |
| Three distinct classes | 15 / 84 | 69 |
| Four distinct classes | 0 / 126 | 126 |

Across all 996 × 9 tests, 636 individual classes were deficient. **All 636 had exactly the same deficiency before pin averaging.** The equality of ranks gives equality of the corresponding row spaces by the containment above. Of the 15 failing triples, each admits a signed sum of zero. The parity-split example in the preceding table shows this pattern depends on the mask family.

For contrast, in the abstract cyclic group of order seven take blocks (1,3), (2,4,5), (2,4,5) at displacement 1. The tuple values cover the group. The raw rank is 4 and the cross-mean rank is 3. For u = (6,0; 3,6,0; 1,0,0), all seven raw scores occur. Yet Hᵀ(Au−4·1) = 0 and at least one pin mass is nonzero modulo 7. This is a nonzero vector in L ∩ V⊥. Some pin counts wrap in this example.

## What happens when the seed potential selects the blocks?

Use the same prime field, group, two slots, and threshold 3. Allow the even and odd affine maps to vary independently in each slot. This yields 3,600 distinct blocks and 6,481,800 unordered block pairs. Write W for the tuple count and C for the sum of squared tuple multiplicities. At this threshold the prescribed mean is μ = (18·3/13)² = 2916/169. Removing the constant and clearing the positive denominator in the potential gives the objective 19·169·C − 5832·W.

Two finite selection rules give different outputs:

| Selection rule | Scaled objective | Occupied tuple values | Rank at nonzero shifts |
|---|---:|---:|---|
| All 144 unordered block pairs minimizing the objective globally | −45,998 | 19 / 19 in every case | Target attained at all 18 shifts for every minimizer |
| Bitwise conditional expectation, with the encoding below | −26,471 | 13 / 19 | 5 / 5 at all nine nonzero sign classes |

For the exhaustive calculation, let k_B(s) count pairs of elements of B with difference s. The energy of the sum of two independently chosen block elements is exactly |B₁||B₂| + 2 ∑ₛ k_B₁(s) k_B₂(s), summing one representative of each nonzero sign class. Grouping by the 700 distinct profiles (|B|, k_B(1), …, k_B(9)) covers every block pair exactly, with no rank inference from profiles: reconstruct all minimizing blocks and test each rank separately. The direct tuple counts also reproduce the minimum objective.

For the bitwise experiment, keep all 156² equally likely map choices per slot, including repeated realizations of the same block. Encode (even multiplier, even offset, odd multiplier, odd offset) as four unsigned four-bit digits, most significant bit first, fixing slot one before slot two; choose zero on ties between nonempty prefixes. The profile identity is bilinear, so the exact conditional expectation uses the profile sums of both remaining seed-completion lists divided by the product of their sizes. This gives block indices (3,4,17) and (4,5,9,13,16). Full tuple support suffices for raw rank but is unnecessary on this instance. The threshold and bit order are explicit choices for this finite model; the asymptotic parameters are outside its range.

## Scope and reproduction

These masks implement the finite coordinate rule. The large-input seed construction has additional quantitative energy requirements, and these examples do not establish a failure on its selected outputs. Four-class success is limited to the finite symmetric family. The complete cross-mean rank theorem and the general deterministic selection problem are unchanged. For the deletion-resilient seed the raw part now follows at every shift. A uniform result still needs modular pin masses and simultaneous rank retention for the averaged equations.

Run `python -B check_displacement_kernels.py` in this directory; see [the transcript](displacement_kernel_verification.txt). The verifier computes ranks over the group-order field and independently enumerates raw tuple incidences for every deficient class and each named witness. It prunes a longer list only when a subset already has full rank. The original fixture totals can be reproduced using `python -B check_short_relations.py`.

Run `python -B check_seed_potential.py` for the exhaustive minimizer and deterministic-prefix comparisons. The [seed transcript](seed_potential_verification.txt) records their results.


[Simultaneous pin loss](SIMULTANEOUS_PIN_LOSS.md) gives a sufficient connectedness condition for one averaged displacement, classifies the extra kernels for blocks missing one group element, and exhibits a shared nonaffine direction invisible at every shift in a family outside the resilient seed regime.
