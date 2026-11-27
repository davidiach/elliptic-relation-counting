# When modular pin averaging loses the same directions at several shifts

This note concerns the two-sided cross-mean matrix R(t) on factor-base columns from [the original system, §2](RESEARCH.md). A shift always means a group displacement. All vector spaces below are over the prime field of the cyclic group's odd order r. The results concern matrix rank before adjoining other equations.

## A shared failure with every ordered-pair pin count below the field order

In C₁₀₁ use seven binary blocks, each ordered as (1,1+dᵢ), with steps (1,2,4,8,16,22,50). Of the 101 full-tuple fibers, 74 have size one and 27 have size two. **Every positive or negative one-pin count at every displacement lies between 64 and 92 as an integer.** The following direction, in the displayed column order, is lost at both 25 and 26:

    u = (0,0 | 78,0 | 9,0 | 1,0 | 2,0 | 1,0 | 0,0).

There are 169 ordered tuple pairs at each of the two shifts. Every one of their signed scores on u is nonzero modulo 101. Several distinct scores occur at each shift, so the raw contrast system also detects u. Yet for every source or target pin at both shifts, the sum of the scores is zero modulo 101. Thus each pin average on u equals zero with an invertible denominator. The raw contrast ranks are both 6, the individual cross-mean ranks are both 5, and the stacked cross-mean rank is still 5. The witness is invisible exactly at ±25 and ±26; each of the other 97 displacements has rank 6.

This rules out the implication from full tuple support, positive unwrapped pin counts, and two distinct nonzero sign classes to full stacked cross-mean rank. It also shows that violation of every individual zero-score test at those shifts can coexist with complete cancellation in every first pin moment. The tuple counts in this example are 1 or 2. The [pointwise seed](POINTWISE_SEEDS.md) has the additional relative concentration d(g) = K(1+O(δ)) at growing K, specific masks, and a deterministic selection rule. Its partial-rank theorem concerns genuine rows obtained from a Gram calculation. This finite calculation does not refute those claims or settle existence of an affordable universal selection rule.

## Two shifts can lose the same direction with full punctured support

On C₁₇ take seven indexed binary blocks (1,1+dᵢ) with steps (1,1,1,3,3,7,7). Every full tuple value occurs 4 through 12 times. After deleting any one block, **every** value occurs 1 through 8 times. At both tested displacements, all the ordered-pair pin masses are units modulo 17. Nevertheless the averaged equations at shifts 2 and 3 miss the same nonaffine direction:

    u = (11,0 | 11,0 | 11,0 | 16,0 | 16,0 | 1,0 | 0,0).

| Measurement over F₁₇ | Shift 2 | Shift 3 | Stacked |
|---|---:|---:|---:|
| Raw contrast rank | 6 | 6 | 6 |
| Cross-mean rank | 2 | 5 | 5 |

The raw rows at both shifts detect u. For the cross means its invisible shifts are exactly {2,3,14,15}. Shift 1 alone has rank 6. The masses at shift 2 have residues in {3,9,14,15}; at shift 3 the residues are in {4,5,7,8,14,15}. This gives a shared loss after modular averaging even when the raw per-shift theorem applies with ε = 0 and **all** means on the tested shifts have defined modular denominators.

A second fixture rules out every fixed pair of nonzero sign classes on C₇ across a family of dilated examples. Use five blocks with steps (1,1,2,2,3) and the same convention (1,1+dᵢ). Full counts lie in [4,5] and all single-deletion counts lie in [1,4]. Raw ranks are 4 at shifts 1 and 2. The two cross-mean ranks are 3 and 2, and the stacked rank is 3. The direction (6,0 | 6,0 | 1,0 | 1,0 | 0,0) is missed exactly at ±1 and ±2. Multiplying every block by a nonzero scalar permutes the displacements and preserves the rank calculations. Multiplication by 2 cycles the three nonzero sign classes, so each of their three possible pairs fails on some dilation. Unlike the C₁₇ case, some individual pin masses at shift 2 vanish.

These are indexed cyclic block fixtures. Their support and count bounds can be checked directly. Neither fixture reproduces the repository's asymptotic seed selection, its vanishing energy error, or its arity schedule. All group labels in the tests are supplied explicitly.

## A sufficient condition at one nonzero shift

For a positive pin in block i selecting A, write gᵢ(x) for the number of representations of the remaining signed slots at displacement x. For candidate table u, let hᵢ(x) be the sum of the remaining scores over these representations. Neither quantity depends on A. At shift s the pin has mass gᵢ(s−A) and its score sum is

    gᵢ(s−A) uᵢ(A) + hᵢ(s−A).

Fix t ≠ 0 and assume some positive or negative pin at t has nonzero mass in the field. For any u in ker R(t), all pin score sums at t equal their masses times a unique scalar c. Reversing each tuple pair yields the same identity at −t with scalar −c, including for individual pins with zero mass.

**Proposition.** In block i join A to A+2t when both belong to the block and gᵢ(−t−A) ≠ 0 modulo r. Write Cᵢ for the number of connected components, counting isolated points. Then

    dim ker R(t) ≤ 1 + Σᵢ Cᵢ.

In particular, if all block graphs are connected, the kernel consists exactly of a common additive slope and one constant per block.

**Proof.** Compare the positive pin at A at shift −t with the positive pin at A+2t at shift t. Their remaining displacement is the same x = −t−A. Subtracting the two equations cancels hᵢ(x). On each edge this gives uᵢ(A+2t)−uᵢ(A) = 2c. There is a unique additive map ℓ from the group to the field with ℓ(t) = c. Thus uᵢ−ℓ is constant on each component. The scalar c and these constants determine u. Every affine table satisfies the original cross means by substitution. ∎

For example, on C₁₁ take blocks (2,4,6) and (2,4,6,8,10,1,3). At t = 1 the indicated consecutive edges have nonzero common masses and the matrix has rank 7 on 10 columns. For sparse irregular blocks the graph can have many components. The proposition gives no general bound on their number for the repository's seed.

## Exact classification for a symmetric block family

Let every one of m ≥ 1 indexed blocks equal F = Fᵣ ∖ {0}, where r ≥ 5. Put h = (r−1)/2 and k = 2m−1 in Fᵣ. There are N = m(r−1) columns. At every nonzero t the raw contrast rank is N−m−1. For m ≥ 2 the tuple sums cover the group, so the full-support lemma applies. For m = 1 the individual transitions A → A+t on the nonzero vertices form a path whose constant increment gives the affine kernel directly.

| Parameter | Excess kernel dimension of one R(t), t ≠ 0, beyond the affine kernel | Effect of stacking |
|---|---:|---|
| k outside {0,1} | 0 | One nonzero displacement attains the raw rank. |
| k = 1 | m−1 | Two displacements s,t with s ≠ ±t attain the raw rank. |
| k = 0 | h−1 | Every nonzero displacement has exactly the same kernel; the zero displacement contributes no rows of positive rank. |

In the last case the common kernel consists exactly of tables uᵢ(A) = f(A)+bᵢ for a **common arbitrary odd function** f and arbitrary block constants. Only the one-dimensional space of additive functions among those odd functions gives affine tables. For instance f(A) = A³ is an explicit nonaffine witness at every shift for r ≥ 5. Taking m = (r+1)/2 gives an infinite family whose excess kernel dimension is (r−3)/2. All individual block deletions leave full sumset support for this choice of m. For these same parameters the mean integer tuple multiplicity (r−1)ᵐ/r exceeds r; this family lies outside the no-wrap hypothesis and the asymptotic arity schedule in [the resilient seed theorem](RESILIENT_SEEDS_AND_VARIANCE.md).

### Modular identities

In the group algebra let J be the sum of all group elements and let 1 be the basis element at the identity. Since J² = 0 in characteristic r and F has indicator J−1, the remaining-slot count is

    gᵢ(x) = k − 1{x=0}.

Indeed the (2m−1)-fold convolution is −1 + (2m−1)J. If k is outside {0,1}, every mass is nonzero. The step-2t cycle with the vertex zero removed is connected, so the proposition proves the first case.

For the exceptional cases, center each table on F by adding the block constant Σₐ uᵢ(a); its new sum on F is zero since |F| = −1 modulo r. Denote the centered table by wᵢ, extend it by wᵢ(0) = 0, and set D(x) = Σⱼ(wⱼ(x)−wⱼ(−x)). Centering leaves all signed tuple scores unchanged. For a positive pin (i,A) at s, with x = s−A, its score sum equals

    (k−1{x=0}) wᵢ(A) + D(x) − wᵢ(x).                         (1)

To see the second term, the unpinned score contributions consist of one copy of each positive centered table except i and the negatives of all reflected centered tables. Their coefficient arrays have sum zero. Convolving any of them with the (2m−2) remaining free-slot indicators acts as the identity: the convolution power is 1−(2m−2)J, and J annihilates every zero-sum array. This proves (1) also for m = 1. Negative pins at s are obtained by negating the positive score sums at −s and retaining their masses.

**When k = 0**, some pin mass is nonzero at each t ≠ 0, namely at A=t. At s in {t,−t} and A ≠ s, the pin mass is zero. Thus (1) forces wᵢ(x) = D(x). Every x ≠ 0 occurs in such a pin for at least one of s = ±t: choose s ≠ x and take A=s−x. Consequently every centered table equals the same odd function f. Conversely D = 2m f = f for any common odd f. Formula (1) then gives the required score sum with scalar f(t) at t and −f(t) at −t. At s = 0 every pin mass is zero, so every cross row vanishes. This proves the stated kernel and its dimension at all shifts.

**When k = 1**, fix t ≠ 0 and subtract the common additive slope ℓ with ℓ(t) = c. This slope is centered on F. Let zᵢ = wᵢ−ℓ. The graph in the proposition is the step-2t cycle with vertex 0 removed and the additional edge from −t to t deleted. In coordinates a = A/t represented by 1,…,r−1, the two path components are the odd and even representatives, each of size h. Thus zᵢ = αᵢ qₜ, where qₜ equals +1 on the odd representatives and −1 on the even ones; centeredness forces these component values to be opposites. The function qₜ is odd. For x ≠ 0,s and s = ±t one has qₜ(s−x) = qₜ(x), by parity of the canonical representatives. Substituting z into (1) shows that the zero-scalar equations at ±t hold exactly when Σᵢ αᵢ = 0. The endpoint x = 0 contributes zero on both sides. This gives precisely m−1 additional degrees of freedom at t.

For two shifts s,t, the slope of a common centered solution is fixed by its block average, since m = 1 in Fᵣ and Σᵢ αᵢ = 0. A nonzero extra direction would require qₛ and qₜ to be proportional. Proportionality means equality up to sign. The sum in Fᵣ of the arguments on which qₜ is positive is t h², and the negative set has sum −t h². Since h² ≠ 0, proportionality requires s = ±t. This proves the two-shift conclusion.

## A separate example with no wrap in the full tuple counts

On C₅ take the blocks (1,4), (1,4), (1,2,3,4). The integer full-tuple counts at values 0 through 4 are (2,4,3,3,4). Thus every fiber is present and every count is strictly below 5. Each displacement also has at least one nonzero modular ordered-pair pin mass. In column order as displayed, take

    u = (4,0 | 4,0 | 1,1,0,0).

Every nonzero displacement has raw contrast rank 4 and pin-mean rank 3. Their combined pin-mean rank is still 3, and all of them miss u. At displacement zero a pin-mean row detects u, bringing the complete rank to 4. Deleting the last block leaves only three of five group values; this example does not satisfy the dense-puncture hypothesis. It distinguishes no wrap in the single-tuple counts from the stronger simultaneous seed requirements. Ordered-pair pin counts can wrap even here.

## Verification and scope

Run `python -B check_simultaneous_pin_loss.py` in this directory. The [verifier](check_simultaneous_pin_loss.py) computes exact integer convolution counts before reduction. It compares every mass and every entry of the positive and negative first moments with literal enumeration on five fixtures. The rank checks exercise all three parameter cases, all nonzero displacements, and every pair of distinct nonzero sign classes in each case. The binary-block fixtures check all single-deletion counts, the two raw ranks, the stacked cross-mean rank, and the complete set of displacements missing each witness. The C₁₇ test checks every pin residue at the four signed shifts. The unequal-block control checks the witness against the individual raw and compressed rows. [Recorded output](simultaneous_pin_loss_verification.txt)

The C₁₀₁ check enumerates all tuple pairs and directly verifies the integer mass range over the full grid, all 338 nonzero raw scores at the two named shifts, every vanishing pinned score sum there, and the exact ranks at all 101 displacements.

The general deterministic selection question for the deletion-resilient seed remains open. The complete shared-kernel obstruction uses a special congruence and permits high integer multiplicities. Its pin masses are supported at a single displacement for each pin, so it also lies outside the earlier all-displacement cross-mean theorem's dense modular pin-support hypothesis. The unequal-block control satisfies the full-tuple bound alone. The single-shift positive result specifies an additional overlap condition that would need a separate argument on the selected seed.
