# Integer-basis extraction and conditioned completion volume

27 November 2026. Follow-up to POINTWISE_SEEDS.md at 3175664ba4f869e8f3a6b6dbd562254ef17c96cf, developed while the branch acquired SIMULTANEOUS_PIN_LOSS.md. The latter is complementary and supplies a concrete obstruction to simultaneous pin averaging. The new argument uses exact integer Gram counts, not noncancellation of pin averages. All asymptotic algorithmic statements use the closed Boolean contraction premise of Liu–Zhang–Ren (LZR), Theorem 2.2, ePrint 2026/3020. No full faster deterministic discrete-logarithm algorithm is proved.

## 1. Near-complete modular rank without a modular extraction oracle

Let A be an integer M-by-D matrix, full column rank over F_r for an odd prime r, with each row having squared Euclidean norm at most s. Rows are implicitly keyed by binary strings of length B. Assume one can compute the exact integer Gram G_I=A_I^T A_I for any row-key prefix/interval I, and extract a row at a singleton key. The matrix need not be listed. In the intended application s=O(m).

**Theorem.** D deterministic integer-basis selections recover a genuine row set whose modular rank is at least

    D - floor(D log_r(s)/2).

The procedure uses O(DB) prefix Gram queries (each query returning D^2 integer entries), and polynomial-in-D,B,entry-bit-length exact integer linear algebra.

Proof. Maintain a rationally independent selected row set W. If its rank is less than D, find a nonzero integer vector u in its rational nullspace by exact elimination and clearing denominators. Since A has full rank even over F_r, it has full rank over Q, hence u^T A^T A u>0. At a prefix the same quadratic form is the exact NONNEGATIVE INTEGER sum of (a dot u)^2 for its rows. Bisect the key tree, preserving strict positivity; at a leaf the selected genuine row has nonzero dot product with u, and increases rational rank. No modular reduction is made in this extraction. After D selections form the square integer matrix B0. Its nonzero determinant satisfies

    |det B0| <= product(row norms) <= s^(D/2).

If its modular rank is D-d, elementary divisors (or elimination by operations invertible over the r-local integers) imply r^d divides det B0. Therefore d<=log_r |det B0|<=D log_r(s)/2. The bit lengths of nullspace lifts, determinants and intermediate exact elimination outputs are polynomial in the matrix dimensions and input bit lengths, by minor bounds and fraction-free elimination. The procedure may output a modularly singular rational basis. It does not presume integer independence implies modular independence.

## 2. Application at the original quarter-power scale

Use LZR Stage I with m=Theta((log b)^(1/4)), b=ceil(log2 p), in the common whole prime-order elliptic-curve group G=<P>, |G|=r. The one-per-block sumset has density >1/2, so the doubled tuples (two points per block) cover all group elements. Use every ordered doubled-tuple pair with val(Y)-val(X)=P, and reduce the signed incidence matrix by deleting one reference column per block. This makes D=N0-m difference variables; each row has known right side one. The true differences solve the system.

Full modular rank follows directly: if every row kills a homogeneous direction, compare source tuples of the same value to a common successor, which exists by coverage. Scores are a function of group value; translation by P leaves that function unchanged. Since P generates G, it is constant. Varying one slot forces the labels to be constant within each block. In reduced coordinates the reference labels are zero, so the direction vanishes. This proof applies to the short pointwise-coverage tuples as well, but does NOT require them.

A doubled-tuple incidence difference has squared norm <=8m (two positive and two negative occurrences per block; merging equal-sign occurrences can only give the worst-case 4+4). Column deletion cannot increase the bound. Each Gram entry is an exact signed linear combination of O(m^2) ordinary counts with at most two pinned occurrences; prefix restrictions and scalar group-law checks are ordinary Boolean arithmetic gadgets. The precision is O(mb) bits for the nonnegative counts. Local membership is the original parity-affine predicate. Hence each count has O(mb) Boolean edges and the degree and entry-access bounds of LZR's closed contractor. Setup, explicit factor enumeration, the polynomial number of Gram entries and queries, exact elimination and storage all fit

    T = 2^{O(b/(log b)^(1/4))},
    S = 2^{O(b/(log b)^(1/4))}.

The recovered modular rank is at least

    D - floor(D log_r(8m)/2).

Thus a deterministic procedure recovers a 1-O(log m/b) fraction of all independent difference equations. This is stronger than the half-rank bound in POINTWISE_SEEDS.md and restores the original quarter-power scale rather than its slower growing-degree seed construction. It also holds for the pointwise short-tuple version with s=2m. This rank statement is not an end-to-end DLP time bound. In particular D is itself subexponentially large, so the absolute number of unresolved directions need not be small.

For target Q, find one doubled tuple of value O and one of value Q by ordinary prefix decomposition. Their signed reduced incidence gives a known linear functional on the true difference table; if it lies in the recovered modular row span, evaluate it from the known right sides and verify by scalar multiplication. Otherwise the rank guarantee alone is insufficient. If the full difference rank is obtained this recovers the target without knowing the blockwise additive constants.

## 3. Conditioned completion, a genuine recursive bound

Let W be k genuine modularly independent rows already selected (therefore rationally independent), d=D-k, and G=A^T A. Define

    Delta(W) = det(W W^T)   (Delta(empty)=1),
    C(W) = [z^k] det(G+z W^T W)
         = (-1)^k det([[G,W^T],[W,0]]).

Cauchy–Binet gives C(W) as the sum of squared D-by-D minors extending the fixed rows W by d distinct row indices of A. This identity holds even for a singular suffix Gram. The bordered determinant identity follows first for invertible G by the matrix determinant lemma, then in general by polynomial identity. All computations are exact integer operations.

Choose an orthonormal basis V of the real nullspace of W. In an orthonormal coordinate system adapted to row(W), the coefficient of z^k equals Delta(W) det(V^T G V), so

    0<C(W)<=Delta(W) Lambda^d,

where Lambda is any upper bound on the largest eigenvalue of G. If each row has at most s nonzero entries, Cauchy–Schwarz gives x^T G x <= s max_j(G_jj) ||x||^2, so that expression is an admissible Lambda. (For the general case distinguish the support bound from the squared norm bound; they can both be bounded by s for the integer rows here.) For an individual rationally independent completion B of W by d genuine rows, orthogonal projection followed by Hadamard gives the sharper

    0<det(B)^2<=Delta(W) s^d.

This latter bound omits the number of candidate rows and is why integer positivity is especially useful. In particular if Delta(W) s^d < r^(2d), any rational completion increases modular rank; if Delta(W) s^d < r^2 it is already a full modular basis. More precisely its new deficiency is at most floor((log_r Delta(W)+d log_r s)/2). The bound is directly checkable using integer comparisons if desired. Delta(W) can be replaced by the squared covolume of the saturated integer lattice in the rational span of W: its index over W is coprime to r, so its reduction has the same modular row span and its genuine right sides can be transported by exact rational elimination with denominators prime to r. Computing a basis of that saturation is standard integer linear algebra. This does not assert that its volume is small.

For the earlier adaptive modular quadratic extractor, if it stalls, let U be integer lifts of a field basis of the remaining kernel. Then GU=0 mod r, U^T G U=0 mod r^2 and WU=0 mod r. Extend U to an integer D-by-D change of coordinates invertible modulo r. In the congruence-transformed matrix G+zW^T W, the first d diagonal-block entries are divisible by r^2, and adjacent off-diagonal blocks by r. Factoring r from each of the first d rows and columns shows **every coefficient**, in particular C(W), is divisible by r^(2d). Since Delta(W)<=s^k, any stall obeys

    d <= D c/(2-a+c),    c=log_r s, a=log_r Lambda,

provided 2-a+c>0. For the pointwise seed Lambda<=2m max G_jj<r eventually, as in POINTWISE_SEEDS.md. This independently strengthens its original modular-extractor stall bound from roughly D/2 to O(D log m/b). It is a true conditioned estimate rather than repeatedly applying the unconditioned determinant bound to a projected matrix with uncontrolled coefficients.

A different exact completion uses C(W) itself. Put v=ord_r C(W), and preserve nonzero squared-minor mass modulo r^(v+1). For any suffix of row keys use its Gram instead of G; difference of two suffix masses gives the mass of completions whose least additional key lies in the corresponding interval. Bisect this quantity, append the selected least row, discard its key and all smaller ones, and repeat. The fixed rows along the search need only be rationally independent, not modularly independent. At the final leaf the mass is det(B)^2, of valuation <=v. This gives a deterministically found completion with new deficiency at most floor(v/2). Retain W and a modularly independent subset of the new rows. With a,c above, the bound on C(W) yields the recurrence

    d_new <= min(d, (a/2)d+(c/2)(D-d)).

If a<1 this contracts until a residual dimension O(D c); it does not yield convergence to zero. The real-completion theorem in section 1 is simpler and improves that worst-case scale's constant. Both methods use polynomially many ordinary exact Gram queries and integer operations. Neither provides a cheap test of whether a remaining nonzero modular direction violates at least one equation. In particular the aggregate completion mass can be divisible by r^2 even if every nonzero individual completion determinant is a unit.

## 4. Exact obstructions and finite verification

The verifier check_integer_completion.py uses exact integer/rational arithmetic (SymPy) on toy matrices, not an implementation of the asymptotic contractor. It checks the bordered mass against direct sums of squared minors, positivity-based completion, modular suffix selection, the individual determinant and aggregate volume bounds, complete fixed-shift cyclic examples, and doubled-tuple examples. See integer_completion_verification.txt.

A genuine fixed-shift row SUBSET on E/F7, y^2=x^3+3, group order 13, has six reduced binary-block coordinates, true differences (3,4,5,8,6,7), and 12 distinct signed rows. It can be realized using blocks {P,[1+lambda_i]P}, with a fixed reference choice for zero slots. The selected W has rank 5 and the full subset rank 6. Every column Gram diagonal is <13; the quadratic mass on its sole residual direction is 845=5*13^2 and the first modular tests vanish. C(W)=169 and Delta(W)=36, with s=4. Thus the old modular quadratic extractor can stall and the aggregate completion mass is divisible by 13^2, but Delta(W)*s=144<169: any rational completion recovers the remaining modular rank. All 12 group-law relations are directly checked. This is not the complete fixed-shift row family of those binary blocks or an example of the asymptotic selected seed. It is a counterexample to inferring success of the old modular tests from small row entries, nonwrapped column counts and genuine known-right-side relations alone.

A separate complete C13 fixed-shift fixture with three three-point blocks has full rank 6 and 55 rows. A displayed rational basis of six actual rows has determinant +/-13 and rank 5 modulo 13. Ordering those six keys first makes a positivity-only rational-basis selection return that basis. This fixture permits wrapped counts and is NOT an obstruction to the specific canonical order or stronger pointwise seed promise. It demonstrates why rational rank recovery alone cannot be claimed to solve the modular problem. Concurrent SIMULTANEOUS_PIN_LOSS.md has a more relevant full-support/no-wrap example for modular pin averages; it is not a counterexample to the conditioned integer bounds.

## 5. Status and next question

The previous recursive half-rank target is substantially improved: the near-complete modular-rank result requires neither high-variance assumptions nor pointwise seed construction. It uses classical exact Gram positivity and Hadamard, rather than a claimed novel general linear-algebra theorem. The application is a conditional resource/rank statement for this implicit row family. Literature priority and independent proof review remain open.

The remaining challenge is to exploit the small relative residual dimension and the arithmetic structure of the complete relation family to obtain a modularly independent completion, or a target-only certificate. A possible next route is to construct low-volume genuine modular bases (or low-dimensional rank-preserving local subfamilies) for which the individual volume criterion proves further progress, without charging an exponential search for them. There is no current theorem guaranteeing the requisite volume. A single unresolved scalar direction may still encode the target discrete logarithm. Consequently no faster uniform deterministic ECDLP exponent is claimed.

References: LZR ePrint 2026/3020, Theorems 2.2 and 3.1 and Proposition 4.1; POINTWISE_SEEDS.md, RESEARCH.md and SIMULTANEOUS_PIN_LOSS.md in PR #4. Cauchy–Binet, Hadamard, the determinant lemma, Smith/Hermite normal forms and integer nullspace elimination are classical.
