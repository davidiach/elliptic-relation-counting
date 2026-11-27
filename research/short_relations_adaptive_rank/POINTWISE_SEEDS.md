# Centered collision moments, pointwise coverage, and guaranteed partial modular rank recovery

Research follow-up to `DISPLACEMENT_KERNELS.md`, `RESIDUAL_KERNELS.md`, and `ANCHORED_TRANSLATIONS.md` in PR #4, based on e9ff1b727c55414f9b3330b53933b0ffa2f0d165. 27 November 2026. Conditional on the closed Boolean tensor contraction theorem stated in Liu–Zhang–Ren (LZR), Theorem 2.2, ePrint 2026/3020, for all algorithmic costs below. The counting/linear-algebra estimates themselves are elementary. This is a pointwise seed construction and a **partial** rank-recovery theorem, not a faster end-to-end deterministic discrete-logarithm algorithm.

## 1. A centered moment bound

Let G be a cyclic group of prime order r, and let the ambient point set be G minus its identity. In each of m blocks choose independent random retention seeds with marginal retention probability theta and h-wise independent indicator variables on distinct ambient points, where h is even. Distinct blocks are independent. On a prime-order elliptic curve over F_p, p>3, these indicators are supplied by the parity-specific degree-h polynomial masks of `RESILIENT_SEEDS_AND_VARIANCE.md`: fix a nonzero leading coefficient and choose the h lower coefficients independently uniform in F_p, separately for the two y parities in each block. Retain a point when the canonical evaluation on its x coordinate lies in [0,T). Thus theta=T/p. The Vandermonde argument gives the required h-wise independence. Each selected block has at most 2hT points. In the large-input regime h<p.

Write d(g) for the number of retained one-per-block tuples with sum g; mu_g=E[d(g)]. Set q'=r theta and mu'=(q')^m/r. Assume q',mu'>=K>=1. Then

    E[(d(g)-mu_g)^h] <= 2^h B_h^m (mu')^h K^(-h/4),

where B_h is the h-th Bell number. This bound holds for **every** target g. Unlike an uncentered collision bound, it improves exponentially in the moment order. The proof does not treat distinct tuple events as independent when they share a point.

### Proof: count only nonisolated tuple patterns

For each ambient tuple X of sum g, let I_X be its retention indicator and J_X=I_X-theta^m. Expand the centered h-th moment as the sum of E[prod J_X] over ordered h-tuples of ambient tuples of sum g. In block i, their point-equality partition pi_i has ell_i classes. Set ell=sum_i ell_i and C=mh-ell. A tuple is isolated if its point is different from every other tuple's point in every block. Its J factor is independent of the others by h-wise independence, so the corresponding expectation is zero. For any pattern, expand the product of centered factors. A subset of retained factors with ell_U distinct point variables has weight theta^(ell_U+m(h-|U|)); since ell_U+m(h-|U|)>=ell, the absolute expectation is at most 2^h theta^ell.

For fixed partitions, let M be the h-by-ell incidence matrix over F_r whose columns are the class indicators in the h tuple positions. The sum equations are Mx=g 1 in abstract group scalar coordinates, used only in the proof. Put R=rank M and D=h-R. There are at most r^(ell-R) assignments, ignoring forbidden identity and distinctness constraints, so the normalized contribution is bounded by

    (q')^ell/r^R / (mu')^h = (mu')^(-D) (q')^(-(C-mD)).

Let S be the union of supports of vectors in ker(M^T), and s=|S|. In every block, any class that meets S contains at least two indices of S: otherwise some kernel vector would have a nonzero sum on that class. The number of mergers in the restriction of pi_i to S is therefore at least s/2, and is also at least D by rank-nullity on the restricted class indicators. Every index outside S must share a class somewhere if the tuple pattern has no isolated tuple. The mergers outside the restriction to S, summed over blocks, number at least (h-s)/2: a class meeting S counts every added outside index, and a pure-outside class of size l>=2 contributes l-1>=l/2. Consequently

    C >= m max(D,s/2) + (h-s)/2.

In particular C>=mD. If D>=h/4, the normalized contribution is <=K^(-h/4). If D<h/4, the displayed inequality gives C-mD >= (h-2D)/2 >= h/4: for s<=2D this follows from the last term, and for s>2D add m(s-2D)/2. Thus the same bound holds. There are at most B_h^m choices of partitions. Summing the absolute centered contributions proves the claim. No factor proportional to the number of targets has been hidden inside a single-target estimate.

## 2. Deterministic seed with pointwise-uniform short fibers

Let b=ceil(log2 p), L=log2 b, and take m=Theta(L^(1/8)), h=1024m (even), K=2^floor(b/(64m)), n=r-1. Define q*=ceil((Kr)^(1/m)), T=ceil(p q*/n), q=nT/p. Asymptotically T<=p and h<p, q=q*+O(1), q^m/r=(1+o(1))K, q'>=K and mu'>=K. Finite initial exceptions can be treated by exhaustive group arithmetic; constants are intentionally crude.

For the exact mean use the number of ambient tuples in (G minus {O})^m: it is [n^m+(r-1)(-1)^m]/r at g=O, and [n^m-(-1)^m]/r elsewhere. This follows by convolution, or character orthogonality. Multiplying by theta^m gives the exact two values of mu_g. They are positive in the asymptotic range, and mu'/mu_g=1+O(m/r). Hence for the nonnegative potential

    Phi = sum_g ((d(g)-mu_g)/mu_g)^h

we have E Phi <= r * 2^h B_h^m (1+O(m/r))^h K^(-h/4) <= 2^(-3b+o(b)), using B_h<=h^h, Hasse, and m^2 log m=o(b). With delta=2^(-floor(b/(1024m))), delta^h>=2^(-b), so E Phi < delta^h eventually. The bitwise conditional-expectation algorithm minimizes this potential over seed completions. It needs exact counts of sum_{g=O} d(g)^j and sum_{g!=O} d(g)^j for j=0,...,h, expanded by the binomial formula with exact rational coefficients. These are ordinary counts of j tuples sharing a sum (and the appropriate identity test), with the seed restricted to its current prefix. For j=0 the sum is just the number of applicable targets times the number of seed completions. No oracle for centered or signed counts is assumed; subtract only after the nonnegative integer counts are obtained. Valid seed completion counts are products of allowed word-interval sizes. Rational denominators and coefficients have polynomial bit length in h,m,b. Choosing the lower conditional expectation over nonempty prefix branches gives Phi<delta^h and therefore

    |d(g)-mu_g| < delta mu_g   for every g.

The rounding errors in q, q* and mu_g are o(delta) relative to K. Thus d(g)=K(1+O(delta)) pointwise; in particular all fibers are nonempty and strictly below r. Put W=prod_i n_i=sum_g d(g), with n_i the selected block sizes. Then W=rK(1+O(delta)); n_i<=2hT and the product lower bound gives log2 min_i n_i >= b/m-O(m log m) (the positive log K/m term may be discarded). These imply K^2/min_i n_i=2^(-Omega(b/m)). This seed uses less redundancy than a seed required to cover after deleting a fixed fraction of blocks.

**Bit budget.** Horner evaluation of each degree-h membership predicate costs O(h) modular arithmetic gadgets, hence O(hb) Boolean edges with degree O(b^(3/4)). At most h tuples and m blocks give O(h^2 mb) total edges and polynomially many vertices. Bounded-fanout equality trees and subdivision preserve this order; canonical seed sharing and coordinate validation are charged. Exact integer count precision is linear in the network edge budget. Local arithmetic tensor entries are uniformly polynomial-time evaluable. The closed LZR theorem therefore gives log bit time O(h^2 mb/sqrt(log(h^2 mb))) per count; O(h^2 mb) prefix/count calls are only a polynomial overhead. Enumerating the selected points by exact prefix queries (and computing the curve order) is charged as in LZR, and the size of the enumerated factor bases is 2^{O(b/m+log b)}. Balancing b/m and m^3 b/sqrt(log b) gives m=Theta(L^(1/8)). Total time is 2^{O(b/L^(1/8))}, and space at most 2^{O(b/L^(1/8))}. The polynomial masks are of growing degree here; the fixed-degree quarter-power estimate from the previous seed theorem cannot be transferred without charging the extra work. No faster multipoint-evaluation compiler is assumed. This is an intermediate **seed** exponent, not a deterministic ECDLP exponent.

## 3. All nonzero displacements have full raw rank and nonzero modular pin masses

For a displacement t consider all ordered pairs X,Y of short tuples with val(Y)-val(X)=t. For a pin fixing A in slot i on one side, its count is a sum of opposite-half fiber counts over W/n_i tuples. Thus, simultaneously for all t and all one-slot pins h,

    n_h(t) = (1+O(delta)) K W/n_i,       0 < n_h(t) < r.

In fact every consistent nonempty collection of pins on one side has positive count below r: the opposite half remains unrestricted and covers G, and the count is bounded by a one-pin count. Total ordered-pair count M_t=(1+O(delta))rK^2. Consequently the modular pin masses never vanish, even at t=0. The two possible causes of rank loss in `DISPLACEMENT_KERNELS.md` are now reduced to the finite-field averaging radical for any prescribed **nonzero** displacement.

To prove the raw rank claim without punctured supports, let A_t have signed incidence rows c_Y-c_X, and suppose A_t u=0 over F_r. Every tuple fiber is nonempty, so comparisons to a common successor imply the score of a tuple depends only on its value g, defining f(g). The homogeneous constraints give f(g+t)=f(g) everywhere. A nonzero t generates the prime-order group; f is constant. Varying one slot forces u constant within each block. Thus A_t has rank N0-m, N0=sum n_i. Similarly if all row scores are the same (the raw contrasts kill u), f(g+t)-f(g) is a constant; as in DISPLACEMENT_KERNELS.md this makes f affine and the raw contrast kernel consists of block constants plus the log direction. Its rank is N0-m-1. The contrast conclusion at t=0 is NOT claimed from full support alone. For actual row equations choose t=P (known RHS one), rather than an arbitrary displacement with unknown RHS.

### The tradeoff with punctured coverage

If the deletion of any one block leaves support >=(1-epsilon)r, then W/n_i >=(1-epsilon)r. For any fixed one-slot pin h in block i, sum_t n_h(t)=W^2/n_i. Its average is therefore >=(1-epsilon)W >=(1-epsilon)^2 r n_i. As soon as (1-epsilon)^2 n_i>=1, some displacement has pin mass >=r. Thus insisting that all one-pin counts be below r at every displacement is incompatible with dense punctured support in the growing-block regime. This does not force a modular zero: a count above r may be a unit modulo r. It explains why the previous punctured-support theorem could not supply the stronger integer no-wrap statement for ordered pairs. The present construction uses full support and a known nonzero shift, not the earlier all-displacement contrast theorem including zero.

## 4. A uniform partial modular-rank theorem

Drop one reference column per block from A_P. Equivalently the unknowns are label differences from these references; block constants cancel from every original signed row. Let A denote the reduced integer matrix, dimension D=N0-m columns. By section 3, A has full column rank over F_r and each true row has RHS one. Its integer Gram G=A^T A is positive definite. Every diagonal entry is at most the sum of the source and target pin counts for that column, hence <=2(1+O(delta))rK^2/min_i n_i. By Hadamard,

    v_r(det G) <= log_r(det G)
                <= D(1-31/(32m)+o(1/m)).

The error uses log r=b+O(1), log K=b/(64m)+O(1), and m^2 log m=o(b). All determinant and matrix operations have polynomial bit cost in D,b. The count network never explicitly lists A: any Gram entry, including restriction to a row-key prefix, is a signed sum of ordinary counts with at most two point pins on the 2m-slot relation. The masks cost O(hmb) edges for this query; the seed step is more expensive.

Apply the adaptive first and r-adic quadratic extractor of RESEARCH.md section 3. For clarity, if W0 is the row space of extracted genuine equations, choose an integer lift U of a field basis of its orthogonal complement. If GU is not zero mod r, use a nonzero sum of a_j(a dot u) mod r and bisect prefixes to extract a genuine row. If U^T G U is nonzero mod r^2, some basis vector or sum of two basis vectors has nonzero quadratic mass mod r^2, by polarization and odd r. Bisect that mass; a legal zero modulo r would have square zero modulo r^2, so the extracted row increases rank. If both tests fail with d=D-dim W0, extend U to an integer square matrix invertible mod r. In the congruent Gram, the leading d-by-d block is divisible by r^2 and both adjacent off-diagonal blocks by r. Factoring r from each of the first d rows and columns gives v_r(det G)>=2d. Thus the algorithm always extracts at least

    (1/2 + 31/(64m) + o(1/m)) D

independent genuine modular equations, in particular at least (1/2+1/(3m))D for sufficiently large b. The estimate is an integer-rank inequality with harmless rounding. It is **not** a guarantee of full rank. All loops and prefix counts are polynomial in N0,b; seed plus extraction fit the 1/8-saving budget. Appending the explicit consistency/normalization rows can help but is not needed for this partial rank bound. If the residual kernel does not vanish, the theorem supplies no faster solution of the DLP.

There is also a weaker guarantee for the direct pin averages. Write H for the source/target pin indicator matrix, n=H^T1, V=H^T A. Since A is a difference of columns of H, rank(V)>=rank(G mod r). If any n_h is nonzero (all are here), span(V) lies in the span of the cross-multiplied conditional means and one anchor row, so their rank is at least rank(G mod r)-1. The rank deficiency of G is at most v_r(det G) (e.g. by Smith form). Consequently direct pin averages have rank >=(31/(32m)+o(1/m))D-1. This quantitative guarantee is not a noncancellation theorem. The radical can still exist; neither integer positive definiteness nor nonzero pin masses makes a Gram nondegenerate over a finite field.

## 5. Finite verification and unresolved frontier

Run `python -B check_centered.py` (NumPy for small exact integer arrays; SymPy for exact determinants). `pointwise_verification.txt` records 4,965 nonisolated partition patterns with 14,895 rank/support checks over three primes; exact centered-moment expansion against 4,608 independent Bernoulli-mask seeds; 150 full-coverage/no-wrap fixed-shift toy fixtures, 40 singular Gram matrices, all 40 rescued by adaptive extraction; and an actual E/F5, y^2=x^3+2x+1, group order 7, with all eleven group-law relations checked. On that curve the Gram rank is 3/4 and adaptive extraction reaches 4/4. No contrast pin-averaging rank loss occurred in the 150 tests. They do not test the enormous asymptotic moment order, implement the contractor, or establish universal full rank. The proof does not infer a universal theorem from the finite results.

An exploratory random search through group order 23 also found no cross-mean counterexample in 10,710 qualifying fixed-shift tests; the subsequent larger sweep was stopped by a runtime limit. This is exploratory evidence, not a completed exhaustive search or an additional formal guarantee. It motivates studying the stricter pointwise/no-wrap regime separately from the punctured-support regime. In particular the next task is to resolve the remaining radical or give a counterexample using the strengthened promises, and to determine whether a uniform joint seed-and-rank theorem is possible. The full deterministic logarithmic-saving exponent remains open. Independent proof review and the upstream contraction-premise audit remain outstanding.
