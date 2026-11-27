# Moment/rank boundary: Cartesian support bounds and sparse sign separation

27 November 2026. Follow-up to [PR #16](https://github.com/davidiach/elliptic-relation-counting/pull/16), head `d9766c7932a063c708b498d3cd0dff5d0657db71`, reviewed against main `2f1fca860a6f02d1bc9c9796b560b02b526e6199` (merged PR #17). PR #16 was subsequently merged at main `444ee2fe541062f15b5f6af39a8439d8a258a186`; the follow-up is submitted against that merge. This note concerns abstract integer row families and exact coordinate moments, not an algorithm or a lower bound for the actual elliptic relation family. The proofs use classical interpolation, coding-theory and graph facts. Independent priority and proof review remain open.

## 1. A stronger necessary condition: number of distinct violating rows

For a multiset of integer vectors H, its moments through degree t are the sums of every coordinate monomial of total degree at most t, including degree zero. Let r be prime. Suppose H and J have the same exact moments through t >= 1 and every H row annihilates a direction u modulo r. Let S be the set of distinct J rows a with a dot u nonzero modulo r and multiplicity n_J(a) nonzero modulo r. Assume S is nonempty. Suppose each coordinate takes values from an alphabet of at most Q distinct residues modulo r, with distinct integer vectors in S giving distinct residue vectors (in particular, this holds for binary or ternary alphabets when r is odd).

Let q=min(Q,r), since the field has only r residues, and let D be the number of coordinates. The case q=1 forces S to be empty by the constant test below, so assume q>=2.

**Theorem 1.** Writing t = k(q-1) + e, 0 <= e < q-1, one has

    |S| >= q^k (e+1).

If t exceeds D(q-1), S must be empty. For binary rows the bound is 2^t. For ternary rows over an odd prime field it is 3^floor(t/2) times (1 + (t mod 2)). In particular, under the no-wrap assumption |J| < r, all violating rows are counted in S. Unlike the earlier score-diversity bound, this is a bound on distinct rows, not distinct scalar scores.

Proof. Set c_a = n_J(a)(a dot u) mod r for a in S. For every polynomial f of degree <= t-1, moment equality applied to (x dot u) f(x) gives

    sum_{a in S} c_a f(a) = 0 mod r.

All c_a are nonzero. The following elementary Cartesian dual-distance lemma proves the claim. A nonzero weighted population on a product of alphabets of size <= q, with nonzero weight at each distinct occupied point, annihilating all polynomials of degree < t has support at least F_q(t)=q^k(e+1), with F_q(t)=1 for t<=0. For positive t its support cannot be a singleton, because the constant polynomial is tested. Pick a coordinate with l >= 2 occupied values. A Lagrange selector of degree l-1 isolates any of its slices, so each nonempty slice annihilates polynomials of degree < t-(l-1). Induct on t, treating nonpositive remaining degree as the vacuous condition. The total support is at least l F_q(t-l+1). This is >= F_q(t): F_q(t) is the minimum product of (a_i+1) over positive integers a_i <= q-1 with sum >= t. To see the minimization formula, reduce the sum to t, then merge two partial parts whenever possible, or replace two partial parts a,b with q-1 and a+b-(q-1) when a+b>q-1. The product does not increase: the relevant differences are ab and (q-1-a)(q-1-b), respectively. Hence at most one part is partial. This proves the induction. Alternatively, a tensor product of one-coordinate Lagrange selectors distinguishes every row when degree D(q-1) is available. QED.

For 1 <= t <= D(q-1), the abstract weighted-support bound is attained on a grid with k alphabets of size q and, if e>0, one more of size e+1. These alphabets exist because q<=r. On those coordinates take the product of one-dimensional finite-difference weight vectors and fix any remaining coordinates. Lagrange denominators or finite differences provide nonzero weights annihilating all degrees smaller than each local size minus one. Every monomial of total degree < t falls below at least one local degree. This construction uses unrestricted field-valued weights. Sharpness for weights n_J(a)(a dot u) from two matched integer populations requires a separate argument. This is the classical dual-distance phenomenon for Cartesian evaluation codes; see [Lopez, Renteria and Villarreal, Affine Cartesian codes](https://arxiv.org/abs/1202.0085). The contribution here is the application to the violating-row population and the explicit separation of count wrapping from row support.

## 2. Positive boundary: nonnegative polynomial rank certificate

Let u be a direction over F_r for an odd prime r. Project onto its nonzero coordinates, ignoring all other coordinates, and write W for an upper bound on the projected integer l1 norm of any row. Assume W >= 2 and r > W. We use rational polynomials evaluated on the integer ball

    X_W = {x in Z^d : sum_i |x_i| <= W}.

**Theorem 2.** There is a rational polynomial P_u of total degree at most 2W-2 with P_u(x) >= 0 on X_W and

    P_u(x) = 0 iff x dot u = 0 mod r.

In particular two nonnegative integer multisets in this ball with equal exact moments through degree 2W-2 have the same status of annihilating u, independently of multiplicity wrapping. The same conclusion holds for nonnegative rational weights. The polynomial depends on u and W; its construction here is not claimed efficient in d,W.

Proof. Partition X_W by absolute-value pattern alpha, where alpha is a nonnegative integer vector with sum alpha_i <= W. Write S=supp(alpha). At each nonzero pattern, a sign vector eps in {-1,1}^S specifies x. Since r is odd and every nonzero magnitude is < r, flipping any one sign changes the linear score by a nonzero residue. Thus at least half the sign patterns violate the equation. Let v_alpha denote their number and B(alpha)=sum_i alpha_i^2. Define the desired function on the alpha pattern to be B(alpha)*2^|S|/v_alpha on violating signs and zero on legal signs. At alpha=0 take zero. Its average on signs is B(alpha), it is nonnegative, and it vanishes precisely on the legal scores. It is invariant under simultaneous reversal of all signs, so all its nonconstant Fourier characters on the Boolean sign cube are indexed by nonempty even subsets R of S. Each such R has |R|>=2.

We show that the function 1_{|x|=alpha} prod_{i in R} sign(x_i) is representable on X_W by a polynomial of degree <= 2W-2. Start with

    G_(alpha,R)(x) = product_{i in S} [ x_i^(2-1_{i in R})
                            product_{j=1}^{alpha_i-1}(x_i^2-j^2) ].

It has degree 2 sum_i alpha_i - |R| <= 2W-2, and vanishes unless |x| >= alpha coordinatewise. At a magnitude pattern beta >= alpha, it is a constant G_(alpha,R)(beta) times the sign character indexed by R. The constant at beta=alpha is nonzero over Q. Proceed in decreasing order of sum beta_i: subtract the already represented sign-character indicators at every strict larger beta, with coefficient G_(alpha,R)(beta), then divide by G_(alpha,R)(alpha). This isolates the desired pattern and does not increase the degree. The induction starts at maximal patterns and is finite.

Now sum the Fourier expansion of the desired function over patterns. The constant terms across all patterns are exactly the single quadratic polynomial sum_i x_i^2. All nonconstant terms use the just-constructed degree-(2W-2) indicators. This constructs P_u. If projecting from a larger ambient dimension, simply ignore the inactive variables in the polynomial. QED.

For W=1, the quadratic sum of active coordinate squares already works, so the uniform bound can be stated as max(2,2W-2). For signed unit rows W is just the maximum number of nonzero active entries, and with larger signed integer coefficients it measures active absolute mass. Without r>W the proof can fail if some nonzero magnitudes reduce to zero. The homogeneous score assumption is essential for the sign symmetry as stated.

**Moment consequence.** A complete degree-(2W-2) exact moment table determines sum P_u(a). The sum is strictly positive exactly when some row of positive weight violates u. For integer populations, if there is a violation the sum is >=1, since P_u(x) on a violating row is at least B(alpha)>=1. This observation is independent of any prime-adic cancellations. It does not provide a cheap way to obtain or evaluate the complete moment table: an explicit polynomial may have up to binomial(d+2W-2,2W-2) coefficients. No target resource bound is inferred. In contrast, arbitrary higher-order moment equality can persist in the binary additive-coset examples of PR #16, whose active weight grows with their arity and bit length.

## 3. Sharpness at W=3, including no-wrap and distinct rows

We make the planar construction of PR #16 more sparse in its active coordinates. Let G be the fan consisting of a path on k vertices 1,...,k and a hub 0 adjacent to each path vertex, embedded with triangle faces (0,i,i+1). Delete the hub column from the oriented primal incidence U and the exterior face column from the oriented dual incidence V. Every spoke has at most one nonzero primal and two nonzero dual entries; every rim edge has at most two nonzero primal and one nonzero dual entry. Thus the rows of A+=[U,V] and A-=[U,-V] have l1 norm <=3. There are n=2k-1 edges and n total columns.

By the cut/cycle identities in PR #16, the matrices are rationally invertible and have the same integer Gram. The reduced primal Laplacian is tridiagonal with off-diagonal entries -1, end diagonals 2 and interior diagonals 3. Its determinant is F_(2k), the even-index Fibonacci number (see [Bogdanowicz, Formulas for the Number of Spanning Trees in a Fan](https://www.m-hikari.com/ams/ams-password-2008/ams-password13-16-2008/bogdanowiczAMS13-16-2008.pdf), 2008; it also follows from the tridiagonal recurrence). Over any prime r dividing this determinant, its nullity is exactly one: each successive coordinate of a kernel vector is determined by the first coordinate, using the nonzero off-diagonal. As in PR #16, both A+ and A- then have modular corank one, and their combined span has full modular rank for odd r.

Pair each row with its negative and append a constant coordinate 1, giving blocks H,J of length 2n with identical moments through degree three. Give the first block tag 0 and the second block tag 1. The two populations are H at both tags versus H at tag 0 and J at tag 1. The tag contributes a new independent direction, since paired rows have the constant coordinate and there are rows with both tag values. Hence the two populations are rationally full rank n+2 and have respective modular ranks n+1 and n+2. Both have distinct rows (rational invertibility makes the paired base blocks distinct). They have equal exact cubic moments, the same known affine solution, and 4n rows. Their deficient direction extends the original null vector by zeros in the constant and tag coordinates; hence W<=3 on its active coordinates.

Choose k=11, r=199. Then F_22=17711=89*199 and 4n=84<199. This gives a degree-three counterexample with full rational rank, no entrywise count wrap, distinct rows, and W=3. Degree four is sufficient by Theorem 2; degree three is not. Six finite examples are verified (for k=11,13,17 and two primes for each). Thus Theorem 2 is sharp at W=3, even with these additional promises. This does not prove sharpness at W>=4. At W=2 the degree-two bound is sharp for unrestricted row families, for example the paired rows +(1,1), -(1,1) versus +(1,-1), -(1,-1), which have equal first moments but different kernels. That last example does not have full rational rank.

## 4. Verification, resources and open scope

Run `python -B verify.py` (Python 3.10+, standard library only). The [transcript](verification.txt) was produced with Python 3.13.5. The verifier checks 5,174 finite Cartesian sub-bound support ranks (exhaustive on the listed small grids, sampled for larger cases), 19 sharp unrestricted-weight product examples, 236 sign-character indicator polynomials on complete integer balls, 201 rational separation certificates at 14,133 points, and six no-wrap fan examples. All arithmetic uses exact integers, rational fractions or modular elimination. These checks are not proof of the asymptotic theorem; the proofs are supplied above. No curve implementation or implicit counting contractor is included.

The primary boundary question remains open for W>=4: can one obtain higher-order exact agreement, no-wrap counts, small active weight and modular rank mismatch at the same time? The support bound in section 1 and the degree cap in section 2 narrow the parameter ranges. They do not exclude other oracle models, key-dependent weights, extra arithmetic information or richer query interfaces. The full-rank distinction is abstract, and none of these results establish a new end-to-end complexity bound for the project's actual relation systems.
