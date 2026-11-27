# Balanced incidence alphabets: dimension-free interpolation and Vandermonde moment trades

27 November 2026. Follow-up to PR #19 (head `322aebc8ef7e1adc68a4a18940d20b365ce88e5d`). The statements concern finite integer populations, including the balanced incidence alphabet of one- and two-copy relation blocks. They are informational theorems about complete exact coordinate moments, not efficient evaluation of those moments on an implicit relation family. No end-to-end ECDLP improvement or lower bound is claimed. The classical ingredients (simplex interpolation, Vandermonde determinants, cycle incidence, tensor products) are not claimed novel; independent proof and priority reviews are outstanding.

## 1. Balanced-block interpolation theorem

For n>=2 and w>=1 let

    B(n,w) = {a in Z^n : sum_i a_i=0, sum_i max(a_i,0)<=w}.

Equivalently a is the difference of two multiset incidence vectors, each of mass w (common filler occurrences may be added). Let T(n,w) be the smallest degree for which evaluation of rational polynomials of total degree <=T surjects onto all rational functions on B(n,w). Write B_unit for its ternary subset a_i in {-1,0,1}, and T_unit for the corresponding interpolation degree.

**Theorem 1.** For all n,w,

    T(n,1)<=3,       T(n,w)<=4w-2 (w>=2),
    T_unit(n,w)<=3w.

Proof. Order all points by decreasing support cardinality and, within equal cardinality, decreasing positive mass, with any tie order. We construct one raw polynomial per point that is nonzero there and vanishes at every point later in the order. The resulting triangular evaluation matrix is invertible over Q. Thus these polynomials span all functions on the alphabet, and back substitution yields every point indicator without increasing the degree.

For a nonzero target a let S=supp(a), s=|S|>=2, alpha_i=|a_i| on S, and k=mass(a). Define

    G_a(x)=prod_{i in S} [x_i prod_{j=1}^{alpha_i-1}(x_i^2-j^2)].

It vanishes unless |x_i|>=alpha_i on S, and its degree is 4k-s. Choose R to be the smaller by cardinality of the positive or negative coordinate sets of a. Put eps_i=sign(a_i), and multiply G_a by prod_{i in R}(alpha_i+eps_i*x_i). The resulting polynomial is nonzero at a. It vanishes at all points with smaller support cardinality, and all points with equal cardinality but different support. For the same support, if it is nonzero at a point of mass <=k, its magnitudes must equal alpha. Every chosen sign must agree with a: an opposite sign gives a zero factor. Since the selected side already accounts for mass k of its sign, balance forces every other sign. Thus a is the only point of its support and mass <=k where the polynomial is nonzero. Its degree is at most 4k-s+floor(s/2)=4k-ceil(s/2). For s>=3 this is <=4w-2. On B_unit, alpha_i=1 and s=2k, so the degree is <=3k<=3w, including s=2.

For the general integer alphabet at s=2, say S={i,j}, every point with exactly that support has x_j=-x_i and x_i in {-w,...,-1,1,...,w}. Use

    x_i*x_j * prod_{t=-w,...,w; t not in {0,a_i}} (x_i-t).

Its degree is 2w+1, at most 4w-2 for w>=2 and equal to 3 for w=1. It vanishes at all later points other than the target and is nonzero there. The zero point is last and uses the constant one. This proves the triangular construction. QED.

All raw polynomials involve at most 2w coordinates and have at most 2^{O(w)} monomials. The integer coefficients, diagonal values, and point evaluations have O(w log(w+1)) bits (up to absolute constant factors). Exact interpolation by rational triangular elimination is polynomial in the number of alphabet points, w, n, and input value bit lengths, by standard exact linear algebra/minor bounds. The number of points is at most binomial(n+w-1,w)^2. This is not a bound polynomial in log n, and the full population can be much larger than the explicit factor table.

**Consequences.** Complete exact moments through T determine every integer or rational row multiplicity: apply the linear functional defined by the moments to the interpolating point indicators. In particular, for any modular linear equation, interpolate its violation indicator (0 on legal, 1 on illegal), obtaining a nonnegative exact rational polynomial on the alphabet. Its sum detects violations of any nonnegative population independently of coefficient characteristic, count wrapping, or multiplicity size. For the homogeneous case on any centrally symmetric alphabet the polynomial can be evenized, so the degree can be rounded down to the greatest even integer <=T; this last statement concerns homogeneous separation, not general reconstruction.

## 2. Sharp and asymptotically sharp lower bounds from alternants

A finite signed population C=(c_a) defines the Laurent polynomial F(z)=sum_a c_a z^a. The coefficients of F(exp(t)) in formal characteristic-zero series are the exact moments divided by multi-index factorials. The smallest nonvanishing moment degree is therefore the order of vanishing at z=(1,...,1). A nonzero signed population annihilating every polynomial of degree <L proves T>=L. Its positive and negative parts have identical moments through L-1. This is the elementary moment/jet correspondence, as in the broader interpolation and Laurent multiplicity literature; see Fabian Gundlach, [Polynomials vanishing at lattice points in a convex set](https://arxiv.org/abs/2107.05353) (2021).

For k>=1, let N=2k+1 and w_k=k(k+1)/2. Consider the normalized Vandermonde

    F_k(z) = (prod_i z_i^(-k)) prod_{i<j}(z_j-z_i).

Its support consists exactly of the N! permutations of (-k,-k+1,...,k). Every coefficient is +1 or -1 according to permutation parity. Every row has coordinate sum zero and positive mass w_k. The first nonzero homogeneous Taylor term at the identity is the Vandermonde in the local coordinates, of degree

    L_k = N(N-1)/2 = k(2k+1) = 4w_k-k.

There is also a direct exact moment proof: for exponent vector e of length N, the moment is det[v_j^(e_i)] with v_j=-k,...,k. When sum e_i<L_k some e_i repeat, so this determinant is zero. For e=(0,1,...,N-1) it is the nonzero ordinary Vandermonde determinant. Thus all lower moments vanish and degree L_k is witnessed exactly. Both populations have N!/2 distinct rows. Counts can be below any prime chosen larger than N!/2, but the two populations are NOT asserted to have different modular rank.

**Corollary 2 (exact small multiplicities).** For n>=3, T(n,1)=3 and T(n,2)=6. For n>=5, T(n,3)=10. The lower bounds are F_1, F_1^2, and F_2, respectively, with unused coordinates set to zero. A Laurent product preserves the balance and adds the mass bounds; its vanishing orders add since the leading homogeneous terms multiply in an integral domain. For n=2 the alphabet is a line of 2w+1 points, so T(2,w)=2w. We do not claim T=10 for n=3 or 4 at w=3.

**Corollary 3 (all weights, sharp leading constant).** For w arbitrary, let k be maximal with w_k<=w and s=w-w_k, so 0<=s<=k. For n>=2k+1 the product F_k F_1^s (F_1 using three of the coordinates) has mass bound w and first nonzero moment degree

    4w-(k+s) >= 4w-2k.

Therefore, for the dimension-unrestricted balanced alphabet, T(w):=sup_n T(n,w) satisfies

    4w-O(sqrt(w)) <= T(w) <= 4w-2   (w>=2),
    lim_{w->infinity} T(w)/w = 4.

If distinct rows are desired for the lower bound, take the s triangle factors on disjoint triples of new coordinates (n>=2k+1+3s). The product then has all coefficients +/-1 and no cancellations; its positive and negative parts have equal size N!6^s/2. The leading-constant lower bound is thus compatible with distinct rows and counts below a sufficiently large prime. Again this is reconstruction, not modular-rank separation.

**Corollary 4 (exact signed-unit bound for every w).** Take w disjoint triangle factors F_1, each on three fresh coordinates. Their product has ternary balanced rows of mass w, coefficients +/-1, and first nonzero moment degree 3w. Hence T_unit(n,w)=3w for n>=3w. Both populations have 6^w/2 distinct rows. The upper and lower bounds hold for every w, not merely a finite list. The contrast with Corollary 3 shows that admitting repeated coordinate occurrences changes the dimension-unrestricted reconstruction leading constant from 3 to 4. It does not imply the same leading constants for the strictly weaker modular-rank task.

## 3. Product alphabets and relevance to one- and two-copy rows

**Theorem 5 (additivity).** For any finite rational point sets A_i in disjoint coordinates, let T_i be their least full interpolation degrees. The full interpolation degree of their Cartesian product is sum_i T_i. For the upper bound tensor the point-indicator polynomials. For the lower bound, in each factor choose a nonzero rational signed measure annihilating all degrees below T_i, by non-surjectivity of that evaluation map (if T_i=0 use any nonzero measure). Tensor the measures. Every monomial of degree <sum T_i has some local degree <T_i, hence its moment vanishes. The tensor measure is nonzero. Clear denominators to obtain integral multiplicities. QED.

For balanced blocks of local multiplicities one and two, with at least three indexed labels per block, exact reconstruction degree is 3m and 6m respectively (m blocks). With mixed w_i in {1,2,3}, and at least five labels for each w_i=3 block, it is exactly sum_i T(w_i). Deleting one reference coordinate per block does not change interpolation degree: the deleted coordinate is the negative sum of the retained coordinates, an invertible affine linear parametrization of the balanced subspace. The results consequently apply to the same reduced row coordinates used by the repository.

The previous generic signed-integer ball rank certificate used W<=2m or 4m, yielding homogeneous degree bounds 4m-2 or 8m-2 when r>W. The balanced block result gives exact full reconstruction with degree 3m or 6m, also a sufficient rank-separation bound; use the smaller applicable bound at small m. In the doubled-copy case this reduces the leading degree from 8m to 6m, and does not need the prime-size assumption because it uses rational interpolation of an arbitrary finite violation indicator. For centrally symmetric domains the homogeneous separator may be evenized as in section 1.

This comparison counts polynomial degree only. Explicit product interpolation or moment tables can still require a product of the local alphabet sizes. Neither the degree reduction nor tensor factorization proves that the actual implicit fixed-shift contraction and rank extraction fit a better asymptotic bit budget. In particular the full population and the relevant coefficient table can exceed the desired scale. No new claim for the prescribed seeds or actual complete group-law relations follows.

## 4. Sharp homogeneous rank separation at active weight four for every odd prime r>4

This is a different lower bound, genuinely about modular rank rather than reconstruction. It strengthens the previous W=4, r=7 repeated-magnitude example by using signed-unit balanced rows, arbitrary odd prime r>=5, distinct rows and full rational rank after reference deletion. It does NOT satisfy no-wrap.

For one block with labels i in Z/rZ, define cyclic edge rows v_i=e_(i+1)-e_i (indices modulo r) and populations C_+={v_i}, C_-={-v_i}. They consist of distinct disjoint rows. With labels u_i=i in F_r their scores are +1 and -1 respectively. Counts agree, first moments vanish by the telescoping edge sum, and second moments agree by sign reversal. The cubic moment x_0^2*x_1 distinguishes them (before reference deletion; one can instead use x_1^2*x_2 afterwards). Thus Delta=C_+-C_- has first nonzero moment degree three.

On two independent blocks form

    H=(C_+ x C_-) union (C_- x C_+),
    J=(C_+ x C_+) union (C_- x C_-).

Then J-H=Delta tensor Delta, so H and J have exact equal moments through degree five and differ at degree six. Each consists of 2r^2 distinct signed-unit rows, with positive mass one in each block and active l1 weight <=4. For z=(u,u), every H row is legal, every J row has score +/-2 and violates. The rank statement is as follows. In either block the v_i span the balanced subspace over every field; their differences span that subspace over Q because sum_i v_i=0 and r is nonzero in Q. In characteristic r the differences span exactly the kernel of the score on that subspace: the quotient is generated by the common image of the v_i, and that image has score one. Fixing the other block in the Cartesian pairs produces these differences in each block independently. Hence both H and J have full rational rank in the reference-deleted dimension 2(r-1), and each has modular corank one. In the quotient by the two block score kernels, their spans are the lines generated by (1,-1) and (1,1), whose intersection is zero for odd r. Thus their combined span is full.

For a literal comparison of a rank-deficient and a full-rank distinct population with equal moment tables, append two inactive coordinates: a binary tag and a constant one. Compare H at tags 0 and 1 versus H at tag 0 and J at tag 1. Tagging preserves moment equality through degree five; the difference is restricted to tag 1. The populations have 4r^2 distinct rows, same known affine solution (the constant coordinate), and full rational rank in dimension 2(r-1)+2=2r after reference deletion. They have respective modular ranks 2r-1 and 2r. Indeed central symmetry in each tagged class lets the tag and constant directions be isolated by adding opposite rows (2 is a unit). The original direction z extended by zeros remains the unique deficient line of the former population and is violated by every J-tagged row. The active weight is still <=4, and it attains four for r>=5. The earlier degree-six upper bound therefore is sharp in this signed-unit balanced abstract model for every odd prime r>4. The row count exceeds r; these examples do not establish no-wrap sharpness at W=4. Nor do they claim to be complete elliptic relation families.

## 5. Verification and limitations

Run `python -B verify.py` with Python 3.10+, standard library only, assertions required. [Recorded output](verification.txt) is from the exact arithmetic verifier. Ten finite balanced alphabets (including multiplicities through four and signed-unit rows) were enumerated; every raw separator has nonzero diagonal and vanishes at all later points, with 887,049 exact zero checks. The alternants for k=1,2,3 have their support and exact nonzero determinant witnesses checked; 2,825 zero moment/determinant checks are exhaustive through the lower degree for k=1,2 and bounded for k=3. The product of two triangles is tested at every degree <=5 and a degree-six witness. Seven cyclic rank fixtures through prime 101 use actual spanning subsets and exact finite-field elimination; full tagged populations and exact moments through degree six are enumerated for r=5,7,11. Rational full rank is certified by full rank modulo the separate prime 1009. These finite checks supplement, but do not replace, the proofs.

Unresolved: exact T(w) for general integer w>=4, rank sharpness for signed-unit W>=5, no-wrap rank sharpness beyond W=3, cheap evaluation of certificates for the actual implicit row family, and any uniform end-to-end exponent improvement. The alternant trades solve a broad reconstruction boundary and do not establish rank gaps merely by having matched high-order moments. An exploratory conjecture T(n,w)<=3w for all integer-balanced rows was false: the normalized five-coordinate Vandermonde at w=3 has first nonzero moment degree ten. This distinction between ternary and repeated-coordinate rows is essential. A preliminary literature scan found general lattice interpolation/Laurent jet results, but did not establish novelty for the precise balanced-block formulas above.
