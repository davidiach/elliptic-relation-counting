# Short relation systems, conditioned means and adaptive prime-adic rank extraction

Research follow-up, 27 November 2026. Parent: PR #3, commit 89211f7652076a4cde9d89bcb9f9e797750268fb. LZR = Liu, Zhang and Ren, ePrint 2026/3020. These results use the closed Boolean contraction premise (LZR Theorem 2.2). The uniform deterministic quarter-power time bound remains OPEN. These are AI-assisted proof candidates, not an independent referee report or a practical timing claim.

## 1. A stronger seed potential and shorter complete relations

Let G=<P> have odd prime order r, b=ceil(log2 p), and m=Theta((log b)^(1/4)). In LZR Stage I use the same independently parity-split field-affine masks on the nonidentity points, with n=r-1. Let K=2^20 m^2, q*=max(Km,ceil((Kr)^(1/m))), T=ceil(pq*/n), q=nT/p, and mu=q^m. Work at sufficiently large b that q*<=n and n>=K; the finite initial range may be handled by exhaustive group arithmetic. The exact constants are not practically optimized.

For each seed, let d(g) count one-point-per-block tuples of sum g, W=prod_i n_i= sum_g d(g), and C=sum_g d(g)^2, where n_i=|F_i|. Choose the potential

    Psi=(rC-W^2+(W-mu)^2)/mu^2
       =(rC-2mu W+mu^2)/mu^2.

Both parts of the first numerator are nonnegative by Cauchy-Schwarz. LZR Lemma 3.2's coincidence-set bound, without fixing the old threshold, gives E[W]=mu and E[C]<=mu+mu^2(1+1/q)^m/n. Thus

    E[Psi]<=r/mu+(r/n)(1+1/q)^m-1<=5/K.

Indeed r/mu<=1/K, r/n<=1+1/K, and (1+1/q)^m<=1+2/K. The same bit-prefix conditional-expectation proof as LZR now uses only the exact W and C totals, plus the valid-seed completion count. All are closed contractions on O(mb) edges. Choose a seed with Psi<=5/K. Then |W/mu-1|<=sqrt(5/K)<1/2, so with p0(g)=d(g)/W,

    chi^2 := r sum_g (p0(g)-1/r)^2 = rC/W^2-1 <= 20/K <= epsilon,
    epsilon=1/(4096m).

Eventually q*=(Kr)^(1/m)+O(1), q=q*+O(1), mu=Theta(Kr), and W=Theta(Kr). Each n_i<=2T; the product lower bound gives n_i>=W/(2T)^(m-1)=2^{b/m-O(m+log b)}. In particular min_i n_i grows faster than any fixed power of m. This uses Hasse and the specified arity. Seed selection, factor-base enumeration and doubled-tuple universal decomposition retain the LZR bit-time scale 2^{O(b/(log b)^(1/4))}. The base size N0=sum_i n_i is 2^{O(b/m+log b)}. No stronger packet interface is used.

**Short relations.** Use independent one-per-block tuples X,Y with displacement val(Y)-val(X)=s. Thus there are only 2m slots per pair. For a fixed displacement the number of ordered pairs is W^2 Z_s, Z_s=sum_g p0(g)p0(g+s). Cauchy-Schwarz on the centered distributions gives (1-epsilon)/r<=Z_s<=(1+epsilon)/r. The count is Theta(K^2 r), less than r^2 eventually. In contrast, the original doubled tuples at LZR's asymptotic density have Theta(r^3) pairs per shift. This reduces arity and exact row counts, not by itself the time exponent.

An explicit polynomial-size small-bias list S of scalar shifts of absolute nontrivial character bias <=1/4 can be constructed without DLP. For nbits=ceil(log2 r), select 16 nbits distinct small primes ell !=r, take L=8 max ell, and include j/ell modulo r for 0<=j<L. For any nonzero character exponent x, at most 2 nbits selected primes divide x(r-x), with x canonical. Otherwise v=x/ell modulo r satisfies 1/ell<=v/r<=1-1/ell; the normalized geometric sum is bounded by ell/(2L)<=1/16. Average bias <=3/16. Prime enumeration is polynomial in nbits. For background see Jalan and Moshkovitz, Near-Optimal Cayley Expanders for Abelian Groups, FSTTCS 2021 (https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.FSTTCS.2021.24).

Include tests with displacement [a]P and known label difference a for a in S, and tests with displacement A_j and label difference y_j for every indexed point, including the dedicated columns P,Q. Their genuine rows have right side a in the first case and right side zero after subtracting e_j in the second. There are O(N+poly(b)) tests.

**Short-tuple constant-gap theorem.** For every incorrect candidate y some test fails on at least 1/8 of its valid pairs. Proof: let u_g be the conditional distribution of label error Lambda_y(X)-log_P(g) in the fiber g, assigning any distribution to an empty fiber. Let ubar=E_g u_g for uniform g. The small-bias averaging operator has nontrivial eigenvalues of absolute value <=1/4, so the average uniform-g success of the known scalar shifts is <=1/4+(3/4)||ubar||_2^2 (Fourier orthogonality and ||u_g||_2<=1).

Call g good when 3/(4r)<=p0(g)<=5/(4r). At most a fraction delta0=16 epsilon<=1/256 is bad (using even the weaker epsilon<=1/4096). Conditioned on displacement s the mass at g is p0(g)p0(g+s)/Z_s, at least 1/(2r) when both g,g+s are good. Consequently actual failure >= (uniform-g failure - 2 delta0)/2. If all actual failures were <1/8, average uniform-g failure for known shifts would be <33/128. This yields ||ubar||_2^2>21/32, so one error e0 has ubar(e0)>21/32. For a wrong column j, the point-shift test fails whenever both errors equal e0: their true required difference is log_P A_j, not y_j. Under uniform g this has probability at least 2 ubar(e0)-1>5/16. The actual failure is then >(5/16-1/128)/2>1/8, a contradiction. Thus the entire short relation system with known right sides has full column rank. This is a completeness theorem, not a construction of a small spanning subset.

## 2. Eliminating unknown displacements by cross-multiplied conditional means

The parent branch proves a full-rank stack of fixed-source and conditioned-class matrices at the prescribed generator (NEIGHBORING_FIBERS.md, Propositions 3-4). Here is a distinct all-displacement certificate using short tuples, with no unknown displacement logarithms on the right side.

For a pin h=(slot, A) in a short pair let n_h(t) count pairs of displacement t satisfying the pin, and let v_h(t) be their integer sum of signed incidence rows a=c_Y-c_X, on the N0 factor-base columns. For any two pins h,k define the homogeneous row

    R_hk(t)= n_k(t) v_h(t) - n_h(t) v_k(t).

The true logarithms lambda annihilate it, because v_h(t) lambda = n_h(t) log_P t. Every row is computable at a specified t using ordinary counts with at most two pins. There is no arbitrary label lookup.

**Dense modular pin support.** For any one pin in either half, the distribution of the remaining displacement is the convolution of a reflected or unreflected complete short-tuple law with another probability distribution. Its chi-square distance to uniform is <=epsilon by the l2 convolution inequality. Thus at most epsilon r displacement values have zero integer count. The total of the counts is W^2/n_i. At most W^2/(n_i r) values have count >=r. Therefore the fraction of t at which a specified pin has zero residue modulo r is <=epsilon+W^2/(n_i r^2)<=2epsilon eventually. Put delta=2epsilon. More particularly the support of each *unpinned remainder* count has density >=1-delta; it is a translate of each pin count's support in the corresponding slot.

**Rank theorem.** Over F_r the kernel of the stack of all R_hk(t), t in G, consists exactly of the tables

    u_{i,A} = kappa log_P A + c_i.

It has rank N0-m-1. Restore the two dedicated columns and the explicit LZR doubled-decomposition and normalization rows; the complete homogeneous kernel is then zero (eventually r>2m-1).

Proof: if u annihilates all cross-multiplied rows, at each t with some nonzero pin mass the ratios v_h(t)u/n_h(t) agree for all nonzero masses. Call the common value f(t); elsewhere define f arbitrarily. Fix a positive slot i. For A in F_i, if the unpinned remainder mass n_i^rem(x) is nonzero modulo r, then

    f(x+A)= h_i(x)/n_i^rem(x)+u_{i,A},

where h_i(x) is the sum of the remaining signed labels at displacement x. This is just separating the pinned contribution. Hence for A,B in F_i,

    f(t+(A-B))-f(t)=u_{i,A}-u_{i,B}

on a set of starting t of density at least 1-delta. The short sumset has density at least 1-epsilon>1/2 by the energy bound. Every displacement s is consequently a sum of m differences A_i-B_i. Telescope any such path: its predicted increment holds on at least (1-m delta)r starting points, by a union bound on translated exceptional sets. Different representations of s have overlapping domains and thus the same increment L(s). For s,t and s+t the domains for the two-step comparison and the direct comparison overlap because 3m delta<1. Therefore L(s+t)=L(s)+L(t). As G has prime order, L(s)=kappa log_P(s), and every within-block label difference has this form (compare the one-step relation to the m-step one). This gives the displayed affine labels. Conversely they give v_h(t)u=n_h(t)kappa log_P(t), so annihilate every cross row. Dimension m+1 follows because the dense sumset ensures at least one nonconstant block. For the explicit rows, homogeneous consistency on every base point forces c_i=2 sum_j c_j, so all c_i=0 when r does not divide 2m-1. Consistency at P and its normalization then give kappa=0; consistency at Q fixes the last coordinate.

**Density of inconsistent conditional means.** For a nonaffine direction u, let bad t mean at least one cross row has nonzero dot product with u. If the bad fraction is eta, define f as above on good t. The one-step proof then holds on density >=1-delta-2eta (both endpoints good). The same telescoping proves affine labels if 3m(delta+2eta)<1. Thus every nonaffine u has bad fraction >=(1/(3m)-delta)/2>=1/(12m) for delta<=1/(6m).

This removes the missing unknown RHS in the parent's all-displacement average by an exact algebraic identity. It does NOT show that an affordable stack of aggregated evaluations preserves rank. Finding enough specific displacements with nonzero rows remains open. Density alone cannot justify deterministic selection. At a supplied t, the full family has polynomial size in N; summing over all r displacements explicitly exceeds the budget.

## 3. Adaptive prime-adic rank extraction

Let A be any integer row family with row keys of O(b) bits and an exact implicit suffix/prefix Gram oracle. Let W denote the row space over F_r of already retained genuine equations. Write U for a matrix of integer lifts of a basis for W^perp, d=dim W^perp, and G=A^T A over Z. r is odd. We do not require G to have unit determinant, and may use several affordable row families or public partitions.

If GU is nonzero modulo r, choose a basis column u and row index j with (Gu)_j nonzero. The sum over rows of a_j(a dot u) is nonzero modulo r. If U^T G U is nonzero modulo r^2, some basis column u or sum of two basis columns has u^T G u nonzero modulo r^2, by polarization and invertibility of 2. Use the sum over rows of (a dot u)^2. In either case bisect the row-key interval using exact prefix/suffix Gram queries. One of the two child sums remains nonzero at the same modulus. At a single key, the row has a dot u nonzero modulo r; the square criterion uses r^2, not r. Append its genuine equation. Repeat until full rank or every test fails.

This characterizes all quadratic lift tests: if GU=0 mod r and U^T G U=0 mod r^2, then for any lift v=Uc+rz of a vector in W^perp, v^T G v=0 mod r^2 by expansion. Conversely the single columns, pair sums and lifts u+r e_j imply both conditions if all such tests vanish. If the process stops, extend U to a square integral matrix L invertible mod r. In L^T G L the leading d-by-d block is divisible by r^2 and both adjacent off-diagonal blocks by r. Extract r from each of the first d rows and columns. Therefore

    ord_r(det G) >= 2d

for every tested family whose determinant is nonzero. The rank-deficiency bound is a certificate, not a universal success claim. A single family of r^2 repeated scalar rows [1] stalls even with full modular rank. More importantly the previous C7 fixed-shift example with blocks (1,4,5),(1,4,5) and target 6, where det G is divisible by 49, is rescued by this method even with one combined family. Partitioning explicit and shift rows also succeeds.

The earlier squared-minor mass extraction is another fast path: for s already selected rows R and suffix Gram G_q, coefficient [z^s] det(G_q+z R^T R) equals the sum of squares of full minors extending R. Compute it exactly by s-th finite differences divided by s!, before modular reduction. Differences of suffix masses select the next least remaining key by bisection. A nonzero initial determinant modulo r^(v+1) yields a basis with 2 ord_r(det B)<=v. The adaptive method can do more even when this criterion does not yield a unit determinant.

## 4. Modular support certificates by higher moments

Let M<r^k (k>=2) be the number of rows in a tested family, let d_i be the integer discrepancy of a candidate, and suppose its nonzero residues mod r occupy at most J distinct values. If any discrepancy is illegal then at least one of

    sum_i d_i^j mod r^k,       k<=j<k(J+1)

is nonzero. In fact the stronger tests mod r^j can be used separately. Neither the actual residues nor their multiplicities need be known.

Proof: over Z/r^k Z the polynomials x^k and (x-a)^k for the distinct nonzero residue representatives a are pairwise comaximal, since all their root differences are units. Polynomial CRT gives a representative Q of degree <k(J+1) that is zero mod x^k and one mod every (x-a)^k. It has no terms of degree <k. At each actual d_i, Q(d_i) is precisely the violation indicator modulo r^k. Its sum is nonzero because 0< number of violations <=M<r^k. Thus one of the displayed moments is nonzero. Prefix bisection maintains the nonzero moment residue until an illegal row is isolated. The modulus r^j is also safe because j>=k and a legal discrepancy has j-th power divisible by r^j.

For the short fixed-shift family use k=2 since M=Theta(K^2 r)<r^2 eventually. A cross-mean family for a fixed pair of pins has at most r rows, so k=2 also applies. A fixed J requires only a constant number of pinned slots/copies. There is no proof that all candidate discrepancies have bounded support. The common-modulus tests alone cannot give such a theorem from the count bound: if h|r-1, take r copies of each of the h nonzero h-th roots of unity mod r. All moments of degree 1 through h-1 vanish modulo r^2 although every row violates, and hr<r^2 for h<r. This is an abstract obstruction to common-modulus tests, not a counterexample to degree-dependent moduli or to actual curve masks.

## 5. Complexity, earlier contributions, and limitations

All ordinary fixed-shift Gram entries expand into at most (4m)^2 signed two-pin counts. A cross-mean product expands into counts of two independent short pairs sharing a displacement; its Gram entries into four such pairs, with specified incidence pins. Prefix restrictions on the canonical coordinates of the shared displacement are ordinary arithmetic constraints. Higher moment j uses O(j) independent pairs and at most j additional point pins (expand the linear discrepancy and apply candidate coefficients outside the graph). The graph has O(jmb) edges, degree O(b^(3/4)) after equality trees and subdivision, and precision O(jmb). Count exact nonnegative integer contributions before reduction. The number of expanded summands is polynomial in N,b for fixed j, and integer coefficients and moments have polynomial bit lengths. Therefore all fixed-order tests, O(b) prefix decisions per selection, at most N retained equations and dense elimination fit 2^{O(b/(log b)^(1/4))} time and space under the closed premise. No full row or displacement array is enumerated in these asymptotic arguments. The finite scripts do enumerate toy arrays only. At growing moment order the number of pins and size of the copied graph must both be charged; no fixed-order guarantee is asserted in that regime.

Two preceding local research rounds supplied additional candidate theorems, recorded here for continuity.

* Original LZR potential implies chi-square distance <1/7 for its one-per-block distribution p0: minimize q x^2+64(x-1)^2 where q=rC/W^2 and x=W/mu. The doubled distribution nu=p0*p0 consequently satisfies 6/(7r)<nu(g)<8/(7r) everywhere by Cauchy-Schwarz on the centered convolution. Fixed-shift conditioning gives density >=9/(14r) on the initial group sum. The small-bias argument in section 1 then gives the 1/8 gap already for the original blocks. Bounded phase score 1-cos(2 pi a(d-c)/r), a from the small-bias list, has some expected score >=3/32 for a wrong candidate. An additive O(1/m) deterministic approximation of every conditional phase score would suffice for prefix extraction, but arbitrary label lookup remains expensive.
* Modal certificate: for each block compute error points [y_A]P-A, its most frequent error point, and L_e=sum of nonmodal fractions. For a pair of doubled tuples, subtract the modal point from each signed occurrence and let B be the set of nonzero deviations. For subsets U of B define f(U) as the 0/1 violation after adding just those deviations to the modal baseline. Let P_t be the sum of Boolean finite differences of f at the empty set of order <=t. Telescoping in an ordering of B bounds |f(B)-P_t|<=2^t binom(|B|,t+1): each omitted remainder is a (t+1)-fold difference with 2^t positive and negative terms in [0,1]. Thus L_t=P_t-2^t binom(|B|,t+1)<=f(B) pointwise. Under the doubled-pair shift law each individual tuple marginal has density ratio <=4/3 relative to independent slots. Cauchy-Schwarz gives E[5^|B|]<=(4/3) exp(48 L_e), and E[2^t binom(|B|,t+1)]<=(1/3)exp(48 L_e)2^(-t). For the test with failure >=1/8, t=min(4m,ceil(70 L_e)+4) gives positive expected L_t, with exactness when t=4m. Its signed integer mass and prefixes are computable by enumerating at most t+1 pins. This is a deterministic separator at the target scale when L_e=O(1), not a general guarantee for candidate tables. Both this and the modular-support test are opportunistic certified fast paths.
* Geometric isolation: injective row key digit vectors u_i in {0,...,D-1}^s induce the convex hull of basis sums. Every exposed vertex is realized by a unique matroid basis: perturb its exposing direction to make all individual row weights distinct, retaining the exposed vertex, and use matroid greedy uniqueness. Independent integer direction coordinates in [1,R], R>=2sN(D-1), isolate a vertex with probability >=1/2 by lower-envelope breakpoints (at most N(D-1) for each variable when the others are fixed). At digit width Theta(b/m), the resulting univariate determinant has degree <=NsR(D-1)=2^{O(b/m+log b)}, so a successful supplied direction can be tested by all degree+1 field points. Direction selection still needs O(b) random bits and no uniform deterministic list has been proved.
* Generic black-box rank condensers cannot simply be imported: a family of h linear maps with t output coordinates preserving every subspace of an ambient M-row space needs ht>=M (common kernel). Chevalley-Warning likewise gives a nonzero vector whose h arbitrary diagonal quadratic forms over F_r simultaneously vanish if M>2h; fix N-1 coordinate rows to extend the obstruction to full-column-rank arrays. These do not establish a barrier for actual elliptic rows. Low-order integer moments alone also fail without a support or concentration promise: interpolate on interlaced integer nodes jr and jr+(r-1)/2, with alternating Lagrange barycentric weights, to obtain two positive normalized distributions agreeing through degree 2k-1, one on legal multiples of r and the other off them. This is a probability-measure obstruction, not a statement about the integer multiplicities of the actual masks.

## 6. Reproducibility and next research

The accompanying finite tests check conditional-mean rank and density on small groups (including an actual prime-order curve), short-tuple robust gap, modular moment support, and the adaptive rank extractor. Use Python 3.10+ and SymPy for the determinant comparison. The scripts and recorded transcripts are in this directory. They are tests of the algebra, not an implementation of LZR's contractor or the stronger seed search at its enormous asymptotic threshold. The local earlier rounds' full notes are retained in the original research packages; this note records their principal proofs and limitations rather than reproducing their entire transcripts.

Highest priority: establish a uniform noncancellation/rank-accumulation theorem on the full cross-mean family, or on the parent's conditioned fixed-generator stack, using their exact convolution structure. A simple sufficient condition for success of our adaptive method is a nonzero permitted modular moment for every remaining candidate/annihilator. The density theorem does not establish this. Neither the old determinant condition nor the improved conditional rank theorem is a complete faster deterministic algorithm. Independent proof review, full source/Lean verification and literature priority checks remain outstanding.

References: LZR https://eprint.iacr.org/2026/3020, especially Theorems 2.2, 3.1, Lemma 3.2 and Proposition 4.1; PR #3 https://github.com/davidiach/elliptic-relation-counting/pull/3 (MAIN_RESULT.md, DETERMINISTIC_POINTS.md, NEIGHBORING_FIBERS.md at the parent above); Gurjar–Thierauf–Vishnoi, Isolating a Vertex via Lattices, https://eccc.weizmann.ac.il/report/2017/127/. No claim that the classical ingredients themselves are new.
