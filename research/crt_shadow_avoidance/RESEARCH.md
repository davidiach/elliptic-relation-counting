# Carry-free CRT avoidance, exact convolution ranks, and common-field shadows

27 November 2026. Abstract finite-group and coding-theory research follow-up to PRs #22, #25, #28 and #30 of `davidiach/elliptic-relation-counting`. Reviewed #25 at `911799c4f2c7d9fde822590db9391c4b32440d98`, #28 at `e6677bbc48dac6ade5f63fd8cbd0e4ce94368e82`, and #30 at `60309a8e34c8dc7ba00e26d7acea82716d658481`. This note removes the independent carry-box overhead of the earlier base-field digit lift. It also gives two constructive exact algorithms, with checkable and complementary hypotheses, for the resulting mixed-modulus problem. A deterministic worst-case sublinear-in-arity algorithm on **all** allowed inputs is not established. No improved end-to-end exponent is claimed. All experiments use abstract finite groups.

## 1. Carry-free mixed-characteristic lifting

Let r be an odd prime, let A_1,...,A_m be nonempty explicitly listed subsets of F_r^2, let n_i bound their sizes, and let M=prod_i n_i. Assume fixed slack M<=r^(1+sigma), 0<=sigma<1. Put b=ceil(log_2 r). Take distinct small primes p_1,...,p_s with product R<=r. Write G_R=(Z/RZ)^2, identified with the product of the two copies of Z/p_jZ by the Chinese remainder theorem. For a plane point a, take canonical integer representatives and define rho(a) by reducing its two coordinates modulo each p_j. Decode every point of the product using CRT representatives in [0,R), which are also legal points of F_r^2. This decoding is a right inverse of rho, not a global inverse on the original plane.

Define the explicit wrap list

    C = {rho(-t r,-u r): 0<=t,u<m},

where rho on integer arguments means ordinary componentwise reduction (not prior reduction modulo r). Then

    rho(A_1+...+A_m mod r) subset rho(A_1)+...+rho(A_m)+C.       (1)

Indeed an unwrapped sum of m canonical coordinates equals its reduced coordinate plus t r with 0<=t<m. Integer addition is coordinatewise addition in CRT, so subtract the corresponding wrap term. This is exact inclusion, including when an original point lies outside the canonical R-square. Consequently any point absent from the right side decodes to an absent original plane point. A probe set not contained in the expanded support also decodes to a successful probe set. It is not necessary to identify the absent member for the abstract probe-list task. Unlike the base-q digit lift in #25, there is **no per-coordinate carry correction**.

### 1.1 Fully charged parameters

For all sufficiently large b there are enough primes in [b,4b] that their product exceeds 2^b (by the classical prime number theorem, in its equivalent theta(x) form). Enumerate them deterministically by trial division or a sieve, and multiply consecutive primes until the next would exceed r. Thus

    r/(4b) < R <= r,
    s=Theta(b/log b), max p_j<=4b.

All integer operations use polynomial work in b; a bounded initial range may be treated directly. The expanded product of list sizes is at most m^2 M, and the ambient size R^2 satisfies

    m^2 M/R^2 <= 16 b^2 m^2 r^(-(1-sigma)) = r^(-1+sigma+o(1))

for growing polylogarithmic m. Each encoded point has O(b) bits. There is only one extra list of size <=m^2, no loss of fixed exponent slack, and no need to list the whole mixed-modulus group. In the balanced regime log n_max=O(b/m+log b), the reduction has polynomial cost in the explicit input length and b. An avoidance or probe-list algorithm on this product of small cyclic groups would retain its sublinear-in-arity complexity after the reduction.

There is a tradeoff: this uses many distinct characteristics, with maximum prime O(b), not the common field of size poly(m) used in #25. In the polylogarithmic-arity regime the product of all distinct primes <=m^c has logarithm at most m^c log(m^c)=o(b) for each fixed c, so a squarefree CRT lift of the required size cannot stay entirely in that alphabet scale. This is **not** a new reduction to the uncolored sparse-linear-span problem over one small field; integer coefficients across different primes do not have that interpretation.

### 1.2 Reverse transfer with a supplied nearby prime

The same argument also runs in reverse for arbitrary composite R and a *supplied* prime p>=R of size at most R poly(log R). Represent each residue mod R canonically in [0,R), then embed the lists in F_p^2. Add the wrap list consisting of (tR,uR) mod p for

    -(m-1) <= t,u <= floor((p-1)/R).

If a prime-plane point x has its reduction mod R in the original sumset, then its canonical coordinate is the unwrapped sum plus one of these integer multiples of R. Thus it belongs to the corrected prime sumset. An absent point decodes by reduction mod R. The additional list has polynomial size, fixed slack persists whenever the original product has fixed slack, and the map is explicit. The statement assumes the nearby prime is supplied: no efficient deterministic construction of arbitrary huge primes is silently charged as polynomial. For the forward lift p=r is already supplied. These two one-sided inclusions identify a convenient alternative formulation, not a proof of general avoidance.

## 2. Exact tensor-train avoidance on product groups

The next result applies to arbitrary product cyclic groups G=prod_{j=1}^d Z/p_jZ, not necessarily with distinct or prime moduli. The p_j are the coordinate group sizes. To avoid confusion, coefficient-field primes below are written ell and can be entirely unrelated to the p_j. The input comprises explicit nonnegative integer tables f_i on G (sparse tables of multiplicities are allowed). Let M=prod_i sum_x f_i(x), assume 0<M<|G|, and put F=f_1*...*f_m for ordinary count convolution. An absent point is exactly a point where F=0.

### 2.1 The rank guard and algorithm

Choose any explicit ordering of coordinates and any binary parenthesization of the list convolutions. A tensor train of a table f over a field is a sequence of slices G_j(x_j), of sizes r_{j-1} by r_j, r_0=r_d=1, with

    f(x_1,...,x_d)=G_1(x_1)...G_d(x_d).

The minimal nonzero-tensor bond ranks are the ranks of the prefix-versus-suffix unfolding matrices. Use rank-one zero cores for the zero tensor. Let B be a supplied positive integer cap. Suppose every intermediate convolution tensor in the chosen parenthesization, including the leaves, has all these unfolding ranks <=B over the rationals (or is zero). This is a sufficient condition; our actual guard tests ranks modulo several small primes, and can succeed even when the rational condition is not met.

**Theorem.** Under this guard, one can uniformly and deterministically construct an absent point in bit time and space polynomial in the explicit table size, d, max p_j, log(M+1), and B. The procedure never constructs the full tensor or enumerates the full convolution support. Without the guard it returns honest failure if a measured bond rank exceeds B. No approximation or uncharged rank oracle is used.

Proof. Take distinct small coefficient primes ell_1,...,ell_t whose product exceeds M. It is sufficient to enumerate the first ceil(log_2(M+1))+1 primes; their magnitudes are polynomial in log(M+1), by the standard prime-counting bound. Work separately over F_ell. Every unfolding rank of an integer tensor after reduction modulo ell is at most its rational rank, because all minors that are zero integrally remain zero. Thus the rational guard guarantees the modular guard.

**Input construction.** A sparse table on n words has an exact train of bond dimension <=n: one path per word, use its multiplicity in the first core and indicators on the remaining cores. This is polynomial in n,d,max p_j. A prefix-trie construction with rank factorizations from right to left, as used by the verifier, avoids materializing many zero entries.

**Convolution.** Given trains for f and g, set

    H_j(z_j) = sum_{x_j+y_j=z_j mod p_j} G_j(x_j) tensorproduct K_j(y_j).

The scalar product of these slices, using the distributive law and the Kronecker product identity, expands to precisely sum_{x+y=z} f(x)g(y). This is exact over every coefficient field; it does not require roots of unity or compatible characteristics. Raw bond dimensions multiply. If both input ranks <=B, every raw bond is <=B^2. Forming the slices takes polynomial work in the sizes, including the local alphabet loops.

**Exact compression.** Use two sweeps of rank factorization, not numerical SVD. Left-to-right flatten each core at its outgoing bond into a matrix A, factor A=X Y with full-column-rank X and full-row-rank Y by Gaussian elimination, replace the core by X and absorb Y into the next core. After this sweep every prefix matrix has full column rank. Right-to-left flatten each core at its incoming bond and again factor A=X Y, replace it by Y and absorb X into the previous core. The processed suffixes now have full row rank. At each right sweep step the full tensor unfolding at that cut is L A R with L full column rank and R full row rank; hence its rank equals rank(A). Absorbing X preserves the left independence at the next relevant cut. This establishes exact minimal bond ranks. A zero factorization represents the zero tensor with the stated convention. The two sweeps have polynomial cost in the dimensions of the cores. Compress every leaf before combining, and every parent before using it again. Check the cap and abort if it is exceeded. No dense unfolding over the ambient group is constructed.

**Selecting an absent point.** The total mass at the root is M<|G|. Given a fixed prefix, the exact count of its completions is obtained by multiplying the fixed core slices and, on the remaining coordinates, the sums of slices. Do this modulo all coefficient primes and recover the unique integer in [0,M] by CRT. Suffix sums can be prepared once. For a prefix with count strictly smaller than its number of completions, scan the next coordinate and compute each child's exact count. Their sum is the parent count, and their completion counts sum to the parent completion count. At least one child remains strictly underfull. Continue along such a child. At the full word its count is a nonnegative integer smaller than one, hence zero. The total number of child counts is at most sum_j p_j. All field arithmetic is on small moduli, and final CRT reconstructions have O(log(M+1)) bits. This proves both correctness and bit complexity. QED.

The result uses ordinary tensor-train algebra, not a new fast numerical tensor decomposition theorem. The contribution is exact modular rank compression of count convolutions with a sufficient rational rank guard, coupled with exact missing-point selection and the carry-free reduction. It is different from #28's capped minimum-cost profiles along a tree of linear interfaces: the ranks here concern count tensors across coordinate cuts, and can be small even without a supplied field-linear direct sum of supports. All intermediate ranks, not only the final rank, must pass the guard. One may also marginalize to any selected subset of coordinates whose product exceeds M and run this procedure there; filling other coordinates arbitrarily preserves absence. This does not prove a cheap universal choice of subset or ordering.

### 2.2 Structured examples beyond low-order coordinate projection

Let the CRT coordinate pairs have alphabets p_j^2. Partition the prime indices into m groups of roughly equal size L (rounding by one is harmless). Choose a common auxiliary prime Q>=max p_j^2 and Q>=L, with Q<=2 max(p_j^2,L)+2. On each group take the evaluations of all polynomials of degree <floor(L/2) at its L distinct points in F_Q. Map each evaluation value surjectively onto the alphabet of the corresponding prime pair by reducing its canonical integer modulo p_j^2 and splitting it into two digits. Keep distinct words. Let the i-th list be this code-image supported on its own group, zero elsewhere. The product of supports across groups is exact; every projection onto at most floor(L/2) prime pairs is full, since interpolation prescribes arbitrary field evaluations on that many coordinates inside each group and the local maps are surjective. Yet the overall number of tuple choices is at most Q^(sum floor(L_i/2))=r^(1+o(1)) in the asymptotic CRT setup, and the ambient size is r^(2-o(1)). Decode the lists as genuine original plane lists and include the wrap list. The expanded support still has every indicated projection full and fixed positive volume slack.

In a coordinate order grouping the code blocks, the main count tensor is a product of the disjoint local indicator tensors; every bond rank is at most the largest local list size, 2^{O(b/m+log b)}. Intermediate products in any order among the main factors satisfy the same rank bound. The final convolution with the at most m^2 wrap vectors is a sum of at most m^2 coordinatewise translates, each preserving unfolding ranks, so its ranks are at most m^2 times this bound. The exact tensor-train algorithm applies with a polynomial-in-explicit-size cap, even though the support of the product need not fit in memory and there is no gap in any of the stated low-order CRT projections. This example is a mathematical list construction, not a claim about selected factor-base masks.

## 3. An exact common-field shadow algorithm

There is a second, complementary way to use a mixed-modulus representation. For now take explicit lists A_i in prod_j Z/p_jZ and an additional explicit correction list C of size L>=1. We treat C as translates, and the A_i as the main summands. Fix an auxiliary prime Q. Regard the canonical residue words of each A_i as integer vectors, and then reduce their coordinates modulo Q. Choose references a_i in each list, put a0=sum_i a_i mod Q, and let V be the F_Q-linear span of all reduced differences a-a_i. Let k=dim V. Thus every unwrapped tuple sum has its reduction modulo Q in a0+V. Compute this affine space directly by elimination; no hidden code coordinates or assumed Reed-Solomon structure are needed.

For coordinate j put

    U_j = sum_i max_{a in A_i} a_j,
    T_j = 1+floor(U_j/p_j),
    beta_j = min(1,T_j ceil(p_j/Q)/p_j).

The allowed residue-shadow set above a canonical residue z is

    D_j(z) = {(z+t p_j) mod Q : 0<=t<T_j, z+t p_j<=U_j}.

The wrap number is nonnegative because the unwrapped sum of canonical input integers is nonnegative. If a target z belongs to the corrected sumset, some c in C and v in a0+V must satisfy

    v_j in D_j((z_j-c_j) mod p_j) for every j.                 (2)

This is only an overapproximation; nonexistence of a compatible pair proves absence.

Choose an information coordinate set I of size k on which projection of V is invertible, via pivot selection on a basis matrix. Let B0=L prod_{j in I} T_j. If

    B0 prod_{j notin I} beta_j < 1,                             (3)

there is an exact deterministic missing-point algorithm costing polynomial in explicit input size, d, Q, max p_j, and B0. The bound B0 can be capped beforehand and (3) tested with exact integer comparisons.

Proof. Fix target coordinates on I arbitrarily, e.g. to zero. For each c in C enumerate all possible choices of v_I in the product of the corresponding D_j. Each such choice determines a unique vector v in a0+V. There are at most B0 candidates (c,v), with duplicates harmless. For any particular candidate and remaining coordinate j, at most T_j ceil(p_j/Q) choices of z_j make (2) hold: for each t, the canonical residue after translating by c_j has to occupy one congruence class modulo Q. Thus among all p_j choices, the average number of surviving candidates is at most beta_j times the current count. Try all p_j residues and choose one minimizing the survivor count. At the end the count is a nonnegative integer <= the left side of (3), hence zero. The resulting target is absent by (2). All trials, affine lifts and comparisons are explicit and exact. The proof also handles k=0, using a singleton affine point. QED.

This is a new checkable condition that is neither a low tensor-rank condition nor a low-degree polynomial separator test. Trying all small auxiliary primes costs a polynomial factor when the prime range is polynomial in b. One can also try different information sets, changing the carry bounds in B0 and (3); no universal guarantee for this optimization is made. The method works at large k if the carry bounds T_j are small: the expensive dimension in the enumeration is k log T rather than k log Q. Arbitrary dictionaries can have full affine rank modulo every auxiliary prime (for instance if the union of differences contains the ordinary coordinate vectors), so this is not a universal theorem under the cardinality promise.

## 4. An explicit all-order exponential-width obstruction that the shadow algorithm solves

We give a single family that sharply separates the two algorithms above. This also rules out a proof that merely chooses a better coordinate ordering for the exact count tensor whenever the CRT support is sparse.

Use the prime-plane/CRT setup in section 1 with p_j in [b,4b], s=Theta(b/log b), d=2s; each p_j occurs twice. Let growing m satisfy m=o(log b), in particular the original polylogarithmic arity m=Theta((log b)^1/4). Choose a prime Q between ceil((max p_j-1)/m)+1 and twice this number (Bertrand). For sufficiently large inputs,

    d<Q<=min p_j,
    m(Q-1)>=max p_j-1,
    C0=1+floor(m(Q-1)/min p_j)=O(1),
    log Q=(1-o(1))log b.

Let k be the least integer with Q^k>=r. Since R<=r and Q<p_j, k>=s; since log R=b-O(log b) and log Q~log b, k=(1+o(1))s<d. Choose d distinct nonzero evaluation parameters in F_Q. Split the k consecutive polynomial coefficients (powers 0,...,k-1) into m nonempty groups of sizes differing by at most one. List i contains the evaluations of all polynomials supported on its group. Interpret each evaluation as its canonical integer in [0,Q), hence a legal residue in the corresponding prime coordinate. These words are distinct within each list since d>=k. Their product size is

    M=Q^k in [r,Q r), n_max<=Q^ceil(k/m)=2^{O(b/m+log b)}.

The lists decode deterministically to ordinary plane inputs. Include the global wrap list C of section 1 (size L<=m^2, containing zero). For any fixed positive slack sigma, the volume guard L M<R^2 holds eventually. The explicit construction cost is polynomial in the list size and b. No assumption about elliptic-curve membership predicates or full original-plane group coverage is made.

### 4.1 Every small coordinate projection is full

Let a=min_i floor(k/m)=floor(k/m), the smallest coefficient-group size. On any set of at most a evaluation positions each group polynomial can take arbitrary values over F_Q, because consecutive powers at distinct nonzero parameters give a scaled Vandermonde matrix of full row rank. Thus each of these coordinates can independently prescribe m canonical integers from [0,Q-1]. Their sum can attain every integer from 0 to m(Q-1), which includes every canonical residue of every coordinate modulus. The main sumset projects surjectively onto **every** set of at most a coordinate positions, and the expanded sumset does also. In particular every projection onto at most floor(a/2) whole prime pairs is full. Here a tends to infinity and is Theta(b/(m log b)). This is not a claim about arbitrary non-coordinate linear maps over the prime plane.

### 4.2 Every coordinate order has an exponentially large bond rank

For any subset J of j coordinate positions, fix a target projection of an uncorrected sum. In each fixed coordinate there are at most C0 possible unwrapped integer totals, since the total is in [0,m(Q-1)] and its residue is fixed. Reduce those integers modulo Q. For any fixed resulting evaluations there are at most Q^max(0,k-j) possible coefficient vectors of a degree-<k polynomial, by the Vandermonde rank. The coefficient vector determines all m chosen tuples. Thus the number of tuples giving any prescribed projection is at most

    C0^j Q^max(0,k-j).                                        (4)

At any balanced cut with s coordinates on each side put B1=C0^s Q^(k-s). For the main count tensor the total mass is M, so each side has at least M/B1 different projected values. Each row and column of the main unfolding has at most B1 nonzero entries (indeed total mass <=B1). For the corrected tensor, sum the L translated main tensors: each row or column has at most L B1 nonzero entries. Because zero is an allowed correction, it still has at least M/B1 nonzero columns as an integer matrix.

Over the rationals, a basis of t rows has union of supports of size at most t L B1; all other rows are combinations of this basis. Thus the unfolding has rank at least M/(L B1^2). More importantly, the obstruction survives the *exact modular algorithm* of section 2. If its coefficient primes have product >L M, each nonzero integer column remains nonzero in at least one such field (some positive entry is <=L M and cannot be divisible by their whole product). With t0 primes, one field therefore has at least M/(B1 t0) nonzero columns. The row-support argument gives rank at least

    M/(L B1^2 t0) = Q^(2s-k)/(L t0 C0^(2s)).                  (5)

For t0 polynomial in b (as in the actual reconstruction algorithm), the logarithm of (5) is b-o(b): s log Q=b-o(b), k log Q=b+O(log b), and s log C0=o(b). Every coordinate order has such a balanced cut. Hence **every** ordering forces an unfolding rank >=2^{b-o(b)} over the rationals, and at least one coefficient prime in an exact polynomial-length reconstruction system does so as well. In particular the exact train algorithm cannot use cap 2^{O(b/m)} on this family. This is a lower bound on one representation width, not on all algorithms, not on approximate numerical tensor methods, and not on the original (non-expanded) plane sumset.

### 4.3 The same family has a polynomial-size common-field shadow certificate

Despite (5), the auxiliary Q-shadow of section 3 has exactly the k-dimensional polynomial-evaluation linear space: each list contains zero and spans the corresponding evaluation columns, and their union spans all coefficient directions. We have T_j<=C0=O(1), and beta_j<=C0(1/Q+1/min p_j)=O(m/b). The guard (3) holds eventually, because

    L C0^k [O(m/b)]^(d-k) = 2^{-b+o(b)} < 1.

The initial candidate count L C0^k is 2^{O(b/log b)}. Since m=o(log b), this is 2^{o(b/m)}. The list and affine-space setup, plus all candidate trials, fit polynomial time in the explicitly enumerated list length and b. Q need not be supplied: trying all primes <=4b and checking the exact guard (and the candidate cap) finds the indicated Q eventually. Therefore this explicit, balanced, fixed-slack class with every small-coordinate projection full and **exponential exact unfolding ranks in every ordering** has a certified missing-point constructor with polynomial explicit-input cost. The method uses a change of coefficient characteristic and a bounded carry *ambiguity* rather than approximating the count tensor. It does not imply that any arbitrary sparse mixed-modulus input has an affordable affine shadow.

## 5. Exact verification and scope

Run `OPENBLAS_NUM_THREADS=1 python -B check_crt_rank.py` with Python and NumPy, assertions enabled. NumPy int64 operations are used only for the small field and matrix sizes of these finite fixtures; theoretical operations use arbitrary-precision exact integers or ordinary finite-field arithmetic. `verification.txt` is the complete transcript. The initial larger mixed-modulus and product-code runs exceeded the runtime limit; the finished verifier uses a smaller mixed fixture and optimized rank-one core convolutions. Two completed runs matched byte-for-byte. There is no inferred asymptotic performance from the wall-clock time of the fixtures.

Verification covers 36 exact CRT plane instances, 854 represented-point inclusions, 6,811 safe absent decodings, and 12,780 reverse-lift checks. The tensor code tests 240 exact tensors, 3,600 point evaluations, 480 independent unfolding-rank comparisons, 80 convolutions and 31 certified missing points, including reconstruction beyond one coefficient prime and an explicit rank-cap failure. The mixed-modulus positive example uses r=211, R=210, four distinct prime moduli, product 7,632<44,100, every single-prime pair projection full, and peak bond 25. A separate 48-coordinate F5 product-code example has 12 explicit lists of 25 entries, ambient 5^48, choices 25^12, peak bond 25, and 240 child counts. All projections onto at most two coordinates are full.

The MDS/carry tests cover four finite mixed fixtures, 52 balanced cuts, 29 full low-coordinate projections, two certified shadow successes (including the separate 12-coordinate positive fixture), three failed sufficient guards retained, and a full-affine-rank negative control. An independent sparse elimination on ten distinct balanced unfoldings of one mixed example gives ranks 1,263 through 2,169 over F101, whereas its auxiliary affine rank is three and eight initial candidates suffice for shadow avoidance. This test uses the actual expanded count tensor and verifies absence against it. Four further large-bitlength integer-parameter guards pass; they use the numeric surrogate 2^b-1 to check the inequalities, not a claim that those numbers are prime or an enumeration of astronomical lists. At b=262,144 the conservative modular rank lower bound already has log2 at least 163,592. All asymptotic statements follow from the proofs rather than those finite guards.

### Literature and limits

CRT residue representation and the absence of interposition carries are classical (NIST DLMF §27.15, https://dlmf.nist.gov/27.15). Exact unfolding-rank tensor trains and core multiplication are classical tensor algebra (Oseledets, *Tensor-Train Decomposition*, SIAM J. Sci. Comput. 33, 2011, https://doi.org/10.1137/090752286). The present application uses finite-field rank factorizations and exact integer reconstruction; no general-purpose new tensor format is claimed. Reed-Solomon interpolation, affine information sets and averaging are classical. We have not established priority of the combined CRT/affine-shadow method or of the all-order width separation. For related limits of sparse linear factorizations see Dvir, Golovnev and Weinstein, *Static Data Structure Lower Bounds Imply Rigidity*, https://golovnev.org/papers/ds.pdf. Volk's *Improved Algorithms for the Remote Point Problem*, https://arxiv.org/abs/2609.25765, still does not supply a general high-rate remote-point algorithm in #25's regime. General circuit range-avoidance results do not give a uniform algorithm for the present arbitrary mixed-modulus dictionaries.

**Research status:** the carry-free mixed-characteristic reduction and both constructive algorithms are unconditional abstract finite-group results. The all-order width lower bound applies to their explicit representation class, not to all avoidance algorithms. The two exact algorithms have additional checkable hypotheses, and their disjunction is not proved universal. The original uniform sublinear-in-arity probe-generator problem remains open. No improved end-to-end exponent, independent proof-chain review, or complete literature-priority review is claimed.
