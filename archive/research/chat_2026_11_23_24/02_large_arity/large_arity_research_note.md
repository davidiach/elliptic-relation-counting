# Large-arity relation search: short ideal projectors, shifted digit blocks, and a taper barrier

23 November 2026. Research note, not an end-to-end speedup or a new asymptotic ECDLP bound. The preceding 8mb compiler gives a candidate numerical coefficient conditional on its packet theorem. The results below remove several auxiliary restrictions and isolate a weaker average-case evaluator target. No new bound below b/(log b)^(1/4) is proved.

## 1. Short ideal projector for arbitrary representatives

Use the preceding polynomial ring R=Z[X]/(X^n-lambda), lambda in {1,2}, with phi:R -> F_p given by X->gamma. Let Lambda=ker(phi). The preceding diagonal-dominant principal ideal construction supplies M in Lambda, with multiplication matrix A_M, D=|det A_M|=pc < p^2, and M0-2 sum_{i>0}|Mi|=Theta(M0). Let e=c*(c^{-1} mod p). In R/(M) the scalar e is the idempotent of the p-primary component. Set J=(M,e) in R. Since [R:(M)]=pc, J/(M) has size p, [R:J]=c, and J intersection Lambda=(M).

Minkowski's cube theorem supplies a nonzero vector in J with infinity norm <=2 c^(1/n); an exact Euclidean SVP call supplies E with infinity norm <=2 sqrt(n)c^(1/n). A nonzero vector in (M) has infinity norm at least M0-2 sum_{i>0}|Mi|, by strict diagonal dominance. As c=2^{o(b)} and M0=2^{(1+o(1))b/n}, E cannot belong to (M). Hence phi(E)!=0. Also E is zero in the complementary c-primary component of R/(M). Consequently, for every F in R,

    phi(F)=0 iff EF belongs to (M) iff there exists Q in R with EF=M Q.

The reverse implication follows by applying phi. Q is unique since A_M is invertible over Q. This proves that arbitrary short representatives can be used; no canonical or designated p-component restriction is required. No coset enumeration is required at inference.

For the balanced composite-degree n=q1 q2=Theta(b/(log b)^2) from the preceding construction, log D <= b+O(b/log b), so log ||E||_infty = O(log b). (The cyclotomic degree deficit O(sqrt(n)) is negligible here.) Thus E is a public polynomial with polynomially bounded coefficients. Multiplying a quartic residual by E adds only O(log b) to the coefficient width, and diagonal dominance bounds Q. A ring modulus of width 4b/d+O(log b) still suffices for EF-MQ=0, where d is the relevant irreducible degree. This is the same leading precision as the old quartic compiler. The short-projector preprocessing is 2^{O(n log n)} by exact SVP, plus polynomial-bit HNF and matrix operations. The test uses LLL for finite examples only.

The graph can therefore test leaf shifts and quartic equations with arbitrary bounded coefficient words. Its leading 8mb budget, subject to the inherited graph-to-contraction transfer, is not increased by this correction.

## 2. Entropy partition and an all-shift sampler

Take a dyadic rectangular coefficient box of total bit length A that contains a representative of every field element. With signed offsets, the evaluation of a word is an affine function of its A bits over F_p. Sort the bits by first-transform owner and partition them into m nearly equal consecutive blocks. For each block j let S_j be the image of its free bits with the other bits zero (after subtracting the fixed affine offset). Then S_1+...+S_m=F_p. Thus

    product_j |S_j| >= p, and max_j |S_j| >= p^(1/m).

Enumerate the m blocks once, choose a largest image S, and reuse it in each slot. Write w=2^ceil(A/m); |S|<=w and log w=b/m+O((A-b)/m+1). Independently choose a_i uniform in F_p. The slot accepts exactly the abscissas in a_i+S. To encode it, fix a short public representative of a_i and allow the chosen block's free bits as additive increments, with all other coefficients fixed. The sum of two short coefficient vectors is short, and the multiplier E removes any p-component restriction. The sampler retains both signs. Enumerate and deduplicate the accepted field abscissas and their rational points. All field elements have exactly the same marginal rho=|S|/p.

The support bit block is contiguous in owner order; hence at most O(1/m+1/V) of V Fourier owners contain variable digits (balanced owner sizes, up to a constant factor). There is no affine pilot or hash output array. If the selected bits have s runs after splitting at coefficient boundaries, |S-S|<=2^s w: each run contributes at most twice its number of assignments to the difference set, and phi is additive. In particular s<=n+O(1) for the owner-sorted block. The general bound |S-S|<=|S|^2 also holds. Put R_S=|S-S|/|S|; then log R_S<=min(log|S|, s+(A-b)/m+O(1)).

## 3. Average target coverage, no affine pilot

Let G=E(F_p), N=|G|. For each slot F_i, define A_i(D)=sum_P Pr(P,P+D in F_i). The identity is excluded from slots. Then A_i(O)=(N-1)rho, sum_D A_i(D)=E|F_i|^2<=4|S|^2. For D!=O, a retained pair has x(P+D)-x(P) in S-S. The rational function x(P+D)-x(P) has two double poles, so at most four preimages of each value. Joint retention is at most rho. Therefore A_i(D)<=4 rho |S-S|.

For uniform T in G, let H(T) count all rational ordered signed relations. Pair counting gives exactly

    E_{a,T} H(T)^2 = (1/N) sum_{D_1+...+D_m=O} product_i A_i(D_i).

Separate the all-zero differences; for every other tuple of differences, choose an index with D_i!=O, bound its A_i pointwise and sum freely over all other indices. With lambda0=N^(m-1)rho^m,

    E H^2 <= lambda0 + m 4^m (p/N)^(2m-1) R_S lambda0^2.

For every nonzero fixed target T, the number of rational affine tuples with sum T and all intermediate prefix sums nonidentity is >=(1-O(m/N))N^(m-1), by forcing a final point and a union bound on bad prefix identities. Their expected retained count X is >=(1-O(m/N))lambda0. Every such tuple supplies a finite abscissa chain, including genuine doubling. If target T is uniform in the nonzero elements of an odd prime subgroup of order r and polynomially bounded cofactor h=N/r, its second moment is bounded by N/(r-1) times the full-group second moment. Paley-Zygmund gives

    Pr(X>0) >= (1-O(m/N))^2 / [(N/(r-1)) (lambda0^-1 + m 4^m (p/N)^(2m-1) R_S)].

Since |S|^m>=p and N=p+O(sqrt(p)), lambda0>=(N/p)^(m-1)=1-o(1) for m=o(sqrt(p)). Thus the logarithm of inverse coverage is O(m+log R_S+log h+log b). This is an average-over-seeds-and-targets statement, not a targetwise guarantee.

## 4. Mixed-twist witnesses and verification instead of torsion pinning

In E(F_{p^2}) define G_+=ker(Frob-1) and G_-=ker(Frob+1), of orders N_+=N and N_-=2p+2-N. Every finite F_p abscissa has lifts in one of the two groups, and a lift lying in both is rational two-torsion. The intersection has size at most 4. A chain of S_3 zeroes can be backtracked from a chosen target T!=O to obtain coherent lifts with sum T: S_3(u,v,w)=0 gives a pair with finite sum having x=w, and simultaneous negation chooses the prescribed sign of that sum. The same polynomial handles tangency. At intermediate two-torsion the prescribed sign is unique.

If 2T!=O, a chain whose recovered leaves are not all rational must have a nonempty proper subset of non-rational (anti-rational) leaves. The sum of that subset belongs to G_+ intersection G_- because the other leaves and T are rational. For a fixed proper subset I of size k and a fixed shared torsion element U, the choices are bounded by N_-^(k-1) N_+^(m-k-1): choose all but one point in each side, then their required sums force the last points. The joint retention probability for their m abscissas is rho^m, by independent shifted masks. Thus, uniformly in T with 2T!=O,

    Pr(any mixed-twist leaf tuple with sum T)
       <= 4 (2^m-2) max(N_+,N_-)^(m-2) rho^m.

In particular this is <=lambda0 * 2^{m-b+O(1)+o(1)} when m=o(sqrt(p)), by Hasse's bound. Also lambda0<=2^{A-b+O(m)+o(1)}. In the polylogarithmic and square-root-bitlength tunings below the mixed-twist probability is 2^{-b+o(b)}. This is negligible compared with the average positive probability. Run the relaxed abscissa-chain finder; if a returned chain does not yield rational leaves and a verified target sum, reject it. No exact torsion pinning or inclusion-exclusion is needed for verified witness search. This does not give an exact signed-count oracle.

## 5. One-pilot size bias: only average evaluator performance is required

For fixed seed a, let kappa(a) be the probability that one query on a uniform nonzero subgroup target returns a verified relation. Let c(a) be its expected bit cost (the cost can correlate with success). Let bar{kappa}=E_a kappa(a)>0 and bar{c}=E_a c(a). All queries use fresh independent targets and internal randomness conditional on a. A fixed factor-base enumeration has B<=2m w labeled columns, retaining labels through cofactor projection.

Draw a fresh seed and run one pilot target; repeat independently until successful. The seed accepted by this procedure has distribution proportional to kappa(a) times the original seed distribution. The pilot's expected total query work is bar{c}/bar{kappa}. If setup is required per trial, add its average cost to bar{c}. Conditional on the accepted seed, one further verified row costs c(a)/kappa(a) in expectation, by the renewal identity. Averaging over the size-biased seed distribution gives at most bar{c}/bar{kappa} per row. This remains valid for zero-success seeds, which are never selected. Obtain B+1 rows, allowing repetitions, and solve for a nonzero linear dependence among factor-base coefficient vectors.

Generate each query as T=[u]P+[v]Q with independent uniform u,v modulo r, excluding T=O. Conditional on T, v is still uniform. After projecting factor-base points by [h], the relation rows give equations in their unknown logarithms relative to [h]P. A nonzero dependence vector chosen using only the coefficient rows (not u,v) yields sum_i c_i u_i + xi sum_i c_i v_i=0. The second sum is uniform, so is nonzero with probability 1-1/r. The resulting logarithm is verified, otherwise repeat. Ordinary Gaussian elimination costs polynomial(B,b) (e.g. O(B^3) field operations); retaining the inherited sparse solve yields O(m B^2) for the numerical comparison.

Therefore, up to polynomial factors in b,m and factor-base setup, the query cost is O(B bar{c}/bar{kappa}). No uniform fixed-target success condition or uniform adaptive-transcript parity oracle is required; only verified witness return, independence of the hidden random representation coefficients, and the charged average ratio are needed. This is weaker than the original fixed-target parity-evaluator target. It is not itself an evaluator.

## 6. Square-root-bitlength preprocessing and graph-size target

To keep all non-evaluator costs within 2^{O(sqrt(b) log b)}, use n an odd prime of order sqrt(b). If n divides p-1 choose gamma of exact order n, lambda=1 and d=n-1; otherwise choose gamma=2^(n^{-1} mod (p-1)), lambda=2 and d=n. The cyclotomic polynomial is 1+X+...+X^(n-1), and the preceding principal-ideal proof applies with d=n or n-1. Exact SVP and the short E construction cost 2^{O(n log n)}. The determinant/coefficient estimates give A=b+O(n log n+b/n)=b+O(sqrt(b) log b), log ||E||=O(log b), and W=4b/d+O(log b).

Prime transform length need not be factored into balanced dimensions. A length-n transform of beta^j A_j uses Bluestein's identity omega^(jk)=omega^(j^2/2) omega^(k^2/2) omega^(-(j-k)^2/2), exponents modulo odd n. Implement the resulting convolution with a public chirp by a padded power-of-two length L_F>=2n-1, using standard balanced two-stage forward and inverse transforms. Public twist/chirp factors are local; the output owners check EF=M Q at each actual frequency k. Choose an auxiliary q=-1 mod 2nL_F (with the inherited effective prime-progression bound), and use its unramified quadratic extension. The base ring has an n-th root of lambda since gcd(n,q-1)=1; its quadratic extension has primitive roots of the required lengths. This construction only claims O(m b) Boolean graph edges, not the 8mb numerical constant. It has polytime local entries and balanced degree O(sqrt(n)W) per constant-size bundle; use constant duplication/equality relays for shared state arrays if needed. This is a qualitative graph bound, not a proof that the specialized finite-packet hypothesis applies at this degree.

At arity m=Theta(sqrt(b)), each slot has log w=O(sqrt(b)), log R_S<=log w=O(sqrt(b)), and log bar{kappa}^{-1}=O(sqrt(b)+log h+log b) for an exact relaxed evaluator, after verification. Preprocessing, factor-base enumeration, ordinary Gaussian elimination, and retries are 2^{O(sqrt(b) log b)}. Thus a verified witness finder for the shifted-slot law with average query cost / average success probability <=2^{O(sqrt(b) polylog(b))} would suffice for the desired scale. All remaining difficulty resides in the evaluator. The generic published Boolean contractor, applied to O(mb)=O(b^{3/2}) edges, does NOT achieve this. No stronger asymptotic upper bound is claimed.

## 7. A barrier to optimizing only the valuation-symmetric Hermite taper

The public quantitative Boolean-contraction proof uses a packet of width P=2^a, deformed polynomial degree at most M(2P-1), and nodes t in [-P,P-1]. Give every nonzero node of 2-adic valuation i the same nonnegative multiplicity q_i, with q_i=0 for i>=u and a>=2u. To reconstruct every polynomial up to the generic degree, the degree condition is sum_i 2^{-i}q_i> M(2-1/P). Define the normalized negative valuation of its Hermite prefactor by

    ell_i=sum_{j=1}^{a-i} (j/2^j)(q_i-q_{i+j}),
    H=max(0,max_{i:q_i>0}ell_i).

The identity v_2(E_t(0))=-(P/2^i)ell_i is the public proof's exact valuation formula. For a>=2u with u->infinity, H >=(6-o(1)) M/u. Proof: inactive levels have ell_i<=0. Since sum_{j=1}^{a-i}j/2^j=2-o(1) uniformly over active levels, backward induction first yields q_i=O(H(u-i)). Extending the sum to j=infinity changes ell_i by o(H), uniformly. Therefore 2q_i-sum_{j>=1}(j/2^j)q_{i+j}<=(1+o(1))H. Define R_0=0 and R_k by 2R_k-sum_{j=1}^{k-1}(j/2^j)R_{k-j}=1. Backward induction gives q_i<=(1+o(1))H R_{u-i}. Its generating function is z(1-z/2)^2 /[(1-z)^2(2-z/2)], so R_k=k/6+O(1). Hence sum_i 2^{-i}q_i<=(1+o(1))H (u/3+O(1)), proving the bound. Conversely q_i=h(u-i)_+ with h~M/u (round upward) meets the degree condition and has ell_i<=6h, making this asymptotically sharp for that loss metric.

This is a barrier to improving the largest normalized prefactor valuation merely by adjusting the symmetric multiplicities, not a lower bound on contraction time, not a bound for nonsymmetric or multivariate interpolations, and not an assertion that the prefactor loss is necessarily incurred by every algorithm. In the current generic proof u=O(log M) and the partition enumeration combines the normalized precision loss with the mode density. Optimizing this taper alone cannot change that accounting's square-root logarithmic gain. The proposed complex cyclic deformation was investigated, but a proof controlling both its interpolation conditioning and crossing rank is missing.

## Evidence and scope

`python short_projector_audit.py` (requires sympy and polynomial_ring_verifier.py) exercises the arbitrary-representative projector on 96- and 128-bit toy lattice fixtures (both binomial branches), checks 120 valid and 120 invalid quartic residuals, and independently compares the energy identity with direct enumeration over all translations and targets on three small curves. It also checks the group-theoretic mixed-twist bound on all 54 eligible targets in those fixtures and the entropy partition on seven toy partitions. These are component regressions only.

Sources: the uploaded v3 paper (especially Sections 2, 4, 8), the preceding `abscissa_chain_note.md`, the repository's fixed-bit supplement, and the public Appendix D quantitative Boolean tensor-contraction theorem of OpenAI's October 25, 2026 manuscript. No novelty or independent peer-review claim is made.
