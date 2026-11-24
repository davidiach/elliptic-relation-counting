# Exterior-power packet deformation and a support-sensitive contractor

24 November 2026. Continuation of `hierarchical_shared_mode_note.md`, `contractor_degree_note.md`, and `seed_local_slice_rank_note.md`.

**Result and scope.** We construct an exact squared-minor deformation of Boolean equality, prove an exact common-factor formula and a uniform dyadic normalization bound, and obtain a first-level contraction algorithm sensitive to Hamming-weight support instead of multilinear degree. Subject to the public quantitative Boolean contractor, it has log time O(M/(log M)^(2/3)) for a class of sparse-support networks whose multilinear degree can be Theta(M). A strengthened slice-rank bound rules out using this parameter saving alone for exact, seed-local ordinary elliptic-relation parity, even with polynomial wire inflation, public bit flips, or a subexponential list of linearly combined networks. The general-prime ECDLP exponent and previous numerical coefficient candidate are **not improved or independently verified here**. The finite checks do not benchmark the asymptotic recursion. No novelty claim against all existing work is made.

The public proof used for packing, geometric crossing expansions, tapered Hermite recovery, and the subsequent ordinary child recursion is *The Exponential Time Hypothesis is False*, 25 October 2026, at `https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex`. The earlier independent numerical finite-packet premise is not used to infer a new leading coefficient. We retain the public contractor's assumptions; the new identities below are direct algebraic statements.

## 1. Squared minors replace affine Boolean selectors

Use packet width L=2^a, grid S={-L,...,L-1}, C(T)=product_(s in S)(T-s), c(T)=C(T)/C'(0), and W_xy(T)=c(T)/(T+x-y). The latter is a polynomial. W(0)=I. The packet row and column assignments are subsets Y and X of the L slots; dummy slots are forced absent. Define the new kernel to be zero if |X|!=|Y| and det(W[Y,X])^2 otherwise (the empty determinant is one). At T=0 this is exactly 1_{X=Y}. Thus one may partition the original contraction by a common occupancy k per packet, and use the product of squared-minor kernels for the corresponding deformed contraction. The tensor entries are lifted canonically at the input precision, as in the public proof.

Fix k, |X|=|Y|=k, and put n_s=#(x in Y,y in X:y-x=s). Cauchy's determinant identity gives, up to a fixed sign for the chosen orderings,

    det W[Y,X] = (Delta(Y) Delta(X) / C'(0)^k)
                   product_(s in S) (T-s)^(k-n_s),

where Delta denotes the Vandermonde product. Since n_s<=min(k,L-|s|), every minor is divisible by

    G_k(T)=product_(s in S) (T-s)^e_s, e_s=max(0,k-L+|s|).

The bound is exact for each root s: choose min(k,L-|s|) distinct matching pairs of positions separated by s, completing X,Y to size k if necessary. Thus this is the monic common factor over Q for the full packet family. Since sum_s e_s=k^2 and sum_s n_s=k^2, the degree after dividing the minor by G_k is exactly 2k(L-k). Also G_k(0)!=0. The normalized squared minor det(W[Y,X])^2 (G_k(0)/G_k(T))^2 has degree 4k(L-k) and still has value 1_{X=Y} at zero. For partial packets the same factor is valid, though it need not be the entire gcd. For a cardinality profile (k_j) the total degree of the product is at most 4 sum_j k_j(L-k_j). In particular, for a near-fully-packed network of M edges its largest possible degree bound is (1+o(1))ML instead of (2+o(1))ML. Near-full packing here is the regime where the public variance/padding errors are o(M), not an assumption about every finite packetization. If sum_j k_j<=H, the degree is at most 4HL, regardless of padding.

## 2. Exact dyadic cost of dividing the common factor

Write alpha=2L. At a grid node t and any unit z in the unramified dyadic ring, define V_t=v2(G_k(t+alpha z)). It is independent of z: the valuation of t+alpha z-s is v2(t-s) for s!=t, and a+1 for s=t. For 1<=j<=a+1 set

    A_j(t)=sum_(s=t mod 2^j) e_s.

Then V_t=sum_j A_j(t). In the cyclic group of order 2L, the function e_s is the circular autocorrelation of the length-k interval {0,...,k-1}, shifted by L. Consequently if q=2^j divides L, and b_x counts that interval in residue x modulo q,

    A_j(t)=sum_(x mod q) b_x b_(x+t).

By Cauchy-Schwarz this is at most sum b_x^2=A_j(0). At the final level j=a+1, A_j(t)=e_t<=k while A_j(0)=0. Since G_k(0) has valuation sum_j A_j(0) (e_0=0), we have proved

    V_t-v2(G_k(0)) <= k.                                      (1)

In particular the additional loss is at most 2k for a squared minor. This is not just a bound on point values. As a polynomial in z,

    G_k(t+alpha z)=2^V_t z^e_t U_t(z),

where U_t has dyadically integral coefficients and unit constant coefficient. W(t+alpha z) has integral polynomial entries by the public Cauchy integrality lemma. The squared minor is divisible by z^(2e_t), so removing this factor preserves integral coefficients. Dividing by U_t(z)^2 uses an integral formal Taylor series. Equation (1) proves that multiplying the normalized squared-minor polynomial by 2^(2k) makes every local Taylor coefficient integral. For several packets write K=sum k_j and

    F_norm(T)=2^(2K) F(T) product_j (G_(k_j)(0)/G_(k_j)(T))^2,

with F the unnormalized squared-minor contraction for that profile. This is a polynomial, locally integral everywhere on the grid, and F_norm(0)=2^(2K) times the desired profile contraction. Read its value modulo 2^(P0+2K) and divide by 2^(2K). To evaluate it near a node, use residue-distinct unit z values in a sufficiently large unramified extension. There are enough such values for full polynomial interpolation, with unit differences. The public scalar multiplier 2^(2K) product_j (G(0)/G(t+alpha z))^2 is integral by (1), so evaluation of the unnormalized contraction at the requested precision suffices. No inversion of a nonunit is performed without accounting for its known valuation. All public products and precisions have polynomial size.

## 3. Bordered determinants and the sign-free group contraction

Use ordinary flat packetization, so at any retained node a packet has at most two group labels. Each crossing mode has one source and one target group. The public geometric expansion at T=t+alpha z, i=v2(t), D_t=L/2^i, gives over the working dyadic ring

    W=A+sum_j 2^w_j u_j v_j^T,

where A is group-local, u_j is supported in the target and v_j in the source, and w_j are positive integers. There are at most c=O(M/D_t) modes of each weight across all packets. One obtains these vectors by extracting the public guaranteed 2^w_j factor from each crossing coefficient before reducing modulo the working precision. All modes of weight >= that precision can be omitted.

For a fixed minor, determinant multilinearity or the Schur complement with an identity block gives

    det W[Y,X] = sum_I (-1)^|I| 2^sum_(j in I)w_j
                   det [ A[Y,X]  U_I[Y,:] ; V_I[:,X]  0 ],                  (2)

where the sum is over subsets of modes and U_I,V_I list their normalized rank-one vectors. To prove (2) without any invertibility condition on A, expand det [ A U ; -V diag(2^w) I ] in the bottom diagonal identity entries, or first use formal rank-one variables and coefficient comparison. The identity is valid over every commutative ring. Terms with |I|>k are zero because the bottom rows have only k nonzero candidate columns. In contrast to the affine-selector expansion, no mode appears twice in a determinant.

For each group g, assign its original rows Y_g and the selected mode rows whose source is g to its row list; assign its original columns X_g and mode columns whose target is g to its column list. The bordered matrix is group-block-diagonal after reordering. If a group has unequal row and column counts, the determinant is zero (every permutation has a zero factor). Otherwise it is the product of the local square determinants times a permutation sign.

There are two kinds of sign. First reorder the original row and column subsets by group, leaving the mode rows and columns after them. This sign depends on X,Y but not on the mode subset. Next reorder mode ports into the groups and interleave them with the original group blocks. The second sign, including (-1)^|I|, depends only on the ordered mode labels and the original per-group cardinalities. Squaring (2), with independently selected subsets I,J, cancels the first sign exactly. Fix the per-group cardinalities at this interpolation evaluation; the remaining sign is public. The local factors are products of two bordered determinants, with zero returned for an incorrect cardinality or a rectangular group block.

The result is an exact sum of products of ordinary child contractions. Each child contains the original group vertices and one packet-group vertex joined to its retained original ports. The new packet entry requires only two determinants of size at most 2L and arithmetic in the working ring. Division-free determinant algorithms give polynomial-time and polynomial-space access over that ring. Dummy ports are constants. All selected mode labels and counts are public in the summand. This does not require splitting original tensors or computing arbitrary multivariate coefficients. The previous shared-mode Newton construction handles arbitrary target fanout; the exterior argument here uses the flat one-source/one-target mode geometry, not arbitrary fanout.

## 4. Sparse cardinality yields a different asymptotic promise

Let M be the Boolean edge count, lambda=log2 M. Assume the public contraction interface with polynomial descriptions, original maximum degree <=M^(4/5), entry time poly(M) exp(O(M^(4/5))) and polynomial space, P0=O(M). Suppose the product of row entries is zero at input precision on every common wire assignment with more than H ones. This promise need not imply that each individual row entry vanishes; it suffices that the original contraction terms with too many ones are zero modulo input precision. We omit their cardinality profiles altogether. For any fixed beta>0, suppose

    H <= M / (lambda^beta log2 lambda).

**Theorem (support-sensitive first level).** Under the public contractor, the contraction can be computed deterministically in polynomial space and time

    poly(M, description length) * 2^O(M/lambda^((1+beta)/2)).

Either bipartition side can supply the support promise. A known public baseline bit string can be removed by flipping those ports on both sides. This result does not require small multilinear degree on either side.

### Packing, degrees and precision

Use r=max(64,ceil(lambda^(beta/2))) groups, a=floor(lambda/(40 log2 r)), u=floor(a/2), L=2^a. These parameters are those of the growing-group argument in `contractor_degree_note.md`. The public second-moment packing construction applies: D_arr=r^(2a)<=M^(1/20), D_arr^2 Delta=o(M), and L D_arr=o(M), so padding is 1+o(1) and each child has O(M/r) wires. L grows faster than any fixed power of lambda, yet L=M^o(1). There are O(M/L) packets. Enumerating their k profiles, and at each node the at most two row and column group counts per packet, costs log time O((M/L) log L), negligible at any fixed logarithmic saving power. Process these branches sequentially.

Let H_*=max(H,ceil(M/lambda^(beta+2))) and use h=ceil(5H_*/u), multiplicities q_t=h(u-v2(t)) at levels v2(t)<u. The total multiplicity exceeds 4H_*L, the degree bound in Section 1. Section 2 ensures local integrality. The usual tapered Hermite calculation gives loss <=6hD_t+q_t(a+1) for the node t: the first term is the public valuation-difference estimate for its cardinal factor, and the reciprocal local series is integral; evaluating the truncated jet at -t/alpha contributes the second loss. Thus sufficient working precision is

    P_t=P0+2K+8hD_t+O(H_* a)=O(H_* D_t/u).

Indeed D_t>=2^(a-u) dominates all fixed powers of lambda. Each local polynomial can be recovered using polynomially many unit-separated nonzero-residue interpolation values. This explicitly pays the scale 2^(2K) and the interpolation readout. All root precisions are M^(1+o(1)).

### Crossing catalogue

For the two selected mode sets, the total selected weight is below P_t and their combined cardinality is at most 2K<=2H_*. For theta>0 and 0<z<=1 their number J satisfies

    J <= exp(theta P_t) z^(-2H_*) product_(w>=1)(1+z exp(-theta w))^(2c).

At z=1, the decreasing-function integral bound sum_(w>=1) ln(1+exp(-theta w))<=pi^2/(12 theta) gives log J=O(sqrt(cP_t)). This leading partition bound for two colors of distinct modes matches the ordinary one-color unrestricted bound; merely replacing powers by subsets is not by itself an asymptotic improvement. For the cardinality-sensitive alternative use ln(1+x)<=x, theta=H_*/P_t, and z=min(1,H_*^2/(cP_t)); together with the previous bound if H_*^2>=cP_t, this gives log J=O(H_* log(2+cP_t/H_*^2)). The case c=0 is trivial. Since c=O(M/D_t), we obtain

    log J=O(min(sqrt(M H_*/u), H_* log(2+M/(H_*u)))).              (3)

All feasible lists can be generated by recursive weighted-budget subset enumeration with polynomial delay and space. Per-packet cardinality caps may also be enforced but are not needed for the upper estimate. The first alternative in (3) suffices for the theorem.

### Charge the child calls

Each child has s=O(M/r) edges, maximum degree at most max(M^(4/5),2L)<=s^(9/10) eventually, and the new packet entries are polynomial-time. Its initial precision is M^(1+o(1)), rather than O(s). This is the same charged high-precision-child variant of the public proof as in `contractor_degree_note.md`: choose the child's own first u'=Theta(log s) with a sufficiently small constant. Its smallest retained D' is s^epsilon for a fixed epsilon>0; the inherited precision is o(sD'/polylog(s)). Hence it is absorbed by the working-precision bound at the child root. Subsequent depths use the cumulative public precision schedule, and directly enumerated leaves need only polynomial-size arithmetic. Reuse an extension with polynomially many residue elements (increase its fixed degree constant if necessary). The ordinary child time is poly(M) exp(O(M/(r sqrt(lambda)))) including the original entry access, which is exp(O(M^(4/5))). Depth-first space remains polynomial. The root enumeration therefore costs in log time

    O(min(sqrt(M H_* log r/lambda), H_* log(2+M/(H_*u)))
      + M/(r sqrt(lambda)) + M^(4/5) + M log L/L).

The claimed bound follows from the first alternative and the chosen r. Only the root needs the support promise.

**Genuine separation of promises.** For example take each row tensor to be the indicator of exactly one occupied port among d ports, with d odd of order lambda^(1/3) log lambda. Over F2 the coefficient of a squarefree monomial of degree j is j mod 2, by Boolean Mobius inversion. Its degree is therefore d, giving total row multilinear degree M. Yet H=M/d. For arbitrary column tensors satisfying the remaining interface this proves the 2/3 saving power while the earlier low-multilinear-degree theorem does not apply. This is a graph class, not an elliptic relation compiler. At typical full entropy the common-factor reduction is a real degree/precision improvement, but no numerical leading constant follows for the earlier separate simplex/finite-packet contractor.

## 5. A stronger seed-local obstruction, now for sparse support

This result uses the paired Fourier diagonal established in `seed_local_slice_rank_note.md`, whose proof applies to the actual random field-affine masks: translations and nonzero dilations of a nonconstant indicator span the zero-sum space on the field; its projection onto the proper rational-abscissa set is surjective. Restrict to symmetric functions on the odd prime subgroup of order r_E, vanishing at O. Over a splitting field of characteristic two the group Fourier transform pairs inverse characters. For a nonidentity target in the prime subgroup, chi(T)+chi(-T)!=0 for every nontrivial character. Ordinary count parity thus restricts to a nonzero diagonal tensor of length s_E=(r_E-1)/2, of slice rank s_E. No random-support approximation is used.

Consider a fixed seed-local Boolean network (each tensor reads at most one slot seed). Permit any public baseline bit string; every nonzero contributing assignment differs from it in at most H bits. For seed i let E_i be its incident edge union, c_i=|E_i|. The actual contraction, not a polynomial lift, is a sum over deviation sets J with |J|<=H, of a product of functions each reading one seed. The factor for seed i depends on the assignment only through J intersect E_i. Also sum c_i<=2M and sum |J intersect E_i|<=2H. These hold even when an edge joins two distinct seed-local tensors. Grouping summands by the first i with |J intersect E_i|<=k_i, if sum(k_i+1)>2H, gives the slice-rank upper bound

    sr(Z) <= sum_i sum_(j<=k_i) binom(c_i,j).

For any fixed i, sr(Z)<=2^c_i without a weight constraint. Both bounds follow directly by pulling the common one-seed factor out of each group. Different public baselines do not change the argument. It applies over the dyadic ring at a chosen precision as well. A fixed mask or public restriction on J is allowed; it only removes terms.

Suppose K such networks combined by a fixed seed-independent linear readout compute exact ordinary relation parity. Write B=log2(s_E/(2Km))>0. By subadditivity of slice rank, some network must have rank greater than s_E/(2K). It has c_i>B for every i. Choose k_i maximal with sum_(j<=k_i) binom(c_i,j)<=2^B, and set t_i=k_i+1. Each t_i<=c_i, and sum t_i<=2H; otherwise the preceding slice bound contradicts the rank. Hence

    B < log2(sum_(j<=t_i) binom(c_i,j))
      <= min(2 sqrt(c_i t_i), t_i log2(e c_i/t_i)).

The elementary binomial bounds follow by weighting the binomial generating function at t_i/c_i (and maximizing log(e x)/sqrt(x) for the square-root version). Summing and applying Cauchy-Schwarz gives

    MH > m^2 B^2/16.                                           (4)

A stronger entropy form follows by the log-sum inequality. If T=sum t_i and C=sum c_i, then sum_i t_i log2(e c_i/t_i)<=T log2(e C/T). Since T<=2H, C<=2M, and T log(e C/T) is increasing for T<=C, we obtain

    2H log2(e M/H) > m B.                                      (5)

Both (4) and (5) hold for the same network. If a packetwise support promise uses distance to the better of the all-zero and all-one packet assignments, enumerate the public baseline choices (one bit per packet); the argument still applies to these restricted assignment sums with K enlarged by at most 2^#packets. For the polylogarithmic r and arities in this checkpoint, #packets=O(M/L)=o(b) when M=O(mb). More generally the conclusion requires log K=o(b). The dyadic scaled-output variant is also valid exactly as in the earlier note: lift the local diagonal restrictions to an unramified extension, with output 2^(P-1)H_T modulo 2^P. The resulting diagonal has nonzero entries, and the diagonal slice-rank induction works over the finite chain ring using free kernel generators from Smith form. Arbitrary nonlinear readouts or seed-dependent coefficients are not covered.

### Both cost estimates face the obstruction

When the subgroup cofactor is polynomial, K=2^o(b), m=polylog(b), M polynomial in b, we have mB=Theta(mb). Equation (4) alone blocks improving the square-root parameter sqrt(MH/u) at u=O(log b). To check the cardinality-refined alternative in (3), put f=M/H>=1. Equation (5) implies H>=c mb/log2(e f). For f<=u^2 (u>=2), the ratio log2(2+f/u)/log2(e f) is Omega(1/log(2+u)), and for f>u^2 it is Omega(1). Thus

    min(sqrt(MH/u), H log2(2+M/(Hu))) = Omega(mb/sqrt(log b))

as a lower bound on this *upper-bound parameter expression* for the obstructed model. Enlarging H to H_* does not lower either term. If such a compiler is merely inflated with additional wires, neither the support-sensitive partition bound nor its cardinality refinement certifies an improved power. This is NOT a time lower bound for ECDLP, for the exterior deformation on particular inputs, or for a contractor that exploits cancellations among the selected mode sets. In particular nonunique verified certificates whose parity is not ordinary relation parity, joint seed processing, and approximate witness recovery escape the proof's hypotheses. The global section-norm candidate may have nonuniform certificate multiplicities; (4) is not an automatic obstruction to every use of it.

## 6. Reproduction and next research gate

The scripts use only Python's standard library and deterministic pseudorandom seeds:

    python -B verify.py
    python -B verify_readout.py
    python -B verify_sparse.py

The first compares exact rational Cauchy minors and the common factor; verifies the tight root orders, autocorrelation valuation identities, real geometric modes modulo powers of two, grouped bordered determinants, and complete finite squared-minor group joins. Its catalogue examples compare counts at P and P/2 by assumption, illustrating a possible degree-driven saving, not actual runtime or the full precision schedule. The second independently builds exact normalized polynomial sums, checks their entire local Taylor coefficient integrality, and reconstructs their origin by the tapered Hermite formula using only modular jet data. It does not claim to run the entire public asymptotic contractor. The third compares actual seed-local sparse contractions with the grouped slice decomposition and checks the strengthened entropy inequalities in finite examples. Retained transcripts give exact counts.

The new deformation achieves a strict asymptotic saving on a specified high-degree subclass, but a complete improved ECDLP algorithm still needs a new step. A concrete gate is to establish that the *actual norm-compiler mode sums*, after the bordered-determinant flow constraints, admit an implicit computation with cost o(mb/sqrt(log b)) at larger arity, or to build a verified-certificate sampler with small effective occupancy and sufficiently many genuine witnesses without requiring exact ordinary parity. A generic sparse recoding of the old local tensors is not sufficient. The best earlier candidate coefficient 11.31486814 remains conditional and has not been independently audited or changed by this work.
