# Mixed-characteristic residue extraction and collision-free scalar norm recovery

Research checkpoint. Continues `collective_trace_note.md`, especially its tagged rational trace and resolvent-only obstruction. **The general-prime asymptotic ECDLP exponent is not improved here.** We remove the separate residue-decoding stage by changing the scalar coefficient ring, and establish a one-call recovery theorem that needs neither random projection nor coordinate tags. We also give a base-ring alternative and a reverse-mode transfer for arithmetic circuits. These results do not provide the missing fast collective norm circuit. The prior numerical coefficient 11.31486814 is not re-audited.

## 1. Setup and the distinction from finite-field readout

Let p be an odd prime and S=∏ S_i be a finite Cartesian grid with S_i subsets of F_p. Set W=|S| and let F∈F_p[U_1,...,U_m]. Write H=|{s∈S:F(s)=0}| and let K≥1 be a promised output cap. Use any integer polynomial lift f of F and any integer representatives of the grid elements. For every grid point, f(s) reduces to F(s) modulo p. The arguments below need only this reduction; they hold even for arbitrary separate lifts of the values and tags. For actual norm evaluation, one uses a succinct polynomial lift, not a table of the W values.

For a polynomial G over a commutative coefficient ring R, write N_R(G)=∏_{s∈S}G(s), using lifted grid coordinates. Whenever R is a Z/p^e-algebra with these grid elements, their pairwise differences within a slot are units. Thus the split grid algebra over R is a product of W copies of R, and this product is exactly the determinant of multiplication by G. No determinant expansion is used in the proofs.

The preceding finite-field readout theorem addressed the information needed to reconstruct a residue from admissible values and jets of a base-field rational Cauchy transform, for arbitrary field-valued histograms. We instead change to mixed characteristic. These are different measurements and do not contradict that theorem.

## 2. One valuation-sensitive rational evaluation suffices for residue readout

**Lemma 1.** Lift a tagged rational trace to Q_p, where Q_p denotes the p-adic field:

    R_P(t)=∑_{s∈S} P(s)/(t+f(s)).

Let π²=p in the ramified quadratic extension and use its valuation normalized by v_π(π)=1. Then

    π R_P(π) is integral,
    π R_P(π) mod π = ∑_{F(s)=0} P(s) mod p.

Proof. If f(s) is a p-adic unit, π/(π+f(s)) has positive valuation and reduces to zero. If p divides f(s), f(s)/π has positive valuation; the ratio is 1/(1+f(s)/π), which reduces to one. The denominator cannot vanish: the valuations of π and a nonzero rational integer f(s) differ, and f(s)=0 causes no problem. Multiply by integral tags and add. ∎

A base-ring version does not even need ramification:

    T_P = ∑_s p P(s)/(f(s)^2+p) ∈ Z_p,
    T_P mod p = ∑_{F(s)=0} P(s) mod p.                  (1)

If f(s)≡0 mod p, the denominator has valuation exactly one; otherwise it is a unit. This proves (1) termwise. One only needs the integral rational trace modulo p (or the unscaled trace to absolute p-adic precision 0). This exact observation avoids reconstruction of a high-degree reduced rational function, but is not a fast collective evaluation of either sum. In particular, computing a trace over the product of W p-adic fields still need not be cheap.

## 3. Recover tagged residues from ordinary scalar norm ratios

**Theorem 2.** Given cap K, compute ordinary scalar norms in Z/p^(K+2)Z:

    N_0 = N(f^2+p),
    N_P = N(f^2+p+p^2 P).

Then detect H>K, or recover H and the tagged witness moment M_P=∑_{F(s)=0}P(s) mod p by

    H = v_p(N_0),
    M_P = ( (N_P/p^H) * (N_0/p^H)^(-1) mod p² - 1 ) / p.     (2)

When the valuation is above K, return overflow instead. Formula (2) means exact division of the represented residues by p^H and reduction of the normalized units modulo p², followed by the indicated integer division. For many tags, the baseline N_0 is reused.

Proof. Every negative factor f(s)^2+p is a unit, and every positive factor has valuation exactly one. Adding p² P does not change this. Thus both valuations equal H, and H≤K leaves at least two p-adic digits after normalization. On a positive address, the ratio of the perturbed to unperturbed factor is 1+pP(s) mod p²; on a negative address it is one mod p². The product of the ratios is 1+pM_P mod p², since every cross term has a factor p². This proves (2). No division by a nonunit is performed inside the grid algebra. ∎

This is a genuine finite scalar readout, with no spectrum reconstruction, no rational function interpolation and no extension field. The only collective primitives are scalar norms of two explicitly perturbed polynomials. The denominator degree doubles rather than becoming p-sized. The bit precision is O((K+2)log p) per coefficient. An evaluator needs to support these perturbed polynomials and the prime-power ring; an evaluator specified solely for F over F_p need not do so.

## 4. Circuit-aware reverse-mode recovery in one pass

**Theorem 3.** Suppose the scalar norm of the polynomial family

    f^2+p(1+∑_j λ_j P_j)

is represented by an arithmetic circuit of S_c gates over R=Z/p^(K+1)Z, as a polynomial identity in the perturbation parameters. Equivalently, a division-free norm algorithm may be used whose identities also hold after adjoining square-zero parameters. Allow unit divisions if their validity persists under the perturbation. Evaluate the circuit at λ=0 and reverse-differentiate its scalar output. In O(S_c) ring operations, plus the cost of constructing its inputs and reading its derivatives, recover H (or overflow) and all tag moments, with no separate norm call per tag:

    D_j = ∂N/∂λ_j at 0,
    M_j = (D_j/p^H) * (N_0/p^H)^(-1) mod p.            (3)

Proof. In Q_p, formal differentiation of the product gives

    D_j=N_0 ∑_s pP_j(s)/(f(s)^2+p).

By (1), the sum is integral, so D_j is divisible by p^H and its normalized ratio reduces to M_j. Since H≤K, precision p^(K+1) suffices. For reverse differentiation, initialize the output adjoint to one, then propagate adjoints backward: addition sends the same adjoint to both inputs, multiplication by input values swaps the multiplier, and inversion of a unit contributes -u^(-2). Each gate takes a constant number of ring operations. This proof uses polynomial identities, not differentiation of an arbitrary function on a finite set. ∎

On a singleton, p² precision in the BASE ring suffices. Take tags U_1,...,U_m; their moments are the witness coordinates. The reverse pass computes all of them with a constant-factor arithmetic overhead relative to the perturbed norm circuit. Forward evaluation must of course have valid symbolic inputs for these directions. For a coefficient-input circuit, supply the perturbations to the coefficient variables; derivative readout then includes the costs of those coefficient directions. General tags might be expensive to represent; the affine-linear tags needed on a singleton are inexpensive.

The reverse implementation stores O(S_c) circuit data/values. We do not assert polynomial space from a polynomial-space streaming forward algorithm. Recomputing/checkpointing may trade time for memory. Nor does a black-box scalar norm oracle automatically expose its arithmetic circuit. Those are material interface restrictions.

An equal-characteristic alternative is to use the family ∏(F+Z(1+∑λ_j P_j)) modulo Z^(K+1). Its z^H output and derivative coefficients satisfy the same normalized identities. A division-free polynomial-identity circuit over F_p lifts to that coefficient algebra, but merely agreeing with a function on the base field does not prove such a lift. The mixed-characteristic proof requires no large formal polynomial degree in its evaluation parameter.

## 5. One scalar norm call recovers all ≤K addresses, with no projection collision

Previous checkpoints used random linear projections and dual-number coordinate tags. A scalar mixed-characteristic ring can encode addresses injectively instead.

Let k=F_(p^m)=F_p[V]/g(V), with g monic irreducible of degree m, and let v denote its generator. The map

    e(s)=∑_{j=1}^m s_j v^(j-1)

is a bijection F_p^m→k. In particular, every distinct grid address has a distinct encoded value. Choose a monic irreducible h(Y) over k of degree K+1, and lift its coefficients to the ring with the lifted g. Define

    R_0=(Z/p^(K+1)Z)[V]/(g_lift),
    R=R_0[Y]/(h_lift).

Write θ for the image of Y. These are finite free modules over Z/p^(K+1)Z with residue fields k and k(θ). In R, a residue nonzero modulo p is a unit (invert modulo p and lift via a finite geometric series). The standard monomial bases make exact division by p^H followed by reduction modulo p well-defined whenever H≤K and the represented element is divisible by p^H. Since θ has degree K+1 over k, θ-e(s) is nonzero in the residue field for every s (K≥1).

**Theorem 4.** Compute the single ordinary scalar norm in R:

    Z=N_R(f^2+p(θ-e_lift(U))).                            (4)

If Z=0 in R, H≥K+1. Otherwise its p-adic valuation in R is H≤K, and

    (Z/p^H) mod p = C ∏_{F(s)=0}(θ-e(s)),
    C=∏_{F(s)≠0}F(s)^2 ∈ F_p^×.                         (5)

Recover the degree-H polynomial C∏(Y-e(s)) simply by reading the coefficients in the basis 1,θ,...,θ^K of the residue extension over k. Divide by its leading coefficient, factor over k, and decode the roots in the v basis. This returns precisely all distinct zero addresses, with no random-projection collision probability. It also handles empty instances and overflow exactly.

Proof. Each negative address has a unit factor in (4), reducing to F(s)^2. Each positive address gives a factor p times a unit, whose divided reduction equals θ-e(s), because f(s)^2/p≡0 mod p. Multiplying proves (5) and exact valuation H until saturation at K+1. Since the degree is ≤K, no reduction by the degree-(K+1) extension polynomial changes the normalized polynomial coefficients. The roots are distinct because e is injective. Ordinary finite-field factorization and basis decoding suffice. ∎

The tower can be constructed by sampling irreducible polynomials and testing them using standard finite-field polynomial arithmetic; this takes expected polynomial time in b,m,K. Factoring the known split squarefree polynomial also takes expected polynomial time, e.g. by randomized gcd splitting. The ring has rank m(K+1) over Z/p^(K+1)Z and each scalar uses O(b m (K+1)^2) bits. Ordinary ring arithmetic is polynomial in these parameters. No formal perturbation variable remains in the norm input. The residue field is large enough to encode addresses collision-free even when p is small relative to the number of potential roots. The ring is still nonreduced: the base prime is nilpotent. We do not claim to have removed all nilpotents or to work solely over a finite field.

**Singleton simplification.** When K=1, use just one extension of degree m+1 modulo p², with generator v, and define e_+(U)=v^m+∑U_j v^(j-1). It is never zero in the residue field. Compute Z=N(f²+p e_+(U)). On H=1, divide by p, reduce mod p, and normalize by the coefficient of v^m. The remaining m coefficients give the witness. On H=0 the norm is a unit, and on H≥2 it vanishes modulo p². This uses one scalar norm over a degree-(m+1) extension, with no extra θ extension and no factorization. If a symbolic circuit is available, Theorem 3 instead avoids the extension completely.

The following comparison is about readout only; the cost of evaluating a norm dominates unless a new evaluator is found.

| Method | Scalar calls | Coefficient ring | Recovery |
|---|---:|---|---|
| Prior projected-support trace | 1 + tags | F_p | High-power trace, projection for multiple fibers |
| Ratio (2) | 1 + tags | Z/p^(K+2)Z | Scalar divisions, Newton/Vandermonde if desired |
| Reverse derivative (3) | one forward and one reverse pass | Z/p^(K+1)Z | All supplied tags, symbolic circuit required |
| Packed scalar (4) | one | unramified tower modulo p^(K+1) | Exact factorization, injective addresses |
| Singleton simplification | one | degree-(m+1) extension modulo p² | Normalize m+1 coefficients |

## 6. Cost at the elliptic parameters and the remaining obstruction

Under the preceding common-translation sampler assumptions, a coefficient lift of the projected summation polynomial of separate degree d=h²2^(m-1) can be produced by prior interpolation at cost 2^O(m²+m log h+log b). Squaring it and adding the public perturbation does not change that scale. For one eliminated interval, the block-product algorithm applies over the new coefficient rings because it uses polynomial multiplication and monic remainder trees. The perturbation has degree at most max(2d,deg P); for coordinate tags and (4) this is 2d. Its circuit construction still costs soft-O(d+sqrt(dw)) coefficient-ring operations, after coefficient setup, with polynomial overhead in m,b,K. The proof of the earlier block-product cost is conditional on fast polynomial arithmetic over the coefficient algebra as stated there. A schoolbook implementation does not establish the fast asymptotic timings.

For K≤2^O(m+log b), polynomial h, balanced log w=b/m+O(log b), all the scalar readout and ring-construction costs above fit 2^O(b/m+m²) (with lower-order log b terms). If a collective norm evaluator with that charged complexity over the required ring were established, at m≈b^(1/3) the prior verified-relation reduction could target the 2^O(b^(2/3)) scale. This is a conditional transfer, not a theorem that such an evaluator exists or a fresh audit of every outer-reduction hypothesis. Recovered algebraic roots must be verified by rational lifts and projected elliptic addition under the preceding coverage assumptions.

Directly evaluating the product over all W addresses still takes W times a polynomial evaluation cost. Eliminating one slot by block products leaves W/w lanes, each with its own norm, so naïvely doing the remaining evaluations costs roughly (W/w) sqrt(dw) coefficient operations, in addition to setup, about 2^(b-b/(2m)+O(m²+log b)) under the earlier balanced conditions. This does not beat the old half-list method, let alone the conditional contractor. Neither the low p-adic precision nor the injective packing compresses those lanes by itself. Valuation of a product is additive, but deciding whether a particular factor is divisible by p still requires the relation test; computing a sum of those decisions is precisely the missing collective operation.

There is no conflict with the prior finite-field resolvent information lower bound. Equations (1)–(5) use mixed-characteristic evaluation, valuation and exact division in a prime-power ring, rather than admissible F_p-extension Cauchy queries. They trade the old global residue readout for a scalar determinant/norm evaluation at finite precision. We have eliminated that readout as a *separate* computational obstacle, not solved the determinant problem or provided a lower bound against arbitrary circuit-aware algorithms.

## 7. Reproduction

Run `python -B verify.py` in this directory. Standard library only. `previous_fiber_verifier.py` is copied unchanged from the previous checkpoint because the new script imports its actual Semaev-polynomial evaluator. The fresh tests do not count the previous verifier's results. The script constructs genuine irreducible polynomials in both residue-field levels, lifts their coefficients to prime-power polynomial rings, computes literal scalar products, performs finite-field splitting, verifies the base-ring ratio identities and reverse-mode arithmetic derivatives, and tests invariance under independent changes of the integer lifts. All grid enumeration is explicit: no faster collective norm is implemented. The reverse-mode test uses a balanced scalar product circuit, not an end-to-end fast determinant circuit.

The first run reached a deliberately checked curve nonsingularity assertion: (p,a,c)=(31,1,1) was singular because 4+27=31. It was replaced by (37,1,1). The final full rerun used the same seed and all assertions passed. See `verification.txt`. No failed mathematical identity was hidden, and this fixture correction does not alter any theorem. The actual Semaev tests cover unprojected algebraic zero addresses, not rational group witness verification for every cofactor class.

Final counts: 200 tower norm readouts, 200 degree-(m+1) singleton readouts, 200 base-ring tagged norm batches, 200 reverse-mode batches, 160 lift-invariance checks, 1199 scalar tag comparisons and 1199 reverse-mode tag comparisons. The grids contain 1292 addresses in total, including 72 actual Semaev grids. There were 40 exact overflow cases, and 104 recovered addresses across the within-cap fixtures, including 29 multiple-address fixtures. Reverse circuits used 33664 total tape nodes. Empty grids of witnesses, whole-grid zeros, and exact cap overflow are included.

## 8. Research frontier

We can read a residue from a single mixed-characteristic rational evaluation, or from ordinary scalar norm products. On a promised singleton, one scalar norm over a small degree extension suffices; a polynomial-identity norm circuit admits an even cheaper base-ring derivative readout. The full tagged problem admits a collision-free, bounded-precision one-call encoding. These strengthen the interface to a potential evaluator and remove the need to compute a rational-function spectrum. They do not beat the general-prime exponent. The open target is computing the scalar norm of the perturbed summation polynomial without expanding the remaining grid. The previous coefficient remains an unverified conditional candidate, not a new claimed bound.
