# Symbolic group splitting, fixed-budget moment ranks, and unit-separated XOR packets

24 November 2026. Continuation of `mode_contraction_note.md` and `section_norm_research_note.md`. The objective is a strict improvement on the scale `b/(log_2 b)^(1/4)` in the logarithm of the general-prime ECDLP time bound.

**Outcome.** We give an exact symbolic composition of two levels of the public packet recursion, with the parent mode indices retained as external arguments of the child tensors. We prove a useful pair-factor sufficient condition for a genuinely improved exponent, including charged factor access. We also exhibit actual packet-group tensors from copy and pin predicates which have exponentially large factorization rank over a finite dyadic ring, even on fixed-weight mode-index blocks. This rules out a uniform small-rank factorization based only on the moment form of packet tensors. Finally, for a different XOR deformation using unit-separated roots in an unramified extension, we prove that every off-coset crossing block has all its elementary divisors of valuation one. A generic interpolation witness shows why this particular deformation does not offer the requisite mode-weight spectrum. The actual unrestricted elliptic group tensors have not been shown to meet the new sufficient condition or to be subject to all these example-based obstructions. The general-prime exponent is not improved in this checkpoint. The earlier 11.31486814 coefficient remains a conditional candidate, not independently verified here.

## Sources and notation

The public source is *The Exponential-Time Hypothesis is False*, 25 October 2026, particularly the Boolean tensor contraction proof at https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-contraction.tex and the interface theorem at https://github.com/openai/math/blob/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/quantitative-reduction.tex . We use its crossing-mode expansion, the child-network identity and tapered Hermite reconstruction as given. The preceding notes supply exact pair-budget encoding and high-precision lifting. The results below on mode moments, ring factorization rank and XOR blocks are algebraic statements independent of the ECDLP reduction. The asymptotic runtime implications are conditional on the stated upstream contractor and reduction, and their indicated access and packing assumptions.

Let M be the original edge count, r the number of groups, u the retained valuation levels, Lp the packet width, i the valuation of the interpolation center, D=Lp/2^i, and P the working precision. Write L=log_2 M in runtime estimates. At the retained nodes the previous estimates give P=O(M D/u). Each raw geometric crossing mode has one source group and one target group, and a source form of the following kind on a tile of columns:

    L_j(X) = sum_t t^j X_t,        weight j+1, j>=0,

with 0^0=1. The tile has at most D/2 column positions. Empty source or target modes can be removed. The row-side coefficients of the forms are divisible by 2^(j+1). Each mode exponent list is truncated by its total weight. For an exact pair budget h_e and maximum number c_e of modes of one weight, the preceding note gives

    log_2 N_e(h_e) <= K sqrt(c_e h_e),
    K = pi sqrt(2/3)/ln 2,
    c_e = O(M/(r^2 D)).

Thus the largest individual pair index has length O(M/(r sqrt u)+1), whereas the sum over all pair indices can be O(M/sqrt u+r^2). The previous treewidth criterion used a weaker incident-degree bound, O(M/sqrt(r u)+r). The difference between pair degree and incident degree will matter below.

## 1. Exact symbolic splitting through two packet levels

The first-level decomposition at an interpolation node is

    F = sum_{beta: wt(beta)<P} product_g Z_g(beta_incident(g)) mod 2^P.

For a fixed beta, each Z_g is the child contraction whose tensor entries include the packet coefficient extractions and source moments specified in the public proof. Treat beta as a formal discrete external input rather than substituting its value into every entry. Each individual mode exponent occurs at exactly one packet vertex in its source group and one packet vertex in its target group. This locality is preserved if its entry is canonically lifted: the lift may depend nonlinearly on its own arguments but does not acquire new external arguments.

Apply one additional packet recursion to each Z_g. Fix a graph-only packing and interpolation schedule, independent of the particular beta. Each retained old vertex belongs to precisely one new subgroup. Thus each occurrence of an original external mode exponent now belongs to exactly one descendant subgroup. New modes have endpoints in descendants of a single original group. This gives an exact incidence description: the first-level modes connect descendants of different original groups, and the second-level modes connect descendants inside the same original group. No materialization of the full Z_g tables is necessary.

### Combining different precisions and Hermite denominators

An issue is that a reconstruction identity for Z_g involves rational coefficients and different local precisions. Let a_g be a common dyadic denominator exponent for the linear reconstruction of Z_g. Odd denominators are units and can be represented at the required precision. For each fixed beta the identity has the form

    A_g(beta) = sum_{xi_g} w_(g,xi_g) C_(g,xi_g)(beta)
              = 2^a_g Z_g(beta) mod 2^(P+a_g).

Here xi_g includes interpolation choices and second-level mode lists, C is the product of the appropriate descendant contractions, and w are scaled reconstruction coefficients. Each individual C may be computed at its own prescribed precision: its reconstruction weight absorbs any difference between that precision and P+a_g. The public Hermite precision guarantee is exactly what is required for this congruence. Canonical lifts of the descendant entries to a common ring are valid, since their changes affect the reconstructed A_g only by multiples of 2^(P+a_g). The description of a descendant entry includes only the external exponents of the retained vertices that it contains, and the ordinary local recursion parameters.

Put A=sum_g a_g. For any choices of valid representatives, the product satisfies

    product_g A_g(beta) = 2^A product_g Z_g(beta) mod 2^(P+A).

Indeed write A_g=2^a_g (Z_g+2^P e_g) for integral e_g and multiply. Hence the entire first-level answer is recovered from the sum over beta of the products of the scaled reconstructions modulo 2^(P+A), by dividing the final residue by 2^A. Expanding the product of sums is now legitimate. There is no division inside an outer mode contraction. The same reasoning works coordinatewise in the unramified coefficient ring. It is essential to retain the lower-precision truncation/creation-time rules in each individual factor; simply interpreting every entry as a fresh high-precision tensor would not be justified.

For fixed interpolation choices the remaining graph is a mode graph on the descendant subgroups. The scaled coefficients are scalar factors. Parent exponent lists are repartitioned by descendant endpoint pairs, while new exponent lists remain within original groups. Enumerate exact budgets for these refined blocks, imposing the original first-level total weight bound and each second-level weight bound. Unrank each block as in the prior note. The number of budget and interpolation choices has logarithm poly(r,r',log M) for r,r' polylogarithmic, and all parameters and dyadic losses are polynomial in M. Here r' bounds the second-level number of groups. The final high precision is removable by the prior digit-lifting theorem when the number of nontrivial tensors is polylogarithmic. This is an exact way of flattening two recursive levels into one join without treating group tensors as indivisible black boxes. It does not imply that the flattened join is cheap.

### Wire accounting

The first-level mode entropy bound is unchanged by repartitioning a mode list into more blocks: for any fixed exact budget allocation, their product of cardinalities is bounded by the same global feasible universe. The rounding overhead increases to O((rr')^2). Summing the second-level entropy bounds costs O(sum_g M_g/sqrt(u')), where M_g is the number of edges in the g-th child. The child edges count the original ports, so sum_g M_g=2M. Consequently the combined edge budget is

    S_flat = O(M/sqrt(u) + M/sqrt(u') + (rr')^2).

This formula requires the stated precision bound on each second-level mode expansion, e.g. a cumulative packet scale large enough to absorb the inherited P. No claim that it holds for arbitrary choices of child widths is intended. Choose compatible polylogarithmic group counts and allocate the available packet scale accordingly. Even granting the most favorable analogous degree balancing and the degree-ratio contractor on the flat mode graph, using u=Theta(L/log r) and u'=Theta(L/log r') in the upper-bound expression yields

    O((M/sqrt L) (sqrt(log r)+sqrt(log r'))/sqrt(log r+log r')).

The ratio lies between 1 and sqrt(2). Thus these standard estimates for a symbolic split do not by themselves improve the order; the extra logarithmic grouping gain is paid for by new modes. This is a comparison of upper-bound expressions, not a lower bound on actual tensors or a proof that a symbolic algorithm cannot improve the exponent. In particular, we have not asserted a maximum-degree guarantee for arbitrary descendant allocations of parent mode incidences. The factorization criterion below avoids requiring such a degree guarantee on the outer graph.

## 2. Pair-factor join: a more direct sufficient condition

For a fixed budget branch, let beta_e denote the index on pair e={g,h}, of alphabet size N_e. Suppose each group tensor has an exact factorization over the coefficient ring of the form

    Z_g(beta_incident(g)) = sum_{a_g in {0,1}^c_g}
        A_g(a_g) product_{h: {g,h} present} F_(g,h)(a_g,beta_{g,h}).

Unused auxiliary indices can be padded with zero. This is a sufficient condition, not a hypothesis established for the elliptic compiler. More general decompositions using internal binary edges of total size c_g can be put in this form by fixing those edges. Let C=sum_g c_g and E_f be the worst charged cost of evaluating the factors and A_g, including any sub-contraction or preprocessing needed for a specified entry. Then the outer contraction is exactly

    sum_a [ product_g A_g(a_g) ]
          product_{e={g,h}} [ sum_{beta_e} F_(g,h)(a_g,beta_e) F_(h,g)(a_h,beta_e) ].

The proof is distribution of finite sums and locality of the two occurrences of beta_e. The logarithm of the time, ignoring polynomial factors but including the exact-budget enumeration, is bounded by

    C + max_e ceil(log_2 N_e) + log_2(1+E_f) + O(r^2 log M).

This bound does not require low treewidth of the catalogue. It uses only the largest individual pair alphabet, not the total incident alphabet. It also permits dense pair interaction graphs. The factors must have a uniform constructive description, and the cost of obtaining them must be charged; an existential tensor decomposition is not enough.

At r=Theta(L^(1/6)), u=Theta(L/log r), the individual pair alphabet bound is O(M sqrt(log L)/L^(2/3)). If factor access is no worse than a charged ordinary child entry, its logarithm is O(M/(r sqrt L))=O(M/L^(2/3)). Therefore the sufficient target on the auxiliary index budget is

    C = O(M sqrt(log L)/L^(2/3)).

With the preceding assumptions, the relation compiler M=O(m b) and inherited outer reduction transfer this to

    log_2 expected ECDLP time = O(b (log log b)^(1/4)/(log b)^(1/3)).

The balance is b/m versus m b sqrt(log log b)/(log b)^(2/3), with m of order (log b)^(1/3)/(log log b)^(1/4). This is a strict asymptotic improvement *if* the pair factors are constructed with the claimed total index and access budgets. Unlike the previous bounded-treewidth sufficient condition, this permits the original mode graph to remain dense and requires only r of order L^(1/6). But proving a sufficiently economical factorization for the actual elliptic group tensors is the missing theorem. We next give a counterexample to any assertion that the source moment form alone guarantees such a factorization.

## 3. Full source-moment rank at small weight

Take one tile with k live Boolean source bits at consecutive positions t=0,...,k-1; all other bits in this tile are pinned to zero. Denote its first k raw moments by L_j(X), j=0,...,k-1. Let

    V_(beta,X) = product_j L_j(X)^beta_j,

where beta ranges over nonnegative exponent lists with sum_j (j+1) beta_j <= k^2 and sum_j beta_j <= k. Let d=(k-1)!.

**Lemma.** There is an integer matrix U, using only these mode indices, such that

    U V = d^k I_(2^k).

**Proof.** The Lagrange basis polynomials at 0,...,k-1 are

    ell_t(z) = product_(s!=t)(z-s)/(t-s).

The denominators divide d. The integer coefficients c_tj of d ell_t(z) satisfy

    sum_j c_tj L_j(X) = d X_t.

For every fixed binary word x, multiply these linear forms for its one positions and the forms d minus these linear forms for its zero positions. The result is d^k times the indicator of X=x. Its ordinary degree is at most k, and each mode has weight at most k, hence all its monomials have weight at most k^2. Their coefficients form the asserted rows of U. QED.

This is stronger than an unrestricted rational-span statement: the exact denominators are explicit and small enough in the dyadic scale. It shows that the weighted mode catalogue can encode every function of k live bits, not just low-order information about their sum.

### Fixed exact weight, without changing the source moments

The outer representation enumerates exact pair budgets rather than a single upper bound. The above rank therefore needs to survive a fixed-budget block. Add a separate source tile in the same endpoint-pair bundle, whose sole live bit is pinned to one at position t=1. Its source forms all have value one. Use two sentinel modes of weights k and 1 (for k=1 one mode suffices). For any original beta of weight w<=k^2, append exponents

    q=floor((k^2-w)/k),  s=(k^2-w)-k q.

The new weight is exactly k^2 and the source value is unchanged. The total mode degree is at most 3k. Choose the tile size and number of target rows to allow this degree. This is a statement about the source-side group tensor on the indexed block; it does not assert that each of these mode exponents has nonzero target coefficient in every possible target tensor. Thus the same matrix U V=d^k I is realized at one exact budget, with the other mode exponents pinned to zero.

## 4. Genuine group tensors with exponentially large rank

Use two such source tiles, in different pair bundles incident to a single group g, and make the k source bits of the two tiles identical by k binary equality predicates in that group. Pin unwanted column bits to zero and the sentinel bits to their required values. At row ports in this group use independent all-ones unary tensors: summing the row selectors l_0+l_1=1 eliminates their local factors; the incoming mode exponents are set to zero. These are actual child-packet source factors and ordinary copy/pin tensors, not arbitrary entries assigned to a finished group table. The geometry is realizable by packets with the two appropriate endpoint labels at the retained coordinate, choosing even-position source tiles at an odd interpolation node. Other coordinate labels and padding may be supplied independently. This construction is a local packet configuration; it is not a claim that a full unrestricted section-norm graph has this configuration at every scale or that it follows every balanced labeling schedule.

The group tensor across the two outgoing pair bundles has submatrix

    G_(beta,gamma) = sum_(X in {0,1}^k) V_(beta,X) V_(gamma,X)
                   = (V V^T)_(beta,gamma).

It follows that

    U G U^T = d^(2k) I_(2^k).

Work modulo 2^P (also valid in an unramified dyadic extension). If P>2k v_2(d), the diagonal scalar is nonzero. The minimal factorization inner dimension of this scaled identity over the finite local principal ideal ring is 2^k: its image is the direct sum of 2^k copies of a nonzero principal ideal, requiring 2^k generators. A submodule of a module generated by q elements over this ring needs at most q generators, by the structure theorem over a principal ideal ring (equivalently lift to a free module over the DVR). Thus if G factored through q<2^k coordinates then so would U G U^T, a contradiction.

Consequently the exact fixed-budget group tensor requires at least 2^k summands in any separated representation across these two pair bundles. This applies to the pair-factor condition in Section 2, which would require c_g>=k in this fixture. The bound is on exact dyadic factorization rank; it is not merely a statement over the rationals where denominators might hide the signal.

For K independent copies the same proof uses tensor products: the rank lower bound is 2^(K k) whenever

    P > 2 K k v_2(d),

and the two fixed pair budgets together cost 2 K k^2. Since v_2(d)<=k-1, choosing k a sufficiently small constant times D/sqrt(u) fits the usual precision and total weight bound when K=O(M/D), P is of order M D/u, and D/sqrt u grows. In a group with K=Theta(M/(rD)) such live blocks, this gives log factorization rank Omega(M/(r sqrt u)). Thus saturation at the natural incident mode-entropy scale is possible already in a very simple packet-generated family. This does not give a lower bound on arbitrary group tensors from actual elliptic norm equations, or on a contractor using nonlinear operations, a different basis jointly with adjacent group factors, or a different precision profile. It shows why a proof of small group factors has to use further structure of the original tensors, not just the fact that their packet source factors are moments of Boolean bits.

## 5. Unit-separated XOR deformation: a full-rank crossing obstruction

The integer XOR diagnostic in `contractor_degree_note.md` measured Walsh eigenvalues but did not prove an asymptotic rank obstruction. Here is a different geometry with exact unit-separated interpolation, and an exact theorem for its first crossing coefficient.

Let q=2^a and let O be an unramified dyadic ring with residue field F_q. Label packet positions by F_q (addition is XOR in a binary basis). Choose residue-distinct representatives t_delta in O, including t_0=0, and put

    C(T)=product_delta (T-t_delta),
    W_(x,y)(T)= C(T)/[C'(t_(x+y))(T-t_(x+y))].

This is the cardinal interpolation polynomial; all root differences are units. At T=t_tau, W is exactly translation by tau and hence flips every set coordinate of tau. It looks attractive for grouping simultaneously by many bits without the integer-grid valuation obstruction.

Group positions into additive cosets of a subspace H of F_q, with the opposite row/column orientation on the selected coordinates as in the ordinary packet construction. A crossing block has row coset A+H and column coset B+H with eta=tau+A+B not in H. At a local perturbation T=t_tau+2z, the crossing matrix vanishes at z=0. Its first Taylor coefficient divided by 2 reduces modulo 2 to

    K_(h,k)=1/(eta+h+k),       h,k in H.

Indeed C mod 2=X^q-X, whose derivative is one in characteristic two. All ratios C'(t_tau)/C'(t_(x+y)) therefore reduce to one; all denominators of the crossing block are units. This is independent of the chosen lifts.

**Theorem.** K is invertible over F_q for every additive subspace H and every eta outside H.

**Proof.** K is the convolution matrix of the function h -> 1/(eta+h) on the elementary abelian two-group H. In characteristic two its group algebra is F_q[Y_1,...,Y_d]/(Y_1^2,...,Y_d^2), by taking Y_j=g_j-1 for a basis of H. An element is a unit precisely when its augmentation (sum of its coefficients) is nonzero. Let P_H(X)=product_(h in H)(X-h). This is an additive polynomial, proved by induction on a basis of H using

    P_(H+<v>)(X)=P_H(X)^2-P_H(v)P_H(X).

It has nonzero constant derivative, since its roots are distinct. Therefore the augmentation is P_H'(eta)/P_H(eta), nonzero. Convolution by a unit is invertible. QED.

For any unit z, the full perturbed crossing block divided by 2 reduces to z K and is invertible. Equivalently all its Smith factors have valuation exactly one. Thus an integral rank-one mode decomposition of that block needs at least |H| modes of weight one: reducing the block divided by 2 modulo 2 gives a full-rank matrix, and modes of weight >=2 disappear. The assertion concerns the usual linear crossing-mode representation with integral factors and each contribution divisible by its claimed power of two. Arbitrary nonlinear or nonintegral decompositions are not covered.

In a fully occupied packet with nontrivial grouping, each source coset has some crossing target coset. The ordinary pair-local catalogue therefore contains at least Lp=|F_q| weight-one modes (and may contain many more). This is the opposite of the low-rank, growing-weight dyadic spectrum that the integer Cauchy contractor exploits.

### A degree-only interpolation precision obstruction

There is also a simple exact readout obstruction for this unit-separated geometry. Generic deformed contractions have degree at most M(q-1). Consider

    f(T)=product_(tau!=0) (1 - T/t_tau)^M.

All t_tau with tau!=0 are units. This polynomial has degree M(q-1), f(0)=1, and integral coefficients. At each nonzero center, f(t_tau+2z) has every coefficient divisible by 2^M, because the factor at that center is (-2z/t_tau)^M and the other factors are integral. Thus a generic polynomial-jet readout (even adaptive and nonlinear) using precision at most M at all nonzero centers cannot distinguish f from zero. At some required center it must obtain precision >M. This is a statement for the generic integral polynomial class of the degree bound, not necessarily for the actual image of the elliptic deformation.

Under full crossing occupancy at that center, the weight-one mode catalogue has Omega(M) independent indices, so the unrestricted weighted exponent enumeration at the required precision contains 2^Omega(M) binary lists (select a small fixed fraction of these indices). They can also respect the total degree cap M. This explains why the simple unit-separated XOR variant does not improve the existing generic readout/enumeration algorithm, despite its ability to flip many label coordinates at once. It is not a lower bound on a different interpolation family or on evaluation methods that never enumerate the list.

## 6. Finite verification

Run `python -B verify_group_tensors.py`. It uses only the Python standard library and exact integer and finite-field arithmetic.

* For k=1,...,6, construct integer indicator expansions directly from Lagrange basis polynomials and evaluate them on all 2^k words. Check their degree and weight bounds. Check sentinel padding to an exact weight and the resulting scaled Gram identity modulo a power of two just above the relevant loss. This is an independent calculation of the source moment matrix and the transformation; it does not run a complete packet recursion or certify the target-side support of all selected mode indices.
* Over F_(2^a) for a=2,...,6, enumerate every coordinate subspace and every disjoint affine coset representative. Calculate the crossing inverse matrices directly, check full rank by finite-field elimination and independently verify the augmentation through the product derivative identity. The proof covers all additive subspaces, not only the coordinate subspaces in the finite tests.
* Test the scaled symbolic product reconstruction with signed integer summands, unequal losses, multiple parent indices, and deliberately noncanonical residue lifts. Compare with direct products modulo the original precision. This tests the precision identity, not an end-to-end nested packet contractor.

The recorded output has 5,460 exact moment indicator identities, 53,000 fixed-budget comparisons, 5,460 transformed Gram comparisons, 965 XOR crossing blocks with 36,455 field matrix entries, and 200 symbolic joins with 65,888 summands and 6,305 changed entry lifts. All assertions passed.

## 7. Status and next research gate

The proposed symbolic splitting is an exact alternative to treating the entire group table as opaque. The pair-factor theorem narrows the missing mathematical requirement: for r of order (log M)^(1/6), construct pair factors with total auxiliary bit count at most O(M sqrt(log log M)/(log M)^(2/3)), with the actual factor access, precision and packing charges. The source-moment examples prove that no such guarantee follows from packet moment algebra alone. The unit-separated XOR alternative also fails the generic linear-mode approach by an exact rank and readout argument. A productive next step would have to exploit additional algebraic properties of the full elliptic constraint tensors *jointly* with adjacent source and target mode factors, or replace the interpolation/expansion algorithm entirely. The general-prime logarithmic exponent remains open. No numerical constant or new end-to-end algorithm is claimed in this checkpoint.
