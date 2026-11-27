# Elliptic-curve relation counting

**Project lead and repository maintainer** · [David Iach (@davidiach)](https://github.com/davidiach)

**Current manuscript · 27 November 2026.** *A polynomial certificate for deterministic elliptic-curve relation systems* studies elliptic curves over prime fields of characteristic greater than three whose rational points form a group of prime order. Assuming Liu–Zhang–Ren (LZR) Theorem 2.2 on closed Boolean tensor contraction, the fixed-shift equations give a matrix polynomial of full generic rank. For field bit length $b$, the resulting Las Vegas algorithm verifies each answer and takes $2^{O(b/(\log_2 b)^{1/4})}$ expected bit time and $O(b)$ expected random bits. With $O(b)$ bits of advice per input length the same time bound holds deterministically. [Manuscript (PDF)](research/fixed_shift_gram/manuscript.pdf) · [Typst source](research/fixed_shift_gram/manuscript.typ) · [Theorem and assumptions](MAIN_RESULT.md)

**What is proved about deterministic evaluation?** Stacking evaluations over the full field grid gives full column rank on every eligible curve. At the prescribed generator, the unweighted shift matrices also retain full shift rank when kept separate by source tuple, or by source group value and one specified slot on each side, under the short list’s parameter conditions. A finite example outside that regime illustrates how rank can fall when classes are added. On the order-19 curve the class matrices have joint rank eight and their sum has rank seven. Whether the proposed short list of weighted evaluations always has full column rank remains open. [Deterministic points](DETERMINISTIC_POINTS.md) · [Neighboring-fiber proofs](NEIGHBORING_FIBERS.md) · [Recorded ranks](evidence/neighboring_fibers.txt)

The manuscript also reduces the number of terms in LZR’s deterministic partial expansion. The residual contraction still gives the uniform deterministic time bound $2^{O(b/\sqrt{\log_2\log_2 b})}$ under LZR’s Appendix D assumptions. [Expansion proof](PARTIAL_EXPANSION.md) · [Bit complexity](COMPLEXITY_AUDIT.md) · [Open questions](OPEN_PROBLEMS.md)

**Common translation covers.** Under the same closed-contraction assumption, a constant-size block of translates gives full coverage with tuple product below the square of the group order at the quarter-power setup cost. The quotient reduction still requires an efficient deterministic algorithm for two-dimensional additive range avoidance to obtain the corresponding end-to-end bound. The finite verifier checks the coefficient-bit choices and a single pair of translations that repairs two different deletion subsets. [Coverage theorem](research/translation_cover/RESEARCH.md) · [Recorded checks](research/translation_cover/verification.txt) · [Quotient reduction](research/universal_quotient_compression/RESEARCH.md)

**What exact moments determine.** For abstract integer row populations, the moment bounds depend on the active weight. The homogeneous degree bound is sharp at weights three, four and five under the respective assumptions. The weight-three example has fewer rows than the modulus. The weight-four and weight-five examples use larger multiplicities. [Moment bounds](research/exact_moment_oracle_barriers/moment_rank_boundary/RESEARCH.md) · [Newton certificates and examples](research/exact_moment_oracle_barriers/moment_rank_boundary/newton_sparse_rank/RESEARCH.md) · [Research ledger](RESEARCH_LEDGER.md)

## Earlier paper and version history

*Counting relations in independently hashed elliptic-curve slots* gives numerical bounds under a separate finite-packet assumption.

[Read the v4 paper (PDF)](paper.pdf) · [Editable Typst source](paper.typ) · [Mathematical supplement](derivations/README.md) · [Recorded experiments](evidence/) · [Historical reports](archive/README.md)

**Version 4 · 24 November 2026.** For general field primes in the paper’s curve and subgroup class, independent field-affine thin-slot masks and a section norm give leading coefficients of **11.31486814 for log expected bit time** and **5.65743407 for log space**. Both multiply $Z_b=b/(\log_2 b)^{1/4}$. The result assumes the finite dyadic packet interface (A.1) and the inherited contraction and reduction estimates. Section 2 proves that the growing composite-degree graph meets A.1, bounds every proposed quotient before reducing modulo the Fourier modulus, and checks the raw field masks for collisions before rational-point filtering. Table 1 retains **24.00246** for the general-prime mixed fixed-bit sampler. [Section-norm derivation and finite checks](archive/research/chat_2026_11_23_24/08_section_norm/)

For a fixed elliptic-curve target $T$ with $2T\ne O$, any adaptive predictor that leaves a whole slot unread has error at least $\lambda-O(\lambda^2)$ at the proposed sparse widths. Joint determinant moments can be computed within the proposed budget. Decoding the relation count from them at growing arity remains open.

**Version 3 · 22 November 2026.** Table 1 compared conditional time and space bounds by field prime and sampler. Within the inherited input class, the reported leading time coefficient was **26.29338** with independent affine hashes and **24.00246** when the algorithm chooses the mixed sampler. The mixed sampler gives **19.59793** for primes with bounded signed radix support and **13.85783** for primes $p=2^b-c$ with $1\le c\le b^C$ for fixed $C$. These rounded coefficients multiply $b/(\log_2 b)^{1/4}$ in the upper bound on log expected time. That revision derived its four rows under finite-packet Assumption A.1 and the inherited reduction estimates. Sections 4–5 give the weak-detection criterion and its finite test. Section 7 covers the target-derivative experiments and the dense three-slot decoder. [Version 3 PDF](https://github.com/davidiach/elliptic-relation-counting/blob/b85f1f0f4969fe84c68be6245c2f8e18e7dfd6ee/paper.pdf) · [Group-equation derivations](derivations/group_equations/results.md) · [Fixed-bit sampler](derivations/strided_slots/results.md) · [Decoding results](derivations/joint_parity/dense_weighted.md)

**Version 2 · 21 November 2026.** The twelfth moment resolves all 168,480 eligible counts in Table 2 over $\mathbb F_{13}$. This revision also made the planting hypotheses explicit and partitioned the fresh hits in both larger-field trials. [Version 2 PDF](https://github.com/davidiach/elliptic-relation-counting/blob/9dcbb4ed9e58fd69d2adf0f39c0deeb25bfbf171/paper.pdf)

### Results

At slightly narrower mixed-sampler widths, almost every successful seed has one signed relation whose address takes $b+O(m)$ bits. A verifier that keeps each seed local faces a communication lower bound. A signed determinant score computed from all public seeds has exponentially small correlation with that address at growing arity. [Sampled-slot results](derivations/sampled_slots/README.md) · [Joint determinant test](derivations/sampled_slots/signed_moment.md)

**False candidates on isolated witnesses, 21 November 2026.** Translating the signed leaves before checking only abscissas inside the chord tree keeps every genuine witness. For each eligible fixed target the mean number of wrong signed addresses on guarded singleton seeds is at most $(1+o(1))/(2m^2)$ at the thinned widths. The direct tree still carries $m-2$ full-field abscissas. Several checked offsets give zero wrong addresses for the same proposed tuple. The conditional end-to-end exponent is unchanged. [Filter, evidence and comparison](derivations/false_candidates/README.md)

**Chord equality without storing either half list, 21 November 2026.** An exact algorithm trades space for time, including on empty targets. Keeping tables of $2^{b/3+o(b)}$ entries costs $2^{2b/3+o(b)}$ time. An isogeny can collapse a whole kernel coset to a single chord-function evaluation, but even adaptively chosen dense cosets cover at most $2^{b/4+o(b)}$ of the $2^{b/2-o(b)}$ distinct half sums on typical successful seeds. A much cheaper predictor meets the average masked-target error budget and loses the bit on isolated witnesses. The conditional general-prime exponent is unchanged. [Algorithms, limits and checks](derivations/implicit_chords/README.md)

**Product differences modulo $p^2$, 21 November 2026.** A base-field algorithm now decides equality of the two lifted determinant products on regular grids, including empty fibers, without constructing either scalar. It recovers any chosen bit on a singleton. Grouping the elliptic chord slopes keeps its cost at the scale of the two half lists. We also show that two calibrated equality comparisons suffice to decide target existence. The conditional end-to-end exponent remains unchanged. [Proofs, costs and finite comparisons](derivations/product_differences/README.md)

**Further parity calculations, 21 November 2026.** A pencil of $p+1$ proposed sections can now be checked in time proportional to one slot width, up to polynomial factors in bit length and arity. We also give an exact parity formula in which every nonwitness factor has quadratic character $+1$, and a nonlinear predictor with correlation at least $1/(4\sqrt K)$ per candidate when a coarse join leaves $K$ pairs. Choosing the pencil, evaluating that character, and constructing the coarse join each have separate costs. The conditional general-prime exponent remains unchanged. [Four approaches and reproducible results](derivations/direct_parity/README.md#further-attempts-at-the-same-exponent)

**Direct witness parity, 21 November 2026.** The signed zero count $M_0-M_{p-1}$ can be computed as one coefficient of a product formed from the slots. On a regular singleton a single elliptic section gives the whole witness. A proposed section can be checked in polynomial time from the seeds, even when the slots are exponentially wide. Finding it cheaply remains open. We derive the first degree that reveals it in the explicit quotient and bounds for scalar powering and sparse Fourier evaluation. The conditional end-to-end exponent gains no improvement from these methods. [Proofs, costs and finite checks](derivations/direct_parity/README.md)

**Nonlinear signed moments.** Over the field of order 13, an exact search found a two-moment table that favors the correct tuple on one target. In six-slot trials over the field of order 65,521, a cheap nonlinear character score ranked the witness above a wrong legal tuple at essentially the chance rate on two independent targets. We prove bounds for random coefficient searches and target-specific sparse lookups, and show how an exact parity decoder can require exponentially many signed determinant values. These moment experiments left the mixed fixed-bit bound at $(24.00246+o(1))b/(\log_2 b)^{1/4}$. [Results, costs and open cases](derivations/nonlinear_moments/README.md)

| Question | Result | Scope |
|---|---|---|
| Can the hash be encoded more cheaply? | Checking the short inverse image reduces the leading arithmetic word budget from 38 to 36 per slot per field bit. The time coefficient falls from 27.0138774237 to 26.2933775156, a 2.6671% reduction. | Conditional on the supplied finite-packet Assumption A.1 and inherited graph and relation-reduction estimates. The scale is $b/(\log_2 b)^{1/4}$ in the log of expected bit time. |
| Is evidence from almost every slot enough? | We derive the exact marginal likelihood under planting. Any adaptive predictor that omits a whole slot has fresh error at least $\lambda-O(\lambda^2)$ at the proposed balanced widths and fixed target with $2T\ne O$. | An access bound for this seed law. Joint moments depend on every slot. |
| Can joint moments be computed on a smaller grid? | Elliptic interpolation computes moments up to order $\lfloor 2^{b/m^2}\rfloor$ within expected log cost $(1+o(1))b/m$ at the proposed growing arity. | Computing a statistic and decoding its zero count have separate costs. |
| Can the moments reveal the exact count? | On the field of order 13, nine pooled moments certify 145,776 of 168,480 eligible counts. Planted calibration cuts nine-moment existence errors from 13,632 to 7,680. The raw twelfth moment determines all 168,480 counts. | Regular three-slot inputs with two sign orbits. On the field of order 1,031 the tested rules missed all 2,929 fresh positives having more than two orbits and regular masks. |

The paper also treats guarded-to-ordinary count conversion, fixed-target isolation, generic packet span, optimal separated approximation, derivative readout, and the time–space tradeoff for implicit resultants. The square-root-in-bitlength objective and a cheaper decoder for typical large pools remain open. These comparisons concern asymptotic formulas and finite counting experiments. This repository has no demonstrated end-to-end cryptographic-scale speedup.

## Reproduce the experiments

From the repository root, run `python -B research/translation_cover/check_translation_cover.py` with Python 3.10 or later and assertions enabled. It uses the standard library. [Expected translation-cover output](research/translation_cover/verification.txt)

Run the three fixed-shift checks from the repository root in a Python 3.12 environment. Each prints a transcript for comparison with the recorded results.

```sh
python -m pip install -r requirements.txt
python -B checks/check_fixed_shift_gram.py
python -B checks/check_fixed_shift_points.py
python -B checks/check_neighboring_fibers.py
```

The checks cover the symbolic certificate on finite examples, the evaluation lists, and the rank of the conditioned sums at the prescribed generator. [Gram transcript](evidence/fixed_shift_gram.txt) · [Point-list transcript](evidence/fixed_shift_points.txt) · [Neighboring-fiber transcript](evidence/neighboring_fibers.txt) · [Research ledger and source checks](RESEARCH_LEDGER.md#reproduction)

### Earlier paper and archive checks

Run `python reproduce.py` for the v4 suite. The recorded runs used Python 3.12.2. The default suite compares 52 outputs with their committed transcripts. To include the four C38 integration checks, pass the directory containing the original `ECDLP_v6_complete.zip` and `simplex-alphabet-ecdl-v6.zip`. The runner checks both archive hashes, loads four reference modules into a temporary directory, and compares all 56 outputs. The historical checkpoint replay has its own inputs and instructions. [Historical checks](archive/README.md#replay-the-historical-checks)

Run the section-norm component check separately with `python -B archive/research/chat_2026_11_23_24/08_section_norm/section_norm_audit.py`. It compares regular norm certificates with direct group sums and exercises lifted residuals on two finite polynomial rings. The growing-degree bound follows from Section 2. [Recorded output](archive/research/chat_2026_11_23_24/08_section_norm/section_norm_audit_output.txt)

```sh
python reproduce.py --source-archives /path/to/local/archives
```

| Script | What it checks |
|---|---|
| [Fixed-bit compiler](checks/strided_compiler.py) · [Compiler checks](checks/check_strided_compiler.py) · [Reciprocal](checks/check_reciprocal.py) · [Sparse quotient](checks/check_sparse.py) · [Weighted transform](checks/check_weighted.py) | C38 integration checks on the original arithmetic fixtures, including rejection of corrupted assignments. |
| [Batching](checks/check_batch.py) · [Chord geometry](checks/check_geometry.py) · [Weak recovery](checks/check_weak_recovery.py) · [Weak gate trials](checks/probe_weak_gate.py) | Finite residual fibers, signed abscissa counts, and the difference between planted and independently marked slots. |
| [Sign identity](checks/check_sign_identity.py) · [Target derivatives](checks/check_target_jets.py) · [Held-out masks](checks/probe_target_jets.py) · [Quartic decoder](checks/check_quartic_decoder.py) · [Opposite targets](checks/check_opposite_targets.py) · [Dense weighted moments](checks/check_dense_weighted.py) · [Seed rank](checks/check_seed_rank.py) | Symbolic identities, exact finite collisions, unseen-mask predictions, and reconstruction from weighted slot measurements. |
| [False-address counts](checks/check_candidate_targets.py) · [Universal offsets](checks/check_candidate_offsets.py) · [Intermediate fingerprints](checks/check_candidate_splices.py) · [Quotient arithmetic](checks/check_candidate_arithmetic.py) · [Cut and additive bounds](checks/check_candidate_fingerprint.py) | Complete fixed-target singleton patterns, joint offset tests, internal completions, and average false addresses for short fingerprints. |
| [Streamed sections](checks/check_chord_sections.py) · [Adaptive cosets](checks/check_quarter_cosets.py) · [Isogeny traces](checks/check_implicit_chord_algebra.py) · [Fixed-target poles](checks/check_implicit_chord_fiber.py) · [Masked errors](checks/check_implicit_chord_projection.py) | Section equality against literal products, subgroup coverage of actual seed states, exact trace and pole identities, and the difference between average and singleton accuracy. |
| [Product-free chord calculation](checks/check_chord_difference.py) · [Local algebra](checks/check_difference_algebra.py) · [Balanced speeds](checks/check_difference_fiber.py) · [Quadratic marker](checks/check_difference_projection.py) | Literal equality modulo $p^2$, recovery of a witness bit, zero versus multiple hits, and an ordinary incidence test for Walsh parity. |
| [Pencil searches](checks/check_section_families.py) · [Tangent readout](checks/check_section_deformation.py) · [Coarse parity score](checks/check_multiscale.py) · [Binary address bit](checks/check_orbit_parity.py) | Batched section sums, deformation characters, all-candidate Walsh guarantees and fixed-target binary ranks. |
| [Direct coefficient](checks/check_projective_coefficient.py) · [Graded kernel](checks/check_direct_algebra.py) · [Scalar join](checks/check_scalar_join.py) | Independent calculations of signed witness parity, polynomial root finding on sections, first kernel degree and intermediate matrix ranks. |
| [Fourier estimators](checks/check_direct_fiber.py) · [Sparse projections](checks/check_direct_projection.py) | Importance-sampling variance, phase-order bounds, exact interpolation and character Gram identities. |
| [Arithmetic and quadrature](checks/check_accounting.py) | Conditional constants, 8,192 exhaustive inverse-hash cases, 160 larger products, offset rounding for 103 targets, and 24 direct interpolation comparisons. |
| [Exact population](checks/check_pooled_zero.py) | All eligible small-curve seed multiplicities, pooled count certificates, planted calibration, and 756 exact marginal-likelihood comparisons. |
| [Fresh-seed trials](checks/probe_pooled_zero.py) | Independent planted training and fresh evaluation on two larger finite fields. Reproduces the raw-moment rank and error measurements. |
| [Generic packet span](checks/verify_span.py) | Literal indicator marginals, endpoint identities and modular span ranks, including full width-64 and a selected width-128 layer. |
| [Mixed-slot isolation and batching](checks/check_strided_isolation.py) | Exact counts on two fields, signed-address pinning, and the failure probability of two simultaneous field checks with one free auxiliary. |
| [Block support](checks/check_block_support.py) · [Seed cuts](checks/check_seed_cuts.py) | Independent collision-energy, complete-table and optimal hidden-slot coverage calculations on 100,352 weighted seeds, including a noncyclic curve group. |
| [Grid values and partial-sum dictionaries](checks/check_grid_values.py) | Leading mixed coefficients on liftable grids, direct elliptic zero checks, and exact finite-seed probability and energy bounds. |
| [Address parity](checks/check_address_parity.py) · [Joint determinant probe](checks/probe_signed_moment.py) | Walsh recovery controls, separated block bounds, and a comparison of the true address with uniform addresses on the same isolated seeds. |
| [Nonlinear moments](checks/probe_nonlinear_moments.py) · [Pinned signs](checks/probe_pinned_moments.py) · [Field residues](checks/probe_dense_residues.py) | Matched comparisons against other legal addresses, including independent six-slot repeats. |
| [Optimal small tables and collision energies](checks/check_nonlinear_readout.py) · [Adaptive equalities](checks/check_nonlinear_branches.py) · [Fourth moments](checks/check_character_fourth.py) | Exact finite checks of nonlinear scores, target averaging, field shifts, and a character-sum diagnostic. |
| [Signed spectrum](checks/check_signed_spectrum.py) · [Private measurement weights](checks/check_hidden_weights.py) | Exact parity recovery by Hankel reconstruction and Fermat, and the effect of refreshing independent field weights. |

The reference code enumerates small curves and tuples to establish ground truth. Its running time does not measure the asymptotic contractor. Fixed pseudorandom seeds reproduce the trials. The empirical confidence statements in the supplement concern the specified finite distributions and model classes.

## Build the paper

From the repository root, the following commands use the Python package `typst==0.15.2` to build the current manuscript and the earlier v4 paper at their published paths.

```sh
python -m pip install typst==0.15.2
python -c "import typst; typst.compile('research/fixed_shift_gram/manuscript.typ', output='research/fixed_shift_gram/manuscript.pdf')"
python -c "import typst; typst.compile('paper.typ', output='paper.pdf')"
```

With the Typst CLI, use `typst compile research/fixed_shift_gram/manuscript.typ` for the current manuscript or `typst compile paper.typ` for v4.

## Attribution and provenance

- **Fixed-shift manuscript.** Shuchang Liu, Bingsheng Zhang and Kui Ren supply the contraction theorem, arithmetic gadgets, and deterministic factor base used in the certificate. [LZR preprint](https://eprint.iacr.org/2026/3020) · [Comparison of assumptions](COMPARATIVE_AUDIT.md)
- **Public baseline for v4.** Akshima, Renas Bacho, Yu-Hsuan Huang, Eike Kiltz, Aysan Nishaburi, Samin Nooripoor, Noah Tittelbach and Emiel Wiedijk, *A Subexponential Time Algorithm for Elliptic Curve Discrete Logarithm*, v0.2, 29 October 2026. The authors give a conditional big-O bound and leave its numerical leading coefficient unspecified. [Authors’ manuscript at the cited revision](https://github.com/bochum-crypto/subexponential-ecdl/blob/1ee44e562b1b9a50ea339f3b8f81112b5a24b96f/subexponential-ecdl.pdf)
- **Supplied intermediate sources.** *Simplex Alphabets and Backward Secants for Elliptic-Curve Discrete Logarithms*, v6.0, and *Compressed Relation Search and Reusable Preprocessing for Elliptic-Curve Discrete Logarithms*, both dated 31 October 2026. The supplied simplex manuscript reports 32.1569516597 for the common project comparison class. The compressed-relation manuscript reports 27.0138774237 and states Assumption A.1. These manuscripts are credited separately from the public baseline and from this project’s contributions. The SHA-256 digests of the supplied ZIPs are, respectively, `b4ed0d2ab9bf5dae521e3e63634ccbfee8a85af234b709053c10abfe60901f6c` and `95f2078edf80cc3ce1090be7a4216fa1f407e91044d376d035c99ebf7b9204c5`.
- **Geometric background.** Qi Cheng, *Hard Problems of Algebraic Geometry Codes*, 2005. [Author’s preprint](https://arxiv.org/abs/cs/0507026)

Published by the **davidiach** repository. The synthesis, derivations and experimental work used AI assistance. The named external researchers are credited for their prior work. This is a working preprint for independent checking.
