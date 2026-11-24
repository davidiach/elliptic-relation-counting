# General-prime ECDLP research checkpoints, 23–24 November 2026

This is a **historical AI-assisted exploratory research archive**, not a revised paper or a claim of a new ECDLP algorithm. It preserves the notes, source files, and recorded experimental outputs of the conversation in chronological order. The original version 3 paper and its established conditional comparison remain unchanged.

## Status and review boundary

- Paper v3 reports a conditional general-prime mixed-sampler time coefficient 24.00246 on the scale `b/(log2 b)^(1/4)`. The conversation subsequently explored conditional coefficient candidates 13.85782672, 12.39481703, and 11.31486814. These candidates have **not received a complete independent proof-chain audit**.
- No strict asymptotic improvement over `b/(log2 b)^(1/4)` is established. Square-root-in-bitlength and other improved exponents in the notes are targets or conditional transfers, not demonstrated algorithms.
- The current frontier is evaluating a collective norm/trace without expanding the remaining Cartesian grid. Local isogeny, resultant, and residue computations do not resolve that cost.
- Claims, assumptions, finite regressions, and unsuccessful experiments are preserved as originally recorded. They may contain gaps or errors; inclusion here is not endorsement. Numerical regressions are not cryptographic-scale performance tests. No independent novelty or priority claim is made.

## Chronological index

| Checkpoint | Original research note | Archived files |
|---|---|---:|
| 00 | [Principal-ideal polynomial-ring compiler (initial verifier)](00_polynomial_ring/polynomial_ring_verifier.py) | 2 |
| 01 | [General-prime coefficient candidate: abscissa-only principal-ideal compiler](01_abscissa_chain/abscissa_chain_note.md) | 5 |
| 02 | [Large-arity relation search: short ideal projectors, shifted digit blocks, and a taper barrier](02_large_arity/large_arity_research_note.md) | 5 |
| 03 | [Field-affine thin-digit samplers and a planted polynomial-readout obstruction](03_affine_digit_sampler/affine_digit_sampler_note.md) | 3 |
| 04 | [A degree-sensitive tensor contractor and two barriers to its ECDLP transfer](04_contractor_degree/contractor_degree_note.md) | 3 |
| 05 | [A slice-rank obstruction to low-degree relation compilers](05_seed_local_slice_rank/seed_local_slice_rank_note.md) | 4 |
| 06 | [Average-case weak relation detection: recovery, robust slice rank, and character readouts](06_average_case_detector/average_case_detector_note.md) | 5 |
| 07 | [Normalized chord residuals: exact box cancellation, dense odd support, and the first informative sketch](07_modular_cancellation/modular_cancellation_note.md) | 5 |
| 08 | [Global elliptic section norms and pairwise compression of dyadic modes](08_section_norm/section_norm_research_note.md) | 7 |
| 09 | [High-precision lifting, degree-aware mode contraction, and the remaining packing gap](09_mode_contraction/mode_contraction_note.md) | 3 |
| 10 | [Symbolic group splitting, fixed-budget moment ranks, and unit-separated XOR packets](10_group_tensor_factorization/group_tensor_factorization_note.md) | 3 |
| 11 | [Hierarchical packet packing and shared-mode Newton contraction](11_hierarchical_shared_mode/hierarchical_shared_mode_note.md) | 5 |
| 12 | [Exterior-power packet deformation and a support-sensitive contractor](12_exterior_packet/exterior_packet_note.md) | 7 |
| 13 | [Direct elliptic evaluation-code route: transfer target and an oracle obstruction](13_prime_exponent/prime_exponent_note.md) | 3 |
| 14 | [Adaptive numerical sketches and full-value moment smoothing on elliptic planted lists](14_numerical_sketch/numerical_sketch_note.md) | 3 |
| 15 | [Exact spectral normal form and higher-order incidence for field-affine elliptic samplers](15_affine_seed_spectrum/affine_seed_spectrum_note.md) | 3 |
| 16 | [Uniform replica mixing and Poisson relation counts for translated elliptic slots](16_replica_mixing/replica_mixing_note.md) | 3 |
| 17 | [Correlated translation samplers, cofactor transfer, and inverse relation queries](17_shared_translation/shared_translation_note.md) | 5 |
| 18 | [Primitive shared-translation relation polynomials, collective resultants, and exact inverse planting](18_shared_norm_aggregate/shared_norm_aggregate_note.md) | 9 |
| 19 | [Block norm elimination, one-call witness extraction, and the grid lifting ideal](19_kronecker_determinant/kronecker_determinant_note.md) | 5 |
| 20 | [Fiber-gcd compression, trace-only witness recovery, and sharp image-certificate degrees](20_fiber_gcd/fiber_gcd_note.md) | 3 |
| 21 | [Projected norms, batched rational traces, and a resolvent readout obstruction](21_collective_trace/collective_trace_note.md) | 5 |
| 22 | [Mixed-characteristic residue extraction and collision-free scalar norm recovery](22_mixed_characteristic_norm/mixed_characteristic_norm_note.md) | 5 |
| 23 | [Coordinate-deformed elliptic resultants: signed recovery and a prime-power evaluator](23_coordinate_deformed_resultant/coordinate_deformed_resultant_note.md) | 6 |
| 24 | [Auxiliary elliptic torsion supports and compact symbolic fiber elimination](24_auxiliary_torsion/auxiliary_torsion_note.md) | 11 |
| 25 | [Smooth auxiliary isogeny towers: exact scalar norms and rational traces](25_smooth_isogeny_norm/smooth_isogeny_norm_note.md) | 4 |

## Integrity and reproduction

`MANIFEST.json` records the SHA-256 of every original ZIP and its individual members, their exact original filenames, and corresponding Git blob IDs. The ZIP members are stored unpacked, byte-for-byte, in their checkpoint directories; original binary ZIP containers are intentionally not duplicated in Git. The initial stand-alone verifier is included separately. The supplied baseline PDF is already the repository's `paper.pdf` (identical Git blob ID).

```sh
python -B archive/research/chat_2026_11_23_24/verify_archive.py
python -m pip install -r requirements.txt
python -B archive/research/chat_2026_11_23_24/25_smooth_isogeny_norm/verify.py
```

The integrity verifier checks member hashes and Python syntax without executing historical code. Checkpoint notes describe individual tests and recorded transcripts. Run historical scripts from their own checkpoint directory, using Python 3.11+ and the NumPy/SymPy dependencies used by individual scripts. The existing repository-wide `reproduce.py` and established tests are unchanged. The whole historical suite is not claimed to have been independently rerun at import time.
