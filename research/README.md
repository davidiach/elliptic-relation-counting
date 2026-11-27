# Research guide

The [fixed-shift paper](fixed_shift_gram/manuscript.pdf) develops the matrix certificate. The notes below study ways of choosing enough equations deterministically. For the scope of the contraction hypothesis and the bit model, start with the [main theorem](../MAIN_RESULT.md) and [complexity audit](../COMPLEXITY_AUDIT.md). The earlier hashed-slot paper has a [separate guide](../derivations/v4-guide.md).

## Fixed-shift equations and prescribed fibers

| Question | Derivation and evidence |
| --- | --- |
| How much rank survives evaluation? | The [evaluation lists](../DETERMINISTIC_POINTS.md) and [neighboring-fiber identities](../NEIGHBORING_FIBERS.md) explain full-grid rank, sparse polynomial consequences and cancellation. The [partial expansion](../PARTIAL_EXPANSION.md) prices the residual contraction. [Independent checks of the argument](../ADVERSARIAL_REVIEW.md) · [Comparison with the v4 section norm](../COMPATIBILITY_MAP.md) |
| Can shorter relations or integer lifts supply independent equations? | [Conditioned means and adaptive rank](short_relations_adaptive_rank/RESEARCH.md), [integral relation lattices](integral_lattice_audit/RESEARCH.md), and [pinning and residual obstructions](deterministic_endgame_barriers/RESEARCH.md) give their hypotheses and counterexamples. |
| What can one learn at a prescribed sum? | [Pointwise obstruction](pointwise_target_obstruction/RESEARCH.md), [two-representation certificates](prescribed_fiber_rigidity/MAIN_RESULT.md), [simultaneous failure of two sufficient criteria](dual_lift_obstruction/RESEARCH.md), and [enumerable prescribed fibers](enumerable_prescribed_fibers/RESEARCH.md). The last result uses a modified seed. |
| How do pinned ranks behave across consecutive targets? | The [reduced-floor examples](sparse_target_runs/RESEARCH.md) compare raw incidence rank, pinned rank and stacked nullity on prime-order curves. [Six-target examples](sparse_target_runs/SIX_TARGETS.md) identify a nonconstant vector shared by all six pinned matrices at orders 109 and 271. [Seeds](sparse_target_runs/witnesses.json) · [Transcript](sparse_target_runs/verification.txt) |

## Additive avoidance

The quotient argument maps the remaining rank question into the plane over the scalar field. When the product of the list sizes is less than the number of plane points, their sumset misses a point. Finding a short list guaranteed to contain such a point would suffice. [Quotient compression](universal_quotient_compression/RESEARCH.md) · [Universal probe lists](additive_probe_lists/RESEARCH.md)

| Result | Scope |
| --- | --- |
| [Quotient-sumset completion](quotient_sumset_completion/RESEARCH.md) and [compact covers](compact_quotient_covers/RESEARCH.md) | Sufficient conditions based on missing projected values, capped supports or explicit covers by affine flats and boxes. |
| [Common translations](translation_cover/RESEARCH.md) | A constant-size translating block gives full coverage and a strict product bound under the closed-contraction hypothesis. [Finite checks](translation_cover/verification.txt) |
| [Avoidance bounds](additive_avoidance_boundaries/RESEARCH.md) and [constructive probes](additive_probe_constructions/RESEARCH.md) | Bounds on prescribed sets that input lists can cover, with explicit constructions under field-size or growth conditions. |
| [Higher-order probes](high_growth_additive_probes/RESEARCH.md) | For two factors of size at most $n$, the construction uses $n^{1+o(1)}$ points when $n\to\infty$ and $n=r^{o(1)}$. The note also constructs high-growth inputs that cover prescribed points. [Finite checks](high_growth_additive_probes/verification.txt) |
| [Carry-aware digit lift](digit_lift_avoidance/RESEARCH.md) and [structured remote syndromes](structured_remote_points/RESEARCH.md) | The lift accounts for carries in digit coordinates. Exact constructions apply to hidden bounded-rank dictionaries or supplied trees with small linear interfaces. The general high-rate case remains open. [Digit checks](digit_lift_avoidance/verification.txt) · [Structured examples](structured_remote_points/verification.txt) |

## What moments determine

For an abstract population of integer rows, sums of coordinate monomials can agree even when reductions of the rows have different ranks. [Moment-oracle counterexamples](exact_moment_oracle_barriers/RESEARCH.md) · [Support and degree bounds](exact_moment_oracle_barriers/moment_rank_boundary/RESEARCH.md)

[Newton certificates](exact_moment_oracle_barriers/moment_rank_boundary/newton_sparse_rank/RESEARCH.md) and [balanced incidence moments](exact_moment_oracle_barriers/moment_rank_boundary/balanced_incidence_moments/RESEARCH.md) give sharp degrees for specified alphabets. For one balanced signed-unit block of positive mass at most $w$ on at least $3w$ labels, degree $3w$ is necessary and sufficient for reconstruction. Sharpness for rank separation with fewer rows than the modulus is a separate question at weights four and five. [Recorded finite checks](exact_moment_oracle_barriers/moment_rank_boundary/balanced_incidence_moments/verification.txt)

## Reproducibility

The [README commands](../README.md#reproduce-the-experiments) cover the core paper. Each research note specifies its own verifier, working directory and recorded output. The [ledger](../RESEARCH_LEDGER.md#reproduction) includes the separate audit of public theorem declarations. The [v4 guide](../derivations/v4-guide.md#reproduce-the-v4-checks) identifies the checks that require the original source archives.
