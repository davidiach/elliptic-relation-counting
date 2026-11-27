# Elliptic-curve relation counting

**[Read the fixed-shift paper (PDF)](research/fixed_shift_gram/manuscript.pdf)** · [Typst source](research/fixed_shift_gram/manuscript.typ) · [Theorem and assumptions](MAIN_RESULT.md) · [Research index](research/README.md)

Can a family of elliptic-curve relations too large to enumerate be used to recover a discrete logarithm? The current preprint, *A polynomial certificate for deterministic elliptic-curve relation systems*, encodes the relations in a matrix whose determinant is a nonzero polynomial. A full-rank evaluation determines the logarithm, and scalar multiplication verifies the answer.

For a nonsingular curve over a prime field of characteristic greater than three whose rational points form a group of prime order, assume the **closed Boolean tensor-contraction bound of Liu, Zhang and Ren (Theorem 2.2)**. With field bit length $b$, the certificate gives expected bit time $2^{O(b/(\log_2 b)^{1/4})}$ and $O(b)$ expected random bits. A deterministic algorithm attains the same bound with $O(b)$ bits of advice per length. Finding suitable evaluations deterministically within that budget for every input remains open. [Bit complexity](COMPLEXITY_AUDIT.md) · [Open questions](OPEN_PROBLEMS.md)

## Research results

| Question | What the repository establishes |
| --- | --- |
| Which evaluations preserve rank? | The full field grid suffices. Under the short-list hypotheses, keeping the source group value and one pin on each side preserves shift rank. Combining those classes can cancel a direction. [Deterministic points](DETERMINISTIC_POINTS.md) · [Neighboring fibers](NEIGHBORING_FIBERS.md) |
| Can a few pinned targets certify rank? | In reduced-floor affine examples of orders 223, 67 and 109, four or five consecutive targets each have excess nullity. Those three stacks recover the forced nullity. At orders 109 and 271, different seeds give six consecutive targets sharing a nonconstant kernel vector. [Four- and five-target data](research/sparse_target_runs/RESEARCH.md) · [Six-target witnesses](research/sparse_target_runs/SIX_TARGETS.md) |
| Can additive probes complete the system? | Common translates give full coverage with fewer than $r^2$ tuples, where $r$ is the group order. Small universal probe lists exist under fixed product slack. Explicit constructions cover restricted inputs, including two small factors. A uniform construction at growing arity remains open. [Translation covers](research/translation_cover/RESEARCH.md) · [Additive probes and structured cases](research/README.md#additive-avoidance) |
| What do the exact moments determine? | Rank separation and reconstruction have different degree thresholds. The results specify the block alphabet, row multiplicities and modulus for each threshold. [Moment bounds and examples](research/README.md#what-moments-determine) |

The [research index](research/README.md) links each result to its derivation and finite checks. The [comparative audit](COMPARATIVE_AUDIT.md) traces the external contraction premise. The [research ledger](RESEARCH_LEDGER.md) records how the arguments and examples developed.

## Reproduce the experiments

From the repository root, use Python 3.12 with assertions enabled. The three manuscript checks exercise the finite rank identities, cancellation examples and deterministic evaluation lists.

```sh
python -m pip install -r requirements.txt
python -B checks/check_fixed_shift_gram.py
python -B checks/check_fixed_shift_points.py
python -B checks/check_neighboring_fibers.py
```

Compare stdout with the [Gram](evidence/fixed_shift_gram.txt), [point-list](evidence/fixed_shift_points.txt) and [neighboring-fiber](evidence/neighboring_fibers.txt) transcripts. To reproduce the sparse-target table, run `python -B -m research.sparse_target_runs.check_four_targets` and compare the [recorded output](research/sparse_target_runs/verification.txt). These scripts use small instances with independently countable relations. The asymptotic time bound comes from the conditional proof.

## Build the paper

From the repository root with the Typst CLI, run `typst compile research/fixed_shift_gram/manuscript.typ`. The Python build uses the same source and output paths.

```sh
python -m pip install typst==0.15.2
python -c "import typst; typst.compile('research/fixed_shift_gram/manuscript.typ', output='research/fixed_shift_gram/manuscript.pdf')"
```

## Earlier paper and version history

*Counting relations in independently hashed elliptic-curve slots* (v4, 24 November 2026) derives numerical coefficients under the separate finite-packet Assumption A.1 and the inherited reduction estimates. For its general-prime field-affine sampler, the leading coefficients on $b/(\log_2 b)^{1/4}$ are 11.31486814 for log expected bit time and 5.65743407 for log space. The v4 guide includes the release history, all 52 default checks, and the four optional C38 checks with their archive requirements. [v4 paper](paper.pdf) · [Typst source](paper.typ) · [v4 guide](derivations/v4-guide.md) · [Archived research](archive/README.md)

## Attribution and provenance

- **Fixed-shift manuscript.** Shuchang Liu, Bingsheng Zhang and Kui Ren supply the contraction theorem, arithmetic gadgets, and deterministic factor base used in the certificate. [LZR preprint](https://eprint.iacr.org/2026/3020) · [Comparison of assumptions](COMPARATIVE_AUDIT.md)
- **Public baseline for v4.** Akshima, Renas Bacho, Yu-Hsuan Huang, Eike Kiltz, Aysan Nishaburi, Samin Nooripoor, Noah Tittelbach and Emiel Wiedijk, *A Subexponential Time Algorithm for Elliptic Curve Discrete Logarithm*, v0.2, 29 October 2026. The authors give a conditional big-O bound and leave its numerical leading coefficient unspecified. [Authors’ manuscript at the cited revision](https://github.com/bochum-crypto/subexponential-ecdl/blob/1ee44e562b1b9a50ea339f3b8f81112b5a24b96f/subexponential-ecdl.pdf)
- **Supplied intermediate sources.** *Simplex Alphabets and Backward Secants for Elliptic-Curve Discrete Logarithms*, v6.0, and *Compressed Relation Search and Reusable Preprocessing for Elliptic-Curve Discrete Logarithms*, both dated 31 October 2026. The supplied simplex manuscript reports 32.1569516597 for the common project comparison class. The compressed-relation manuscript reports 27.0138774237 and states Assumption A.1. These manuscripts are credited separately from the public baseline and from this project’s contributions. The SHA-256 digests of the supplied ZIPs are, respectively, `b4ed0d2ab9bf5dae521e3e63634ccbfee8a85af234b709053c10abfe60901f6c` and `95f2078edf80cc3ce1090be7a4216fa1f407e91044d376d035c99ebf7b9204c5`.
- **Geometric background.** Qi Cheng, *Hard Problems of Algebraic Geometry Codes*, 2005. [Author’s preprint](https://arxiv.org/abs/cs/0507026)

**Project lead and repository maintainer** · [David Iach (@davidiach)](https://github.com/davidiach). The synthesis, derivations and experimental work used AI assistance. The named external researchers are credited for their prior work. This is a working preprint for independent checking.
