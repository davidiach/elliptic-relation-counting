# Elliptic-curve relation counting

How much information about a discrete logarithm survives when many elliptic-curve relations are compressed into a single matrix? This repository contains two research reports on that question, their mathematical derivations, and Python programs that check the results on finite examples.

## Research results

**[A polynomial certificate for deterministic elliptic-curve relation systems](research/fixed_shift_gram/manuscript.pdf)** · 27 November 2026

Summing relations can lose rank. The fixed-shift paper assigns a polynomial weight to each relation and proves that the resulting matrix has a nonzero determinant. For a nonsingular elliptic curve over a prime field of characteristic greater than three whose rational points form a group of prime order, a random evaluation yields a discrete-log algorithm with expected bit time $2^{O(b/(\log_2 b)^{1/4})}$, where $b$ is the field bit length. The bound assumes the closed Boolean tensor-contraction theorem of Liu, Zhang and Ren. Finding suitable evaluations deterministically within the same time bound for every input remains open. [Theorem and assumptions](MAIN_RESULT.md) · [Source](research/fixed_shift_gram/manuscript.typ) · [Open questions](OPEN_PROBLEMS.md)

Akshima et al. previously obtained the same conditional expected-time scale using randomly hashed factor bases. The fixed-shift paper uses a deterministic relation system and a single shared evaluation, with $O(b)$ expected fair bits for the whole computation. [Akshima et al., v0.2](https://github.com/bochum-crypto/subexponential-ecdl/blob/1ee44e562b1b9a50ea339f3b8f81112b5a24b96f/subexponential-ecdl.pdf) · [Weighted-count audit](COMPARATIVE_AUDIT.md#weighted-counts-against-the-closed-theorem)

### Earlier paper and version history

**[Counting relations in independently hashed elliptic-curve slots](paper.pdf)** · Version 4 · 24 November 2026

The earlier paper studies relations whose slots are sampled independently. For its field-affine sampler, the leading coefficients on $b/(\log_2 b)^{1/4}$ are 11.31486814 for log expected bit time and 5.65743407 for log space in the stated input class. These bounds assume the finite-packet interface A.1 and the inherited reduction estimates. [Source](paper.typ) · [Derivations and version history](derivations/v4-guide.md)

The companion note *[Target certificates from prescribed sum fibers in a group of prime order](research/prescribed_fiber_rigidity/manuscript.pdf)* studies what can be certified from representations of a single target. The [research index](research/README.md) connects this and other follow-up questions to their proofs and examples.

## Attribution and provenance

[David Iach](https://github.com/davidiach) is the project lead and repository maintainer. AI agents helped develop mathematical arguments and computational experiments for both reports. OpenAI Codex edited the manuscripts and reran the checks described below. The project has no recorded independent human verification of the general proofs. [Questions and corrections](https://github.com/davidiach/elliptic-relation-counting/issues)

The fixed-shift construction uses the contraction theorem, arithmetic gadgets, and factor base of Shuchang Liu, Bingsheng Zhang and Kui Ren. [Liu–Zhang–Ren preprint](https://eprint.iacr.org/2026/3020) · [Comparison of assumptions](COMPARATIVE_AUDIT.md)

The v4 analysis builds on the public discrete-logarithm preprint by Akshima and coauthors, two subsequent working manuscripts supplied to the project, and Qi Cheng’s elliptic evaluation-code construction. The v4 guide identifies which source supplies each assumption and numerical baseline. [Sources and archive hashes](derivations/v4-guide.md#sources-and-archive-hashes)

## Reproduce the experiments

From the repository root, use Python 3.12 with assertions enabled. To reproduce the fixed-shift results, run

```sh
python -m pip install -r requirements.txt
python -B checks/check_fixed_shift_gram.py
python -B checks/check_fixed_shift_points.py
python -B checks/check_neighboring_fibers.py
python -B -m research.sparse_target_runs.check_four_targets
```

The first three outputs matched the recorded counts and ranks. The fourth reproduced every row of the paper’s sparse-target table. [Gram matrices](evidence/fixed_shift_gram.txt) · [Evaluation lists](evidence/fixed_shift_points.txt) · [Neighboring fibers](evidence/neighboring_fibers.txt) · [Sparse targets](research/sparse_target_runs/verification.txt)

For v4, `python -B reproduce.py` runs 52 standalone checks. With the original archives, `python -B reproduce.py --source-archives DIR` runs all 56, including four C38 integration checks. All 56 outputs matched their transcripts. A separate section-norm audit also matched its transcript. [v4 checks and inputs](derivations/v4-guide.md#reproduce-the-v4-checks) · [Section-norm audit](archive/research/chat_2026_11_23_24/08_section_norm/section_norm_audit_output.txt)

For the companion note, run `verify.py`, `verify_restrictions.py`, and `verify_height.py` with Python in the [note’s directory](research/prescribed_fiber_rigidity/). Each output matched the corresponding transcript there.

## Build the paper

From the repository root, build the two main reports with Typst 0.15.2.

```sh
typst compile research/fixed_shift_gram/manuscript.typ
typst compile paper.typ
```

With Python, install `typst==0.15.2` and call `typst.compile(source, output=pdf)` with the corresponding `.typ` and `.pdf` paths.

The [companion note](research/prescribed_fiber_rigidity/manuscript.tex) uses LaTeX. Build it with `tectonic research/prescribed_fiber_rigidity/manuscript.tex` from the repository root.

## How to cite

The [BibTeX citation file](CITATION.bib) gives a title-first reference and a permanent PDF link for each report and the companion note. For example, the citation for v4 is

> *Counting relations in independently hashed elliptic-curve slots*. Version 4. Elliptic Relation Counting Project, 24 November 2026.

When citing a computation, also give the script path and the Git commit you ran.
