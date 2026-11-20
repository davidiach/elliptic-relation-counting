# Counting relations in independently hashed elliptic-curve slots

[Read the paper (PDF)](paper.pdf) · [Editable Typst source](paper.typ) · [Index of mathematical derivations](derivations/README.md) · [Recorded experiments](evidence/)

This AI-assisted research preprint studies a fixed-target counting problem behind elliptic-curve relation search. Each slot retains points using an independent hash of the abscissa. How much computation, and how much joint information about the slots, is needed to decide whether any tuple sums to a prescribed target?

## Results

| Question | Result | Scope |
|---|---|---|
| Can the hash be encoded more cheaply? | Checking the short inverse image reduces the leading arithmetic word budget from 38 to 36 per slot per field bit. The time coefficient falls from 27.0138774237 to 26.2933775156, a 2.6671% reduction. | Conditional on the supplied finite-packet Assumption A.1 and inherited graph and relation-reduction estimates. The scale is $b/(\log_2 b)^{1/4}$ in the log of expected bit time. |
| Is evidence from almost every slot enough? | We derive the exact marginal likelihood under planting. Any adaptive predictor that omits a whole slot has fresh error at least $\lambda-O(\lambda^2)$ at the proposed balanced widths and fixed target. | An access bound for this seed law. Joint moments depend on every slot. |
| Can joint moments be computed on a smaller grid? | Elliptic interpolation computes moments up to order $\lfloor 2^{b/m^2}\rfloor$ within expected log cost $(1+o(1))b/m$ at the proposed growing arity. | Computing a statistic and decoding its zero count have separate costs. |
| Can the moments reveal the exact count? | Nine pooled moments certify 145,776 of 168,480 eligible seed triples over the field of order 13. Planted calibration cuts existence errors from 13,632 to 7,680. | Regular three-slot inputs with two sign orbits. On the larger field of order 1,031 the tested rules missed all 2,929 fresh positives having more than two orbits and regular masks. |

The paper also treats guarded-to-ordinary count conversion, fixed-target isolation, generic packet span, optimal separated approximation, derivative readout, and the time–space tradeoff for implicit resultants. The square-root-in-bitlength objective and a cheaper decoder for typical large pools remain open. These comparisons concern asymptotic formulas and finite counting experiments. This repository has no demonstrated end-to-end cryptographic-scale speedup.

## Reproduce the experiments

The recorded runs used Python 3.12.2. In a Python 3.12 environment run the following from the repository root.

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

The runner checks all four outputs against the committed transcripts. The finite-rank check is usually the slowest. Each script can also be run independently.

| Script | What it checks |
|---|---|
| [Arithmetic and quadrature](checks/check_accounting.py) | Conditional constants, 8,192 exhaustive inverse-hash cases, 160 larger products, offset rounding for 103 targets, and 24 direct interpolation comparisons. |
| [Exact population](checks/check_pooled_zero.py) | All eligible small-curve seed multiplicities, pooled count certificates, planted calibration, and 756 exact marginal-likelihood comparisons. |
| [Fresh-seed trials](checks/probe_pooled_zero.py) | Independent planted training and fresh evaluation on two larger finite fields. Reproduces the raw-moment rank and error measurements. |
| [Generic packet span](checks/verify_span.py) | Literal indicator marginals, endpoint identities and modular span ranks, including full width-64 and a selected width-128 layer. |

The reference code enumerates small curves and tuples to establish ground truth. Its running time does not measure the asymptotic contractor. Fixed pseudorandom seeds reproduce the trials. The empirical confidence statements in the supplement concern the specified finite distributions and model classes.

## Build the paper

The PDF was built with Typst 0.15.2. With the Typst CLI, run `typst compile paper.typ paper.pdf`. Alternatively, install the Python package `typst==0.15.2` and run `python -c "import typst; typst.compile('paper.typ', output='paper.pdf')"`.

## Attribution and provenance

- **Public baseline.** Akshima, Renas Bacho, Yu-Hsuan Huang, Eike Kiltz, Aysan Nishaburi, Samin Nooripoor, Noah Tittelbach and Emiel Wiedijk, *A Subexponential Time Algorithm for Elliptic Curve Discrete Logarithm*, v0.2, 29 October 2026. The authors give a conditional big-O bound and leave its numerical leading coefficient unspecified. [Authors’ manuscript at the cited revision](https://github.com/bochum-crypto/subexponential-ecdl/blob/1ee44e562b1b9a50ea339f3b8f81112b5a24b96f/subexponential-ecdl.pdf)
- **Supplied intermediate sources.** *Simplex Alphabets and Backward Secants for Elliptic-Curve Discrete Logarithms*, v6.0, and *Compressed Relation Search and Reusable Preprocessing for Elliptic-Curve Discrete Logarithms*, both dated 31 October 2026. The supplied simplex manuscript reports 32.1569516597 for the common project comparison class. The compressed-relation manuscript reports 27.0138774237 and states Assumption A.1. These manuscripts are credited separately from the public baseline and from this project’s contributions. The SHA-256 digests of the supplied ZIPs are, respectively, `b4ed0d2ab9bf5dae521e3e63634ccbfee8a85af234b709053c10abfe60901f6c` and `95f2078edf80cc3ce1090be7a4216fa1f407e91044d376d035c99ebf7b9204c5`.
- **Geometric background.** Qi Cheng, *Hard Problems of Algebraic Geometry Codes*, 2005. [Author’s preprint](https://arxiv.org/abs/cs/0507026)

Published by the **davidiach** repository. The synthesis, derivations and experimental work used AI assistance. The named external researchers are credited for their prior work. This is a working preprint for independent checking.
