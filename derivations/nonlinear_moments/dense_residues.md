# Dense readouts of the raw first and second moments

This experiment uses the actual residues of the joint signed moments $M_1,M_2\in\mathbb F_p$. The determinant and mixed-discriminant formulas make both moments affordable in one challenge at the arities of interest. We tested whether the canonical ordering of their residues, field division, and a learned two-dimensional table could identify the witness among the other legal signed tuples on the same root. [Moment formulas and paired test design](probe.md)

## Readouts

Let $C(0)=0$, $C(z)=1$ for representatives $1\le z\le (p-1)/2$, and $C(z)=-1$ otherwise. The fixed tests use $C(M_2+M_1^2)$, $C(M_2+M_1^2+\lfloor p/4\rfloor)$ and $C(M_2/(M_1^2-d))$, with arithmetic modulo $p$. A public quadratic nonresidue $d$ ensures the denominator is nonzero for every challenge. The finite fixture uses the first nonresidue and a table of inverses. A general query needs only one inverse, and a suitable public $d$ can be selected by random sampling with a field exponentiation check in expected polynomial bit time.

The learned test reserves a separate bin for a zero moment and partitions the nonzero canonical representatives into eight contiguous bins of almost equal width. The pair of bin indices selects one of 81 entries. Two calibration objectives maximize signed contrast or seek higher squared magnitude at the true address. The latter uses eight public-seed-reproducible starting points for coordinate ascent. A separate interaction score replaces the fitted power-table entry $\theta_{ij}$ by

$$
(\theta_{ij}-\theta_{i0}-\theta_{0j}+\theta_{00})/4 .
$$

This remains bounded by one, removes the contributions of arbitrary additive one-coordinate functions, and vanishes whenever either moment is zero. The full joint table and all three fixed scores are evaluated separately. With a given 81-entry table and $d$, every readout costs $(w+2^m)\operatorname{poly}(b,m)$ expected bit operations per challenge, including the moments. The code uses full spectra and truth labels to train and measure the finite models. No asymptotic calibration or recovery bound follows from their constant description size.

## Independent paired comparisons

For each field, arity and fixed target, 1,200 conditional roots trained the table and 1,600 independently generated roots evaluated it. Every root in these experiments has a unique ordinary tuple that passes the addition guards. Each legal non-witness address on that same root contributes equally to its control. The held-out streams differ from those used for the quadratic-character experiments. In the table, rank means the probability that the marked coefficient has larger magnitude than a uniformly chosen wrong legal coefficient, counting ties half. SE is the sample standard error over roots.

| $p$ | $m$ | Target | Both raw moments nonzero | Power-table rank (SE) | Interaction rank (SE) | Interaction magnitude gap (SE) |
|---:|---:|---|---:|---:|---:|---:|
| 1,031 | 4 | $G$ | 61.90% | 0.4904 (0.0073) | 0.5055 (0.0071) | +0.00003 (0.00033) |
| 1,031 | 4 | $2G$ | 61.93% | 0.4832 (0.0072) | 0.4975 (0.0071) | −0.00010 (0.00033) |
| 65,521 | 6 | $G$ | 54.58% | 0.5049 (0.0071) | 0.4947 (0.0072) | −0.00005 (0.00026) |
| 65,521 | 6 | $2G$ | 54.59% | 0.4929 (0.0072) | 0.4985 (0.0072) | +0.00010 (0.00022) |

The held-out challenges populated all 81 cells at $p=1031$ and 79 cells at $p=65521$. For the fixed ratio readout the ranks were $0.5089$, $0.5067$, $0.5014$ and $0.4887$ in the same row order. The fixed shifted and unshifted additive readouts gave ranks between $0.4915$ and $0.5060$. The measurements therefore include variation across the raw residues even when both moments are nonzero. They give no replicated improvement across arities and targets.

For scale, Hoeffding's inequality and a union bound over the four interaction ranks give simultaneous 95% intervals of half-width $\sqrt{\ln(160)/3200}<0.040$, conditional on the fitted tables. These bounds concern these four finite comparisons. A small correlation everywhere, or a different nonlinear readout that concentrates on larger arities, remains possible. [Executable test](../../checks/probe_dense_residues.py) · [Recorded output for all six readouts](../../evidence/probe_dense_residues.txt)
