# Does varying the target help decode pooled determinant moments?

The organizing idea from [Levent Alpöge’s remark, pp. 1 and 4](https://alpo.ge/remark.pdf) is to gain multiplicity while charging the cost of the operation that produces it. The experiment here uses only elementary finite-field algebra. Its other methodological influence is the [independently stated Lean challenge](https://github.com/Deicyde/riemann-hypothesis/blob/main/src/Challenge.lean). The statements and tests below require no Riemann hypothesis.

## What can be computed

Keep the hash seeds fixed. In the three-slot case, let $X=x(T)$ and $Q(X)=Z_2(T)$, the sum of squared determinants. In the sign-paired setting $Q$ is a polynomial of degree at most four. This degree bound and the general Gram construction are already in the project's [incidence note](target_incidence.tex). Here we test a different readout from that construction. The five coefficients of $Q(x(T)+s)$ give the value and its four higher Hasse derivatives. A Hasse derivative here means a Taylor coefficient over the field. At our fixed targets $y(T)\ne0$, the abscissa difference is also a local parameter.

For one abscissa $x_i$ in each slot write $x_4=X,\ y_4=y(T)$. Let $V_i$ be the three-variable Vandermonde product on the other three abscissas. Then the exact specification is

$$
\sum_{\epsilon\in\{\pm1\}^3}
\det[\,1,\ x,\ x^2,\ y\,]_{\epsilon,-T}^{\,2}
=8\sum_{i=1}^4 y_i^2 V_i^2 .
$$

The rows have last coordinates $\epsilon_1y_1,\epsilon_2y_2,\epsilon_3y_3,-y_4$. Expand in the last column. Each mixed product is odd in some sign and sums to zero. The remaining cofactors are signed Vandermonde products. Substituting $y_i^2=x_i^3+a x_i+c$ shows the target polynomial has degree at most four. Summing this identity over abscissa triples is valid in the full paired case. The independent symbolic check expands the determinant by permutations and the right side by polynomial multiplication over the integers. The two polynomials agree in all 76 nonzero monomials before imposing a curve equation. This check covers the identity in arbitrary characteristic by reduction. The learning claims below are separate.

For fast evaluation use the five sums $\sum_{P\in F_i}x(P)^j$ for $0\le j\le4$ and $\sum_{P\in F_i} y(P)^2$ in each slot. The Gram matrix of any subset of slots is block diagonal, with a $3\times3$ even block $B$ and odd scalar $g$. For $u=(1,X,X^2)^\top$, sum

$$
(-1)^{3-|S|}\left[g\,u^\top\operatorname{adj}(B)u+
              (X^3+aX+c)\det B\right]
$$

over the nonempty subsets $S$. The mixed determinant identity gives exactly $Q$. The readout takes $O(\sum_i |F_i|+1)$ field operations from the three slot lists, and slot construction is charged separately. It needs no tuple grid. In the general $m$-slot setting the existing even/odd Gram formula already gives the entire $Z_2(T)$ as a polynomial in $x(T)$ of degree at most $m+1$. Extracting its coefficients from the Gram matrices adds polynomial work in $m$ to $2^m$ subset operations. Thus obtaining all these derivatives costs the same asymptotic exponent as obtaining $Z_2$ in the proposed regime. Nothing in that calculation makes the derivatives a zero-count decoder.

## Exact information test

For fixed $T=(0,1)$, enumerate the true mask law, including hash-seed multiplicities. Retain all regular nonempty three-slot inputs. Each observable includes the ordered slot sizes. For each feature value $v$ the best possible yes/no table makes $\min(n_0(v),n_1(v))$ errors, where these are the two exact fresh-seed masses in counts. This is an oracle comparison within this population. It uses the labels of every possible mask triple. It is a lower bound on the mistakes any decoder of that observable can make there.

| Field and common width | Eligible seed triples | Positive | Errors with five raw moments | Errors with five Taylor coefficients | Errors with both |
|---|---:|---:|---:|---:|---:|
| 13, width 2 | 367,104 | 258,288 | 1,200 | 0 | 0 |
| 13, width 4 | 80,064 | 71,040 | 0 | 0 | 0 |
| 23, width 2 | 25,884,768 | 11,769,600 | 43,584 | 3,840 | 0 |

For the last row, the subpopulation with at least three sign orbits had mass 3,951,360. The corresponding error counts were 2,448, 768 and zero. The direct determinant reference evaluated the polynomial at every affine target on each curve. The five coefficients computed from the slot-level Gram formula also agreed with explicit tuple sums for all ordered mask triples, including masks outside the scoring population.

## An independent test of learnability

At each field, fit bin decisions using 30,000 planted samples and 30,000 independent, unlabeled fresh samples. A bin is called positive when $2\mu\,\widehat Q(v)>\widehat P(v)$, with $\mu=\mathbb E H$ and equal sample sizes. This score uses the planted count bias and can differ from the Bayes existence score when multiple hits are common. Test on 20,000 newly drawn fresh seed triples including empty and irregular masks. The decoder defaults to absence on those cases. Calibration uses sorted slot sizes and the five polynomial coefficients, an invertible coordinate change from the Taylor coefficients at a fixed target. The comparator memorizes the unordered triple of complete masks using exactly the same draws, score, and fallback.

| Field, width | Fresh positives | Taylor lookup errors | Mask memorization errors | True positives found by Taylor lookup on previously unseen eligible triples |
|---|---:|---:|---:|---:|
| 23, width 2 | 3,165 | 1,976 | 1,975 | 0 of 254 |
| 101, width 4 | 4,763 | 4,744 | 4,744 | 0 of 3,284 |

An unseen mask triple is one whose unordered triple of exact retained sets was absent from both training samples. There were 756 such eligible test draws over the field of order 23 and 9,540 over the field of order 101. This control exposes why a nearly perfect finite lookup table alone can be misleading. A transferable arithmetic function of the coefficients could behave differently.

## Multiplicity and the cost of an auxiliary function

There is a useful test for proposed multiplicity constructions. Fix the seeds on a prime-order curve, then choose one nonzero rational function on the target curve independently of which target will be queried. If its pole divisor has degree $D$, the number of targets where it vanishes to order at least $s$ is at most $\lfloor D/s\rfloor$. This follows from equality of zero and pole divisor degrees, including multiplicities. For uniform nonidentity targets the accepted mass is consequently at most $D/(s(r-1))$, with $r$ the group order. A scheme that accepts only on such a zero and achieves fresh error $O(\lambda^2)$ for every prescribed target must have $\mathbb E D\ge s(r-1)(\lambda-O(\lambda^2))$, since the averaged true positive mass is $\lambda+O(\lambda^2)$. The same accounting applies to a union of such tests by adding their degree budgets divided by their orders. This prices pole degree for this particular representation. A high-degree function may have a short evaluation circuit, and arbitrary nonlinear readouts of the target coefficients fall outside the argument.

The unresolved target is a structural arithmetic rule that recognizes relations on unseen masks while keeping fresh false positives within the error budget. No end-to-end exponent improvement follows from these experiments.

From the project root run the following with Python, NumPy and SymPy.

~~~sh
python -B checks/check_sign_identity.py
python -B checks/check_target_jets.py
python -B checks/probe_target_jets.py
~~~
