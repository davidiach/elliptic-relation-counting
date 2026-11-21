# A smaller conditional bound from fixed abscissa bits

Let $b=\lceil\log_2 p\rceil$ and $Z_b=b/(\log_2 b)^{1/4}$. On the input class of the supplied compressed relation-search manuscript, under its finite-packet Assumption A.1 and the same contraction and verified-solver estimates, we obtain

$$
\log_2 \mathbb E T\le (24.002459963814+o(1))Z_b,\qquad
\log_2 \mathrm{space}\le(12.001229981907+o(1))Z_b .
$$

The time coefficient in the project's previous construction was 26.293377515591. This is an 8.7129% decrease. The change requires a new sampling distribution for all but one slot. The bound for a caller who supplies independent affine hashes in every slot remains open to further improvement. The square-root-in-bitlength goal also remains open. [Prior coefficient and assumptions](../short_hash_count/count.tex) · [Source for the one-pilot theorem and C38 ledger](../../README.md#attribution-and-provenance)

## Sampling and target coverage

Work on a smooth short Weierstrass curve in characteristic $p>3$. Write $N=\#E(\mathbb F_p)$, $M=2^b$, $w=2^\ell$, $\rho=w/M$, and $\lambda=N^{m-1}\rho^m$. Keep the ordinary odd affine hash of width $w$ in slot 1. For each of slots $2,\ldots,m$, choose $\ell$ public free bit positions and fill the other $b-\ell$ positions with independent fair bits. The free positions are deterministic and may be identical across slots. The random fixed bits are independent between slots and independent of slot 1. Scatter an $\ell$-bit counter into the free positions, discard abscissas at least $p$, and retain every affine point above the others. This takes $w\,\mathrm{poly}(b)$ expected bit operations per slot with the inherited verified square-root routine. Each slot contains at most $2w$ points. Every affine point has marginal acceptance probability $\rho$.

Fix a nonidentity target $T$. Write $H$ for the ordered signed tuple count, $X$ for its generic part with the original addition guards, and $F_j$ for the sampled point set in slot $j$. The original deterministic guard count and independent one-point probabilities give $\mathbb E X\ge(1-3m/N)\lambda$ for $N\ge30m$. The pair calculation needs a different argument.

Let the free positions in a fixed-bit slot form at most $s$ runs. Write $S$ for all words supported on these positions. A run of $a$ bits has at most $2^{a+1}-1$ differences, up to its place value. Multiplying over runs gives $|S-S|\le 2^s w$ for integer differences, hence also modulo $p$. For $D\ne O$, the rational function $x(P+D)-x(P)$ has two distinct poles, each of order two, at $-D$ and $O$. Its degree is four. Each field value has at most four affine preimages. This holds for nonzero two-torsion as well. [Pole order of the Weierstrass abscissa, Sutherland's Lecture 24](https://math.mit.edu/classes/18.783/2019/LectureNotes24.pdf)

Set $A_j(D)=\sum_P \Pr(P,P+D\in F_j)$, with identities contributing zero. Two retained abscissas differ by an element of $S-S$, since they have the same fixed bits. Each admissible pair is retained with probability at most $\rho$. Therefore for $j\ge2$,

$$
A_j(D)\le4\rho\,2^s w\ (D\ne O),\qquad
A_j(O)=(N-1)\rho,\qquad
\sum_D A_j(D)=\mathbb E|F_j|^2\le(2w)^2 .
$$

To bound $\mathbb EH^2$, choose both tuples in slots $2,\ldots,m$. Their completions in slot 1 are determined. If their abscissas differ, their probability of joint acceptance in slot 1 is at most $\rho^2$ by the affine pair bound. This contributes at most $\rho^2(2w)^{2m-2}$. Identity completions contribute zero.

Otherwise the two completions are equal or opposite. In the equal case the remaining differences sum to $O$. In the opposite case negate the first remaining tuple, using the sign symmetry of every slot. The resulting differences sum to $2T$. For either sum, the all-zero differences contribute at most $\rho[(N-1)\rho]^{m-1}\le\lambda$. In any other case choose an index $j$ of a nonzero difference. After all other differences have been chosen, its value is forced. Bound its mass pointwise by $4\rho 2^s w$, sum all other masses, then multiply by the at-most-$\rho$ completion probability. Summing over $j$ and both signs gives

$$
\mathbb EH^2\le \rho^2(2w)^{2m-2}
       \left(1+\frac{2(m-1)2^s}{w}\right)+2\lambda.
$$

Since $X\le H$, Cauchy–Schwarz proves the uniform lower bound

$$
\Pr(X>0)\ge
\frac{(1-3m/N)^2}
{(2M/N)^{2m-2}\bigl(1+2(m-1)2^s/w\bigr)+2/\lambda}.
\tag{1}
$$

## Choosing the free bits

The main Fourier transform has length $L_F=L_1V=\Theta(b/(\log b)^2)$, with $L_1\le V\le2L_1$, limb width $k=\Theta((\log b)^2)$, and length-$W$ boundary words with $L_FW=(4+o(1))b$. The source row owning bit $i$ is $\lfloor i/k\rfloor\bmod V$. Sort the $b$ bit positions first by owner and then by bit position, and take the first $\ell$. Every owner holds at least $k\lfloor b/(kV)\rfloor$ bits. Thus the number $h$ of rows containing free bits is $O(1+\ell V/b)$. Inside each block of $kV$ consecutive positions the free bits form an initial segment, possibly empty. Consequently $s\le\lceil b/(kV)\rceil$.

Choose the smallest dyadic $w$ with $N\rho\ge\max(2N^{1/m},20m)$, as in the supplied manuscript. At $m\asymp(\log b)^{1/4}$, we have $\ell=b/m+O(1)$, $\lambda\ge2^m$, $h/V=O(1/m+1/V)$, and $s=O(1+\sqrt b/\log b)=o(\ell)$. Hence $m2^s/w=o(1)$. Hasse's bound gives $M/N=O(1)$. Equation (1) therefore gives inverse coverage at most $2^{O(m)}$, uniformly in $T$. Its logarithm is lower order than $Z_b$.

## Specializing the graph

For a fixed-bit slot, rows outside the $h$ active owners know every abscissa limb, including padding. Its first-stage Fourier symbols are public constants at both endpoints. Delete those edges and substitute their canonical symbols. Each active source row already recovers its actual limbs by local inversion. It checks the prescribed bits there. Retain the range, distinctness, and carry tests needed by the arithmetic. Remove the affine hash output and its hash-only carries for each changed slot. A literal implementation can specialize these hash expressions to zero, then pin their unique outputs and their canonical signed-zero carries. Given the surviving symbols, expansion of the public values is unique. This preserves the count for the new acceptance predicate, including after a consistent pin. No factor table is read by a graph vertex.

Let $t_m=\lfloor(m-1)/2\rfloor$. Before specialization, the supplied graph contains exactly $8m-7+3t_m$ complete transform arrays. After removing $m-1$ hash arrays and the inactive source rows of $m-1$ abscissa arrays, the number of full-array equivalents is

$$
6m-5+3t_m+(m-1)h/V=(7.5+o(1))m.
$$

This uses the old full affine hash for slot 1. Its cost is included. Multiplication by $L_FW$ gives $(30+o(1))mb$ Boolean word edges. The previous short-hash construction used $(36+o(1))mb$. Auxiliary binary edges remain $O(mb/\log b)$. Alphabet rounding, graph size and polynomial-time local entry bounds meet the same contraction hypotheses. Deleting edges does not increase maximum degree or parallel multiplicity. Local public precomputation depends only on the sampled patterns and costs polynomial time and space.

Three classes now each cost $(10+o(1))mb$ Boolean word edges. They are the remaining coordinates and slopes, the shifted quotient inputs, and the arithmetic outputs. The membership hash and the surviving boundary words for the sampled abscissas together cost $o(mb)$. Further savings confined to those latter words cannot improve the leading coefficient of this construction.

## Charging the complete reduction

With $K=2\pi\sqrt{14\log_2 90}/(9\ln2)$, the supplied conditional contractor and its pinning reduction give a verified witness attempt of cost

$$
\log_2 C_m\le (15K+o(1))mb/\sqrt{\log_2 b}.
$$

The deterministic column bound is $B=2mw$. Enumeration costs $B\,\mathrm{poly}(b,m)$, with $\log_2 B=b/m+O(\log b+\log m)$. As in the original reduction, retain slot and point labels through cofactor projection. Training and freshly masked queries have the same target law and use the same witness-selection rule. For a fixed seed the successful row law therefore agrees, and (1) bounds the reciprocal of the seed-averaged success $\kappa$. The supplied one-pilot theorem with twice the row-space dimension in training charges $O((B+BC_m)/\kappa)$ up to polynomial factors, the sparse solve costs $O(mB^2)$ field operations, and each subsequent query costs $O(C_m/\kappa)$ in expectation up to polynomial factors. The original fallback for small subgroup order is also covered. Retries reuse space. The inherited sequential contraction, factor table and sparse solver give the previous space balance with the new $C_m$.

For one target the leading log-time terms are $b/m+15Kmb/\sqrt{\log_2 b}$ and $2b/m$. Rounding $m=(15K)^{-1/2}(\log_2 b)^{1/4}(1+o(1))$ balances them. This gives the coefficients $2\sqrt{15K}$ in time and $\sqrt{15K}$ in space. All bounds are asymptotic. The calculation supplies neither a measured speedup nor independent validation of Assumption A.1.

Relative to the supplied intermediate manuscript at 27.0138774237, the new time coefficient is 11.1477% lower. The original public Akshima et al. manuscript states a conditional big-O bound and gives no numerical coefficient for this comparison. [Baseline provenance and comparison class](../../README.md)

## Reproduction and scope

The [reproduction instructions](../../README.md#reproduce-the-experiments) run the compiler check with its pinned C38 reference modules. It exhausts 100,352 weighted seed combinations across three tiny configurations, compares the targetwise first and second moments against the formulas, and checks guarded coverage at every nonidentity target. The order-20 example includes two-torsion. The geometric check also runs on a curve of order 103. Two actual C38 predicate fixtures compare pins, accept valid witnesses, and reject corruptions. Their edge counts are 4,688 to 3,988 and 8,122 to 7,202. The asymptotic coefficients come from the proof. [Source specialization](../../checks/strided_compiler.py) · [Verification](../../checks/check_strided_compiler.py)

A separate analysis gives a sufficient weak statistical condition on the original all-affine sampling law. The coefficient above is independent of that condition. [Weak planted detection and witness recovery](weak_detection.md)
