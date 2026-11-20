# Sparse mixed slots isolate one signed address

Fix a smooth short Weierstrass curve over $\mathbb F_p$, $p>3$, and a target $T$ with $2T\ne O$. Write $N=\#E(\mathbb F_p)$, $b=\lceil\log_2p\rceil$, $M=2^b$, $w=2^\ell$, $\rho=w/M$, and $\lambda=N^{m-1}\rho^m$. Let $H$ count all ordered signed affine tuples with sum $T$, and $X$ count those that pass the usual intermediate-sum guards. Sampling uses an odd affine abscissa hash in slot 1. In each other slot, fix the locations of $\ell$ free bits and sample the remaining bits uniformly. Different slots are independent. Retain both affine signs above each legal abscissa. Let $s$ bound the number of runs of free bits in each fixed-bit slot.

It is possible to choose the common width so that $1/\lambda\le 2^{O(m)}m^2$ and almost every positive seed has exactly one guarded tuple. For such a seed and target the whole relation predicate, including its restrictions to pinned addresses, admits an $O(b+m)$-bit description. The unknown address is advice specific to both the seed and the target. Computing it from the sampled masks remains open at the required cost.

## A factorial moment bound for this sampling law

For $N\ge30m$, the deterministic guard count on the full group gives $\mathbb EX\ge(1-3m/N)\lambda$, since each affine tuple survives with probability $\rho^m$. Also $\mathbb EH\le\lambda$ by forcing the last point. We claim

$$
\mathbb E[H(H-1)]\le
\lambda^2(2M/N)^{2m-2}\left(1+2(m-1)2^s/w\right). \tag{1}
$$

Choose the two candidate tuples in slots $2,\ldots,m$. Both slot-1 completions are forced. If their abscissas differ, the affine hash accepts them jointly with probability at most $\rho^2$. Indeed their hashed difference is uniform on residues of the same dyadic valuation. Conditional on acceptance of the first abscissa, the second has acceptance probability $\rho$ or zero according to whether this valuation is below $\ell$. Each remaining slot has at most $2w$ points. Thus the unequal-abscissa contribution is at most $\rho^2(2w)^{2m-2}$.

Equal abscissas mean equal points or opposite points. In the equal case the differences of the remaining two tuples sum to $O$. All differences zero gives the exact diagonal, which the factorial moment excludes. In the opposite case negate the first tuple in the remaining slots. Their sets are sign symmetric and the differences now sum to $2T$. All differences zero is impossible. Double-counting a two-torsion completion in these upper bounds is harmless.

For a remaining slot with random set $F$, put $A(D)=\sum_P\Pr(P,P+D\in F)$. We have $\sum_D A(D)=\mathbb E|F|^2\le(2w)^2$. The difference set of the free-bit words has cardinality at most $2^s w$: a run of $a$ bits has at most $2^{a+1}-1$ differences, scaled by its place value. For $D\ne O$, the rational function $x(P+D)-x(P)$ has distinct double poles at $-D$ and $O$, hence degree four. Each field value has at most four affine preimages. Two accepted abscissas must have difference in the free-bit difference set. Their joint retention has probability at most $\rho$. Thus $A(D)\le4\rho\,2^s w$ for $D\ne O$. [Pole order of the abscissa, Sutherland, Example 24.8](https://math.mit.edu/classes/18.783/2019/LectureNotes24.pdf)

Select a nonzero difference among the $m-1$ remaining slots. It is forced by the other differences and their prescribed sum. Bound its mass pointwise, sum all the other masses, and multiply by the slot-1 completion probability at most $\rho$. The two signs contribute at most $2(m-1)\rho(4\rho\,2^s w)(2w)^{2m-4}$. Add the unequal-abscissa contribution and divide by $\lambda^2$ to obtain (1).

## Width, success and retries

Let $D_0=2M/N$ and use the largest dyadic $w$ for which $\lambda\le\lambda_*=(4m^2 D_0^{2m-2})^{-1}$. In the growing-arity regime of interest this gives

$$
  2^{-m}\lambda_*<\lambda\le\lambda_*,\qquad
  \ell=b/m+O(1),\qquad
  \lambda^{-1}\le 2^{O(m)}m^2.
$$

Here Hasse's bound gives $N\asymp M$. Suppose $m2^s/w=o(1)$, as with contiguous free bits, or the strided layout with $s=O(1+\sqrt b/\log b)$ and $m=O(\sqrt b)$. From (1), eventually $F_2=\mathbb E[H(H-1)]\le\lambda/(2m^2)$. Since $H1_{\{H\ge2\}}\le H(H-1)$,

$$
\Pr(H=X=1)\ \ge\ \mathbb EX-F_2\
\ge\lambda(1-3m/N-1/(2m^2)).
$$

The chance of two or more tuples is at most $F_2/2$. Conditional on $H>0$, the chance of ambiguity is $O(m^{-2})$. A unique unguarded tuple has unconditional probability at most $\mathbb E(H-X)\le3m\lambda/N$. These estimates are targetwise and permit empty slots. They impose no assumption that the sets for different targets are independent.

A procedure allowed to skip ambiguous cases would therefore have inverse mean seed success at most $2^{O(m)}m^2$ **if** it could obtain the appropriate advice from the seed and the target. This factor charges $O(m+\log m)$ extra bits of log work in the inherited one-pilot reduction when training and querying use the same rule. The number of retained labels is still at most $2mw$. This argument bounds the loss of sampling probability. It supplies no bound on the cost of identifying a singleton.

## What is in the representation

In slot 1 a label consists of its $\ell$-bit hash address and one sign bit. In each other slot, use the actual free bits and a sign bit. Give the smaller affine ordinate sign 0 and the larger one sign 1, using only sign 0 at rational two-torsion. Concatenate the labels into a bit string of length $L=m(\ell+1)=b+O(m)$. Both the curve point and the slot position can be recovered, including when different slots contain the same point.

For a fixed seed and prescribed target the guarded address set has exactly $X$ elements. Listing them gives an exact comparison predicate of $O(1+LX)$ Boolean size and description bits, with $O(1+L\lambda)$ expected description bits. A bounded-fan-in comparison tree per address suffices, with fanout implemented by equality copies if required. On a negative guarded fiber its predicate is constantly false. With one guarded address it is equality to that one literal string. With two or more, use a list for exactness or deliberately reject the seed-target pair for the singleton-only finder. The analysis above counts $X=H=1$ toward the latter finder's success; accepting additional cases with $X=1<H$ is optional. Setting any subset of address bits simply restricts these same predicates. In the singleton case, a pin either contradicts the literal or leaves a unique compatible address.

Building this explicit list requires the signed addresses of the solutions. Giving a verifier the singleton literal would already give a witness to its caller. Giving a decision procedure only a bit test on complete addresses can still require searching the address space. A prefix-*existence* oracle could recover the bits in at most $L$ calls once existence is known. In particular, an advice string prepared for one target cannot be assumed to work for the next masked target. Thus instance-wise representation size alone cannot imply a uniform $o(mb)$ construction or an improved end-to-end exponent. [Earlier all-affine isolation argument](../mixed_coefficient/mixed.tex)

The distinction matters even for seed-only preprocessing. For any set $\mathcal T$ of targets with $2T\ne O$, linearity of expectation puts the expected number of guarded-singleton targets in a seed at least $|\mathcal T|\lambda(1-3m/N-1/(2m^2))$. A literal directory storing one entry per such target must therefore charge this many entries in expectation. A capped directory of $J$ target-address pairs can answer a fresh target uniform on $\mathcal T$ with probability at most $J/|\mathcal T|$. For a subgroup of order $2^{b-O(1)}$ and $m=O(\sqrt b)$, matching a constant fraction of the isolation probability by literal table hits requires $J\ge 2^{b-O(m+\log m)}$. The argument permits a different, implicit representation shared by all targets.

## Public batching on a sparse address domain

There is also a constructive reduction in the number of scalar comparisons. Suppose a fixed seed and target have a publicly specified residual vector $R(a,z)$ over $\mathbb F_p$, where $a$ ranges over the signed address domain and $z$ over at most $r$ freely guessed field scalars. Deterministically derived intermediates may occur in $R$. Any further guessed bits or encodings with multiple possible completions must be included in the domain count. Assume the zero vector is the desired acceptance predicate and that any separate eligibility tests are kept. The address count is at most

$$
 D=(2w)^m=N\lambda D_0^m,\qquad
 D/p\le (N/p)D_0^{\,2-m}/(4m^2)=O(m^{-2})
$$

at the widths above. For the last inequality, $N/p=1+o(1)$ and $D_0\ge 1$ for sufficiently large $p$. Choose $r+1$ independent uniform linear functionals on the residual space, independently of the fixed masks and target. For any nonzero residual the chance all checks return zero is $p^{-r-1}$. There are at most $Dp^r$ candidate pairs. A union bound gives probability at least $1-D/p$ that the projected system is correct *simultaneously on every address and auxiliary assignment*, even on a fixed positive seed. All true solutions survive every challenge. When the checks are used for pinning, the same good challenge works for every restriction. The terminal candidate still needs exact verification. For a fresh independent challenge, a guarded singleton seed consequently yields a verified witness with probability at least $1-D/p$ given an exact finder for the projected system.

For a residual computed canonically from the address alone this uses one scalar test. More generally $r=o(m)$ permits $o(m)$ tests. The theorem supplies the random coefficients by sampling them once. It does not supply such a low-wire evaluation of $R$, and each coefficient, its multiplication, and all transport must be charged in a proposed graph. In particular, the alternating chord system has $e=2m-2+\lfloor(m-1)/2\rfloor$ equations and $e-1$ free non-leaf field variables after leaf abscissas have been fixed. Here the lemma gives $e$ tests, saving none. Counting signed addresses only increases the domain used in this upper bound. The earlier random-batching obstruction for freely varying chord witnesses and this positive low-auxiliary statement therefore address different domains.

## Finite verification

The checker exhausts 8,192 mixed seeds for all 18 nonidentity targets on the order-19 curve over $\mathbb F_{13}$. It tests the factorial bound, the signed-address labels and single-bit pins of every guarded singleton. For the order-126 curve $y^2=x^3+x+1$ over $\mathbb F_{127}$, take $T=(0,1)$, $m=3$, and width one. Of 2,097,152 equiprobable abscissa patterns, 2,082,019 are empty, 14,410 give a unique guarded tuple, 359 give a unique unguarded tuple, and 364 give at least two tuples. The proved lower bound on guarded-singleton probability is $0.0021467$, and the observed exact probability is $14410/2097152$. On one of these singletons a three-component residual uses the computed curve-sum coordinates and one freely guessed field checksum. Classifying kernels by their rank gives an exact two-row failure probability of about $0.00794$ over all $127^6$ challenges, below $8/127$. The fixture computes the group sum directly, so this batching experiment establishes no circuit saving. [Reproduction script](../../checks/check_strided_isolation.py)
