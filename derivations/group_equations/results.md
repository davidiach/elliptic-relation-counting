# Removing the quotient input from the curve equations

The curve equations can be verified with one output transform per equation and zero input transforms for the integer quotients. For arbitrary primes, our exact construction uses wider words and ties the project's previous conditional time coefficient of 24.00246. For primes with bounded signed radix support the coefficient falls to 19.59793. For $p=2^b-c$ with positive $c=\mathrm{poly}(b)$, a weighted cyclic transform lowers it further to 13.85783. All three bounds have the same scale $Z_b=b/(\log_2 b)^{1/4}$, where $b$ is the field bit length. We obtained no $o(Z_b)$ bound.

| Construction using the fixed-bit slots | Leading Boolean word-edge count | Coefficient of $Z_b$ in log expected time | Scope |
|---|---:|---:|---|
| Previous shifted-quotient graph | $30mb$ | 24.00246 | All primes in the inherited input class |
| Public reciprocal and one output | $30mb$ | 24.00246 | Same class |
| Local radix division and one output | $20mb$ | 19.59793 | Bounded signed radix support |
| Weighted cyclic arithmetic and one output | $10mb$ | 13.85783 | $p=2^b-c,\ 0<c=\mathrm{poly}(b)$ |

The bit bounds use finite-packet Assumption A.1 and the inherited contraction and verified-solver estimates. They include the factor bases, training, sparse solve, and retries. The savings below do not require a new hypothesis about coverage. [Previous bound and sampling proof](../strided_slots/results.md) · [Source for the supplied C38 compilation and contraction](../../README.md#attribution-and-provenance)

## An exact general-prime predicate

Let $D$ be one of C38's integer quadratic residuals on canonical field representatives. We know $|D|<8p^2$. Set $E=D+8p^2$, so $0<E<16p^2$. With the inherited radix $B=2^k$ and power-of-two transform length $L$, use $B^L-1>16p^2$. The existing padding suffices. Choose $\epsilon=+1$ unless $p\mid B^L-1$, in which case take $\epsilon=-1$. Since $p$ is odd, $R=B^L-\epsilon$ is coprime to $p$. Compute the public constant $U=p^{-1}\bmod R$.

Let $z$ be the least nonnegative residue of $UE$ modulo $R$. Then

$$
                             p\mid D \quad\Longleftrightarrow\quad z<16p .
$$

Indeed, a divisible $E$ has $E/p$ in that interval. Conversely, $pz\equiv E\bmod R$, and if $z<16p$, both integers are in $[0,16p^2)\subset[0,R)$. They must be equal. The bound on the lifted residual is essential to this test.

Form the radix polynomial $E(Z)$ from the same quadratic operations on the source limbs and let $U(Z)$ contain the public digits of $U$. Compute

$$
                 C(Z)=U(Z)E(Z)\pmod{Z^L-\epsilon}.
$$

Evaluation at $B$ preserves this congruence modulo $R$. At each Fourier point, the middle processor evaluates the quadratic expression, multiplies by the precomputed public value of $U(Z)$, and returns one output. For the negative sign, twist each input coefficient at degree $i$ by $\tau^i$ for a primitive $2L$-th root $\tau$. Undo that twist in the source predicates before checking digits and fixed bits, and at closure before checking coefficients. A prime $1\bmod 2L$ and a power of that prime provide the same root and invertibility conditions as in the original compiler. Twisting uses local public constants, so the previous fixed rows of the sampled abscissas are still fixed. The separate ordinary modulo-$2^b$ hash on $x_1$ may use its original graph and be linked by $O(b)$ binary relays, a lower-order cost. [Classical cyclic and negacyclic convolution](https://members.loria.fr/PZimmermann/papers/fft.pdf)

Write $c_i$ for the centered, exactly recovered coefficients of $C$. At closure compute the radix digit $d_i$ and carry $h_{i+1}$ by Euclidean division:

$$
             c_i+h_i=d_i+Bh_{i+1},\qquad
             0\le d_i<B,\qquad h_0=\epsilon h_L .
$$

The local comparison tests $\sum_i d_iB^i<16p$, using initialized radix borrows. The digits are computed where the coefficients are recovered. No full copy of $z$ is shipped. Summing the recurrence gives $C(B)-\sum_i d_iB^i=Rh_L$. An accepted digit vector therefore represents the unique $z$ in the stated interval. Given that $z$, $h_L$, $h_0$, and all intervening carries are unique. This excludes the extra all-ones representation in the cyclic case. In the negacyclic case the sole residue $B^L$ without an $L$-digit representative is outside the acceptance range. Existing guards on abscissas, slopes, and canonical representatives still apply. An accepted assignment of the remaining coordinates has exactly one extension, so pinning to recover a witness uses the same reduction.

### Where the general-prime saving goes

Each raw quadratic coefficient is $O(LB^2)$. Multiplication by the $L$ public radix digits of $U$ bounds a folded coefficient by $O(L^2 B^3)$. Exact centered reconstruction under this bound uses a Fourier modulus of width $3k+O(\log b)$, compared with $2k+O(\log b)$ before. The wrapped carry has $2k+O(\log b)$ bits. Routing $h_L$ to initialize $h_0$ costs one such word per equation. The number of carry boundaries and the auxiliary alphabet charges remain lower order under the inherited merging schedule.

Since $Lk=(2+o(1))b$, each remaining full transform costs $(6+o(1))b$ word bits. There are $(2.5+o(1))m$ full-array equivalents for the sources after fixed-row specialization and $(2.5+o(1))m$ arithmetic outputs. The total is $(30+o(1))mb$. Charging only the deleted quotient arrays would incorrectly predict a general-prime improvement. A tighter construction or coefficient-recovery bound could still change this calculation.

## When the prime has few nonzero signed radix digits

If $p=\sum_{d\in S}p_dB^d$ with a bounded number of public terms and odd $p_0$, the closure can infer the quotient one radix digit at a time. Each inferred digit depends on one carry and a bounded set of earlier quotient digits. Grouping coefficient owners so that few of those dependencies cross between groups makes the communication per equation $o(b)$, even after charging relays and alphabet recoding. The ordinary $2k+O(\log b)$-bit plane now suffices for both sources and raw arithmetic outputs. This gives $(20+o(1))mb$ leading word edges. The sparse assumption concerns the radix used by the compiler. For example, it covers $p=2^b-c$ for positive $c=\mathrm{poly}(b)$ at the stated limb width. It gives no guarantee for a general dense prime. [Sparse recurrence, routing proof and checks](sparse.md)

For $K=2\pi\sqrt{14\log_2 90}/(9\ln 2)$, the general and sparse budgets give contraction coefficients $15K$ and $10K$, respectively, in $mb/\sqrt{\log_2 b}$. The outer time terms are $b/m+A mb/\sqrt{\log_2 b}$ and $2b/m$. Optimizing $m$ gives $2\sqrt A Z_b$. Thus the sparse time coefficient is $2\sqrt{10K}=19.59792649$, an 18.35% reduction relative to 24.00245996 on qualifying inputs. Its corresponding space coefficient is $\sqrt{10K}=9.79896325$.

## When multiplication can fold directly modulo the prime

For $p=2^b-c,\ 0<c=\mathrm{poly}(b)$, mixed-width binary digits and a weighted cyclic transform reduce a quadratic product modulo the curve prime before it leaves the Fourier plane. This cuts the number of words per complete array in half. Conjugate frequency pairs transmit one two-component value, and self-paired frequencies transmit one base word. After charging these components, an array costs $(2+o(1))b$ bits. A cyclic carry verifies exact divisibility without the quotient array. The time coefficient becomes $2\sqrt{5K}\approx13.85782672$, with space coefficient 6.92891336. This cuts the time coefficient by 42.26% relative to the previous 24.00245996 on qualifying inputs. The benefit depends on the size of $c$: for a general $b$-bit prime, the folded coefficient bound can require an extra $\Theta(b)$ bits per word. This construction includes arbitrary alignment of the mixed digits and the original first-slot hash. [Weighted transform and routing proof](weighted.md)

## What limits the exponent

The previous graph already carries $(10+o(1))mb$ word bits in coordinates and slopes. The weighted construction halves that cost on its subclass. Both still scale as $mb$. To make $2b/m=o(Z_b)$ would require $m/(\log_2b)^{1/4}\to\infty$. At that arity the inherited $mb/\sqrt{\log_2 b}$ contraction estimate already exceeds $Z_b$ in order, even if checking all equations were free. This is an accounting limit of the retained representation and analysis. It establishes no time lower bound for other representations or contractors.

Two other approaches give more specific information.

* A random linear combination can give zero on many *different* invalid assignments. On two negative finite fixtures, all 64 sampled batches accepted at two of five equations and four of seven equations, respectively. A backward-fiber argument explains the effect for growing fields. It works directly with curve residuals. [Batching proof and finite experiment](batching.md)
* Eliminating intermediate variables into an exact algebraic relation incurs degree. A generic fixed-target relation has exponentially many possible final abscissas as the signs of the earlier factors vary. This gives a Bézout constraint on lifts with bounded-degree equations. Its scope is algebraic correctness on a dense open over the algebraic closure. Proofs for a restricted finite input distribution need a separate analysis. [Degree of the chord relation](chord_geometry.md)

## Verification

The four checks against the original C38 implementation use the local reference modules described in the [reproduction instructions](../../README.md#reproduce-the-experiments). Run `python reproduce.py --source-archives /path/to/local/archives` from the repository root to exercise these checks with their pinned inputs. The individual scripts are listed below. The reciprocal test exhausts 36,160 bounded scalar inputs and both wrap signs, including an odd prime dividing $B^L-1$. It also runs both Fourier stages of the modified arithmetic and checks every recovered coefficient against direct integer polynomial convolution. Among 660 checks of C38 equations on three curves, 48 residuals vanished and 612 did not. The sparse and batching tests cover separate claims. These checks do not run the conditional packet contractor end to end.

~~~text
python -B checks/check_reciprocal.py
python -B checks/check_sparse.py
python -B checks/check_weighted.py
python -B checks/check_batch.py
python -B checks/check_geometry.py
~~~
