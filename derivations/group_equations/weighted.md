# A half-length arithmetic transform when the prime is close to a power of two

Suppose $p=2^b-c$ is prime and $1\le c\le b^C$ for a fixed $C$. We obtain a $(10+o(1))mb$ Boolean word-edge count on the fixed-bit slot construction. Under the project's inherited assumptions this gives time coefficient $2\sqrt{5K}=13.85782672$ and space coefficient $\sqrt{5K}=6.92891336$ in $Z_b=b/(\log_2b)^{1/4}$. No divisibility of $b$ by the transform length is required. Weighted transforms for nonaligned digits are classical. The work here is constructing the exact root system, verifying the curve equations with local carries, and showing how many words cross the actual graph boundaries. [Weighted transform precedent](https://doi.org/10.1016/j.jco.2016.03.001) · [Split-radix finite Galois transforms](https://www.reed.edu/physics/faculty/crandall/papers/confgt.pdf)

## Mixed digits and folded products

Choose a power of two $n=\Theta(b/(\log b)^2)$. Set $u_i=\lceil ib/n\rceil$ for $0\le i\le n$, $k=\lceil b/n\rceil$, and $B_i=2^{u_{i+1}-u_i}$. Each source is a canonical integer written $\sum_{i<n}a_i2^{u_i}$, where $0\le a_i<B_i$. The limb widths differ by at most one bit. For $r=i+j-tn$ and $t=\lfloor(i+j)/n\rfloor$, the product of digits at $i,j$ contributes to position $r$ with the integer factor

$$
              2^{u_i+u_j-u_r-tb}c^t .
$$

The exponent of 2 is zero or one by the ceiling identity. The folding preserves evaluation modulo $p$, since $2^b\equiv c\bmod p$. Any C38 quadratic residual thus has a folded signed coefficient vector $g$ with $\sum_i g_i2^{u_i}\equiv D\bmod p$. Each $|g_i|$ is $O(n(1+c)2^{2k})$, using the constant number of summands and bounded integer multipliers per residual. Public field coefficients are supplied in their own digit representation.

## An exact weighted Fourier transform

Choose a prime $q\equiv-1\bmod 8nc$. The effective least-prime bound used by the supplied manuscript gives $q\le (8nc)^{O(1)}=\mathrm{poly}(b)$, so deterministic search is polynomial in $b$. Take the least $q^h$ greater than twice a safe bound for the folded coefficients. Then each base-ring word in $A=\mathbb Z/q^h\mathbb Z$ uses $W=2k+O(\log b)$ bits.

Both 2 and $c$ are squares modulo $q$. For 2 use $q\equiv7\bmod8$. For each odd prime $l\mid c$, reciprocity gives $(l/q)=1$ since $q\equiv-1\bmod l$ and $q\equiv3\bmod4$. These square units lift modulo $q^h$. Since the square subgroup has odd order $q^{h-1}(q-1)/2$, raising to the power $n$ permutes it. Compute $\beta\in A^\times$ such that $\beta^n=c\,2^{-b}$, by inverting $n$ in that odd group order. No factoring of $c$ is needed by the construction.

Work in $A[I]/(I^2+1)$ and use conjugation $I\mapsto-I$. The residue field is $\mathbb F_{q^2}$. Its norm-one subgroup has order $q+1$, divisible by $n$. Enumerate elements in the polynomial-size residue field to find a root of exact order $n$, then lift the simple root of $X^n-1$ to the extension modulo $q^h$. Its lift $\omega$ satisfies $\bar\omega=\omega^{-1}$ by uniqueness of lifting. Its distinct powers have unit differences. This suffices for exact invertible transforms and polynomial local construction.

Define real, meaning conjugation-fixed, weights $t_i=2^{u_i}\beta^i\in A^\times$. Apply the cyclic transform to $t_i a_i$ for each source, evaluate the original quadratic expression at each frequency, and invert. Untwist coefficient $r$ by $t_r^{-1}$. Indeed $t_i t_j/t_r=2^{u_i+u_j-u_r-tb}c^t$, giving exactly the folded coefficients above after centering. Linear terms and digit expansions of public constants are transformed with the same weights. There is no public reciprocal convolution and no shifted quotient.

## Counting what is transmitted

Factor $n=L_1 V$ with $L_1\le V\le2L_1$. Source row $j$ owns indices $i=aV+j$. Send the row transform before its interstage twiddle,

$$
                 U_{j,s}=\sum_a t_{aV+j}a_{aV+j}(\omega^V)^{as}.
$$

It satisfies $U_{j,-s}=\overline{U_{j,s}}$. At $s=0,L_1/2$ the values lie in the base ring. At each other pair store a single extension value, which costs two base words. Thus a row needs exactly $2+2(L_1-2)/2=L_1$ base words. One middle processor handles each pair of frequencies, supplies the known factors $\omega^{js}$, and performs the second-stage transforms, pointwise arithmetic and second-stage inverses. After reversing the twiddle it transmits the same compressed row format for each output. Overall frequency negation sends $s+L_1 t$ to the other member's frequency block (or to the same block at a self-paired $s$). All coefficient operations and weights are real, so the omitted frequencies are determined on both boundaries. Source predicates invert their own row, untwist its digits and check range and fixed bits. Closure predicates do the same inverse row operation to recover their signed coefficients. Each remaining full array therefore uses **$nW=(2+o(1))b$** Boolean word bits. Charging a single base word for an extension value would be wrong. Keeping both members of a conjugate pair would double-charge them.

Closure imposes $g_i+h_i=B_i h_{i+1}$ with $h_0=c h_n$. Multiplying by $2^{u_i}$ and summing gives $\sum_i g_i2^{u_i}=p h_n$. Conversely divisibility fixes $h_n$, hence $h_0$, and recursively forces every intermediate integral carry. The absolute value of $h_n$ is at most $O(\max_i |g_i|/(2^{k-1}-1))$ for large $b$, by summing the geometric bound on the place weights and using $p\sim2^b$. The extra factor $c$ at the wrap and the intermediate recurrence still fit $O(k+\log b)$ bits. Consecutive owner merging as in C38 makes the carry budget lower order. Routing one wrap between the last and first owner costs a single additional signed word per equation.

Choose the free abscissa bits by the new mixed-limb owners, in the same order as in the fixed-bit proof. Each owner has $L_1$ limbs, all of width $k$ or $k-1$. The active row fraction is $O(1/m+1/V)$, and the free bits form at most $L_1$ runs in normal bit order. This satisfies the same difference-set bound in the target-coverage proof. The ranges on all arithmetic variables still use the original prime $p$. The first slot's affine hash is modulo $2^b$. Keep that one ordinary hash circuit, and bind its copy of $x_1$ to the mixed-radix source by canonical binary digits. Per-bit equality relays cost $O(b)$ edges and avoid large parallel multiplicities. Each owner sees $O(b/V)$ such bits up to constant factors. Even after the binary alphabet charge $O(\log\log b)$, this first-slot cost is $o(mb)$ at the optimized arity. [Fixed-bit coverage and outer costs](../strided_slots/results.md)

The graph respects the inherited degree and multiplicity bounds: pairing middle processors costs at most a constant factor, products in the quadratic extension are local, and all words between a given pair of arithmetic processors are counted in the ledger. The $(2.5+o(1))m$ full-array equivalents of coordinates/slopes and $(2.5+o(1))m$ arithmetic output arrays therefore cost $(10+o(1))mb$ in total. No additional full group-coordinate or quotient array is needed. The contraction coefficient is $5K$. Balancing against enumeration and the sparse solve gives the stated time and space coefficients. Inverse target coverage, witness recovery by pinning, and retries obey the previous accounting.

The same arithmetic argument applies if an odd positive $v=\mathrm{poly}(b)$ and an $s=b+O(\log b)$ satisfy $vp=2^s-c$ with positive $c=\mathrm{poly}(b)$: check the residual $vD$ against the latter modulus, since $vp\mid vD$ exactly when $p\mid D$. The upper $s-b$ bits of each canonical field source are public zeros, and the factor $v$ adds only $O(\log b)$ to the coefficient width. This observation provides no such sparse multiple for arbitrary primes.

The finite regression exercises a genuine quadratic extension over a prime power and compares the compressed transforms, signed integer folding, and closure against exact finite-field arithmetic. It checks all 320 pairs of basis vectors in two nonaligned configurations. It then evaluates the actual C38 quadratic list on two fields at arities three, four and five. All 44 equations of the valid witnesses passed and all 18 deliberately perturbed witnesses failed the joint check. These are component tests. The coefficient includes the full outer reduction by proof. [Reproduction script](../../checks/check_weighted.py)
