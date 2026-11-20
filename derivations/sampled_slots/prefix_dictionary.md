# A fixed-target cost for exact partial-sum dictionaries

On a typical successful sampled instance, the middle stages of a left-to-right addition chain cannot all replace their full group elements by short addresses in prefix-dependent dictionaries. The result allows the dictionary maker unlimited computation, knowledge of the prescribed target, and a limited message from outside the prefix. It charges an address that identifies **one exact partial sum**. The proof averages over fresh seeds at one fixed target. The companion block-support argument examines the full support of a realized block. [Block support](block_support.md)

## Collision energy of a sampled block

Let $G=E(\mathbb F_p)$, $N=|G|$, $M=2^b$, and let slot $j$ have dyadic width $w_j$, marginal affine-point inclusion probability $\rho_j=w_j/M$, and $\nu_j=N\rho_j\ge1$. Slots are independent. Each is an odd affine abscissa hash or a fixed-bit mask with at most $s$ runs of free bits. Write $F_j$ for its signed point set, so $|F_j|\le2w_j$. For a block $I$ of $k$ slots let $\nu_I=\prod_{j\in I}\nu_j$, and let $H_I(U)$ count its ordered tuples of sum $U$. For an absolute constant $C$, for all sufficiently large $p$,

$$
\mathcal E_I:=\mathbb E\sum_{U\in G} H_I(U)^2
 \ \le\ \nu_I + C^k k\,2^s\,\nu_I^2/N.
\tag{1}
$$

Here an identity block sum is included. To see (1), define the mean difference multiplicity $A_j(D)=\sum_P\Pr(P,P+D\in F_j)$. We have $A_j(0)\le\nu_j$ and $\sum_D A_j(D)=\mathbb E|F_j|^2\le4w_j^2$. In an affine hash, distinct abscissas have pair probability at most $\rho_j^2$. For $D\ne O$, equal abscissas require $[2]P=-D$, with at most four solutions. Hence $A_j(D)\le N\rho_j^2+4\rho_j\le5\nu_j^2/N$.

In a fixed-bit slot, two accepted canonical abscissas differ by an element of the difference set of the free-bit words. If the run lengths are $a_1,\dots,a_s$, this difference set has size at most $\prod_i(2^{a_i+1}-1)\le 2^s w_j$ (use the actual number of runs if smaller). For any nonzero group element $D$, the function $P\mapsto x(P+D)-x(P)$ has distinct double poles at $-D$ and $O$, hence at most four preimages per field value. Each admissible pair is retained with probability at most $\rho_j$. This gives $A_j(D)\le4\rho_j2^s w_j=O(2^s\nu_j^2/N)$, using $M/N=O(1)$ by Hasse.

Independence gives $\mathcal E_I=\sum_{\sum_{j\in I}D_j=0}\prod_{j\in I}A_j(D_j)$. The all-zero difference contributes at most $\nu_I$. For every other choice, mark one nonzero $D_j$. It is forced after all other differences are specified. Bound its factor pointwise and sum all remaining factors by their total masses. Summing this upper bound over the possible marked positions proves (1). This avoids any second-moment estimate at a specified target or any exception for two-torsion targets.

## A fixed target and an adaptive dictionary

Write $\lambda=(\prod_{j=1}^m\nu_j)/N$, and suppose $1\le k<m$. The target $T$ may be any fixed element of $G$. A dictionary $\mathcal D_I$ is a set of at most $L$ group elements, constructed from $T$, the block seeds and independent randomness. Let $E_I$ be the event that some complete tuple sums to $T$ with its block sum in the dictionary. Additional guards only reduce this event. Then

$$
\Pr(E_I)\ \le\ \frac{\lambda}{\nu_I}\sqrt{L\mathcal E_I}.
\tag{2}
$$

Condition on everything used to make the dictionary. Every block tuple has at most $N^{m-k-1}$ possible completing tuples in the group, even with the last complementary point forced. Any particular affine completion is present with probability $\prod_{j\notin I}\rho_j$. The expected number of completions of all listed block tuples is bounded by $(\lambda/\nu_I)\sum_{U\in\mathcal D_I}H_I(U)$. Cauchy–Schwarz over the dictionary and then over the block seeds gives (2). The argument permits unbounded preprocessing of the block seeds. If the maker also receives at most $a$ fixed bits from outside the block, take the union of its lists over all $2^a$ messages before conditioning. This replaces $L$ by $2^aL$ without an independence assumption on the actual message. For variable-length messages of maximum length $a$, use $2^{a+1}$.

Suppose $m\to\infty$, $m+s=o(b)$, and $\sum_j|\log_2\nu_j-(\log_2 N)/m|=o(b)$. These conditions include both the balanced paid widths and widths reduced by a constant factor in each slot. They imply $\lambda=2^{o(b)}$. For $k/m$ between $1/4$ and $1/2$, $\nu_I=N^{k/m}2^{o(b)}$, so (1) gives $\mathcal E_I\le(1+o(1))\nu_I$ uniformly. Any procedure retaining fixed-target success at least $2^{-o(b)}$ through such a listed state must therefore satisfy

$$
\log_2 L+a\ \ge\ k(\log_2 N)/m-o(b).
\tag{3}
$$

In particular, take all middle prefixes of a sequential chain. Even if list sizes vary from seed to seed, restrict each list to its first $2^{\lfloor b/8\rfloor}$ entries in any deterministic order. On an instance with a smaller original list, this keeps it entire. With $o(b)$ outside advice bits per stage, (2) and a union bound give probability at most $2^{-b/16+o(b)}$ that a successful instance has even one middle prefix with a shorter list. Conditional on success at least $2^{-o(b)}$, all those prefixes therefore have more than $2^{\lfloor b/8\rfloor}$ entries except with exponentially small probability. A representation transmitting a fixed-width exact-state address at every one of those cuts has total width $\Omega(mb)$ on these instances. If the address itself has a variable-length code, the same truncation applies to the set of states having short codes, as long as the corresponding codebook is block-local as specified above.

The conclusion concerns a chain with separately decodable partial sums. A balanced tree has fewer macroscopic cuts. A verifier can also distribute the computation without exposing any full partial sum, or generate a joint description that already incorporates substantial information about the complementary slots. Equation (2) does not bound these models or the running time of arbitrary joint verification. [Finite energy and dictionary checks](../../checks/check_grid_values.py)
