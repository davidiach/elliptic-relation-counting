# How many bits must a one-sided chord fingerprint carry?

Write $G=E(\mathbb F_p)$, $N=|G|$, $b=\lceil\log_2 p\rceil$, and let $\sigma$ contain all the sampled slot seeds. Each legal signed address $a$ gives one tuple of affine points with sum $S_\sigma(a)$. Write $U$ for a regular seed and target with exactly one ordinary tuple, which passes the addition guards. On $U$, let $F$ be the number of **distinct wrong legal addresses** accepted by a test that always accepts the witness. Auxiliary certificates for the same address are counted only once. Here are precision and communication costs when $\mathbb E[F\mid U]=\epsilon<1/3$. They apply regardless of whether the mistakes for different addresses are independent.

## A reference fingerprint

Suppose a candidate's group sum is available canonically. Encode an affine sum by an infinity flag, its abscissa, and its sign, with one fixed code for infinity. Apply the same independent uniform binary $r$-row linear hash to the sum's code and the target's code, and compare the outputs. On a singleton seed with $D=\prod_i|F_i|$ legal addresses, no wrong address has sum $T$. Each wrong address collides with probability exactly $2^{-r}$, irrespective of collisions among the sums of wrong addresses. Thus the witness always survives and the expected number of wrong addresses for this seed is $(D-1)2^{-r}$. Choosing $r=\lceil\log_2(3D)\rceil$ gives a strict bound below $1/3$. For a balanced join one can instead apply the same hash to the two signed half sums $P$ and $T-Q$, with exactly the same expectation. Exact group labels can be compared directly once they are computed.

This comparison has a short description and is inexpensive on one supplied tuple. Obtaining the two group sums for candidate addresses inside the counting graph, or joining the two explicit half lists, retains its previous cost. An error budget on candidate addresses does not price that work. [Exact balanced calculation](../product_differences/chords.md)

## Public seeds and a masked target

Let $T$ be uniform on a seed-independent set of $n$ eligible targets. All parties can have all of $\sigma$ and arbitrary shared preprocessing of it. Pick any fixed nonempty block $I$ of $k$ address slots. The input to one side of a cut is the block address $u$. The input to the other is the outside address $c$ and the target. Suppose the accepted pairs are covered exactly by at most $K$ rectangles of the form $A_j\times B_j$ on these two inputs. Overlaps are allowed. For instance, a deterministic binary protocol of worst-case cost $d$ has at most $2^d$ accepting transcript rectangles. A Boolean feasibility graph separated by $d$ bits has at most $2^d$ such rectangles after existentially quantifying local certificates. Even revealing $c$ to the first side for free does not affect the following bound, since the argument fixes it throughout.

Fix the seed, any random choices, and one $c$. Let $L_c$ be the number of targets with event $U$ whose witness extends this $c$. Their $L_c$ witness block addresses are distinct, as $(u,c)$ has only one group sum. Assign each witness to one accepting rectangle covering it. If $s_j$ witnesses are assigned to rectangle $j$, its cross pairs give $s_j(s_j-1)$ wrong accepted addresses, counted over those targets. These cross pairs are distinct even if the original rectangles overlap. Writing $L=\sum_c L_c$ and $D_{\rm out}=\prod_{i\notin I}|F_i|$, Cauchy–Schwarz twice yields

$$
\sum_{T:\,U}(F(\sigma,T)+1)
\ \ge\ \frac1K\sum_c L_c^2
\ \ge\ \frac{L^2}{K D_{\rm out}}.
$$

When $D_{\rm out}=0$, define the right side to be zero. The outside addresses partition all the wrong addresses counted in this proof. Another application of Cauchy–Schwarz over the seed law gives

$$
K(1+\epsilon)\ \ge\ \frac{\mathbb E L}{\mathbb E D_{\rm out}}
       =\frac{n\,\Pr(U)}{\mathbb E D_{\rm out}}. \tag{1}
$$

Public randomness can be fixed to a common tape for all inputs before taking expectations. Probability-one completeness suffices on these finite input sets. If $d'$ bits describing the target are preloaded on the block side, partition the targets by this advice before applying the argument and include the resulting factor $2^{d'}$ in $K$. Target-dependent random challenges delivered free to that side would also carry information and require this accounting.

Under the independent thinned mixed sampler each of the $N-1$ affine points survives with probability $\rho=w/2^b$. Thus $\mathbb E D_{\rm out}=((N-1)\rho)^{m-k}$. Let $\nu=N\rho$. The targetwise isolation result gives $\Pr(U)\asymp\lambda=\nu^m/N$ in the growing-arity regime; the regularity restriction loses a vanishing fraction of singleton successes. Consequently (1) yields, uniformly in the block,

$$
\log_2 K \ \ge\ k\log_2\nu+\log_2(n/N)-O(1).
$$

In particular a mask of size $n=\Theta(N)$ and a cut separating **the whole address** from the target require $b-O(m+\log m)$ bits, even for a target decoder optimized jointly for all the seed's singleton targets. A chain whose prefix blocks of sizes $1,\ldots,m-1$ reach the target side through distinct interfaces must carry $\Omega(mb)$ bits in total: sum the cut bounds and use $\log_2\nu=b/m+O(1)$. In a balanced binary tree with the target at the root, summing the bounds over the disjoint interfaces at each level gives $\Omega(b\log m)$ total for a power-of-two arity. Each level accounts for all $m$ slots, while the constant per-edge loss sums to $O(m)$. These estimates permit arbitrary messages, interactive protocols and false candidates. No interface is required to name an exact group sum. They concern a reusable, target-independent representation on the block side; a compiler that uses the full target freely inside every block can escape (1). [Isolation probability and fresh masks](../sampled_slots/isolation.md) · [Conditional regularity](../direct_parity/scalar_join.md)

## A fixed target with private block seeds

A different bound holds when $T$ is available to both sides, but each side sees only the seeds of its own slots. Suppose there are at most $K$ accepting rectangles across the cut on the inputs consisting of each side's local seeds and its local address. Pick an accepting rectangle uniformly by its index among $K$ slots, padding the list if necessary. Each side selects a uniformly random legal local address from its side of the rectangle, or fails if none exists. No communication is needed for this trial.

On a singleton there is a rectangle containing the true address. In such a rectangle the product of the two numbers of compatible legal addresses is at most $1+F$, since every crossed pair accepts. The trial therefore outputs the witness with probability at least $1/[K(1+F)]$ on that seed. Averaging and using convexity gives an unconditional witness probability at least $\Pr(U)/[K(1+\epsilon)]$. Conditional on the rectangle and public coins, the two outputs use independent seeds. For a block of $j$ slots, the chance of representing any preassigned group sum is at most $N^{j-1}\rho^j$, by choosing $j-1$ points, forcing the last, and applying independence. Condition on the other side's output and use the better of the two block bounds to obtain, with $1\le k<m$ and $\nu\ge1$,

$$
 K(1+\epsilon)\ \ge\ \frac{N\,\Pr(U)}{\nu^{\min(k,m-k)}}. \tag{2}
$$

For one private slot this requires $\log_2K\ge b-b/m-O(m+\log m)$ bits. The sum over disjoint single-slot regions charges each crossing edge at most twice, again requiring $\Omega(mb)$ total Boolean edges. The full target can be public in (2). Shared preprocessing that reveals other slots' seeds can invalidate the independence step. Together (1) and (2) price two common ways of delivering the inputs to a test. They leave open computations with all public seeds and target inside the same arithmetic predicates. [Earlier cut for always-sound verification](../sampled_slots/seed_cuts.md)

For local tags combined by addition alone, approximate additivity rules out almost all pruning regardless of the tag width over inappropriate output groups. [Additive tags on singleton seeds](additive_tags.md) · [Finite check](../../checks/check_candidate_fingerprint.py) · [Recorded output](../../evidence/check_candidate_fingerprint.txt)
