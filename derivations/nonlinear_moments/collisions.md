# A two-challenge test for nonlinear moment readouts

Two independent queries on a marked seed can measure whether a vector of signed moments carries a recoverable address signal after a random binary encoding. For the raw field values of low-order moments, we can prove that this signal is too small on the masked-target distribution. The bound also applies to sparse lookup tables picked separately for each seed and target. Dense arithmetic readouts remain possible.

## What two challenges measure

Let $v$ be the binary address of a legal tuple, $r$ a uniform binary challenge of the same length, and $\chi_v(r)=(-1)^{r\cdot v}$. Fix the seed and target. Let $Y(r)$ be any deterministically computed vector of moment values or bins. For each observed response $y$, put

$$
a_y(v)=\mathbb E_r[\chi_v(r)1_{\{Y(r)=y\}}],\qquad C(v)=\sum_y a_y(v)^2 .
$$

Draw independent $r,r'$. Expanding the square gives

$$
C(v)=\mathbb E_{r,r'}\!\left[\chi_v(r+r')\,1_{\{Y(r)=Y(r')\}}\right].
$$

If the two challenges collide on the response, their relative parity counts with its sign. An ordinary collision probability would omit that sign and measure something else. The energy is nonnegative. Parseval applied to each response indicator gives $\sum_v C(v)=1$, where this sum includes *every* binary address.

Give the response an injective binary encoding and return the score $h(Y)=(-1)^{u\cdot\mathrm{enc}(Y)+c}$ with a fresh uniform key $(u,c)$ held fixed across the challenges of one recovery attempt. Distinct response labels have uncorrelated signs. Hence $\mathbb E_h|\widehat{h\circ Y}(v)|^2=C(v)$. No response table has to be enumerated to compute a query. At the other extreme, a seed- and target-specific lookup knowing the witness, bounded by one and depending on the challenge only through $Y$, can attain absolute coefficient $\sum_y|a_y(v)|$. With at most $L$ distinct responses this upper bound is at most $\sqrt{L C(v)}$. Thus random hashing and an optimally chosen lookup can have very different advantages when the response alphabet is large.

## Transferring a measured signal to fresh roots

Fix a target $T$ with $2T\ne O$. Write $X$ for the count of guarded tuples and $\mu=\mathbb E X$. A guarded plant draws a seed and one marked guarded tuple by giving each valid seed–tuple pair mass proportional to the seed probability. This law can be sampled by drawing a uniform ordered affine tuple of sum $T$ satisfying the guards, then sampling each slot conditioned on its point. Independence and the constant marginal point-inclusion probability give the stated weights. This procedure knows the marked address. It can estimate $\Gamma_T=\mathbb E_{\rm plant} C(v)$ directly with the two-challenge test, without having to establish uniqueness on the sampled root. An additive-$\epsilon$ estimate at any fixed confidence costs $O(\epsilon^{-2})$ independent bounded trials, plus the cost of sampling each pair and evaluating both responses.

Suppose a known lower bound $\Gamma_T\ge\gamma>0$ holds for each needed target. Under guarded planting and an independent score key the expected squared coefficient of the mark is at least $\gamma$. Its squared coefficient is at least $\gamma/2$ with probability at least $\gamma/2$. On fresh seeds and keys let $W$ count guarded tuples above that threshold. The planting identity gives $\mathbb E W\ge\mu\gamma/2$. Since $W\le X$, Cauchy–Schwarz gives $\Pr(W>0)\ge\mu^2\gamma^2/(4\mathbb E X^2)$. The thinned mixed sampler has $\mu\asymp\lambda$ and $\mathbb E X^2=O(\lambda)$. The fresh probability is thus $\Omega(\lambda\gamma^2)$. Walsh list decoding with threshold $\sqrt{\gamma/2}$ and terminal verification recovers a tuple with constant conditional probability on that event. This permits multiple valid tuples and charges every retry. [Isolation and second moments](../sampled_slots/isolation.md) · [Verified Walsh recovery and call costs](../sampled_slots/address_parity.md)

For example, if $m\asymp\sqrt b$, a response of cost $2^{O(\sqrt b)}$ and planted energy $\gamma\ge 2^{-O(\sqrt b)}$ would suffice for the project's square-root-in-bitlength goal under the inherited reduction. Forming the response, encoding it, charging seed setup and applying the same selection rule in training and queries are part of the requirement. A positive energy on one small curve establishes none of these growing-parameter guarantees.

## A bound for raw field values on the actual curve

Let $T$ be uniform on a fixed set of $n$ nonidentity targets, independently of the seed. For the outer comparison on a large subgroup $\log_2 n=b-O(1)$. Let $U$ be the unique guarded event with no additional ordinary hit and $q=\Pr(U)$. On $U$ let $v_*(T)$ denote the witness address. The following bound is unconditional over the seeds and targets:

$$
\mathbb E[1_U C_{\sigma,T}(v_*)]\ \le\ (1+D)/n.
$$

Here $D$ bounds the number of extra target collisions of any two formally different response vectors at fixed seed and challenge pair. For the vector of **all raw signed moments through order $K$**, one may take $D=(m+1)K$, independently of the number of included orders. Indeed, the determinant is linear in its endpoint evaluation row. The $j$-th signed moment, viewed as a function of $T$, has only a pole at the identity, of order at most $j(m+1)$. If two function-valued vectors differ, at least one component difference is nonzero, so it vanishes on at most $D$ targets. This argument permits seed- and challenge-dependent point weights independent of $T$. The same proof works for other publicly computed response vectors satisfying the collision bound. [Target pole calculation](adaptive_equalities.md)

For the proof, fix a seed and group challenges by equality of their entire vectors as functions on the curve. Compute the corresponding energy $C_{\rm gen}(v)$ using these formal labels. Parseval gives $\sum_v C_{\rm gen}(v)=1$. Distinct singleton targets have distinct addresses, since the points of any fixed address have one sum. Thus $\sum_{T:U}C_{\rm gen}(v_*(T))\le1$. Specializing to a target can merge distinct labels. For every $v$ the absolute change in its energy is at most the probability over $r,r'$ of an extra collision at that target. Summing over all targets costs at most $D$ by the pairwise pole bound. Divide by $n$ and average over seeds. The proof even permits choosing one valid address separately at each represented target in place of the singleton.

Three consequences explain which kinds of nonlinear decoding this rules out.

* **Hashing the entire vector to a bit.** The conditional mean squared witness coefficient averaged over independent keys is at most $(1+D)/(nq)$. At the isolating widths and $D=2^{o(b)}$, the mean absolute coefficient is at most $2^{-b/2+o(b)}$. This is a bound on the randomized scheme. A purpose-built deterministic map can exploit arithmetic that a random encoding loses.
* **A target-specific dictionary of exceptions.** Allow an arbitrary constant default and changes at up to $L$ raw response values. Choose the table using the full seed and target, even with knowledge of the witness. The deviations have absolute value at most two, so the conditional mean absolute coefficient is at most $1/(nq)+2\sqrt{L(1+D)/(nq)}$. The constant default contributes only at the all-zero address, which can be a witness for at most one target per seed. The dictionary and default are fixed across challenges. This rules out sparse lookup even with uncharged table construction when $L,D=2^{o(b)}$.
* **Combining separately hashed moments nonlinearly.** Give each of $k$ component functions an independent random binary hash with pairwise independent outputs on distinct values. The readout may compute *any* bounded function of those $k$ bits, the seed and the challenge. Allow target access through the hashed values. For two challenges, the joint hash bits have the same law as hashes on formal function labels unless a component has an extra target collision. If the component degree bounds are $D_i$, their total charge is $\sum_i D_i$. The generic score for any fixed formal hash assignment is target-independent and has total squared Walsh energy at most one. The difference of the two expected products of scores has magnitude at most two at each exceptional target. Consequently the unconditional squared coefficient is at most $(1+2\sum_i D_i)/n$, averaged over random keys. The readout accesses keys only through the queried hash outputs. Revealing raw field values or selecting keys using labeled training changes the model.

The singleton upper bounds price a typical success in the specified masked-target experiment. They do not assert targetwise bounds for every $T$, and the planted lower-bound criterion concerns a different, explicitly weighted law. In particular, a low-order scalar expression acquired by repeated squaring can yield a character test of very high pole degree. A dense, cheaply evaluated field map remains a possible route.

## Finite checks

On the order-103 curve over $\mathbb F_{101}$ the checker evaluates the joint first and second moments for all 64 challenges and all 102 nonidentity targets on one mixed seed. There are 57 distinct function-valued responses. A pair of formally distinct responses agrees at up to eight targets, attaining the degree bound. The sum of the 22 singleton energies is $1444/4096$, below the theoretical upper bound of nine.

On the curve over $\mathbb F_{13}$ it separately enumerates all 8,192 weighted seeds at each of two targets, all 64 challenges, and all 256 binary scoring tables up to a global sign on the nine joint quadratic-character bins. This gives the exact optimum mean absolute coefficient even over bounded real lookups on those bins, since the absolute-coefficient objective is convex and a maximum occurs at a cube vertex. The separate witness-minus-wrong optimization reported by the script is over *binary* tables only. The larger gain occurred at $2G$: the best joint mean absolute witness coefficient was $0.144509$, against $0.119782$ for the best lookup using either character alone and $0.103856$ for the same joint rule on wrong legal addresses. At $G$ these three values were $0.102992$, $0.092420$ and $0.103175$. These are finite target-specific optima. [Executable population and collision checks](../../checks/check_nonlinear_readout.py)
