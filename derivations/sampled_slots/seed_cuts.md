# The cost of keeping each slot seed local

Fix a target $T$ in a finite abelian group. Independent seeds $\sigma_i$ select point sets $F_i$. Assume each prescribed point $P$ satisfies $\Pr(P\in F_i)\le\rho_i$. In the mixed sampler, $\rho_i=w/2^b$. For a fixed target, a sound verifier must produce points $P_i\in F_i$ with $\sum_iP_i=T$. It may deliberately skip roots. Write $\kappa>0$ for its fresh probability of accepting. Under the following locality conditions, each hidden slot forces a communication cost even when the verifier has arbitrary topology and unlimited local computation.

## A single cut

View a Boolean feasibility network as a conjunction of local predicates. An assignment to its edge labels and private variables is a certificate. Fix all seeds except $\sigma_i$, and fix any public randomness independent of the seeds. Separate a region containing the private slot input and its point choice from the rest of the network. The region can contain many vertices. Assume the following.

1. Outside the region, the graph, its interface syntax, local tests, and point decoders can have at most $2^{a_i}$ different descriptions as $\sigma_i$ varies. A description may be selected by arbitrarily expensive joint preprocessing. For each description, at most $2^{c_i}$ boundary assignments are possible.
2. The outside assignment, the outside data and the boundary labels determine the signed points in all the other slots. For a fixed description and label, any feasible inside and outside assignments splice to give a globally feasible assignment.
3. Every such global assignment decodes a legal relation. This is soundness of certificates. Only existence on the accepted inputs is needed for completeness.

Then

$$
 \kappa\ \le\ \min(1,2^{a_i+c_i}\rho_i).
$$

For every possible outside description and boundary label choose one feasible outside assignment in a fixed order, skipping empty cases. Every chosen assignment gives a candidate $T-\sum_{j\ne i}P_j$. The list of candidates is fixed **before** drawing the hidden seed, even if the realized description will depend on that seed. If a global certificate exists, its inside part splices with the canonical outside part at its description and label. Soundness forces the corresponding candidate to belong to $F_i$. The union bound over the list proves the inequality, followed by averaging over the frozen data. If an outside label determines only abscissas, there are at most $2^{m-1}$ sign combinations for its candidate. This adds at most $m-1$ bits to the right side of the logarithmic budget.

## Summing the cuts

Take disjoint regions, one per slot. For fixed-width interfaces let $c_i$ be the sum of actual incident edge widths in bits, and $S$ the total edge width. Each edge is counted in at most two cuts. For fixed interface widths we obtain

$$
 2S\ \ge\ \sum_i\bigl(\log_2(\kappa/\rho_i)-a_i\bigr)_+ .
$$

For a varying graph, summing worst-case cut widths could mix different inputs. Instead, for each cut consider certificates on instances with an interface of at most $C$ bits. Provided its syntax is included in the outside description, at most $2^{C+1}$ binary strings of length at most $C$ can occur per description. The same proof gives

$$
 \Pr(\text{accepted and }c_i\le C)\le 2^{a_i+C+1}\rho_i .
$$

Choose $C_i=\lfloor\log_2(\kappa/\rho_i)-a_i-\log_2(4m)\rfloor$ and take a union bound over the nonnegative $C_i$. On at least half the accepted probability mass, all those cuts simultaneously exceed $C_i$. For the paid widths $\ell=b/m+O(1)$, with $m\to\infty$, $m=o(b)$, $\log_2(1/\kappa)=o(b)$, and $\max_i a_i=o(b)$, both conclusions give a budget of at least $(1/2-o(1))mb$ for the corresponding networks. With variable widths this is a statement about a constant fraction of successes. The possible sign loss above is lower order in this regime.

## Why public seeds change the conclusion

All hash seeds are public in the project. A node can use them in its local tests, in the choice of graph, or in the decoder for its labels. For example, pad the seed words to a common length and publish their bitwise XOR. The complement of each slot, which already knows the other words, recovers the entire missing word. One shared word of length $O(b)$ can thus reveal $\Theta(b)$ bits at *each* cut. Charging that length just once as though it bounded $\sum_i a_i$ would produce a false lower bound. The XOR example supplies no curve-equation check.

The cut theorem applies to seed-local sound certificates with locally decoded point choices. A fully public joint arithmetic test may fall outside all these cuts. Contractions using cancellation or a global search that chooses among unsound assignments also need separate analysis. The result does not imply a general lower bound on the time to find a relation.

[Exact finite cut check](../../checks/check_seed_cuts.py) · [Sparse sampling and the address of a unique relation](isolation.md)
