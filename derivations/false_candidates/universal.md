# Exact signed incidence from a short list of abscissa tests

The [single-offset tree](target.md) meets the requested average false-address bound. A stronger address-level statement holds at the paid arity. For each sufficiently large prime curve one can prepare $m+2$ offsets such that the intersection of their tests is exact for **all** signed point tuples and **all** targets simultaneously. The same offsets can be reused for new hash seeds. This result also holds for nongeneric tuples, including repeated points or identities, with the projective conventions of the single-offset test.

For an offset $U$, the test on proposed signed points $P_i$ asks if the abscissas of $P_i+U$ can be combined to give the abscissa of $T+[m]U$. Equivalently some signs satisfy $\sum_i\epsilon_i P_i=T+[2j]U$, with $j$ the number of negative signs. Fix a wrong tuple, meaning $\sum_iP_i\ne T$. Thus every successful probe has $1\le j\le m$. There is no independence assumption about the signs in different probes.

## A checked common sequence

Set $d=m+1$, $h=\lceil d^{d/2}\rceil$, $L=2mh+1$, and
$$C=2mh\sum_{a=0}^{d} L^a.$$
Choose a rational point $R$ of order greater than $C$. Use $U_a=[L^a]R$ for $0\le a\le d$. All true tuples pass every test by taking positive signs.

Suppose a wrong tuple passes every test. For each probe form the integer row $v_a=(\epsilon_{a1},\ldots,\epsilon_{am},-1)$, with corresponding $j_a\ge1$. There are $d+1$ rows of length $d$. A minimally dependent subset gives an integer relation $\sum_a c_a v_a=0$ with $0<\max_a|c_a|\le h$. Indeed, restrict to enough independent columns to have full column rank on this circuit of rows, and take its signed maximal minors. Hadamard's bound on their absolute values is at most $d^{d/2}$. Coefficients off the subset are zero. Adding the group equations with coefficients $c_a$ now gives
$$[k]R=O,\qquad k=\sum_{a=0}^d (2j_a c_a)L^a.$$
The lowest nonzero coefficient is nonzero modulo $L$ because its absolute value is at most $2mh=L-1$. Hence $k\ne0$, while $|k|\le C$. The order of $R$ excludes this possibility.

This proof certifies distinct **addresses**. The exact intersection might still admit several internal abscissa assignments on a true address. Feasibility with pinned address bits suffices to read a unique address. Counting internal assignments and declaring a unique address only when that count is one would require an extra argument.

## Uniform preparation and cost

One can check $[j]R\ne O$ for $1\le j\le C$ by successive additions. Over $\mathbb F_p$ at most $\sum_{j\le C}j^2=O(C^3)$ rational points fail, since the kernel of multiplication by $j$ has at most $j^2$ geometric points. At $\log C=O(m^2\log(m+1))=o(b)$, random rational points pass with probability $1-o(1)$. Once checked, the answer for every subsequent input is deterministic. At $m=\Theta((\log b)^{1/4})$, even this simple preprocessing costs only $b^{o(1)}$ group operations, besides sampling and the bit cost of each operation. No factorization of the curve order is needed.

The direct graph with separate tests has $(m+2)(m-2)$ internal full-field abscissas before its leaf and arithmetic costs. A common-address pin still needs joint feasibility over all tests. Separate existence calls on the probes could return different addresses. Thus zero false addresses, even uniformly over all seeds, does not supply a new end-to-end exponent.

If arithmetic over $\mathbb F_{p^2}$ is allowed, a single probe also gives exact signed incidence. Choose the [checked quadratic-twist marker](../product_differences/projection.md) with $[2j]R\notin G$ for $1\le j\le m$. A branch with $j$ negative signs would require $[2j]R=\sum_i\epsilon_iP_i-T\in G$, a contradiction. The extension-field abscissas and the cost of their arithmetic are then part of the representation. The rational construction above gives a universal choice while keeping every probe in the base field.

The checker exhausts all point tuples and targets on an order-19 elliptic curve for a four-offset fixture. At arity three it tests all 16,807 possible joint nontrivial sign patterns by linear consistency. On the order-1,019 curve over the field with 1,031 elements, the five multiples $1,97,238,668,599$ of the generator already give an exact test for all triples and targets. This small instance is certified directly by the rank checks. The same check also runs on an abstract cyclic group large enough for the general order bound. [Checker](../../checks/check_candidate_offsets.py) · [Transcript](../../evidence/check_candidate_offsets.txt)
