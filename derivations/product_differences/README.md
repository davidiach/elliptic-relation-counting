# Equality of the products modulo $p^2$

**21 November 2026.** We can decide equality from the sampled seeds without constructing either lifted determinant product. On a regular singleton, moving one slot according to one address bit makes the products equal for bit zero and unequal for bit one. A balanced calculation now decides the literal equality on regular empty and multiple-witness grids as well. Its log cost is at most $b/2+O(b/m+m+\log b)$, where $b=\lceil\log_2p\rceil$. This improves the full-grid evaluation of the equality statistic. It gives no improvement over the project's best conditional general-prime log-time bound $(24.00245996+o(1))b/(\log_2 b)^{1/4}$. [Previous precision-two construction](../direct_parity/deformation.md) · [End-to-end accounting](../sampled_slots/README.md)

## What the equality answers

The [subsequent half-list analysis](../implicit_chords/README.md) gives a bounded-storage version of this exact calculation and determines which poles survive cancellation for a raw abscissa bit.

Let $h$ count the ordinary target tuples and let the endpoint be $-T$. Use the same original lifts in both products. Speeds are measured in the invariant tangent direction, and the speed at each row is the difference between the two calls. On a regular grid:

| Fiber | Exact equality condition | Computation |
|---|---|---|
| Empty | A sum of logarithmic derivatives is zero in $\mathbb F_p$. Empty inputs can give either answer. | Split the rows into two blocks and sum elliptic chord slopes using multipoint polynomial evaluation. |
| Singleton | The row speeds on the unique tuple sum to zero in $\mathbb F_p$. | A signed group join locates the unique match. With just one bit moving, its speed is that bit. |
| Two or more | Both products vanish modulo $p^2$. | Detect a second match. |

The balanced calculation takes $(2w)^{\lceil m/2\rceil}\operatorname{poly}(m,b)$ expected bit operations from the seed words, including expansion of the slots. For a raw abscissa bit, the calculation also identifies a rational function in the target abscissa whose pole residue on a regular singleton is twice the requested bit. On typical bit-one roots under the stated probability bound, explicitly writing its reduced denominator requires nearly $2^b$ coefficients in a dense representation. A single-target calculation need not construct this denominator. [Proof and product-free evaluator](chords.md)

## Three other approaches

* **How hard is an equality oracle?** At the project widths, target existence on a nonempty regular mixed grid with $2T\ne O$ reduces deterministically to at most two equality comparisons with speeds constant within each slot. The proof uses $m\ge2$, $n_i<p$ and $p\nmid(m+1)$. The individual slot data choose the shifts in $O(mw)$ field operations. When the initial coefficients all vanish, pinning one free bit of a fixed-bit slot yields two usable comparisons. This prices the preparation and leaves the comparisons themselves to be evaluated. [Calibration and exceptional examples](fiber.md)
* **One ordinary incidence for a full parity query.** Over $\mathbb F_{p^2}$, a checked marker point can encode the local address bits so that exactly the even-parity relations survive one ordinary target test. It keeps the same target and number of slots, and enlarges just one slot by a factor of at most $\lfloor m/2\rfloor+1$. Eager preparation of the transformed slots costs $O(mw\,\operatorname{poly}(m,b))$. They are translated lists, so the inherited bit-mask contraction bound does not automatically apply. An individual bit can already be selected by pinning its original slot over the base field. [Marker construction and costs](projection.md)
* **A scalar over the base field.** Exact equality on a singleton is membership of a specified polynomial in the ideal cut out by the slot equations. A one-dimensional quotient suffices as output. With a consistent scalar evaluator for that quotient, $m$ evaluations on compact inputs recover the witness section with failure probability at most $1/p$ in the stated degree regime. Explicit reduction already involves a space of dimension $\prod_{i\ne j}|F_i|$ on omitting slot $j$. [Top-degree quotient and tests](local_algebra.md)

The open computation is a target-specific test of the restricted group-sum distribution, or an equivalent implicit evaluation of the graded scalar. An algorithm may use the full arithmetic of every seed. The current balanced algorithm still pays for the individual half tuples; typical successful roots contain $2^{b/2-o(b)}$ of them.

## Evidence

The [chord check](../../evidence/check_chord_difference.txt) makes 470 literal comparisons over two primes with three or six sampled slots, including agreements and disagreements on empty targets. It checks the full normalized field digit on empty and singleton inputs, and exceptional abscissas on a third curve. The [graded check](../../evidence/check_difference_algebra.txt) makes 130 literal comparisons and tests the empty/singleton pairing. The [marker check](../../evidence/check_difference_projection.txt) tests 640 determinants over quadratic fields, including a base curve with fully rational two-torsion. The [calibration check](../../evidence/check_difference_fiber.txt) exercises the exceptional fixed-bit split on grids with zero, one and two witnesses. The asymptotic bounds use the proofs above; these finite fixtures check the identities.
