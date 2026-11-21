# Balanced shifts and an exact existence comparison

Move the sampled points in each slot at one constant invariant tangent speed and the evaluated endpoint at the negative sum of those speeds. On a regular grid, the change of the full determinant product modulo $p^2$ factors as its unperturbed value modulo $p$ times an elementary sum over the individual slots, multiplied by $p$. That sum can always be made nonzero using at most two subgrids. Consequently the equality-only oracle for these constant shifts already tests target existence. Computing that last equality directly at the desired cost remains open.

## An elliptic alternant identity

Write $d=m+1$ and label the evaluated rows by $P_1,\ldots,P_m,E$, where $E=-T$. Let $\Delta$ be their determinant in a basis of $L(dO)$. The regularity assumption makes all row abscissas distinct. For each row let $D_i$ be differentiation in the invariant tangent direction $(2y,3x^2+A)$, dual to $\omega=dx/(2y)$. Define

$$
K(P,Q)=\frac{y_P+y_Q}{x_P-x_Q},\qquad K_i=\sum_{j\ne i}K(P_i,P_j).
$$

The following rational identity also holds at all regular determinant zeros:

$$
(D_i-D_j)\Delta=(K_i-K_j)\Delta .
$$

Here is an algebraic proof in characteristic $p>3$. Fix the other $d-2$ generic points $R$ and the generic sum $C=P_i+P_j$. Vary $P=P_i$ and write $Q=C-P$. The alternant is nonzero generically by the divisor-sum criterion. It vanishes at $P=R$, at $Q=R$, and at all four solutions of $P=Q$; multiplication by two is separable of degree four. These conditions give $2d$ distinct zeros. Only $P=O$ and $Q=O$ can be poles, of order at most $d$ each, so both orders are exactly $d$ and all the zeros are simple.

The differential $(K_i-K_j)\omega$ is

$$
\left[2K(P,Q)+\sum_R(K(P,R)-K(Q,R))\right]\omega.
$$

It has the same simple-pole residues as $d\Delta/\Delta$. At $P=R$, $K(P,R)\omega$ has residue one. At $Q=R$ use $\omega_Q=-\omega_P$. At $P=Q$, with a local invariant parameter for $P$, the numerator of $K(P,Q)$ begins with $2y_P$ while the derivative of its denominator is $4y_P$, giving residue one for the doubled term. At $O$, $K(P,R)\omega$ has residue $-1$ by the pole expansions of $x,y$. The doubled term and the $d-2$ single terms yield $-d$ there and also at $Q=O$. The apparent pole with $P=-R$ is removable. Thus the difference of the two differentials is holomorphic and is a constant multiple of $\omega$. It is invariant under exchanging $P,Q$, whereas $\omega$ changes sign. The constant vanishes. Clearing the denominators extends the rational identity to all regular tuples, including specializations with rational two-torsion. This logarithmic identity is also consistent with the classical Frobenius–Stickelberger alternant. [Frobenius and Stickelberger, translated original article](https://arxiv.org/abs/2603.27466)

## The part of a general comparison that comes from individual slots

Use the same curve-point lifts in both products. Let $c_i(P)$ and $e$ be the differences of the slot and endpoint speeds between the calls. For a tuple $a$ the first-order change of its determinant, divided by $p$ and reduced modulo $p$, is

$$
\dot\Delta_a=\Delta_a\sum_i c_i(P_i)(K_i-K_E)
                  +\big(e+\sum_i c_i(P_i)\big)D_E\Delta_a .
$$

With $V=\prod_a\Delta_a\bmod p$, the difference of the products divided by $p$ is therefore

$$
V C_{\rm loc}+
\sum_a\big(e+\sum_i c_i(P_i)\big)D_E\Delta_a\prod_{b\ne a}\Delta_b,\qquad
C_{\rm loc}=\sum_a\sum_i c_i(P_i)(K_i-K_E).
$$

This polynomial product rule is valid regardless of the number of zero determinants. For an empty slot both lifted products are the empty product, so consider nonempty slots of sizes $n_i<p$. Let $D=\prod_i n_i$ and $L_j(z)=n_j^{-1}\sum_{Q\in F_j}(z-x_Q)^{-1}$. Inversion symmetry in the slots cancels the $y_Q$ terms when summing $K(P,Q)$. Consequently

$$
\frac{C_{\rm loc}}{D}
=\sum_i\frac1{n_i}\sum_{P\in F_i}c_i(P)
 \left[2K(P,E)+\sum_{j\ne i}\big(y_P L_j(x_P)-y_E L_j(x_E)\big)\right].
$$

This costs $O(m^2w^2)$ elementary field operations directly. Alternatively, build $g_j(X)=\prod_{Q\in F_j}(X-x_Q)$ with multiplicity, so $L_j=g_j'/(n_jg_j)$. Product trees and multipoint evaluation at the other slot abscissas give $\widetilde O(m^2w)$ field operations. Regularity keeps all evaluations away from poles. This computes the local part for arbitrary address-dependent speeds, leaving just the endpoint derivative weighted by the sum of all speeds.

## Choosing the constant shifts on every nonempty regular grid

Set $c_i(P)=s_i$ constant and choose the endpoint speed $e=-\sum_i s_i$. The collective term vanishes on every tuple. Put

$$
a_i=\frac1{n_i}\sum_{P\in F_i}(x_P-x_E)^{-1},\qquad S=\sum_i a_i.
$$

The unweighted mean of $K(P_i,P_j)$ across two distinct slots vanishes by simultaneous sign reversal. The unweighted mean of $K(P_i,E)$ equals $y_Ea_i$. Thus

$$
C_{\rm loc}=D y_E\sum_i(s_i-e)a_i=D y_E\sum_i s_i(a_i+S).
$$

If $2T\ne O$ and $p\nmid m+1$, some basis direction with one $s_i=1$, the others zero and $e=-1$ gives a nonzero coefficient unless every $a_i=0$. Indeed, if all $a_i+S=0$, summing gives $(m+1)S=0$. The coefficients can be computed explicitly in $O(mw)$ field operations after forming the slots.

In the exceptional case choose **one fixed-bit slot** in the mixed sampler ($m\ge2$). Its canonical abscissa has the integer expansion $x_P=C+\sum_j2^{b_j}u_j(P)$, where $C$ is the fixed part and $u_j$ are its remaining free address bits. Summing $(x_P-x_E)/(x_P-x_E)$ gives

$$
0\ne n_i=(C-x_E)\sum_{P\in F_i}(x_P-x_E)^{-1}
 +\sum_j 2^{b_j}\sum_{P\in F_i:\,u_j(P)=1}(x_P-x_E)^{-1}\pmod p .
$$

The first term is zero. Some single free bit therefore has a nonzero reciprocal sum on its one branch; the zero branch has the opposite sum. Pin that bit to zero and to one in two subqueries. Both are ordinary fixed-bit masks, both subgrids are nonempty, and each now has a nonzero coefficient for some basis direction. A fully pinned nonempty slot has just one abscissa and cannot have $a_i=0$. To find the bit in $O(w)$ field operations, put the reciprocal contributions at the leaves of a binary tree indexed by free addresses. Compute all subtree sums, and at each level add the right-child sums to get the total for that bit. The project permits pinning any free bit; a model restricted to contiguous prefix intervals would require that extra capability.

For a calibrated subgrid the equality of the shifted and original products modulo $p^2$ is precisely the event $V=0$, equivalently the existence of an ordinary target tuple. The equivalence covers empty witness fibers, singletons and arbitrarily many witnesses. Testing one or two calibrated subgrids decides the original existence question. On a regular guarded singleton, successive sign-symmetric abscissa pins therefore require $O(m\log w)$ such comparisons to identify all abscissas, followed by at most $2^m$ sign checks. All comparisons use the same original curve lift and constant speeds at each slot and endpoint. The expensive primitive is still equality of the two full products. The computable local sum only calibrates the comparison.

For fixed nonempty slots and uniform independent target the exceptional case itself is sparse. For any one slot the rational sum $\sum_P (x_P-z)^{-1}$ has one simple pole at each distinct abscissa with nonzero residue of magnitude one or two. Its numerator has degree at most $w-1$ and is nonzero. Thus at most $2(w-1)$ affine target points off the poles make all $a_i=0$. Exceptional fixed targets do occur. The finite check includes regular width-two grids with zero, one and several witnesses for which **all** balanced coefficients vanish before the split. [Check](../../checks/check_difference_fiber.py) · [Output](../../evidence/check_difference_fiber.txt)
