# An exact section-space computation and its first informative degree

There is a direct geometric calculation of $M_0-M_{p-1}$ that uses no scalar moment recurrence. It counts projective sections through the target. On a regular isolated input the desired section is a single reduced point. Its coordinate ring admits a one-dimensional kernel that recovers every address sign. We determine the **first degree** at which this homogeneous kernel appears, including on positive inputs. Its ambient dimension there is already the product of all but one slot size.

## Dualizing the zero test

Write $d=m+1$, take the usual evaluation row $e(P)$ for $L(dO)$ with constant coordinate first, and let $F_i$ be the signed points retained in slot $i$. The nonidentity target is $T$. Sections through $-T$ up to scaling form $\mathbb P^{m-1}$. In coordinates $z$ on the underlying $m$-dimensional vector space, evaluation at $P$ gives a linear form $\ell_P(z)$. One choice has coefficients $e_j(P)-e_j(-T)$ at each nonconstant basis function $j$. For a tuple, the determinant of its $m$ linear rows differs from the original determinant by at most a sign.

Give each slot point its Walsh sign $\epsilon_i(P)$ for the challenge. Over $\mathbb F_p$,

$$
M_0-M_{p-1}
=\sum_{[z]\in\mathbb P^{m-1}(\mathbb F_p)}
 \prod_{i=1}^m \ \sum_{P\in F_i:\,\ell_P(z)=0}\epsilon_i(P).
$$

Indeed, expand the product. A tuple with nullity $a>0$ is counted at $1+p+\cdots+p^{a-1}$ projective points, congruent to one. A full-rank tuple contributes zero. This identity handles repeated rows and empty slots as well. It exchanges the tuple grid for the projective space and by itself gives no fast algorithm. [Determinant and address convention](../nonlinear_moments/signed_spectrum.md)

## The kernel for a regular singleton

Assume the slots are nonempty, have disjoint abscissas between slots and avoid $x(T)$. Set $n_i=|F_i|$ and $Q_i(z)=\prod_{P\in F_i}\ell_P(z)$. The common projective points of all $Q_i$ are in bijection with full relations. In fact, a nonzero section in $L(dO)$ through $-T$ and one point in each slot vanishes at $d$ distinct points. They exhaust its zeros and sum to the identity. Conversely such a tuple determines that section up to scalar. In particular it cannot vanish at any other point from these slots.

Keep the first $m-1$ slot equations in the homogeneous coordinate ring

$$
R=k[z_1,\ldots,z_m]/(Q_1,\ldots,Q_{m-1}),
\qquad s=\sum_{i<m}(n_i-1),\qquad L=\prod_{i<m} n_i.
$$

Here $k$ may be the algebraic closure for the proof. The ranks of the finite-field linear systems used below are unchanged by scalar extension. Every partial choice together with $-T$ prescribes $m$ distinct zeros. The remaining divisor has degree one, so on a genus-one curve Riemann–Roch gives a one-dimensional space of sections. Thus the chosen $m-1$ linear conditions in section space are independent and intersect in one projective point. Taking their union over all partial choices shows that the projective support cut out by the $Q_i$ is finite. These $m-1$ forms have height $m-1$ in a polynomial ring and form a regular sequence. Their quotient has dimension one and no contribution supported only at the irrelevant ideal. Partial choices that yield the same point may cause nonreduced structure. The length counted with multiplicity is $L$, since

$$
\operatorname{Hilb}_R(t)=
\frac{\prod_{i<m}(1+t+\cdots+t^{n_i-1})}{1-t}.
$$

Suppose now that there is precisely one full tuple, with projective point $z_*$. In each slot exactly its selected factor $\ell_{i,*}$ vanishes there. In a neighborhood of $z_*$, the first $m-1$ equations are these independent linear equations times units. Consequently the local ring on the projective intersection is a copy of $k$. At any other partial section $Q_m$ is a unit locally. In an affine chart containing the whole finite support, the kernel of multiplication by the dehomogenized $Q_m$ is precisely the one-dimensional summand supported at $z_*$. This also describes the kernel of $R_j\xrightarrow{\cdot Q_m}R_{j+n_m}$ restricted to degree $j$, since a linear chart coordinate avoiding the finite support is a nonzerodivisor in $R$.

The first $j$ with nonzero kernel is **exactly $s$**. For existence at $s$ take the product $S$ of all unselected factors of the first $m-1$ slots. Its value at $z_*$ is nonzero. Every $\ell_{i,*} S$ is zero in $R$ since it is a multiple of $Q_i$. The independent $\ell_{i,*}$ generate the point ideal, so $Q_m S=0$ in $R$.

For minimality, let $K$ be any nonzero homogeneous kernel element and let $u$ be the chart coordinate above. Since $K$ is supported at the reduced witness, all $\ell_{i,*}K=0$ in $R$. In $A=R/(u)$ the images of these linear forms generate the homogeneous maximal ideal, as $u(z_*)\ne0$. The image of $K$, if nonzero, must be in the socle, the subspace annihilated by that ideal. The Artinian complete intersection $A$ has its socle precisely in degree $s$. For a direct degree calculation, view $A$ as the quotient of a polynomial ring on $m-1$ variables by the regular sequence of degrees $n_i$. Its top Tor with $k$ is one-dimensional in degree $\sum_{i<m}n_i$ by the Koszul resolution of the sequence. Computing the same Tor by the Koszul resolution on the variables yields the socle with grading shifted up by $m-1$. If $\deg K<s$, we deduce $K=uK'$ in $R$. Since $u$ is regular, $K'$ is a nonzero kernel element in one smaller degree. Iterating gives a contradiction. This proof permits nonreduced partial sections away from the witness. [Complete intersections and the separator principle](https://archive.mpim-bonn.mpg.de/id/eprint/3773/1/preprint_1993_71.pdf)

At degree $s$, $\dim R_s=L$ and the kernel has dimension one. Compute any nonzero kernel vector $K$ there. For each slot, the sole retained point satisfying $\ell_P K=0$ in $R_{s+1}$ is its chosen point. Multiplying their Walsh signs returns $M_0-M_{p-1}$. The same kernel can answer every challenge. If there is no relation then $Q_m$ is everywhere a unit in the affine chart and all graded kernels are zero. Equivalently, the Hilbert function of the full homogeneous ideal on a singleton agrees with the formula for an empty projective intersection in every degree strictly below $s+n_m$. This follows by subtracting the rank of the multiplication map from the Hilbert function of $R$ in each degree.

## Actual sizes

The growth statement uses the thinned mixed sampler, with $\nu=Nw/2^b$, $\nu^m=N\lambda$, $\log_2(1/\lambda)=O(m+\log m)$, and $m\to\infty$, $m=O(\sqrt b)$ in its isolating regime. Conditional on a guarded singleton, with probability $1-O(m^{-2})$ the root is regular and every slot has at least $\nu/m^3$ points. The bound applies simultaneously, hence also after selecting which slot to drop. It gives

$$
\log_2 L \ge (1-1/m)b-O(m\log m).
$$

Representing a full vector or basis of the first informative graded quotient therefore exceeds the targeted budget on these roots. A separator at this degree has a product representation of size $O(\sum_i n_i)$ **given the chosen factors**. Recovering or testing it implicitly is a distinct open problem. The degree and dimension calculation is not a lower bound on arbitrary symbolic algorithms. [Conditional slot-size estimate](../nonlinear_moments/signed_spectrum.md)

### Thinning before forming a quotient

One might try to make many smaller quotients and stop at the first hit. Here is a limited accounting bound on that strategy. Fix the slots, and sample $T$ uniformly and independently from a set of $n$ possible targets. Suppose subgrid selection has access to a target-dependent label with at most $J$ values. Beyond this label, its only target-dependent feedback is a bit saying whether the queried subgrid contains a tuple summing to $T$. Its random coins are independent of $T$. For each query with nonempty subslots $G_i\subseteq F_i$ and a freely chosen omitted index $j$, charge $L'=\prod_{i\ne j}|G_i|$, the graded quotient width. If $h=\max_i |F_i|$, then the probability of seeing a first hit at cumulative charge at most $B$ is at most $JhB/n$.

To see this, fix the coins and one label. Follow the path on which all preceding answers are negative. A queried grid on this path contains at most $hL'$ tuples and hence covers at most $hL'$ targets. Every possible first hit under this label occurs on this path. Unioning these targets for all queries through budget $B$ and all $J$ labels gives the bound. Conditioning on any original-root singleton event of probability $q>0$ multiplies the upper bound by at most $1/q$. The bound charges the indicated width on every query. It places no restriction on a method that computes a cheaper invariant on negative grids, or selects grids using richer target-dependent arithmetic.

The finite check constructs the graded rings by field row reduction. In the three-slot fixture the first two slots both have size four. Two different positive targets give respectively 15 and 16 distinct partial sections in length 16. In both cases the unique kernel first appears in degree six and identifies the address. An empty control and a four-slot positive fixture provide additional comparisons. [Executable check](../../checks/check_direct_algebra.py)
