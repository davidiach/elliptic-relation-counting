# A binary group-algebra calculation for one witness bit

An isolated signed address has a parity bit even though the two signs $+1$ and $-1$ coincide in characteristic two. The missing bit can be stored as a derivative. This gives an exact simultaneous calculation from the sampled slots. We determine how much information one class of composable binary-algebra implementations must retain, including for even-order curve groups.

## The address derivative

Put $G=E(\mathbb F_p)$, $N=|G|$, and work in $\mathbb F_2[G][u]/(u^2)$, writing $[P][Q]=[P+Q]$. If a challenge assigns the bit $c_i(P)=r_i\cdot a_i(P)$ to a signed point in slot $F_i$, form

$$
\prod_i\sum_{P\in F_i}(1+u c_i(P))[P]=A+uB .
$$

At any target $T$, $A(T)$ is the parity of the number of tuples summing to $T$. The bit $B(T)$ is the XOR of the challenge bits of all those tuples, since the linear coefficient of a product of dual numbers is the sum of the individual perturbations. On a singleton with witness $a_*$, this gives $B(T)=r\cdot a_*$, and on a target-regular singleton $(-1)^{B(T)}$ is exactly $M_0-M_{p-1}$. On general inputs it computes an unguarded group-incidence quantity. [Singleton and regularity conventions](../nonlinear_moments/signed_spectrum.md)

Even one fixed challenge yields a dimension requirement for an exact implementation of this identity. Take the high bit of the **raw affine anchor address** and zero bits in all other slots. If the width is $w=2^\ell\ge2$, this bit selects the upper half $[w/2,w)$ of the retained interval. Here $B$ is simply the group-algebra product of that restricted anchor and the ordinary indicators of all remaining slots.

## Which binary vectors occur in the span of the actual products?

For each abscissa word $x$ among the $M=2^b$ words, let $e_x$ be the sum of its rational affine point symbols, or zero at a nonliftable or out-of-field word. Write $W$ for the span of all distinct pairs $[P]+[-P]$ with $2P\ne O$. These pairs have disjoint support, so $\dim W=(N-h)/2$, $h=|G[2]|\le4$. Their span consists of inversion-symmetric functions on $G$ which vanish on $G[2]$. It is closed under convolution: at an inversion-fixed output, pairs of nonfixed input points contributing to a product match their simultaneous negatives and cancel.

Consider only the permitted affine multiplier $1$. Varying its translation supplies every cyclic interval sum $I_c=\sum_{j=0}^{v-1}e_{c+j}$ of length $v=w/2$, where word addresses wrap modulo $M$. Consecutive sums obey $I_c+I_{c+1}=e_c+e_{c+v}$. If **every residue modulo $v$ contains some word with $e_x=0$**, the telescope along each residue recovers every $e_x$ in the binary span of these masks. Without any missing-word assumption the codimension of the cyclic interval span in the unprojected word space is at most $v-1$: in its orthogonal space the word is $v$-periodic and its sum on one period vanishes.

The missing-word condition holds uniformly for sufficiently large odd primes at widths $v=o(\sqrt p)$. If a residue contains no out-of-field word, its canonical abscissas give an arithmetic progression of length at least $\lfloor p/v\rfloor$. Write $f(X)=X^3+aX+c$ for the squarefree curve cubic and $\chi$ for the quadratic character. In any progression of step $v\ne0\bmod p$, a run of $L$ consecutive values with $\chi(f(x))=1$ has $L=O(\sqrt p)$. To see this, put $a_j=\chi(f(x_0+vj))$, take $H=\lfloor\sqrt p\rfloor$, and sum the squares of all cyclic length-$H$ block sums over the full field. The diagonal contributes at most $Hp$, and each off-diagonal correlation is bounded in absolute value by $5\sqrt p$ by the complete quadratic Weil bound. Indeed the two translates of the squarefree cubic have different root sets: a nonzero translation has order $p>3$ and cannot stabilize a set of three geometric roots. Their product is therefore nonsquare. A positive run of length $L\ge H$ contributes at least $(L-H+1)H^2$, giving $L\le H-1+p/H+5\sqrt p$. There are at most three zeros of $f$ in the original progression, so if it has no nonsquare it splits into at most four positive runs. This contradicts $p/v\gg\sqrt p$. The project's growing-arity widths have $w=2^{b/m+O(1)}=o(\sqrt p)$. [Weil, *On some exponential sums*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1079093/) · [Sampling widths](../sampled_slots/isolation.md)

For each fixed-bit slot, summing its indicators over all fixed-bit patterns gives $J+[O]$ where $J=\sum_{P\in G}[P]$, because these masks partition all rational affine points. Every $f\in W$ has augmentation zero, so $fJ=0$ and $f(J+[O])=f$. Expand the product over the independent choices of all the fixed-bit slots. Every $f\in W$ is thereby a binary linear combination of actual $m$-slot outputs $B$, with the stipulated anchor challenge and any prescribed free-bit layouts. The argument concerns the span over all seeds. Individual seed outputs need not belong to $W$ when rational two-torsion is present.

## Fixed-target dimension of a composable sketch

Suppose the calculation applies an $\mathbb F_2$-algebra homomorphism $\pi$ to the group distributions, multiplies there, then uses one $\mathbb F_2$-linear functional of $\pi(B)$ to recover $B(T)$ exactly on every seed. Both the homomorphism and functional may depend on the fixed target, layouts and curve, and must be the same for all seeds under examination. The domain need only include the algebra generated by the relevant masks. It contains $W$ by the span argument. Linearity shows the readout agrees with $f\mapsto f(T)$ on all of $W$.

For $f\in W\cap\ker\pi$ and every $g\in W$ we now have $(fg)(T)=0$. The radical of the bilinear form $\beta_T(f,g)=(fg)(T)$ can be calculated without semisimplicity. Testing a pair generator $g=[P]+[-P]$ gives $f(T-P)=f(T+P)$. This condition holds automatically at fixed $P$, so it is precisely reflection invariance about $T$ on the full group. Together with inversion invariance it is equivalent to invariance under translation by $2T$.

Write $H=\langle 2T\rangle$, $K=G/H$, $k_2=|K[2]|$, and $h'=|\operatorname{im}(G[2]\to K)|$. Radical vectors are exactly the pullbacks of inversion-symmetric functions on $K$ that vanish on the image of $G[2]$. Their dimension is $(|K|+k_2)/2-h'$. It follows that

$$
\dim\pi(W)\ \ge\ \operatorname{rank}\beta_T
 =\frac{N-h-|K|-k_2}{2}+h'.
$$

In fact $k_2\le2h'$. Sending $P+H\in K[2]$ to $2P+2H\in H/2H$ has kernel the image of $G[2]$, and the cyclic group $H$ has $|H/2H|\le2$. Thus for **every $2T\ne O$**, the dimension bound is at least $N/4-2$. For odd-order $G$ and target order $t$ the exact rank is $(N-N/t)/2$. In particular, for a nonzero target on a prime-order curve the sketch must be injective on all $(N-1)/2$ independent pair coordinates.

This rules out reducing the **dimension of an explicit composable binary-algebra state** below $2^{b-O(1)}$ by selecting a few Fourier factors, bundling conjugates, or using a nilpotent quotient with a linear readout. The dimension does not bound the time to evaluate one coordinate of an implicitly represented state.

## A binary rank bound across the seed split

For the same fixed target and high address bit, arrange the answers $B(T)$ into a zero-one matrix across any nonempty partition of the slots. The side containing the anchor spans at least $W$ after summing its other slot masks. On the other side choose one fixed-bit slot. Its nonempty bins partition the affine points. Discard bins containing rational two-torsion and let $s$ be the number remaining, so $s\ge(N-1)/(2w)-(h-1)$. Their indicator vectors have disjoint supports, belong to $W$, and span a space $S$ of dimension $s$. The other slots on this side can again be summed out. The binary rank of the seed matrix is therefore at least the rank of $\beta_T$ restricted to $W\times S$, namely $s-\dim(S\cap\operatorname{rad}\beta_T)$.

Write $D=2T\ne O$. A vector in the latter intersection assigns a bit to each surviving bin, assigns zero on the exterior, and must be invariant under translation by $D$. Link bins whenever some $P$ and $P+D$ belong to them, including an exterior vertex when one of the two lies outside. Labels must agree on every connected component and vanish on the exterior component. A one-bin component independent of the exterior requires that bin to be stable under addition by $D$. For any point in such a bin the two canonical abscissas have the same fixed bits. Their field difference lies in the difference set of the $w$ free-bit words, containing at most $w^2$ possibilities. The function $x(P+D)-x(P)$ has distinct double poles at $-D$ and $O$ and thus at most four points per finite value. There are at most $4w^2$ one-bin components independent of the exterior. All other free components use at least two bins, hence the dimension of the intersection is at most $(s+4w^2)/2$. [Abscissa difference bound](../sampled_slots/block_support.md)

It follows that the seed matrix has binary rank at least $(s-4w^2)/2=\Omega(N/w)$ at the growing-arity widths $w=o(N^{1/3})$. Thus any exact XOR of products of arbitrary binary functions of the two halves requires $2^{b-b/m-O(1)}$ terms or more. This assertion allows arbitrary computation within either half. It applies to an all-seed **linear combination of separated terms**, and it does not give a lower bound on the size of a circuit with arbitrary cross-half gates.

Both dimension arguments use the entire seed family, including seeds with ambiguous, irregular or empty fibers and the particular multiplier $1$. They give no conditional error bound for a predictor on fresh guarded singletons. [Finite verification](../../checks/check_orbit_parity.py) · [Recorded output](../../evidence/check_orbit_parity.txt)
