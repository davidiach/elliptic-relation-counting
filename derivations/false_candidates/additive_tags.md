# When always keeping the witness forces an additive test to accept almost every address

Consider the prime-order subclass $G=E(\mathbb F_p)$, $N=|G|$ an odd prime, with the usual independent thinned mixed slots, and mask the target uniformly over $G\setminus\{O\}$. Let $U$ denote a regular guarded singleton. Take any finite abelian group $A$ with $N\nmid |A|$ and any local functions. On a legal address the proposed equality test is

$$
                  \sum_i f_i(\sigma_i,P_i)=g(T)\quad\text{in }A.
$$

The functions can differ between slots, see their entire local seed, and have arbitrary lookup or computation cost. Shared random coins independent of the inputs are allowed. Suppose the test always accepts the witness on $U$. At growing arity it then accepts an expected **fraction at least $1-O(1/m)$** of the legal addresses on $U$. In fact, writing $D=\prod_i|F_i|$ and $F$ for the number of wrong addresses accepted,

$$
 \mathbb E[F\mid U]\ \ge\ (1-O(1/m))\,\mathbb E D-1
                      =2^{b-o(b)}\quad\text{as a lower bound}.
$$

This covers local additive checksums even with wide tags in a group of order coprime to $N$. It applies when the displayed equality alone is the candidate test. Extra filters using joint seed or target arithmetic, and tags with an arbitrary jointly seed-dependent target value, require a different argument.

## Change to the planted law

Let each $P_i$ be uniform on the $N-1$ affine points. Independently at each slot sample its seed conditional on containing that point, and set $T=\sum_i P_i$. Each affine point has identical marginal inclusion probability $\rho$. Thus each pair consisting of a seed and a legal tuple gets mass $\Pr(\sigma)/\mathbb E D$, with $\mathbb E D=((N-1)\rho)^m$. On $U$ this is exactly the original independent masked law conditioned on $U$, because the target singles out one tuple. In particular the new measure conditions the points of the true tuple without requiring independence between competing tuples.

Under this planted law $U$ fails with probability $O(m^{-2})$. Indeed its probability is $(N-1)\Pr_{\rm mask}(U)/\mathbb E D$. The [targetwise factorial bound](../sampled_slots/isolation.md) gives $\Pr_{\rm mask}(H=X=1)\ge\lambda(1-3m/N-1/(2m^2))$, where $\lambda=(N\rho)^m/N$. Cross-slot and endpoint irregularity have probability $2^{-b+o(b)}$ before conditioning, and $\lambda^{-1}=2^{o(b)}$ in the stated regime. Finally $(N-1)\lambda/\mathbb E D=(N/(N-1))^{m-1}$. One may replace the affine point law by uniform independent $X_i\in G$, defining the missing local tag at $O$ arbitrarily and choosing an arbitrary $g(O)$. This costs at most $m/N$ in total variation. Fixing the public tape, the local tags $Z_i$ then satisfy $\sum_i Z_i=g(\sum_i X_i)$ with probability at least $1-\delta$, for $\delta=O(m^{-2})$. Each $Z_i$ is drawn from its own conditional law given $X_i$.

## Correcting the common increments

For uniform independent $h,x\in G$, let the increment of slot $i$ be the tag at $x+h$ minus the tag at $x$, sampling fresh local seeds given both point values. For distinct $i,j$, shift the two points of a random tuple by $+h,-h$ and refresh just those two seeds. Both tuples have the same sum and each has the original distribution. Consequently, except with probability $2\delta$, the increment of $i$ equals an increment of $j$ at an independent basepoint with the same shift $h$. Conditional on $h$ the two increments are independent.

Let $d(h)$ be a most likely slot-one increment. Agreement of slots one and two shows that a slot-one increment equals $d(h)$ except with average probability $2\delta$: the largest atom of its conditional distribution is at least the conditional agreement probability. Comparing each other slot with slot one gives error at most $4\delta$. For uniform $h,k,x$, take independent tags at the three positions $x,x+h,x+h+k$ for slot one, reusing the middle tag when adding increments. Each of the three increments has the preceding marginal. It follows that $d(h+k)=d(h)+d(k)$ except with probability $e\le6\delta$.

For completeness, an elementary stability argument corrects $d$. For fixed $v$ and uniform independent $y,z$, the differences $d(v+y)-d(y)$ and $d(v+y+z)-d(y+z)$ agree except with probability $2e$ by applying the approximate equation twice. Because $y,y+z$ are independent, the modal difference $h_0(v)$ has probability at least $1-2e$, for every $v$. For fixed $v,w$ the differences along the path $y,y+v,y+v+w$ agree with their respective modes, including the direct difference, with positive probability when $6e<1$. Hence $h_0(v+w)=h_0(v)+h_0(w)$ for every $v,w$. Also $d(v)=h_0(v)$ except with probability $3e$ over uniform $v$.

Finally, take two independent uniform points with fresh local tags in slot $i$. Their difference is uniform. The tags minus $h_0$ of the corresponding point agree except with probability at most $4\delta+3e\le22\delta$. Another modal-atom bound thus produces $c_i$ with $\Pr(Z_i=h_0(X_i)+c_i)\ge1-22\delta$. The group $G$ has prime order $N\nmid |A|$, so $h_0=0$. The original approximate equation and a union bound also give $g(\sum_iX_i)=\sum_i c_i$ except with probability $(22m+1)\delta$.

## Count the wrong candidates

Return to the affine planted law. Draw a second legal address uniformly on the same seed. Its joint marginal with the seed equals that of the first legal address: planting weights the seed by $D$, and both choices conditional on the seed are uniform on its Cartesian grid. In particular the alternative's local tags equal $c_i$ except with total probability $O(m\delta+m/N)$. The target value from the original tuple equals $\sum_i c_i$ with the same error order. This comparison uses no independence between the alternative and the original tuple. Conditioning on $U$, of planted probability $1-O(m^{-2})$, proves the expected accepted fraction $\mathbb E[\alpha\mid U]\ge1-O(1/m)$ for $\alpha=(1+F)/D$.

The planted measure satisfies $\mathbb E_{\rm plant}[1/D]=\Pr_{\rm seed}(D>0)/\mathbb E_{\rm seed}D$. By Cauchy–Schwarz,

$$
\mathbb E[1+F\mid U]\ \ge\
\frac{(\mathbb E[\sqrt\alpha\mid U])^2}{\mathbb E[1/D\mid U]}
\ \ge\ (1-O(1/m))\,\mathbb E_{\rm seed} D .
$$

The constants in these probability estimates are uniform over the fixed public tape. The [finite checker](../../checks/check_candidate_fingerprint.py) tests the planting identity on the mixed seed law and the corrected-increment step on small cyclic groups. On the order-19 curve it also fits **all** additive binary checks with the local table indexed by the labeled legal mask and chosen point. A true singleton yields one linear equation in all the table bits and the target bits. A wrong row lies in the linear span of these equations exactly when every complete check must accept it. The constraint rank is 231, and all 103,620 weighted wrong pairs lie in that span. Thus allowing an unrestricted number of these complete checks simultaneously prunes zero wrong addresses in this fixture. [Recorded output](../../evidence/check_candidate_fingerprint.txt)
