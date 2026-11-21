# False addresses across a cut and the cost of quotient digits

Write $U$ for the fixed-target event of one ordinary relation passing the addition guards, with probability $q_U$. Legal addresses decode to ordered signed point tuples. A candidate test is one-sided if it accepts the witness on every $U$ instance. Let $F$ count distinct accepted wrong legal addresses on $U$, irrespective of how many auxiliary assignments certify each address.

## Seed-local interfaces with false candidates

Hide slot $i$, whose marginal inclusion probability for a prescribed point is at most $\rho_i$. Use the [seed-local cut model](../sampled_slots/seed_cuts.md). Before reading its seed, the rest of the graph has at most $2^a$ possible descriptions and $2^c$ possible boundary strings per description. Outside assignments determine the other signed addresses. Feasible inside and outside assignments with a common description and boundary can be spliced, and the resulting addresses are legal. The test can accept arbitrary false addresses on negative instances. On $U$ we have

$$
q_U\,\mathbb E[(1+F)^{-1}\mid U]\le 2^{a+c}\rho_i,
\qquad
\mathbb E[F\mid U]\ge\frac{q_U}{2^{a+c}\rho_i}-1.
$$

Freeze the other seeds and public test coins. For each description and boundary, order its distinct feasible outside addresses by an independent uniform random permutation and take the first. The selected outside address forces one candidate point in the hidden slot for the prescribed target. These at most $2^{a+c}$ points are chosen independently of the hidden seed. Their probability of hitting its slot is at most $2^{a+c}\rho_i$. On a singleton choose a witness certificate independently of the auxiliary permutations. At its description and boundary at most $1+F$ outside addresses can occur, since splicing each one to that certificate's inside assignment gives a distinct accepted address. The chance the selected outside address is the witness's outside address is at least $1/(1+F)$. Averaging proves the first inequality. Convexity of $1/(1+x)$ gives the second.

Thus mean false count below $1/3$ implies $2^{a+c}\rho_i>3q_U/4$. On the thinned mixed sampler $\log_2(q_U/\rho_i)=b-b/m-o(b)$. For disjoint slot regions, summing the cut widths counts each edge at most twice. If all outside descriptions have $a=o(b)$ bits of variation, a one-sided test meeting the mean bound still requires $(1/2-o(1))mb$ total edge bits. Public joint use of the hash seeds can violate the premise on descriptions. The argument also requires legal splicing; an oracle that discards nonlegal outputs globally must include that check in its network model.

## When a short residue test works exactly

Consider a quadratic integer residual $f$ on a fixed number $t$ of candidate coordinates with public centers and at most $w$ possible values each. Suppose $|f|<Cp^2$ and the quotients of divisible values are known to lie in an interval $J$ of $H$ consecutive integers with $|z|=O(p)$ throughout. For example, for $f=xy-A$ on two contiguous abscissa windows inside $[0,p)$, both of width $w$, the variation of $xy$ is less than $2pw$, hence $H\le2w+1$. This example applies to the local product and does not assume contiguous support for internal group sums.

Choose a prime $\ell>H$ different from $p$. Test whether $f\equiv pz\bmod\ell$ for some $z\in J$. This always keeps a true zero modulo $p$. If $p\nmid f$, every $f-pz$ is a nonzero integer of magnitude $O(p^2)$, with $O(b)$ distinct prime factors. Over all wrong candidate inputs and possible $z$, at most $O(bHw^t)$ primes are bad. One can find a good prime with $O(\log b+\log H+t\log w)$ bits by testing successive primes above $H$; the elementary bound on the $s$th prime is $O(s\log(2s))$. Certify a proposed prime by comparing the two predicates on all local input choices. For fixed $t$ and $H=O(w)$ this construction and its validation take $2^{O(b/m)}\operatorname{poly}(b,m)$ time at the project widths. With that prime, the test accepts zero wrong *local* inputs. Both the evaluation of $f$ and the interval comparison can use the short modulus, precomputed reduced constants, and short offsets. The interval check compares the residue of $p^{-1}f-\min J$ with $H$.

Conversely, if a residual's only quotient constraint is one interval of $H$ integers and a congruence of coprime modulus $L\le H$, every integer residual satisfies the test. Multiplication by $p$ permutes the residue classes, and the interval contains all of them. Multiple congruences sharing the same quotient have the same obstruction at their least common multiple. Thus an unrestricted $H=\Theta(p)$ quotient in this construction requires a modulus of $\Omega(b)$ bits even to exclude one input.

## Quotients of products with a prefix abscissa

A product in the alternating chord equations is $uv$, where $u=x(S_I)$ is the canonical abscissa of a running sum and $v=x(P_j)$ is the canonical abscissa of the next leaf. Could one keep, on actual successful seeds, a short public list of its possible integer quotients $\lfloor uv/p\rfloor$? The quotient itself need not be materialized by every verifier. The following result concerns tests that restrict it this way.

Take a prefix $I$ of $k$ slots followed by $j$, and leave at least one further independent slot. A list of up to $H$ possible quotient values may depend on the target, all seeds in $I$, the entire seed of slot $j$, and the actual selected $P_j$. It can also use independent coins. Use $\nu_i=N\rho_i$, $\nu_I=\prod_{i\in I}\nu_i$ and $\lambda=\prod_i\nu_i/N$, with $N$ the curve group order. Write $\beta_I=C^k k 2^s\nu_I^2/N$ for the [off-diagonal collision bound](../sampled_slots/prefix_dictionary.md). The probability of a guarded singleton whose true quotient belongs to its selected list is at most

$$
\frac{\lambda}{\nu_I}\bigl(O(H(1+\log p))+\beta_I\bigr).
$$

For $v>0$ any one quotient restricts integer $u$ to at most $\lceil p/v\rceil$ possible abscissas, with at most two group points per abscissa. For $v=0$ allow all $N$ points. Over all possible $P_j$ the sum of these bounds is $O(Hp(1+\log p)+N)$, because at most two leaf points have any particular $v$. Consequently the expected sum of the list sizes *in prefix group points* over sampled $P_j$ is $O(\nu_j H(1+\log p))$ by the uniform point marginal and $N\asymp p$. This bound holds even with list choice depending on all allowed seeds.

If $h(z)$ is the prefix multiplicity at a group point, the number of prefix tuples in any set $D$ is at most $|D|+\sum_z h(z)(h(z)-1)$. The second term has mean at most $\beta_I$ and is independent of slot $j$. Since $\mathbb E|F_j|\le\nu_j$, the expected number of prefix/leaf choices surviving the list is at most $\nu_j[O(H(1+\log p))+\beta_I]$. For each such choice the remaining independent slots give at most $\lambda/(\nu_I\nu_j)$ expected completions at the prescribed target, by forcing their last member. This proves the bound. A guarded prefix has finite $S_I$; allowing extra states in this upper estimate causes no loss.

At balanced $k=m/2+O(1)$, $\nu_I=2^{b/2+o(b)}$, $\beta_I\le2^{o(b)}$, and $q_U\asymp\lambda$. A list of $2^{o(b)}$ quotient values therefore covers the true prefix quotient on at most $2^{-b/2+o(b)}$ of the singleton seeds. Even $2^{(1/2-\epsilon)b}$ values have vanishing conditional coverage for every fixed $\epsilon>0$. Up to $a$ bits of outside advice replace $H$ by $2^aH$ by taking the union over advice strings. Full access to the remaining seeds, jointly encoded quotients, and predicates evaluating divisibility without exposing such a quotient lie outside this result.

Both limitations apply to mathematical representations of the address test. Neither gives an unrestricted lower bound on evaluating its accepted count. [Finite checker](../../checks/check_candidate_arithmetic.py) · [Recorded output](../../evidence/check_candidate_arithmetic.txt)
