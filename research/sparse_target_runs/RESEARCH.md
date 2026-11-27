# Four and five consecutive targets with excess pinned nullity

## Four targets at order 223

On $y^2=x^3+x+1$ over $\mathbb F_{211}$, the rational group has prime order
$r=223$. For the generator $P=(204,101)$, a sparse affine seed with $m=5$
has

$$
\operatorname{nullity}_{\mathbb F_{223}} M_{kP}=6
\qquad (k=0,1,2,3).
$$

Here $M_t$ has one column per tagged block element. An ordered tuple uses
two slots per block. Row $(i,x)$ sums the incidence vectors of the tuples
with total $t$ whose first slot in block $i$ equals $x$.
The result concerns accepting one of the four matrices individually. Their
stack has nullity $4=m-1$, and the next individual matrix $M_{4P}$ has
nullity $5=m$.

## Seed and potential

The reduced-floor experiment uses the least $q$ with $q^m\ge16r$, then
$T=\lceil pq/(r-1)\rceil=6$. The two pairs for each block in the fixture
give the affine maps for even and odd $y$, respectively. A finite point is
selected when the residue of $\alpha x+\beta\pmod p$ is less than $T$.
The five blocks have sizes $2,8,8,8,10$. [Exact affine seed](witnesses.json)

Let $h(s)$ count ordered half-tuples with one element from each block and
sum $s$. All 223 residues occur. For $W=\sum_s h(s)=10{,}240$,
$C=\sum_s h(s)^2=473{,}400$, and $\mu=((r-1)T/p)^m$,

$$
\Phi=\frac{r C+64(W-\mu)^2}{\mu^2}
=\frac{2372513226408939176812671967937}{2197618292246305904679118307328}
<\frac{64}{57}.
$$

Thus the potential bound alone allows a run of four failures in the
reduced-floor model. The question for the canonical bit-fixed seed with
the original parameter floor remains open.

## Why the failures occur in a progression

Write group elements as multiples of $B=(0,1)$. Start with a seed whose
blocks are invariant under negation. Negating every slot permutes the
pin rows and the tagged columns and takes target $u$ to $-u$. The
source seed had excess nullity at $dB$ and $3dB$ for $d=60$, so it
also had excess nullity at their negatives.

Translate just the two-point first block by $90B=(3d/2)B$. This translates
the half-sum histogram by $90B$, preserving $C,W,\Phi$. Each full
tuple uses that block twice, so its target moves by $180B=3dB$.
Relabeling that block in the pins and columns preserves the rank.
Consequently the four targets with coefficients $-3d,-d,d,3d$ become
$0,P,2P,3P$ for $P=2dB=120B$. The source first block used the
map $(97,24)$ for both parities. The verifier reconstructs the translated
block from the two affine maps in the fixture and checks exact equality
of the point sets. This last check matters because the affine-mask family
need not realize an arbitrary group translation.

## Exact verification

| Target | Coefficient of $B$ | Ordered tuples | Incidence rank | Pinned rank |
|---|---:|---:|---:|---:|
| $0$ | 0 | 469,952 | 31 | 30 |
| $P$ | 120 | 469,728 | 31 | 30 |
| $2P$ | 17 | 469,728 | 31 | 30 |
| $3P$ | 137 | 469,952 | 31 | 30 |

There are 36 tagged columns. For each target $t$, let $A_t$ have one
incidence row per ordered tuple. Its kernel contains the $m-1$ block
contrasts and the log vector shifted by $t/(2m)$ in $B$-coordinates.
The verifier checks independence and exhibits enough genuine rows to reach
the matching upper bound on rank. Exchanging the two slots of each block
gives $A_t^\mathsf{T}A_t=2M_t$. The one-dimensional rank loss occurs
when relations are summed modulo 223. Every pin count is nonzero modulo
223 at all four targets.

## Five targets at order 67

On $y^2=x^3+x+13$ over $\mathbb F_{59}$, take $B=(1,29)$ and
$P=29B=(4,50)$. The reduced-floor affine seed has four blocks of sizes
$2,10,10,10$. All five targets $kP$ for $0\le k\le4$ have nullity 5
over $\mathbb F_{67}$. The remaining 62 targets have nullity 4. Stacking
the five matrices gives nullity 3, the forced block-contrast dimension.

At threshold $T=6$, the 2,000 half-tuples cover every group element.
Their collision count is $C=60{,}876$ and the exact potential is
$151752917869129250809/151182490735070429184<64/57$.
The ordered tuple counts at the five targets, in order, are 59,566,
60,419, 60,876, 60,419, and 59,566. Each raw incidence matrix has rank
28 on 32 columns. Each pinned matrix has rank 27.

For a five-term run, start with negation-symmetric blocks and check the
three source targets $0,P,2P$. A failure at each gives failures at
$-2P,-P,0,P,2P$ by symmetry. Translating the first block by $P$ shifts
every full target by $2P$ and permutes the half-sum histogram. For this
witness the first block used $(43,34)$ for both parities before translation.
The translated mask uses $(24,46)$ on even $y$ and $(22,4)$ on odd $y$.
The verifier checks equality of the translated point set and the affine mask.

### A second five-target witness

On $y^2=x^3+x+46$ over $\mathbb F_{101}$, the group has order 109 and
$P=(36,51)$. The five blocks have sizes $2,6,6,6,10$ and
$\Phi=1.0211280659\ldots$. At $0,P,2P,3P,4P$ the raw rank is 25 and
the pinned rank is 24 on 30 columns, giving nullity 6 at all five.
Every other target has nullity 5, and the five matrices stacked have
nullity 4. The fixture records the exact rational potential.

## Reproduction

Run from the repository root with NumPy and python-flint installed.

```sh
python -B -m research.sparse_target_runs.check_four_targets
```

The verifier computes all target ranks by counting the remaining slots.
At each target in a claimed run it also joins the half-tuples, counts every
pinned matrix entry directly, and checks equality over the integers before
reducing modulo $r$. [Verifier](check_four_targets.py) · [Recorded output](verification.txt)
