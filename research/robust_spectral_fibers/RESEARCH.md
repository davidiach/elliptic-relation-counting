# Spectral dispersion of robust fiber kernels and uniform target mixing

**Research follow-up, 27 November 2026.** This note extends PR #43 (head
`fd8fc3f9125b45fb3fe566369ab359d7cb830835`) rather than making a claim
about all affine-mask outputs. The new statements are (i) a uniform
single-multiplier circular-variance exclusion, (ii) a strict separation
from the earlier modal/trimmed-width exclusions, and (iii) a conditional
Fourier mixing theorem at every target, with a new unweighted randomized
rank-completion branch. The structural results are abstract finite-group
statements; the cost of producing and enumerating the arithmetic fibers
uses the inherited LZR closed counting premise. A uniform deterministic
end-to-end saving remains open. No deployed cryptographic tests are used.

## 1. Stronger good-fiber tolerance at the same conditional scale

Use the actual parity-specific affine masks and notation of #43. Write
`r` for the prime-order group size, `b` for the field-prime bit length,
`L=log_2 b`, and `m` for the block count. Let `h>=256` be an integer,
`h+k=o(m)`, with `k` the retained-set deletion allowance. All asymptotic
claims concern the sufficiently large parameter range of #43. Its seed
potential gives, for the retained sets I of size at least m-k,

    sum_I barPhi_I < 2^(-a)/256,
    3/4 <= W_I/q^|I| <= 5/4,
    a=floor((h+k+1)b/(128m)).

Replace #43's goodness threshold by

    tau = (100m)^(-h-4),   delta=8 tau.

Then delta <= (100m)^(-h-3). All of its marginal and fiber-mass
consequences continue to hold, since this is a stricter threshold. In
particular a good fiber with uniform tuple law F satisfies, for every
J of at most h distinct blocks and every pinned assignment A,

    (1-delta) prod_(i in J) 1/n_i
        <= Pr_F[X_J=A] <=
    (1+delta) prod_(i in J) 1/n_i.                         (1)

Its main tuple count is at most 2 mu_I and its normalized group mass
p_I(t)=d_I(t)/W_I is at least 1/(2r). The latter follows directly from
the empty-pin bound relative to W_I/r and tau tending to zero, not
from the looser count lower bound relative to mu_I.

The fraction of bad group values for any I is at most
`beta=2^(-a)/(256 tau^2)` by Markov. Also `chi_I^2=r sum_g(p_I(g)-1/r)^2
<=4 barPhi_I` by #43. The exponent a dominates `(h+4)log(100m)`, so
beta and chi_I^2 tend to zero exponentially on the b/m scale. For a
batch of at most 2^(a/4) targets, the adaptive doubled-translate potential
of #43 has expectation at most `2*2^(a/4) barPhi_I < tau^2` eventually.
The same exact nonnegative counts select a translate with all main
fibers good at the stronger tolerance. No new moment-copy cost is
introduced: changing the rational threshold only changes the numerical
guard. For a fixed h,k, the conditional cost is still

    log_2 T = O((h+k+1)b/m + mb/sqrt(L)),

which gives the quarter-power component at m asymptotic to L^(1/4).
The previous fixed-translation-list version can likewise use the
stronger tau because log(1/tau)=o(a), without changing its resource
envelope. For growing h this is the actual h-dependent envelope of
#43, not a fixed-h complexity estimate applied outside its range.

## 2. A uniform circular-variance exclusion

This lemma does not assume any elliptic-curve or mask structure. Let U
be the independent uniform tuple law on m nonempty finite blocks. Let
F be ANY law on their Cartesian product satisfying (1), where
256<=h<=m and delta<=(100m)^(-h-3). Let u_i(A) be labels in F_r, r an
odd prime. Assume sum_i u_i(X_i) is constant modulo r on the support of
F. For a block label distribution define V_i as the least possible
variance of the normalized real lifts u_i(A)/r in an interval of
length one (the block is uniformly distributed). The minimum is over
all circular cuts; offsets do not matter. Set V=sum_i V_i.

**Theorem 1.** If some block label is nonconstant, then

    V > h/60000.                                        (2)

This holds for every nonzero scaled version of any surviving kernel
direction, and for every nonconstant linear combination of such
directions, including any available true affine direction. It is a
single-multiplier conclusion. It does not require small modal deficit,
a bound on the sum of individual widths, or a favorable multiplier
distribution.

**Proof.** Choose minimum-variance cuts (there are only finitely many
different lifts up to common shifts in each block). Let X_i,X_i' be
two independent samples of the real lifts under U, and put
D_i=X_i-X_i', D=sum_i D_i. They are independent, symmetric, satisfy
|D_i|<=1, and E_U D_i^2=2V_i. Use the same differences under two
independent samples of F, denoting that pair law by F^2. By the
constant-residue assumption, D is always an integer under F^2.

Write M_i=E_U exp(2 pi i X_i). Choose a circular cut centered at the
argument of M_i (any center if M_i=0). For lifts relative to that
center, x in [-1/2,1/2], the elementary inequality
`1-cos(2 pi x)>=8x^2` gives mean square at most `(1-|M_i|)/8`.
The minimum variance V_i is no larger than that mean square. Hence

    1-E_U cos(2 pi D) = 1-prod_i |M_i|^2
      >= 1-exp(-16V) > (3/4) min(V,1)     if V>0.       (3)

For the final inequality, concavity on [0,1] and monotonicity above
one suffice. Nonconstant labels guarantee V>0. We will obtain a
contradiction if V<=h/60000.

Put H=4 floor(h/4), and define the real polynomial

    P_H(x)=sum_(j=1..H/2) (-1)^(j+1)(2 pi x)^(2j)/(2j)!.

For every real x, `P_H(x)<=1-cos(2 pi x)` with error at most
`|2 pi x|^(H+2)/(H+2)!`. These global alternating Taylor bounds follow
by integrating twice the alternating-sign error for the preceding
cosine truncation, starting with `1-cos>=0`; they are not an assumption
that |x| is small. In particular P_H(D)<=0 on F^2.

For the independent pair law, symmetry and |D_i|<=1 give

    E_U exp(8D_i) = E_U cosh(8D_i)
      <= 1+2(cosh(8)-1)V_i <= exp(4000 V_i).

Here `2(cosh(8)-1)<4000`, for example by e<11/4. Multiplying the
independent moment generating functions, the Taylor remainder has
expectation at most

    (pi/4)^(H+2) E_U[cosh(8D)-1]
      <= 4000 V exp(4000V) (4/5)^(H+2)
      <= 5000 V exp(-2h/15)
      < (1/4) min(V,1).                               (4)

In the penultimate inequality use H+2>=h-1, log(4/5)<=-1/5, and
4000V<=h/15. For the last, V<=h/60000<=h and
`5000 h exp(-2h/15)<1/4` for h>=256: it decreases there and is already
below 1/4 at 256 (even the weaker bound using 2^floor(2h/15) in the
denominator suffices). All constants are deliberately conservative.

It remains to transfer the polynomial expectation to the paired fiber.
For any at most h block indices, the paired fiber marginal has density
within 3 delta of the paired independent marginal, multiplicatively,
by (1) and independent sampling of the two fiber tuples. Each monomial
of D^(2j) uses at most 2j<=h block indices. The multinomial expansion
and absolute values therefore give

    |E_(F^2) D^(2j) - E_U D^(2j)|
      <= 3 delta E_U (sum_i |D_i|)^(2j)
      <= 6 delta m^(2j-1) V.

The second inequality uses sum |D_i|<=m and Cauchy-Schwarz for its
square. Consequently, using 2 pi<7, H<=h<=m, V<=m, and the stated delta,

    |E_(F^2) P_H(D) - E_U P_H(D)|
      <= 6 delta m^2 (7m)^h min(V,1)
      < (1/4) min(V,1).                               (5)

Equations (3)--(5) give positive E_(F^2) P_H(D), contradicting its
pointwise nonpositivity. This proves (2). QED.

The strengthened tolerance is essential to this proof: summing the
absolute monomial errors is not justified using only a constant
multiplicative h-wise error. Conversely the required tolerance costs
nothing new at the asymptotic arithmetic scale of Section 1.

### 2.1 Spectral consequence for the entire residual quotient

From the phase-centered cut inequality in the proof,

    prod_i |E_A exp(2 pi i u_i(A)/r)|
       <= exp(-8V) < rho := exp(-h/7500).              (6)

Use reduced incidence coordinates, with one reference per block. Let
K be the kernel of the complete zero-value contrast stack containing
one good fiber, so that all u in K are constant in score on that
fiber. The true log-difference direction lambda belongs to K. Every
nonzero reduced direction is nonconstant on some block. Equation (6)
is simultaneous for all nonzero elements of K. If one chooses a basis
of K including lambda, the distribution of its score vector on a
uniform independent short tuple is rho-biased in every nontrivial
additive character of the corresponding quotient. This is a
consequence of the structural lemma; no Fourier evaluation or true
logarithm lookup is needed to obtain the good fiber.

For deletion or pinning of at most k blocks the sum of the minimum
variances falls by at most k/4, since an individual lift interval has
length one. If k<=h/30000, the remaining variance is >h/120000, so the
same product bias bound with rho'=exp(-h/15000) applies to every
nonzero original kernel direction on every permitted retained subset.
Fixed pins contribute only a constant phase. This uses #43's deletion
energy guarantees separately for the group marginal.

For comparison with the older Hadamard bound, a simple pigeonhole
corollary is `r^dim(K)<=B^D`, where D is the reduced coordinate
count and B=ceil(sqrt(60000 m/h)). Indeed two kernel vectors whose
canonical normalized coordinates lie in the same B-adic boxes have
a difference with lifts of width <=2/B and total variance <=m/B^2;
(2) forbids a nonzero such difference. This is an angular packing
bound, not a stronger rank fraction than the separate exponential
quotient-compression result of #18.

## 3. Strict separation from modal and trimmed-width exclusions

The new criterion is not just a restatement of #43's earlier
criteria. The following construction works for candidate labels on
ANY fixed sufficiently large blocks obeying #43's size bounds, with
any specified group scalar labels lambda. It is a separation of
necessary tests, not a claim that the constructed direction belongs
to the actual kernel (Theorem 1 proves it cannot).

Let `ell=floor(r/(10 sqrt(m)))`, assuming the asymptotic regime of
#43, so ell tends to infinity. Independently assign each indexed
point a label u uniformly in [0,ell-1]. Subtract the reference label
in each block if reduced coordinates are required. Such translations
have no effect on the tests. For any pair (a,c) in F_r^2 with c!=0,
any block, and any cyclic interval J of cardinality <=floor(ell/2),

the independent labels a lambda(A)+c u(A) hit J with probability at
most 1/2 individually. This is because the image of the original
ell-point interval under a nonzero dilation and translation is
injective, and J has at most ell/2 points. Hoeffding therefore bounds
the probability that at least 3/4 of a block lies in J by exp(-n_i/8).
Union over m blocks, fewer than r^2 pairs and r starting points (the
longest interval suffices) gives probability at most
`m r^3 exp(-n_min/8)=o(1)`, since n_min grows as 2^Theta(b/m).

Choose labels satisfying this event. Every nonzero combination with
c!=0 has modal deficit at least m/4. In a trim with p pinned blocks
and sum of excluded fractions E, p+E<=h/16, at most h/4 blocks are
pinned or lose more than 1/4 of their labels. In every other block,
no retained set of at least 3/4 of its labels fits in an interval of
cardinality floor(ell/2). Its integer lift range is at least
floor(ell/2). Thus the sum of retained ranges is at least

    (m-h/4) floor(ell/2) > r

for sufficiently large m with h=o(m). This excludes all earlier
low-loss trimmed-width tests on ALL scalar multiples and affine
combinations with c!=0, not just a particular cut. Combinations with
c=0 are true affine directions; they cannot be used as nonaffine
obstructions. Yet the unscaled u has

    V_min(u) <= m (ell-1)^2/(4r^2) <= 1/400 < h/60000

for h>=256, by the elementary bounded-range variance inequality.
Theorem 1 excludes it. This is a genuine new uniform restriction on
high-modal-deficit, wide-lift directions. No random choice of the
actual arithmetic seed is needed for the theorem; randomness here is
only an existence proof separating the numerical exclusion tests.

## 4. All-target conditional spectral mixing

We can turn the angular theorem into a stronger quantitative
statement about the actual remaining directions. Let C_F be the matrix
of all within-fiber contrasts of a good main fiber F. As in Section 2.1,
require span(lambda) to be contained in K and K to be contained in
ker(C_F). Thus every score in K is constant on F. Pick any complement
basis of lambda in K, and use its scores on a uniform independent short
tuple as a vector S in F_r^d (d=dim(K)-1). Write G for the scalar group value,
and mu for the joint law of (G,S). If d=0 there is no nonaffine
residual direction, and the following statement is vacuous. The
basis and scalar group coordinates are for the proof, not an oracle
for their hidden values. Let p be the marginal of G, and suppose
`chi^2=r sum_g(p(g)-1/r)^2<=epsilon<1`.

For every nonzero character index v of the score space and every
scalar character index a of G, the weighted labels are a lambda plus
a nontrivial complementary combination. They are a nonconstant
kernel direction, so (6) yields `|mu_hat(a,v)|<=rho`. Parseval in the
first coordinate gives

    sum_a |mu_hat(a,v)|^2
      = r sum_g |sum_s mu(g,s) exp(2 pi i <v,s>/r)|^2
      <= r sum_g p(g)^2 <= 1+epsilon.                 (7)

For ell>=2 copies, Cauchy-Schwarz (or simply (7) for identical
copies) and the maximum bound on the remaining factors imply

    r |sum_s mu^{*ell}(t,s) exp(2 pi i <v,s>/r)|
       <= (1+epsilon) rho^(ell-2).

Also `r p^{*ell}(t)` lies between `1-epsilon` and `1+epsilon`, since
`sum_(a!=0)|p_hat(a)|^2<=epsilon`, and further Fourier factors have
absolute value at most one. Hence we obtain the following theorem.

**Theorem 2 (every prescribed target).** The score law of ell
independent short tuples conditional on their group total being any
specified t has all its nontrivial character biases at most

    eta_ell = [(1+epsilon)/(1-epsilon)] rho^(ell-2).    (8)

In particular the probability of ANY specified proper affine
hyperplane of the residual score space is at most `1/r+eta_ell`, by
the character expansion of its indicator. The same upper bound
applies to every proper affine subspace, since it is contained in a
hyperplane. This holds uniformly in t and every nonaffine kernel
direction. For ell=3 and h growing, the failure mass is exponentially
small in h (with the explicit conservative constant above). Earlier
complete triple-fiber rank or constant-gap results do not by themselves
give this every-target bias estimate.

The proof also works with different permitted retained sets or fixed
pins in each copy when the corresponding marginal energy bounds are
at most epsilon: use Cauchy-Schwarz for the first two factors in (7)
and the common maximum rho' of Section 2.1 for the rest. In the group
marginal calculation use the same Cauchy-Schwarz on nontrivial
frequencies. Fixed pins only shift the group total and score by
constants. This provides a quantitative pin-robust version, not a
guarantee under arbitrarily many pins.

## 5. A count-only randomized rank increment, and the remaining gap

This is a constructive consequence of the new analytical theorem,
but it is NOT a new uniform deterministic exponent. Keep the same
arithmetic seed and the initial good-fiber contrast stack. Let the
current zero-value row space W contain that stack, and suppose its
kernel has some nonaffine direction. Fix any target t, and a genuine
reference ell-copy tuple of value t. For ell>=3 a reference exists:
fix ell-3 short tuples and use #43's adaptive good translated fiber
for the remaining total. No extra unknown labels are needed.

Draw T as ell-1 independent uniform short tuples. If the main fiber
at t-val(T) is bad, reject it. Otherwise enumerate that whole good
fiber (the adaptive algorithm is not used to change T during this
random test), and compare the resulting ell-copy incidence vectors
with the reference. Each difference has zero group value, and
ordinary modular elimination tests whether it lies in W. Append any
independent differences. If u is a fixed nonaffine member of the
current kernel, failure of this test on a good T implies that every
main completion has the reference score in direction u.

Write beta for the uniform fraction of bad main values. The law of
val(T) has each point mass <=(1+epsilon)/r by the marginal convolution
bound. Thus bad T has probability <=(1+epsilon) beta. A T which is
good has main count probability >=1/(2r). Relative to the uniform
T law, its weight in a truly uniform ell-copy target fiber is
`p(t-val T)/p^{*ell}(t)>=1/(2(1+epsilon))`. Therefore the failure
probability of the rank-increment test is bounded by

    (1+epsilon) beta
       + 2(1+epsilon)(1/r+eta_ell).                   (9)

This comparison uses counts (not distinct projected support sizes).
It is valid for every current nonaffine direction and every fixed
target, including after earlier successful rank increases. It also
holds with the permitted pinning variant if the main fiber is good
and its mass bound holds. It does not require a candidate-dependent
weight in the tensor contraction or a high-arity lookup.

For example fix h>=256, ell=256, and take sufficiently large original
parameters with epsilon<=1/4, beta<=1/16 and r>=31. The crude
bound `rho^(ell-2)<=2^-8` follows from h(ell-2)/7500>8 and e>2. Then
(9) is strictly below 1/2. With growing h, ell=3 suffices eventually
and (9) is O(beta+1/r+exp(-h/7500)). No practical efficiency claim is
made about these constants. Repeating independent tests yields
expected O(D) attempts for the abstract zero-value rank-completion
stage: each success adds an independent genuine relation, and the
true affine direction gives the rank ceiling D-1. At this point the
zero-value row rank is exactly D-1. Sampling the explicit point indices
costs O(b) expected random bits per attempt for fixed ell; the
ordinary goodness counts and complete good-fiber enumeration have
the conditional #43 budget, with constants depending on ell,h,k.
The total expected random-bit bound here is O(Db), NOT the O(b)
expected total of the earlier weighted polynomial trial. The value
of this branch is its uniform all-target guarantee and its use of
only ordinary exact unweighted counts plus explicit arithmetic. It
is not a better total-randomness or bit-time exponent claim.

A nonuniform fixed-list consequence is also immediate. For an
initial quotient of dimension d+1 including the true slope and a
uniform upper bound q0<1 on (9), choose s independently sampled
translation tuples. For any fixed nonaffine direction and any target,
the probability all s fail to produce a contrast nonzero on that
direction is at most q0^s (use a fixed genuine reference at each
target). A union bound over at most r^(d+1) directions and r targets
shows that a single list of length O((d+2) log r / log(1/q0)) exists
which works for all targets and directions. This is a size/existence
claim. A uniform inexpensive constructor is not obtained by the
union bound; it is not a deterministic guarantee. Density of the
successful translations, even exponentially high density in h,
does not itself give a small deterministic hitting list. The
remaining problem is to select them within budget without assuming
an arbitrary weighted score predicate can be contracted cheaply.

## 6. Verification, novelty, and scope

Run `python -B check_spectral_fibers.py` (standard Python, assertions
enabled). The verifier performs exact integer/rational controls of
circular-variance cuts, the strengthened tolerance and proof constants,
wide-lift separation, h-wise modular conditioning and moments, and
finite Fourier mixing / random-translation failure bounds. Certified
rational trigonometric intervals are used for small-field Fourier
checks, not floating-point acceptance. The verifier rejects optimized
Python, checks all 1,010 distinct projective classes, and compares the
paired moments through degree 256 with direct integer powers. The
perturbation fixture checks the exact changes from adding a small atom
to the uniform law. An additional three-copy example over F_11 has two
positive failure counts that can be enumerated by hand. The transcript
records the bad-fiber and constant-score counts separately.
Synthetic modular arrays need not be full arithmetic
fibers; the high-degree modular-conditioning controls deliberately
illustrate that exact marginal uniformity alone is not full modular
rank. None implements the asymptotic contraction oracle or deployed
cryptographic systems.

The Taylor-minorant technique and conditional marginal bounds have
ancestors in the repository's PR #4 (`RESILIENT_SEEDS_AND_VARIANCE.md`)
and PR #43. The novelty proposed here is the uniform single-multiplier,
relative-error variance exclusion and the subsequent all-target
Fourier estimate. Limited independence and product-test estimates
are classical subjects; related background is *Bounded Independence
Plus Noise Fools Products*, CCC 2017,
https://doi.org/10.4230/LIPIcs.CCC.2017.14 , and Gopalan, Kane, Meka's
work on pseudorandom generators for Fourier shapes. None of those
results is assumed as a derandomizer of the conditioned translation
predicate. No independent priority review or proof of the upstream
closed contraction premise has been completed. The main deterministic
selection lemma remains missing; the full mandatory breakthrough
checklist does not pass.
