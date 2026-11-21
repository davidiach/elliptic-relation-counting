HILLCLIMB 005 -- BALANCED PACKET MIXING AND COEFFICIENT EXTRACTION
31 October 2026

Outcome
-------
Built two exactly balanced finite networks containing an actual restricted
C38 first tensor, with explicit pinned padding. Their selected packets
are fully occupied. This closes the singleton-packet gap in our earlier
finite compression experiments, but not the asymptotic scheduling gap.

The short target generating formula survives:

  product_pairs (1-u_a-u_b+2*u_a*u_b) * product_pinned_rows (1-u_x),

where each u_x is a linear combination of incoming modes. Its factors
now share modes. Repeated powers of a mode matter: restricting every
mode exponent to zero or one changes the full join, even with all
exterior source moments included.

The first fixture has 648 edges, 81 full packets and 70 checked target
entries. The second has 324 edges and a selected packet whose two label
pairs are both unequal; it has 330 checked target entries. Both complete
node joins match independent 64-term sums of original assignments. The
negative-control join differs at jet degree two, with 2-adic valuations
10 and 11 respectively. There are 6400 source-assignment evaluations and
6900 checks that padding packets contribute one.

Two specific routes rejected
----------------------------
1. Bounded mode scopes inherited from the small original C38 factors.
   One row of a size-L packet depends on all L/2 incoming modes. Their
   exact leading valuations are a+1+4k-2*s_2(k), L=2^a. All survive with
   precision 2L+1. This was proved and checked for L=8 through 128.
   It does not imply a lower bound against auxiliary variables, sparse
   polynomial representations, or all arithmetic circuits.

2. Apply the Cauchy permanent determinant identity after normalization.
   A concrete normalized 3-by-3 matrix violates the identity by
   -8/495990495. The actual finite packet entry equals its appropriately
   scaled permanent, and the identity defect survives at jet degree six
   with valuation 49 modulo 2^64. A raw Cauchy matrix is a passing positive
   control. Different algorithms for structured permanents remain open.

Keeping the raw Cauchy basis instead does not preserve the old enumeration
bound: Theta(M) modes then have only the same common guaranteed weight,
and all their binary exponent vectors fit the asymptotic budget. That
explicit candidate universe is exponential in M. This does not prove
all candidates contribute on C38 or exclude an implicit joining method.

Why this matters for the goal
----------------------------
We cannot infer cheap coefficient extraction or joining from a small
arithmetic circuit, or carry singleton separation through real mixing.
The next representation must handle the complete source-weighted
coefficient pairing. These counterexamples narrow that search; they do
not improve the complete ECDLP runtime or invalidate the conditional
27.0138774 baseline from iteration 004. The goal is active and unachieved.

The next concrete spectrum experiment will test nonuniform dyadic packet
coordinates. They can increase diagonal mode valuations, but may also
increase the difference-grid size and destroy local integrality. Those
costs must be measured together before accepting any apparent entropy
gain. next_experiment.txt records the proposed falsification test.

Reproduce
---------
    python -B research/hillclimb_005/run_all.py

Standard library only. The launcher refuses disabled assertions and saves
stdout, stderr, source hashes and exact result records in a new runs/
directory. latest_run.json identifies the checkpoint.

balanced_projection.py     Actual balanced C38 fixtures and full node joins.
normalized_permanent.py    Rational and finite-ring determinant counterexample;
                          dense-scope valuation checks.
closure.tex               Proofs, explicit limitations and basis tradeoff.
frontier.json             Accepted evidence and next research gate.

The diagnostic geometry is finite and has three groups. Padding extends
the construction to nine groups algebraically, but a complete asymptotic
schedule and the full unpinned C38 network are not executed here. All
original archives, published-source modules and earlier checkpoints are
unchanged, apart from the explicit frontier continuation link.
