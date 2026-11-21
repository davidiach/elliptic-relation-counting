HILLCLIMB 007 -- BALANCED CLUSTERS, WITH INTERPOLATION FULLY CHARGED
31 October 2026

Outcome
-------
Derived exact formulas for the balanced cluster family

  X = {2^u*x+r : x in {0,1,2,4,...,2^(n-2)}, 0 <= r < 2^u}.

It has equal low-bit orientation populations, addressing a limitation
of the previous geometric example. The formulas account jointly for
local integrality, normalized mode weights, interpolation degree,
Hermite multiplicities and the readout's precision loss.

With triangular multiplicities, the standard full Hermite readout, and
unrestricted mode enumeration, neither of the two tested node policies
improves the old asymptotic enumeration scale:

1. Retain every coarse node at each active valuation level.
   Integrality forces n <= 3u+3, and leaves roughly n/3 low bits unused.
   The feasible catalogue contains
     2^(Omega(M * min(1,sqrt(n/v))))
   candidates, where v <= u is the number of active levels and M is
   the edge count. Since u <= log2 M, this is no better than the old
   M/sqrt(log M) exponent scale. Growing n makes the bound worse.

2. Retain only positive central nodes.
   Their constants are integral, but fewer nodes require larger
   multiplicities for the unchanged deformation degree. The same
   obstruction remains for n<v; for n>=v the feasible catalogue
   already contains at least 2^(M/128) candidates under the stated
   rounding convention.

These are statements about the specified method's feasible candidate
universe. They do not say that every candidate contributes. They do not
rule out another interpolation functional, coefficient-dependent
precision, support-aware pruning, selective coarse nodes, or a different
representation. No complete ECDLP exponent improved in this iteration.

What was checked
----------------
Thirty finite fixtures; 4640 direct derivative valuations; 1395 exact
mode weights; nine independent rational Gaussian eliminations; full
interpolation identities on 75 monomials; and 74 large-parameter integer
certificates. A negative control detects the coarse-node constant
valuation -16 that checking only central nodes would miss.

The asymptotic conclusions come from the proofs in cluster_bounds.tex.
The large-parameter checks evaluate closed formulas only; they do not
construct enormous packets or claim benchmark speedups. An initial
diagnostic accidentally tried to construct a large coarse grid; that
process was stopped, and these checks now use the proved closed forms.

Next gate
---------
The current full Hermite functional reconstructs a polynomial space
whose degree can exceed the degree of the actual deformation. Its
highest-jet coefficient can therefore be an unnecessarily expensive
choice of readout. This is an explicit limitation of the present bound.

The next experiment will optimize the exact linear readout on the actual
degree-limited space, using a 2-adic module calculation. Any improvement
must have an exact rational reconstruction certificate and a matching
minimality certificate for the stated model. Then its precision savings
must be propagated through mode enumeration and recursion. Lower
precision at a few finite nodes alone will not establish the goal.

Reproduce
---------
    python -B research/hillclimb_007/run_all.py

The launcher retains source hashes, result JSON, stdout and stderr in a
fresh runs/ directory. latest_run.json identifies the checkpoint.

The square-root-in-b ECDLP goal remains active and unachieved. The latest
complete conditional baseline is still 27.0138774*b/(log2 b)^(1/4).
Original archives, vendor modules and previous executable checkpoints
are unchanged.
