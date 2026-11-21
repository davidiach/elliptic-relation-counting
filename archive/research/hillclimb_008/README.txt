HILLCLIMB 008 -- DEGREE-AWARE LINEAR READOUT
31 October 2026

Outcome
-------
The standard full Hermite functional can waste precision when its jet
set reconstructs a larger polynomial space than the actual degree bound.
This iteration computes and certifies the optimum exact linear readout
on the smaller degree-limited space.

The savings are real. Examples for consecutive coordinates:

  packet L   levels u   mass M   degree   full loss   optimal loss
     4           2         4      28         45            27
     8           3         8     120        128           102

These numbers are denominator-loss exponents in bits, not complete
working precisions or algorithm runtimes. In these two cases they are
also the optimum scale for budgets floor(lambda / 2^v2(tau)). The masses
are compatible with full packets, but these are polynomial readout
experiments, not complete C38 contractions. Smaller mass fixtures are
explicitly marked as polynomial tests.

The saving does not change the asymptotic precision scale of this model.
An explicit integer polynomial proves that ANY exact linear readout on
the full degree-limited space, with consecutive centers and the stated
geometric precision profile, requires

  lambda >= M*L/(8*u), for M >= 4*u.

The standard readout supplies O(M*L/u), so optimizing the functional
alone cannot change that order. The lower witness remains valid when
more jets are supplied or some centers are omitted.

Analogous lower bounds for the balanced-cluster grids of iteration 007
are M*n^2*B/(32*v) for all coarse nodes and M*n^2*B/(16*v) for positive
central nodes, under the conditions stated in readout.tex. These are
operator-precision bounds. They do not say every node must be computed,
every feasible mode contributes, or every possible precision profile
has the same contraction cost.

What was constructed
--------------------
optimal_readout.py builds the exact integer jet matrix and solves its
dyadic row-module membership problem. Each optimum has two certificates:

* Rational weights that exactly reconstruct the constant coefficient.
* An integer polynomial proving that a smaller denominator cannot work.

The solver also handles unequal node budgets and optimizes the geometric
profile. A proved inverse-minor coefficient bound controls its finite
working precision. Fraction-free elimination lifts modular answers to
exact rational identities. There is no guessed truncation precision.

verify_certificate.py independently reconstructs jet coefficients by
polynomial multiplication and directly checks both certificates, without
using the solver's matrix constructor or module diagonalization.

Validation
----------
The completed run passed:
  45 common-budget optimality certificates
  19 geometric-profile optimality certificates
  37 direct confluent-minor determinant checks
  656 inverse-minor entry checks
  45 augmented-module index checks
  4 controls where the unique full Hermite weights must be recovered
  3 corrupted-certificate negative controls
  7278 expanded polynomial-jet checks
  90 large-parameter formula certificates

The large-parameter checks evaluate closed formulas only. They do not
construct exponentially large grids or execute a contractor.

An initial audit was deliberately interrupted when a conservative
determinant precision and rational elimination made the degree-120
certificate slow. Its partial output and interrupted status are retained.
After proving a tighter inverse-entry precision bound and replacing the
lift with fraction-free elimination, the entire suite completed in about
46 seconds. The latest_run.json checkpoint identifies the successful run.

Scope and next gate
-------------------
This model requires an exact linear identity on ALL rational polynomials
of the stated degree. Actual deformation polynomials may occupy a smaller
space or satisfy stronger coefficient constraints. Computing only a
constant modulo two on promised inputs is also a different problem.
Those possibilities are not ruled out by this result.

The next experiment will construct the actual polynomial space of small
complete packets from their tensor definitions, then test whether its
structure changes the readout optimum. It must not infer a universal
factor merely from the W entries: selector polynomials have constant
terms as well. Any useful restriction must be proved for the complete
deformation and extended beyond finite fixtures before it is used in a
runtime claim. See next_experiment.txt.

Reproduce
---------
    python -B research/hillclimb_008/run_all.py

The launcher retains source hashes, exact certificates, stdout, stderr,
and an execution summary in a fresh runs/ directory. Assertions are
mandatory. The successful source hashes were stable throughout the run.

The ultimate 2^(tilde O(sqrt(b))) general-prime-field ECDLP goal remains
active and unachieved. The complete conditional baseline remains
27.0138774*b/(log2 b)^(1/4). No supplied archive, PDF, vendor module, or
previous executable checkpoint was modified.
