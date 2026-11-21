HILLCLIMB 001 -- CONDITIONAL COMPOSITION AND EXECUTABLE INTERFACES
31 October 2026

Outcome
-------
The first iteration combines C38 alternating-anchor arithmetic, root-fixed
joint simplex words, and the sharper adaptive nine-group contraction rule.
It also replaces conservative matrix parameters with stronger published
ones after checking the support inequalities again.

Claim status: conditional derivation under the supplied foundational
contraction, encoding, uniformity, and ECDLP-reduction proofs. Those
foundations have NOT been independently certified in this iteration.

On the common general-prime, supplied-order, polynomial-coprime-cofactor
class, the derived bounds have scale b/(log2 b)^(1/4):

                                      Time coefficient   Space coefficient
Supplied simplex branch                   32.15695166         16.07847583
Composition, original sharper parameters   27.15403884         13.57701942
Composition, improved published inputs     27.13326883         13.56663441
Additional recent rectangular preprint     27.01387742         13.50693871

The published-input row is 15.62238513% below the supplied simplex leading
time coefficient. This is a comparison of asymptotic upper-bound formulas,
not a measured speedup. The logarithmic saving power remains 1/4. No
2^(~O(sqrt(b))) algorithm, stronger power, or practical crossover is claimed.

What changed
------------
1. The C38 Fourier words have shared endpoints and support canonical
   simplex recoding. An explicit auxiliary-wire estimate survives the
   growing alphabet penalty. The effective relation coefficient becomes
   19K instead of 22K from the C44 simplex application.
2. The new pair-sector join implementation chooses among the 280 balanced
   partitions and three internal choices. It performs exact algebra using
   an explicit Strassen backend and plugs into the supplied simplex engine.
3. The proof extends the sharper schedule uniformly to t growing as
   loglog(b), charging factorial normalization at every call, coefficient
   dimension, auxiliary wires, rounding, degree sectors and input access.
4. Exact support certificates admit alpha=0.32133 and Omega=2.37134 from
   published matrix results. A separate row uses the recent preprint.

A failed shortcut was retained as a counterexample: changing matrix
parameters can move a support-function maximum into the interval interior.
The old endpoint-only monotonicity argument then fails. The new certificate
maximizes on the full interval and records an exact violating interior point.

Run
---
Python 3.10+ and its standard library are sufficient. From this folder:

    python run_all.py

Every run gets a fresh runs/<id>/ directory with stdout, stderr, exact
records, source hashes, and an execution summary. latest_run.json identifies
the most recent run. Assertions are mandatory; standalone scripts refuse
Python -O. No network access or external package is needed to rerun checks.

Files
-----
proof.tex                     Derivation, named premises and resource audit.
c38_simplex.py                Adapter for actual C38 words and binary states.
adaptive_pair_join.py         New exact pair-sector join implementation.
verify_c38_simplex.py         Arithmetic recoding and corruption checks.
verify_adaptive_simplex.py    Independent sums and real recursive integration.
certify_bounds.py             Rational support, schedule and coefficient checks.
source_manifest.json          Pinned input identities and proof-source references.
vendor/                      Five unchanged source modules from supplied archives.
frontier.json                 Results, rejected route and next research priorities.

Evidence and boundaries
-----------------------
The arithmetic suite exercises 60 C38 instances, 36,608 Fourier words,
every reconstructed local predicate, tuple recovery, 285 corruptions, and
182 exhaustive small-word cases. It does not search for a relation by
contracting the entire C38 graph.

The join suite compares 48 finite mode sums with independent Cartesian
enumeration, checks 50 realized partition objectives, tests replicated
coordinate families and an actual nontrivial Strassen product over a
nilpotent ring, and checks all 70 active calls in two simplex integrations.

The certificate scans 30 matrix-parameter pairs and checks 60 finite
schedule cases, including integers adjacent to powers of two. The finite
schedule checks are precision/coordinate diagnostics; they do not enumerate
the asymptotic packing search. Whole-interval support inequalities come
from the algebraic argument and rational certificates, not sampled grids.

The executable join uses Strassen, while the smaller asymptotic constants
use published matrix-formula existence and the inherited charged catalog
construction. These are different evidence levels.

Next target
-----------
The next implementation milestone is to compile an actual C38 network
into executable local descriptors and connect contraction plus witness
recovery, first on explicitly pinned small instances. The next mathematical
milestone is a stronger saving power, requiring a change to the normalized
mode/precision spectrum or a representation that avoids full intermediate
tables. Alphabet optimization inside the existing selector model alone
cannot beat its limiting factor 1/2.

Original project archives and PDF are unchanged. Nothing was published.
