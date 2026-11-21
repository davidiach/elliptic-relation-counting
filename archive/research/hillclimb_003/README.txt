HILLCLIMB 003 -- EXACT FACTORS AND ARITY-INDEPENDENT ROOT DEGREE
31 October 2026

Outcome
-------
Implemented an exact refactoring of the full general-prime C38 graph.
First-stage tests split into disjoint one- or two-source factors, closures
split by equation, and middle equations get explicit equality-copy wires.
Each equation uses at most eight sources; each source feeds at most six
equations. The new maximum-degree bound is ~O(sqrt(b)), without the old
factor m. Every added wire is counted and every old witness has exactly
one extension.

This removes the specific large-arity root-degree obstacle found in the
previous iteration. It does NOT improve the runtime exponent or certify
the other reduction assumptions at larger arities. The extra wires can
worsen the existing contraction constant; this is a structural research
step toward a different representation or contractor.

What the compression experiments found
-------------------------------------
- Natural Fourier cuts of actual C38 first tensors contain large identity
  minors. Their exact rank can be large even though the predicates have
  short descriptions. Six finite minors were evaluated, up to size 64.
- Actual projected children of a restricted C38 first tensor have the
  exact product formula with per-pair coefficient matrix

      [[1, -c(Z)], [-c(Z), 2*c(Z)^2]].

  Fixed-degree slices have provably exponential exact rank over the
  finite-precision jet ring, up to polynomial factors. Nevertheless the
  representation has only O(r) factors for r independent pairs.
- Including both exterior source tables gives the full join

      (4 - 4*c(Z) + 2*c(Z)^2)^r.

  Dropping those source moments gives the wrong answer. The corrected
  formula was checked against actual child contractions and, for small
  cases, direct deformation sums.

The rank bounds concern specified group/Fourier cuts. They do not exclude
better orderings, bases, circuits, or all compact representations. The
projected examples use explicit finite singleton packets, not the balanced
asymptotic packing schedule. They do not establish closure for the full
C38 recursion.

Evidence
--------
Compression probe: 340 actual target projected entries, 680 exterior child
contractions, six root minors, protected finite-precision rank slices.
Full refactoring: four graphs, 324 local-product comparisons, 48 complete
graph mutation comparisons, 32 rejected inconsistent copy branches.
Dependency counts checked at 99 arities, backed by the source-occurrence
proof in compression.tex.

Example (m=8, one finite fixture):
                        Original       Factored
Maximum degree             4,013            288
Edges                     13,754         24,890

Run
---
Python 3.10+ and the standard library are sufficient:

    python run_all.py

Each run is retained under runs/<id> with result JSON, stdout/stderr and
source hashes. latest_run.json identifies the latest checkpoint.

Files
-----
compression.tex          Proofs, scope qualifications and next hypothesis.
c38_slices.py            Actual first-tensor restrictions and Fourier minors.
probe_compression.py     Projected-table, finite-precision and factored-join probes.
factor_network.py        Full graph refactoring with all new equality wires.
verify_factoring.py      Equivalence and resource-ledger checks.
frontier.json            Research decisions and unresolved work.

Next gate
---------
Find a structural invariant for the refactored C38 factors that stays
small through mode projection, coefficient extraction and exact joins.
Pay for construction and access. A compact formula at one node or small
observed numerical rank is insufficient. Apply the complete arity balance
before claiming progress to the square-root-in-b ECDLP exponent.

The 27.13327 coefficient remains a conditional baseline. The goal remains
active and unachieved. Original archives and PDF are unchanged.
