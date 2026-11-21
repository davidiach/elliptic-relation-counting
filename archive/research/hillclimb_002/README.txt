HILLCLIMB 002 -- EXPONENT TARGET, ENUMERATION BARRIER, EXECUTABLE C38 GRAPH
31 October 2026

Outcome
-------
The ECDLP goal remains 2^(tilde O(sqrt(b))) expected classical bit time.
No improved time coefficient or logarithmic saving power is claimed here.
The previous 27.13327 coefficient remains conditional on inherited proofs.

This iteration changes the research direction in two concrete ways:

1. The paid relation reduction does not turn a square-root-in-network-size
   contractor into a square-root-in-field-bitlength ECDLP algorithm.
   A contractor with exponent S^theta and graph size S ~ m*b gives an
   optimized ECDLP exponent with power 2*theta/(1+theta), subject to a
   uniform arity extension. theta=1/2 gives 2/3; theta=1/3 gives 1/2.
   The latter needs m ~ sqrt(b), outside the existing low-degree regime.
   A different reduction or specialized relation algorithm may avoid this.

2. An explicit family of feasible mode tuples has cardinality at least
   2^(M/(4*sqrt(u))), even with the total-degree cap. Thus unrestricted
   candidate enumeration retains the existing power. A stronger saving
   must exploit actual support, avoid enumeration, or change the linked
   spectrum/precision/geometry bounds. This is NOT a lower bound on ECDLP
   or on actual nonzero C38 coefficients.

Implemented and checked
-----------------------
- A general-prime C38 graph compiler using public data only, with real
  edge IDs, incident-only tensor entry access, and binary equality relays.
- 18 compiled graphs: 89,892 edges and 10,518 relays across fixtures.
  54 split-relay faults and 90 coherent Fourier mutations rejected.
- 12 exact recursive contractions and conditional witness recoveries,
  with all but one or two edges explicitly pinned. Recovered curve tuples
  were verified. These are interface checks, not an unpinned relation search.
- Exhaustive source searches for 12 m=2 instances over F_13 with hashes
  chosen independently of any witness: 26,364 source triples and 3,468
  independent point-pair checks. All 67 graph witnesses matched; three
  target instances had no decomposition.
- 60 independent exact rational Smith computations agree with the claimed
  finite scalar spectrum. Six large exact candidate-subset certificates
  verify the enumeration barrier's construction.

The foundation audit also exposes the scope of the inherited theorem:
as stated, it implies subexponential exact counting of arbitrary width-3
CNF formulas via a constant-degree tensor network. audit.tex gives the
reduction. This is an implication to scrutinize, not a disproof or a claim
that the inherited theorem has now been independently certified.

Reproduce
---------
Python 3.10+ and the standard library are sufficient. From this folder:

    python run_all.py

Each run records stdout/stderr, result JSON, code/proof and dependency
hashes in runs/<id>/. latest_run.json selects the latest checkpoint.
Root-level *_results.json files are the most recent standalone outputs.

Files
-----
audit.tex                       Mathematical audit and next hypotheses.
c38_network.py                  Public-data graph compiler and pinning API.
verify_network.py               Incidence, relay, contraction and recovery checks.
verify_source_enumeration.py     Complete small source search with curve oracle.
spectrum_audit.py               Independent Smith and enumeration certificates.
run_all.py                      Reproducible checkpoint runner.
frontier.json                   Accepted evidence, rejected routes, next actions.

Next research gate
------------------
Investigate exact representations of the actual C38 projected-child tables
that remain compact under recursion. A candidate needs proved construction,
access and composition costs; small observed rank is insufficient. Continue
auditing the broad contraction theorem independently. Any proposed stronger
contractor must pass the full ECDLP arity balance before being presented as
progress to the square-root target.

Original archives and PDF are unchanged. No files were published.
